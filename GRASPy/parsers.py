###############################################################################
# Date: 17/01/22
# Author: Sebastian Porras
# Aims: Any scripts related to translating one form of data into another.
# Most functions convert to and from JSON format into the various classes
# and also any string formating functions.
###############################################################################

from typing import Tuple, Union, Optional
from . import pog_tree
from . import pog_graph
import numpy as np
import pandas as pd
from . import sequence


def find_p(s: str) -> Tuple[int, int]:
    """Locates the first index positions for the top most
    embedding within a nwk string. This allows substrings
    to be located
    For example:

    (A:1,(B:2)C:3)D:4; -> 0, 13
    """

    start = s.find('(')

    end = s[::-1].find(')')

    if start == -1 and end == -1:
        return start, end

    true_end = len(s) - end - 1

    return start, true_end


def find_comma(s: str, level: int = 0) -> list[int]:
    """
    Finds all commas within a particular embedding level for
    a nwk string. This essentially locates where children are
    to allow subsetting of the string.

    Parameters:
        s(str): nwk string
        level(int): current embedding level

    Returns:
        list: contains indices for each comma within
        the embedding
    """

    mylevel = 0
    coms = []

    for i in range(len(s)):

        if s[i] == '(':
            mylevel += 1

        elif s[i] == ')':
            mylevel -= 1

        elif s[i] == ',' and mylevel == level:
            coms.append(i)

    return coms


def make_label(child: str, idx: int, count: int) -> str:
    """
    Used within the nwk_split method to create labels for nodes.

    Parameters:
        child(str): current child
        idx(int): marks where the child label begins
        count(int): the current ancestor node count

    Returns:
        str: The label of the node
    """

    return (str(count) + ":" + child[idx+1:].split(":")[1])


def nwk_split(n: str, data=None) -> dict:
    """Performs a depth first search (DFS) of a nwk string and records
    the order of sequences and the parent of each node. Note the
    recursive nature of the function.

    Parameters:
        n(str): The nwk string
        data(None): will record tree info during DFS
        isLabeled(bool): marks if internal nodes have
        been assigned labels or not.

    Returns:
        dict: contains order and parents for each node
    """

    if data is None:
        data = dict()
        data["order"] = []
        data["count"] = 0
        data["Parents"] = {}

    # locates the first embedding in a nwk string
    start, end = find_p(n)

    # Assigns an internal node a number starting with 0 at the root

    parent = make_label(n, end, data["count"])

    data["order"].append(parent)

    # Grabs all children of this ancestor
    children = n[start+1: end]

    # finds all children for current embedding
    coms = find_comma(children)

    sub_strings = []

    prev_com = 0

    # save each child and any children it has
    for i in coms:

        sub_strings.append(children[prev_com: i])
        prev_com = i + 1

    sub_strings.append(children[coms[-1]+1:])

    # now check these children for their own children
    for child in sub_strings:

        # look for the next embedding within "(...)"
        p_start, p_end = find_p(child)

        if p_start != -1 and p_end != -1:

            # implies another ancestor has been found
            data["count"] += 1

            data["Parents"][make_label(child, p_end, data["count"])] = parent

            # repeat on the next ancestor
            nwk_split(child, data)

            continue

        data["order"].append(child)

        # record the parent of this current child
        data["Parents"][child] = parent

    return data


def nwkToJSON(nwk: str) -> dict:
    """Parses a nwk string into a JSON format.
    This method assumes that a tree either has no names for internal
    nodes OR that the internal nodes have been named according to
    GRASP's naming system i.e. N0, N1 etc.

    Parameters:
        nwk(str): A nwk string

    Returns:
        dict: Contains a representation of an IdxTree in
        JSON format.
    """

    # find the first ')' to locate the root
    end = len(nwk) - nwk[::-1].find(")")

    # standard nwk format e.g ...);
    if nwk[end:] == ";":

        # remove ;  and add root distance
        job = nwk[:-1] + ":0"

    # for nwk output from a GRASP job e.g ...)N0;
    elif len(nwk[end:].split(":")) == 1:

        job = nwk[:-1] + ":0"

    # for nwk output from GRASPy e.g ...)I:0.2;
    elif len(nwk[end:].split(":")) == 2:
        # job = nwk[:-1]
        job = nwk[:end] + ":0"

    else:
        raise RuntimeError("nwk in unsupported or incorrect format")

    raw_tree = nwk_split(job)

    # record what the indices that will be assigned to each node
    idxs = {name: i for i, name in enumerate(raw_tree["order"])}

    Parents = []
    Labels = []
    Distances = []

    # add information to idx tree in a DFS order
    for i, node in enumerate(raw_tree["order"]):

        # grab name for node
        Labels.append(node.split(":")[0])
        # grab distance of that node
        Distances.append(float(node.split(":")[1]))

        # record index of the parent, -1 for root parent
        p_idx = -1

        if i != 0:
            p_idx = idxs[(raw_tree["Parents"][node])]

        Parents.append(p_idx)

    json_idx = dict()
    json_idx["Parents"] = Parents
    json_idx["Labels"] = Labels
    json_idx["Distances"] = Distances
    json_idx["Branchpoints"] = len(Labels)

    return json_idx

def alnToJSON(file_name: str, data_type: Optional[str] = None) -> dict:
    """
    Creates a dictionary where seq ids are the key
    and alignment is the value.

    Parameters:
        file_name(str): path to aln file
        data_type(str): user must specify what alignment letters are
        e.g DNA or Protein otherwise it will guess

    Returns:
        list: alignment seqs in JSON format

    """

    sequences = sequence.readFastaFile(file_name)

    j_seqs = []

    for seq in sequences:

        tmp = {}

        tmp["Name"] = seq.name

        # remove any gap characters
        # TODO: Note from Sam - this treats gaps as missing data rather than gaps. Should use gappy alphabet instead (see methods for modes)?
        tmp["Seq"] = [None if s == "-" else s for s in seq.sequence]

        j_seqs.append(tmp)

    alignments = dict()
    alignments["Sequences"] = j_seqs

    # will guess based on sequence what the alphabet is
    if data_type is None:
        data_type = sequences[0].alphabet.name

    alignments["Datatype"] = {"Predef": data_type}

    return alignments


def make_anc_label(jlabels: dict, i: int) -> str:

    lab = jlabels[i]

    if lab.isdigit():
        lab = "N" + lab

    return lab


def record_children(PIdx: int, parents: list, nBranches: int):

    curr_children = [CIdx for CIdx in range(
        nBranches) if parents[CIdx] == PIdx]

    # Leaves will have no children which is recorded accordingly
    if len(curr_children) == 0:
        curr_children = [None]

    return curr_children


def TreeFromJSON(serial: dict) -> dict:
    """Creates a IdxTree in a dict from a JSON file

    Parameters:
        serial: JSON form of an IdxTree

    Returns:
        dict
    """

    if ("Distances" or "Branchpoints" or "Labels" or "Parents") \
            not in serial.keys():

        raise RuntimeError("JSON in incorrect format")

    jdists = serial["Distances"]
    nBranches = serial["Branchpoints"]
    jlabels = serial["Labels"]
    jparents = serial["Parents"]

    # map the names of b_points to their index in the tree
    indices = {make_anc_label(jlabels, i): i for i in range(nBranches)}

    # keep track of the children of each b_point
    children = [record_children(PIdx, jparents, nBranches)
                for PIdx in range(nBranches)]

    bpoints = {}

    # branch points represent nodes on a tree
    for BIdx in range(nBranches):

        branch_name = make_anc_label(jlabels, BIdx)

        branch_children = [None if i is None else make_anc_label(
            jlabels, i) for i in children[BIdx]]

        branch_PIdx = jparents[BIdx]

        parent_name = None

        # -1 is for the root node which has None as its parent
        if branch_PIdx != -1:
            parent_name = make_anc_label(jlabels, branch_PIdx)

        bpoints[branch_name] = pog_tree.BranchPoint(id=branch_name,
                                                    parent=parent_name,
                                                    dist=jdists[BIdx],
                                                    children=branch_children)

    tree = {}
    tree['nBranches'] = nBranches
    tree['branchpoints'] = bpoints
    tree['parents'] = jparents
    tree['children'] = children
    tree['indices'] = indices
    tree['distances'] = jdists

    return tree


def makeEdges(indices, adj, i) -> list[pog_graph.Edge]:
    '''Creates a list of Edge objects for a particular node'''

    # set [] to -999 to signify end of seq
    if len(adj[i]) == 0:
        return [pog_graph.Edge(start=indices[i], end=-999)]

    # add all adjacent edges, note that ancestors can have multiple adjacent edges
    return [pog_graph.Edge(start=indices[i], end=adj[i][j]) for j in range(len(adj[i]))]


def locateEdgeIndex(edges, i, indices) -> int:
    '''Finds location of the node where the edge begins'''

    # check if edge starts at -1 -> virtual start
    if edges[i][0] == -1:
        # this edge will be added to first real node
        return 0

    # edges are stored in random order, find index where edge is
    return np.where(indices == edges[i][0])[0][0]


def removeDuplicateEdges(edge: list[int], node: pog_graph.SymNode):
    '''Ancestors have extra edges that aren't simply adjacent. Some 
    of the extra edges are identical to the adjacent edges and so 
    this redundancy is removed'''

    # check if there's already an edge with the same start and end
    dupeedges = [j for j in range(len(node.edges)) if
                 (edge[0] == node.edges[j].start and
                  edge[1] == node.edges[j].end)]

    [node.edges.pop(j) for j in dupeedges]


def addAncestralEdges(edgeIndices: list[list[int]], edge_info: dict,
                      edgeType: str, indices, nodes: list[pog_graph.SymNode]):

    START = 0
    END = 1

    for i in range(len(edgeIndices)):

        # find position of node where edge begins
        edge_location = locateEdgeIndex(edgeIndices, i, indices)

        ancestor_node = nodes[edge_location]

        # some edges are identical to adj edges and can be replaced
        removeDuplicateEdges(edgeIndices[i], ancestor_node)

        edge = pog_graph.Edge(start=edgeIndices[i][START],
                              end=edgeIndices[i][END],
                              edgeType=edgeType,
                              recip=edge_info[i]["Recip"],
                              backward=edge_info[i]["Backward"],
                              forward=edge_info[i]["Forward"],
                              weight=edge_info[i]["Weight"])

        ancestor_node.edges.append(edge)


def POGraphFromJSON(jpog: dict, isAncestor: bool = False) -> pog_graph.POGraph:
    """Takes JSON format of a partial order graph and transcribes
    this information into a POGraph object.

    Parameters:
        jpog(dict): serialised JSON format of a POG for a sequence

        isAncestor(bool): Only ancestors have multiple edges unlike
        extants which only have adjacent edges.
    """

    # each position in a sequence is given an index based on the alignment
    indices = np.array(jpog["Indices"])

    # the index is the start of the edge and the value is the index of the end
    adj = jpog["Adjacent"]

    # contains the most likely residue at that position
    node_vals = jpog["Nodes"]

    # store node value and all adjacent edges with this node
    nodes = [pog_graph.SymNode(name=indices[i],
                               symbol=node_vals[i]["Value"],
                               edges=makeEdges(indices, adj, i)) for i in range(len(indices))]

    if isAncestor:

        addAncestralEdges(jpog["Edgeindices"], jpog["Edges"],
                          jpog["Edgetype"], indices, nodes)

    # adds ancestor identifier "N"
    name = jpog["Name"]

    if name.isdigit():
        name = "N" + name

    return pog_graph.POGraph(version=jpog["GRASP_version"], indices=indices,
                             start=jpog["Starts"], end=jpog["Ends"], size=jpog["Size"],
                             terminated=jpog["Terminated"], directed=jpog["Directed"],
                             name=name, isAncestor=isAncestor, nodes=nodes)


def POGTreeFromJointReconstruction(nwk: Union[str, dict], POG_graphs: dict) -> pog_tree.POGTree:
    """Creates an instance of the POGTree data structure. A nwk
    file OR output from g_requests.requestPOGTree() can be used
    to create tree topology with the second option also creating
    POGraphs for extants.

    Parameters:
        nwk (str or dict): Users can input a nwk file path
        or can provide the output from g_requests.requestPOGTree().


        POG_graphs(dict): The POGraphs for ancestors generated from
        output from g_requests.requestJointReconstruction().

    Returns:
        POGTree
    """
    # this will hold all POGraphs
    graphs = dict()

    # case for using a nwk file for IdxTree
    if isinstance(nwk, str):

        with open(nwk, 'r') as f:
            tree_parts = ""
            for line in f:
                tree_parts += line.strip()

            j_tree = nwkToJSON(tree_parts)

        tree = TreeFromJSON(j_tree)

    # case when using output from server for IdxTree
    elif isinstance(nwk, dict):

        tree = TreeFromJSON(nwk["Result"]["Tree"])

        for e in nwk["Result"]["Extants"]:

            ex = POGraphFromJSON(e, isAncestor=False)

            graphs[ex.name] = ex

    else:
        raise RuntimeError("Nwk tree in unsupported format")

    # Translate the POGraphs JSON and save to the tree
    ancestors = POG_graphs["Result"]["Ancestors"]

    for a in ancestors:

        g = POGraphFromJSON(a, isAncestor=True)

        graphs[g.name] = g

    return pog_tree.POGTree(nBranches=tree['nBranches'],
                            branchpoints=tree['branchpoints'],
                            parents=tree['parents'],
                            children=tree['children'],
                            indices=tree['indices'],
                            distances=tree['distances'],
                            POGraphs=graphs)


def csvDataToJSON(file_name: str) -> dict[str, list]:
    """Reads in datafile and formats in the correct format for use
    in LearnLatentDistributions and MarginaliseDistOnAncestor.

    ONLY implemented for continuous values not discrete multinomial
    symbols.

    Bnkit requires that data must be in a square matrix, therefore
    the number of observations must be uniform and null values are
    in the place of "missing" observations.

    Column names in CSV MUST be called "Headers" & "Data"

    If there are multiple observations for an annotation,
    separate the observations with whitespace within a cell
    e.g. "7.0 3.3"
    """

    # None get converted to null for JSON
    raw = pd.read_csv(file_name).replace(pd.NA, None)

    j_data = dict()

    # Headers represent names of sequences
    j_data["Headers"] = raw["Headers"].tolist()

    # place observations within lists and convert multiple
    # observations into floats
    formatted = [[float(obs) for obs in annot.split()] if isinstance(annot, str)
                 else [annot] for annot in raw["Data"]]

    # find the largest list which represents the max number of datapoints
    max_len = len(max(formatted, key=lambda x: len(x)))

    # add null values to ensure that there is uniform length
    for obs in formatted:
        if len(obs) < max_len:
            while len(obs) < max_len:
                obs.append(None)

    j_data["Data"] = formatted

    return j_data


def motifDistribJSON(node_names : list[str],
                     n_latent : int,
                     alpha : str = 'Protein with gap'
                     ):
    """
    Construct the input JSON dict for motif mode training

    NOTE: Specific to motifs as data type for a single mode - need to generalise this further
    """

    distrib = dict()

    # Mapping of each real node to the set of target latent nodes
    # NOTE: In this specific case, all mapped to the single latent node
    targets = [[0] for node in node_names]
    distrib["Targets"] = targets

    # Properties of latent states (modes)
    # NOTE: Hard-coded for a single mode, with n_latent possible states
    state_vals = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    modetypes = [{"Size" : n_latent, "Values" : [state_vals[i] for i in range(n_latent)], "Datatype" : "Character"}]
    distrib["Modetypes"] = modetypes

    # Real nodes, indexed by node_names
    # NOTE: Hard-coded as CPTs for representing characters within a motif
    nodes = [{"Condition" : [],
              "Pr" : [],
              "Variable" : {"Domain" : {"Predef" : alpha}, "Name" : f"N0__{name}"},
              "Nodetype" : "CPT",
              "Index" : []
        } for name in node_names]
    distrib["Nodes"] = nodes

    # Distrib name matches prefix for node names - hard-coded as 'N0' here
    distrib["Name"] = "N0"

    return distrib


def modeDataToJSON(file_name : str) -> dict:  # TODO: update with specific dict data type
    """ Reads in data for evolutionary modes. 'Data' (observations) are indexed by 'item' (node label), and feature
    (real plate nodes associated with a given latent node).

    NOTE:
    Currently this outputs datasets in the format used for training and inferring evolutionary 'modes'. It is also
    only implemented for single observations per feature per item. TODO: I still don't understand how it works with multiple

    First column headers must be "Items", subsequent column headers specify Feature names
    """

    raw = pd.read_csv(file_name, index_col=0).replace(pd.NA, None)

    j_data = dict()

    # Items (node labels most likely)
    j_data["Items"] = raw.index.tolist()

    # All other cols are features
    j_data["Features"] = raw.columns.tolist()

    # TODO: currently assumes one observation per feature per item
    data = [[[raw[feature][item] for feature in j_data["Features"]] for item in j_data["Items"]]]
    j_data["Data"] = data

    return j_data
