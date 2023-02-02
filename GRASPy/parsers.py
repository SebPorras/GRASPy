###############################################################################
# Date: 17/01/22
# Author: Sebastian Porras
# Aims: Any scripts related to translating one form of data into another.
# Most functions convert to and from JSON format into the various Classes
# and also any string formating functions.
###############################################################################

from typing import Tuple, Union, Optional
from . import pog_tree
from . import pog_graph
import numpy as np
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

        # this is the base case for an extant sequence
        if p_start == -1 and p_end == -1:

            data["order"].append(child)

            # record the parent of this current child
            data["Parents"][child] = parent

        else:

            # implies another ancestor has been found
            data["count"] += 1

            label = make_label(child, p_end, data["count"])

            data["Parents"][label] = parent

            # repeat on the next ancestor
            nwk_split(child, data)

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

        raw_tree = nwk_split(job)

    # for nwk output from a GRASP job e.g ...)N0;
    elif len(nwk[end:].split(":")) == 1:

        job = nwk[:-1] + ":0"

        raw_tree = nwk_split(job)

    # for nwk output from GRASPy e.g ...)N0:0.0;
    elif len(nwk[end:].split(":")) == 2:

        job = nwk[:-1]

        raw_tree = nwk_split(job)

    else:
        raise RuntimeError("nwk in unsupported or incorrect format")

    idxs = {name: i for i, name in enumerate(raw_tree["order"])}

    Parents = []
    Labels = []
    Distances = []

    for i, node in enumerate(raw_tree["order"]):

        Labels.append(node.split(":")[0])
        Distances.append(float(node.split(":")[1]))

        if i == 0:
            p_idx = -1
        else:
            # grab the index of the parent for this node
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

    seqs = []

    for seq in sequences:

        tmp = {}

        tmp["Name"] = seq.name

        # remove any gap characters
        tmp["Seq"] = [None if s == "-" else s for s in seq.sequence]

        seqs.append(tmp)

    alignments = dict()
    alignments["Sequences"] = seqs

    # will guess based on sequence what the alphabet is
    if data_type is None:
        alignments["Datatype"] = {"Predef": sequences[0].alphabet.name}
    else:
        alignments["Datatype"] = {"Predef": data_type}

    return alignments


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

    # index by child branch point index and map their parent index
    parents = [jparents[i] for i in range(nBranches)]

    # using same index as branch point, map distances
    distances = [jdists[i] for i in range(nBranches)]

    indices = dict()

    for i in range(nBranches):

        # index by sequence id and map branch point index
        lab = jlabels[i]

        if lab.isdigit():
            lab = "N" + lab

        indices[lab] = i

    # Next step is to record children of each parent
    children = []

    for PIdx in range(nBranches):

        curr_children = []

        # reference parent array and check if it matches current parent index
        for CIdx in range(nBranches):

            if (parents[CIdx] == PIdx):
                curr_children.append(CIdx)

        # Leaves will have no children which is recorded accordingly
        if len(curr_children) == 0:
            curr_children = [None]

        children.append(curr_children)

    # branch points represent sequences on the tree
    bpoints = {}

    for BIdx in range(nBranches):

        # marks ancestor labels with an 'N' to avoid confusion
        name = jlabels[BIdx]
        if name.isdigit():

            name = "N" + name

        # convert parent indexes to seq names
        c_lab = []

        if None in children[BIdx]:
            c_lab.append(None)

        else:
            for i in children[BIdx]:

                lab = jlabels[i]

                if lab.isdigit():

                    lab = "N" + lab

                c_lab.append(lab)

        # convert parent indexes to seq names
        p_idx = parents[BIdx]

        if p_idx == -1:
            p = None

        else:
            p = "N" + jlabels[p_idx]

        bp = pog_tree.BranchPoint(id=name,
                                  parent=p,
                                  dist=distances[BIdx],
                                  children=c_lab)

        bpoints[name] = bp

    tree = {}
    tree['nBranches'] = nBranches
    tree['branchpoints'] = bpoints
    tree['parents'] = parents
    tree['children'] = children
    tree['indices'] = indices
    tree['distances'] = distances

    return tree


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

    # adjacent edges start at each index and end at value stored there
    adj = jpog["Adjacent"]

    if len(indices) != len(adj):
        raise RuntimeError("JSON is incorrectly formatted")

    # contains outgoing edges and the most likely residue
    node_vals = jpog["Nodes"]
    nodes = []

    for i in range(len(indices)):

        es = []

        # set [] to -999 to signify end of seq
        if len(adj[i]) == 0:
            edge = pog_graph.Edge(start=indices[i], end=-999)
            es.append(edge)
        else:
            # there can be multiple adjacent edges for each sequence position
            for j in range(len(adj[i])):

                edge = pog_graph.Edge(start=indices[i], end=adj[i][j])
                es.append(edge)

        node = pog_graph.SymNode(name=indices[i],
                                 symbol=node_vals[i]["Value"], edges=es)

        nodes.append(node)

    if isAncestor:

        # tuple of start and end of each edge
        edgeInd = jpog["Edgeindices"]

        # extra annotations about each edge
        edges = jpog["Edges"]

        for i in range(len(jpog["Edgeindices"])):

            edge_info = edges[i]

            # virtual start is always referenced as -1
            if edgeInd[i][0] == -1:

                # virtual start edges will be assigned to the first "real" node
                cor_node = nodes[0]

                edge = pog_graph.Edge(start=edgeInd[i][0], end=edgeInd[i][1],
                                      edgeType=jpog["Edgetype"],
                                      recip=edge_info["Recip"],
                                      backward=edge_info["Backward"],
                                      forward=edge_info["Forward"],
                                      weight=edge_info["Weight"])

                cor_node.edges.append(edge)

            else:

                # edges are stored in random order, need to find correct start

                # returns a tuple with an array containg correct index
                node_loc = np.where(indices == edgeInd[i][0])[0][0]

                cor_node = nodes[node_loc]

                # some edges are identical to adj edges and can be replaced
                dupeedges = []

                for j in range(len(cor_node.edges)):

                    adjacent_edge = cor_node.edges[j]

                    if edgeInd[i][0] == adjacent_edge.start and\
                            edgeInd[i][1] == adjacent_edge.end:

                        dupeedges.append(j)

                # remove any identified duplicates
                [cor_node.edges.pop(j) for j in dupeedges]

                edge = pog_graph.Edge(start=edgeInd[i][0],
                                      end=edgeInd[i][1],
                                      edgeType=jpog["Edgetype"],
                                      recip=edge_info["Recip"],
                                      backward=edge_info["Backward"],
                                      forward=edge_info["Forward"],
                                      weight=edge_info["Weight"])

                cor_node.edges.append(edge)

    # adds ancestor identifier "N"
    name = jpog["Name"]

    if name.isdigit():
        name = "N" + name

    return pog_graph.POGraph(version=jpog["GRASP_version"], indices=indices,
                             start=jpog["Starts"], end=jpog["Ends"], size=jpog["Size"],
                             terminated=jpog["Terminated"], directed=jpog["Directed"],
                             name=name, isAncestor=isAncestor, nodes=nodes)


def POGTreeFromJSON(nwk: Union[str, dict], POG_graphs: dict) -> pog_tree.POGTree:
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

            e = POGraphFromJSON(e, isAncestor=False)

            graphs[e.name] = e

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
