###############################################################################
# Date: 17/01/22
# Author: Sebastian Porras
# Aims: These functions perform formatting on all input file for a request
# to the GRASP server. Visit to commands.py to see their usage.
###############################################################################
from typing import Tuple


def _find_p(s: str) -> Tuple[int, int]:
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


def _find_comma(s: str, level: int = 0) -> list[int]:
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


def _make_label(child: str, isLabeled: bool, idx: int, count: int) -> str:
    """
    Used within the _nwk_split method to create labels for nodes.

    Parameters:
        child(str): current child 
        isLabeled(bool): if the child has a label or not 
        idx(int): marks where the child label begins
        count(int): the current ancestor node count 

    Returns:
        str: The label of the node 
    """

    if isLabeled:
        label = str(count) + \
            ":" + child[idx+1:].split(":")[1]

    else:
        label = str(count) + child[idx+1:]

    return label


def _nwk_split(n: str, data=None, isLabeled=False) -> dict:
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
    start, end = _find_p(n)

    # Assigns an internal node a number starting with 0 at the root

    parent = _make_label(n, isLabeled, end, data["count"])

    data["order"].append(parent)

    # Grabs all children of this ancestor
    children = n[start+1: end]

    # finds all children for current embedding
    coms = _find_comma(children)

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
        p_start, p_end = _find_p(child)

        # this is the base case for an extant sequence
        if p_start == -1 and p_end == -1:

            data["order"].append(child)

            # record the parent of this current child
            data["Parents"][child] = parent

        else:

            # implies another ancestor has been found
            data["count"] += 1

            label = _make_label(child, isLabeled, p_end, data["count"])

            data["Parents"][label] = parent

            # repeat on the next ancestor
            _nwk_split(child, data, isLabeled)

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

    # no internal node labels or root distance
    if nwk[end:] == ";":

        # remove ;  and add root distance
        job = nwk[:-1] + ":0"

        raw_tree = _nwk_split(job)

    # root has dist but no internal nodes e.g. ...):0.123;
    elif nwk[end:].split(":")[0] == "":

        job = nwk[:end] + nwk[end:-1]

        raw_tree = _nwk_split(job)

    # for nwk output from a GRASP job e.g ...)N0;
    elif nwk[end:].split(":")[0] != "":

        job = nwk[:-1]

        raw_tree = _nwk_split(job, isLabeled=True)

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


def alnToJSON(file_name: str, data_type: str) -> dict:
    """
    Creates a dictionary where seq ids are the key
    and alignment is the value.

    Parameters:
        file_name(str): path to aln file
        data_type(str): user must specify what alignment letters are.
        E.g DNA or Protein

    Returns:
        list: alignment seqs in JSON format

    """

    seqs = []

    with open(file_name, "r") as fa:

        line = fa.readline()

        while line:

            if line[0] == ">":

                tmp = {}

                # this would save extra info in fasta file
                # key = "_".join(line.split())

                # currently, just saves seq id and removes '>'
                key = line.split()[0][1:]

                aln = ""

                line = fa.readline()

                while line[0] != ">":

                    aln += line.strip()
                    line = fa.readline()

                    if not line:
                        break

                tmp["Name"] = key

                sequence = []

                for letter in aln:
                    if letter == "-":
                        sequence.append(None)
                    else:
                        sequence.append(letter)

                tmp["Seq"] = sequence

                seqs.append(tmp)
            else:
                line = fa.readline()

    alignments = dict()
    alignments["Sequences"] = seqs
    alignments["Datatype"] = {"Predef": data_type}

    return alignments
