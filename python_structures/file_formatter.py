###############################################################################
# Date: 17/01/22
# Author: Sebastian Porras
# Aims: These functions perform formatting on all input file for a request
# to the GRASP server. Visit to commands.py to see their usage.
###############################################################################


def _find_p(s: str):

    start = s.find('(')

    end = s[::-1].find(')')

    if start == -1 and end == -1:
        return start, end

    true_end = len(s) - end - 1

    return start, true_end


def _find_comma(s: str, level: int = 0) -> list[int]:

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


def _make_label(child, idx, count):

    if child[idx+1:].split(":")[0] == "":
        label = "N" + str(count) + child[idx+1:]

    else:
        label = child

    return label


def _nwk_split(n: str, visited=None, count=None, anc_count=None) -> dict:

    # keep track of order of visited nodes in DFS
    if visited is None:
        visited = dict()
        visited["order"] = []

    if count is None:
        count = 0

    # ancestor count which allows labeling
    if anc_count is None:
        anc_count = []

    start, end = _find_p(n)

    # base case for extant
    if start == -1 and end == -1:
        parent = n

    # case for ancestor node
    else:
        parent = "N" + str(count) + n[end+1:]
        anc_count.append(count)
        count += 1

    visited["order"].append(parent)

    children = n[start+1: end]

    # Splits up subtree into children

    coms = _find_comma(children)

    sub_strings = []

    prev_com = 0

    # record children of this parent
    if len(coms) == 0:

        return visited

    for i in coms:

        sub_strings.append(children[prev_com: i])
        prev_com = i + 1

    sub_strings.append(children[coms[-1]+1:])
    ##################################

    # now check these children for their children
    for child in sub_strings:

        _, p_end = _find_p(child)

        if count not in anc_count:

            label = _make_label(child, p_end, count)

            visited[label] = parent

            _nwk_split(child, visited, count, anc_count)

        # this loop here avoids count reverting back
        # to its value further up the tree
        else:
            while (count in anc_count):
                count += 1

            label = _make_label(child, p_end, count)

            visited[label] = parent

            _nwk_split(child, visited, count, anc_count)

    return visited


def _nwk_split_labels(n: str, visited=None, count=None, anc_count=None) -> dict:

    # keep track of order of visited nodes in DFS
    if visited is None:
        visited = dict()
        visited["order"] = []

    start, end = _find_p(n)

    # base case for extant
    if start == -1 and end == -1:
        parent = n

    # case for ancestor node
    else:
        parent = n[end+1:]

    visited["order"].append(parent)

    children = n[start+1: end]

    # Splits up subtree into children

    coms = _find_comma(children)

    sub_strings = []

    prev_com = 0

    # record children of this parent
    if len(coms) == 0:

        return visited

    for i in coms:

        sub_strings.append(children[prev_com: i])
        prev_com = i + 1

    sub_strings.append(children[coms[-1]+1:])
    ##################################

    # now check these children for their children
    for child in sub_strings:

        _, p_end = _find_p(child)

        label = child[p_end+1:]

        visited[label] = parent

        _nwk_split_labels(child, visited, count, anc_count)

    return visited


def nwkToJSON(nwk: str) -> dict:
    """Assumes that internal nodes have not been added"""

    # find the first ')' to locate the root
    end = len(nwk) - nwk[::-1].find(")")

    # no internal node labels or root distance
    if nwk[end:] == ";":

        # remove ; and add root distance
        job = nwk[:-1] + ":0"

        raw_tree = _nwk_split(job)

    # root has dist but no internal nodes e.g. ...):0.123;
    elif nwk[end:].split(":")[0] == "":

        job = nwk[:end] + nwk[end:-1]

        raw_tree = _nwk_split(job)

    # for nwk output from a GRASP job e.g ...)N0;
    elif nwk[end:].split(":")[0] != "N0":

        job = nwk[:-1] + ":0"

        raw_tree = _nwk_split_labels(job)

    else:
        raise RuntimeError("nwk in unsupported or incorrect format")

    idxs = {name: i for i, name in enumerate(raw_tree["order"])}
    Parents = []
    Labels = []
    Distances = []

    for name in raw_tree["order"]:

        Labels.append(name.split(":")[0])

        Distances.append(float(name.split(":")[1]))

        if name.split(":")[0] == "N0":
            Parents.append(-1)
        else:
            Parents.append(idxs[raw_tree[name]])

    json_idx = dict()
    json_idx["Parents"] = Parents
    json_idx["Labels"] = Labels
    json_idx["Distances"] = Distances
    json_idx["Branchpoints"] = len(Labels)

    return json_idx


def readAln(file_name: str, data_type: str) -> dict:
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
                tmp["Seq"] = [letter for letter in aln]

                seqs.append(tmp)
            else:
                line = fa.readline()

    alignments = dict()
    alignments["Sequences"] = seqs
    alignments["Datatype"] = {"Predef": data_type}

    return alignments
