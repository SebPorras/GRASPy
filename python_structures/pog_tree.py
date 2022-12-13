###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: The POG tree will have two components, 1) the Idx tree which encodes
# the topology of the tree. 2) The POG graph that will be stored according
# to the corresponding index on the IDX tree so that both structures
# are in agreement with each other.
#################################################################################

import json
from pog_graph import *
from idx_tree import *


class POGTree(object):
    """Contains an IdxTree and a POGraph for each branch point.
    All graphs and branchpoints are assigned the same index for
    continuity
    """

    def __init__(self, idxTree: IdxTree, POGraphs: dict) -> None:
        """Constructs instance of POGTree

        Parameters:
            idxTree(IdxTree)
            POGraph(POGraph)
        """

        self._idxTree = idxTree
        self._graphs = POGraphs

    def writeNwk(self, file_name: str) -> None:
        """Writes the tree including ancestors in the 
        Newick Standard (nwk) format  

        Parameters: 
            file_name(str): specify the file name 
            and file path for output
        """

        if file_name[-4:] != ".nwk":
            file_name += ".nwk"

        nwk = self._idxTree._createNwk()

        with open(file_name, "w") as out:
            out.write(nwk)


def POGTreeFromJSON(json_path: str) -> POGTree:
    """Instantiates a POGraph object from a JSON file.

    Parameters:
        json_path (json): path to JSON file
    """

    with open(json_path, "r") as file:
        data = json.load(file)

    tree = IdxTreeFromJSON(data)

    extants = data["Input"]["Extants"]
    ancestors = data["Ancestors"]

    graphs = dict()

    for e in extants:

        g = POGraphFromJSON(e, isAncestor=False)

        idx = tree.getIndexOf(g._name)

        graphs[idx] = g

    for a in ancestors:

        g = POGraphFromJSON(a, isAncestor=True)

        idx = tree.getIndexOf(g._name)

        graphs[idx] = g

    return POGTree(idxTree=tree, POGraphs=graphs)


if __name__ == "__main__":

    poggers = POGTreeFromJSON("./python_structures/ASR_big.json")

    test = poggers.writeNwk("tester.nwk")

    # for n in poggers._graphs[0]._nodes:
    #     if len(n._edges) >= 2:
    #         for e in n._edges:
    #             print(e)
