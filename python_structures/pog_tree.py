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
    """The Partial Order Graph Tree (POGTree), is a phylogenetic tree that 
    contains two data structures. 1) An IdxTree object which holds
    information about the topology and information about the sequence at 
    each branchpoint of the tree. 2) A POGraph object which describes
    the graph of the sequence at that branchpoint.
    """

    def __init__(self, idxTree: IdxTree, POGraphs: dict[str, POGraph]) -> None:
        """Constructs instance of POGTree

        Parameters:
            idxTree(IdxTree)
            POGraph(POGraph)
        """

        self.idxTree = idxTree
        self.graphs = POGraphs

    def pogToNwk(self, root: str = "N0",) -> str:
        """Writes the tree including ancestors in the
        Newick Standard (nwk) format. Root is based on
        sequence ID.

        Parameters:
            root(str): Default set to N0 at the "root" ancestor
            but can be changed to create subtrees if desired. 

        Returns:
            str: The POGTree in nwk format 
        """

        nwk = ""

        if None in self.idxTree.branchpoints[root].children:

            nwk += f"{self.idxTree.branchpoints[root].id}:{self.idxTree.branchpoints[root].dist}"

        else:

            for c in self.idxTree.branchpoints[root].children:
                nwk += self.pogToNwk(self.idxTree.branchpoints[c].id) + ','

            # slicing removes final ',' at end of subtree
            nwk = "(" + nwk[:-1] + \
                f"){self.idxTree.branchpoints[root].id}:{self.idxTree.branchpoints[root].dist}"

        return nwk

    def writeNwk(self, file_name: str, root: str = "N0") -> str:
        """Writes a nwk string of the tree to a file
        
        Parameters: 

            file_name(str): name of nwk file

            root(str): Default set to N0 at the "root" ancestor
            but can be changed to create subtrees if desired. 
        """

        nwk = self.pogToNwk(root) + ';'

        with open(file_name, 'w') as f:
            f.write(nwk)

        return nwk


def POGTreeFromJSON(json_path: str) -> POGTree:
    """Instantiates a POGTree object from a JSON file.

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

        graphs[g.name] = g

    for a in ancestors:

        g = POGraphFromJSON(a, isAncestor=True)

        graphs[g.name] = g

    return POGTree(idxTree=tree, POGraphs=graphs)


if __name__ == "__main__":

    poggers = POGTreeFromJSON("./test_data/small_test_data/ASR.json")
    print(poggers.idxTree.branchpoints["N0"].children)
    # print(poggers.pogToNwk())

    # print(poggers.graphs['XP_004050792.2'].version)
    # for i in poggers.idxTree.branchpoints:
    #     # if None in i.children:
    #     print(i)
    #     print()

  