###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: The POG tree will have two components, 1) the Idx tree which encodes
# the topology of the tree. 2) The POG graph that will be stored according
# to the corresponding index on the IDX tree so that both structures
# are in agreement with each other.
#################################################################################


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

    def _POGTreeToNwk(self, root: str = "N0",) -> str:
        """Converts the POGTree to nwk form Newick Standard (nwk) format.
        Root can be changed if the user wishes to create subtrees.

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
                nwk += self._POGTreeToNwk(
                    self.idxTree.branchpoints[c].id) + ','

            # slicing removes final ',' at end of subtree
            nwk = "(" + nwk[:-1] + \
                f"){self.idxTree.branchpoints[root].id}:{self.idxTree.branchpoints[root].dist}"

        return nwk

    def writeToNwk(self, file_name: str, root: str = "N0") -> str:
        """Writes a nwk string of the tree to a file

        Parameters: 

            file_name(str): name of nwk file

            root(str): Default set to N0 at the "root" ancestor
            but can be changed to create subtrees if desired. 
        """

        nwk = self._POGTreeToNwk(root) + ';'

        with open(file_name, 'w') as f:
            f.write(nwk)

        return nwk
