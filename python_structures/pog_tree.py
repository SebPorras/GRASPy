###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: The POG tree will have two components, 1) the Idx tree which encodes
# the topology of the tree. 2) The POG graph that will be stored according
# to the corresponding index on the IDX tree so that both structures
# are in agreement with each other.
#################################################################################


from pog_graph import *
from typing import Union


class BranchPoint(object):
    """Represents a branchpoint on a phylogenetic tree. Can contain
    information about the parents or children of that point and how
    long that branch point is.

    Future use -> can associate annotations with this branch point
    but must be in JSON format
    """

    def __init__(self, id: str, parent: Union[str, None], dist: float,
                 children: list[str]) -> None:
        """ Constructs instance of a branchpoint.

        Parameters:
            id(str): Sequence ID

            parent(int): Index of parent branchpoint

            dist(float): Distance to parent

            children(np.array): all children of branchpoint
        """

        self.id = id
        self.parent = parent
        self.dist = dist
        self.children = children
        self.seq = ""

    def __str__(self) -> str:

        return (f"Name: {self.id}\nParent: {self.parent}\nDistance To Parent {self.dist}\nChildren: {self.children}")


class POGTree(object):
    """The Partial Order Graph Tree (POGTree), is a phylogenetic tree made up 
    of branchpoints which represent nodes on the tree.
    Each branchpoint is assigned an index and a BranchPoint
    object, allowing easy access of information via the sequence name of
    an extant or an ancestor. 2) A POGraph object which describes
    the graph of the sequence at that branchpoint.
    """

    def __init__(self, nBranches: int, branchpoints: dict[str, BranchPoint],
                 parents: list[int], children: list[list[Union[int, None]]],
                 indices: dict[str, int], distances: list[float],
                 POGraphs: dict[str, POGraph]) -> None:
        """Constructs instance of POGTree

        Parameters:
            nBranchs(int): number of branch points in the tree

            branchpoints(dict[str, BranchPoint]): Contains BranchPoint objects

            parents(list[int]): maps the index of the child to the index
            of the parent

            children(list[list[Union[int, None]]]): maps the index of the
            parent to an array containing the indexes the child/children

            indices(dict[str, int]): Maps the sequence ID to the index on the
            tree

            distances(list[float]): Maps the branchpoint to the distance to
            its parent

            POGraph(POGraph): Instance of a POGraph object
        """

        self.nBranches = nBranches
        self.branchpoints = branchpoints
        self.parents = parents
        self.children = children
        self.indices = indices
        self.distances = distances
        self.graphs = POGraphs

        # Annotate branchpoints with any sequences
        for key, value in self.graphs.items():

            self.branchpoints[key].seq = \
                ''.join([s.symbol for s in value.nodes])

    def __str__(self) -> str:
        return f"Number of branchpoints: {self.nBranches}\nParents: {self.parents}\nChildren: {self.children}\nIndices: {self.indices}\nDistances: {self.distances}"

    def _parseToNwk(self, root: str = "N0",) -> str:
        """Converts the POGTree to nwk form Newick Standard (nwk) format.
        Root can be changed if the user wishes to create subtrees.

        Parameters:
            root(str): Default set to N0 at the "root" ancestor
            but can be changed to create subtrees if desired. 

        Returns:
            str: The POGTree in nwk format 
        """

        nwk = ""

        if None in self.branchpoints[root].children:

            nwk += f"{self.branchpoints[root].id}:{self.branchpoints[root].dist}"

        else:

            for c in self.branchpoints[root].children:
                nwk += self._parseToNwk(
                    self.branchpoints[c].id) + ','

            # slicing removes final ',' at end of subtree
            nwk = "(" + nwk[:-1] + \
                f"){self.branchpoints[root].id}:{self.branchpoints[root].dist}"

        return nwk

    def writeToNwk(self, file_name: str, root: str = "N0") -> str:
        """Writes a nwk string of the tree to a file

        Parameters: 

            file_name(str): name of nwk file

            root(str): Default set to N0 at the "root" ancestor
            but can be changed to create subtrees if desired. 
        """

        nwk = self._parseToNwk(root) + ';'

        with open(file_name, 'w') as f:
            f.write(nwk)

        return nwk
