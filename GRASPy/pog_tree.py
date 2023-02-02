###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: The POG tree will have two components, 1) the Idx tree which encodes
# the topology of the tree. 2) The POG graph that will be stored according
# to the corresponding index on the IDX tree so that both structures
# are in agreement with each other.
#################################################################################


from . import pog_graph
from . import sequence
from typing import Union, Optional


class BranchPoint(object):
    """Represents a branchpoint on a phylogenetic tree. Can contain
    information about the parents or children of that point and how
    long that branch point is.
    """

    def __init__(self, id: str, parent: Union[str, None], dist: float,
                 children: list[str], seq: Optional[sequence.Sequence] = None) -> None:
        """ Constructs instance of a branchpoint.

        Parameters:
            id(str): Sequence ID

            parent(str or None]): ID of parent

            dist(float): Distance to parent

            children(list): IDs children of current BranchPoint

            seq(Sequence): the sequence based on a joint reconstruction
            if the BranchPoint is an ancestor otherwise it is just
            the sequence of an extant. Contains a Sequence object. 
        """

        self.id = id
        self.parent = parent
        self.dist = dist
        self.children = children
        self.seq = seq

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
                 POGraphs: dict[str, pog_graph.POGraph]) -> None:
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

        # Annotate branchpoints with Sequences
        for key, value in self.graphs.items():

            self.branchpoints[key].seq = \
                sequence.Sequence(
                    ''.join([s.symbol for s in value.nodes]), name=key)

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

    def writeToFasta(self, file_name: str) -> None:
        """Writes all sequences of the tree to file.
        Sequence for ancestors are based on a joint 
        reconstruction and each symbol is the most likely 
        at each position. 

        Parameters:

            file_name(str): name of fasta file
        """

        seqs = [bp.seq for bp in self.branchpoints.values()
                if bp.seq is not None]

        sequence.writeFastaFile(file_name, seqs)
