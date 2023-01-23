###############################################################################
# Date: 20/1/23
# Author: Sebastian Porras
# Aims: IdxTree is condensed phylogenetic tree representation within a POGTree
# object. Each branchpoint is assigned an index and a BranchPoint object,
# allowing easy access of information via the sequence name of an extant
# or an ancestor.
#################################################################################

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
        self.annotations = dict()

    def __str__(self) -> str:

        return (f"Name: {self.id}\nParent: {self.parent}\nDistance To Parent {self.dist}\nChildren: {self.children}")


class IdxTree(object):
    """IdxTree is condensed phylogenetic tree representation within a 
    POGTree object. Each branchpoint is assigned an index and a BranchPoint
    object, allowing easy access of information via the sequence name of
    an extant or an ancestor. 
    """

    def __init__(self, nBranches: int, branchpoints: dict[str, BranchPoint],
                 parents: list[int], children: list[list[Union[int, None]]],
                 indices: dict[str, int], distances: list[float]) -> None:
        """
        Parameters:
            nBranchs(int): number of branch points in the tree

            branchpoints(dict[str, BranchPoint]): Contains BranchPoint objects

            parents(list[int]): maps the index of the child to the index
            of the parent

            children(list[list[Union[int, None]]]): maps the index of the parent to an array
            containing the indexes the child/children

            indices(dict[str, int]): Maps the sequence ID to the index on the tree

            distances(list[float]): Maps the branchpoint to the distance to its parent
        """

        self.nBranches = nBranches
        self.branchpoints = branchpoints
        self.parents = parents
        self.children = children
        self.indices = indices
        self.distances = distances

    def __str__(self) -> str:
        return f"Number of branchpoints: {self.nBranches}\nParents: {self.parents}\nChildren: {self.children}\nIndices: {self.indices}\nDistances: {self.distances}"
