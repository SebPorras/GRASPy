###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: This is an implementation of the IDx tree. Can be initialised from a
# JSON object and created.Currently based on output from a joint reconstruction
# Has very basic structures that allow most variables to be accessed or set.
# To do:
# 1) Make a method that can parse this information into a nwk file so that
# it could be read by other tree viewers.
# 2) Investigate about adding in Branchpoint objects that store particular
# info about what parents are and how far away they are etc.
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

            branchpoints(Array): Contains BranchPoint objects

            parents(Array): maps the index of the child to the index
            of the parent

            children(Array): maps the index of the parent to an array
            containing the indexes the child/children

            indices(dict): Maps the sequence ID to the index on the tree

            distances(Array): Maps the branchpoint to the distance to its parent
        """

        self.nBranches = nBranches
        self.branchpoints = branchpoints
        self.parents = parents
        self.children = children
        self.indices = indices
        self.distances = distances

    def __str__(self) -> str:
        return f"Number of branchpoints: {self.nBranches}\nParents: {self.parents}\nChildren: {self.children}\nIndices: {self.indices}\nDistances: {self.distances}"


def IdxTreeFromJSON(serial: dict) -> IdxTree:
    """Instantiates a IdxTree object from a Json file

    Parameters:
        json_file (os.PathLike): path to json file

    Returns:
        IdxTree
    """

    # Currently, assuming that output will contain all these fields

    if ("Distances" or "Branchpoints" or "Labels" or "Parents") \
            not in serial["Input"]["Tree"].keys():
        raise RuntimeError("JSON in incorrect format")

    jdists = serial["Input"]["Tree"]["Distances"]
    nBranches = serial["Input"]["Tree"]["Branchpoints"]
    jlabels = serial["Input"]["Tree"]["Labels"]
    jparents = serial["Input"]["Tree"]["Parents"]

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

    # branch points mirror information in the tree but in human readble form
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

        bp = BranchPoint(id=name,
                         parent=p,
                         dist=distances[BIdx],
                         children=c_lab)

        bpoints[name] = bp

    return IdxTree(nBranches=nBranches, branchpoints=bpoints, parents=parents,
                   children=children, indices=indices, distances=distances)
