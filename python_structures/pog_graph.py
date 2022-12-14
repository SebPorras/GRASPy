###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: Trying to implement python version of the POG graph.
# To do:
# 1) Create a sym node class for annotating
# 2) Create a POG graph class
###############################################################################

import numpy as np
#import nptyping as npt
import numpy.typing as npt
from typing import Optional


class Edge(object):
    """Creates instance of an edge between two positions in a sequence.
    Currently only works for bidirectional edges.
    """

    def __init__(self, start: int, end: int, edgeType: Optional[str] = None,
                 recip: Optional[bool] = None, backward: Optional[bool] = None,
                 forward: Optional[bool] = None, weight: Optional[float] = None) -> None:
        """Constructs instance of an edge.

        Parameters:
            start(int): position of beginning of edge

            end(int): position of end of edge

            edgeType(str): Currently only supports bidirectional edge

            recip(bool): ASK ABOUT THIS

            backward(bool): Direction of edge

            forward(bool): Direction of edge

            weight(float): Support of the edge
        """

        self._start = start
        self._end = end
        self._edgeType = edgeType
        self._recip = recip
        self._backward = backward
        self._forward = forward
        self._weight = weight

    def __str__(self) -> str:

        return (f"Start: {self._start}\nEnd: {self._end}\nType: {self._edgeType}\nWeight: {self._weight}\nBackward: {self._backward}\nForward: {self._forward}\nRecip: {self._recip}")


class SymNode(object):
    """Only implemented for output from joint reconstruction.
    Stores the most likely character at a given sequence position
    and all of the edges at this position.
    """

    def __init__(self, name: int, value: str, edges: list) -> None:
        """Constructs instance of the object.

        Parameters:
            name(int): index position in sequence

            value(str): Most likely amino acid based on joint reconstruction

            edges(list): Contains all outgoing edges at this position
        """

        self._name = name
        self._value = value
        self._edges = edges

    def __str__(self) -> str:
        return (f"Name: {self._name}\n Value: {self._value}\n# of edges: {len(self._edges)}")


class POGraph(object):
    """Representation of a sequence as a partial order graph (POG).
    Each sequence position is assigned a SymNode with Edges.
    """

    def __init__(self, version: str, indices: npt.NDArray, nodes: list,
                 start: int, end: int, size: int, terminated: bool,
                 directed: bool, name: str, isAncestor: bool) -> None:
        """Constructs instance of POGraph.

        Parameters:
            version(str): Version of GRASP used for inference

            indices(np.array): Contains all sequence position indices

            nodes(np.array): contains SymNode objects for each position

            start(int): Start index of sequence

            end(int): End index of sequence

            size(int): Size of sequence

            terminated(bool): If the sequence has an end??

            directed(bool): If graph has order??

            name(str): Sequence ID

            isAncestor(bool): Identifier for extant or ancestor
        """

        self._version = version
        self._indices = indices
        self._nodes = nodes
        self._start = start
        self._end = end
        self._size = size
        self._terminated = terminated
        self._directed = directed
        self._name = name
        self._isAncestor = isAncestor

    def __str__(self) -> str:
        return (f"Sequence ID: {self._name}\nSize: {self._size}\nStart: {self._start}\nEnd: {self._end}")


def POGraphFromJSON(jpog: dict, isAncestor: bool = False) -> POGraph:
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
            edge = Edge(start=indices[i], end=-999)
            es.append(edge)
        else:
            # there can be multiple adjacent edges for each sequence position
            for j in range(len(adj[i])):

                edge = Edge(start=indices[i], end=adj[i][j])
                es.append(edge)

        node = SymNode(name=indices[i],
                       value=node_vals[i]["Value"], edges=es)

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

                edge = Edge(start=edgeInd[i][0], end=edgeInd[i][1],
                            edgeType=jpog["Edgetype"],
                            recip=edge_info["Recip"],
                            backward=edge_info["Backward"],
                            forward=edge_info["Forward"],
                            weight=edge_info["Weight"])

                cor_node._edges.append(edge)

            else:

                # edges are stored in random order, need to find correct start

                # returns a tuple with an array containg correct index
                node_loc = np.where(indices == edgeInd[i][0])[0][0]

                cor_node = nodes[node_loc]

                # some edges are identical to adj edges and can be replaced
                dupe_edges = []

                for j in range(len(cor_node._edges)):

                    adjacent_edge = cor_node._edges[j]

                    if edgeInd[i][0] == adjacent_edge._start and\
                            edgeInd[i][1] == adjacent_edge._end:

                        dupe_edges.append(j)

                # remove any identified duplicates
                [cor_node._edges.pop(j) for j in dupe_edges]

                edge = Edge(start=edgeInd[i][0],
                            end=edgeInd[i][1],
                            edgeType=jpog["Edgetype"],
                            recip=edge_info["Recip"],
                            backward=edge_info["Backward"],
                            forward=edge_info["Forward"],
                            weight=edge_info["Weight"])

                cor_node._edges.append(edge)

    # adds ancestor identifier "N"
    name = jpog["Name"]

    if name.isdigit():
        name = "N" + name

    return POGraph(version=jpog["GRASP_version"], indices=indices,
                   start=jpog["Starts"], end=jpog["Ends"], size=jpog["Size"],
                   terminated=jpog["Terminated"], directed=jpog["Directed"],
                   name=name, isAncestor=isAncestor, nodes=nodes)
