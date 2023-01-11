###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: Trying to implement python version of the POG graph.
# To do:
# 1) Create a sym node class for annotating
# 2) Create a POG graph class
###############################################################################

import numpy as np
from numpy.typing import NDArray
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

        self.start = start
        self.end = end
        self.edgeType = edgeType
        self.recip = recip
        self.backward = backward
        self.forward = forward
        self.weight = weight

    def __str__(self) -> str:

        return (f"Start: {self.start}\nEnd: {self.end}\nType: {self.edgeType}\nWeight: {self.weight}\nBackward: {self.backward}\nForward: {self.forward}\nRecip: {self.recip}")


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

        self.name = name
        self.value = value
        self.edges = edges

    def __str__(self) -> str:
        return (f"Name: {self.name}\nValue: {self.value}\n# of edges: {len(self.edges)}")


class POGraph(object):
    """Representation of a sequence as a partial order graph (POG).
    Each sequence position is assigned a SymNode with Edges.
    """

    def __init__(self, version: str, indices: NDArray, nodes: list[SymNode],
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

        self.version = version
        self.indices = indices
        self.nodes = nodes
        self.start = start
        self.end = end
        self.size = size
        self.terminated = terminated
        self.directed = directed
        self.name = name
        self.isAncestor = isAncestor

    def __str__(self) -> str:
        return (f"Sequence ID: {self.name}\nSize: {self.size}\nStart: {self.start}\nEnd: {self.end}")


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

                cor_node.edges.append(edge)

            else:

                # edges are stored in random order, need to find correct start

                # returns a tuple with an array containg correct index
                node_loc = np.where(indices == edgeInd[i][0])[0][0]

                cor_node = nodes[node_loc]

                # some edges are identical to adj edges and can be replaced
                dupeedges = []

                for j in range(len(cor_node.edges)):

                    adjacent_edge = cor_node.edges[j]

                    if edgeInd[i][0] == adjacent_edge.start and\
                            edgeInd[i][1] == adjacent_edge.end:

                        dupeedges.append(j)

                # remove any identified duplicates
                [cor_node.edges.pop(j) for j in dupeedges]

                edge = Edge(start=edgeInd[i][0],
                            end=edgeInd[i][1],
                            edgeType=jpog["Edgetype"],
                            recip=edge_info["Recip"],
                            backward=edge_info["Backward"],
                            forward=edge_info["Forward"],
                            weight=edge_info["Weight"])

                cor_node.edges.append(edge)

    # adds ancestor identifier "N"
    name = jpog["Name"]

    if name.isdigit():
        name = "N" + name

    return POGraph(version=jpog["GRASP_version"], indices=indices,
                   start=jpog["Starts"], end=jpog["Ends"], size=jpog["Size"],
                   terminated=jpog["Terminated"], directed=jpog["Directed"],
                   name=name, isAncestor=isAncestor, nodes=nodes)
