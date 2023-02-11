###############################################################################
# Date: 13/12/22
# Author: Sebastian Porras
# Aims: Contains all the data structures for a POG graph.
# A POGraph contains SymNodes that represent characters within a sequence.
# Within these nodes there are Edges which contain information about
# how these SymNodes are connected to one another.
###############################################################################

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

    def __init__(self, name: int, symbol: str, edges: list[Edge]) -> None:
        """Constructs instance of the object.

        Parameters:
            name(int): index position in sequence

            symbol(str): Most likely amino acid based on joint reconstruction

            edges(list): Contains all outgoing edges at this position
        """

        self.name = name
        self.symbol = symbol
        self.edges = edges

    def __str__(self) -> str:
        return (f"Name: {self.name}\nsymbol: {self.symbol}\n# of edges: {len(self.edges)}")


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
