from typing import List, Tuple
import numpy as np  # type: ignore
from numba.experimental import jitclass  # type: ignore
from numba import typed, types  # type: ignore
from .alias_method import alias_draw
from .build_alias import build_alias_edges, build_alias_nodes
from .graph_types import (
    numba_vector_nodes_type,
    numba_edges_type,
    alias_list_type,
    edges_type_list
)


@jitclass([
    ('_destinations', numba_vector_nodes_type),
    ('_sources', numba_vector_nodes_type),
    ('_neighbors', types.ListType(edges_type_list)),
    ('_nodes_alias', types.ListType(alias_list_type)),
    ('_edges_alias', types.ListType(alias_list_type)),
    ('has_traps', types.boolean),
    ('uniform', types.boolean),
    ('random_walk_preprocessing', types.boolean),
])
class NumbaGraph:

    def __init__(
        self,
        nodes_number: int,
        sources: List[int],
        destinations: List[int],
        node_types: List[np.uint16] = None,
        edge_types: List[np.uint16] = None,
        weights: List[float] = None,
        default_weight: float = 1.0,
        return_weight: float = 1.0,
        explore_weight: float = 1.0,
        change_node_type_weight: float = 1.0,
        change_edge_type_weight: float = 1.0,
    ):
        """Crate a new instance of a undirected graph with given edges.

        Parameters
        -------------------------
        nodes_number: int,
            Number of nodes in the graph.
        sources: List[int],
            List of the source nodes of the graph.
        destinations: List[int],
            List of the destination nodes of the graph.
        node_types: List[np.uint16] = None,
            The node types for each node.
            This is an optional parameter to make the graph behave as if it
            is colored within the walk.
        edge_types: List[np.uint16],
            The edge types for each source and sink.
            This is an optional parameter to make the graph behave as if it
            is a multigraph within the walk.
        weights: List[float] = None,
            The weights for each source and sink. By default None. If you want
            to specify a single value to be used for every weight, use the
            default_weight parameter.
        default_weight: int = 1.0,
            The default weight to use when no weight is provided.
        return_weight : float = 1.0,
            Weight on the probability of returning to node coming from
            Having this higher tends the walks to be
            more like a Breadth-First Search.
            Having this very high  (> 2) makes search very local.
            Equal to the inverse of p in the Node2Vec paper.
        explore_weight : float = 1.0,
            Weight on the probability of visiting a neighbor node
            to the one we're coming from in the random walk
            Having this higher tends the walks to be
            more like a Depth-First Search.
            Having this very high makes search more outward.
            Having this very low makes search very local.
            Equal to the inverse of q in the Node2Vec paper.
        change_node_type_weight: float = 1.0,
            Weight on the probability of visiting a neighbor node of a
            different type than the previous node. This only applies to
            colored graphs, otherwise it has no impact.
        change_edge_type_weight: float = 1.0,
            Weight on the probability of visiting a neighbor edge of a
            different type than the previous edge. This only applies to
            multigraphs, otherwise it has no impact.
        
        Raises
        -------------------------
        ValueError,
            If given node types has not the same length of given nodes list.
        ValueError,
            If given sources length does not match destinations length.
        ValueError,
            If return_weight is not a strictly positive real number.
        ValueError,
            If explore_weight is not a strictly positive real number.
        ValueError,      
            If change_node_type_weight is not a strictly positive real number.
        ValueError,
            If change_edge_type_weight is not a strictly positive real number.

        """

        if len(sources) != len(destinations):
            raise ValueError(
                "Given sources length does not match destinations length."
            )
        if edge_types is not None and len(edge_types) != len(destinations):
            raise ValueError(
                "Given edge types length does not match destinations length."
            )
        if weights is not None and len(weights) != len(destinations):
            raise ValueError(
                "Given weights length does not match destinations length."
            )
        if node_types is not None and len(nodes) != len(node_types):
            raise ValueError(
                "Given node types has not the same length of given nodes list."
            )
        if return_weight <= 0:
            raise ValueError("Given return weight is not a positive number")
        if explore_weight <= 0:
            raise ValueError("Given explore weight is not a positive number")
        if change_node_type_weight <= 0:
            raise ValueError(
                "Given change_node_type_weigh is not a positive number"
            )
        if change_edge_type_weight <= 0:
            raise ValueError(
                "Given change_edge_type_weight is not a positive number"
            )
        
        self._destinations = destinations
        self._sources = sources

        # Each node has a list of neighbors.
        # These lists are initialized as empty.
        self._neighbors = typed.List.empty_list(edges_type_list)
        for _ in range(len(nodes_number)):
            self._neighbors.append(
                typed.List.empty_list(numba_edges_type)
            )

        # Allocating the vectors of the mappings
        edges_set = set()
        self._destinations = np.empty(len(edges), dtype=numpy_nodes_type)
        self._sources = np.empty(len(edges), dtype=numpy_nodes_type)

        # The following proceedure ASSUMES that the edges only appear
        # in a single direction. This must be handled in the preprocessing
        # of the graph parsing proceedure.
        for i, (source, destination) in enumerate(edges):
            # Create the sources numeric ID
            src = self._nodes_mapping[str(source)]
            # Create the destinations numeric ID
            dst = self._nodes_mapping[str(destination)]
            # Store the destinations into the destinations vector
            self._destinations[i] = dst
            # Store the sources into the sources vector
            self._sources[i] = src
            # If the preprocessing is required we compute the neighbours
            if self.random_walk_preprocessing:
                # Appending outbound edge ID to SRC list.
                self._neighbors[src].append(i)
                # Storing the edges mapping.
                edges_set.add((src, dst))

        if not self.random_walk_preprocessing:
            return

        # Creating the node alias list, which contains tuples composed of
        # the list of indices of the opposite extraction events and the list
        # of probabilities for the extraction of edges neighbouring the nodes.
        self._nodes_alias = build_alias_nodes(self._neighbors, weights)

        # Creating the edges alias list, which contains tuples composed of
        # the list of indices of the opposite extraction events and the list
        # of probabilities for the extraction of edges neighbouring the edges.
        self._edges_alias = build_alias_edges(
            edges_set, self._neighbors,
            node_types, edge_types, weights,
            self._sources, self._destinations,
            return_weight, explore_weight,
            change_edge_type_weight, change_edge_type_weight
        )

        # To verify if this graph has some walker traps, meaning some nodes
        # that do not have any neighbors, we have to iterate on the list of
        # neighbors and to check if at least a node has no neighbors.
        # If such a condition is met, we cannot anymore do the simple random
        # walk assuming that all the walks have the same length, but we need
        # to create a random walk with variable length, hence a list of lists.

        self.has_traps = False
        for src in range(len(self._nodes_alias)):
            if self.is_node_trap(src):
                self.has_traps = True
                break

    @property
    def nodes_number(self) -> int:
        """Return the total number of nodes in the graph.

        Returns
        -------
        The total number of nodes in the graph.
        """
        return len(self._nodes_alias)

    def is_node_trap(self, node: int) -> bool:
        """Return boolean representing if node is a dead end.

        Parameters
        ----------
        node: int,
            Node numeric ID.

        Returns
        -----------------
        Boolean True if node is a trap.
        """
        return len(self._neighbors[node]) == 0

    def is_edge_trap(self, edge: int) -> bool:
        """Return boolean representing if edge is a dead end.

        Parameters
        ----------
        edge: int,
            Edge numeric ID.

        Returns
        -----------------
        Boolean True if edge is a trap.
        """
        return self.is_node_trap(self._destinations[edge])

    def neighbors(self, node: str) -> List:
        """Return neighbors of given node.

        Parameters
        ---------------------
        node: str,
            The node whose neigbours are to be identified.

        Returns
        ---------------------
        List of neigbours of given node.
        """
        # Get the numeric ID of the node
        node_id = self._nodes_mapping[node]
        # We get the node neighbours
        neighbours = self._neighbors[node_id]
        # And translate them back.
        return [
            # For each neighbor edge, we need to get first the destination ID
            # and then remap the destination ID to the original node name
            self._reverse_nodes_mapping[self._destinations[neighbor]]
            for neighbor in neighbours
        ]

    def degree(self, node: str) -> int:
        """Return degree of given node.

        Parameters
        ---------------------
        node: str,
            The node whose neigbours are to be identified.

        Returns
        ---------------------
        Number of neighbors of given node.
        """
        return len(self.neighbors(node))

    def extract_transition_informations(self, src: int, j: np.ndarray, q: np.ndarray) -> Tuple[int, int]:
        # TODO! Add docstring!
        neighbor_index = alias_draw(j, q)
        neighbours = self._neighbors[src]
        edge = neighbours[neighbor_index]
        # Get the destination of the chosen edge.
        dst = self._destinations[edge]
        # Return the obtained tuple
        return dst, edge

    def extract_random_node_neighbour(self, src: int) -> Tuple[int, int]:
        """Return a random adiacent node to the one associated to node.
        The destination is extracted by using the normalized weights
        of the edges as probability distribution.

        Parameters
        ----------
        src: int
            The index of the source node that is to be considered.

        Returns
        -------
        A tuple containing the index of a random adiacent node to given
        source node and the ID of th edge used for the transition between
        the two nodes.
        """
        # !TODO! UPDATE THIS METHOD!
        # Get the information relative to the source node, composed of a tuple:
        # - The numpy array of opposite events for the alias method (j)
        # - The probabilities for the extractions for the alias method (q)
        j, q = self._nodes_alias[src]

        return self.extract_transition_informations(src, j, q)

    def extract_random_edge_neighbour(self, edge: int) -> Tuple[int, int]:
        """Return a random adiacent edge to the one associated to edge.
        The Random is extracted by using the normalized weights of the edges
        as probability distribution.

        Parameters
        ----------
        edge: int
            The index of the egde that is to be considered.

        Returns
        -------
        The index of a random adiacent edge to edge.
        """
        # !TODO! UPDATE THIS METHOD!
        # We retrieve the destination of edge currently used, which can be
        # considered the source of the edge we are looking for now.
        src = self._destinations[edge]
        # Get the information relative to the source node, composed of a tuple:
        # - The numpy array of opposite events for the alias method (j)
        # - The probabilities for the extractions for the alias method (q)
        j, q = self._edges_alias[edge]

        return self.extract_transition_informations(src, j, q)
