from embiggen.graph import Graph, GraphFactory
from unittest import TestCase
import pytest
import numpy as np


class TestGraph(TestCase):

    def setUp(self):
        self._paths = [
            'tests/data/unweighted_small_graph.txt',
            'tests/data/small_het_graph_edges.tsv',
            'tests/data/small_graph.txt',
        ]
        self._legacy_space_paths = [
            "tests/data/karate.train",
            "tests/data/karate.test",
        ]
        self._legacy_paths = [
            "tests/data/small_graph_LEGACY.txt",
            "tests/data/small_g2d_test.txt",
            "tests/data/rand_100nodes_5000edges.graph",
            "tests/data/ppt_train.txt",
            "tests/data/ppt_test.txt",
            *[
                f"tests/data/ppismall_with_validation/{filename}"
                for filename in (
                    "neg_test_edges_max_comp_graph",
                    "neg_train_edges_max_comp_graph",
                    "neg_validation_edges_max_comp_graph",
                    "pos_test_edges_max_comp_graph",
                    "pos_train_edges_max_comp_graph",
                    "pos_validation_edges_max_comp_graph"
                )
            ],
            *[
                f"tests/data/ppismall/{filename}"
                for filename in (
                    "neg_test_edges",
                    "neg_train_edges",
                    "pos_test_edges",
                    "pos_train_edges",
                )
            ],
            *[
                f"tests/data/karate/{filename}"
                for filename in (
                    "neg_test_edges",
                    "neg_train_edges",
                    "neg_validation_edges",
                    "pos_test_edges",
                    "pos_train_edges",
                    "pos_validation_edges",
                )
            ]
        ]

        self._factory = GraphFactory(Graph)

    def test_setup_from_dataframe(self):
        graph = self._factory.read_csv('tests/data/unweighted_small_graph.txt')
        assert len(graph.nodes_indices) > 0

        # Testing illegal weights arguments
        with pytest.raises(ValueError):
            graph = self._factory.read_csv(
                'tests/data/unweighted_small_graph.txt',
                weights=0
            )
        with pytest.raises(ValueError):
            graph = self._factory.read_csv(
                'tests/data/unweighted_small_graph.txt',
                weights=[0]
            )

        # Testing illegal edge types parameter
        with pytest.raises(ValueError):
            graph = self._factory.read_csv(
                'tests/data/unweighted_small_graph.txt',
                edge_types=[0]
            )

        # Testing illegal node types parameter
        with pytest.raises(ValueError):
            graph = self._factory.read_csv(
                'tests/data/unweighted_small_graph.txt',
                node_types=[0]
            )

    def test_normalize_graph(self):
        """Testing that the normalization process actually works."""
        for path in self._paths:
            graph = self._factory.read_csv(path, normalize_weights=True)
            assert all(
                np.isclose(neighbours_weights.sum(), 1)
                for neighbours_weights in graph._neighbours_weights
            )
            assert all(
                graph.has_edge(edge)
                for edge in graph.edges_indices
            )
            assert all(
                not graph.has_edge((-src-1, -dst-1))
                for src, dst in graph.edges_indices
            )

    def test_not_normalize_graph(self):
        """Testing that the normalization process actually works."""
        for path in self._paths:
            graph = self._factory.read_csv(path, normalize_weights=False)
            assert not all(
                np.isclose(neighbours_weights.sum(), 1)
                for neighbours_weights in graph._neighbours_weights
            )
            assert all(
                graph.has_edge(edge)
                for edge in graph.edges_indices
            )
            assert all(
                not graph.has_edge((-src-1, -dst-1))
                for src, dst in graph.edges_indices
            )

    def test_normalize_graph_legacy(self):
        """Testing that the normalization process actually works."""
        for path in self._legacy_paths:
            graph = self._factory.read_csv(
                path,
                normalize_weights=True,
                edge_has_header=False,
                start_nodes_column=0,
                end_nodes_column=1,
                weights_column=2
            )
            assert len(graph.edges_indices) == graph.edges_number
            assert len(graph.nodes_indices) == graph.nodes_number
            assert all(
                np.isclose(neighbours_weights.sum(), 1)
                for neighbours_weights in graph._neighbours_weights
            )
            assert all(
                graph.has_edge(edge)
                for edge in graph.edges_indices
            )
            assert all(
                not graph.has_edge((-src-1, -dst-1))
                for src, dst in graph.edges_indices
            )

        for path in self._legacy_space_paths:
            graph = self._factory.read_csv(
                path,
                edge_sep=" ",
                normalize_weights=True,
                edge_has_header=False,
                start_nodes_column=0,
                end_nodes_column=1,
                weights_column=2
            )
            assert all(
                np.isclose(neighbours_weights.sum(), 1)
                for neighbours_weights in graph._neighbours_weights
            )
            assert all(
                graph.has_edge(edge)
                for edge in graph.edges_indices
            )
            assert all(
                not graph.has_edge((-src-1, -dst-1))
                for src, dst in graph.edges_indices
            )

    def test_setup_from_custom_dataframe(self):
        # TODO: integrate all other remaining columns
        graph = self._factory.read_csv(
            "tests/data/small_9606.protein.actions.txt",
            normalize_weights=False,
            start_nodes_column="item_id_a",
            end_nodes_column="item_id_b",
            weights_column="score"
        )
        assert not all(
            np.isclose(neighbours_weights.sum(), 1)
            for neighbours_weights in graph._neighbours_weights
        )
        graph = self._factory.read_csv(
            "tests/data/small_9606.protein.actions.txt",
            normalize_weights=True,
            start_nodes_column="item_id_a",
            end_nodes_column="item_id_b",
            weights_column="score"
        )
        assert all(
            np.isclose(neighbours_weights.sum(), 1)
            for neighbours_weights in graph._neighbours_weights
        )