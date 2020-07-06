from unittest import TestCase

import os.path
from embiggen import CSFGraph
from embiggen.hetnode2vec import N2vGraph
from tests.utils.utils import calculate_total_probs

from embiggen.utils import serialize, deserialize


class TestHetGraph(TestCase):

    def setUp(self):
        data_dir = os.path.join(os.path.dirname(__file__), 'data')

        edge_file = os.path.join(data_dir, 'small_graph_all_genes_edges.tsv')
        node_file = os.path.join(data_dir, 'small_graph_all_genes_nodes.tsv')

        g = CSFGraph(edge_file=edge_file, node_file=node_file)
        self.graph = g
        self.nodes = g.nodes()

    def testGraphNodeCounts(self):
        """
        #We expect 6 gene nodes
        :return:
        """
        g = self.graph
        n = g.node_count()
        self.assertEqual(6, n)
