from unittest import TestCase

import os.path
from embiggen import CSFGraph
from embiggen.hetnode2vec_edgetype import N2vGraph
from tests.utils.utils import calculate_total_probs


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

    def testGraphEdgeCounts(self):
        """
         # We expect 18 edges
         # Note that currently we have 2 directed edges for each undirected edge. This
         # means that edge_count() returns . This is an implementation detail that
         # may change in the future.
        :return:
        """
        g = self.graph
        m = g.edge_count()
        self.assertEqual(18, m)

    def test_raw_probs_1(self):
        p = 1
        q = 1
        gamma = 1
        g = N2vGraph(self.graph, p, q, gamma, doxn2v=True)
        src = 'g2'
        dst = 'g6'
        [j_alias, q_alias] = g.get_alias_edge_xn2v_edgetype(src, dst)
        self.assertEqual(len(j_alias), len(q_alias))
        # outgoing edges from g6: g2
        self.assertEqual(1, len(j_alias))
        # recreate the original probabilities.
        original_probs = calculate_total_probs(j_alias, q_alias)
        self.assertAlmostEqual(1.0, original_probs[0])


    def test_raw_probs_2(self):
        p = 1
        q = 1
        gamma = 1
        g = N2vGraph(self.graph, p, q, gamma, doxn2v=True)
        src = 'g2'
        dst = 'g4'
        [j_alias, q_alias] = g.get_alias_edge_xn2v_edgetype(src, dst)
        self.assertEqual(len(j_alias), len(q_alias))
        # outgoing edges from g4: g2, g1 both of type 2. Weight g2 <->g4 =9, weight g4 <->g1 =10
        #0.5 * 10 = 5, 0.5 * 9 = 4.5
        #5/9.5 = 0.526, 4.5/9.5 = 0.473
        self.assertEqual(2, len(j_alias))
        # recreate the original probabilities.
        original_probs = calculate_total_probs(j_alias, q_alias)
        self.assertAlmostEqual(5/9.5, original_probs[0])
        self.assertAlmostEqual(4.5/9.5, original_probs[1])

    def test_raw_probs_3(self):
        p = 1
        q = 1
        gamma = 1
        g = N2vGraph(self.graph, p, q, gamma, doxn2v=True)
        src = 'g1'
        dst = 'g3'
        [j_alias, q_alias] = g.get_alias_edge_xn2v_edgetype(src, dst)
        self.assertEqual(len(j_alias), len(q_alias))
        # outgoing edges from g3: g2, g1, g5. g3 <-> g1 of type 1 and weight 10, g3 <->g2 of type 1 and weight 14,
        # g3 <-> g5 of type 3 and weight 10

        self.assertEqual(3, len(j_alias))
        # recreate the original probabilities.
        original_probs = calculate_total_probs(j_alias, q_alias)
        self.assertAlmostEqual(2.5/12, original_probs[0])
        self.assertAlmostEqual(3.5/12, original_probs[1])
        self.assertAlmostEqual(6/12, original_probs[2])