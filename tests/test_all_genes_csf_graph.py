from unittest import TestCase
import os.path
from embiggen import CSFGraph
from embiggen import Edge
from embiggen.csf_graph.csf_graph import CSFGraphNoSubjectColumnError, \
    CSFGraphNoObjectColumnError


class TestCSFGraph(TestCase):
    def setUp(self):
        data_dir = os.path.join(os.path.dirname( __file__), 'data')

        # files for canonical test graph
        self.edge_file = os.path.join(data_dir, 'small_graph_all_genes_edges.txt')
        self.node_file = os.path.join(data_dir, 'small_graph_all_genes_nodes.txt')
        g = CSFGraph(edge_file=self.edge_file)
        self.g = g
        str(g)

    def test_notnull(self):
        self.assertIsNotNone(self.g)

        #
        # check maps
        #

    def test_get_node_to_index_map(self):
        self.assertIsNotNone(self.g.get_node_to_index_map())

    def test_get_index_to_node_map(self):
        self.assertIsNotNone(self.g.get_index_to_node_map())

    def test_edgetype2count_dictionary(self):
        self.assertIsInstance(self.g.edgetype2count_dictionary, dict)
        self.assertEqual(self.g.edgetype2count_dictionary['1'], 3)
        self.assertEqual(self.g.edgetype2count_dictionary['2'], 3)
        self.assertEqual(self.g.edgetype2count_dictionary['3'], 3)

    #def test_nodetype2count_dictionary(self):
        #het_g = CSFGraph(edge_file=self.edge_file, node_file=self.node_file)
        #self.assertIsInstance(self.g.nodetype2count_dictionary, dict)
        #self.assertEqual(het_g.nodetype2count_dictionary['gene'], 9)

    # check nodetype to index map
    def test_csfgraph_makes_nodetype_to_index_map(self):
        self.assertIsInstance(self.g.nodetype_to_index_map, dict)

    def test_csfgraph_assigns_default_nodetype_to_nodetype_to_index_map(self):
        self.assertIsInstance(self.g.nodetype_to_index_map, dict)
        self.assertEqual(self.g.nodetype_to_index_map[self.g.default_node_type],
                         list(range(self.g.node_count())))

    def test_csfgraph_populates_edgetype_to_index_map(self):
        self.assertCountEqual(self.g.edgetype_to_index_map.keys(),
                              ['1','2','3'])
        self.assertEqual(6, len(self.g.edgetype_to_index_map['1']))
        self.assertEqual(6, len(self.g.edgetype_to_index_map['2']))

    def test_csfgraph_constructor_makes_index_to_edgetype_map(self):
        self.assertIsInstance(self.g.index_to_edgetype_map, dict)

    def test_csfgraph_populates_index_to_edgetype_map(self):
        self.assertEqual(18, len(self.g.index_to_edgetype_map))
        self.assertEqual(self.g.index_to_edgetype_map[0], '1')
        self.assertEqual(self.g.index_to_edgetype_map[7], '2')
        self.assertEqual(self.g.index_to_edgetype_map[16], '3')

    def test_count_edges(self):
        # Note that the graph is transformed into an undirected graph by
        # adding the inverse of each edge
        # 9 edges are included in the file, thus we expect 2*9=18 edges
        expected_num = 18
        self.assertEqual(expected_num, self.g.edge_count())
        return None

    def test_get_weight(self):
        expected = 10  # for g1<->g2
        weight = self.g.weight('g1', 'g2')
        self.assertEqual(expected, weight)
        expected = 15  # for g3<->g5
        weight = self.g.weight('g3', 'g5')
        self.assertEqual(expected, weight)

        return None

    def test_get_neighbors1(self):
        # the neighbors of g4 are g1 and g2 (sorted neighbors)
        nbrs = self.g.neighbors('g4')
        # print(nbrs)
        self.assertEqual(['g1', 'g2'], nbrs)

        return None

    def test_get_neighbors_as_ints_1(self):
        # The neighbors of g4 are g1, g2,  and their indices are 0,1
        g4_idx = self.g.node_to_index_map['g4']
        nbrs = self.g.neighbors_as_ints(g4_idx)
        self.assertEqual([0,1], nbrs)

    def test_get_neighbors2(self):
        # the neighbors of g2 are g1, g3, g4, g5,g6
        nbrs = self.g.neighbors('g2')
        # print(nbrs)
        self.assertEqual(['g1', 'g3', 'g4', 'g5', 'g6'], nbrs)
        return None

    def test_degree(self):
        #test the degrees of nodes
        node_1 = "g1"
        node_2 = "g2"
        node_3 = "g4"
        self.assertEqual(4, self.g.node_degree(node_1))
        self.assertEqual(5, self.g.node_degree(node_2))
        self.assertEqual(2, self.g.node_degree(node_3))

    def test_edge_type_1(self):
        node_1 = "g1"
        node_2 = "g3"
        edge_type = self.g.edgetype(node_1, node_2)
        expected = "1"
        self.assertEqual(expected,edge_type)

    def test_edge_type_2(self):
        node_1 = "g1"
        node_2 = "g5"
        edge_type = self.g.edgetype(node_1, node_2)
        expected = "3"
        self.assertEqual(expected,edge_type)

    def test_edge_type_3(self):
        node_1 = "g2"
        node_2 = "g6"
        edge_type = self.g.edgetype(node_1, node_2)
        expected = "3"
        self.assertEqual(expected,edge_type)