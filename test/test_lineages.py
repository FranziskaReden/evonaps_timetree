from Bio import Phylo
import unittest
import sys
import os
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))
from utils import *

class TestLineageRetrieval(unittest.TestCase):

    def test_replace_tax_id(self):

        data = TimeTree('./')
        data.read_lineages('lineages1.json', 'reverse_lineage.json')

        tax_id = 'C1'
        children = data.lineages_dict[tax_id]['neighbours']
        ancestor = data.reversed_dict[tax_id]

        self.assertEqual(children, ['C'])
        self.assertEqual(data.reversed_dict['C'], tax_id)
        self.assertEqual(ancestor, "ABC")

        data.replace_neighbour('C1', children[0])

        self.assertEqual(data.lineages_dict[ancestor]['neighbours'], ['C', 'AB2'])
    
    def test_combine_neighbours(self):

        data = TimeTree('./')
        data.read_lineages('lineages1.json', 'reverse_lineage.json')

        id = 'AB1'
        ancestor = 'ABC'

        data.combine_neighbours(id, ancestor)
        self.assertEqual(data.lineages_dict['ABC']['neighbours'], ['C1', 'AB2', 'A', 'B'])

    def test_remove_taxon(self):

        data = TimeTree('./')
        data.read_lineages('lineages1.json', 'reverse_lineage.json')

        id = 'AB2'
        ancestor = 'ABC'

        data.remove_from_neighbours(id, ancestor)
        self.assertEqual(data.lineages_dict['ABC']['neighbours'], ['C1'])

    def test_check_ancestry(self):

        data = TimeTree('./')
        data.read_lineages('lineages1.json', 'reverse_lineage.json')


        id1 = 'AB1'
        id2 = 'C1'

        mrca, lin1, lin2 = data.check_ancestry(id1, id2)

        self.assertEqual(mrca, 'ABC')
        self.assertEqual(lin1, ['AB1', 'AB2'])
        self.assertEqual(lin2, ['C1'])

    def test_clean_up_lineages(self):

        data = TimeTree('./')
        data.read_lineages('lineages1.json', 'reverse_lineage.json')

        id1 = 'AB1'
        id2 = 'C1'

        mrca, lin1, lin2 = data.check_ancestry(id1, id2)

        data.clean_up_lineages(lin1, mrca)
        self.assertEqual(data.lineages_dict['ABC']['neighbours'], ['C1', 'A', 'B'])

        data.clean_up_lineages(lin2, mrca)
        self.assertEqual(data.lineages_dict['ABC']['neighbours'], ['A', 'B', 'C'])

    def test_get_leaf_age(self):

        data = TimeTree('./')
        data.tree, data.leaves = read_timetree('tree1.nwk')
        data.read_lineages('lineages1.json', 'reverse_lineage.json')

        pot_leaf = data.get_leaf('ABC')

        self.assertEqual('C', pot_leaf)

        age = data.get_distance("*ABC", 'C')
        self.assertEqual(age, 2)

        age = data.get_distance("*ABC", 'A')
        self.assertEqual(age, 2)

        age = data.get_distance("*AB", 'A')
        self.assertEqual(age, 1)