"""Some tests to make sure the filters work as intended"""
import unittest
from pkg_resources import resource_filename

from locals import SourceCatalog

 
class TestSourceCatalog(unittest.TestCase):
    """Tests for Filter class"""
    def setUp(self):
        self.cat_path = resource_filename('locals', 'data/fake/')

    def test_Catalog(self):
        """Test if a SourceCatalog object is created properly"""
        cat = SourceCatalog(self.cat_path, color_cut='brown_dwarfs')

        self.assertEqual(len(cat.results), 5)