import unittest
from collections import namedtuple

from geo2dtree import Geo2dTree

GeoLocation = namedtuple('GeoLocation', 'lat lon')

class TestGeo2dTree(unittest.TestCase):

    # Example from https://en.wikipedia.org/wiki/K-d_tree
    # where X is longitude and Y is latitude.
    def make_points(self):
        coordinates = [(3,2), (4,5), (6,9), (7,4), (1,8), (2,7)]
        return [GeoLocation(lat, lon) for lat, lon in coordinates]
    
    def assertLatLon(self, loc, lat, lon):
        self.assertAlmostEqual(lat, loc.lat)
        self.assertAlmostEqual(lon, loc.lon)

    def test_build_tree(self):
        tree = Geo2dTree(self.make_points())
        self.assertEqual(6, len(tree))
        root = tree.root
        self.assertLatLon(root, 2, 7)

        node54 = root.left
        self.assertLatLon(node54, 4, 5)

        node96 = root.right
        self.assertLatLon(node96, 6, 9)
        self.assertIsNone(node96.right)

        node81 = node96.left
        self.assertLatLon(node81, 1, 8)
        self.assertIsNone(node81.left)
        self.assertIsNone(node81.right)

        node23 = node54.left
        self.assertLatLon(node23, 3, 2)
        self.assertIsNone(node23.left)
        self.assertIsNone(node23.right)

        node47 = node54.right
        self.assertLatLon(node47, 7, 4)
        self.assertIsNone(node47.left)
        self.assertIsNone(node47.right)

    def test_duplicate(self):
        points = self.make_points()
        # Add a duplicate point
        points.append(GeoLocation(1, 8))

        tree = Geo2dTree(points)
        root = tree.root
        self.assertLatLon(root, 2, 7)

        node81 = root.right;
        self.assertLatLon(node81, 1, 8)

        node81dup = node81.left;
        self.assertLatLon(node81dup, 1, 8)

        node96 = node81.right;
        self.assertLatLon(node96, 6, 9)

    def test_nearest_neighbor(self):
        tree = Geo2dTree(self.make_points())

        # The query is exactly in the tree
        d, best = tree.nearest_neighbor(GeoLocation(3, 2))
        self.assertLatLon(best, 3, 2)
        self.assertEqual(d, 0)

        # The best is on the same side
        d, best = tree.nearest_neighbor(GeoLocation(3.1, 2.1))
        self.assertLatLon(best, 3, 2)
        self.assertAlmostEqual(d, 15714, 0)

        # The best is on the other side
        d, best = tree.nearest_neighbor(GeoLocation(0, 6.9))
        self.assertLatLon(best, 1, 8)

        # The best is the root
        d, best = tree.nearest_neighbor(GeoLocation(3, 6.9))
        self.assertLatLon(best, 2, 7)
        
    def test_wrap_at_180_meridian(self):
        locations = [GeoLocation(0, 0), GeoLocation(0, 178), GeoLocation(0,-179)]
        tree = Geo2dTree(locations)
        root = tree.root

        self.assertLatLon(root, 0, 0);
        self.assertLatLon(root.left, 0, -179)
        self.assertLatLon(root.right, 0, 178)

        # The best is in the same side.
        d, best = tree.nearest_neighbor(GeoLocation(0, 179.4))
        self.assertLatLon(best, 0, 178)

        # The best is in the other side across the 180 degree meridian.
        d, best = tree.nearest_neighbor(GeoLocation(0, 179.6))
        self.assertLatLon(best, 0, -179)
        
    def test_nearest_neighbors(self):
        tree = Geo2dTree(self.make_points())

        neighbors = tree.nearest_neighbors(GeoLocation(2, 6), max_size=3)
        neighbors.sort()
        self.assertAlmostEqual(neighbors[0][0], 111127, 0) #111Km
        self.assertLatLon(neighbors[0][1], 2, 7)
        self.assertAlmostEqual(neighbors[1][0], 248569, 0) #248Km
        self.assertLatLon(neighbors[1][1], 1, 8) 
        self.assertAlmostEqual(neighbors[2][0], 248569, 0) #248Km
        self.assertLatLon(neighbors[2][1], 4, 5) 
        
        neighbors = tree.nearest_neighbors(GeoLocation(2, 6), max_size=3, max_dist=250000)
        neighbors.sort()
        self.assertEqual(3, len(neighbors))
        self.assertLatLon(neighbors[0][1], 2, 7) #111Km
        self.assertLatLon(neighbors[1][1], 1, 8) #248Km
        self.assertLatLon(neighbors[2][1], 4, 5) #248Km

        neighbors = tree.nearest_neighbors(GeoLocation(2, 6), max_size=3, max_dist=150000)
        self.assertEqual(1, len(neighbors))
        self.assertLatLon(neighbors[0][1], 2, 7) #111Km
        
        neighbors = tree.nearest_neighbors(GeoLocation(2, 6), max_size=3, max_dist=100000)
        self.assertEqual(0, len(neighbors))
        
if __name__ == '__main__':
    unittest.main()