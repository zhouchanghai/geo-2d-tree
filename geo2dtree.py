import math, operator, heapq

__all__ = ('distance', 'Geo2dTree')

def _distance(lat1, lon1, lat2, lon2):
    R = 6371000
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    deltaPhi = phi2 - phi1
    deltaLambda = math.radians(lon2 - lon1)

    sinPhi = math.sin(deltaPhi * 0.5)
    sinLambda = math.sin(deltaLambda * 0.5)
    a = sinPhi * sinPhi + math.cos(phi1) * math.cos(phi2) * sinLambda * sinLambda
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c


def distance(*args):
    '''Distance in meters calculated by Haversine formula.
    See https://www.movable-type.co.uk/scripts/latlong.html
    
    The arguments can be lat1, lon1, lat2, lon2 
    or loc1, loc2 with lat and lon attributes.
    '''
    if len(args) == 4:
        lat1, lon1, lat2, lon2 = args
    elif len(args) == 2:
        lat1, lon1 = args[0].lat, args[0].lon 
        lat2, lon2 = args[1].lat, args[1].lon 
    else:
        raise TypeError
    return _distance(lat1, lon1, lat2, lon2)


class _Node:
    '''A node in Geo2dTree. Internal use only.'''
    
    def __init__(self, location):
        self.loc = location
        self.left = None
        self.right = None

    def coordinate(self, axis: int):
        # Longitude first
        if axis == 0:
            return self.loc.lon
        elif axis == 1:
            return self.loc.lat
        else:
            raise IndexError

    @property
    def lat(self):
        return self.loc.lat
        
    @property
    def lon(self):
        return self.loc.lon

class Geo2dTree:
    '''A 2-dimensional k-d tree divided by longitude and latitude.
    The input locations must have lat and lon attributes.
    
    See https://en.wikipedia.org/wiki/K-d_tree
    The code is adapted from https://github.com/gisalgs/indexing
    '''
    
    def __init__(self, locations):
        '''Build the tree from a list of locations where 
        each location has lat and lon attributes.
        '''
        self.nodes = [_Node(loc) for loc in locations]
        sortedByLat = self.nodes
        sortedByLon = list(self.nodes)
        sortedByLat.sort(key = operator.attrgetter('lat'))
        sortedByLon.sort(key = operator.attrgetter('lon'))
        # First split is along meridian to handle 180 degree meridian wrap.
        self.root = self._build_tree(sortedByLon, sortedByLat)

    def __len__(self):
        return len(self.nodes)

    def nearest_neighbors(self, query, max_dist = math.inf, max_size = 1):
        '''Returns an unsorted list of (distance, location) of the nearest
        neighbors in the tree for the given query. The query argument must 
        have lat and lon attributes. The distances are in meters.
        '''
        results = []
        query = _Node(query)
        self._find_nearest(self.root, query, 0, results, max_dist, max_size)
        return [(-d, loc) for d, loc in results]
        
    def nearest_neighbor(self, query):
        '''Returns (distance, location) of the nearest neighbor in the tree 
        for the given query. The query argument must have lat and lon attributes. 
        The distances are in meters.
        '''
        results = self.nearest_neighbors(query, max_size=1)
        return results[0]

    def _build_tree(self, sortedByMainAxis, sortedBy2ndAxis):
        if not sortedByMainAxis:
            return None

        n = len(sortedByMainAxis)
        midIndex = n // 2
        midNode = sortedByMainAxis[midIndex]
        leftByMain = sortedByMainAxis[ : midIndex]
        rightByMain = sortedByMainAxis[midIndex + 1 : ]

        leftSet = set(leftByMain)
        leftBy2nd = []
        rightBy2nd = []
        for node in sortedBy2ndAxis:
            if node in leftSet:
                leftBy2nd.append(node)
            elif node != midNode:
                rightBy2nd.append(node)

        # Build sub tree on the other axis
        midNode.left = self._build_tree(leftBy2nd, leftByMain)
        midNode.right = self._build_tree(rightBy2nd, rightByMain)
        return midNode

    def _find_nearest(self, root, query, depth, results, max_dist, max_size):
        if not root:
            return
            
        if not root.left and not root.right:
            self._update_nearest(query, root, results, max_dist, max_size)
            return
            
        nearerTree = None
        fartherTree = None
        axis = depth % 2
        if query.coordinate(axis) < root.coordinate(axis):
            nearerTree = root.left
            fartherTree = root.right
        else:
            nearerTree = root.right
            fartherTree = root.left
            
        self._find_nearest(nearerTree, query, depth + 1, results, max_dist, max_size)
        self._update_nearest(query, root, results, max_dist, max_size)
        
        if axis == 0:
            # Splitted by meridian. Find the distance to
            # the other half plane along the parallel.
            distToFartherTree = distance(query.lat, query.lon, query.lat, root.lon)
            if depth == 0:
                # Special case to handle meridian at 180E & 180W
                distTo180Meridian = distance(query.lat, query.lon, query.lat, 180)
                distToFartherTree = min(distToFartherTree, distTo180Meridian)
        else:
            # Splitted by parallel. Find the distance to
            # the other half plane along the meridian.
            distToFartherTree = distance(query.lat, query.lon, root.lat, query.lon)
        
        if (distToFartherTree < max_dist and 
            (len(results) < max_size or distToFartherTree < -results[0][0])):
            self._find_nearest(fartherTree, query, depth + 1, results, max_dist, max_size)

    def _update_nearest(self, query, candidate, results, max_dist, max_size):
        d = distance(query, candidate)
        if d > max_dist:
            return

        pair = (-d, candidate.loc)
        if max_size == 1:
            # Optimization for base case
            if not results:
                results.append(pair)
            elif d < -results[0][0]:
                results[0] = pair
        elif len(results) < max_size:
            heapq.heappush(results, pair)
        elif d < -results[0][0]:
            heapq.heappushpop(results, pair)


