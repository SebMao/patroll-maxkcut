from test_tsp_obstacle import VisibilityGraph, compute_all_pairs_paths
import math
import networkx as nx

class DistanceCalculator:
    def __init__(self, method="euclidean", polygons=None):
        self.method = method
        if method == "visibility" and polygons is not None:
            self.vg = VisibilityGraph(polygons)
            self.path_cache, self.graph = {}, None
            self.dist_matrix = None
        elif method == "euclidean":
            pass
        else:
            raise ValueError("Unsupported distance method or missing polygon data.")

    def compute_all_pairs(self, points):
        if self.method == "euclidean":
            G = nx.Graph()
            for i, p1 in enumerate(points):
                for j, p2 in enumerate(points):
                    if i < j:
                        d = math.dist(p1, p2)
                        G.add_edge(tuple(p1), tuple(p2), weight=d)
            return {}, G
        elif self.method == "visibility":
            self.path_cache, self.graph, self.dist_matrix = compute_all_pairs_paths(self.vg, points)
            return self.path_cache, self.graph, self.dist_matrix

    def get_distance(self, p1, p2):
        if self.method == "euclidean":
            return math.dist(p1, p2)
        elif self.method == "visibility":
            path = self.path_cache.get((tuple(p1), tuple(p2)))
            return sum(math.dist(path[i], path[i + 1]) for i in range(len(path) - 1)) if path else float("inf")

    def get_path(self, p1, p2):
        if self.method == "euclidean":
            return [p1, p2]
        elif self.method == "visibility":
            return self.path_cache.get((tuple(p1), tuple(p2)), [])
