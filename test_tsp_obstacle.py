
import math
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import nearest_points
import heapq
import networkx as nx
import matplotlib.pyplot as plt
import copy


class VisibilityGraph:
    def __init__(self, polygons):
        self.polygons = polygons
        self.edges = self._extract_edges()
        self.graph = self.build_visibility_graph()

    def _extract_edges(self):
        edges = []
        for poly in self.polygons:
            coords = list(poly.exterior.coords)
            for i in range(len(coords) - 1):
                edges.append(LineString([coords[i], coords[i + 1]]))
        return edges

    def is_visible(self, p1, p2):
        segment = LineString([p1, p2])
        for edge in self.edges:
            if segment.crosses(edge):
                return False
            if segment.touches(edge):
                edge_points = [Point(edge.coords[0]), Point(edge.coords[1])]
                if not (Point(p1) in edge_points or Point(p2) in edge_points):
                    return False

        # Also ensure the midpoint lies within free space
        midpoint = Point((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
        for poly in self.polygons:
            if poly.contains(midpoint):
                return False
        return True

    def build_visibility_graph(self):
        points = [tuple(p) for poly in self.polygons for p in poly.exterior.coords[:-1]]
        graph = {p: {} for p in points}
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                if self.is_visible(points[i], points[j]):
                    dist = math.dist(points[i], points[j])
                    graph[points[i]][points[j]] = dist
                    graph[points[j]][points[i]] = dist
        return graph

    def add_point_to_graph(self, point):
        point_t = tuple(point)
        if point_t in self.graph:
            return  # already added
        self.graph[point_t] = {}
        for v in self.graph.keys():
            if v == point_t:
                continue
            if self.is_visible(point, v):
                dist = math.dist(point, v)
                self.graph[point_t][v] = dist
                self.graph[v][point_t] = dist

    def dijkstra(self, start, goal):
        start = tuple(start)
        goal = tuple(goal)

        # Add start and goal into the graph
        added_nodes = []
        if start not in self.graph:
            self.add_point_to_graph(start)
            added_nodes.append(start)
        if goal not in self.graph:
            self.add_point_to_graph(goal)
            added_nodes.append(goal)

        # Dijkstra's algorithm
        queue = [(0, start)]
        distances = {start: 0}
        previous = {}

        while queue:
            cost, u = heapq.heappop(queue)
            if u == goal:
                break
            for v, w in self.graph[u].items():
                alt = cost + w
                if v not in distances or alt < distances[v]:
                    distances[v] = alt
                    previous[v] = u
                    heapq.heappush(queue, (alt, v))

        # Reconstruct path
        path = []
        u = goal
        while u in previous:
            path.append(u)
            u = previous[u]
        path.append(start)
        path.reverse()

        # Clean up added nodes
        for node in added_nodes:
            for neighbor in self.graph[node]:
                if node in self.graph[neighbor]:
                    del self.graph[neighbor][node]
            del self.graph[node]
        return path



def plot_visibility_graph(vg: VisibilityGraph, path=None, start=None, goal=None):
    fig, ax = plt.subplots(figsize=(8, 8))

    # 绘制障碍物
    for poly in vg.polygons:
        x, y = poly.exterior.xy
        ax.fill(x, y, alpha=0.4, fc='gray', ec='black', linewidth=1.5, label="Obstacle")

    # 绘制可视图中的边
    for u in vg.graph:
        for v in vg.graph[u]:
            if u < v:  # 避免重复绘制对称边
                line = LineString([u, v])
                ax.plot(*line.xy, color='lightblue', linewidth=0.5, zorder=1)

    # 绘制最短路径（高亮）
    if path:
        px, py = zip(*path)
        ax.plot(px, py, color='red', linewidth=2, label='Shortest Path', zorder=3)

    # 绘制起点和终点
    if start:
        ax.plot(start[0], start[1], 'go', markersize=10, label='Start', zorder=4)
    if goal:
        ax.plot(goal[0], goal[1], 'ro', markersize=10, label='Goal', zorder=4)

    ax.set_title('Visibility Graph with Shortest Path')
    ax.set_aspect('equal')
    ax.legend(loc='best')
    plt.grid(True)
    plt.savefig('visibility_graph.png', dpi=300)
    plt.show()

def compute_all_pairs_paths(vg, points):
    """Returns a dict of shortest paths and distances between all pairs of points"""
    n = len(points)
    path_dict = {}
    dist_graph = nx.Graph()

    for i in range(n):
        for j in range(i + 1, n):
            p1, p2 = tuple(points[i]), tuple(points[j])
            path = vg.dijkstra(p1, p2)
            if path is None:
                continue
            dist = sum(math.dist(path[k], path[k+1]) for k in range(len(path)-1))
            path_dict[(p1, p2)] = path
            path_dict[(p2, p1)] = list(reversed(path))
            dist_graph.add_edge(p1, p2, weight=dist)
    return path_dict, dist_graph

def solve_tsp_path(dist_graph):
    tsp_path = nx.approximation.traveling_salesman_problem(dist_graph, cycle=False)
    return tsp_path

def reconstruct_full_path(tsp_path, path_dict):
    full_path = []
    for i in range(len(tsp_path)):
        p1 = tsp_path[i]
        p2 = tsp_path[(i + 1) % len(tsp_path)]
        segment = path_dict.get((p1, p2))
        if segment:
            if full_path and segment[0] == full_path[-1]:
                full_path.extend(segment[1:])
            else:
                full_path.extend(segment)
    return full_path

def visualize(polygons, points, path):
    fig, ax = plt.subplots()
    for poly in polygons:
        xs, ys = poly.exterior.xy
        ax.plot(xs, ys, 'k-')
    for x, y in points:
        ax.plot(x, y, 'ro')
    if path:
        xs, ys = zip(*path)
        ax.plot(xs, ys, 'b--')
    plt.axis('equal')
    plt.savefig('tsp_obstacles.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    polygons = [
        Polygon([(1, 1), (1, 3), (2, 3), (2, 1)]),
        Polygon([(4, 4), (5, 4), (5, 5), (4, 5)]),
    ]
    target_points = [(0, 0), (2, 4), (6, 6), (3, 3), (3, 1)]

    vg = VisibilityGraph(polygons)
    path_dict, dist_graph = compute_all_pairs_paths(vg, target_points)
    tsp_path = solve_tsp_path(dist_graph)
    full_path = reconstruct_full_path(tsp_path, path_dict)
    print("Shortest path:", full_path)
    print("tsp_path:", tsp_path)
    visualize(polygons, target_points, full_path)

