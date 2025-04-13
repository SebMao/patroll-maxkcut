import math
import heapq
from shapely.geometry import LineString, Polygon, Point
import matplotlib.pyplot as plt


class VisibilityGraph:
    def __init__(self, polygons, epsilon=1e-9):
        self.polygons = polygons
        self.vertices = self.extract_vertices()
        self.edges = self.extract_edges()
        self.epsilon = epsilon
        self.graph = self.build_visibility_graph()

    def extract_vertices(self):
        vertices = []
        for poly in self.polygons:
            coords = list(poly.exterior.coords)[:-1]  # Remove duplicate last point
            vertices.extend([Point(p) for p in coords])
        return vertices

    def extract_edges(self):
        edges = set()
        for poly in self.polygons:
            coords = list(poly.exterior.coords)
            for i in range(len(coords) - 1):
                edges.add(LineString([coords[i], coords[i + 1]]))
        return edges

    def is_visible(self, p1, p2):
        segment = LineString([p1, p2])

        # Check if the segment intersects any polygon edge (excluding shared endpoints)
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
        graph = {tuple(v.coords)[0]: {} for v in self.vertices}
        points = list(graph.keys())
        n = len(points)

        for i in range(n):
            for j in range(i + 1, n):
                if self.is_visible(points[i], points[j]):
                    dist = math.dist(points[i], points[j])
                    graph[points[i]][points[j]] = dist
                    graph[points[j]][points[i]] = dist  # Symmetric
        return graph

    def add_point_to_graph(self, point):
        """Adds start or goal point to visibility graph."""
        point_t = tuple(point)
        self.graph[point_t] = {}
        for v in self.graph.keys():
            if v == point_t:
                continue
            if self.is_visible(point, v):
                dist = math.dist(point, v)
                self.graph[point_t][v] = dist
                self.graph[v][point_t] = dist

    def shortest_path(self, start, goal):
        start = tuple(start)
        goal = tuple(goal)

        # Add start and goal into the graph
        self.add_point_to_graph(start)
        self.add_point_to_graph(goal)

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

if __name__ == "__main__":
    # 定义障碍物
    polygons = [
        Polygon([(1, 1), (1, 3), (3, 3), (3, 1)]),
        Polygon([(4, 4), (4, 6), (6, 6), (6, 4)]),
    ]

    vg = VisibilityGraph(polygons)
    start = (0, 0)
    goal = (7, 7)

    path = vg.shortest_path(start, goal)
    plot_visibility_graph(vg, path=path, start=start, goal=goal)
    print("Shortest path:", path)
