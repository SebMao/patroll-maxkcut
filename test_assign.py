from ortools.graph.python import min_cost_flow
from scipy.spatial.distance import pdist, squareform, euclidean
import numpy as np
from scipy.optimize import linear_sum_assignment

# def allocate_robots(distances, n_i, scaled_factor = 1_000_000):
#     smcf = min_cost_flow.SimpleMinCostFlow()
#     m = len(distances)
#     k = len(n_i)
    
#     # 节点编号：0源，1~m机器人，m+1~m+k子图，m+k+1汇
#     source = 0
#     sink = m + k + 1
#     robot_nodes = list(range(1, m+1))
#     subgraph_nodes = list(range(m+1, m+k+1))
#     scaled_distances = [
#         [int(d * scaled_factor) for d in row] 
#         for row in distances
#     ]

#     # 源到机器人
#     for j in robot_nodes:
#         smcf.add_arc_with_capacity_and_unit_cost(source, j, 1, 0)
    
#     # 机器人到子图
#     for idx, j in enumerate(robot_nodes):
#         for i in range(k):
#             smcf.add_arc_with_capacity_and_unit_cost(j, subgraph_nodes[i], 1, scaled_distances[idx][i])
    
#     # 子图到汇
#     for i in range(k):
#         smcf.add_arc_with_capacity_and_unit_cost(subgraph_nodes[i], sink, n_i[i], 0)
    
#     # 设置供需
#     smcf.set_node_supply(source, m)
#     smcf.set_node_supply(sink, -m)
#     for node in robot_nodes + subgraph_nodes:
#         smcf.set_node_supply(node, 0)
    
#     status = smcf.solve()
#     if status != smcf.OPTIMAL:
#         return None
    
#     # 解析结果
#     allocation = {}
#     for arc in range(smcf.num_arcs()):
#         if smcf.tail(arc) in robot_nodes and smcf.flow(arc) > 0:
#             robot = smcf.tail(arc) - 1
#             subgraph = smcf.head(arc) - (m+1)
#             allocation[robot] = subgraph
#     return allocation

def allocate_robots(distances, n_i):
    # 重复子图节点以满足 n_i 需求（例如 n_i=[2,1] → 子图0出现2次）
    expanded_indices = []
    for i, count in enumerate(n_i):
        expanded_indices.extend([i] * count)
    
    # 构建代价矩阵（机器人 x 展开后的子图需求）
    cost_matrix = []
    for robot_dists in distances:
        cost_matrix.append([robot_dists[i] for i in expanded_indices])
    
    # 求解
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # 映射回原始子图
    allocation = {}
    for robot, expanded_idx in zip(row_ind, col_ind):
        subgraph = expanded_indices[expanded_idx]
        allocation[robot] = subgraph
    
    return allocation
# 示例用法

coords = [(4, 1), (2, 3), (5, 5), (8, 8), (12, 3), (6, 9), (14, 10), (7, 2),(10,5)]
raw_robot_positions = {0: (8, 4), 1: (7, 5), 2: (5, 6)}
centroids = {}
partition = {
    0: [0,1,2,7],
    1: [3,4,5,6,8]
}
for subgraph_id, nodes in partition.items():
    sub_coords = np.array([coords[i] for i in nodes])
    centroids[subgraph_id] = np.mean(sub_coords, axis=0)  # 计算质心
# 计算机器人到子图质心的距离
distances = []
for robot_id, pos in raw_robot_positions.items():
    distances.append([euclidean(pos, centroids[subgraph_id]) for subgraph_id in partition.keys()])
n_i = [1, 2]  
print(allocate_robots(distances, n_i))  # 输出：{0: 0, 1: 1}