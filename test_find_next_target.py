import numpy as np

def euclidean_distance(p1, p2):
    return np.linalg.norm(np.array(p1) - np.array(p2))

def find_next_target(robot_pos, tsp_path, node_coords):
    n = len(tsp_path)
    
    # 找到机器人是否恰好在某个节点上
    for i, node in enumerate(node_coords):
        if np.allclose(robot_pos, node, atol=1e-6):  # 机器人刚好在节点上
            node_index = tsp_path.index(i + 1)  # 在TSP路径中的索引
            return node_coords[tsp_path[(node_index) % n] - 1]  # 下一个目标点
    
    # 机器人不在节点上，找到所在的路径段
    for i in range(len(tsp_path) - 1):  # 遍历TSP路径
        idx1, idx2 = tsp_path[i] - 1, tsp_path[i + 1] - 1
        A, B = np.array(node_coords[idx1]), np.array(node_coords[idx2])
        
        # 检查机器人是否在线段 AB 上
        AB = B - A
        AP = np.array(robot_pos) - A
        proj = np.dot(AP, AB) / np.dot(AB, AB)  # 计算投影比例
        
        cross_product = np.cross(AB, AP)

        if np.allclose(cross_product, atol=1e-3) and 0 <= proj <= 1:  # 机器人确实在AB线段上
            return tuple(node_coords[idx2])  # 返回 B 作为下一个目标
    
    return None  # 机器人不在任何TSP路径上

# 示例数据
tsp_path = [1, 2, 3, 4, 1]
node_coords = [(1, 2), (3, 4), (5, 6), (7, 8)]
robot_pos = (1, 2)

next_target = find_next_target(robot_pos, tsp_path, node_coords)
print("机器人下一个目标:", next_target)
