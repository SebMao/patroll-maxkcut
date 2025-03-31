from max_k_cut import *
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform, euclidean
from itertools import combinations
from collections import defaultdict
from ortools.graph.python import min_cost_flow
from scipy.optimize import linear_sum_assignment


def create_graph(coordinates):
    """
    根据节点的二维坐标创建加权无向图，边的权重为节点间的欧氏距离。
    """
    distances = squareform(pdist(coordinates, 'euclidean'))
    G = nx.Graph()
    for i, j in combinations(range(len(coordinates)), 2):
        G.add_edge(i, j, weight=distances[i, j])
        # 将节点坐标存储为节点属性
    for i, (x, y) in enumerate(coordinates):
        G.nodes[i]['pos'] = (x, y)

    return G

def maxKcut(G, k):
    my_instance = Instance(G)
    my_instance.Params.Num_Partitions = k
    my_instance.solve()
    partition = {t: [] for t in range(k)}
    nodes = list(G.nodes())

    for i in range(len(nodes)):
        partition[my_instance.graph.nodes[i]["partition"]].append(i)
            
    return partition

def tsp_length_and_path(subset, distance_matrix):
    """ 计算子图的TSP路径长度 """
    if len(subset) < 2:
        return 0
    subgraph = distance_matrix[np.ix_(subset, subset)]
    G = nx.complete_graph(len(subset))
    for i, j in G.edges():
        G[i][j]['weight'] = subgraph[i, j]
    tsp_path_local = nx.approximation.traveling_salesman_problem(G, cycle=True)
    # print("局部编号的tsp路径:", tsp_path_local)
    # 将局部编号转换为全局编号
    index_mapping = {i: subset[i] for i in range(len(subset))}
    tsp_path_global = [index_mapping[i] for i in tsp_path_local]
    # print("全局编号的tsp路径:", tsp_path_global)
    path_length = sum(subgraph[tsp_path_local[i], tsp_path_local[i+1]] for i in range(len(tsp_path_local)-1))
    return path_length, tsp_path_global

def solve_robot_distribution_brute_force(partition, weights, m, distance_matrix):
    """
    通过遍历所有可能的机器人分配方案，找到最小的最大刷新时间。

    :param partition: 子图分割方案 (字典 {子图编号: 节点索引列表})
    :param weights: 权重列表
    :param m: 机器人总数 (假设 m ≤ n)
    :param distance_matrix: 距离矩阵
    :return: (best_robot_assignment, min_max_refresh_time)
    """
    k = len(partition)  # 子图数量
    tsp_lengths = {}
    tsp_paths = {}
    
    for j, nodes in partition.items():
        tsp_lengths[j], tsp_paths[j] = tsp_length_and_path(nodes, distance_matrix)
    # tsp_lengths = {j: tsp_length(nodes, distance_matrix) for j, nodes in partition.items()}
    
    best_time = float('inf')
    best_assignment = None
    best_partition_refresh_times = None
    best_phi1 = {}

    for assignment in generate_valid_assignments(m, k):  # 遍历所有合法分配方案
        max_time = 0
        partition_refresh_times = {}
        for j, robots in enumerate(assignment):
            nodes = partition[j]
            sorted_weights = sorted([weights[i] for i in nodes], reverse=True)  # 按权重降序排序
            print("sorted_weights:", sorted_weights)
            
            # 计算 1/w 之和 (考虑前 `robots` 个权重)
            weight_sum_inv = sum(1 / sorted_weights[t] for t in range(min(robots, len(sorted_weights))))
            print("weight_sum_inv:", weight_sum_inv)
            print("tsp_lengths[j]:", tsp_lengths[j])
            if weight_sum_inv == 0:
                continue
            refresh_time = tsp_lengths[j] / weight_sum_inv  # 计算刷新时间
            partition_refresh_times[j] = refresh_time
            print("refresh_time:", refresh_time)
            
            max_time = max(max_time, refresh_time)  # 记录最大刷新时间
        
        # 更新最优方案
        if max_time < best_time:
            best_time = max_time
            best_assignment = assignment
            best_partition_refresh_times = partition_refresh_times
            # 计算停留时间
            best_dwell_times = {i: 0 for i in range(len(weights))}  # 初始化所有节点停留时间
            for j, robots in enumerate(best_assignment):
                nodes = partition[j]
                if robots == 0:
                    continue

                sorted_nodes = sorted(nodes, key=lambda i: weights[i], reverse=True)  # 按权重排序
                phi_max_nodes = sorted_nodes[:robots]  # 取前 `robots` 个点
                phi_min_nodes = sorted_nodes[robots:]  # 其余点

                if not phi_min_nodes:
                    continue

                phi_1 = max(weights[i] for i in phi_min_nodes)  # 计算 `φ_1`
                best_phi1[j] = phi_1

                for i in phi_max_nodes:
                    phi_alpha = weights[i]
                    best_dwell_times[i] = (phi_alpha - phi_1) / (phi_alpha * phi_1) * best_time  # 计算 `δ_α`

    return best_assignment, best_time, best_dwell_times, best_partition_refresh_times, best_phi1



def generate_valid_assignments(m, k):
    """ 生成所有合法的机器人分配方案 (每个子图至少分配 1 个机器人，总数等于 m) """
    base = [1] * k  # 先分配 1 个机器人给每个子图
    remaining = m - k  # 剩余的机器人数量
    for extra in itertools.combinations_with_replacement(range(k), remaining):
        assignment = base[:]
        for i in extra:
            assignment[i] += 1
        yield tuple(assignment)


def find_best_patrolling_plan(coords, weights, m):
    """
    计算全局最优的机器人巡逻方案。

    :param coords: n个点的二维坐标
    :param weights: n个点的权重
    :param m: 机器人的个数
    :return: (best_partition, best_assignment, best_time)
    """
    n = len(coords)
    distance_matrix = np.array([[euclidean(coords[i], coords[j]) for j in range(n)] for i in range(n)])

    best_time = float('inf')
    best_partition = None
    best_assignment = None
    best_dwell_times = None
    best_refresh_times = None
    best_phi1 = None

    for k in range(2, m + 1):  # 2个子图到m个子图的情况
        G = create_graph(coords)
        partition = maxKcut(G, k)
        assignment, max_refresh_time, dwell_times, refresh_times, phi1_dict = solve_robot_distribution_brute_force(partition, weights, m, distance_matrix)
        print(max_refresh_time, k)
        if max_refresh_time < best_time:
            best_time = max_refresh_time
            best_partition = partition
            best_assignment = assignment
            best_dwell_times = dwell_times
            best_refresh_times = refresh_times
            best_phi1 = phi1_dict
    
    # 单独处理不分割的情况
    partition = {0: list(range(n))}
    assignment, max_refresh_time, dwell_times_nocut, refresh_times_nocut, phi1_dict_nocut = solve_robot_distribution_brute_force(partition, weights, m, distance_matrix)
    if max_refresh_time < best_time:
        best_time = max_refresh_time
        best_partition = partition
        best_assignment = assignment
        best_dwell_times = dwell_times_nocut
        best_refresh_times = refresh_times_nocut
        best_phi1 = phi1_dict_nocut
    
    return best_partition, best_assignment, best_time, best_dwell_times, best_refresh_times, best_phi1

def calculate_trajectory_for_leader_robot(tsp_path, node_coords, dwell_times, distance_matrix, raw_robot_positions, robot_list, speed=1):
    """
    计算0号机器人的完整周期运动轨迹，并返回分段点。
    
    :param tsp_path: 0号机器人TSP路径，包含访问节点的顺序
    :param dwell_times: 每个节点的停留时间
    :param distance_matrix: 距离矩阵
    :param raw_robot_positions: 机器人原始位置字典 {机器人编号: 二维坐标}
    :param speed: 机器人匀速移动的速度，默认为1
    :param robot_list: 子图包含机器人的编号列表
    :return: 分段点列表 [(time, position)]
    """
    trajectory = []
    current_time = 0
    n = len(tsp_path)

    # 计算列表中机器人的初始位置与tsp_path中包含节点的位置的距离矩阵
    distance_matrix_ = np.array([[euclidean(raw_robot_positions[i], node_coords[j]) for j in tsp_path[:-1]] for i in robot_list])
    # 找到距离矩阵中最小值的索引
    min_index = np.unravel_index(np.argmin(distance_matrix_), distance_matrix_.shape)
    # 获取最小值的行和列索引
    min_row_index, min_col_index = min_index
    # 获取最小值对应的机器人编号和节点编号
    leader_index = robot_list[min_row_index]
    # 将tsp_path进行旋转，以leader_node_index为起点，[1,2,3,4,1] -> [2,3,4,1,2]
    tsp_path = shift_to_index(tsp_path, min_col_index)
    for i in range(n):
        node_current = tsp_path[i]
        node_next = tsp_path[(i + 1) % n]  # 下一个节点（循环）
        
        # 获取当前节点和下一个节点的坐标
        coord_current = coords[node_current]        
        # 计算从node_current到node_next的移动时间
        distance = distance_matrix[node_current][node_next]
        travel_time = distance / speed
        
        # 记录当前节点和对应的停留时间
        trajectory.append((current_time, coord_current))  # 当前节点的时间戳和位置
        current_time += dwell_times[node_current]  # 停留时间
        trajectory.append((current_time, coord_current))  # 停留后的时间戳
        
        # 移动到下一个节点
        current_time += travel_time  # 机器人移动的时间
    return trajectory, leader_index

def calculate_trajectory_for_robot(tsp_path, node_coords, dwell_times, distance_matrix, initial_position, speed=1):
    """
    计算跟随机器人的完整周期运动轨迹，并返回分段点。
    
    :param tsp_path: 路径，包含访问节点的顺序
    :param node_coords: 节点坐标
    :param dwell_times: 每个节点的停留时间
    :param distance_matrix: 距离矩阵
    :param initial_position: 机器人的初始位置
    :param speed: 机器人匀速移动的速度，默认为1
    :return: 分段点列表 [(time, position)]
    """
    trajectory = []
    current_time = 0.0
    n = len(tsp_path)
    # 寻找tsp_path中下一个目标位置
    next_target, target_node_id = find_next_target(initial_position, tsp_path, node_coords)
    # 计算初始位置到下一个目标位置的距离
    initial_distance = euclidean(initial_position, next_target)
    trajectory.append((current_time, initial_position))  # 初始位置的时间戳和位置
    current_time += initial_distance / speed  # 移动到下一个目标位置的时间
    shifted_tsp_path = shift_to_index(tsp_path, target_node_id)
    for i in range(n-1):
        node_current = shifted_tsp_path[i]
        node_next = shifted_tsp_path[(i + 1) % n]  # 下一个节点（循环）
        # 获取当前节点和下一个节点的坐标
        coord_current = node_coords[node_current]
        coord_next = node_coords[node_next]
        # 计算从node_current到node_next的移动时间
        distance = distance_matrix[node_current][node_next]
        travel_time = distance / speed
        # 记录当前节点和对应的停留时间
        trajectory.append((current_time, coord_current))  # 当前节点的时间戳和位置
        current_time += dwell_times[node_current]  # 停留时间
        trajectory.append((current_time, coord_current))  # 停留后的时间戳
        # 移动到下一个节点
        current_time += travel_time  # 机器人移动的时间

    # 计算返回初始位置的时间
    distance_back = euclidean(coord_current, initial_position)
    travel_time_back = distance_back / speed
    current_time += travel_time_back  # 移动到初始位置的时间
    trajectory.append((current_time, initial_position))  # 当前节点的时间戳和位置
    return trajectory

def shift_to_index(arr, index):
    if not arr:
        return arr
    n = len(arr)
    truncated_arr = arr[:-1]
    shifted_arr = truncated_arr[index:] + truncated_arr[:index]
    shifted = shifted_arr + [shifted_arr[0]]
    return shifted

def interpolate_position(trajectory, time):
    """
    通过插值计算在给定时间点的位置。
    
    :param trajectory: 分段轨迹，按时间排序
    :param time: 给定时间
    :return: 插值后的位置信息
    """
    max_time = trajectory[-1][0]
    if time >= max_time:
        time = time % max_time
    for i in range(1, len(trajectory)):
        t_j, p_j = trajectory[i - 1]
        t_j1, p_j1 = trajectory[i]
        
        # 如果时间位于这一段区间内，进行线性插值
        if t_j <= time < t_j1:
            # 计算在区间 [t_j, t_j1] 内的比例
            ratio = (time - t_j) / (t_j1 - t_j)
            # 插值位置
            interpolated_position = (
                p_j[0] + ratio * (p_j1[0] - p_j[0]),  # x坐标插值
                p_j[1] + ratio * (p_j1[1] - p_j[1])   # y坐标插值
            )
            return interpolated_position
    
    return trajectory[-1][1]  # 如果时间超出轨迹范围，返回最后一个位置

def get_robot_initial_position(trajectory, m, T, phi_1, i):
    """
    计算机器人 i 在周期 T 内的初始位置。
    
    :param trajectory: 0号机器人的完整轨迹
    :param m: 总机器人数
    :param T: 子图的刷新时间
    :param phi_1: 子图中的第m+1大权重
    :param i: 机器人编号
    :return: 机器人 i 在周期 T 内的初始位置
    """
    # 计算时间点
    time_i = (m - i) * T / phi_1
    
    # 查找时间点位于哪一段
    return interpolate_position(trajectory, time_i)

def compute_subgraph_centroids(coords, partition):
    """ 计算每个子图的几何中心 """
    centroids = {}
    for subgraph_id, nodes in partition.items():
        sub_coords = np.array([coords[i] for i in nodes])
        centroids[subgraph_id] = np.mean(sub_coords, axis=0)  # 计算质心
    return centroids

def optimal_robot_assignment_hungarian(raw_robot_positions, best_assignment, partition, coords):
    """ 使用匈牙利算法找到最优机器人分配 """
    centroids = compute_subgraph_centroids(coords, partition)
    subgraph_ids = list(centroids.keys())
    robot_ids = list(raw_robot_positions.keys())
    # 计算机器人到子图质心的距离
    distances = []
    for robot_id, robot_pos in raw_robot_positions.items():
        distances.append([euclidean(robot_pos, centroids[subgraph_id]) for subgraph_id in best_partition.keys()])

    # 重复子图节点以满足 n_i 需求（例如 n_i=[2,1] → 子图0出现2次）
    expanded_indices = []
    for i, count in enumerate(best_assignment):
        expanded_indices.extend([i] * count)
    
    # 构建代价矩阵（机器人 x 展开后的子图需求）
    cost_matrix = []
    for robot_dists in distances:
        cost_matrix.append([robot_dists[i] for i in expanded_indices])
    
    # 使用匈牙利算法求解最优分配
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # 生成最终的 robot_assignment
    robot_assignment = {subgraph_id: [] for subgraph_id in subgraph_ids}
    for robot_idx, expanded_idx in zip(row_ind, col_ind):
        subgraph_idx = expanded_indices[expanded_idx]
        robot_assignment[subgraph_idx].append(robot_idx)    
    
    return robot_assignment

def optimal_robot_assignment_min_cost_flow(raw_robot_positions, best_assignment, best_partition, coords):
    """
    使用最小费用流算法找到最优的机器人分配方案。
    
    :param raw_robot_positions: 机器人原始坐标，字典{机器人编号: (x, y)}
    :param best_assignment: 每个子图所需的机器人数量，字典{子图编号: 机器人数量}
    :param best_partition: 每个子图包含的节点标号，字典{子图编号: 节点索引列表}
    :param coords: 所有节点的坐标，列表[(x1, y1), (x2, y2), ...]
    :return: 最优的机器人分配方案，字典{子图编号: 机器人索引列表}
    """

    scale_factor = 1_000_000
    centroids = compute_subgraph_centroids(coords, best_partition)
        
    # 计算机器人到子图质心的距离
    distances = []
    for robot_id, robot_pos in raw_robot_positions.items():
        distances.append([euclidean(robot_pos, centroids[subgraph_id]) for subgraph_id in best_partition.keys()])
    scaled_distances = [[int(d * scale_factor) for d in row] for row in distances]
    smcf = min_cost_flow.SimpleMinCostFlow()
    m = len(distances)
    k = len(best_assignment)  # 子图数量
    n_i = best_assignment  # 每个子图所需的机器人数量
    # 节点编号：0源，1~m机器人，m+1~m+k子图，m+k+1汇
    source = 0
    sink = m + k + 1
    robot_nodes = list(range(1, m+1))
    subgraph_nodes = list(range(m+1, m+k+1))

        # 源到机器人
    for j in robot_nodes:
        smcf.add_arc_with_capacity_and_unit_cost(source, j, 1, 0)
    
    # 机器人到子图
    for idx, j in enumerate(robot_nodes):
        for i in range(k):
            smcf.add_arc_with_capacity_and_unit_cost(j, subgraph_nodes[i], 1, scaled_distances[idx][i])
    
    # 子图到汇
    for i in range(k):
        smcf.add_arc_with_capacity_and_unit_cost(subgraph_nodes[i], sink, n_i[i], 0)
    
    # 设置供需
    smcf.set_node_supply(source, m)
    smcf.set_node_supply(sink, -m)
    for node in robot_nodes + subgraph_nodes:
        smcf.set_node_supply(node, 0)
    
    status = smcf.solve()
    if status != smcf.OPTIMAL:
        return None
    
    # 解析结果
    subgraph_ids = list(centroids.keys())
    robot_assignment = {subgraph_id: [] for subgraph_id in subgraph_ids}
    for arc in range(smcf.num_arcs()):
        if smcf.tail(arc) in robot_nodes and smcf.flow(arc) > 0:
            robot = smcf.tail(arc) - 1
            subgraph = smcf.head(arc) - (m+1)
            robot_assignment[subgraph].append(robot)
    return robot_assignment

def solve_intial_position_(partition, coords, robot_assignment, global_dwell_times, refresh_times, phi1_dict):
    """
    计算每个机器人的初始位置。
    
    :param partition: 子图分割方案 (字典 {子图编号: 节点索引列表})
    :param coords: n个点的二维坐标
    :param robot_assignment: 每个子图的机器人的编号列表
    :param global_dwell_times: 每个节点的停留时间
    :param refresh_times: 每个子图的刷新时间
    :param phi1_dict: 每个子图的第m+1大权重
    :return: 机器人初始位置字典 {机器人编号: 二维坐标}
    """

    initial_positions = {}
    tsp_lengths = {}
    tsp_paths = {}
    n = len(coords)
    distance_matrix_ = np.array([[euclidean(coords[i], coords[j]) for j in range(n)] for i in range(n)])

    for j, nodes in partition.items():
        tsp_lengths[j], tsp_paths[j] = tsp_length_and_path(nodes, distance_matrix_)
    for subgraph_id, robots in robot_assignment.items():
        nodes = partition[subgraph_id]
        leader_robot_traj, leader_index = calculate_trajectory_for_leader_robot(tsp_path=tsp_paths[subgraph_id], node_coords=coords, dwell_times=global_dwell_times, distance_matrix=distance_matrix_, raw_robot_positions=raw_robot_positions, robot_list=robots)
        t = 1
        for robot_id in robots:
            if robot_id == leader_index:
                initial_positions[robot_id] = leader_robot_traj[0][1]
            else:
                initial_positions[robot_id] = get_robot_initial_position(leader_robot_traj, len(robots), refresh_times[subgraph_id], phi1_dict[subgraph_id], t)
                t += 1
    return initial_positions
def assign_robots_to_subgraphs(raw_robot_positions, best_assignment, best_partition, coords):
    """
    将机器人分配到最近的子图，考虑机器人到子图节点的最短距离，返回最优的机器人分配方案。
    
    :param raw_robot_positions: 机器人原始坐标，字典{机器人编号: (x, y)}
    :param best_assignment: 每个子图所需的机器人数量，字典{子图编号: 机器人数量}
    :param best_partition: 每个子图包含的节点标号，字典{子图编号: 节点索引列表}
    :param coords: 所有节点的坐标，列表[(x1, y1), (x2, y2), ...]
    
    :return: 最优的机器人分配方案，字典{子图编号: 机器人索引列表}
    """
    # 计算每个子图的节点坐标列表
    subgraph_node_coords = {}
    for subgraph_id, nodes in best_partition.items():
        sub_coords = np.array([coords[i] for i in nodes])
        subgraph_node_coords[subgraph_id] = sub_coords
    
    # 计算每个机器人到子图中节点的最短距离
    robot_distances = defaultdict(list)
    for robot_id, position in raw_robot_positions.items():
        for subgraph_id, sub_coords in subgraph_node_coords.items():
            # 计算机器人到每个子图节点的距离
            distances_to_nodes = np.linalg.norm(sub_coords - position, axis=1)
            min_distance = np.min(distances_to_nodes)
            robot_distances[subgraph_id].append((min_distance, robot_id))
    
    # 根据距离对每个子图的机器人进行排序
    for subgraph_id in robot_distances:
        robot_distances[subgraph_id].sort()

    # 分配机器人到子图，避免机器人重复分配
    robot_assignment = defaultdict(list)
    assigned_robots = set()  # 用一个集合记录已分配的机器人

        # 假设 best_assignment = (2, 3)  # 代表子图 0 需要 2 个机器人，子图 1 需要 3 个机器人
    for subgraph_id, num_robots in enumerate(best_assignment):  # 按索引访问
        assigned_robots_for_subgraph = 0
        for _, robot_id in robot_distances[subgraph_id]:
            # 如果机器人没有被分配到任何子图，则分配给当前子图
            if robot_id not in assigned_robots and assigned_robots_for_subgraph < num_robots:
                robot_assignment[subgraph_id].append(robot_id)
                assigned_robots.add(robot_id)  # 记录该机器人已经分配
                assigned_robots_for_subgraph += 1
            if assigned_robots_for_subgraph >= num_robots:
                break  # 当前子图的机器人数量已经足够


    return robot_assignment


def plot_patrolling_plan(coords, partition, robot_assignment, raw_robot_positions, initial_positions):
    """
    可视化最优的巡逻方案。
    
    :param coords: n个点的二维坐标
    :param partition: {子图编号: 节点索引列表}
    :param robot_assignment: {子图编号: 机器人索引列表}
    :param raw_robot_positions: 机器人原始位置字典 {机器人编号: 二维坐标}
    :param initial_positions: 机器人任务初始位置字典 {机器人编号: 二维坐标}
    """
    plt.figure(figsize=(10, 8))
    colors = plt.cm.get_cmap("tab10", len(partition))

    for subgraph_id, nodes in partition.items():
        sub_coords = np.array([coords[i] for i in nodes])
        
        # 绘制节点
        plt.scatter(sub_coords[:, 0], sub_coords[:, 1], color=colors(subgraph_id), label=f"Subgraph {subgraph_id}")
        if len(nodes) == 1:
            continue
        
        # 计算TSP路径并绘制
        G = nx.complete_graph(len(nodes))
        for i, j in G.edges():
            G[i][j]['weight'] = np.linalg.norm(sub_coords[i] - sub_coords[j])
        tsp_path = nx.approximation.traveling_salesman_problem(G, cycle=True)
        for i in range(len(tsp_path) - 1):
            plt.plot([sub_coords[tsp_path[i]][0], sub_coords[tsp_path[i + 1]][0]],
                     [sub_coords[tsp_path[i]][1], sub_coords[tsp_path[i + 1]][1]],
                     color=colors(subgraph_id), linestyle='--', alpha=0.7)
        
        # 标注点索引
        for idx, (x, y) in zip(nodes, sub_coords):
            plt.text(x, y, str(idx), fontsize=12, ha='right', color='black')
        
        # 在子图中心标注机器人数量
        centroid = np.mean(sub_coords, axis=0)
        num_robots = len(robot_assignment.get(subgraph_id, []))
        plt.text(centroid[0], centroid[1], f"{num_robots} robots", fontsize=12,
                 ha='center', va='center', color='red', bbox=dict(facecolor='white', alpha=0.6))
    
    # 可视化机器人原始位置
    for robot_id, position in raw_robot_positions.items():
        assigned_subgraph = None
        for subgraph_id, robots in robot_assignment.items():
            if robot_id in robots:
                assigned_subgraph = subgraph_id
                break
        color = colors(assigned_subgraph) if assigned_subgraph is not None else 'black'
        plt.scatter(position[0], position[1], color=color, marker='x', s=100, label=f"Robot {robot_id} (raw)")
        plt.text(position[0], position[1], f"R{robot_id}", fontsize=10, ha='left', color=color)
    
    # 可视化机器人初始位置
    for robot_id, position in initial_positions.items():
        plt.scatter(position[0], position[1], color='black', marker='o', s=100, label=f"Robot {robot_id} initial position")
        plt.text(position[0], position[1], f"Robot {robot_id}", fontsize=12, ha='left', color='black')

    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.title("Patrolling Plan")
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig("patrolling_plan.png")

def get_robots_trajectories(robot_assignment, initial_positions, coords, dwell_times, partition, tsp_paths, distance_matrix):
    """
    计算每个机器人的完整周期运动轨迹，并返回分段点。
    
    :param robot_assignment: {子图编号: 机器人索引列表}
    :param initial_positions: 机器人初始位置字典 {机器人编号: 二维坐标}
    :param coords: n个点的二维坐标
    :param dwell_times: 每个节点的停留时间
    :param partition: {子图编号: 节点索引列表}
    :param tsp_paths: {子图编号: TSP路径}
    :param distance_matrix: 距离矩阵
    
    :return: 分段点列表 [(time, position)]
    """
    trajectories = {}
    for subgraph_id, robots in robot_assignment.items():
        # 将子图中所有节点的位置组合成一个字典
        sub_coords = {}
        for node in partition[subgraph_id]:
            sub_coords[node] = coords[node]
        for robot_id in robots:
            robot_traj = calculate_trajectory_for_robot(
                tsp_path=tsp_paths[subgraph_id],
                node_coords=sub_coords,
                dwell_times=dwell_times,
                distance_matrix=distance_matrix,
                initial_position=initial_positions[robot_id],
                speed=1
            )
            trajectories[robot_id] = robot_traj
    return trajectories

def find_next_target(robot_pos, tsp_path, node_coords):
    n = len(tsp_path)
    
    # 找到机器人是否恰好在某个节点上
    for id, node_pos in node_coords.items():
        if np.allclose(robot_pos, node_pos, atol=1e-6):  # 机器人刚好在节点上
            node_index = tsp_path.index(id)  # 在TSP路径中的索引
            return node_coords[tsp_path[(node_index) % n]], node_index  # 下一个目标点
    
    # 机器人不在节点上，找到所在的路径段
    for i in range(len(tsp_path) - 1):  # 遍历TSP路径
        idx1, idx2 = tsp_path[i], tsp_path[i + 1]
        A, B = np.array(node_coords[idx1]), np.array(node_coords[idx2])
        
        # 检查机器人是否在线段 AB 上
        AB = B - A
        AP = np.array(robot_pos) - A
        proj = np.dot(AP, AB) / np.dot(AB, AB)  # 计算投影比例
        cross_product = np.cross(AB, AP)  # 计算叉积
        if np.linalg.norm(cross_product) > 1e-6:  # 机器人不在AB线段上
            continue
        if 0 <= proj <= 1:  # 机器人确实在AB线段上
            return tuple(node_coords[idx2]), tsp_path.index(idx2)  # 返回 B 作为下一个目标
    
    return None  # 机器人不在任何TSP路径上

# coords = [(4, 1), (2, 3), (5, 5), (8, 8), (12, 3), (6, 9), (14, 10), (7, 2),(10,5), (14, 6)]
coords = [(4, 1), (2, 3), (5, 5), (8, 8), (12, 3), (6, 9), (14, 10), (7, 2),(10,5)]

# raw_robot_positions = {0: (1, 2), 1: (3, 4), 2: (5, 6), 3: (7,8)}
raw_robot_positions = {0: (8, 4), 1: (7, 5), 2: (5, 6)}

# weights = [1.5, 2.0, 1.2, 1.9, 2.5, 1.8, 2.2, 3.0, 6.0, 4.0]
weights = [1.5, 2.0, 1.2, 1.9, 2.5, 1.8, 2.2, 3.0, 6.0]

# 机器人数量
m = max(raw_robot_positions.keys()) + 1

best_partition, best_assignment, best_time, best_dwell_times, best_refresh_times, best_phi1_dict = find_best_patrolling_plan(coords, weights, m)
robot_assignment = optimal_robot_assignment_min_cost_flow(raw_robot_positions, best_assignment, best_partition, coords)
initial_positions = solve_intial_position_(best_partition, coords, robot_assignment, best_dwell_times, best_refresh_times, best_phi1_dict)

# 每个机器人的轨迹
tsp_paths = {}
distance_matrix_ = np.array([[euclidean(coords[i], coords[j]) for j in range(len(coords))] for i in range(len(coords))])

for j, nodes in best_partition.items():
    _, tsp_paths[j] = tsp_length_and_path(nodes, distance_matrix_)
robot_trajectories = get_robots_trajectories(robot_assignment, initial_positions, coords, best_dwell_times, best_partition, tsp_paths, distance_matrix_)
plot_patrolling_plan(coords, best_partition, robot_assignment, raw_robot_positions, initial_positions)

print("最优分割方案:", best_partition)
print("最优机器人分配方案 (子图编号: 机器人数量):")
for i, robots in enumerate(best_assignment):
    print(f" - 子图 {i}: {robots} 个机器人")
print("机器人分配方案:", robot_assignment)

for i in range(m):
    initial_position = initial_positions[i]
    print(f"机器人 {i} 的initial position: {initial_position}")
print("最小的最大刷新时间:", best_time)
for i, time in best_dwell_times.items():
    print(f"节点 {i} 的停留时间: {time}")
for i, traj in robot_trajectories.items():
    print(f"机器人 {i} 的完整周期运动轨迹: {traj}")
    print("\n")


