from max_k_cut import *
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform, euclidean
from itertools import combinations

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


def plot_patrolling_plan(coords, partition, robot_assignment):
    """
    可视化最优的巡逻方案。
    
    :param coords: n个点的二维坐标
    :param partition: {子图编号: 节点索引列表}
    :param robot_assignment: 每个子图的机器人数量
    """
    plt.figure(figsize=(10, 8))
    colors = plt.cm.get_cmap("tab10", len(partition))

    for subgraph_id, nodes in partition.items():
        sub_coords = np.array([coords[i] for i in nodes])
        
        # 绘制节点
        plt.scatter(sub_coords[:, 0], sub_coords[:, 1], color=colors(subgraph_id), label=f"{subgraph_id}")
        if(len(nodes) == 1):
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
        plt.text(centroid[0], centroid[1], f"{robot_assignment[subgraph_id]} robot", fontsize=12,
                 ha='center', va='center', color='red', bbox=dict(facecolor='white', alpha=0.6))
    
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.title("patrolling plan")
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig("patrolling_plan.png")



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

def calculate_trajectory_for_robot_0(tsp_path, dwell_times, distance_matrix, speed=1):
    """
    计算0号机器人的完整周期运动轨迹，并返回分段点。
    
    :param tsp_path: 0号机器人TSP路径，包含访问节点的顺序
    :param dwell_times: 每个节点的停留时间
    :param distance_matrix: 距离矩阵
    :param speed: 机器人匀速移动的速度，默认为1
    :return: 分段点列表 [(time, position)]
    """
    trajectory = []
    current_time = 0
    n = len(tsp_path)
    
    # 0号机器人从路径的第一个节点开始
     # 0号机器人从路径的第一个节点开始
    for i in range(n):
        node_current = tsp_path[i]
        node_next = tsp_path[(i + 1) % n]  # 下一个节点（循环）
        
        # 获取当前节点和下一个节点的坐标
        coord_current = coords[node_current]
        coord_next = coords[node_next]
        
        # 计算从node_current到node_next的移动时间
        distance = distance_matrix[node_current][node_next]
        travel_time = distance / speed
        
        # 记录当前节点和对应的停留时间
        trajectory.append((current_time, coord_current))  # 当前节点的时间戳和位置
        current_time += dwell_times[node_current]  # 停留时间
        trajectory.append((current_time, coord_current))  # 停留后的时间戳
        
        # 移动到下一个节点
        current_time += travel_time  # 机器人移动的时间
        # trajectory.append((current_time, coord_next))  # 记录到达下一个节点的时间戳
    # node_current = tsp_path[0]
    # coord_current = coords[node_current]
    # trajectory.append((current_time, coord_current))
    
    return trajectory

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

def solve_intial_position(partition, coords, assignment, global_dwell_times, refresh_times, phi1_dict):
    initial_positions = {}
    tsp_lengths = {}
    tsp_paths = {}
    n = len(coords)
    distance_matrix_ = np.array([[euclidean(coords[i], coords[j]) for j in range(n)] for i in range(n)])

    for j, nodes in partition.items():
        tsp_lengths[j], tsp_paths[j] = tsp_length_and_path(nodes, distance_matrix_)
    global_index = 0
    for i, robots in enumerate(assignment):
        nodes = partition[i]
        # dwell_times_par = [global_dwell_times[i] for i in nodes]
        robot_0_traj = calculate_trajectory_for_robot_0(tsp_path=tsp_paths[i], dwell_times=global_dwell_times, distance_matrix = distance_matrix_ )
        # initial_positions[global_index] = robot_initial_position
        for t in range(robots):
            initial_positions[global_index+t] = get_robot_initial_position(robot_0_traj, robots, refresh_times[i], phi1_dict[i],t)
        global_index += robots
    return initial_positions




coords = [(4, 1), (2, 3), (5, 5), (8, 8), (12, 3), (6, 9), (14, 10), (7, 2),(10,5)]
weights = [1.5, 2.0, 1.2, 1.9, 2.5, 1.8, 2.2, 3.0, 6.0]
m = 3

best_partition, best_assignment, best_time, best_dwell_times, best_refresh_times, best_phi1_dict = find_best_patrolling_plan(coords, weights, m)
plot_patrolling_plan(coords, best_partition, {i: robots for i, robots in enumerate(best_assignment)})
initial_positions = solve_intial_position(best_partition, coords, best_assignment, best_dwell_times, best_refresh_times, best_phi1_dict)

print("最优分割方案:", best_partition)
print("最优机器人分配方案 (子图编号: 机器人数量):")
for i, robots in enumerate(best_assignment):
    print(f" - 子图 {i}: {robots} 个机器人")

for i in range(m):
    initial_position = initial_positions[i]
    print(f"机器人 {i} 的initial position: {initial_position}")
print("最小的最大刷新时间:", best_time)
for i, time in best_dwell_times.items():
    print(f"节点 {i} 的停留时间: {time}")



