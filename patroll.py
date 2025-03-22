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

def tsp_length(subset, distance_matrix):
    """ 计算子图的TSP路径长度 """
    if len(subset) < 2:
        return 0
    subgraph = distance_matrix[np.ix_(subset, subset)]
    G = nx.complete_graph(len(subset))
    for i, j in G.edges():
        G[i][j]['weight'] = subgraph[i, j]
    tsp_path = nx.approximation.traveling_salesman_problem(G, cycle=True)
    return sum(subgraph[tsp_path[i], tsp_path[i+1]] for i in range(len(tsp_path)-1))

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
    tsp_lengths = {j: tsp_length(nodes, distance_matrix) for j, nodes in partition.items()}
    
    best_time = float('inf')
    best_assignment = None

    for assignment in generate_valid_assignments(m, k):  # 遍历所有合法分配方案
        max_time = 0
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
            print("refresh_time:", refresh_time)
            
            max_time = max(max_time, refresh_time)  # 记录最大刷新时间
        
        # 更新最优方案
        if max_time < best_time:
            best_time = max_time
            best_assignment = assignment

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

                for i in phi_max_nodes:
                    phi_alpha = weights[i]
                    best_dwell_times[i] = (phi_alpha - phi_1) / (phi_alpha * phi_1) * best_time  # 计算 `δ_α`


    return best_assignment, best_time, best_dwell_times



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

    for k in range(2, m + 1):  # 2个子图到m个子图的情况
        G = create_graph(coords)
        partition = maxKcut(G, k)
        assignment, max_refresh_time, dwell_times = solve_robot_distribution_brute_force(partition, weights, m, distance_matrix)
        print(max_refresh_time, k)
        if max_refresh_time < best_time:
            best_time = max_refresh_time
            best_partition = partition
            best_assignment = assignment
            best_dwell_times = dwell_times
    
    # 单独处理不分割的情况
    partition = {0: list(range(n))}
    assignment, max_refresh_time, dwell_times_nocut = solve_robot_distribution_brute_force(partition, weights, m, distance_matrix)
    if max_refresh_time < best_time:
        best_time = max_refresh_time
        best_partition = partition
        best_assignment = assignment
        best_dwell_times = dwell_times_nocut
    
    return best_partition, best_assignment, best_time, best_dwell_times

coords = [(4, 1), (2, 3), (5, 5), (8, 8), (12, 3), (6, 9), (14, 10), (7, 2),(10,5)]
weights = [1.5, 2.0, 1.2, 1.9, 2.5, 1.8, 2.2, 3.0, 6.0]
m = 3

best_partition, best_assignment, best_time, best_dwell_times = find_best_patrolling_plan(coords, weights, m)
print("最优分割方案:", best_partition)
print("最优机器人分配方案 (子图编号: 机器人数量):")
for i, robots in enumerate(best_assignment):
    print(f" - 子图 {i}: {robots} 个机器人")
print("最小的最大刷新时间:", best_time)
for i, time in best_dwell_times.items():
    print(f"节点 {i} 的停留时间: {time}")
plot_patrolling_plan(coords, best_partition, {i: robots for i, robots in enumerate(best_assignment)})


# print(f'节点的分区分配：{partition_assignment}')
# print(f'最大 k 切割值：{max_cut_value}')
