from max_k_cut import *
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform, euclidean
from itertools import combinations, combinations_with_replacement
import itertools
import matplotlib.pyplot as plt
from collections import defaultdict
from ortools.graph.python import min_cost_flow
from scipy.optimize import linear_sum_assignment
from matplotlib.animation import FuncAnimation, PillowWriter
from distance_caculation import DistanceCalculator, compute_all_pairs_paths
from shapely.geometry import Polygon


class GraphCreator:
    """
    用于创建和处理图结构的类
    """
    @staticmethod
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

    @staticmethod
    def maxKcut(G, k):
        my_instance = Instance(G)
        my_instance.Params.Num_Partitions = k
        my_instance.solve()
        partition = {t: [] for t in range(k)}
        nodes = list(G.nodes())

        for i in range(len(nodes)):
            partition[my_instance.graph.nodes[i]["partition"]].append(i)
                
        return partition

    @staticmethod
    def tsp_length_and_path(subset, distance_matrix):
        """ 计算子图的TSP路径长度 """
        if len(subset) < 2:
            return 0, subset
        subgraph = distance_matrix[np.ix_(subset, subset)]
        G = nx.complete_graph(len(subset))
        for i, j in G.edges():
            G[i][j]['weight'] = subgraph[i, j]
        tsp_path_local = nx.approximation.traveling_salesman_problem(G, cycle=True)
        # 将局部编号转换为全局编号
        index_mapping = {i: subset[i] for i in range(len(subset))}
        tsp_path_global = [index_mapping[i] for i in tsp_path_local]
        path_length = sum(subgraph[tsp_path_local[i], tsp_path_local[i+1]] for i in range(len(tsp_path_local)-1))
        return path_length, tsp_path_global


class PatrolPlanner:
    """
    巡逻规划器类，用于计算最优的巡逻方案
    """
    def __init__(self, coords, weights, robot_positions, distance_matrix, path_dict=None, polygons=None):
        """
        初始化巡逻规划器
        
        :param coords: 节点坐标列表 [(x1, y1), (x2, y2), ...]
        :param weights: 节点权重列表 [w1, w2, ...]
        :param robot_positions: 机器人初始位置字典 {robot_id: (x, y)}
        :param distance_matrix: 距离矩阵
        :param path_dict: 路径字典 {(start_coord, end_coord): path}
        :param polygons: 障碍物多边形列表
        """
        self.coords = coords
        self.weights = weights
        self.raw_robot_positions = robot_positions
        self.m = max(robot_positions.keys()) + 1  # 机器人数量
        self.distance_matrix = distance_matrix
        self.path_dict = path_dict
        self.polygons = polygons
        
        # 结果存储
        self.best_partition = None
        self.best_assignment = None
        self.best_time = float('inf')
        self.best_dwell_times = None
        self.best_refresh_times = None
        self.best_phi1 = None
        self.robot_assignment = None
        self.initial_positions = None
        self.tsp_paths = {}
        self.robot_trajectories = None
        
    def generate_valid_assignments(self, m, k):
        """ 生成所有合法的机器人分配方案 (每个子图至少分配 1 个机器人，总数等于 m) """
        base = [1] * k  # 先分配 1 个机器人给每个子图
        remaining = m - k  # 剩余的机器人数量
        for extra in combinations_with_replacement(range(k), remaining):
            assignment = base[:]
            for i in extra:
                assignment[i] += 1
            yield tuple(assignment)
            
    def solve_robot_distribution_brute_force(self, partition, weights, m, distance_matrix):
        """
        通过遍历所有可能的机器人分配方案，找到最小的最大刷新时间。
        """
        k = len(partition)  # 子图数量
        tsp_lengths = {}
        tsp_paths = {}
        
        for j, nodes in partition.items():
            tsp_lengths[j], tsp_paths[j] = GraphCreator.tsp_length_and_path(nodes, distance_matrix)
        
        best_time = float('inf')
        best_assignment = None
        best_partition_refresh_times = None
        best_phi1 = {}

        for assignment in self.generate_valid_assignments(m, k):  # 遍历所有合法分配方案
            max_time = 0
            partition_refresh_times = {}
            for j, robots in enumerate(assignment):
                nodes = partition[j]
                sorted_weights = sorted([weights[i] for i in nodes], reverse=True)  # 按权重降序排序
                
                # 计算 1/w 之和 (考虑前 `robots` 个权重)
                weight_sum_inv = sum(1 / sorted_weights[t] for t in range(min(robots, len(sorted_weights))))
                if weight_sum_inv == 0:
                    continue
                refresh_time = tsp_lengths[j] / weight_sum_inv  # 计算刷新时间
                partition_refresh_times[j] = refresh_time
                
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

        return best_assignment, best_time, best_dwell_times, best_partition_refresh_times, best_phi1, tsp_paths
    
    def find_best_patrolling_plan(self):
        """
        计算全局最优的机器人巡逻方案。
        """
        n = len(self.coords)

        for k in range(2, self.m + 1):  # 2个子图到m个子图的情况
            G = GraphCreator.create_graph(self.coords)
            partition = GraphCreator.maxKcut(G, k)
            assignment, max_refresh_time, dwell_times, refresh_times, phi1_dict, tsp_paths = self.solve_robot_distribution_brute_force(
                partition, self.weights, self.m, self.distance_matrix
            )
            if max_refresh_time < self.best_time:
                self.best_time = max_refresh_time
                self.best_partition = partition
                self.best_assignment = assignment
                self.best_dwell_times = dwell_times
                self.best_refresh_times = refresh_times
                self.best_phi1 = phi1_dict
                self.tsp_paths = tsp_paths
        
        # 单独处理不分割的情况
        partition = {0: list(range(n))}
        assignment, max_refresh_time, dwell_times_nocut, refresh_times_nocut, phi1_dict_nocut, tsp_paths_nocut = self.solve_robot_distribution_brute_force(
            partition, self.weights, self.m, self.distance_matrix
        )
        if max_refresh_time < self.best_time:
            self.best_time = max_refresh_time
            self.best_partition = partition
            self.best_assignment = assignment
            self.best_dwell_times = dwell_times_nocut
            self.best_refresh_times = refresh_times_nocut
            self.best_phi1 = phi1_dict_nocut
            self.tsp_paths = tsp_paths_nocut
        
        return self.best_partition
    
    def compute_subgraph_centroids(self):
        """ 计算每个子图的几何中心 """
        centroids = {}
        for subgraph_id, nodes in self.best_partition.items():
            sub_coords = np.array([self.coords[i] for i in nodes])
            centroids[subgraph_id] = np.mean(sub_coords, axis=0)  # 计算质心
        return centroids
    
    def assign_robots_to_subgraphs(self):
        """
        使用匈牙利算法为每个子图分配最优的机器人
        """
        centroids = self.compute_subgraph_centroids()
        
        # 计算机器人到子图质心的距离
        distances = []
        for robot_id, robot_pos in self.raw_robot_positions.items():
            distances.append([euclidean(robot_pos, centroids[subgraph_id]) for subgraph_id in self.best_partition.keys()])

        # 重复子图节点以满足 n_i 需求（例如 n_i=[2,1] → 子图0出现2次）
        expanded_indices = []
        for i, count in enumerate(self.best_assignment):
            expanded_indices.extend([i] * count)
        
        # 构建代价矩阵（机器人 x 展开后的子图需求）
        cost_matrix = []
        for robot_dists in distances:
            cost_matrix.append([robot_dists[i] for i in expanded_indices])
        
        # 使用匈牙利算法求解最优分配
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        # 生成最终的 robot_assignment
        self.robot_assignment = {subgraph_id: [] for subgraph_id in self.best_partition.keys()}
        for robot_idx, expanded_idx in zip(row_ind, col_ind):
            subgraph_idx = expanded_indices[expanded_idx]
            self.robot_assignment[subgraph_idx].append(robot_idx)
        
        return self.robot_assignment
    
    def find_next_target(self, position, tsp_path, node_coords):
        """
        找到机器人当前位置最近的目标节点
        
        :param position: 机器人当前位置
        :param tsp_path: TSP路径
        :param node_coords: 节点坐标
        :return: (目标节点坐标, 目标节点在tsp_path中的索引)
        """
        min_dist = float('inf')
        target_node = None
        target_idx = None
        
        for i, node_id in enumerate(tsp_path[:-1]):  # 不包括最后一个重复的节点
            node_coord = node_coords[node_id]
            dist = euclidean(position, node_coord)
            if dist < min_dist:
                min_dist = dist
                target_node = node_coord
                target_idx = i
        
        return target_node, target_idx
    
    def shift_to_index(self, arr, index):
        """
        将数组循环移位，使指定索引位置的元素成为第一个元素
        
        :param arr: 输入数组
        :param index: 目标索引
        :return: 移位后的数组
        """
        if not arr:
            return arr
        n = len(arr)
        truncated_arr = arr[:-1]  # 去掉最后一个重复的元素
        shifted_arr = truncated_arr[index:] + truncated_arr[:index]
        shifted = shifted_arr + [shifted_arr[0]]  # 添加循环元素
        return shifted
    
    def calculate_trajectory_for_leader_robot(self, subgraph_id):
        """
        计算子图中领导机器人的完整周期运动轨迹
        
        :param subgraph_id: 子图ID
        :return: (轨迹, 领导机器人ID)
        """
        tsp_path = self.tsp_paths[subgraph_id]
        robot_list = self.robot_assignment[subgraph_id]
        
        trajectory = []
        current_time = 0
        n = len(tsp_path)

        # 计算列表中机器人的初始位置与tsp_path中包含节点的位置的距离矩阵
        distance_matrix_ = np.array([[euclidean(self.raw_robot_positions[i], self.coords[j]) 
                                     for j in tsp_path[:-1]] for i in robot_list])
        
        # 找到距离矩阵中最小值的索引
        min_index = np.unravel_index(np.argmin(distance_matrix_), distance_matrix_.shape)
        # 获取最小值的行和列索引
        min_row_index, min_col_index = min_index
        # 获取最小值对应的机器人编号
        leader_index = robot_list[min_row_index]
        # 将tsp_path进行旋转，以最近节点为起点
        tsp_path = self.shift_to_index(tsp_path, min_col_index)
        
        trajectory.append((0, self.coords[tsp_path[0]]))
        for i in range(n-1):
            node_current = tsp_path[i]
            coord_current = self.coords[node_current]        
            current_time += self.best_dwell_times[node_current]  # 停留时间
            trajectory.append((current_time, coord_current))  # 停留后的时间戳

            node_next = tsp_path[(i + 1) % n]  # 下一个节点（循环）
            coord_next = self.coords[node_next]
            
            if self.path_dict:
                path_to_next = self.path_dict.get((coord_current, coord_next), 
                                                 [(coord_current, coord_next)])
                for j in range(len(path_to_next) - 1):
                    # 计算路径上的每个点
                    path_length = euclidean(path_to_next[j], path_to_next[j + 1])
                    travel_time = path_length
                    current_time += travel_time
                    trajectory.append((current_time, path_to_next[j + 1]))
            else:
                # 如果没有路径字典，直接计算欧氏距离
                distance = euclidean(coord_current, coord_next)
                current_time += distance  # 假设速度为1
                trajectory.append((current_time, coord_next))
        
        return trajectory, leader_index, tsp_path
    
    def calculate_trajectory_for_follower_robot(self, subgraph_id, tsp_path, initial_position, robot_id):
        """
        计算跟随机器人的完整周期运动轨迹
        
        :param subgraph_id: 子图ID
        :param tsp_path: TSP路径
        :param initial_position: 机器人初始位置
        :param robot_id: 机器人ID
        :return: 轨迹
        """
        trajectory = []
        current_time = 0.0
        n = len(tsp_path)
        
        # 寻找tsp_path中下一个目标位置
        next_target, target_node_id = self.find_next_target(initial_position, tsp_path, self.coords)
        
        # 计算初始位置到下一个目标位置的距离
        initial_distance = euclidean(initial_position, next_target)
        trajectory.append((current_time, initial_position))  # 初始位置的时间戳和位置
        current_time += initial_distance  # 移动到下一个目标位置的时间（假设速度为1）
        
        shifted_tsp_path = self.shift_to_index(tsp_path, target_node_id)
        
        for i in range(n-1):
            node_current = shifted_tsp_path[i]
            coord_current = self.coords[node_current]
            trajectory.append((current_time, coord_current))  # 当前节点的时间戳和位置
            current_time += self.best_dwell_times[node_current]  # 停留时间
            trajectory.append((current_time, coord_current))  # 停留后的时间戳
            
            # 移动到下一个节点
            node_next = shifted_tsp_path[(i + 1) % n]  # 下一个节点（循环）
            coord_next = self.coords[node_next]

            if self.path_dict:
                path_to_next = self.path_dict.get((coord_current, coord_next), 
                                                 [(coord_current, coord_next)])
                for j in range(len(path_to_next) - 1):
                    # 计算路径上的每个点
                    path_length = euclidean(path_to_next[j], path_to_next[j + 1])
                    travel_time = path_length
                    current_time += travel_time
                    trajectory.append((current_time, path_to_next[j + 1]))
            else:
                # 如果没有路径字典，直接计算欧氏距离
                distance = euclidean(coord_current, coord_next)
                current_time += distance  # 假设速度为1
                trajectory.append((current_time, coord_next))
        
        return trajectory
    
    def interpolate_position(self, trajectory, time):
        """
        通过插值计算在给定时间点的位置
        
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
    
    def get_robot_initial_position(self, trajectory, robots_count, T, phi_1, robot_idx):
        """
        计算机器人在周期T内的初始位置
        
        :param trajectory: 领导机器人的完整轨迹
        :param robots_count: 子图中的机器人数量
        :param T: 子图的刷新时间
        :param phi_1: 子图中的第m+1大权重
        :param robot_idx: 机器人在子图中的索引
        :return: 机器人在周期T内的初始位置
        """
        # 计算时间点
        time_i = (robots_count - robot_idx) * T / phi_1
        
        # 查找时间点位于哪一段
        return self.interpolate_position(trajectory, time_i)
    
    def calculate_all_trajectories(self):
        """
        计算所有机器人的轨迹
        
        :return: 机器人轨迹字典 {robot_id: trajectory}
        """
        self.robot_trajectories = {}
        
        # 为每个子图计算轨迹
        for subgraph_id, robots in self.robot_assignment.items():
            if not robots:
                continue
                
            # 计算领导机器人的轨迹
            leader_trajectory, leader_id, tsp_path = self.calculate_trajectory_for_leader_robot(subgraph_id)
            self.robot_trajectories[leader_id] = leader_trajectory
            
            # 计算跟随机器人的轨迹
            robots_count = len(robots)
            phi_1 = self.best_phi1.get(subgraph_id, 1.0)  # 默认为1.0
            refresh_time = self.best_refresh_times.get(subgraph_id, 0.0)
            
            for i, robot_id in enumerate(robots):
                if robot_id == leader_id:
                    continue  # 跳过领导机器人
                
                # 计算初始位置
                initial_position = self.get_robot_initial_position(
                    leader_trajectory, robots_count, refresh_time, phi_1, i
                )
                
                # 计算轨迹
                follower_trajectory = self.calculate_trajectory_for_follower_robot(
                    subgraph_id, tsp_path, initial_position, robot_id
                )
                
                self.robot_trajectories[robot_id] = follower_trajectory
        
        return self.robot_trajectories
    
    def plan(self):
        """
        执行完整的巡逻规划过程
        
        :return: (best_partition, robot_assignment, robot_trajectories)
        """
        # 1. 找到最佳分区和机器人分配
        self.find_best_patrolling_plan()
        
        # 2. 将机器人分配到子图
        self.assign_robots_to_subgraphs()
        
        # 3. 计算所有机器人的轨迹
        self.calculate_all_trajectories()
        
        return self.best_partition, self.robot_assignment, self.robot_trajectories
    
    def visualize_partition(self, ax=None):
        """
        可视化子图分区
        
        :param ax: matplotlib轴对象，如果为None则创建新的
        :return: matplotlib轴对象
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
            
        colors = plt.cm.tab10(np.linspace(0, 1, len(self.best_partition)))
        
        # 绘制节点和分区
        for subgraph_id, nodes in self.best_partition.items():
            color = colors[subgraph_id]
            
            # 绘制TSP路径连线
            tsp_path = self.tsp_paths.get(subgraph_id, [])
            if tsp_path:
                for i in range(len(tsp_path)-1):
                    node_i, node_j = tsp_path[i], tsp_path[i+1]
                    coord_i, coord_j = self.coords[node_i], self.coords[node_j]
                    
                    # 如果有路径字典，使用实际路径
                    if self.path_dict and (coord_i, coord_j) in self.path_dict:
                        path = self.path_dict[(coord_i, coord_j)]
                        for p in range(len(path)-1):
                            ax.plot([path[p][0], path[p+1][0]], 
                                   [path[p][1], path[p+1][1]],
                                   color=color, alpha=0.5, linestyle='-', linewidth=1.5)
                    else:
                        # 直接连线
                        ax.plot([coord_i[0], coord_j[0]], 
                               [coord_i[1], coord_j[1]],
                               color=color, alpha=0.5, linestyle='-', linewidth=1.5)
            
            # 绘制节点
            for node_id in nodes:
                x, y = self.coords[node_id]
                ax.scatter(x, y, color=color, s=100, alpha=0.7)
                ax.text(x, y, str(node_id), fontsize=12)
        
        # 绘制障碍物
        if self.polygons:
            for poly in self.polygons:
                x, y = poly.exterior.xy
                ax.fill(x, y, alpha=0.3, fc='gray', ec='black')
        
        # 绘制机器人初始位置
        for robot_id, pos in self.raw_robot_positions.items():
            ax.scatter(pos[0], pos[1], marker='^', color='red', s=150)
            ax.text(pos[0], pos[1], f'R{robot_id}', fontsize=12)
            
        ax.set_title('子图分区和机器人初始位置')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True)
        
        return ax
    
    def visualize_trajectories(self, duration=10, interval=200, save_path=None, frames=200):
        """
        可视化机器人轨迹动画
        
        :param duration: 动画持续时间（秒）
        :param interval: 帧间隔（毫秒），值越大动画越慢
        :param save_path: 保存路径，如果为None则不保存
        :param frames: 动画帧数，值越大动画越平滑
        :return: 动画对象
        """
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # 绘制节点
        for i, (x, y) in enumerate(self.coords):
            ax.scatter(x, y, color='blue', s=100, alpha=0.5)
            ax.text(x, y, str(i), fontsize=12)
        
        # 绘制TSP路径连线，与机器人实际运动路径一致
        for subgraph_id, nodes in self.best_partition.items():
            color = plt.cm.tab10(subgraph_id / len(self.best_partition))
            tsp_path = self.tsp_paths.get(subgraph_id, [])
            if tsp_path:
                for i in range(len(tsp_path)-1):
                    node_i, node_j = tsp_path[i], tsp_path[i+1]
                    coord_i, coord_j = self.coords[node_i], self.coords[node_j]
                    
                    # 如果有路径字典，使用实际路径
                    if self.path_dict and (coord_i, coord_j) in self.path_dict:
                        path = self.path_dict[(coord_i, coord_j)]
                        for p in range(len(path)-1):
                            ax.plot([path[p][0], path[p+1][0]], 
                                   [path[p][1], path[p+1][1]],
                                   color=color, alpha=0.3, linestyle='--', linewidth=1)
                    else:
                        # 直接连线
                        ax.plot([coord_i[0], coord_j[0]], 
                               [coord_i[1], coord_j[1]],
                               color=color, alpha=0.3, linestyle='--', linewidth=1)
        
        # 绘制障碍物
        if self.polygons:
            for poly in self.polygons:
                x, y = poly.exterior.xy
                ax.fill(x, y, alpha=0.3, fc='gray', ec='black')
        
        # 为每个机器人创建散点对象和轨迹线
        robot_dots = {}
        robot_trails = {}
        colors = plt.cm.rainbow(np.linspace(0, 1, len(self.robot_trajectories)))
        
        for i, robot_id in enumerate(self.robot_trajectories.keys()):
            color = colors[i]
            dot, = ax.plot([], [], 'o', markersize=10, color=color, label=f'Robot {robot_id}')
            trail, = ax.plot([], [], '-', linewidth=1, color=color, alpha=0.5)
            robot_dots[robot_id] = dot
            robot_trails[robot_id] = trail
        
        # 计算所有轨迹的最大时间
        max_times = [traj[-1][0] for traj in self.robot_trajectories.values()]
        max_time = max(max_times) if max_times else 1.0
        
        def init():
            for dot in robot_dots.values():
                dot.set_data([], [])
            for trail in robot_trails.values():
                trail.set_data([], [])
            return list(robot_dots.values()) + list(robot_trails.values())
        
        # 存储轨迹历史
        trail_history = {robot_id: ([], []) for robot_id in self.robot_trajectories.keys()}
        
        def update(frame):
            # 计算当前时间
            current_time = (frame / frames) * max_time
            
            # 更新每个机器人的位置和轨迹
            for robot_id, dot in robot_dots.items():
                trajectory = self.robot_trajectories[robot_id]
                position = self.interpolate_position(trajectory, current_time)
                
                # 更新当前位置
                dot.set_data([position[0]], [position[1]])
                
                # 更新轨迹历史
                x_history, y_history = trail_history[robot_id]
                x_history.append(position[0])
                y_history.append(position[1])
                
                # 只保留最近的100个点，避免轨迹过长
                if len(x_history) > 100:
                    x_history = x_history[-100:]
                    y_history = y_history[-100:]
                
                trail_history[robot_id] = (x_history, y_history)
                robot_trails[robot_id].set_data(x_history, y_history)
            
            return list(robot_dots.values()) + list(robot_trails.values())
        
        ani = FuncAnimation(fig, update, frames=frames, init_func=init, 
                           interval=interval, blit=True)
        
        ax.legend()
        ax.set_title('机器人巡逻轨迹')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True)
        
        # 保存动画
        if save_path:
            writer = PillowWriter(fps=15)  # 降低帧率使动画更容易观察
            ani.save(save_path, writer=writer)
        
        return ani

