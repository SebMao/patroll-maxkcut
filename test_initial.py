import numpy as np
from scipy.spatial.distance import pdist, squareform, euclidean

def calculate_trajectory_for_robot_0(tsp_path, dwell_times, distance_matrix, coords, speed=1):
    """
    计算0号机器人的完整周期运动轨迹，并返回分段点。
    
    :param tsp_path: 0号机器人TSP路径，包含访问节点的顺序
    :param dwell_times: 每个节点的停留时间
    :param distance_matrix: 距离矩阵
    :param coords: 节点坐标 [(x1, y1), (x2, y2), ...]
    :param speed: 机器人匀速移动的速度，默认为1
    :return: 分段点列表 [(time, position)]，其中position为(x, y)坐标
    """
    trajectory = []
    current_time = 0
    n = len(tsp_path)

    
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

# 下面是一个简单的例子如何使用这些函数
# 假设tsp_path为0号机器人TSP路径，dwell_times为停留时间，distance_matrix为距离矩阵

tsp_path = [0, 1, 2, 3, 0]  # 示例TSP路径
dwell_times = [0, 10.0, 15.0, 5.0]  # 示例停留时间
# distance_matrix = np.array([[0, 10, 15, 20], [10, 0, 10, 15], [15, 10, 0, 10], [20, 15, 10, 0]])  # 示例距离矩阵
coords = [(4, 1), (4, 5), (5, 5), (1, 5)]  # 示例节点坐标
n = 4
distance_matrix = np.array([[euclidean(coords[i], coords[j]) for j in range(n)] for i in range(n)])


# 计算0号机器人的轨迹
trajectory_robot_0 = calculate_trajectory_for_robot_0(tsp_path, dwell_times, distance_matrix, coords)
# 打印轨迹信息
print("0号机器人完整周期的轨迹:")
for time, position in trajectory_robot_0:
    print(f"时间: {time:.2f}, 位置: {position}")
# 获取其他机器人的初始位置
T = 41  # 假设刷新时间
phi_1 = 2.5  # 假设第m+1大权重
m = 3  # 假设有3个机器人
for i in range(1,m):
    initial_position = get_robot_initial_position(trajectory_robot_0, m, T, phi_1, i)
    print(f"机器人{i}的初始位置: {initial_position}")
