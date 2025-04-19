import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from patroll_new import PatrolPlanner
from distance_caculation import DistanceCalculator

def main():
    # 测试数据
    coords = [(4, 1), (2, 3), (5, 5), (8, 8), (12, 3), (6, 9), (14, 10), (7, 2), (10, 5)]
    weights = [1.5, 2.0, 1.2, 1.9, 2.5, 1.8, 2.2, 3.0, 6.0]
    raw_robot_positions = {0: (8, 4), 1: (7, 5), 2: (5, 6)}
    
    # 障碍物
    polygons = [
        Polygon([(12, 6), (12, 8), (14, 8), (14, 6)]),
        Polygon([(3, 4), (3, 5), (4, 5), (4, 4)]),
        Polygon([(8, 6), (8, 7), (10, 7), (10, 6)]),
    ]
    
    # 距离计算
    distance_calculator = DistanceCalculator(method="visibility", polygons=polygons)
    path_dict, dist_graph, global_distance_matrix = distance_calculator.compute_all_pairs(coords)
    
    # 创建巡逻规划器
    planner = PatrolPlanner(
        coords=coords,
        weights=weights,
        robot_positions=raw_robot_positions,
        distance_matrix=global_distance_matrix,
        path_dict=path_dict,
        polygons=polygons
    )
    
    # 执行规划
    best_partition, robot_assignment, robot_trajectories = planner.plan()
    
    # 打印结果
    print("最优分割方案:", best_partition)
    print("\n最优机器人分配方案 (子图编号: 机器人列表):")
    for subgraph_id, robots in robot_assignment.items():
        print(f" - 子图 {subgraph_id}: {robots} 个机器人")
    
    print(f"\n最小的最大刷新时间: {planner.best_time}")
    
    print("\n节点停留时间:")
    for i, time in planner.best_dwell_times.items():
        if time > 0:
            print(f" - 节点 {i}: {time:.2f}")
    
    # 可视化分区
    plt.figure(figsize=(12, 10))
    ax = planner.visualize_partition()
    plt.savefig('partition_visualization.png')
    plt.close()
    
    # 可视化轨迹动画
    ani = planner.visualize_trajectories(interval=100, save_path='robot_trajectories.gif')
    plt.show()
    
    # 打印每个机器人的轨迹信息
    print("\n机器人轨迹信息:")
    for robot_id, trajectory in robot_trajectories.items():
        print(f" - 机器人 {robot_id} 的轨迹长度: {len(trajectory)} 个点, 周期时间: {trajectory[-1][0]:.2f}")

if __name__ == "__main__":
    main()