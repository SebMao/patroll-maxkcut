import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def interpolate_position(traj, t):
    if not traj:
        return (0, 0)  # 空轨迹兜底

    # 获取整个周期总时长（最后一个点的时间）
    total_time = traj[-1][0]
    
    # 周期性处理：将 t 映射回 [0, total_time)
    t = t % total_time

    # 查找 t 所在的时间段
    for i in range(len(traj) - 1):
        t0, p0 = traj[i]
        t1, p1 = traj[i + 1]
        if t0 <= t <= t1 and t1 > t0:
            ratio = (t - t0) / (t1 - t0)
            x = p0[0] + ratio * (p1[0] - p0[0])
            y = p0[1] + ratio * (p1[1] - p0[1])
            return (x, y)

    # 如果 t 恰好在最后一个点（理论上不会发生），返回最后点
    return traj[-1][1] if traj else (0, 0)


# 示例轨迹（你实际中替换为 robot_trajectories 字典）
robot_trajectories = {
    0: [(0.0, (10, 5)), (0.0, (10, 5)), (12.3581, (10, 5)), (18.015, (6, 9)), 
        (18.015, (6, 9)), (20.251, (8, 8)), (20.251, (8, 8)), (26.5756, (14, 10)),
        (26.5756, (14, 10)), (33.8557, (12, 3)), (36.1973, (12, 3)), (41.8541, (10, 5))],
    
    1: [(0.0, (7.34, 8.33)), (0.7382, (8, 8)), (0.7382, (8, 8)), (7.0628, (14, 10)),
        (7.0628, (14, 10)), (14.3429, (12, 3)), (16.6844, (12, 3)), (19.5128, (10, 5)), 
        (31.871, (10, 5)), (37.5278, (6, 9)), (37.5278, (6, 9)), (41.2618, (7.34, 8.33))]
}

# 统一的采样时间间隔（秒）
time_interval = 0.5

# 获取整个周期的最大时间
max_time = max(max([pt[0] for pt in traj]) for traj in robot_trajectories.values())
times = np.arange(0, max_time, time_interval)

# 所有时间点的快照 [{id: (x,y)}]
snapshots = []
for t in times:
    snapshot = {}
    for robot_id, traj in robot_trajectories.items():
        pos = interpolate_position(traj, t)
        if not isinstance(pos, tuple):
            print(f"Warning: robot {robot_id} at time {t} -> bad position {pos}")
        snapshot[robot_id] = pos
    snapshots.append(snapshot)

def update(frame):
    snap = snapshots[frame]
    for i, robot_id in enumerate(robot_trajectories.keys()):
        pos = snap.get(robot_id, (0, 0))  # 避免 KeyError
        scatters[i].set_data([pos[0]], [pos[1]])  # 注意：加中括号！
    return scatters


# 可视化+动画
fig, ax = plt.subplots()
colors = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
scatters = [ax.plot([], [], 'o', color=colors[i % len(colors)], label=f'Robot {i}')[0] 
            for i in range(len(robot_trajectories))]
ax.set_xlim(0, 20)
ax.set_ylim(0, 20)
ax.legend()



ani = animation.FuncAnimation(fig, update, frames=len(snapshots), interval=200)

# 保存为GIF
ani.save("robot_trajectory.gif", writer='pillow')
