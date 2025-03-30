from ortools.graph.python import min_cost_flow

def allocate_robots(distances, n_i):
    smcf = min_cost_flow.SimpleMinCostFlow()
    m = len(distances)
    k = len(n_i)
    
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
            smcf.add_arc_with_capacity_and_unit_cost(j, subgraph_nodes[i], 1, distances[idx][i])
    
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
    allocation = {}
    for arc in range(smcf.num_arcs()):
        if smcf.tail(arc) in robot_nodes and smcf.flow(arc) > 0:
            robot = smcf.tail(arc) - 1
            subgraph = smcf.head(arc) - (m+1)
            allocation[robot] = subgraph
    return allocation

# 示例用法
distances = [
    [1, 3],  # 机器人0到子图0和1的距离
    [2, 2],  # 机器人1到子图0和1的距离
]
n_i = [1, 1]  # 每个子图需要1个机器人
print(allocate_robots(distances, n_i))  # 输出：{0: 0, 1: 1}