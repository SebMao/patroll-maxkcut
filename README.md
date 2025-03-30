## 相关库的安装
```
conda create -n patroll python==3.8
conda activate patroll
pip3 install -r requirements.txt
```

## 参考资料
最大K割问题：https://github.com/QCOL-LU/MaxKcut

机器人巡逻：
Cooperative Patrolling via Weighted Tours: Performance Analysis and Distributed Algorithms, TRO, 2012、


## 协议
- 输入：节点经纬度、艇经纬度、节点权重（间隔要求）、节点上次访问时间、任务区域、禁止区域、（艇速）
- 输出：艇的路线、节点停留时间


## 算法设计
给定n个节点，m个机器人
1. 对n个节点进行最大k割划分，共有m种划分方法(k=1,...,m)。
2. 对每种划分，求解最小整体刷新时间的分配方案，找到最优的划分。
3. 按照分配方案，对机器人进行指派。
4. 根据指派方案，求解初始阶段机器人的移动方案。


