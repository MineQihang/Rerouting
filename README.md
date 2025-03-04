# Rerouting in an optical network
[华为算法精英实战营第十五期-Rerouting in an optical network](https://competition.huaweicloud.com/information/1000042138/introduction)

赛题简介：一个光网络可以表示为图的形式。图中的每条边是光缆（光纤）。信息通过特定波长（频率）的光进行传输。光纤可以同时传输多个不同波长的信号。随着时间的推移，网络中的某些元件（光纤）可能会发生故障，我们需要为中断的服务（连接）找到新的路径。此外，网络中可能会发生随机的故障序列。任务的目标是尽可能提供更多服务的正常运行。

最终名次：第一名

## 代码说明

主要代码逻辑见`src/rerouting.cpp`。

## 解题思路

对未分配的业务依次用bfs找最短路径（不同频率都找），然后对找到的路径按照一个权值排序，选择权值最大的路径作为该业务的候选路。

权值初始为0，首先考虑这个路径与未分配业务的可以找到路径之间的重合边，多一条重合边权值减1；然后考虑这个路径与分配好的业务在不使用当前路径后再找路径的重合边，多一条重合边权值减1。（这个权值就是在考虑这条路径对当前其他未分配业务找路径和未来断边后找路径的影响）

- 优化点1: 分配好的业务会很多，所以只需要挑出路径长度最大的$k$个即可。（$k$需要调参）

- 优化点2: 考虑未分配业务的顺序很重要，因此可以迭代t次，每次随机shuffle一次顺序。（$t$根据运行时间调整）

- 优化点3：为了能剔除肯定找不到路径的业务，从而加快求解速度，我使用了Floyd-Warshall算法做预处理。

详细解题思路见[文档](./docs/华为算法精英实战营第十五期作品思路-mine_qihang.pdf)。
