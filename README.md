# 1、PiecewiseQuadFit（恒定步长+KKT调整）

步长10个点，两段之间重合一点，每段先用分方向二次多项式拟合，再加上约束条件用KKT求出最终结果

约束条件：段之间位置和速度不变


$$
J(\theta) = \frac{1}{2} \left\| W^{1/2} (A\theta - b) \right\|^2
$$


W：权重矩阵

A：设计矩阵

b：观测值


$$
\begin{gather*}
\min_{\theta} J(\theta) \quad \text{s.t.} \quad C\theta = 0 \\
\\
L(\theta, \lambda) = \frac{1}{2} \theta^T H \theta - g^T \theta + \lambda^T C \theta
\end{gather*}
$$


C：节点连续性约束

H = A^TWA

g = A*TWb

# 2、WaldDpFit（Wald检验+动态规划+KKT调整）

1. wald检验：由位置、速度和加速度驱动，选择左右局部轨迹差异较大的点作为切分点的候选点

   

2. 动态规划：

   ​	
$$
   	dp[k][i][j] = \min_{h \lt i} \left\{ dp[k - 1][h][i] + \mathrm{SSR}_{i,j} + \mu_p \left( \frac{\left\| p_{h,i}^R - p_{i,j}^L \right\|}{s_p} \right)^2 + \mu_v \left( \frac{\left\| v_{h,i}^R - v_{i,j}^L \right\|}{s_v} \right)^2 \right\}
$$
   

   dp[k-1]\[h][i]：前(k-1)段的最小代价

   SSR<sub>i,j</sub>：当前段独立拟合的残差平方和

   位置惩罚：上一段右端与当前段左端的位置差

   速度惩罚：上一段右端与当前段左端的速度差

3. 5%边际规则选分割点
