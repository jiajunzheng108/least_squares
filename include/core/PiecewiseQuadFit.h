#ifndef PIECEWISE_QUAD_FIT_H
#define PIECEWISE_QUAD_FIT_H

#include <vector>

// 单个分段的拟合结果。
// 段内局部时间 u 统一归一化到 [-1, 1]，坐标为 q(u)=a*u^2+b*u+c。
struct PiecewiseSegmentResult {
    int index;        // 段号(1-based)
    int startPoint;   // 起始点号(1-based)
    int pointCount;   // 点数(普通段 10，最后一段 12)
    double halfSpan;  // 局部半宽 h(秒)
    double coeffX[3]; // X 的 [a, b, c]
    double coeffY[3]; // Y 的 [a, b, c]
    double rmseX;     // 段内 X 的 RMSE(含重合点)
    double rmseY;     // 段内 Y 的 RMSE(含重合点)
};

// 节点(相邻段重合点)的拟合结果。
// 位置/速度/加速度取自左侧段右端；da 为右侧段加速度减左侧段加速度。
struct PiecewiseNodeResult {
    int index;     // 节点号(1-based)
    int point;     // 节点所在点号(1-based)
    double t;      // 时间(秒)
    double x, y;   // 位置
    double vx, vy; // 速度
    double ax, ay; // 左侧段加速度
    double dax;    // 加速度跳变 X 分量
    double day;    // 加速度跳变 Y 分量
    double ep;     // 位置跳变模长
    double ev;     // 速度跳变模长
    double ea;     // 加速度跳变模长
};

// 分段二次拟合的完整结果。
struct PiecewiseQuadFitResult {
    int pointCount; // 原始点数
    double dt;      // 时间步长(秒)
    int window;     // 普通段点数
    int step;       // 相邻段起点间隔
    int lastWindow; // 最后一段点数

    std::vector<PiecewiseSegmentResult> segments;
    std::vector<PiecewiseNodeResult> nodes;

    std::vector<double> xFit;       // 各点 X 拟合值(每个点只填一次)
    std::vector<double> yFit;       // 各点 Y 拟合值
    std::vector<double> xResidual;  // X 残差 = 观测 - 拟合
    std::vector<double> yResidual;  // Y 残差
    std::vector<int> pointSegment;  // 各点所属段号(1-based)
    std::vector<int> pointType;     // 0=普通点，1=节点

    double sseX;             // X 残差平方和(每个点只计一次)
    double sseY;             // Y 残差平方和
    double rmseX;            // sqrt(sseX / pointCount)
    double rmseY;            // sqrt(sseY / pointCount)
    double maxAbsResidualX;  // X 最大绝对残差
    double maxAbsResidualY;  // Y 最大绝对残差
};

// 分段二次最小二乘拟合。
// 前 K-1 段取 window 个点、步长 step；最后一段从同一节点开始，取 lastWindow 个点。
// 相邻段在重合点强制位置连续和速度连续(C1)，不强制加速度连续。
// 默认 window=10, step=9, lastWindow=12，对应 111 点：前 11 段各 10 点，
// 第 12 段覆盖第 100..111 点(12 点)。
PiecewiseQuadFitResult fitPiecewiseQuadratic(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double dt,
    int window = 10,
    int step = 9,
    int lastWindow = 12);

#endif
