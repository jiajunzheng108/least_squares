#ifndef WALD_DP_FIT_H
#define WALD_DP_FIT_H

#include <vector>

// Wald 候选点筛选 + DP 最优分段 + KKT 硬约束联合拟合的参数。
struct WaldDpConfig {
    double dt = 0.001;          // 时间步长(秒)
    int windowHalfWidth = 7;    // Wald 窗口半宽 m
    int minSegmentLength = 6;   // 最小段长 L_min(点数)
    int maxSegmentLength = 0;   // 最大段长 L_max(点数)，0 表示不限制
    double alpha = 0.05;        // 全局显著性水平
    double muP = 1.0e4;         // DP 中位置连续项权重(归一化后)
    double muV = 1.0e4;         // DP 中速度连续项权重(归一化后)
};

struct WaldCandidateInfo {
    int point;    // 候选点号(1-based)
    double wald;  // Wald 统计量 W_i
    double pValue;
};

struct WaldSegmentInfo {
    int startPoint;  // 1-based
    int endPoint;    // 1-based
    int pointCount;
    double halfSpan; // 秒
    double coeffX[3];  // 局部 u 下的系数
    double coeffY[3];
    double sseX;       // 段内未加权 SSE(含重合点)
    double sseY;
};

struct WaldNodeInfo {
    int point;     // 节点所在点号(1-based)
    double t;      // 秒
    double x, y;
    double vx, vy;
    double ax, ay; // 左侧段加速度
    double dax, day; // 加速度跳变(右-左)
    double ep, ev, ea;
};

struct WaldKPathResult {
    int K = 0;
    std::vector<int> nodes; // 1-based，含首尾
    std::vector<WaldSegmentInfo> segments;
    std::vector<WaldNodeInfo> nodeInfo;
    double sseX = 0.0;
    double sseY = 0.0;
    double sseTotal = 0.0;
    double rmseX = 0.0;
    double rmseY = 0.0;
    double rmseTotal = 0.0;
    double epMax = 0.0;
    double evMax = 0.0;
    double eaMax = 0.0;
    double dpCost = 0.0;
};

struct WaldDpFitResult {
    int pointCount = 0;
    WaldDpConfig config;

    int rawCandidateCount = 0; // NMS 前的候选数
    std::vector<WaldCandidateInfo> candidates; // NMS 后，不含首尾
    double waldCritical = 0.0;
    double positionScale = 0.0;
    double velocityScale = 0.0;

    // 下标为 K；某些 K 无合法路径时为空。
    std::vector<WaldKPathResult> dpPaths;
    std::vector<WaldKPathResult> bestPaths;

    int chosenK = 0;
    WaldKPathResult chosen;
    bool enumerationTruncated = false;
};

// 使用全部系数(A,B,C)的 Wald 检验筛选候选点，然后 DP + KKT 拟合。
WaldDpFitResult fitWaldDp(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const WaldDpConfig& config);

#endif
