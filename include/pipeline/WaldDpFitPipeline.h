#ifndef WALD_DP_FIT_PIPELINE_H
#define WALD_DP_FIT_PIPELINE_H

#include "core/WaldDpFit.h"

#include <string>

// 从 inputPath 的 inputSheet 中读取 x/y 两列，执行：
//   Wald 检验(全部系数) -> DP 最优分段 -> KKT 硬约束联合拟合。
// 打印结果，并把 DP 得到的分段结果追加到“分段拟合”表。
WaldDpFitResult runWaldDpFitPipeline(
    const std::string& inputPath,
    const std::string& inputSheet,
    const std::string& xHeader,
    const std::string& yHeader,
    double dt,
    int windowHalfWidth = 7,
    int minSegmentLength = 6,
    int maxSegmentLength = 0,
    double alpha = 0.05,
    double muP = 1.0e4,
    double muV = 1.0e4);

#endif
