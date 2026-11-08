#ifndef PIECEWISE_QUAD_FIT_PIPELINE_H
#define PIECEWISE_QUAD_FIT_PIPELINE_H

#include "core/PiecewiseQuadFit.h"

#include <string>

// 从 inputPath 的 inputSheet 中按表头名读取 x/y 两列，做分段二次拟合
// (相邻段 C1 连续)，并把分段系数、整体精度和节点连续性打印到控制台。
// 不写回 xlsx。默认 10 点/段、步长 9、最后一段 12 点。
PiecewiseQuadFitResult runPiecewiseQuadFitPipeline(
    const std::string& inputPath,
    const std::string& inputSheet,
    const std::string& xHeader,
    const std::string& yHeader,
    double dt,
    int window = 10,
    int step = 9,
    int lastWindow = 12);

#endif
