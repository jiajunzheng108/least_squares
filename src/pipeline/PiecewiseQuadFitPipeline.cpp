// 分段二次拟合的操作封装函数
#include "pipeline/PiecewiseQuadFitPipeline.h"

#include "core/TimeQuadFit.h"
#include "utils/XlsxReader.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double roundTo6(double value) {
    if (std::isnan(value))
        return value;
    return std::round(value * 1.0e6) / 1.0e6;
}

void printSegmentHeader() {
    std::cout << "  seg  points      n   h(ms)"
                 "          aX           bX           cX"
                 "          aY           bY           cY"
                 "      RMSE_X      RMSE_Y\n";
}

}  // namespace

PiecewiseQuadFitResult runPiecewiseQuadFitPipeline(
    const std::string& inputPath, const std::string& inputSheet,
    const std::string& xHeader, const std::string& yHeader, double dt,
    int window, int step, int lastWindow) {
    std::vector<double> x;
    std::vector<double> y;
    readPointsFromXlsx(inputPath, inputSheet, xHeader, yHeader, x, y);

    const PiecewiseQuadFitResult result =
        fitPiecewiseQuadratic(x, y, dt, window, step, lastWindow);

    std::cout << "=== PiecewiseQuadFit: " << inputPath << " / " << inputSheet
              << " (" << xHeader << ", " << yHeader << ") ===\n";
    std::cout << "points = " << result.pointCount << "  dt = " << result.dt
              << " s  window = " << result.window
              << "  step = " << result.step
              << "  lastWindow = " << result.lastWindow << "\n";
    std::cout << "segments = " << result.segments.size()
              << "  nodes = " << result.nodes.size() << "\n\n";

    std::cout << "[segments] local u in [-1,1], q(u)=a*u^2+b*u+c, "
                 "h = half span\n";
    printSegmentHeader();
    std::cout << std::fixed << std::setprecision(6);
    for (const PiecewiseSegmentResult& segment : result.segments) {
        const int endPoint =
            segment.startPoint + segment.pointCount - 1;
        std::cout << std::setw(5) << segment.index << "  " << std::setw(5)
                  << segment.startPoint << "-" << std::setw(4) << endPoint
                  << "  " << std::setw(5) << segment.pointCount << "  "
                  << std::setw(6) << std::setprecision(2)
                  << segment.halfSpan * 1000.0 << "  " << std::setprecision(6)
                  << std::setw(12) << segment.coeffX[0] << "  "
                  << std::setw(12) << segment.coeffX[1] << "  "
                  << std::setw(12) << segment.coeffX[2] << "  "
                  << std::setw(12) << segment.coeffY[0] << "  "
                  << std::setw(12) << segment.coeffY[1] << "  "
                  << std::setw(12) << segment.coeffY[2] << "  "
                  << std::setw(10) << segment.rmseX << "  "
                  << std::setw(10) << segment.rmseY << "\n";
    }

    std::cout << "\n[overall]\n";
    std::cout << std::setprecision(6);
    std::cout << "  X: SSE=" << result.sseX << "  RMSE=" << result.rmseX
              << "  max|res|=" << result.maxAbsResidualX << "\n";
    std::cout << "  Y: SSE=" << result.sseY << "  RMSE=" << result.rmseY
              << "  max|res|=" << result.maxAbsResidualY << "\n";

    std::cout << "\n[nodes] position/velocity/acceleration from left segment; "
                 "da = right - left\n";
    std::cout << " node point   t(ms)            X            Y"
                 "           vX           vY           aX           aY"
                 "          daX          daY          e_a          e_p"
                 "          e_v\n";
    for (const PiecewiseNodeResult& node : result.nodes) {
        std::cout << std::setw(5) << node.index << "  " << std::setw(5)
                  << node.point << "  " << std::setw(6)
                  << std::setprecision(3) << node.t * 1000.0 << "  "
                  << std::setprecision(6) << std::setw(12) << node.x << "  "
                  << std::setw(12) << node.y << "  " << std::setw(12)
                  << node.vx << "  " << std::setw(12) << node.vy << "  "
                  << std::setw(12) << node.ax << "  " << std::setw(12)
                  << node.ay << "  " << std::setw(12) << node.dax << "  "
                  << std::setw(12) << node.day << "  " << std::setw(12)
                  << node.ea << "  " << std::setw(12) << node.ep << "  "
                  << std::setw(12) << node.ev << "\n";
    }

    // 写入“分段拟合”表：每段一行 + 最后一行总体 SSE。
    // 系数使用段内局部时间 t(从每段第一个点 t=0 开始，单位秒)：
    // q(t) = a*t^2 + b*t + c。
    const std::string outputSheet = "分段拟合";
    const std::size_t rowCount = result.segments.size() + 1;
    const std::size_t columnCount = 10;
    std::vector<std::vector<double>> columns(
        columnCount, std::vector<double>(rowCount, 0.0));
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t k = 0; k < result.segments.size(); ++k) {
        const PiecewiseSegmentResult& segment = result.segments[k];
        const double h = segment.halfSpan;

        // 局部 u 系数 -> 段内时间 t(秒) 系数。
        const double aX = segment.coeffX[0] / (h * h);
        const double bX =
            (-2.0 * segment.coeffX[0] + segment.coeffX[1]) / h;
        const double cX =
            segment.coeffX[0] - segment.coeffX[1] + segment.coeffX[2];

        const double aY = segment.coeffY[0] / (h * h);
        const double bY =
            (-2.0 * segment.coeffY[0] + segment.coeffY[1]) / h;
        const double cY =
            segment.coeffY[0] - segment.coeffY[1] + segment.coeffY[2];

        // 段内未加权 SSE，包含重合点。
        double sseX = 0.0;
        double sseY = 0.0;
        const int start = segment.startPoint - 1;
        for (int j = 0; j < segment.pointCount; ++j) {
            const std::size_t index = static_cast<std::size_t>(start + j);
            sseX += result.xResidual[index] * result.xResidual[index];
            sseY += result.yResidual[index] * result.yResidual[index];
        }

        columns[0][k] = static_cast<double>(segment.index);
        columns[1][k] = roundTo6(aX);
        columns[2][k] = roundTo6(bX);
        columns[3][k] = roundTo6(cX);
        columns[4][k] = roundTo6(aY);
        columns[5][k] = roundTo6(bY);
        columns[6][k] = roundTo6(cY);
        columns[7][k] = roundTo6(sseX);
        columns[8][k] = roundTo6(sseY);
        columns[9][k] = roundTo6(sseX + sseY);
    }

    // 合计行：系数列留空(NaN)，SSE 用 111 个点各算一次的结果。
    const std::size_t totalRow = result.segments.size();
    columns[0][totalRow] = 0.0;  // 0 表示合计行
    for (std::size_t c = 1; c <= 6; ++c)
        columns[c][totalRow] = nan;
    columns[7][totalRow] = roundTo6(result.sseX);
    columns[8][totalRow] = roundTo6(result.sseY);
    columns[9][totalRow] = roundTo6(result.sseX + result.sseY);

    const std::vector<std::string> headers = {
        "段号(0=合计)",
        "X二次项aX(t^2)",
        "X一次项bX(t)",
        "X常数项cX",
        "Y二次项aY(t^2)",
        "Y一次项bY(t)",
        "Y常数项cY",
        "X残差平方和",
        "Y残差平方和",
        "残差平方和合计"};

    writeColumnsToXlsx(inputPath, outputSheet, headers, columns, true);
    std::cout << "\n[write] " << inputPath << " / " << outputSheet << ": "
              << rowCount << " rows written (last row = overall SSE)\n";

    // 全局二次拟合：同样使用 L1 的 X0/Y0，时间从第一个点 t=0 开始，单位秒。
    // 结果追加到分段表下方：第 16 行表头，第 17 行数据。
    const TimeQuadFitResult global = fitTimeQuadratic({x, y}, dt);
    if (global.axes.size() < 2)
        throw std::runtime_error(
            "global quadratic fit did not return two axes");

    const std::vector<std::string> globalHeaders = {
        "全局二次",
        "X二次项aX(t^2)",
        "X一次项bX(t)",
        "X常数项cX",
        "Y二次项aY(t^2)",
        "Y一次项bY(t)",
        "Y常数项cY",
        "X残差平方和",
        "Y残差平方和",
        "残差平方和合计"};

    std::vector<std::vector<double>> globalColumns(
        globalHeaders.size(), std::vector<double>(1, nan));
    globalColumns[0][0] = nan;  // 模型名放在表头，数据行第一列留空
    globalColumns[1][0] = roundTo6(global.axes[0].coeff[0]);
    globalColumns[2][0] = roundTo6(global.axes[0].coeff[1]);
    globalColumns[3][0] = roundTo6(global.axes[0].coeff[2]);
    globalColumns[4][0] = roundTo6(global.axes[1].coeff[0]);
    globalColumns[5][0] = roundTo6(global.axes[1].coeff[1]);
    globalColumns[6][0] = roundTo6(global.axes[1].coeff[2]);
    globalColumns[7][0] = roundTo6(global.axes[0].sse);
    globalColumns[8][0] = roundTo6(global.axes[1].sse);
    globalColumns[9][0] = roundTo6(global.sse);

    writeColumnsToXlsx(inputPath, outputSheet, globalHeaders, globalColumns,
                       false, 16);
    std::cout << "[write] " << inputPath << " / " << outputSheet
              << ": global quadratic appended at rows 16-17\n";

    return result;
}
