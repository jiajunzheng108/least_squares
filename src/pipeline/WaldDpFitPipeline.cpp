// Wald + DP + KKT 分段拟合的操作封装
#include "pipeline/WaldDpFitPipeline.h"

#include "utils/XlsxReader.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

double roundTo6(double value) {
    if (std::isnan(value))
        return value;
    return std::round(value * 1.0e6) / 1.0e6;
}

std::string joinNodes(const std::vector<int>& nodes) {
    std::ostringstream os;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (i > 0)
            os << ',';
        os << nodes[i];
    }
    return os.str();
}

void printTableHeader() {
    std::cout << "  K | DP路径 | KKT复核路径 | SSE_X | SSE_Y | SSE合计 | "
                 "RMSE | e_p | e_v | e_a\n";
}

}  // namespace

WaldDpFitResult runWaldDpFitPipeline(
    const std::string& inputPath, const std::string& inputSheet,
    const std::string& xHeader, const std::string& yHeader, double dt,
    int windowHalfWidth, int minSegmentLength, int maxSegmentLength,
    double alpha, double muP, double muV) {
    std::vector<double> x;
    std::vector<double> y;
    readPointsFromXlsx(inputPath, inputSheet, xHeader, yHeader, x, y);

    WaldDpConfig config;
    config.dt = dt;
    config.windowHalfWidth = windowHalfWidth;
    config.minSegmentLength = minSegmentLength;
    config.maxSegmentLength = maxSegmentLength;
    config.alpha = alpha;
    config.muP = muP;
    config.muV = muV;

    const WaldDpFitResult result = fitWaldDp(x, y, config);

    std::cout << "=== Wald + DP + KKT: " << inputPath << " / " << inputSheet
              << " (" << xHeader << ", " << yHeader << ") ===\n";
    std::cout << "points = " << result.pointCount << "  dt = " << result.config.dt
              << " s  m = " << result.config.windowHalfWidth
              << "  L_min = " << result.config.minSegmentLength
              << "  L_max = "
              << (result.config.maxSegmentLength == 0
                      ? std::string("unlimited")
                      : std::to_string(result.config.maxSegmentLength))
              << "\n";
    std::cout << "alpha = " << result.config.alpha
              << "  alpha' = " << result.config.alpha / result.pointCount
              << "  chi2(6) critical = " << result.waldCritical << "\n";
    std::cout << "scales: s_p = " << result.positionScale
              << "  s_v = " << result.velocityScale
              << "  mu_p = " << result.config.muP
              << "  mu_v = " << result.config.muV << "\n";
    std::cout << "Wald raw candidates = " << result.rawCandidateCount
              << "  kept after NMS = " << result.candidates.size() << "\n";
    for (const WaldCandidateInfo& c : result.candidates) {
        std::cout << "  point " << std::setw(4) << c.point
                  << "  W = " << std::setw(12) << c.wald
                  << "  p = " << std::setw(12) << c.pValue << "\n";
    }

    std::cout << "\n[K comparison] DP 路径 vs KKT 复核最优路径\n";
    printTableHeader();
    std::cout << std::fixed << std::setprecision(6);
    const int maxK = static_cast<int>(result.bestPaths.size()) - 1;
    for (int k = 1; k <= maxK; ++k) {
        const WaldKPathResult& best =
            result.bestPaths[static_cast<std::size_t>(k)];
        if (best.K == 0)
            continue;
        const WaldKPathResult& dp =
            result.dpPaths[static_cast<std::size_t>(k)];
        std::cout << std::setw(4) << k << " | "
                  << (dp.K == 0 ? std::string("-") : joinNodes(dp.nodes))
                  << " | " << joinNodes(best.nodes) << " | "
                  << std::setw(10) << best.sseX << " | " << std::setw(10)
                  << best.sseY << " | " << std::setw(10) << best.sseTotal
                  << " | " << std::setw(10) << best.rmseTotal << " | "
                  << std::setw(10) << best.epMax << " | " << std::setw(10)
                  << best.evMax << " | " << std::setw(10) << best.eaMax
                  << "\n";
    }

    std::cout << "\n[chosen] K = " << result.chosenK << "  nodes = "
              << joinNodes(result.chosen.nodes) << "\n";
    std::cout << "  SSE_X = " << result.chosen.sseX
              << "  SSE_Y = " << result.chosen.sseY
              << "  SSE_total = " << result.chosen.sseTotal << "\n";
    std::cout << "  RMSE_X = " << result.chosen.rmseX
              << "  RMSE_Y = " << result.chosen.rmseY
              << "  RMSE_total = " << result.chosen.rmseTotal << "\n";
    std::cout << "  e_p_max = " << result.chosen.epMax
              << "  e_v_max = " << result.chosen.evMax
              << "  e_a_max = " << result.chosen.eaMax << "\n";

    std::cout << "\n[segments]\n";
    std::cout << "  seg  points    n   SSE_X      SSE_Y\n";
    for (std::size_t k = 0; k < result.chosen.segments.size(); ++k) {
        const WaldSegmentInfo& s = result.chosen.segments[k];
        std::cout << "  " << std::setw(3) << (k + 1) << "  " << std::setw(4)
                  << s.startPoint << "-" << std::setw(4) << s.endPoint
                  << "  " << std::setw(4) << s.pointCount << "  "
                  << std::setw(10) << s.sseX << "  " << std::setw(10)
                  << s.sseY << "\n";
    }

    std::cout << "\n[nodes]\n";
    std::cout << "  point   e_p          e_v          e_a          daX"
                 "          daY\n";
    for (const WaldNodeInfo& node : result.chosen.nodeInfo) {
        std::cout << "  " << std::setw(4) << node.point << "  "
                  << std::setw(10) << node.ep << "  " << std::setw(10)
                  << node.ev << "  " << std::setw(10) << node.ea << "  "
                  << std::setw(10) << node.dax << "  " << std::setw(10)
                  << node.day << "\n";
    }

    // 追加写入“分段拟合”表：第 19 行表头，20-21 段，22 合计。
    if (result.chosen.K > 0 && !result.chosen.segments.empty()) {
        const int dpBlockStartRow = 19;
        const std::size_t segCount = result.chosen.segments.size();
        const std::size_t dataRows = segCount + 1;
        std::vector<std::vector<double>> columns(
            10, std::vector<double>(dataRows, 0.0));
        const double nan = std::numeric_limits<double>::quiet_NaN();

        for (std::size_t k = 0; k < segCount; ++k) {
            const WaldSegmentInfo& s = result.chosen.segments[k];
            const double h = s.halfSpan;
            const double aX = s.coeffX[0] / (h * h);
            const double bX = (-2.0 * s.coeffX[0] + s.coeffX[1]) / h;
            const double cX = s.coeffX[0] - s.coeffX[1] + s.coeffX[2];
            const double aY = s.coeffY[0] / (h * h);
            const double bY = (-2.0 * s.coeffY[0] + s.coeffY[1]) / h;
            const double cY = s.coeffY[0] - s.coeffY[1] + s.coeffY[2];

            columns[0][k] = static_cast<double>(k + 1);
            columns[1][k] = roundTo6(aX);
            columns[2][k] = roundTo6(bX);
            columns[3][k] = roundTo6(cX);
            columns[4][k] = roundTo6(aY);
            columns[5][k] = roundTo6(bY);
            columns[6][k] = roundTo6(cY);
            columns[7][k] = roundTo6(s.sseX);
            columns[8][k] = roundTo6(s.sseY);
            columns[9][k] = roundTo6(s.sseX + s.sseY);
        }

        const std::size_t total = segCount;
        columns[0][total] = 0.0;  // 0 表示合计行
        for (std::size_t c = 1; c <= 6; ++c)
            columns[c][total] = nan;
        columns[7][total] = roundTo6(result.chosen.sseX);
        columns[8][total] = roundTo6(result.chosen.sseY);
        columns[9][total] = roundTo6(result.chosen.sseTotal);

        const std::vector<std::string> headers = {
            "DP分段(0=合计)",
            "X二次项aX(t^2)",
            "X一次项bX(t)",
            "X常数项cX",
            "Y二次项aY(t^2)",
            "Y一次项bY(t)",
            "Y常数项cY",
            "X残差平方和",
            "Y残差平方和",
            "残差平方和合计"};

        writeColumnsToXlsx(inputPath, "分段拟合", headers, columns, false,
                           dpBlockStartRow);
        std::cout << "\n[write] " << inputPath << " / 分段拟合: DP block at rows "
                  << dpBlockStartRow << "-" << dpBlockStartRow + dataRows
                  << "\n";
    } else {
        std::cout << "\n[warn] no chosen DP path, nothing written.\n";
    }

    if (result.enumerationTruncated)
        std::cout << "\n[warn] candidate path enumeration was truncated.\n";

    return result;
}
