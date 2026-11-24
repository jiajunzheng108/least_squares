// 分段二次最小二乘拟合(相邻段 C1 连续，支持任意重合点数)
#include "core/PiecewiseQuadFit.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

struct SegmentLayout {
    int start;        // 0-based 起始点
    int count;        // 点数
    double halfSpan;  // 局部半宽(秒)
};

// 段内局部坐标 u：第 j 个点(j=0..count-1)映射到 [-1, 1]。
double localU(int j, int count) {
    return (2.0 * static_cast<double>(j) -
            static_cast<double>(count - 1)) /
           static_cast<double>(count - 1);
}

// 第 k 个结合点(两段重合区中间)的全局点号(0-based)。
int junctionPoint(const std::vector<SegmentLayout>& layouts, int k,
                  int overlap) {
    return layouts[static_cast<std::size_t>(k + 1)].start + overlap / 2;
}

// 按前 K-1 段取 window 点、最后一段取 lastWindow 点的规则分段。
std::vector<SegmentLayout> buildSegments(int pointCount, double dt, int window,
                                         int step, int lastWindow) {
    if (window < 3)
        throw std::invalid_argument("window must be at least 3");
    if (lastWindow < 3)
        throw std::invalid_argument("lastWindow must be at least 3");
    if (step < 1)
        throw std::invalid_argument("step must be positive");
    const int overlap = window - step;
    if (overlap < 1 || overlap >= window)
        throw std::invalid_argument(
            "overlap = window - step must be in [1, window-1]");
    if (pointCount < lastWindow)
        throw std::invalid_argument(
            "point count is smaller than lastWindow");
    if (window > lastWindow)
        throw std::invalid_argument(
            "window must not be larger than lastWindow in this layout");
    if ((pointCount - lastWindow) % step != 0)
        throw std::invalid_argument(
            "point count is not compatible with step/lastWindow");

    const int segmentCount = (pointCount - lastWindow) / step + 1;
    std::vector<SegmentLayout> segments;
    segments.reserve(static_cast<std::size_t>(segmentCount));

    for (int k = 0; k < segmentCount - 1; ++k) {
        const int start = k * step;
        if (start + window > pointCount)
            throw std::invalid_argument("segment exceeds point range");
        segments.push_back(
            {start, window, (static_cast<double>(window) - 1.0) * dt / 2.0});
    }

    const int lastStart = (segmentCount - 1) * step;
    if (lastStart + lastWindow != pointCount)
        throw std::invalid_argument("last segment does not end at last point");
    segments.push_back({lastStart, lastWindow,
                        (static_cast<double>(lastWindow) - 1.0) * dt / 2.0});
    return segments;
}

// 原地求解 A*x=b。A 为 row-major n*n，b 为 n 维；部分主元高斯消元。
void solveLinearSystem(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double best = std::fabs(A[static_cast<std::size_t>(col) * n + col]);
        for (int row = col + 1; row < n; ++row) {
            const double v =
                std::fabs(A[static_cast<std::size_t>(row) * n + col]);
            if (v > best) {
                best = v;
                pivot = row;
            }
        }
        if (best < 1e-12)
            throw std::runtime_error(
                "piecewise quad fit: singular KKT matrix");

        if (pivot != col) {
            for (int c = col; c < n; ++c)
                std::swap(A[static_cast<std::size_t>(col) * n + c],
                          A[static_cast<std::size_t>(pivot) * n + c]);
            std::swap(b[col], b[pivot]);
        }

        const double piv = A[static_cast<std::size_t>(col) * n + col];
        for (int row = col + 1; row < n; ++row) {
            const double factor =
                A[static_cast<std::size_t>(row) * n + col] / piv;
            if (factor == 0.0)
                continue;
            A[static_cast<std::size_t>(row) * n + col] = 0.0;
            for (int c = col + 1; c < n; ++c)
                A[static_cast<std::size_t>(row) * n + c] -=
                    factor * A[static_cast<std::size_t>(col) * n + c];
            b[row] -= factor * b[col];
        }
    }

    for (int row = n - 1; row >= 0; --row) {
        double sum = b[row];
        for (int c = row + 1; c < n; ++c)
            sum -= A[static_cast<std::size_t>(row) * n + c] * b[c];
        b[row] = sum / A[static_cast<std::size_t>(row) * n + row];
    }
}

}  // namespace

PiecewiseQuadFitResult fitPiecewiseQuadratic(
    const std::vector<double>& x, const std::vector<double>& y, double dt,
    int window, int step, int lastWindow) {
    if (x.size() != y.size())
        throw std::invalid_argument("x and y must have the same size");
    const int pointCount = static_cast<int>(x.size());
    if (pointCount < 3)
        throw std::invalid_argument("at least three points are required");
    if (dt <= 0.0)
        throw std::invalid_argument("dt must be positive");

    const std::vector<SegmentLayout> layouts =
        buildSegments(pointCount, dt, window, step, lastWindow);
    const int segmentCount = static_cast<int>(layouts.size());
    const int overlap = window - step;
    const int variableCount = 3 * segmentCount;
    const int constraintCount = 2 * (segmentCount - 1);

    int rowCount = 0;
    for (const SegmentLayout& layout : layouts)
        rowCount += layout.count;

    // 每个全局点被多少个段覆盖；权重取 1/覆盖次数，保证每个观测总权重为 1。
    std::vector<int> cover(static_cast<std::size_t>(pointCount), 0);
    for (const SegmentLayout& layout : layouts)
        for (int j = 0; j < layout.count; ++j)
            ++cover[static_cast<std::size_t>(layout.start + j)];

    // 设计矩阵 A：每段 [u^2, u, 1] 三列。
    std::vector<double> A(static_cast<std::size_t>(rowCount) * variableCount,
                          0.0);
    std::vector<double> weight(static_cast<std::size_t>(rowCount), 1.0);
    std::vector<int> rowPoint(static_cast<std::size_t>(rowCount), 0);

    int row = 0;
    for (int k = 0; k < segmentCount; ++k) {
        const SegmentLayout& layout = layouts[k];
        for (int j = 0; j < layout.count; ++j) {
            const double u = localU(j, layout.count);
            A[static_cast<std::size_t>(row) * variableCount + 3 * k + 0] =
                u * u;
            A[static_cast<std::size_t>(row) * variableCount + 3 * k + 1] = u;
            A[static_cast<std::size_t>(row) * variableCount + 3 * k + 2] =
                1.0;

            rowPoint[row] = layout.start + j;
            weight[row] = 1.0 / static_cast<double>(
                                    cover[static_cast<std::size_t>(
                                        rowPoint[row])]);
            ++row;
        }
    }

    // 法方程 H = A^T W A。
    std::vector<double> H(static_cast<std::size_t>(variableCount) *
                              variableCount,
                          0.0);
    for (int i = 0; i < rowCount; ++i) {
        const double wi = weight[i];
        const std::size_t baseI = static_cast<std::size_t>(i) * variableCount;
        for (int p = 0; p < variableCount; ++p) {
            const double aip = A[baseI + p];
            if (aip == 0.0)
                continue;
            const double wa = wi * aip;
            for (int q = 0; q < variableCount; ++q)
                H[static_cast<std::size_t>(p) * variableCount + q] +=
                    wa * A[baseI + q];
        }
    }

    // 约束 C*theta=0：在重合区中间点处位置连续 + 速度连续。
    std::vector<double> C(static_cast<std::size_t>(constraintCount) *
                              variableCount,
                          0.0);
    int con = 0;
    for (int k = 0; k + 1 < segmentCount; ++k) {
        const int baseK = 3 * k;
        const int baseNext = 3 * (k + 1);

        const int junction = junctionPoint(layouts, k, overlap);
        const double uL =
            localU(junction - layouts[k].start, layouts[k].count);
        const double uR =
            localU(junction - layouts[k + 1].start, layouts[k + 1].count);

        // 位置：p_k(uL) - p_{k+1}(uR) = 0
        const std::size_t posRow =
            static_cast<std::size_t>(con) * variableCount;
        C[posRow + baseK + 0] = uL * uL;
        C[posRow + baseK + 1] = uL;
        C[posRow + baseK + 2] = 1.0;
        C[posRow + baseNext + 0] = -uR * uR;
        C[posRow + baseNext + 1] = -uR;
        C[posRow + baseNext + 2] = -1.0;
        ++con;

        // 速度：h_{k+1}(2a_k uL+b_k) - h_k(2a_{k+1} uR+b_{k+1}) = 0
        // 这里用毫秒作为长度单位，仅对约束行做统一缩放，不改变解。
        const double hk = layouts[k].halfSpan * 1000.0;
        const double hkNext = layouts[k + 1].halfSpan * 1000.0;
        const std::size_t velRow =
            static_cast<std::size_t>(con) * variableCount;
        C[velRow + baseK + 0] = 2.0 * hkNext * uL;
        C[velRow + baseK + 1] = hkNext;
        C[velRow + baseNext + 0] = -2.0 * hk * uR;
        C[velRow + baseNext + 1] = -hk;
        ++con;
    }

    // KKT 矩阵：
    // [ H  C^T ]
    // [ C   0  ]
    const int total = variableCount + constraintCount;
    std::vector<double> kkt(static_cast<std::size_t>(total) * total, 0.0);
    for (int p = 0; p < variableCount; ++p)
        for (int q = 0; q < variableCount; ++q)
            kkt[static_cast<std::size_t>(p) * total + q] =
                H[static_cast<std::size_t>(p) * variableCount + q];
    for (int c = 0; c < constraintCount; ++c) {
        for (int p = 0; p < variableCount; ++p) {
            const double value =
                C[static_cast<std::size_t>(c) * variableCount + p];
            kkt[static_cast<std::size_t>(p) * total + variableCount + c] =
                value;
            kkt[static_cast<std::size_t>(variableCount + c) * total + p] =
                value;
        }
    }

    // 对 X、Y 各解一次；两者共用同一个 KKT 矩阵。
    auto solveAxis = [&](const std::vector<double>& values) {
        std::vector<double> gradient(static_cast<std::size_t>(variableCount),
                                     0.0);
        for (int i = 0; i < rowCount; ++i) {
            const double wi = weight[i];
            if (wi == 0.0)
                continue;
            const double bi = values[rowPoint[i]];
            const std::size_t baseI =
                static_cast<std::size_t>(i) * variableCount;
            for (int p = 0; p < variableCount; ++p)
                gradient[p] += wi * A[baseI + p] * bi;
        }

        std::vector<double> rhs(static_cast<std::size_t>(total), 0.0);
        for (int p = 0; p < variableCount; ++p)
            rhs[p] = gradient[p];

        std::vector<double> matrix = kkt;
        solveLinearSystem(matrix, rhs, total);
        rhs.resize(static_cast<std::size_t>(variableCount));
        return rhs;
    };

    const std::vector<double> thetaX = solveAxis(x);
    const std::vector<double> thetaY = solveAxis(y);

    PiecewiseQuadFitResult result;
    result.pointCount = pointCount;
    result.dt = dt;
    result.window = window;
    result.step = step;
    result.lastWindow = lastWindow;
    result.overlap = overlap;
    result.segments.resize(static_cast<std::size_t>(segmentCount));
    result.nodes.resize(static_cast<std::size_t>(segmentCount - 1));

    result.xFit.assign(static_cast<std::size_t>(pointCount), 0.0);
    result.yFit.assign(static_cast<std::size_t>(pointCount), 0.0);
    result.xResidual.assign(static_cast<std::size_t>(pointCount), 0.0);
    result.yResidual.assign(static_cast<std::size_t>(pointCount), 0.0);
    result.pointSegment.assign(static_cast<std::size_t>(pointCount), 0);
    result.pointType.assign(static_cast<std::size_t>(pointCount), 0);

    for (int k = 0; k < segmentCount; ++k) {
        const SegmentLayout& layout = layouts[k];
        PiecewiseSegmentResult& segment =
            result.segments[static_cast<std::size_t>(k)];
        segment.index = k + 1;
        segment.startPoint = layout.start + 1;
        segment.pointCount = layout.count;
        segment.halfSpan = layout.halfSpan;
        for (int c = 0; c < 3; ++c) {
            segment.coeffX[c] =
                thetaX[static_cast<std::size_t>(3 * k + c)];
            segment.coeffY[c] =
                thetaY[static_cast<std::size_t>(3 * k + c)];
        }

        double sseX = 0.0;
        double sseY = 0.0;
        for (int j = 0; j < layout.count; ++j) {
            const int index = layout.start + j;
            const double u = localU(j, layout.count);
            const double fx =
                segment.coeffX[0] * u * u + segment.coeffX[1] * u +
                segment.coeffX[2];
            const double fy =
                segment.coeffY[0] * u * u + segment.coeffY[1] * u +
                segment.coeffY[2];
            const double rx = x[static_cast<std::size_t>(index)] - fx;
            const double ry = y[static_cast<std::size_t>(index)] - fy;
            sseX += rx * rx;
            sseY += ry * ry;
        }
        segment.rmseX =
            std::sqrt(sseX / static_cast<double>(layout.count));
        segment.rmseY =
            std::sqrt(sseY / static_cast<double>(layout.count));
    }

    // 点级结果：以结合点 g 为界，i <= g 归左段、i > g 归右段。
    std::vector<int> junctions(static_cast<std::size_t>(segmentCount - 1), 0);
    for (int k = 0; k + 1 < segmentCount; ++k)
        junctions[static_cast<std::size_t>(k)] =
            junctionPoint(layouts, k, overlap);

    for (int i = 0; i < pointCount; ++i) {
        const int owner = static_cast<int>(
            std::lower_bound(junctions.begin(), junctions.end(), i) -
            junctions.begin());
        if (owner < 0 || owner >= segmentCount)
            throw std::runtime_error(
                "piecewise quad fit: some points are not covered");
        const SegmentLayout& layout =
            layouts[static_cast<std::size_t>(owner)];
        const PiecewiseSegmentResult& segment =
            result.segments[static_cast<std::size_t>(owner)];
        const double u = localU(i - layout.start, layout.count);
        result.xFit[static_cast<std::size_t>(i)] =
            segment.coeffX[0] * u * u + segment.coeffX[1] * u +
            segment.coeffX[2];
        result.yFit[static_cast<std::size_t>(i)] =
            segment.coeffY[0] * u * u + segment.coeffY[1] * u +
            segment.coeffY[2];
        result.pointSegment[static_cast<std::size_t>(i)] = owner + 1;
    }

    for (int k = 0; k + 1 < segmentCount; ++k)
        result.pointType[static_cast<std::size_t>(
            junctions[static_cast<std::size_t>(k)])] = 1;

    result.sseX = 0.0;
    result.sseY = 0.0;
    result.maxAbsResidualX = 0.0;
    result.maxAbsResidualY = 0.0;
    for (int i = 0; i < pointCount; ++i) {
        const double rx =
            x[static_cast<std::size_t>(i)] -
            result.xFit[static_cast<std::size_t>(i)];
        const double ry =
            y[static_cast<std::size_t>(i)] -
            result.yFit[static_cast<std::size_t>(i)];
        result.xResidual[static_cast<std::size_t>(i)] = rx;
        result.yResidual[static_cast<std::size_t>(i)] = ry;
        result.sseX += rx * rx;
        result.sseY += ry * ry;
        result.maxAbsResidualX =
            std::max(result.maxAbsResidualX, std::fabs(rx));
        result.maxAbsResidualY =
            std::max(result.maxAbsResidualY, std::fabs(ry));
    }
    result.rmseX =
        std::sqrt(result.sseX / static_cast<double>(pointCount));
    result.rmseY =
        std::sqrt(result.sseY / static_cast<double>(pointCount));

    for (int k = 0; k + 1 < segmentCount; ++k) {
        PiecewiseNodeResult& node = result.nodes[static_cast<std::size_t>(k)];
        const int junction = junctionPoint(layouts, k, overlap);
        node.index = k + 1;
        node.point = junction + 1;
        node.t = static_cast<double>(junction) * dt;

        const double hLeft = layouts[k].halfSpan;
        const double hRight = layouts[k + 1].halfSpan;
        const double uL =
            localU(junction - layouts[k].start, layouts[k].count);
        const double uR =
            localU(junction - layouts[k + 1].start, layouts[k + 1].count);

        const double* leftX = &thetaX[static_cast<std::size_t>(3 * k)];
        const double* leftY = &thetaY[static_cast<std::size_t>(3 * k)];
        const double* rightX =
            &thetaX[static_cast<std::size_t>(3 * (k + 1))];
        const double* rightY =
            &thetaY[static_cast<std::size_t>(3 * (k + 1))];

        // 左侧段在 uL 处
        const double xLeft =
            leftX[0] * uL * uL + leftX[1] * uL + leftX[2];
        const double yLeft =
            leftY[0] * uL * uL + leftY[1] * uL + leftY[2];
        const double vxLeft = (2.0 * leftX[0] * uL + leftX[1]) / hLeft;
        const double vyLeft = (2.0 * leftY[0] * uL + leftY[1]) / hLeft;
        const double axLeft = 2.0 * leftX[0] / (hLeft * hLeft);
        const double ayLeft = 2.0 * leftY[0] / (hLeft * hLeft);

        // 右侧段在 uR 处
        const double xRight =
            rightX[0] * uR * uR + rightX[1] * uR + rightX[2];
        const double yRight =
            rightY[0] * uR * uR + rightY[1] * uR + rightY[2];
        const double vxRight = (2.0 * rightX[0] * uR + rightX[1]) / hRight;
        const double vyRight = (2.0 * rightY[0] * uR + rightY[1]) / hRight;
        const double axRight = 2.0 * rightX[0] / (hRight * hRight);
        const double ayRight = 2.0 * rightY[0] / (hRight * hRight);

        node.x = xLeft;
        node.y = yLeft;
        node.vx = vxLeft;
        node.vy = vyLeft;
        node.ax = axLeft;
        node.ay = ayLeft;
        node.dax = axRight - axLeft;
        node.day = ayRight - ayLeft;
        node.ep = std::hypot(xRight - xLeft, yRight - yLeft);
        node.ev = std::hypot(vxRight - vxLeft, vyRight - vyLeft);
        node.ea = std::hypot(node.dax, node.day);
    }

    return result;
}
