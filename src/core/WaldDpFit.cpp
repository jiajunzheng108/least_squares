// Wald 候选点筛选 + DP 最优分段 + KKT 硬约束拟合(二维轨迹)
#include "core/WaldDpFit.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// ---------------- 小型线性代数 ----------------

// 原地求解 A*x=b。A 为 row-major n*n；部分主元高斯消元。
void solveLinearSystem(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double best = std::fabs(A[static_cast<std::size_t>(col) * n + col]);
        for (int r = col + 1; r < n; ++r) {
            const double v =
                std::fabs(A[static_cast<std::size_t>(r) * n + col]);
            if (v > best) {
                best = v;
                pivot = r;
            }
        }
        if (best < 1e-14)
            throw std::runtime_error("singular linear system");
        if (pivot != col) {
            for (int c = col; c < n; ++c)
                std::swap(A[static_cast<std::size_t>(col) * n + c],
                          A[static_cast<std::size_t>(pivot) * n + c]);
            std::swap(b[col], b[pivot]);
        }
        const double piv = A[static_cast<std::size_t>(col) * n + col];
        for (int r = col + 1; r < n; ++r) {
            const double f =
                A[static_cast<std::size_t>(r) * n + col] / piv;
            if (f == 0.0)
                continue;
            A[static_cast<std::size_t>(r) * n + col] = 0.0;
            for (int c = col + 1; c < n; ++c)
                A[static_cast<std::size_t>(r) * n + c] -=
                    f * A[static_cast<std::size_t>(col) * n + c];
            b[r] -= f * b[col];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = b[r];
        for (int c = r + 1; c < n; ++c)
            s -= A[static_cast<std::size_t>(r) * n + c] * b[c];
        b[r] = s / A[static_cast<std::size_t>(r) * n + r];
    }
}

bool invertMatrix(const std::vector<double>& A, int n,
                  std::vector<double>& inv) {
    inv.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int c = 0; c < n; ++c) {
        std::vector<double> M = A;
        std::vector<double> e(static_cast<std::size_t>(n), 0.0);
        e[static_cast<std::size_t>(c)] = 1.0;
        try {
            solveLinearSystem(M, e, n);
        } catch (...) {
            return false;
        }
        for (int r = 0; r < n; ++r)
            inv[static_cast<std::size_t>(r) * n + c] = e[r];
    }
    return true;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        s += a[i] * b[i];
    return s;
}

// ---------------- 卡方分布 ----------------

double lnGamma(double x) {
    static const double cof[6] = {
        76.18009172947146,      -86.50532032941677, 24.01409824083091,
        -1.231739572450155,     0.1208650973866179e-2, -0.5395239384953e-5};
    double y = x;
    double tmp = x + 5.5;
    tmp -= (x + 0.5) * std::log(tmp);
    double ser = 1.000000000190015;
    for (int j = 0; j < 6; ++j)
        ser += cof[j] / ++y;
    return -tmp + std::log(2.5066282746310005 * ser / x);
}

double gammaQ(double a, double x);

double gammaP(double a, double x) {
    if (x <= 0.0)
        return 0.0;
    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n = 1; n < 300; ++n) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * 1e-15)
                break;
        }
        return sum * std::exp(-x + a * std::log(x) - lnGamma(a));
    }
    return 1.0 - gammaQ(a, x);
}

double gammaQ(double a, double x) {
    if (x < a + 1.0)
        return 1.0 - gammaP(a, x);
    const double tiny = 1e-300;
    double b = x + 1.0 - a;
    double c = 1.0 / tiny;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i < 300; ++i) {
        const double an = -static_cast<double>(i) * (i - a);
        b += 2.0;
        d = an * d + b;
        if (std::fabs(d) < tiny)
            d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny)
            c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-15)
            break;
    }
    return std::exp(-x + a * std::log(x) - lnGamma(a)) * h;
}

double chiSquareCDF(double x, int df) {
    if (x <= 0.0)
        return 0.0;
    return gammaP(df / 2.0, x / 2.0);
}

double chiSquareSF(double x, int df) {
    if (x <= 0.0)
        return 1.0;
    return gammaQ(df / 2.0, x / 2.0);
}

// 返回上尾概率为 tail 的卡方分位数：P(X > q) = tail。
double chiSquareUpperQuantile(double tail, int df) {
    double lo = 0.0;
    double hi = 1.0;
    const double target = 1.0 - tail;
    while (chiSquareCDF(hi, df) < target && hi < 1e6)
        hi *= 2.0;
    for (int it = 0; it < 120; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (chiSquareCDF(mid, df) < target)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

// ---------------- Wald 窗口拟合 ----------------

struct WindowPair {
    double coefL[3][2] = {};
    double coefR[3][2] = {};
    double covL[2][9] = {};
    double covR[2][9] = {};
    bool ok = false;
};

// 对一组索引做二维二次拟合，返回系数和系数协方差。
bool fitOneWindow(const std::vector<double>& x, const std::vector<double>& y,
                  int first, int last, int center, int m,
                  double coef[3][2], double cov[2][9]) {
    const int n = last - first + 1;
    if (n < 4)
        return false;

    std::vector<double> X(static_cast<std::size_t>(n) * 3, 0.0);
    for (int j = 0; j < n; ++j) {
        const double u =
            static_cast<double>(first + j - center) / static_cast<double>(m);
        X[static_cast<std::size_t>(j) * 3 + 0] = u * u;
        X[static_cast<std::size_t>(j) * 3 + 1] = u;
        X[static_cast<std::size_t>(j) * 3 + 2] = 1.0;
    }

    // X^T X
    std::vector<double> XtX(9, 0.0);
    for (int j = 0; j < n; ++j)
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
                XtX[static_cast<std::size_t>(p) * 3 + q] +=
                    X[static_cast<std::size_t>(j) * 3 + p] *
                    X[static_cast<std::size_t>(j) * 3 + q];

    std::vector<double> invXtX;
    if (!invertMatrix(XtX, 3, invXtX))
        return false;

    for (int axis = 0; axis < 2; ++axis) {
        const std::vector<double>& values = (axis == 0) ? x : y;
        std::vector<double> rhs(3, 0.0);
        for (int j = 0; j < n; ++j) {
            const double v = values[static_cast<std::size_t>(first + j)];
            for (int p = 0; p < 3; ++p)
                rhs[p] += X[static_cast<std::size_t>(j) * 3 + p] * v;
        }
        std::vector<double> M = XtX;
        try {
            solveLinearSystem(M, rhs, 3);
        } catch (...) {
            return false;
        }
        for (int p = 0; p < 3; ++p)
            coef[p][axis] = rhs[p];

        double sse = 0.0;
        for (int j = 0; j < n; ++j) {
            const double u =
                static_cast<double>(first + j - center) / static_cast<double>(m);
            const double fit =
                coef[0][axis] * u * u + coef[1][axis] * u + coef[2][axis];
            const double e = values[static_cast<std::size_t>(first + j)] - fit;
            sse += e * e;
        }
        const double sigma2 = sse / static_cast<double>(n - 3);
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
                cov[axis][static_cast<std::size_t>(p) * 3 + q] =
                    sigma2 * invXtX[static_cast<std::size_t>(p) * 3 + q];
    }
    return true;
}

// 左窗口 [center-m, center]，右窗口 [center, center+m]，共用归一化时间。
WindowPair fitWindowPair(const std::vector<double>& x,
                         const std::vector<double>& y, int center, int m) {
    WindowPair out;
    const bool leftOk =
        fitOneWindow(x, y, center - m, center, center, m, out.coefL, out.covL);
    const bool rightOk =
        fitOneWindow(x, y, center, center + m, center, m, out.coefR, out.covR);
    out.ok = leftOk && rightOk;
    return out;
}

// ---------------- 段拟合 ----------------

struct SegFit {
    bool valid = false;
    int a = 0;  // 1-based
    int b = 0;
    int n = 0;
    double h = 0.0;
    double sseX = 0.0;
    double sseY = 0.0;
    double pL[2] = {};
    double vL[2] = {};
    double pR[2] = {};
    double vR[2] = {};
    double coeffX[3] = {};
    double coeffY[3] = {};
};

SegFit fitSegment(const std::vector<double>& x, const std::vector<double>& y,
                  int a, int b, double dt) {
    SegFit s;
    s.a = a;
    s.b = b;
    s.n = b - a + 1;
    if (s.n < 3)
        return s;
    s.h = (static_cast<double>(s.n) - 1.0) * dt / 2.0;

    std::vector<double> X(static_cast<std::size_t>(s.n) * 3, 0.0);
    for (int j = 0; j < s.n; ++j) {
        const double u =
            (2.0 * static_cast<double>(j) - static_cast<double>(s.n - 1)) /
            static_cast<double>(s.n - 1);
        X[static_cast<std::size_t>(j) * 3 + 0] = u * u;
        X[static_cast<std::size_t>(j) * 3 + 1] = u;
        X[static_cast<std::size_t>(j) * 3 + 2] = 1.0;
    }

    std::vector<double> XtX(9, 0.0);
    for (int j = 0; j < s.n; ++j)
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
                XtX[static_cast<std::size_t>(p) * 3 + q] +=
                    X[static_cast<std::size_t>(j) * 3 + p] *
                    X[static_cast<std::size_t>(j) * 3 + q];

    for (int axis = 0; axis < 2; ++axis) {
        const std::vector<double>& values = (axis == 0) ? x : y;
        double* coef = (axis == 0) ? s.coeffX : s.coeffY;
        std::vector<double> rhs(3, 0.0);
        for (int j = 0; j < s.n; ++j) {
            const double v = values[static_cast<std::size_t>(a - 1 + j)];
            for (int p = 0; p < 3; ++p)
                rhs[p] += X[static_cast<std::size_t>(j) * 3 + p] * v;
        }
        std::vector<double> M = XtX;
        try {
            solveLinearSystem(M, rhs, 3);
        } catch (...) {
            return s;
        }
        for (int p = 0; p < 3; ++p)
            coef[p] = rhs[p];

        double sse = 0.0;
        for (int j = 0; j < s.n; ++j) {
            const double u =
                (2.0 * static_cast<double>(j) -
                 static_cast<double>(s.n - 1)) /
                static_cast<double>(s.n - 1);
            const double fit = coef[0] * u * u + coef[1] * u + coef[2];
            const double e =
                values[static_cast<std::size_t>(a - 1 + j)] - fit;
            sse += e * e;
        }
        if (axis == 0)
            s.sseX = sse;
        else
            s.sseY = sse;
    }

    s.pL[0] = s.coeffX[0] - s.coeffX[1] + s.coeffX[2];
    s.pL[1] = s.coeffY[0] - s.coeffY[1] + s.coeffY[2];
    s.pR[0] = s.coeffX[0] + s.coeffX[1] + s.coeffX[2];
    s.pR[1] = s.coeffY[0] + s.coeffY[1] + s.coeffY[2];
    s.vL[0] = (-2.0 * s.coeffX[0] + s.coeffX[1]) / s.h;
    s.vL[1] = (-2.0 * s.coeffY[0] + s.coeffY[1]) / s.h;
    s.vR[0] = (2.0 * s.coeffX[0] + s.coeffX[1]) / s.h;
    s.vR[1] = (2.0 * s.coeffY[0] + s.coeffY[1]) / s.h;
    s.valid = true;
    return s;
}

// ---------------- KKT 联合拟合 ----------------

WaldKPathResult fitKkt(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::vector<int>& nodes, double dt) {
    WaldKPathResult result;
    const int K = static_cast<int>(nodes.size()) - 1;
    result.K = K;
    result.nodes = nodes;

    struct KSeg {
        int a;
        int n;
        double h;
        std::vector<double> u;
    };
    std::vector<KSeg> segs(static_cast<std::size_t>(K));
    int rowCount = 0;
    for (int k = 0; k < K; ++k) {
        const int a = nodes[static_cast<std::size_t>(k)];
        const int b = nodes[static_cast<std::size_t>(k + 1)];
        KSeg s;
        s.a = a;
        s.n = b - a + 1;
        s.h = (static_cast<double>(s.n) - 1.0) * dt / 2.0;
        s.u.resize(static_cast<std::size_t>(s.n));
        for (int j = 0; j < s.n; ++j)
            s.u[static_cast<std::size_t>(j)] =
                (2.0 * static_cast<double>(j) -
                 static_cast<double>(s.n - 1)) /
                static_cast<double>(s.n - 1);
        segs[static_cast<std::size_t>(k)] = std::move(s);
        rowCount += segs[static_cast<std::size_t>(k)].n;
    }

    const int nvar = 3 * K;
    const int ncon = 2 * (K - 1);
    std::vector<double> A(static_cast<std::size_t>(rowCount) * nvar, 0.0);
    std::vector<double> w(static_cast<std::size_t>(rowCount), 1.0);
    std::vector<int> rowPoint(static_cast<std::size_t>(rowCount), 0);

    int row = 0;
    for (int k = 0; k < K; ++k) {
        const KSeg& s = segs[static_cast<std::size_t>(k)];
        for (int j = 0; j < s.n; ++j) {
            const double u = s.u[static_cast<std::size_t>(j)];
            A[static_cast<std::size_t>(row) * nvar + 3 * k + 0] = u * u;
            A[static_cast<std::size_t>(row) * nvar + 3 * k + 1] = u;
            A[static_cast<std::size_t>(row) * nvar + 3 * k + 2] = 1.0;
            const bool sharedPrev = (k > 0 && j == 0);
            const bool sharedNext = (k + 1 < K && j == s.n - 1);
            if (sharedPrev || sharedNext)
                w[static_cast<std::size_t>(row)] = 0.5;
            rowPoint[static_cast<std::size_t>(row)] = s.a - 1 + j;
            ++row;
        }
    }

    std::vector<double> H(static_cast<std::size_t>(nvar) * nvar, 0.0);
    for (int i = 0; i < rowCount; ++i) {
        const double wi = w[static_cast<std::size_t>(i)];
        const std::size_t base = static_cast<std::size_t>(i) * nvar;
        for (int p = 0; p < nvar; ++p) {
            const double aip = A[base + p];
            if (aip == 0.0)
                continue;
            const double wa = wi * aip;
            for (int q = 0; q < nvar; ++q)
                H[static_cast<std::size_t>(p) * nvar + q] +=
                    wa * A[base + q];
        }
    }

    std::vector<double> C(static_cast<std::size_t>(ncon) * nvar, 0.0);
    int con = 0;
    for (int k = 0; k + 1 < K; ++k) {
        const int b0 = 3 * k;
        const int b1 = 3 * (k + 1);
        C[static_cast<std::size_t>(con) * nvar + b0 + 0] = 1.0;
        C[static_cast<std::size_t>(con) * nvar + b0 + 1] = 1.0;
        C[static_cast<std::size_t>(con) * nvar + b0 + 2] = 1.0;
        C[static_cast<std::size_t>(con) * nvar + b1 + 0] = -1.0;
        C[static_cast<std::size_t>(con) * nvar + b1 + 1] = 1.0;
        C[static_cast<std::size_t>(con) * nvar + b1 + 2] = -1.0;
        ++con;

        const double hk = segs[static_cast<std::size_t>(k)].h * 1000.0;
        const double hk1 =
            segs[static_cast<std::size_t>(k + 1)].h * 1000.0;
        C[static_cast<std::size_t>(con) * nvar + b0 + 0] = 2.0 * hk1;
        C[static_cast<std::size_t>(con) * nvar + b0 + 1] = hk1;
        C[static_cast<std::size_t>(con) * nvar + b1 + 0] = 2.0 * hk;
        C[static_cast<std::size_t>(con) * nvar + b1 + 1] = -hk;
        ++con;
    }

    const int total = nvar + ncon;
    std::vector<double> kkt(static_cast<std::size_t>(total) * total, 0.0);
    for (int p = 0; p < nvar; ++p)
        for (int q = 0; q < nvar; ++q)
            kkt[static_cast<std::size_t>(p) * total + q] =
                H[static_cast<std::size_t>(p) * nvar + q];
    for (int c = 0; c < ncon; ++c) {
        for (int p = 0; p < nvar; ++p) {
            const double v =
                C[static_cast<std::size_t>(c) * nvar + p];
            kkt[static_cast<std::size_t>(p) * total + nvar + c] = v;
            kkt[static_cast<std::size_t>(nvar + c) * total + p] = v;
        }
    }

    auto solveAxis = [&](const std::vector<double>& values) {
        std::vector<double> g(static_cast<std::size_t>(nvar), 0.0);
        for (int i = 0; i < rowCount; ++i) {
            const double wi = w[static_cast<std::size_t>(i)];
            if (wi == 0.0)
                continue;
            const double bi = values[static_cast<std::size_t>(rowPoint[i])];
            const std::size_t base = static_cast<std::size_t>(i) * nvar;
            for (int p = 0; p < nvar; ++p)
                g[p] += wi * A[base + p] * bi;
        }
        std::vector<double> rhs(static_cast<std::size_t>(total), 0.0);
        for (int p = 0; p < nvar; ++p)
            rhs[p] = g[p];
        std::vector<double> M = kkt;
        solveLinearSystem(M, rhs, total);
        rhs.resize(static_cast<std::size_t>(nvar));
        return rhs;
    };

    const std::vector<double> cx = solveAxis(x);
    const std::vector<double> cy = solveAxis(y);

    const int N = static_cast<int>(x.size());
    std::vector<double> fitX(static_cast<std::size_t>(N), 0.0);
    std::vector<double> fitY(static_cast<std::size_t>(N), 0.0);
    std::vector<bool> assigned(static_cast<std::size_t>(N), false);

    result.segments.resize(static_cast<std::size_t>(K));
    for (int k = 0; k < K; ++k) {
        const KSeg& s = segs[static_cast<std::size_t>(k)];
        WaldSegmentInfo info;
        info.startPoint = nodes[static_cast<std::size_t>(k)];
        info.endPoint = nodes[static_cast<std::size_t>(k + 1)];
        info.pointCount = s.n;
        info.halfSpan = s.h;
        for (int p = 0; p < 3; ++p) {
            info.coeffX[p] = cx[static_cast<std::size_t>(3 * k + p)];
            info.coeffY[p] = cy[static_cast<std::size_t>(3 * k + p)];
        }
        double sseX = 0.0;
        double sseY = 0.0;
        for (int j = 0; j < s.n; ++j) {
            const double u = s.u[static_cast<std::size_t>(j)];
            const double fx =
                info.coeffX[0] * u * u + info.coeffX[1] * u +
                info.coeffX[2];
            const double fy =
                info.coeffY[0] * u * u + info.coeffY[1] * u +
                info.coeffY[2];
            const int idx = s.a - 1 + j;
            const double rx = x[static_cast<std::size_t>(idx)] - fx;
            const double ry = y[static_cast<std::size_t>(idx)] - fy;
            sseX += rx * rx;
            sseY += ry * ry;
            if (!assigned[static_cast<std::size_t>(idx)]) {
                assigned[static_cast<std::size_t>(idx)] = true;
                fitX[static_cast<std::size_t>(idx)] = fx;
                fitY[static_cast<std::size_t>(idx)] = fy;
            }
        }
        info.sseX = sseX;
        info.sseY = sseY;
        result.segments[static_cast<std::size_t>(k)] = info;
    }

    result.sseX = 0.0;
    result.sseY = 0.0;
    for (int i = 0; i < N; ++i) {
        const double rx = x[static_cast<std::size_t>(i)] -
                          fitX[static_cast<std::size_t>(i)];
        const double ry = y[static_cast<std::size_t>(i)] -
                          fitY[static_cast<std::size_t>(i)];
        result.sseX += rx * rx;
        result.sseY += ry * ry;
    }
    result.sseTotal = result.sseX + result.sseY;
    result.rmseX = std::sqrt(result.sseX / static_cast<double>(N));
    result.rmseY = std::sqrt(result.sseY / static_cast<double>(N));
    result.rmseTotal =
        std::sqrt(result.sseTotal / static_cast<double>(N));

    result.nodeInfo.resize(static_cast<std::size_t>(K - 1));
    for (int k = 0; k + 1 < K; ++k) {
        WaldNodeInfo node;
        node.point = nodes[static_cast<std::size_t>(k + 1)];
        node.t = static_cast<double>(node.point - 1) * dt;

        const double hL = segs[static_cast<std::size_t>(k)].h;
        const double hR = segs[static_cast<std::size_t>(k + 1)].h;
        const double* lx = &cx[static_cast<std::size_t>(3 * k)];
        const double* ly = &cy[static_cast<std::size_t>(3 * k)];
        const double* rx =
            &cx[static_cast<std::size_t>(3 * (k + 1))];
        const double* ry =
            &cy[static_cast<std::size_t>(3 * (k + 1))];

        const double xL = lx[0] + lx[1] + lx[2];
        const double yL = ly[0] + ly[1] + ly[2];
        const double vxL = (2.0 * lx[0] + lx[1]) / hL;
        const double vyL = (2.0 * ly[0] + ly[1]) / hL;
        const double axL = 2.0 * lx[0] / (hL * hL);
        const double ayL = 2.0 * ly[0] / (hL * hL);

        const double xR = rx[0] - rx[1] + rx[2];
        const double yR = ry[0] - ry[1] + ry[2];
        const double vxR = (-2.0 * rx[0] + rx[1]) / hR;
        const double vyR = (-2.0 * ry[0] + ry[1]) / hR;
        const double axR = 2.0 * rx[0] / (hR * hR);
        const double ayR = 2.0 * ry[0] / (hR * hR);

        node.x = xL;
        node.y = yL;
        node.vx = vxL;
        node.vy = vyL;
        node.ax = axL;
        node.ay = ayL;
        node.dax = axR - axL;
        node.day = ayR - ayL;
        node.ep = std::hypot(xR - xL, yR - yL);
        node.ev = std::hypot(vxR - vxL, vyR - vyL);
        node.ea = std::hypot(node.dax, node.day);

        result.epMax = std::max(result.epMax, node.ep);
        result.evMax = std::max(result.evMax, node.ev);
        result.eaMax = std::max(result.eaMax, node.ea);
        result.nodeInfo[static_cast<std::size_t>(k)] = node;
    }

    return result;
}

}  // namespace

WaldDpFitResult fitWaldDp(const std::vector<double>& x,
                          const std::vector<double>& y,
                          const WaldDpConfig& config) {
    if (x.size() != y.size())
        throw std::invalid_argument("x and y must have the same size");
    const int N = static_cast<int>(x.size());
    if (N < 8)
        throw std::invalid_argument("not enough points");
    if (config.dt <= 0.0)
        throw std::invalid_argument("dt must be positive");
    if (config.windowHalfWidth < 2)
        throw std::invalid_argument("windowHalfWidth must be at least 2");
    if (config.minSegmentLength < 3)
        throw std::invalid_argument("minSegmentLength must be at least 3");
    if (config.minSegmentLength > N)
        throw std::invalid_argument("minSegmentLength exceeds point count");
    if (N <= 2 * config.windowHalfWidth)
        throw std::invalid_argument("not enough points for Wald windows");

    WaldDpFitResult result;
    result.pointCount = N;
    result.config = config;

    // 位置/速度尺度，用于 DP 中无量纲化。
    double meanX = 0.0, meanY = 0.0;
    for (int i = 0; i < N; ++i) {
        meanX += x[static_cast<std::size_t>(i)];
        meanY += y[static_cast<std::size_t>(i)];
    }
    meanX /= N;
    meanY /= N;
    double varX = 0.0, varY = 0.0;
    for (int i = 0; i < N; ++i) {
        const double dx = x[static_cast<std::size_t>(i)] - meanX;
        const double dy = y[static_cast<std::size_t>(i)] - meanY;
        varX += dx * dx;
        varY += dy * dy;
    }
    varX /= N;
    varY /= N;
    result.positionScale = std::sqrt(0.5 * (varX + varY));
    if (result.positionScale <= 0.0)
        result.positionScale = 1.0;

    double meanVX = 0.0, meanVY = 0.0;
    const int nv = N - 1;
    std::vector<double> vx(static_cast<std::size_t>(nv), 0.0);
    std::vector<double> vy(static_cast<std::size_t>(nv), 0.0);
    for (int i = 0; i < nv; ++i) {
        vx[static_cast<std::size_t>(i)] =
            (x[static_cast<std::size_t>(i + 1)] -
             x[static_cast<std::size_t>(i)]) /
            config.dt;
        vy[static_cast<std::size_t>(i)] =
            (y[static_cast<std::size_t>(i + 1)] -
             y[static_cast<std::size_t>(i)]) /
            config.dt;
        meanVX += vx[static_cast<std::size_t>(i)];
        meanVY += vy[static_cast<std::size_t>(i)];
    }
    meanVX /= nv;
    meanVY /= nv;
    double varVX = 0.0, varVY = 0.0;
    for (int i = 0; i < nv; ++i) {
        const double dx = vx[static_cast<std::size_t>(i)] - meanVX;
        const double dy = vy[static_cast<std::size_t>(i)] - meanVY;
        varVX += dx * dx;
        varVY += dy * dy;
    }
    varVX /= nv;
    varVY /= nv;
    result.velocityScale = std::sqrt(0.5 * (varVX + varVY));
    if (result.velocityScale <= 0.0)
        result.velocityScale = 1.0;

    // ---- Wald 检验(全部 6 个系数：X/Y 各 3 个) ----
    const double alphaPrime = config.alpha / static_cast<double>(N);
    result.waldCritical = chiSquareUpperQuantile(alphaPrime, 6);

    struct RawCandidate {
        int point;
        double wald;
        double p;
    };
    std::vector<RawCandidate> raw;
    const int m = config.windowHalfWidth;
    for (int center = m; center < N - m; ++center) {
        const WindowPair wp = fitWindowPair(x, y, center, m);
        if (!wp.ok)
            continue;

        double wald = 0.0;
        bool ok = true;
        for (int axis = 0; axis < 2; ++axis) {
            std::vector<double> d(3, 0.0);
            std::vector<double> V(9, 0.0);
            for (int p = 0; p < 3; ++p)
                d[static_cast<std::size_t>(p)] =
                    wp.coefL[p][axis] - wp.coefR[p][axis];
            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                    V[static_cast<std::size_t>(p) * 3 + q] =
                        wp.covL[axis][static_cast<std::size_t>(p) * 3 + q] +
                        wp.covR[axis][static_cast<std::size_t>(p) * 3 + q];
            std::vector<double> M = V;
            std::vector<double> z = d;
            try {
                solveLinearSystem(M, z, 3);
            } catch (...) {
                ok = false;
                break;
            }
            wald += dot(d, z);
        }
        if (!ok)
            continue;

        RawCandidate rc;
        rc.point = center + 1;
        rc.wald = wald;
        rc.p = chiSquareSF(wald, 6);
        raw.push_back(rc);
    }

    // NMS：按 W 降序，保留彼此距离 >= L_min 的点。
    std::vector<RawCandidate> kept;
    for (const RawCandidate& rc : raw)
        if (rc.wald > result.waldCritical)
            ++result.rawCandidateCount;
    std::vector<RawCandidate> rawSorted;
    for (const RawCandidate& rc : raw)
        if (rc.wald > result.waldCritical)
            rawSorted.push_back(rc);
    std::sort(rawSorted.begin(), rawSorted.end(),
              [](const RawCandidate& a, const RawCandidate& b) {
                  return a.wald > b.wald;
              });
    for (const RawCandidate& rc : rawSorted) {
        bool tooClose = false;
        for (const RawCandidate& k : kept)
            if (std::abs(rc.point - k.point) < config.minSegmentLength) {
                tooClose = true;
                break;
            }
        if (!tooClose)
            kept.push_back(rc);
    }
    std::sort(kept.begin(), kept.end(),
              [](const RawCandidate& a, const RawCandidate& b) {
                  return a.point < b.point;
              });
    result.candidates.resize(kept.size());
    for (std::size_t i = 0; i < kept.size(); ++i) {
        result.candidates[i].point = kept[i].point;
        result.candidates[i].wald = kept[i].wald;
        result.candidates[i].pValue = kept[i].p;
    }

    // 候选点集合，含首尾。
    std::vector<int> C;
    C.push_back(1);
    for (const WaldCandidateInfo& c : result.candidates)
        C.push_back(c.point);
    C.push_back(N);
    std::sort(C.begin(), C.end());
    C.erase(std::unique(C.begin(), C.end()), C.end());

    const int M = static_cast<int>(C.size());
    const int maxK = M - 1;
    result.dpPaths.assign(static_cast<std::size_t>(maxK + 1),
                          WaldKPathResult{});
    result.bestPaths.assign(static_cast<std::size_t>(maxK + 1),
                            WaldKPathResult{});

    // 预计算候选段。
    std::vector<std::vector<SegFit>> seg(
        static_cast<std::size_t>(M), std::vector<SegFit>(M));
    for (int i = 0; i < M; ++i) {
        for (int j = i + 1; j < M; ++j) {
            const int n = C[static_cast<std::size_t>(j)] -
                          C[static_cast<std::size_t>(i)] + 1;
            if (n < config.minSegmentLength)
                continue;
            if (config.maxSegmentLength > 0 && n > config.maxSegmentLength)
                continue;
            seg[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                fitSegment(x, y, C[static_cast<std::size_t>(i)],
                           C[static_cast<std::size_t>(j)], config.dt);
        }
    }

    // ---- DP：对每个 K 求软约束最优路径 ----
    const double sp2 = result.positionScale * result.positionScale;
    const double sv2 = result.velocityScale * result.velocityScale;

    std::vector<std::vector<std::vector<double>>> dp(
        static_cast<std::size_t>(maxK + 1),
        std::vector<std::vector<double>>(
            static_cast<std::size_t>(M),
            std::vector<double>(static_cast<std::size_t>(M), kInf)));
    std::vector<std::vector<std::vector<int>>> prev(
        static_cast<std::size_t>(maxK + 1),
        std::vector<std::vector<int>>(
            static_cast<std::size_t>(M),
            std::vector<int>(static_cast<std::size_t>(M), -1)));

    for (int j = 1; j < M; ++j) {
        const SegFit& s = seg[0][static_cast<std::size_t>(j)];
        if (!s.valid)
            continue;
        dp[1][0][static_cast<std::size_t>(j)] = s.sseX + s.sseY;
        prev[1][0][static_cast<std::size_t>(j)] = -2;
    }

    for (int k = 2; k <= maxK; ++k) {
        for (int i = 1; i < M - 1; ++i) {
            for (int j = i + 1; j < M; ++j) {
                const SegFit& cur = seg[static_cast<std::size_t>(i)]
                                       [static_cast<std::size_t>(j)];
                if (!cur.valid)
                    continue;
                double best = kInf;
                int bestH = -1;
                for (int h = 0; h < i; ++h) {
                    const double base =
                        dp[static_cast<std::size_t>(k - 1)]
                          [static_cast<std::size_t>(h)]
                          [static_cast<std::size_t>(i)];
                    if (!std::isfinite(base))
                        continue;
                    const SegFit& prv = seg[static_cast<std::size_t>(h)]
                                           [static_cast<std::size_t>(i)];
                    if (!prv.valid)
                        continue;
                    const double dpx = cur.pL[0] - prv.pR[0];
                    const double dpy = cur.pL[1] - prv.pR[1];
                    const double dvx = cur.vL[0] - prv.vR[0];
                    const double dvy = cur.vL[1] - prv.vR[1];
                    const double pen =
                        config.muP * (dpx * dpx + dpy * dpy) / sp2 +
                        config.muV * (dvx * dvx + dvy * dvy) / sv2;
                    const double val =
                        base + cur.sseX + cur.sseY + pen;
                    if (val < best) {
                        best = val;
                        bestH = h;
                    }
                }
                dp[static_cast<std::size_t>(k)]
                  [static_cast<std::size_t>(i)]
                  [static_cast<std::size_t>(j)] = best;
                prev[static_cast<std::size_t>(k)]
                    [static_cast<std::size_t>(i)]
                    [static_cast<std::size_t>(j)] = bestH;
            }
        }
    }

    for (int k = 1; k <= maxK; ++k) {
        double best = kInf;
        int bestI = -1;
        for (int i = 0; i < M - 1; ++i) {
            const double v = dp[static_cast<std::size_t>(k)]
                               [static_cast<std::size_t>(i)]
                               [static_cast<std::size_t>(M - 1)];
            if (v < best) {
                best = v;
                bestI = i;
            }
        }
        if (bestI < 0 || !std::isfinite(best))
            continue;

        std::vector<int> nodes;
        nodes.push_back(C[static_cast<std::size_t>(M - 1)]);
        int i = bestI;
        int j = M - 1;
        for (int kk = k; kk >= 1; --kk) {
            nodes.push_back(C[static_cast<std::size_t>(i)]);
            const int h = prev[static_cast<std::size_t>(kk)]
                              [static_cast<std::size_t>(i)]
                              [static_cast<std::size_t>(j)];
            j = i;
            i = h;
        }
        std::reverse(nodes.begin(), nodes.end());

        WaldKPathResult path = fitKkt(x, y, nodes, config.dt);
        path.dpCost = best;
        result.dpPaths[static_cast<std::size_t>(k)] = path;
    }

    // ---- 枚举所有合法候选路径，做 KKT 复核，取每个 K 的 KKT-SSE 最小 ----
    int pathCount = 0;
    const int kPathLimit = 50000;
    std::vector<int> current;
    current.push_back(C[0]);

    std::function<void(int)> dfs = [&](int idx) {
        if (pathCount >= kPathLimit) {
            result.enumerationTruncated = true;
            return;
        }
        if (idx == M - 1) {
            if (current.size() < 2)
                return;
            ++pathCount;
            const WaldKPathResult path = fitKkt(x, y, current, config.dt);
            WaldKPathResult& best =
                result.bestPaths[static_cast<std::size_t>(path.K)];
            if (best.K == 0 || path.sseTotal < best.sseTotal)
                best = path;
            return;
        }
        for (int j = idx + 1; j < M; ++j) {
            if (!seg[static_cast<std::size_t>(idx)]
                    [static_cast<std::size_t>(j)]
                    .valid)
                continue;
            current.push_back(C[static_cast<std::size_t>(j)]);
            dfs(j);
            current.pop_back();
        }
    };
    dfs(0);

    // ---- 5% 边际收益规则 ----
    int chosen = 0;
    for (int k = 1; k <= maxK; ++k) {
        if (result.bestPaths[static_cast<std::size_t>(k)].K == 0)
            continue;
        int next = k + 1;
        while (next <= maxK &&
               result.bestPaths[static_cast<std::size_t>(next)].K == 0)
            ++next;
        if (next > maxK) {
            chosen = k;
            break;
        }
        const double cur =
            result.bestPaths[static_cast<std::size_t>(k)].sseTotal;
        const double nxt =
            result.bestPaths[static_cast<std::size_t>(next)].sseTotal;
        const double rel = (cur > 0.0) ? (cur - nxt) / cur : 0.0;
        if (rel < 0.05) {
            chosen = k;
            break;
        }
    }
    if (chosen == 0) {
        for (int k = maxK; k >= 1; --k)
            if (result.bestPaths[static_cast<std::size_t>(k)].K != 0) {
                chosen = k;
                break;
            }
    }
    result.chosenK = chosen;
    if (chosen >= 1 && chosen <= maxK)
        result.chosen = result.bestPaths[static_cast<std::size_t>(chosen)];

    return result;
}
