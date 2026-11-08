// 联合优化和分布优化操作封装函数
#include "pipeline/LineFitPipeline.h"

#include "utils/XlsxReader.h"
#include "utils/SortUtils.h"

#include <iostream>
#include <stdexcept>
#include <vector>

LineFitSummary runLineFitPipeline(const std::string& inputPath,
                                  const std::string& inputSheet,
                                  const std::string& inputXHeader,
                                  const std::string& inputYHeader,
                                  const std::string& outputSheet,
                                  const std::vector<std::string>& outputHeaders) {
    std::vector<double> x, y;
    readPointsFromXlsx(inputPath, inputSheet, inputXHeader, inputYHeader, x, y);
    
    // 按 x 升序排序，y 列随 x 列一起重排
    // sortutils::sortPoints(x, &y);

    LineFitSummary summary;
    summary.pointCount = x.size();
    summary.uniform = uniformLineFit(x, y);
    // summary.orthogonal = orthogonalLineFit(x, y);
    // summary.orthogonal = orthogonalFitResidual(x, y);

    std::cout << "=== LineFit: " << inputPath << " / " << inputSheet << " ("
              << inputXHeader << ", " << inputYHeader << ") ===\n";
    std::cout << "points = " << summary.pointCount << "\n";
    std::cout << "[uniform]    k=" << summary.uniform.k
              << "  b=" << summary.uniform.b << "  Xb=" << summary.uniform.Xb
              << "  d=" << summary.uniform.d << "  sse=" << summary.uniform.sse
              << "\n";
    for (std::size_t i = 0; i < summary.uniform.xFit.size(); ++i) {
        std::cout << "  [" << i << "] xFit=" << summary.uniform.xFit[i]
                  << "  yFit=" << summary.uniform.yFit[i]
                  << "  residual=" << summary.uniform.residual[i] << "\n";
    }          

    // const std::vector<std::vector<double>> columns = {
    // summary.uniform.xFit, summary.uniform.yFit};
    // writeColumnsToXlsx(inputPath, outputSheet, outputHeaders, columns);

    return summary;
}
