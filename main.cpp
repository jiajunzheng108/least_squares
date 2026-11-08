#include "pipeline/LineFitPipeline.h"
#include "pipeline/CurveFitPipeline.h"
#include "pipeline/TimeQuadFitPipeline.h"
#include "pipeline/PiecewiseQuadFitPipeline.h"
#include "pipeline/WaldDpFitPipeline.h"

#include <exception>
#include <iostream>

int main() {
    try {
        // runLineFitPipeline("label.xlsx", "L1", "X0", "Y0");
        // runCurveFitPipeline("lab110.xlsx", "L1", "X0", "Y0");
        // runTimeQuadFitPipeline("lab110.xlsx", "L1", "X0", "Y0", 0.001);
        // runPiecewiseQuadFitPipeline("lab110.xlsx", "L1", "X0", "Y0", 0.001,
        //                             10, 9, 12);
        runWaldDpFitPipeline("lab110.xlsx", "L1", "X0", "Y0", 0.001,
                             7, 6, 0, 0.05, 1.0e4, 1.0e4);
    } catch (const std::exception& e) {
        std::cerr << "失败: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
