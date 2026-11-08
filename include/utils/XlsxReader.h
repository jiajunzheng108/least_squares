#ifndef XLSX_READER_H
#define XLSX_READER_H

#include <string>
#include <vector>

// 读取 xlsx 中指定工作表的点数据。
// 先在该工作表中按表头名 xHeader / yHeader 定位两列，再读取其下方所有数据行，
// 分别作为 x / y 坐标。
void readPointsFromXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::string& xHeader,
                        const std::string& yHeader,
                        std::vector<double>& x,
                        std::vector<double>& y);

// 按列名向指定工作表写数据：
// - overwrite=false：在表头行按 headers 中的名称定位列，找不到的列在右侧新建；
//   数据追加到已有数据末尾(多列共用同一公共起始行以保持对齐)。
// - overwrite=true：清空目标工作表，重新写表头(第 1 行)和数据(第 2 行起)。
// - startRow>0：在指定行(startRow)写表头，数据从 startRow+1 行开始，
//   保留该工作表其它单元格；此时忽略 overwrite。
// 写入时 columns 中的 NaN 会被跳过，对应单元格留空。
// 保留其它工作表，原地修改 xlsx 文件。headers 与 columns 数量须一致，
// 各数据列长度须一致。
void writeColumnsToXlsx(const std::string& path,
                        const std::string& sheetName,
                        const std::vector<std::string>& headers,
                        const std::vector<std::vector<double>>& columns,
                        bool overwrite = false,
                        int startRow = 0);

#endif
