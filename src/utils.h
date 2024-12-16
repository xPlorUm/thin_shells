//
// Created by hugues on 11/12/24.
//

#pragma once

#define PRINT_SHAPE(matrix) std::cout << #matrix << " : " << matrix.rows() << " " << matrix.cols() << std::endl;
#define CHECK_VECTOR_SIZE(vec, size) assert((vec).size() == (size))
#define CHECK_MATRIX_SIZE(matrix, rows, cols) assert((matrix).rows() == (rows) && (matrix).cols() == (cols))
