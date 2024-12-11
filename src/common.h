//
// Created by hugues on 18/11/24.
//

#ifndef THINSHELLS_COMMON_H
#define THINSHELLS_COMMON_H

#include <Eigen/Core>
#include <iostream>

void pretty_print_vector(const Eigen::VectorXd& v) {
    std::cout << "Vector : " << std::endl;
    for (int i = 0; i < v.size(); i++) {
        std::cout << v(i) << " ";
    }
    std::cout << std::endl;
}

void pretty_print_vectorX3d(const Eigen::MatrixX3d& v) {
    std::cout << "Vector : " << std::endl;
    for (int i = 0; i < v.rows(); i++) {
        std::cout << v(i, 0) << " " << v(i, 1) << " " << v(i, 2) << std::endl;
    }
    std::cout << std::endl;
}


void extendMatrix(const Eigen::SparseMatrix<double> *M, Eigen::SparseMatrix<double> &M_extended) {
    for (int i = 0; i < M->rows(); i++) {
        for (int j = 0; j < M->cols(); j++) {
            M_extended.coeffRef(3 * i, 3 * j) = M->coeff(i, j);
            M_extended.coeffRef(3 * i + 1, 3 * j + 1) = M->coeff(i, j);
            M_extended.coeffRef(3 * i + 2, 3 * j + 2) = M->coeff(i, j);
        }
    }
}

#endif //THINSHELLS_COMMON_H
