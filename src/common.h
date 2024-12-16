//
// Created by hugues on 18/11/24.
//

#ifndef THINSHELLS_COMMON_H
#define THINSHELLS_COMMON_H

#include <Eigen/Core>
#include <iostream>

#define PRINT_SHAPE(matrix) std::cout << #matrix << " : " << matrix.rows() << " " << matrix.cols() << std::endl;
#define CHECK_VECTOR_SIZE(vec, size) assert((vec).size() == (size))
#define PRINT_FLAT_VECTOR(vec) std::cout << #vec << "\n" << deflatten_vector(vec) << std::endl;
#define PRINT_MATRIX(vec) std::cout << #vec << "\n" << vec << std::endl;
// I hate typing std::cout and std::endl
#define LOG(x) std::cout << x << std::endl;


inline void pretty_print_vector(const Eigen::VectorXd& v) {
    std::cout << "Vector : " << std::endl;
    for (int i = 0; i < v.size(); i++) {
        std::cout << v(i) << " ";
    }
    std::cout << std::endl;
}

inline void pretty_print_vectorX3d(const Eigen::MatrixX3d& v) {
    std::cout << "Vector : " << std::endl;
    for (int i = 0; i < v.rows(); i++) {
        std::cout << v(i, 0) << " " << v(i, 1) << " " << v(i, 2) << std::endl;
    }
    std::cout << std::endl;
}


inline void extendMatrix(const Eigen::SparseMatrix<double> *M, Eigen::SparseMatrix<double> &M_extended) {
    for (int i = 0; i < M->rows(); i++) {
        for (int j = 0; j < M->cols(); j++) {
            M_extended.coeffRef(3 * i, 3 * j) = M->coeff(i, j);
            M_extended.coeffRef(3 * i + 1, 3 * j + 1) = M->coeff(i, j);
            M_extended.coeffRef(3 * i + 2, 3 * j + 2) = M->coeff(i, j);
        }
    }
}

/**
 * Returns a flattened version of the matrix, where each row is concatenated to the previous one.
 *
 * Eg :
 * 1 2 3
 * 4 5 6
 * 7 8 9
 * Will be
 * 1 2 3 4 5 6 7 8 9
 */
inline Eigen::VectorXd flatten_matrix(Eigen::MatrixXd &m) {
    Eigen::MatrixXd m_transposed = m.transpose();
    return Eigen::Map<Eigen::VectorXd>(m_transposed.data(), m_transposed.cols() * m_transposed.rows());
}

/**
 * "deflattens" a vector into a matrix, where each row is a part of the vector.
 *
 * Eg :
 * 1 2 3 4 5 6 7 8 9
 * will be
 * 1 2 3
 * 4 5 6
 * 7 8 9
 */
inline Eigen::MatrixXd deflatten_vector(const Eigen::VectorXd &v) {
    // TODO : this is garbage, find a better way
    Eigen::MatrixXd m = Eigen::MatrixXd(v.size() / 3, 3);
    for (int i = 0; i < v.size(); i += 3) {
        m.row(i / 3) = v.segment<3>(i).transpose();
    }
    return m;
}

/**
 * What you would expect it to do.
 */
template <typename Derived>
inline auto computeArea(const Eigen::MatrixBase<Derived> &v0,
                        const Eigen::MatrixBase<Derived> &v1,
                        const Eigen::MatrixBase<Derived> &v2)
-> typename Derived::Scalar {
    using T = typename Derived::Scalar;
    Eigen::Matrix<T,3,1> e1 = (v1 - v0).transpose();
    Eigen::Matrix<T,3,1> e2 = (v2 - v0).transpose();
    return T(0.5) * (e1.cross(e2)).norm();
}

#endif //THINSHELLS_COMMON_H
