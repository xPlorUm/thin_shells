//
// Created by hugues on 09/12/24.
//

#ifndef THINSHELLS_SOLVER_H
#define THINSHELLS_SOLVER_H

#include <Eigen/Core>
#include <Eigen/Sparse>


class DiscreteShell;

class Solver {
public:
    Solver() {}
    Solver(double beta, double gamma, double dt, const Eigen::SparseMatrix<double> *M): m_beta(beta), m_gamma(gamma), m_dt(dt){
        // Initialize the mass matrix, identity matrix of size 3nx3n from the original mass matrix
        m_M_extended = new Eigen::SparseMatrix<double>(3 * M->rows(), 3 * M->cols());
        for (int i = 0; i < M->rows(); i++) {
            for (int j = 0; j < M->cols(); j++) {
                m_M_extended->coeffRef(3 * i, 3 * j) = M->coeff(i, j);
                m_M_extended->coeffRef(3 * i + 1, 3 * j + 1) = M->coeff(i, j);
                m_M_extended->coeffRef(3 * i + 2, 3 * j + 2) = M->coeff(i, j);
            }
        }
    };

    void solve(Eigen::MatrixXd *Position_i, Eigen::MatrixXd *Velocity_i, Eigen::MatrixXd *Acceleration_i);

    /**
     * Sets the discrete shell
     */
    void setDiscreteShell(DiscreteShell *discreteShell) {
        m_discreteshell = discreteShell;
    }

private :
    double m_beta;
    double m_gamma;
    double m_dt;
    // Extended mass matrix (of size 3n x 3n)
    Eigen::SparseMatrix<double> *m_M_extended;
    DiscreteShell *m_discreteshell;

    void
    linearSolve(const Eigen::SparseMatrix<double> &systemMatrix, const Eigen::VectorXd &rhs, Eigen::VectorXd &delta_x);
};


#endif //THINSHELLS_SOLVER_H
