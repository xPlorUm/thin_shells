//
// Created by hugues on 09/12/24.
//

#include <Eigen/SparseCore>
#include <iostream>
#include "Solver.h"
#include "DiscreteShell.h"
#include "common.h"
#include "constants.h"

bool Solver::solve(Eigen::MatrixXd &Position_solution, Eigen::MatrixXd *Velocity_i, Eigen::MatrixXd *Acceleration_i,
                   Eigen::MatrixXd *Position_i) {
    // TODO : pass those as constants, or handle better the flattened versions
    Eigen::MatrixXd C = 0.5 * Eigen::MatrixXd::Identity(m_discreteshell->N_VERTICES * 3, m_discreteshell->N_VERTICES * 3);

    // THe hessian of the bending energy. This matrix is made global to avoid repetitive calls from
    // addBendingForcesAndHessianTo.
    Eigen::SparseMatrix<double> HessianBending = Eigen::SparseMatrix<double>(m_discreteshell->N_VERTICES * 3, m_discreteshell->N_VERTICES * 3);

    // ---- RESIDUAL FUNCTION ----
    auto R = [this, &HessianBending](
                    const Eigen::VectorXd &u_new,
                    const Eigen::VectorXd &u_n,
                    const Eigen::VectorXd &v_n,
                    const Eigen::VectorXd &a_n,
                    const Eigen::VectorXd &f_ext_flatten) -> Eigen::VectorXd {
        // Step 1: Compute acceleration at n+1 using Newmark's formula
        // a_new = \ddot{x}_{i + 1}
        Eigen::VectorXd a_new =
                (u_new - u_n - m_dt * v_n - (m_dt * m_dt) * (0.5 - m_beta) * a_n) / (m_beta * m_dt * m_dt);
        // Step 2: Compute velocity at n+1 using Newmark's formula
        Eigen::VectorXd v_new = v_n + m_dt * ((1.0 - m_gamma) * a_n + m_gamma * a_new);
        // Step 3: Map u_new to Position matrix for force computation
        // Step 4: Compute internal forces based on the updated positions
        Eigen::MatrixXd f_int_matrix = Eigen::MatrixX3d::Zero(m_discreteshell->N_VERTICES, 3);
        // Convert back u_n to Eigen::MatrixXd and takes its pointer
        Eigen::MatrixXd u_new_matrix = deflatten_vector(u_new);
        m_discreteshell->addStrechingForcesTo(f_int_matrix, &u_new_matrix); // Use updated positions
        m_discreteshell->addBendingForcesAndHessianTo(f_int_matrix, HessianBending, &u_new_matrix);
        Eigen::VectorXd f_int = flatten_matrix(f_int_matrix);
        // Step 5: Compute residual: R = M * a_new + C * v_new + f_int - f_ext
        // As of now, damping matrix is a simple identity matrix SET TO ZERO
        Eigen::VectorXd residual = *m_M_extended * a_new - f_int - f_ext_flatten;
        return residual;
    };

    Eigen::VectorXd Position_flatten = flatten_matrix(*Position_i);
    Eigen::VectorXd Velocity_flatten = flatten_matrix(*Velocity_i);
    Eigen::VectorXd Acceleration_flatten = flatten_matrix(*Acceleration_i);
    Eigen::VectorXd u_new = Position_flatten;

    Eigen::MatrixXd f_ext_matrix = Eigen::MatrixXd::Zero(m_discreteshell->N_VERTICES, 3);
    // m_discreteshell->add_F_ext(f_ext_matrix);
    Eigen::VectorXd f_ext_flatten = flatten_matrix(f_ext_matrix);
    // ---- SYSTEM MATRIX FUNCTION ----
    auto S = [this, C, Position_i, &HessianBending](Eigen::VectorXd &Position, Eigen::VectorXd &Velocity, Eigen::VectorXd &Acceleration) {
        std::vector<Eigen::Triplet<double>> triplets = std::vector<Eigen::Triplet<double>>();
        Eigen::SparseMatrix<double> systemMatrix = Eigen::SparseMatrix<double>(Position_i->rows() * 3, Position_i->rows() * 3);

        m_discreteshell->addStretchingHessianTo(triplets, Position_i);
        systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
        // Bending hessian is computed in the residual function. See comment of the variable
        systemMatrix += HessianBending;

        systemMatrix += (*m_M_extended / (m_beta * (m_dt * m_dt))); //+ C * (m_gamma / (m_beta * m_dt));
        return systemMatrix;
    };

    // Reference : https://chatgpt.com/c/67585038-9a24-800f-a9cd-a44d464ad723
    std::cout << "-----------------------------" << std::endl;
    for(int iter = 0; iter < MAX_ITERATIONS_NEWTON; iter++) {
        Eigen::VectorXd residual = R(u_new, Position_flatten, Velocity_flatten, Acceleration_flatten, f_ext_flatten);
        double residual_norm = residual.norm();
        if (residual_norm < TOLERANCE_NEWTON) {
            Eigen::MatrixXd sol = deflatten_vector(u_new);
            Position_solution = sol;
            std::cout << "Converged in " << iter << " iterations,  Δx : " << (u_new - Position_flatten).norm() << std::endl;
            return true;
        }
        Eigen::SparseMatrix<double> systemMatrix = S(u_new, Velocity_flatten, Acceleration_flatten);
        // Sove delta_x = K^-1 * -R (The system matrix is the jacobian of R, so J = K)
        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(Position_i->rows() * 3);
        if (!linearSolve(systemMatrix, residual, delta_x)) {
            std::cerr << "Linear solve failed!" << std::endl;
            break;
        }
        u_new -= delta_x;
        std::cout << "Iteration " << iter << ": delta_x norm = " << delta_x.norm() <<  " residual=" << residual_norm << std::endl;
    }
    std::cerr << "Did not converge!" << std::endl;
    return false;
}



bool Solver::linearSolve(const Eigen::SparseMatrix<double> &systemMatrix,
                         const Eigen::VectorXd &rhs,
                         Eigen::VectorXd &delta_x) {
    // We can use a direct sparse solver from Eigen, for example SparseLU or SparseCholesky
    // depending on the system properties. Here we'll use SparseLU for generality.

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(systemMatrix);
    solver.factorize(systemMatrix);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed!" << std::endl;
        return false;
    }

    delta_x = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
        return false;
    }
    return true;
}
