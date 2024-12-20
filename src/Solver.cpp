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
    // Same but for stretching
    Eigen::SparseMatrix<double> HessianStretching = Eigen::SparseMatrix<double>(m_discreteshell->N_VERTICES * 3, m_discreteshell->N_VERTICES * 3);
    // Same, but for area preservation
    Eigen::SparseMatrix<double> HessianAreaPreserving = Eigen::SparseMatrix<double>(m_discreteshell->N_VERTICES * 3, m_discreteshell->N_VERTICES * 3);

    // ---- RESIDUAL FUNCTION ----
    auto R = [this, &HessianBending, &HessianStretching, &HessianAreaPreserving](
                    const Eigen::VectorXd &u_new,
                    const Eigen::VectorXd &u_n,
                    const Eigen::VectorXd &v_n,
                    const Eigen::VectorXd &a_n,
                    const Eigen::VectorXd &f_ext_flatten) -> Eigen::VectorXd {
        // Step 1: Compute acceleration at n+1 using Newmark's formula
        Eigen::VectorXd a_new;
        if (SIMULATION_NEWMARK_BETA == 0) {
            a_new = a_n;
        } else {
            a_new = (u_new - u_n - m_dt * v_n - (m_dt * m_dt) * (0.5 - m_beta) * a_n) / (m_beta * m_dt * m_dt);
        }
        // Step 2: Compute velocity at n+1 using Newmark's formula
        Eigen::VectorXd v_new = v_n + m_dt * ((1.0 - m_gamma) * a_n + m_gamma * a_new);
        // Step 3: Map u_new to Position matrix for force computation
        // Step 4: Compute internal forces based on the updated positions
        Eigen::MatrixXd f_int_matrix = Eigen::MatrixX3d::Zero(m_discreteshell->N_VERTICES, 3);
        // Convert back u_n to Eigen::MatrixXd and takes its pointer
        Eigen::MatrixXd u_new_matrix = deflatten_vector(u_new);
        if (ENABLE_STRETCHING_FORCES)
            m_discreteshell->addStretchingForcesAndHessianTo_AD(f_int_matrix, HessianStretching, &u_new_matrix);
        if (ENABLE_AREA_PRESERVATION_FORCES)
            m_discreteshell->addAreaPreserationForcesAndHessianTo(f_int_matrix, HessianAreaPreserving, &u_new_matrix);
        if (ENABLE_BENDING_FORCES)
            m_discreteshell->addBendingForcesAndHessianTo(f_int_matrix, HessianBending, &u_new_matrix);
        Eigen::VectorXd f_int = flatten_matrix(f_int_matrix);
        // Step 5: Compute residual: R = M * a_new + C * v_new + f_int - f_ext
        // As of now, damping matrix is a simple identity matrix SET TO ZERO
        Eigen::VectorXd residual = *m_M_extended * a_new + f_int - f_ext_flatten;
        return residual;
    };

    const Eigen::VectorXd Position_i_flatten = flatten_matrix(*Position_i);
    const Eigen::VectorXd Velocity_i_flatten = flatten_matrix(*Velocity_i);
    const Eigen::VectorXd Acceleration_i_flatten = flatten_matrix(*Acceleration_i);
    // Initial guess. u_new is the iterative guessed variable.
    Eigen::VectorXd u_new = Position_i_flatten;

    Eigen::MatrixXd f_ext_matrix = Eigen::MatrixXd::Zero(m_discreteshell->N_VERTICES, 3);
    m_discreteshell->add_F_ext(f_ext_matrix, -1);
    Eigen::VectorXd f_ext_flatten = flatten_matrix(f_ext_matrix);
    // ---- SYSTEM MATRIX FUNCTION ----
    auto S = [this, C, &HessianBending, &HessianStretching, &HessianAreaPreserving](const Eigen::VectorXd &u_new) {
        std::vector<Eigen::Triplet<double>> triplets = std::vector<Eigen::Triplet<double>>();
        Eigen::SparseMatrix<double> systemMatrix = Eigen::SparseMatrix<double>(m_discreteshell->N_VERTICES * 3, m_discreteshell->N_VERTICES * 3);
        Eigen::MatrixXd u_new_matrix = deflatten_vector(u_new);
        if (ENABLE_STRETCHING_FORCES) {
            systemMatrix += HessianStretching;
        }
        if (ENABLE_AREA_PRESERVATION_FORCES) {
            systemMatrix += HessianAreaPreserving;
        }
        // Bending hessian is computed in the residual function. See comment of the variable
        if (ENABLE_BENDING_FORCES) {
            systemMatrix += HessianBending;
        }
        systemMatrix += (*m_M_extended / (m_beta * (m_dt * m_dt))); //+ C * (m_gamma / (m_beta * m_dt));
        return systemMatrix;
    };

    // Reference : https://chatgpt.com/c/67585038-9a24-800f-a9cd-a44d464ad723
    std::cout << "-----------------------------" << std::endl;
    for(int iter = 0; iter < MAX_ITERATIONS_NEWTON; iter++) {
        Eigen::VectorXd residual = R(u_new, Position_i_flatten, Velocity_i_flatten, Acceleration_i_flatten, f_ext_flatten);
        double residual_squared_norm = residual.squaredNorm();
        if (residual_squared_norm < TOLERANCE_NEWTON) {
            Eigen::MatrixXd sol = deflatten_vector(u_new);
            Position_solution = sol;
            std::cout << "Converged in " << iter << " iterations,  Δx (correction) : " << (u_new - Position_i_flatten).norm() << std::endl;
            return true;
        }
        Eigen::SparseMatrix<double> systemMatrix = S(u_new);
        // Sove delta_x = K^-1 * -R (The system matrix is the jacobian of R, so J = K)
        Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(Position_i->rows() * 3);
        if (!linearSolve(systemMatrix, residual, delta_x)) {
            std::cerr << "Linear solve failed!" << std::endl;
            break;
        }

#ifdef USE_AMIJO
        const double alpha_initial = 4.0;      // Initial step size
        const double alpha_min = 0.25;         // Minimum step size
        const double rho = 0.8;                // Step size reduction factor
        const double c = 0.1;                  // Armijo condition constant
        // Initialize step size
        double alpha = alpha_initial;
        Eigen::VectorXd u_trial;
        // Backtracking loop
        int i = 10;
        // ref : https://math.stackexchange.com/questions/972890/how-to-find-the-gradient-of-norm-square
        Eigen::VectorXd gradient_squared_norm_residual = 2 * systemMatrix * residual;
        while (alpha >= alpha_min || i--) {
            u_trial = u_new - alpha * delta_x;
            // Compute new residual with u_trial
            Eigen::VectorXd residual_trial = R(u_trial, Position_i_flatten, Velocity_i_flatten, Acceleration_i_flatten, f_ext_flatten);
            double lhs = residual_squared_norm - residual_trial.squaredNorm();
            double rhs = alpha * c * gradient_squared_norm_residual.dot(delta_x);
            // PRINT_VAR(gradient_squared_norm_residual.dot(delta_x))
            // PRINT_VAR(lhs)
            // PRINT_VAR(rhs)
            if (lhs >= 0) {
                break; // Sufficient decrease achieved
            }
            alpha *= rho; // Reduce step size
        }

        if (alpha < alpha_min) {
            std::cerr << "Armijo condition not satisfied. Line search failed at iteration " << iter << "." << std::endl;
            alpha = 1;
        }
#endif
        // Update the guess with the accepted step
#ifndef USE_AMIJO
        const double alpha = 1;
        // delta_x(0) = 0;
        // delta_x(1) = 0;
        // delta_x(2) = 0;
        u_new -=  alpha * delta_x;
        std::cout << "Iteration " << iter << ": Step size = " << alpha << ", Update norm = " << (alpha * delta_x).norm()
                  << " residual squared norm : " << residual_squared_norm << std::endl;
#endif
    }
    std::cerr << "Did not converge!" << std::endl;
    return false;
}

bool Solver::linearSolve(const Eigen::SparseMatrix<double> &systemMatrix,
                         const Eigen::VectorXd &rhs,
                         Eigen::VectorXd &delta_x) {
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