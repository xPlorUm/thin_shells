//
// Created by hugues on 09/12/24.
//

#include <Eigen/SparseCore>
#include <iostream>
#include "Solver.h"
#include "DiscreteShell.h"

void Solver::solve(Eigen::MatrixXd *Position_i, Eigen::MatrixXd *Velocity_i, Eigen::MatrixXd *Acceleration_i) {
    // TODO : pass those as constants, or handle better the flattened versions
    Eigen::MatrixXd C = 0.5 * Eigen::MatrixXd::Identity(Position_i->rows() * 3, Position_i->rows() * 3);
    auto R = [this, C](Eigen::VectorXd &Position, Eigen::VectorXd &Velocity, Eigen::VectorXd &Acceleration) {
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(Position.rows());
        // reference : eq (4) from Chang paper
        residual += (*m_M_extended / (m_beta * (m_dt * m_dt))) *
                    (Position + m_dt * Velocity + (0.5 - m_beta) * m_dt * m_dt * Acceleration);
        residual += -C * (Velocity + (1 - m_gamma) * (m_dt) * Acceleration);
        residual += C * (m_gamma / (m_beta * m_dt)) * (Position + m_dt * Velocity);
        residual += (0.5 - m_beta) * (m_dt * m_dt) * Acceleration;
        return residual;
    };

    // Create new vectors with flattened values, copied
    Eigen::VectorXd Position_flatten = Eigen::Map<const Eigen::VectorXd>(Position_i->data(), Position_i->rows() * 3);
    Eigen::VectorXd Velocity_flatten = Eigen::Map<const Eigen::VectorXd>(Velocity_i->data(), Velocity_i->rows() * 3);
    Eigen::VectorXd Acceleration_flatten = Eigen::Map<const Eigen::VectorXd>(Acceleration_i->data(), Acceleration_i->rows() * 3);
    Eigen::VectorXd residual = R(Position_flatten, Velocity_flatten, Acceleration_flatten);
    // Need to add forces
    Eigen::MatrixX3d f_i = Eigen::MatrixX3d::Zero(Position_i->rows(), 3);
    // SEGFAULT HERE
    m_discreteshell->addStrechingForcesTo(f_i, Position_i);
    Eigen::VectorXd f_i_flatten = Eigen::Map<const Eigen::VectorXd>(f_i.data(), f_i.rows() * 3);
    residual += f_i_flatten;

    // Compute System matrix of the system
    auto S = [this, C](Eigen::VectorXd &Position, Eigen::VectorXd &Velocity, Eigen::VectorXd &Acceleration) {
        Eigen::SparseMatrix<double> K = Eigen::SparseMatrix<double>(Position.rows(), Position.rows());
        K = (*m_M_extended / (m_beta * (m_dt * m_dt))) + C * (m_gamma / (m_beta * m_dt));
        return K;
    };

    std::vector<Eigen::Triplet<double>> triplets = std::vector<Eigen::Triplet<double>>();
    Eigen::SparseMatrix<double> systemMatrix = Eigen::SparseMatrix<double>(Position_i->rows() * 3, Position_i->rows() * 3);
    m_discreteshell->addStretchingHessianTo(triplets, Position_i);
    systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
    systemMatrix += S(Position_flatten, Velocity_flatten, Acceleration_flatten);
    // K delta_x = -R : solving that system
    Eigen::VectorXd delta_x = Eigen::VectorXd::Zero(Position_i->rows() * 3);
    linearSolve(systemMatrix, -residual, delta_x);
}



void Solver::linearSolve(const Eigen::SparseMatrix<double> &systemMatrix,
                         const Eigen::VectorXd &rhs,
                         Eigen::VectorXd &delta_x) {
    // We can use a direct sparse solver from Eigen, for example SparseLU or SparseCholesky
    // depending on the system properties. Here we'll use SparseLU for generality.

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(systemMatrix);
    solver.factorize(systemMatrix);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Factorization failed!" << std::endl;
        return;
    }

    delta_x = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
    }
}
