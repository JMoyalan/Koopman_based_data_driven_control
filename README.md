# Koopman operator-based data-driven optimal control
The code in this respository is concerned with the data-driven optimal control of nonlinear systems. We present a
convex formulation of the optimal control problem with a discounted cost function. We consider
optimal control problems with both positive and negative discount factors. The convex approach
relies on lifting nonlinear system dynamics in the space of densities using the linear Perron–
Frobenius operator. This lifting leads to an infinite-dimensional convex optimization formulation of
the optimal control problem. The data-driven approximation of the optimization problem relies on
the approximation of the Koopman operator and its dual: the Perron–Frobenius operator, using a
polynomial basis function. We write the approximate finite-dimensional optimization problem as a
polynomial optimization which is then solved efficiently using a sum-of-squares-based optimization
framework.
