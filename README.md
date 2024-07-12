# Koopman operator-based data-driven optimal control
The code in this repository is concerned with the data-driven optimal control of nonlinear systems. We present a
convex formulation of the optimal control problem with a discounted cost function. We consider
optimal control problems with both positive and negative discount factors. The convex approach
relies on lifting nonlinear system dynamics in the space of densities using the linear Perron–
Frobenius operator. This lifting leads to an infinite-dimensional convex optimization formulation of
the optimal control problem. The data-driven approximation of the optimization problem relies on
the approximation of the Koopman operator and its dual: the Perron–Frobenius operator, using a
polynomial basis function. We write the approximate finite-dimensional optimization problem as a
polynomial optimization which is then solved efficiently using a sum-of-squares-based optimization
framework.

# Koopman and Perron-Frobenius (P-F) operators
 Koopman operator is a linear operator in the function space. If **f**(**x**) is the vector field of the dynamics and the $`\psi`$(**x**) is an observable function or lifting function, the Koopman operator is defined w.r.t to the **f**(**x**) and is given as composition of the function $`\psi`$(**x**) with **f**(**x**). The linearity is a simple consequence of the composition $`\psi \ocircle`$ **f**(**x**) being a linear function  of the lifting function. Similarly, the P-F operator is a linear operator and gives you the evolution of the density of system trajectories in the function space. The Koopman and P-F operators are implemented for discrete-time systems. For continuous-time systems, we use Koopman and P-F generators.

