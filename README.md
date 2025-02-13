# Schur-Q-Plethysm
Computes Schur's Q-functions as polynomials in the $q_n$'s, and computes the plethysm $Q_\lambda\circ Q_\mu$. Also gives conversions between Schur's $Q$-functions and power sums.

# How to Use

We recommend the [online Macaulay2 editor](https://www.unimelb-macaulay2.cloud.edu.au/#editor).
First, run all of the code in the file [schur-Q-plethysm.m2](https://github.com/j-graf/Schur-Q-Plethysm/blob/main/schur-Q_plethysm.m2) (you may copy-paste the contents into the editor on the left, then click the orange "Run all editor code" button).
Then, you may run individual lines or blocks of code.

# Examples

## Schur's $Q$-function Basics

1. Compute $\ell(\lambda)$ and $|\lambda|$, where $\lambda=(6,5,2)$:
```
lam = {6,5,2}
#lam
sum lam
```

2. Compute $Q_{\lambda}$, where $\lambda=(6,5,2)$:
```
lam = {6,5,2}
Q lam
```

3. Compute the matrices $M(\lambda)$ and $M(\lambda,(0))$, and verify that $Q_{\lambda/0}=Q_\lambda$, where $\lambda=(6,5)$:
```
lam = {6,5}
M lam
M(lam,{0})
Q lam == skewQ(lam,{0})
```
