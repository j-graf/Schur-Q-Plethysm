# Schur-Q-Plethysm
Computes Schur's Q-functions as polynomials in the $q_n$'s with the Pfaffian formula, and computes the plethysm $Q_\lambda\circ Q_\mu$. Also gives conversions between Schur's $Q$-functions and power sums.

# How to Use

We recommend the [online Macaulay2 editor](https://www.unimelb-macaulay2.cloud.edu.au/#editor).
First, run all of the code in the file [schur-Q-plethysm.m2](https://github.com/j-graf/Schur-Q-Plethysm/blob/main/schur-Q_plethysm.m2) (you may copy-paste the contents into the editor on the left, then click the orange "Run all editor code" button).
Then, you may run individual lines or blocks of code.

# Examples

## Schur's $Q$-function Computations

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

3. Compute the matrices $M(\lambda)$, $M(\lambda0)$, and $M(\lambda00)$, and verify that $Q_\lambda=Q_{\lambda0}=Q_{\lambda00}$, where $\lambda=(6,5)$:
```
M {6,5}
M {6,5,0}
M {6,5,0,0}
Q {6,5} == Q {6,5,0} and Q {6,5} == Q {6,5,0,0}
```

4. Compute the matrices $M(\lambda)$ and $M(\lambda,0)$, and verify that $Q_\lambda=Q_{\lambda/0}$, where $\lambda=(6,5)$:
```
lam = {6,5}
M lam
skewM(lam,{0})
Q lam == skewQ(lam,{0})
```

5. Compute the matrix $M(p\lambda)$ and the function $Q_{p\lambda}$, where $p=-2$ and $\lambda=(5,2)$:
```
M {-2,5,2}
Q {-2,5,2}
```

6. Verify the example $Q_{(5,3,1)}=(q_2q_1-q_3)Q_{(5,3,2,1)/(5)}-q_1Q_{(5,3,2,1)/(3)}+Q_{(5,3,2,1)/(2)}$:
```
Q {5,3,1} == (q_2*q_1-q_3)*(skewQ({5,3,2,1},{5}))-q_1*(skewQ({5,3,2,1},{3}))+(skewQ({5,3,2,1},{2}))
```

## Plethysm Computations

7. Compute the plethysm $F\circ G$ for $F,G\in\mathbb{Q}[q_1,q_2,\ldots]$, where $F=Q_{(3,1)}$ and $G=q_4$:
```
pleth(Q {3,1}, Q {4})
```

## Power Sum Conversions

8. Write $F\in\mathbb{Q}[q_1,q_2,\ldots]$ as a polynomial in the $p_n$'s, where $F=Q_{(6,5)}+Q_{(3,1)}$:
```
TOp(Q {6,5} + Q {3,1})
```

9. Write $F\in\mathbb{Q}[p_1,p_2,\ldots]$ as a polynomial in the $q_n$'s, where $F=p_1p_3^2$:
```
TOq(p_1*p_3^2)
```
