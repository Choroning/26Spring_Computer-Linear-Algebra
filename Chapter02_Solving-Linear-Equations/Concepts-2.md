# Chapter 2 Lecture — Solving Linear Equations — Part 2/2

> **📑 This document is split into 2 parts.**
>
> [Part 1](Concepts.md) · **Part 2**

---

<br>

## Table of Contents

- [5. Permutations and Transpose (2.4)](#5-permutations-and-transpose-24)
  - [5.1 Permutation Matrices](#51-permutation-matrices)
  - [5.2 Properties of Permutation Matrices](#52-properties-of-permutation-matrices)
  - [5.3 The PA = LU Factorization](#53-the-pa--lu-factorization)
  - [5.4 Partial Pivoting](#54-partial-pivoting)
  - [5.5 PAQ: Row and Column Permutations](#55-paq-row-and-column-permutations)
  - [5.6 The Transpose of A](#56-the-transpose-of-a)
  - [5.7 Inner Products and the Transpose](#57-inner-products-and-the-transpose)
  - [5.8 Symmetric Matrices](#58-symmetric-matrices)
  - [5.9 Symmetric Products and LDL^T](#59-symmetric-products-and-ldlt)
- [6. Derivatives and Finite Difference Matrices (2.5)](#6-derivatives-and-finite-difference-matrices-25)
  - [6.1 Taylor Series and Approximations](#61-taylor-series-and-approximations)
  - [6.2 Derivatives from Differences](#62-derivatives-from-differences)
  - [6.3 Second Difference Matrices K, T, B](#63-second-difference-matrices-k-t-b)
  - [6.4 Properties of K](#64-properties-of-k)
  - [6.5 Free-Fixed Matrices T](#65-free-fixed-matrices-t)
  - [6.6 Free-Free Matrices B](#66-free-free-matrices-b)
- [Summary](#summary)

---

<br>

## 5. Permutations and Transpose (2.4)

### 5.1 Permutation Matrices

A **permutation matrix** $`P`$ has the same rows as $`I \in \mathbb{R}^{n \times n}`$.

There are $`n!`$ different orders.

**Example:** $`P \in \mathbb{R}^{3 \times 3}`$: 3 rows. At 1st row: 3 cases; at 2nd row: 2 cases; at 3rd row: 1 case $`\implies 3! = 6`$ orders.

$`P`$ times $`\mathbf{x}`$ puts the components $`x_1`$ to $`x_n`$ in that new order.

And $`P^T`$ equals $`P^{-1}`$.

**Example:**

```math
P\mathbf{x} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \end{pmatrix}
```

```math
P^T = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}
```

$`P^{-1} = ?`$ Using Gauss-Jordan method $`(A|I) \Rightarrow (I|A^{-1})`$:

```math
\begin{pmatrix} 0 & 0 & 1 & | & 1 & 0 & 0 \\ 1 & 0 & 0 & | & 0 & 1 & 0 \\ 0 & 1 & 0 & | & 0 & 0 & 1 \end{pmatrix} \xrightarrow{\text{row exchange}} \begin{pmatrix} 1 & 0 & 0 & | & 0 & 1 & 0 \\ 0 & 1 & 0 & | & 0 & 0 & 1 \\ 0 & 0 & 1 & | & 1 & 0 & 0 \end{pmatrix}
```

```math
P^{-1} = P^T
```

**Transpose properties:**

- Columns of $`A`$ are rows of $`A^T`$.
- The transposes of $`A\mathbf{x}`$ and $`AB`$ are $`\mathbf{x}^T A^T`$ and $`B^T A^T`$.

**Inner product property:**

```math
A\mathbf{x} \cdot \mathbf{y} = \mathbf{x} \cdot A^T\mathbf{y}
```

because $`(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T A^T \mathbf{y} = \mathbf{x}^T(A^T\mathbf{y})`$

**Symmetric matrix:** $`S^T = S`$. The products $`A^T A`$ and $`AA^T`$ are always symmetric.

---

### 5.2 Properties of Permutation Matrices

Permutation matrices have a 1 in every row and a 1 in every column.

When we multiply $`P`$ with a vector $`\mathbf{x}`$, it changes the order of its components:

```math
P\mathbf{x} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \end{pmatrix}
```

$`P`$ shifts $`x_1`$ to second position.

```math
PP\mathbf{x} = P^2\mathbf{x} = \begin{pmatrix} x_2 \\ x_3 \\ x_1 \end{pmatrix}
```

$`P^2`$ shifts $`x_1`$ to third position.

```math
PPP\mathbf{x} = P^3\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = I\mathbf{x}
```

$`P^3`$ recovers $`x_1`$ to its original position. $`P^3 = I`$.

**Consider 4x4 permutation matrices:**

**(a)** $`P`$ reverses the order of $`\mathbf{x}`$:

```math
\begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_4 \\ x_3 \\ x_2 \\ x_1 \end{pmatrix}
```

$`P(P\mathbf{x}) = \mathbf{x} \implies P^2 = I`$

**(b)** $`P`$ does not change $`x_4`$ position:

```math
\begin{pmatrix} 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \\ x_4 \end{pmatrix}
```

$`P(P(P\mathbf{x})) = \mathbf{x} \implies P^3 = I`$

**(c)** $`P`$ cyclically shifts the elements:

```math
\begin{pmatrix} 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_4 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix}
```

$`PPPP\mathbf{x} = P^4\mathbf{x} = I\mathbf{x} \implies P^4 = I`$

**(d)** Even-odd separation:

```math
\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_0 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_0 \\ x_2 \\ x_1 \\ x_3 \end{pmatrix} \quad \text{(even, odd separation)}
```

$`P^2 = I`$

This extends to $`\mathbf{x} \in \mathbb{R}^8`$ vectors using an 8x8 permutation matrix that separates even-indexed and odd-indexed entries.

**Proof that $`P^T P = I`$:**

The rows of any $`P`$ are the columns of $`P^{-1} = P^T`$.

```math
P^T P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}\begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}
```

```math
= \mathbf{h}_1\mathbf{h}_1^T + \mathbf{h}_2\mathbf{h}_2^T + \mathbf{h}_3\mathbf{h}_3^T
```

Since $`\mathbf{h}_i`$ are standard basis vectors (canonical unit vectors):

```math
= \begin{pmatrix} 1 & & \\ & 1 & \\ & & 1 \end{pmatrix} = I
```

**Properties of Permutation Matrices:**

1. A permutation matrix $`P`$ has exactly a 1 in each row and exactly a 1 in each column.

2. The columns of $`P`$ are **orthogonal**. (Dot products between columns are all zero.)

3. The product $`P_1 P_2`$ of permutations is a permutation. So is the inverse of $`P`$.

4. If $`A`$ is invertible, then there exists a permutation $`P`$ to order its rows in advance, so that elimination on $`PA`$ meets no zeros in the pivot positions:

```math
PA = LU
```

---

### 5.3 The PA = LU Factorization

**Row exchanges from $`P`$:**

Consider a matrix $`A`$:

```math
A = \begin{pmatrix} 1 & 2 & a \\ 2 & 4 & b \\ 3 & 7 & c \end{pmatrix}
```

Take 1 as pivot in row 1:

```math
R_2 - 2R_1, \quad R_3 - 3R_1
```

```math
EA = \begin{pmatrix} 1 & 2 & a \\ 0 & 0 & b - 2a \\ 0 & 1 & c - 3a \end{pmatrix}
```

Due to zero pivot, swap row 2 and row 3:

```math
PEA = \begin{pmatrix} 1 & 2 & a \\ 0 & 1 & c - 3a \\ 0 & 0 & b - 2a \end{pmatrix} = U
```

$`A`$ is invertible iff $`b - 2a \neq 0`$. When $`b = 2a`$, then $`\text{rank}(A) = 2 < 3`$, $`A`$ is not invertible.

**We can exchange rows 2 and 3 first:**

```math
PA = \begin{pmatrix} 1 & 2 & a \\ 3 & 7 & c \\ 2 & 4 & b \end{pmatrix}
```

```math
\xrightarrow{R_2 - 3R_1, R_3 - 2R_1}
```

```math
EPA = \begin{pmatrix} 1 & 2 & a \\ 0 & 1 & c - 3a \\ 0 & 0 & b - 2a \end{pmatrix} = U
```

```math
\therefore PA = E^{-1}U = LU
```

**Daniel Drucker's method for tracking $`P`$:**

Augment the matrix with a column tracking row indices:

```math
\begin{pmatrix} 1 & 2 & a & | & 1 \\ 2 & 4 & b & | & 2 \\ 3 & 7 & c & | & 3 \end{pmatrix} \xrightarrow{\text{elimination}} \begin{pmatrix} 1 & 2 & a & | & 1 \\ 0 & 0 & b-2a & | & 2 \\ 0 & 1 & c-3a & | & 3 \end{pmatrix} \xrightarrow{\text{swap}} \begin{pmatrix} 1 & 2 & a & | & 1 \\ 0 & 1 & c-3a & | & 3 \\ 0 & 0 & b-2a & | & 2 \end{pmatrix}
```

The final column gives $`P_{132}`$:

```math
P_{132} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}
```

---

### 5.4 Partial Pivoting

**"Partial Pivoting"** to reduce roundoff errors.

The computation is more stable if we exchange rows to produce the **largest possible number** in the pivot.

**Example:**

```math
\begin{pmatrix} 1 & 2 & a & | & 1 \\ 2 & 4 & b & | & 2 \\ 3 & 7 & c & | & 3 \end{pmatrix}
```

```math
\xrightarrow{R_3 \leftrightarrow R_1} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 2 & 4 & b & | & 2 \\ 1 & 2 & a & | & 1 \end{pmatrix}
```

```math
\xrightarrow{R_2 - \frac{2}{3}R_1, R_3 - \frac{1}{3}R_1} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 0 & -\frac{2}{3} & b - \frac{2}{3}c & | & 2 \\ 0 & -\frac{1}{3} & a - \frac{1}{3}c & | & 1 \end{pmatrix}
```

```math
\xrightarrow{R_3' - R_2'(\frac{1}{2})} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 0 & -\frac{2}{3} & b - \frac{2}{3}c & | & 2 \\ 0 & 0 & a - \frac{1}{2}b & | & 1 \end{pmatrix} = U
```

All entries of $`L`$ are $`\leq 1`$ when we make each pivot larger than all numbers below it:

```math
L = \begin{pmatrix} 1 & 0 & 0 \\ 2/3 & 1 & 0 \\ 1/3 & 1/2 & 1 \end{pmatrix}
```

```math
P_{321} = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}
```

---

### 5.5 PAQ: Row and Column Permutations

$`PAQ`$ has row permutation $`P`$ and column permutation $`Q`$.

Start with $`A \in \mathbb{R}^{3 \times 3}`$. Reorder its rows by $`P \in \mathbb{R}^{3 \times 3}`$.

```math
P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}, \quad PA = \begin{pmatrix} a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \\ a_{11} & a_{12} & a_{13} \end{pmatrix}
```

Multiply $`Q = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}`$ from the right:

```math
(PA)Q = \begin{pmatrix} a_{23} & a_{22} & a_{21} \\ a_{33} & a_{32} & a_{31} \\ a_{13} & a_{12} & a_{11} \end{pmatrix}
```

Column permutation $`Q`$ reorders those columns.

Is the column space of $`A`$ equal to the column space of $`PA`$?

```math
C(A) \stackrel{?}{=} C(PA)
```

**Yes**, $`P`$ does not change the linear relationships. Thus $`C(A) = C(PA)`$.

**Q:** A matrix $`A`$ has 9 numbers. How many different ways can you arrange the 9 numbers in $`A`$? **A:** $`9!`$

$`P, Q`$ can arrange the numbers $`6 \times 6 = 36`$ ways for $`PAQ`$. $`PAQ`$ is very special, satisfying $`C(A) = C(PAQ)`$.

---

### 5.6 The Transpose of A

The **transpose** of $`A`$, denoted by $`A^T`$. The columns of $`A^T`$ are the rows of $`A`$.

```math
A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & \vdots & \ddots & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}_{n \times m}
```

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}_{m \times n}
```

**Example:**

```math
A = \begin{pmatrix} 1 & 2 & 3 \\ 0 & 0 & 4 \end{pmatrix}, \quad A^T = \begin{pmatrix} 1 & 0 \\ 2 & 0 \\ 3 & 4 \end{pmatrix}
```

The matrix "flips over" its main diagonal.

```math
(A^T)_{ij} = A_{ji}
```

**Rules for transpose:**

- **Sum:** $`(A + B)^T = A^T + B^T`$
- **Product:** $`(AB)^T = B^T A^T`$ (reverse order)
- **Inverse:** $`(A^{-1})^T = (A^T)^{-1}`$

**Proof of $`(A\mathbf{x})^T = \mathbf{x}^T A^T`$:**

```math
A\mathbf{x} = \begin{pmatrix} \sum_{j=1}^n a_{1j}x_j \\ \sum_{j=1}^n a_{2j}x_j \\ \vdots \\ \sum_{j=1}^n a_{mj}x_j \end{pmatrix}
```

```math
(A\mathbf{x})^T = \left(\sum_{j=1}^n a_{1j}x_j \quad \sum_{j=1}^n a_{2j}x_j \quad \cdots \quad \sum_{j=1}^n a_{mj}x_j\right)
```

```math
= (x_1, x_2, \ldots, x_n)\begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}_{n \times m}
```

```math
= \mathbf{x}^T A^T
```

**Proof of $`(AB)^T = B^T A^T`$:**

Interpret $`B = \begin{pmatrix} | & | & & | \\ \mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_p \\ | & | & & | \end{pmatrix}_{n \times p}`$

```math
(AB)^T = (A\mathbf{x}_1 \quad A\mathbf{x}_2 \quad \cdots \quad A\mathbf{x}_p)^T = \begin{pmatrix} \mathbf{x}_1^T A^T \\ \mathbf{x}_2^T A^T \\ \vdots \\ \mathbf{x}_p^T A^T \end{pmatrix} = \begin{pmatrix} \mathbf{x}_1^T \\ \mathbf{x}_2^T \\ \vdots \\ \mathbf{x}_p^T \end{pmatrix} A^T = B^T A^T
```

**Example:**

```math
AB = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 5 & 0 \\ 4 & 1 \end{pmatrix} = \begin{pmatrix} 5 & 0 \\ 9 & 1 \end{pmatrix}
```

```math
B^T A^T = \begin{pmatrix} 5 & 4 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 5 & 9 \\ 0 & 1 \end{pmatrix} = (AB)^T
```

Also: $`(ABC)^T = (AB \cdot C)^T = C^T(AB)^T = C^T B^T A^T`$. The reverse order rule still holds true.

**Proof of $`(A^{-1})^T = (A^T)^{-1}`$:**

Consider $`A^{-1}A = I`$:

```math
\Rightarrow (A^{-1}A)^T = I^T = I
```
```math
\Leftrightarrow A^T(A^{-1})^T = I
```

$`\therefore (A^{-1})^T`$ is the inverse of $`A^T`$:

```math
(A^{-1})^T = (A^T)^{-1}
```

This implies that $`A^T`$ is invertible exactly when $`A`$ is invertible.

**Example:**

```math
A = \begin{pmatrix} 1 & 0 \\ 6 & 1 \end{pmatrix}, \quad A^{-1} = \begin{pmatrix} 1 & 0 \\ -6 & 1 \end{pmatrix} \implies (A^{-1})^T = \begin{pmatrix} 1 & -6 \\ 0 & 1 \end{pmatrix}
```

```math
A^T = \begin{pmatrix} 1 & 6 \\ 0 & 1 \end{pmatrix} \implies (A^T)^{-1} = \begin{pmatrix} 1 & -6 \\ 0 & 1 \end{pmatrix}
```

They agree: $`(A^{-1})^T = (A^T)^{-1}`$.

---

### 5.7 Inner Products and the Transpose

**The meaning of inner products.**

The dot product of $`\mathbf{x}`$ and $`\mathbf{y}`$ (inner product) is the sum of numbers $`x_i y_i`$:

```math
\mathbf{x} \cdot \mathbf{y} = \sum x_i y_i
```

Let $`\mathbf{x}, \mathbf{y} \in \mathbb{R}^n`$:

```math
\mathbf{x} \cdot \mathbf{y} = \sum_{i=1}^n x_i y_i = \mathbf{x}^T \mathbf{y} \quad (1 \times n)(n \times 1) = (1 \times 1)
```

```math
\mathbf{x}\mathbf{y}^T = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix}(y_1, y_2, \ldots, y_n) = \mathbf{x} \otimes \mathbf{y} \quad (n \times 1)(1 \times n) = (n \times n)
```

This is the **outer product** (rank one product, rank one matrix).

**Examples of dot products:**

```math
\text{work} = \text{force} \cdot \text{distance} = \mathbf{f}^T \mathbf{d} \quad [J] = [N] \cdot [m]
```

```math
\text{income} = \text{quantities} \cdot \text{prices} = \mathbf{q}^T \mathbf{p}
```

**A better definition of $`A^T`$:**

We define $`A^T`$ by flipping the matrix $`A`$ across its main diagonal, but a better definition of $`A^T`$ is that $`A^T`$ is the matrix that makes the following inner products equal:

```math
(A\mathbf{x}) \cdot \mathbf{y} = \mathbf{x} \cdot (A^T\mathbf{y})
```

i.e.,

```math
(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T A^T \mathbf{y}
```

**Example 1:**

```math
A = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}, \quad \mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix} y_1 \\ y_2 \end{pmatrix}
```

```math
A\mathbf{x} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \end{pmatrix}
```

```math
(A\mathbf{x})^T\mathbf{y} = (x_2 - x_1)y_1 + (x_3 - x_2)y_2
```

```math
A^T\mathbf{y} = \begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} -y_1 \\ y_1 - y_2 \\ y_2 \end{pmatrix}
```

```math
\mathbf{x}^T(A^T\mathbf{y}) = (x_1, x_2, x_3)\begin{pmatrix} -y_1 \\ y_1 - y_2 \\ y_2 \end{pmatrix} = -x_1y_1 + x_2(y_1 - y_2) + x_3y_2
```

Both are equal. $`\checkmark`$

**Example 2: Inner product of functions.**

```math
\mathbf{x}^T\mathbf{y} = x_1 y_1 + x_2 y_2 + \cdots + x_n y_n
```

In the continuous world:

```math
(x, y) := \int_{-\infty}^{\infty} x(t) \, y(t) \, dt
```

Similarly, $`(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T(A^T\mathbf{y})`$:

```math
(Ax, y) = (x, A^T y)
```

Let $`A = \frac{d}{dt}`$, $`A^T = -\frac{d}{dt} = -A`$.

```math
\left(\frac{dx}{dt}, y\right) = \left(x, -\frac{dy}{dt}\right)
```

```math
\int_{-\infty}^{\infty} \frac{dx}{dt} y \, dt = \int_{-\infty}^{\infty} x \left(-\frac{dy}{dt}\right) dt
```

**Proof via integration by parts (IBP):**

```math
(f(t)g(t))' = f'(t)g(t) + f(t)g'(t)
```

```math
\int f'g \, dt = \int (fg)' \, dt - \int fg' \, dt = (fg)\Big|_{-\infty}^{\infty} - \int fg' \, dt
```

Assuming $`f(\infty) = f(-\infty) = 0`$:

```math
\boxed{\int_{-\infty}^{\infty} f'g \, dt = -\int_{-\infty}^{\infty} fg' \, dt}
```

The derivative is **anti-symmetric**: $`A^T = -A`$.

Symmetric matrices have $`A^T = A`$.

---

### 5.8 Symmetric Matrices

A **symmetric matrix** has $`S^T = S`$:

```math
\implies (S^T)_{ij} = S_{ij} = S_{ji}
```

**Examples:**

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix} = S^T, \quad D = \begin{pmatrix} 1 & 0 \\ 0 & 10 \end{pmatrix} = D^T
```

**The inverse of a symmetric matrix is a symmetric matrix:**

```math
(S^{-1})^T = (S^T)^{-1} = S^{-1}
```

$`\Rightarrow`$ When $`S`$ is invertible, $`S^{-1}`$ is symmetric.

**Example:**

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}, \quad S^{-1} = \begin{pmatrix} 5 & -2 \\ -2 & 1 \end{pmatrix} = (S^{-1})^T
```

---

### 5.9 Symmetric Products and LDL^T

**Symmetric Products $`A^T A`$, $`AA^T`$, $`LDL^T`$:**

For $`A \in \mathbb{R}^{m \times n}`$:

- $`A^T A \in \mathbb{R}^{n \times n}`$: $`(A^T A)^T = A^T A \implies`$ **symmetric**
- $`AA^T \in \mathbb{R}^{m \times m}`$: $`(AA^T)^T = AA^T \implies`$ **symmetric**

**Example 3:**

```math
A = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}
```

```math
AA^T = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}
```

```math
A^T A = \begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}
```

**Symmetric matrices in elimination:** $`S^T = S`$ makes elimination **twice as fast**.

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 7 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix} = U
```

```math
L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}
```

```math
S = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}
```

We do not see the symmetry of $`S`$ in $`LU`$ decomposition. For symmetric matrix $`S`$, we can further decompose $`U`$ into $`D`$ and $`L^T`$:

```math
U = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix} = \begin{pmatrix} d_1 & \\ & d_2 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & \\ & 3 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}
```

```math
S = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 3 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = LDL^T
```

**Example 4: Saddle-point matrix.**

For a rectangular $`A \in \mathbb{R}^{m \times n}`$, the **saddle-point matrix** $`S`$ is symmetric and important:

```math
S = \begin{pmatrix} I_{m \times m} & A_{m \times n} \\ A^T_{n \times m} & 0_{n \times n} \end{pmatrix} = S^T \quad (m+n) \times (m+n)
```

Block elimination: $`R_2 - A^T R_1`$:

```math
ES = \begin{pmatrix} I & A \\ 0 & -A^T A \end{pmatrix} = U
```

**Block factorization:**

```math
S = \begin{pmatrix} I & 0 \\ A^T & I \end{pmatrix}\begin{pmatrix} I & 0 \\ 0 & -A^T A \end{pmatrix}\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} = LDL^T
```

$`S`$ is invertible $`\iff`$ $`A^T A`$ is invertible $`\iff`$ $`\text{rank}(A^T A) = n`$ $`\iff`$ $`A\mathbf{x} \neq \mathbf{0}`$ whenever $`\mathbf{x} \neq \mathbf{0}`$ $`\iff`$ columns of $`A`$ are linearly independent.

---

<br>

## 6. Derivatives and Finite Difference Matrices (2.5)

### 6.1 Taylor Series and Approximations

The matrices imitate **derivatives**. Derivatives tell us what is happening at one point $`x`$ of space or at one moment $`t`$ in time.

**Example:** $`y = x^2 + 2`$

```math
\frac{dy}{dx} = 2x \xrightarrow{x=1} 2 > 0
```

```math
\frac{d^2y}{dx^2} = 2 > 0
```

The graph of $`y`$ bends upward (slope is increasing).

**Consider** $`\Delta y = y(x+h) - y(x)`$:

1. The difference is approximately: $`\Delta y \approx h \cdot \frac{dy}{dx}\Big|_{x=x_1}`$

2. For better approximation: $`\Delta y \approx h \cdot \frac{dy}{dx}\Big|_{x=x_1} + \frac{h^2}{2}\frac{d^2y}{dx^2}\Big|_{x=x_1}`$ (tangent line + parabola)

3. The exact $`\Delta y`$ is the integral: $`\Delta y = y(x+h) - y(x) = \int_x^{x+h} \frac{dy}{dx} \, dx`$

The accuracy of $`\Delta y`$ increases by adding the derivative terms.

**Taylor Series:**

```math
y(x+h) = y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + \cdots + \frac{h^n}{n!}\frac{d^{(n)}y}{dx^n} + \cdots
```

**Example:** $`e^x`$ (entire analytic function):

```math
e^{x+h} = e^x + h \cdot e^x + \frac{h^2}{2} \cdot e^x + \cdots + \frac{h^n}{n!} e^x + \cdots
```

```math
= e^x\left(1 + h + \frac{h^2}{2} + \cdots + \frac{h^n}{n!} + \cdots\right)
```

```math
\therefore e^h = 1 + h + \frac{h^2}{2} + \cdots + \frac{h^n}{n!} + \cdots
```

---

### 6.2 Derivatives from Differences

**Turning the Formulas Around: Derivatives from Differences**

Start with a tangent parabola:

```math
y(x+h) \approx y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2}
```

**Q:** If we know $`y(x)`$ and $`y(x+h)`$, then how do we estimate $`\frac{dy}{dx}`$?

**Q:** If we know $`y(x-h)`$, then can we have an estimation of $`\frac{dy}{dx}`$, $`\frac{d^2y}{dx^2}`$?

Using **Finite difference method**, we can approximate the derivatives.

**Forward Difference:**

```math
y(x+h) \approx y(x) + h\frac{dy}{dx}
```

```math
\implies \frac{dy}{dx} \approx \frac{y(x+h) - y(x)}{h} \quad \text{(Forward Difference)}
```

This is the **first order** approximation:

```math
y(x+h) = y(x) + h\frac{dy}{dx} + O(h^2)
```

```math
\implies \frac{dy}{dx} = \frac{y(x+h) - y(x)}{h} + O(h) \quad \text{(truncation error)}
```

**Backward Difference:**

```math
y(x-h) \approx y(x) - h\frac{dy}{dx}
```

```math
\implies \frac{dy}{dx} \approx \frac{y(x) - y(x-h)}{h}
```

**Q:** Can we increase the accuracy of the approximation to $`\frac{dy}{dx}`$?

**Centered Difference (2nd order accuracy):**

```math
y(x+h) = y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + O(h^3)
```

```math
y(x-h) = y(x) - h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + O(h^3)
```

Subtracting:

```math
y(x+h) - y(x-h) = 2h\frac{dy}{dx} + O(h^3)
```

```math
\implies \frac{dy}{dx} = \frac{y(x+h) - y(x-h)}{2h} + O(h^2) \approx \frac{y(x+h) - y(x-h)}{2h}
```

This centered difference formula has the **2nd order accuracy**.

**Second Difference (approximation to $`\frac{d^2y}{dx^2}`$):**

Adding the two Taylor expansions:

```math
y(x+h) + y(x-h) = 2y(x) + h^2\frac{d^2y}{dx^2} + O(h^4)
```

```math
\implies \frac{d^2y}{dx^2} = \frac{y(x+h) - 2y(x) + y(x-h)}{h^2} + O(h^2)
```

```math
\approx \frac{y(x+h) - 2y(x) + y(x-h)}{h^2} \quad \text{(Second difference)}
```

```math
\frac{d^2y}{dx^2} \approx \frac{1}{h^2}(1 \quad -2 \quad 1)\begin{pmatrix} y(x-h) \\ y(x) \\ y(x+h) \end{pmatrix}
```

---

### 6.3 Second Difference Matrices K, T, B

**Consider a 1D domain:**

```math
x_0, x_1, x_2, \ldots, x_{N-1}, x_N, x_{N+1}
```

We decompose the whole domain into $`N+1`$ non-overlapping elements with uniform spacing $`h`$.

Here $`x_0 = 0`$, $`x_{N+1} = 1`$ are boundary points.

```math
h = \frac{1}{N+1} \implies x_i = ih = \frac{i}{N+1}
```

**Discretize** the equation $`-\frac{d^2u}{dx^2} = f(x)`$ with $`u(0) = u(1) = 0`$ (boundary conditions).

From BC: $`u_0 = u_{N+1} = 0`$.

**Q:** What are $`u_1, u_2, \ldots, u_N`$?

```math
\left.\frac{d^2u}{dx^2}\right|_{x=x_1} \approx \frac{u_0 - 2u_1 + u_2}{h^2} = -f(x_1)
```

```math
\left.\frac{d^2u}{dx^2}\right|_{x=x_2} \approx \frac{u_1 - 2u_2 + u_3}{h^2} = -f(x_2)
```

```math
\vdots
```

```math
\left.\frac{d^2u}{dx^2}\right|_{x=x_N} \approx \frac{u_{N-1} - 2u_N + u_{N+1}}{h^2} = -f(x_N)
```

In matrix form:

```math
\frac{1}{h^2}\begin{pmatrix} 2 & -1 & 0 & \cdots & 0 \\ -1 & 2 & -1 & \cdots & 0 \\ 0 & -1 & 2 & \cdots & 0 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & -1 & 2 \end{pmatrix}\begin{pmatrix} u_1 \\ u_2 \\ u_3 \\ \vdots \\ u_{N-1} \\ u_N \end{pmatrix} = \begin{pmatrix} f_1 \\ f_2 \\ f_3 \\ \vdots \\ f_{N-1} \\ f_N \end{pmatrix}
```

where $`f_i := f(x_i)`$.

```math
\boxed{\frac{1}{h^2}K\mathbf{u} = \mathbf{f}}
```

The matrix $`K`$ gives a natural approximation to $`-\frac{d^2u}{dx^2}`$ with **fixed-fixed BC**: $`u_0 = 0`$ and $`u_{N+1} = 0`$.

---

### 6.4 Properties of K

Let $`N = 4`$:

```math
K_4 = \begin{pmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}
```

**Properties:**

**1. $`K`$ is symmetric.** $`K^T = K`$.

**2. $`K`$ is banded.** All the nonzeros in $`K`$ lie in a "band" around the main diagonal. A matrix with a narrow band is **sparse** (mostly zeros). It is a **tridiagonal matrix**.

Example: $`N = 100`$:
- Number of 2's: 100
- Number of $`-1`$'s: $`99 + 99 = 198`$
- Number of nonzeros $`\approx 300`$. Out of $`10000`$ entries: $`\frac{300}{10000} = 3\%`$.

**3. $`K`$ has constant diagonals,** which is related to Fourier transforms, filters, convolution matrices, Toeplitz matrix. $`K`$ is **shift-invariant** because $`(-1, 2, -1)`$ pattern appears in each row.

**4. $`K`$ is invertible.** We can check it by elimination. If the resulting upper triangular matrix $`U`$ has no zero pivot, then $`K`$ is invertible.

**5. The symmetric matrices $`K_n`$ are positive definite:**

```math
\mathbf{x}^T K \mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}
```

**Example:** $`K = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}`$

```math
K\mathbf{x} = \begin{pmatrix} 2x - y \\ -x + 2y \end{pmatrix}
```

```math
\mathbf{x}^T(K\mathbf{x}) = 2x^2 - 2xy + 2y^2 = x^2 + (x-y)^2 + y^2 > 0
```

When $`\mathbf{x}^T K \mathbf{x} \geq 0 \; \forall \; \mathbf{x} \neq \mathbf{0}`$, we say $`K`$ is **positive semi-definite**.

**Pivots:**
- An invertible matrix has $`n`$ nonzero pivots.
- A positive definite symmetric matrix has $`n`$ **positive** pivots.
- A positive semi-definite symmetric matrix has $`n`$ **nonnegative** pivots.

---

### 6.5 Free-Fixed Matrices T

```math
-\frac{d^2u}{dx^2} = f(x), \quad \text{with } \frac{du}{dx} = 0 \text{ at } x = 0, \quad u(1) = 0
```

The Neumann BC at $`x = 0`$: $`\frac{u_1 - u_0}{h} = 0 \implies \frac{1}{h^2}(u_0 - u_1) = 0`$. Here $`u_0`$ is **unknown**.

This gives the matrix $`T`$:

```math
\frac{1}{h^2}\begin{pmatrix} 1 & -1 & 0 & 0 & \cdots & 0 \\ -1 & 2 & -1 & 0 & \cdots & 0 \\ 0 & -1 & 2 & -1 & \cdots & 0 \\ \vdots & & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & -1 & 2 \end{pmatrix}\begin{pmatrix} u_0 \\ u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N \end{pmatrix} = \begin{pmatrix} 0 \\ f_1 \\ f_2 \\ \vdots \\ f_{N-1} \\ f_N \end{pmatrix}
```

For $`N = 4`$:

```math
T_4 = \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}
```

**Elimination of $`T_4`$:**

```math
\xrightarrow{R_2 + R_1} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix} \xrightarrow{R_3' + R_2'} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix} \xrightarrow{R_4'' + R_3''} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 1 \end{pmatrix} = U = L^T
```

```math
\therefore T_4 = LL^T
```

```math
L = \begin{pmatrix} 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}
```

**What is $`L^{-1}`$?** Use Gauss-Jordan method: $`(L|I) \Rightarrow (I|L^{-1})`$

```math
\begin{pmatrix} 1 & 0 & 0 & 0 & | & 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 & | & 0 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 & | & 0 & 0 & 1 & 0 \\ 0 & 0 & -1 & 1 & | & 0 & 0 & 0 & 1 \end{pmatrix} \implies \begin{pmatrix} 1 & 0 & 0 & 0 & | & 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & | & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 & | & 1 & 1 & 1 & 0 \\ 0 & 0 & 0 & 1 & | & 1 & 1 & 1 & 1 \end{pmatrix}
```

```math
L^{-1} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 1 & 1 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix}
```

```math
T_4^{-1} = (LL^T)^{-1} = L^{-T}L^{-1} = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 1 & 1 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix} = \begin{pmatrix} 4 & 3 & 2 & 1 \\ 3 & 3 & 2 & 1 \\ 2 & 2 & 2 & 1 \\ 1 & 1 & 1 & 1 \end{pmatrix}
```

---

### 6.6 Free-Free Matrices B

**Free-Free matrices $`B`$ are singular.**

- Not invertible
- There exists nonzero $`\mathbf{x}`$ s.t. $`B\mathbf{x} = \mathbf{0}`$

For the equation $`-\frac{d^2u}{dx^2} = f(x)`$ with:

```math
\frac{du}{dx} = 0 \text{ at } x = 0, \quad \frac{du}{dx} = 0 \text{ at } x = 1
```

```math
\frac{u_1 - u_0}{h} = 0 \implies \frac{1}{h^2}(u_0 - u_1) = 0 \quad (u_0 \text{ is unknown})
```

```math
\frac{u_{N+1} - u_N}{h} = 0 \implies \frac{1}{h^2}(-u_N + u_{N+1}) = 0 \quad (u_{N+1} \text{ is unknown})
```

```math
\frac{1}{h^2}\begin{pmatrix} 1 & -1 & 0 & 0 & \cdots & 0 & 0 \\ -1 & 2 & -1 & 0 & \cdots & 0 & 0 \\ 0 & -1 & 2 & -1 & \cdots & 0 & 0 \\ \vdots & & & \ddots & & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 & 0 \\ 0 & 0 & \cdots & 0 & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} u_0 \\ u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N \\ u_{N+1} \end{pmatrix} = \begin{pmatrix} 0 \\ f_1 \\ f_2 \\ \vdots \\ f_{N-1} \\ f_N \\ 0 \end{pmatrix}
```

**Consider $`B_3`$:**

```math
B_3 = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}
```

```math
B_3 \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`\mathbf{x} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$ is in the **null space** of $`B_3`$.

$`\mathbf{x} = c\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$ vectors are in the null space of $`B_3`$.

$`\Rightarrow`$ $`B_3`$ **cannot be invertible.**

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| $`A\mathbf{x} = \mathbf{b}`$ | System of $`n`$ linear equations with $`n`$ unknowns |
| Elimination | Subtract $`l_{ij}`$ times row $`j`$ from row $`i`$ to produce zeros below pivots |
| Back substitution | Solve $`U\mathbf{x} = \mathbf{c}`$ from bottom row upward |
| Three cases | Unique solution ($`\text{rank} = n`$), no solution, infinitely many solutions |
| Augmented matrix | $`(A \mid \mathbf{b})`$ tracks both sides during elimination |
| Homogeneous system | $`A\mathbf{x} = \mathbf{0}`$; nontrivial solutions exist when $`\text{rank}(A) < n`$ |
| Pivots | Diagonal entries of $`U`$; must be nonzero for invertibility |
| Row exchange | Swap rows when zero appears in pivot position |
| Elimination matrix $`E_{ij}`$ | Identity matrix with $`-l_{ij}`$ in position $`(i,j)`$ |
| $`EA = U`$ | Product of all elimination matrices transforms $`A`$ to upper triangular $`U`$ |
| $`A = LU`$ | $`L = E^{-1}`$ is lower triangular with multipliers $`l_{ij}`$ in correct positions |
| Inverse $`A^{-1}`$ | Exists iff $`A`$ has $`n`$ independent columns ($`\text{rank}(A) = n`$) |
| Uniqueness of inverse | Left inverse equals right inverse; $`BA = I`$ and $`AC = I \implies B = C`$ |
| $`(AB)^{-1} = B^{-1}A^{-1}`$ | Inverses come in reverse order |
| $`2 \times 2`$ inverse | $`A^{-1} = \frac{1}{ad-bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}`$; requires $`\det(A) \neq 0`$ |
| Gauss-Jordan | $`(A \mid I) \Rightarrow (I \mid A^{-1})`$ to find inverse explicitly |
| Cost of elimination | $`A \to U`$: $`\frac{1}{3}n^3`$; $`\mathbf{b} \to \mathbf{c} \to \mathbf{x}`$: $`n^2`$ |
| Second proof of $`A = LU`$ | Column times row: $`A = \sum \mathbf{l}_k \mathbf{u}_k`$ |
| $`A = LU`$ without row exchange | All upper-left $`k \times k`$ submatrices must be invertible |
| Permutation matrix $`P`$ | Rows of $`I`$ in different order; $`P^{-1} = P^T`$; $`n!`$ possible |
| $`PA = LU`$ | General factorization with row exchanges recorded in $`P`$ |
| Partial pivoting | Exchange rows to make largest element the pivot; reduces roundoff errors |
| $`PAQ`$ | Row permutation $`P`$, column permutation $`Q`$; $`C(A) = C(PA)`$ |
| Transpose $`A^T`$ | $`(A^T)_{ij} = A_{ji}`$; columns of $`A^T`$ are rows of $`A`$ |
| Transpose rules | $`(A+B)^T = A^T + B^T`$; $`(AB)^T = B^T A^T`$; $`(A^{-1})^T = (A^T)^{-1}`$ |
| Inner product | $`\mathbf{x} \cdot \mathbf{y} = \mathbf{x}^T\mathbf{y}`$; $`(A\mathbf{x}) \cdot \mathbf{y} = \mathbf{x} \cdot (A^T\mathbf{y})`$ |
| Outer product | $`\mathbf{x}\mathbf{y}^T`$ is a rank-one $`n \times n`$ matrix |
| Symmetric matrix $`S`$ | $`S^T = S`$; $`S^{-1}`$ is also symmetric |
| $`A^T A`$ and $`AA^T`$ | Always symmetric; $`A^T A`$ invertible iff columns of $`A`$ are independent |
| $`S = LDL^T`$ | Symmetric factorization; reveals symmetry that $`LU`$ does not |
| IBP and transpose | $`A = d/dt \implies A^T = -d/dt`$ (anti-symmetric); $`(Ax, y) = (x, A^T y)`$ |
| Forward difference | $`\frac{dy}{dx} \approx \frac{y(x+h)-y(x)}{h}`$; 1st order accuracy $`O(h)`$ |
| Centered difference | $`\frac{dy}{dx} \approx \frac{y(x+h)-y(x-h)}{2h}`$; 2nd order accuracy $`O(h^2)`$ |
| Second difference | $`\frac{d^2y}{dx^2} \approx \frac{y(x+h)-2y(x)+y(x-h)}{h^2}`$; 2nd order accuracy |
| Matrix $`K`$ (fixed-fixed) | Tridiagonal $`(-1, 2, -1)`$; symmetric, banded, invertible, positive definite |
| Matrix $`T`$ (free-fixed) | First row is $`(1, -1, 0, \ldots)`$; $`T = LL^T`$; invertible |
| Matrix $`B`$ (free-free) | Singular; $`B\mathbf{1} = \mathbf{0}`$; constant vector is in the null space |
| Positive definite | $`\mathbf{x}^T K\mathbf{x} > 0`$ for all $`\mathbf{x} \neq \mathbf{0}`$; all pivots are positive |
| Positive semi-definite | $`\mathbf{x}^T K\mathbf{x} \geq 0`$ for all $`\mathbf{x} \neq \mathbf{0}`$; all pivots are nonnegative |

---
