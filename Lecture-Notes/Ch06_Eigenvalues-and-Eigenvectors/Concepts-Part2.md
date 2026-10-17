# Chapter 6 Lecture — Eigenvalues and Eigenvectors — Part 2/3

> **📑 This document is split into 3 parts.**
>
> [Part 1](Concepts-Part1.md) · **Part 2** · [Part 3](Concepts-Part3.md)

---

<br>

## Table of Contents

- [3. Symmetric Positive Definite Matrices (6.3)](#3-symmetric-positive-definite-matrices-63)
  - [3.1 Symmetric Matrices: Key Properties](#31-symmetric-matrices-key-properties)
  - [3.2 Spectral Theorem](#32-spectral-theorem)
  - [3.3 Proof: Symmetric Matrices Have Orthonormal Eigenbasis](#33-proof-symmetric-matrices-have-orthonormal-eigenbasis)
  - [3.4 Positive Definite Matrices: Definitions](#34-positive-definite-matrices-definitions)
  - [3.5 Properties of Positive Definite Matrices](#35-properties-of-positive-definite-matrices)
  - [3.6 How to Check if a Matrix is Positive Definite](#36-how-to-check-if-a-matrix-is-positive-definite)
  - [3.7 Worked Examples: Positive Definite and Semidefinite](#37-worked-examples-positive-definite-and-semidefinite)
  - [3.8 The Ellipse and Quadratic Forms](#38-the-ellipse-and-quadratic-forms)
  - [3.9 Positive Definite Matrices and Minimum Problems](#39-positive-definite-matrices-and-minimum-problems)
  - [3.10 Positive Semidefinite Matrices](#310-positive-semidefinite-matrices)
  - [3.11 Congruent Matrices](#311-congruent-matrices)
  - [3.12 Optimization and Machine Learning](#312-optimization-and-machine-learning)

---

<br>

## 3. Symmetric Positive Definite Matrices (6.3)

### 3.1 Symmetric Matrices: Key Properties

**(1)** A symmetric matrix $`S`$ has $`n`$ **real** eigenvalues $`\lambda_i`$, and $`n`$ **orthonormal** eigenvectors $`\mathbf{q}_i`$.

**(2)** $`S`$ is diagonalized by an orthogonal eigenvector matrix $`Q`$:

```math
S = Q\Lambda Q^{-1} = Q\Lambda Q^T
```

**(3a)** **Positive definite:** all $`\lambda_i > 0`$ and all pivots $`> 0`$ and all upper-left determinants $`> 0`$.

**(3b)** The energy test is $`\mathbf{x}^T S\mathbf{x} > 0`$ for all $`\mathbf{x} \neq \mathbf{0}`$. Then $`S = A^T A`$ with independent columns in $`A`$.

**(4)** **Positive semidefinite** allows $`\lambda = 0`$, pivot = 0, determinant = 0, energy $`\mathbf{x}^T S\mathbf{x} = 0`$, any $`A`$.

**Symmetric matrices ($`S = S^T`$) are special:**

1. All $`n`$ eigenvalues $`\lambda`$ are **real** numbers.
2. The $`n`$ eigenvectors $`\mathbf{q}`$ can be chosen **orthogonal**.

e.g., $`I = I^T = \begin{pmatrix} 1 & & \\ & 1 & \\ & & \ddots \end{pmatrix}`$, $`\lambda = 1, 1, \ldots, 1`$. $`I\mathbf{x} = 1 \cdot \mathbf{x}`$ — every nonzero vector $`\mathbf{x}`$ is an eigenvector.

We can choose them to be **orthogonal**, and we can rescale them to be **unit vectors** $`\implies`$ **orthonormal**.

### 3.2 Spectral Theorem

Let $`Q = (\mathbf{q}_1 \; \mathbf{q}_2 \; \cdots \; \mathbf{q}_n)`$ where $`\|\mathbf{q}_i\| = 1`$, and $`\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 1 & i = j \\ 0 & i \neq j \end{cases}`$

```math
Q^T Q = \begin{pmatrix} - \mathbf{q}_1^T - \\ - \mathbf{q}_2^T - \\ \vdots \\ - \mathbf{q}_n^T - \end{pmatrix}\begin{pmatrix} | & | & & | \\ \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_n \\ | & | & & | \end{pmatrix} = I
```

For a square matrix $`Q`$: $`Q^T Q = I \implies Q^T = Q^{-1}`$.

Recall that $`AX = X\Lambda`$. For a symmetric matrix $`S`$, we use orthonormal $`Q`$ instead of $`X`$:

```math
SQ = Q\Lambda
```

```math
\boxed{SQQ^T = Q\Lambda Q^{-1} = Q\Lambda Q^T}
```

**Spectral Theorem:** Every real symmetric matrix $`S`$ has the form

```math
S = Q\Lambda Q^T
```

Proof of symmetry: $`S^T = (Q\Lambda Q^T)^T = Q\Lambda^T Q^T = Q\Lambda Q^T = S`$.

If $`S`$ has orthonormal eigenvectors, then $`S`$ is symmetric.

### 3.3 Proof: Symmetric Matrices Have Orthonormal Eigenbasis

**Statement:** If $`S = S^T`$, then $`S`$ has an orthonormal eigenbasis.

**Proof:** Consider two eigenvectors $`\mathbf{u}, \mathbf{v}`$.

i) $`\mathbf{v} \cdot (S\mathbf{u}) = \mathbf{v}^T S\mathbf{u} = \mathbf{v}^T S^T\mathbf{u}`$ (by symmetry) $`= (S\mathbf{v})^T\mathbf{u} = (S\mathbf{v}) \cdot \mathbf{u}`$

ii) Let $`\alpha, \beta`$ be the corresponding eigenvalues: $`S\mathbf{u} = \alpha\mathbf{u}`$, $`S\mathbf{v} = \beta\mathbf{v}`$.

From $`\mathbf{v} \cdot (S\mathbf{u}) = (S\mathbf{v}) \cdot \mathbf{u}`$, we have:

```math
\mathbf{v} \cdot (\alpha\mathbf{u}) = (\beta\mathbf{v}) \cdot \mathbf{u}
```

```math
(\alpha - \beta)\mathbf{v} \cdot \mathbf{u} = 0
```

When $`\alpha \neq \beta`$: $`\mathbf{v} \cdot \mathbf{u} = 0`$.

$`\therefore \mathbf{v} \perp \mathbf{u}`$ — **orthogonal**.

Rescale $`\mathbf{u}, \mathbf{v}`$ to $`\|\mathbf{u}\| = \|\mathbf{v}\| = 1`$. Then $`\mathbf{u}, \mathbf{v}`$ are **orthonormal** eigenvectors. $`\square`$

### 3.4 Positive Definite Matrices: Definitions

**Def.** $`n \times n`$ symmetric real matrix $`S`$ is **positive definite** if

```math
\mathbf{x}^T S\mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}, \; \mathbf{x} \in \mathbb{R}^n
```

(equivalently, $`\forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$)

**Def.** $`n \times n`$ symmetric real matrix $`S`$ is **positive semidefinite** (= non-negative definite) if

```math
\mathbf{x}^T S\mathbf{x} \geq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n
```

**Def.** $`n \times n`$ symmetric real matrix $`S`$ is **negative definite** if

```math
\mathbf{x}^T S\mathbf{x} < 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace
```

**Def.** $`n \times n`$ symmetric real matrix $`S`$ is **negative semidefinite** if

```math
\mathbf{x}^T S\mathbf{x} \leq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n
```

### 3.5 Properties of Positive Definite Matrices

**Property 1.** Positive definite matrix $`S`$ has **all positive eigenvalues**.

**Proof:** $`S`$ is symmetric $`\implies S = Q\Lambda Q^T`$, $`S\mathbf{x} = Q\Lambda Q^T\mathbf{x}`$.

```math
\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\Lambda Q^T\mathbf{x} = (Q^T\mathbf{x})^T \Lambda (Q^T\mathbf{x})
```

Let $`\mathbf{y} = Q^T\mathbf{x}`$:

```math
= \mathbf{y}^T \Lambda \mathbf{y} = \mathbf{y}^T\begin{pmatrix} \lambda_1 y_1 \\ \lambda_2 y_2 \\ \vdots \\ \lambda_n y_n \end{pmatrix} = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2
```

Since $`S`$ is positive definite, $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$.

Choose $`\mathbf{x} = Q\mathbf{e}_i`$ where $`\mathbf{e}_i`$ is the $`i`$-th column of $`I`$. Then $`\mathbf{y} = Q^T\mathbf{x} = Q^T Q\mathbf{e}_i = \mathbf{e}_i`$, which leads to $`\mathbf{x}^T S\mathbf{x} = \lambda_i > 0`$.

Therefore all $`\lambda_1, \lambda_2, \ldots, \lambda_n > 0`$. $`\square`$

**Remark:** All positive eigenvalues does NOT imply the matrix is positive definite.

e.g., $`A = \begin{pmatrix} 1 & -3 \\ 0 & 1 \end{pmatrix} \to \lambda_1 = \lambda_2 = 1 > 0`$, but $`\mathbf{x}^T A\mathbf{x} = (x_1, x_2)\begin{pmatrix} x_1 - 3x_2 \\ x_2 \end{pmatrix} = x_1^2 - 3x_1 x_2 + x_2^2 < 0`$ when $`x_1 = x_2 = 1`$.

(The matrix must be **symmetric** for the equivalence to hold.)

**Positive energy is closely connected to positive eigenvalues ($`\lambda > 0`$):**

If $`S\mathbf{x} = \lambda\mathbf{x}`$, then $`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T\lambda\mathbf{x} = \lambda\mathbf{x}^T\mathbf{x} = \lambda\|\mathbf{x}\|^2`$.

So $`\lambda > 0 \implies \mathbf{x}^T S\mathbf{x} > 0`$ for $`\mathbf{x} \neq \mathbf{0}`$.

**Statement:** If $`\mathbf{x}^T S\mathbf{x} > 0`$ for the eigenvectors of $`S`$, then $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \neq \mathbf{0}`$.

**Proof:** Let $`\mathbf{x} = Q\mathbf{c} = c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n`$ where $`\mathbf{q}_i`$ is the $`i`$-th eigenvector of $`S`$ and an orthonormal vector.

```math
\mathbf{x}^T S\mathbf{x} = (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T S(c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)
```

```math
= (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T(c_1\lambda_1\mathbf{q}_1 + c_2\lambda_2\mathbf{q}_2 + \cdots + c_n\lambda_n\mathbf{q}_n)
```

```math
= c_1^2\lambda_1\mathbf{q}_1^T\mathbf{q}_1 + c_2^2\lambda_2\mathbf{q}_2^T\mathbf{q}_2 + \cdots + c_n^2\lambda_n\mathbf{q}_n^T\mathbf{q}_n > 0
```

if $`\lambda_1, \lambda_2, \ldots, \lambda_n > 0`$. $`\square`$

**Property 2.** If $`S`$ is positive definite, then it is **invertible**, $`\det(S) > 0`$, and $`S^{-1}`$ is positive definite.

**Proof:**

i) $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$

$`\implies S\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$

$`\implies \mathcal{N}(S) = \lbrace\mathbf{0}\rbrace`$

$`\implies S`$ has full rank $`\iff \dim C(S) = n`$

$`\therefore S`$ is **invertible**.

ii) Since $`S`$ is positive definite, $`S`$ has all positive eigenvalues.

Therefore $`0 < \lambda_1\lambda_2\cdots\lambda_n = \det(S)`$.

iii) $`(S^{-1})^T = (S^T)^{-1} = S^{-1}`$, so $`S^{-1}`$ is **symmetric**.

$`\mathbf{x}^T S^{-1}\mathbf{x} = \mathbf{x}^T S^{-1}SS^{-1}\mathbf{x} = (S^{-T}\mathbf{x})^T S(S^{-1}\mathbf{x}) = (S^{-1}\mathbf{x})^T S(S^{-1}\mathbf{x})`$

Let $`\mathbf{z} = S^{-1}\mathbf{x}`$: $`= \mathbf{z}^T S\mathbf{z} > 0`$. $`\square`$

**Property 3.** If $`S`$ is positive definite, then $`S = A^T A`$ with independent columns of $`A`$.

**Proof:** Symmetric matrix $`S = Q\Lambda Q^T`$. Since $`S`$ is positive definite, all eigenvalues are positive. The diagonal matrix $`\Lambda`$ has a square root $`\sqrt{\Lambda}`$:

```math
\Lambda = \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix} = \begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix}\begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix} = \sqrt{\Lambda}\sqrt{\Lambda}
```

Let $`A = Q\sqrt{\Lambda}Q^T`$. Then:

```math
A^T A = (Q\sqrt{\Lambda}Q^T)(Q\sqrt{\Lambda}Q^T) = Q\sqrt{\Lambda}\underbrace{Q^T Q}_{I}\sqrt{\Lambda}Q^T = Q\sqrt{\Lambda}\sqrt{\Lambda}Q^T = Q\Lambda Q^T = S
```

**Remark:** Energy $`= \mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2`$.

$`\|A\mathbf{x}\| > 0`$ if $`A\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$ $`\iff`$ $`A`$ has full rank.

### 3.6 How to Check if a Matrix is Positive Definite

**(1) Positive eigenvalues:** A positive definite matrix $`S`$ has all positive eigenvalues.

e.g., $`S = \mathbf{u}\mathbf{v}^T`$, rank 1 matrix is NOT positive definite.

$`A = \begin{pmatrix} a \\ b \end{pmatrix}(1 \; k) = \begin{pmatrix} a & ka \\ b & kb \end{pmatrix}`$

Rank 1 matrix $`\to`$ 1 LI vector $`\to`$ $`(n-1)`$ dependent vectors.

$`A\mathbf{u} = \mathbf{u}\mathbf{v}^T\mathbf{u} = (\mathbf{v}^T\mathbf{u})\mathbf{u}`$, so $`\lambda = \mathbf{v}^T\mathbf{u} = (1, k)\begin{pmatrix} a \\ b \end{pmatrix} = a + kb`$.

Let $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_{n-1}`$ be basis vectors of $`\mathcal{N}(A)`$. $`A\mathbf{w}_i = \mathbf{0} = 0 \cdot \mathbf{w}_i`$, so $`\lambda = 0`$ for $`i = 1, 2, \ldots, n-1`$.

$`\therefore (n-1)`$ zero eigenvalues.

**(2) Positive energy test:** $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$.

**(3) Positive energy test for $`S = A^T A`$:**

```math
\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2
```

e.g., $`S = A^T A`$ where $`A = \begin{pmatrix} 1 & 1 & 2 \\ 1 & 2 & 3 \end{pmatrix}`$ with columns $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3`$.

Note: $`\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2`$, so $`\text{rank}(A) = 2 < 3`$, $`\dim C(A^T) = 2`$, $`\dim \mathcal{N}(A) = 1`$.

Since $`\dim \mathcal{N}(A) \neq 0`$, $`\mathbf{x}^T S\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0`$. $`S`$ is NOT positive definite. $`S`$ is **positive semidefinite**.

**(4) Determinant test:** Check if $`\det(S) > 0`$ — but more precisely, check all **upper-left determinants**.

```math
S = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & -1 & 2 & -1 \\ & & -1 & 2 \end{pmatrix}
```

$`D_1 = |2| = 2 > 0`$

$`D_2 = \begin{vmatrix} 2 & -1 \\ -1 & 2 \end{vmatrix} = 3 > 0`$

$`D_3 = \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{vmatrix} = 8 - 4 = 4 > 0`$

$`D_4 = |S| = 2D_3 + \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & -1 & -1 \end{vmatrix} = 8 - 4 + 1 = 5 > 0`$

All upper-left determinants are positive, so $`S`$ is positive definite.

**(5) Pivot test:** Check if the **pivots** are positive.

For $`S = \begin{pmatrix} 2 & -1 \\ -1 & 2 & -1 \\ & -1 & 2 \end{pmatrix}`$:

After elimination: pivots are $`2, \frac{3}{2}, \frac{4}{3}`$ — all positive.

$`S = LU`$, $`S = LDL^T`$, $`S = A^T A`$.

Leading determinants are closely related to pivots: $`D_2/D_1`$, $`D_3/D_1`$, etc.

**For a $`2 \times 2`$ SPD matrix** $`S = \begin{pmatrix} a & b \\ b & d \end{pmatrix}`$:

- Determinant test: $`D_1 = a > 0`$, $`D_2 = ad - b^2 > 0`$.
- Pivot test: $`d_1 = a > 0`$, $`d_2 = \frac{ad - b^2}{a} > 0`$.
- Eigenvalues: $`\lambda_1 > 0`$, $`\lambda_2 > 0`$.
- Energy: $`(x \; y)\begin{pmatrix} a & b \\ b & d \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = ax^2 + bxy + byx + dy^2 = ax^2 + 2bxy + dy^2 > 0`$.

### 3.7 Worked Examples: Positive Definite and Semidefinite

**Example:** $`S = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix} \implies \lambda_1 = 2 > 0, \lambda_2 = 6 > 0`$

$`S\mathbf{x} = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix}`$

$`\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix} = 2x_1^2 + 6x_2^2 > 0`$. $`\therefore S`$ is positive definite.

**Example:** $`S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}`$

$`|S - \lambda I| = \begin{vmatrix} 5 - \lambda & 4 \\ 4 & 5 - \lambda \end{vmatrix} = \lambda^2 - 10\lambda + 25 - 16 = (\lambda - 9)(\lambda - 1) = 0`$

$`\therefore \lambda_1 = 9 > 0`$, $`\lambda_2 = 1 > 0`$.

i) $`(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = x_2`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 + x_2 = 0`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) $`Q = (\mathbf{x}_1 \; \mathbf{x}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}`$, $`Q^{-1} = Q = Q^T`$

iv) $`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T(Q\Lambda Q^T)\mathbf{x}`$, let $`\mathbf{y} = Q^T\mathbf{x}`$:

$`= \lambda_1 y_1^2 + \lambda_2 y_2^2 = 9y_1^2 + y_2^2 > 0`$

$`\therefore S`$ is positive definite.

**Alternative approach (energy test):**

$`\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = (x_1, x_2)\begin{pmatrix} 5x_1 + 4x_2 \\ 4x_1 + 5x_2 \end{pmatrix}`$

$`= 5x_1^2 + 4x_1x_2 + 4x_1x_2 + 5x_2^2 = x_1^2 + 4x_1^2 + 8x_1x_2 + 4x_2^2 + x_1^2`$

$`= x_1^2 + \underbrace{4(x_1^2 + 2x_1x_2 + x_2^2)}_{4(x_1 + x_2)^2} > 0`$

**Example:** $`S = \begin{pmatrix} 4 & 5 \\ 5 & 4 \end{pmatrix}`$

$`|S - \lambda I| = \lambda^2 - 8\lambda + 16 - 25 = \lambda^2 - 8\lambda - 9 = (\lambda - 9)(\lambda + 1) = 0`$

$`\therefore \lambda_1 = 9`$ and $`\lambda_2 = -1`$.

Since $`\lambda_2 < 0`$, $`S`$ is **NOT** positive definite.

### 3.8 The Ellipse and Quadratic Forms

**The Ellipse** $`ax^2 + 2bxy + cy^2 = 1`$.

**Example 1:** Consider an ellipse $`5x^2 + 8xy + 5y^2 = 1`$.

```math
(x \; y)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 1 \implies \mathbf{x}^T S\mathbf{x} = 1
```

$`S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}`$, $`\lambda^2 - 10\lambda + 9 = (\lambda - 9)(\lambda - 1) = 0`$, so $`\lambda = 9, 1`$.

i) $`(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) $`\mathbf{q}_1 = \frac{1}{\sqrt{2}}\mathbf{x}_1`$, $`\mathbf{q}_2 = \frac{1}{\sqrt{2}}\mathbf{x}_2`$

$`Q = (\mathbf{q}_1 \; \mathbf{q}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}`$, $`Q^{-1} = Q^T = Q`$

$`S = Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T`$

iv) $`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T\mathbf{x}`$

Let $`\mathbf{z} = Q^T\mathbf{x}`$: $`= \mathbf{z}^T\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}\mathbf{z} = (z_1, z_2)\begin{pmatrix} 9z_1 \\ z_2 \end{pmatrix} = 9z_1^2 + z_2^2`$

$`\mathbf{z} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix} x_1 + x_2 \\ x_1 - x_2 \end{pmatrix}`$

$`\therefore z_1 = \frac{x_1 + x_2}{\sqrt{2}}, \; z_2 = \frac{x_1 - x_2}{\sqrt{2}}`$

$`\mathbf{x}^T S\mathbf{x} = 9z_1^2 + z_2^2 = 1`$

When $`z_2 = 0`$: $`z_1^2 = \frac{1}{9}`$, so $`z_1 = \pm \frac{1}{3}`$ (semi-axis along $`\mathbf{q}_1`$).

When $`z_1 = 0`$: $`z_2^2 = 1`$, so $`z_2 = \pm 1`$ (semi-axis along $`\mathbf{q}_2`$).

**Example 2:** Positive semidefinite

```math
T = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}\begin{pmatrix} 3 & 1 \end{pmatrix} = A^T A
```

Determinants: $`D_1 = 9 > 0`$, $`D_2 = 9 - 9 = 0 \to`$ No inverse.

$`\lambda^2 - \text{trace}(T)\lambda + \det(T) = \lambda^2 - 10\lambda = 0`$, so $`\lambda_1 = 10, \lambda_2 = 0`$.

$`(T - 10I)\mathbf{x}_1 = \begin{pmatrix} -1 & 3 \\ 3 & -9 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 3 \\ 1 \end{pmatrix}`$

$`(T - 0I)\mathbf{x}_2 = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -3 \end{pmatrix}`$

Energy: $`ax^2 + 2bxy + dy^2 = 9x^2 + 6xy + y^2 = (3x + y)^2 \geq 0`$.

$`E(x,y) = 1`$ is a band: $`3x + y = \pm 1`$. The graph of energy $`E(x,y)`$ is a valley. Axes of the band along eigenvectors of $`T`$.

### 3.9 Positive Definite Matrices and Minimum Problems

Consider energy $`E = \mathbf{x}^T S\mathbf{x} = \begin{pmatrix} x \\ y \end{pmatrix}^T\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 5x^2 + 8xy + 5y^2 > 0`$

(Bowl shape, **convex**)

$`E(x,y) = 0`$ when $`x = y = 0`$. This connects to **minimum problems**.

The matrix of second derivatives is positive definite at all points.

1st derivatives: $`\frac{\partial E}{\partial x} = 10x + 8y`$, $`\frac{\partial E}{\partial y} = 8x + 10y`$.

2nd derivatives: $`\frac{\partial^2 E}{\partial x^2} = 10 > 0`$, $`\frac{\partial^2 E}{\partial x \partial y} = 8 > 0`$, $`\frac{\partial^2 E}{\partial y^2} = 10 > 0`$.

```math
\nabla E = \left(\frac{\partial E}{\partial x}, \frac{\partial E}{\partial y}\right) \quad \xrightarrow{(x,y) = (0,0)} \quad \nabla E = (0, 0)
```

```math
J(\nabla E^T) = \begin{pmatrix} \frac{\partial^2 E}{\partial x^2} & \frac{\partial^2 E}{\partial y \partial x} \\ \frac{\partial^2 E}{\partial x \partial y} & \frac{\partial^2 E}{\partial y^2} \end{pmatrix} = H \quad \text{(Hessian matrix)}
```

$`H_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}`$

If we define the energy by $`\frac{1}{2}\mathbf{x}^T S\mathbf{x}`$, then $`H = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix} = S`$.

e.g., $`E = x^2 + y^2`$ has minimum at $`x = y = 0`$.

What about $`f = e^{x^2 + y^2}`$?

$`\nabla f = \left(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}\right)`$ where $`\frac{\partial f}{\partial x} = e^{x^2+y^2} \cdot 2x`$, $`\frac{\partial f}{\partial y} = e^{x^2+y^2} \cdot 2y`$.

$`J(\nabla f)^T = \begin{pmatrix} (2 + 4x^2)e^{x^2+y^2} & 4xy \, e^{x^2+y^2} \\ 4xy \, e^{x^2+y^2} & (2 + 4y^2)e^{x^2+y^2} \end{pmatrix} = S`$

$`S`$ is positive definite, so $`f`$ is **strictly convex**.

**Example 1 (Positive Definite):**

$`S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix}`$

Determinants: $`D_1 = 9 > 0`$, $`D_2 = 27 - 9 = 18 > 0`$.

Pivots: $`\begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} \xrightarrow{R_2 - \frac{1}{3}R_1} \begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix}`$. Pivots: 9, 2 (both positive).

Energy: $`E(x,y) = 9x^2 + 6xy + 3y^2 = (3x + y)^2 + 2y^2 > 0`$. (Strictly convex function.)

Trace and determinants: $`\text{trace}(S) = 12`$, $`\det(S) = 18`$.

$`\lambda^2 - 12\lambda + 18 = 0 \implies \lambda = 6 \pm \sqrt{36 - 18} = 6 \pm 3\sqrt{2}`$

i) $`\lambda_1 = 6 + 3\sqrt{2}`$: $`(S - \lambda_1 I)\mathbf{x} = \begin{pmatrix} 3 - 3\sqrt{2} & 3 \\ 3 & 3 - 3\sqrt{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_2 = (-1 + \sqrt{2})x_1`$, choose $`x_1 = 1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ -1 + \sqrt{2} \end{pmatrix}`$

ii) $`\lambda_2 = 6 - 3\sqrt{2}`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 - \sqrt{2} \end{pmatrix}`$

**Decomposition:**

$`S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} = LU = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix} = LDL^T = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & \frac{1}{3} \\ 0 & 1 \end{pmatrix}`$

$`= (L\sqrt{D})(\sqrt{D}L^T) = A^T A`$ where $`A = \begin{pmatrix} 3 & 1 \\ 0 & \sqrt{2} \end{pmatrix}`$

### 3.10 Positive Semidefinite Matrices

$`\mathbf{x}^T S\mathbf{x} \geq 0`$

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}
```

$`\det(S) = 4 - 4 = 0`$ — singular matrix.

$`\text{trace}(S) = 1 + 4 = 5`$. $`\lambda^2 - 5\lambda + 0 = \lambda(\lambda - 5) = 0`$, so $`\lambda_1 = 5, \lambda_2 = 0`$.

i) $`(S - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_2 = 2x_1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$

ii) $`S\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 + 2x_2 = 0`$: $`\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}`$

$`E_{21}S = U`$: $`S = E_{21}^{-1}U = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = LDL^T`$

$`= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}`$

$`= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}`$

$`S`$ is decomposed into $`A^T A`$ with **dependent** columns in $`A`$.

$`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0`$

### 3.11 Congruent Matrices

If $`S`$ is positive semidefinite, so is every matrix $`A^T SA`$:

```math
\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} \geq 0 \quad \forall \mathbf{x}
```

where $`\mathbf{y} = A\mathbf{x}`$. $`\square`$

Suppose $`\mathbf{x}^T S\mathbf{x} > 0`$ and $`A\mathbf{x} \neq \mathbf{0}`$ for all $`\mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$. Then $`A^T SA`$ is **positive definite**.

```math
\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} > 0 \quad \forall \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace \quad \square
```

The matrix $`A^T SA`$ is called **"congruent"** to $`S`$.

If $`S^T = S`$ has $`P`$ positive eigenvalues, $`N`$ negative eigenvalues, and $`Z`$ zero eigenvalues, then the same is true for $`A^T SA`$ provided $`A`$ is invertible.

**Proof of skew-symmetric eigenvalue property:**

$`A = -\bar{A}^T`$. $`A\mathbf{x} = \lambda\mathbf{x}`$, $`\overline{A\mathbf{x}} = \bar{\lambda}\bar{\mathbf{x}}`$.

$`\bar{\mathbf{x}}^T A\mathbf{x} = \bar{\mathbf{x}}^T(-\bar{A}^T)\mathbf{x} = -(\bar{A}\bar{\mathbf{x}})^T\mathbf{x}`$

$`= \bar{\mathbf{x}}^T\lambda\mathbf{x}`$ and $`= -(\bar{\lambda}\bar{\mathbf{x}})^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}`$

$`\lambda\bar{\mathbf{x}}^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}`$

$`(\lambda + \bar{\lambda})\bar{\mathbf{x}}^T\mathbf{x} = 0`$

$`\therefore \lambda + \bar{\lambda} = (r + is) + (r - is) = 2r = 0`$

Real part of $`\lambda`$ is **zero**. Antisymmetric matrix has zero or **pure imaginary** eigenvalues.

### 3.12 Optimization and Machine Learning

**Gradient descent** to minimize $`f(x)`$.

i) $`f(x) = x^2 + 4x + 4 \implies f'(x) = 2x + 4`$

ii) $`x_0 = 10`$, $`\eta = 0.1`$ (learning rate), stop criterion $`|f'(x)| < 0.01`$

iii) $`x_{k+1} = x_k - \eta f'(x_k)`$ — **steepest direction**.

Iterate the approximation until $`x_k \to x^*`$ (the minimum point).

For $`f(\mathbf{x})`$: $`\mathbf{x}_{k+1} = \mathbf{x}_k - \eta \nabla f(\mathbf{x}_k)`$

If $`\mathbf{x}_k \to \mathbf{x}^*`$, then $`\nabla f`$ at $`\mathbf{x}^*`$ is zero $`\implies \mathbf{x}_{k+1} \approx \mathbf{x}_k`$.

If $`J(\nabla f)^T = S`$ is positive definite, then $`f`$ is **convex** — easy to find $`\mathbf{x}^*`$.

---

<br>
