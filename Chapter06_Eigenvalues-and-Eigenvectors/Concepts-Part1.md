# Chapter 6 Lecture — Eigenvalues and Eigenvectors

> **Last Updated:** 2026-06-06
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 6

> **Prerequisites**: [Linear Algebra] Determinants, matrix operations (Ch 1-5).
>
> **Learning Objectives**:
> 1. Define eigenvalues and eigenvectors and compute them
> 2. Diagonalize matrices using eigendecomposition
> 3. Analyze stability and dynamics using eigenvalues

> **📑 This document is split into 3 parts.**
>
> **Part 1** · [Part 2](Concepts-Part2.md) · [Part 3](Concepts-Part3.md)

---

<br>

## Table of Contents

- [1. Introduction to Eigenvalues (6.1)](#1-introduction-to-eigenvalues-61)
  - [1.1 Definition of Eigenvalues and Eigenvectors](#11-definition-of-eigenvalues-and-eigenvectors)
  - [1.2 Powers of A and Eigenvalues](#12-powers-of-a-and-eigenvalues)
  - [1.3 Properties of Eigenvalues](#13-properties-of-eigenvalues)
  - [1.4 The Equation for Eigenvalues](#14-the-equation-for-eigenvalues)
  - [1.5 Determinant and Trace](#15-determinant-and-trace)
  - [1.6 Worked Examples](#16-worked-examples)
  - [1.7 Imaginary Eigenvalues](#17-imaginary-eigenvalues)
  - [1.8 Rotation Matrix Eigenvalues](#18-rotation-matrix-eigenvalues)
  - [1.9 Eigenvalues of AB and A+B](#19-eigenvalues-of-ab-and-ab)
- [2. Diagonalizing a Matrix (6.2)](#2-diagonalizing-a-matrix-62)
  - [2.1 Key Facts of Diagonalization](#21-key-facts-of-diagonalization)
  - [2.2 Diagonalization Procedure](#22-diagonalization-procedure)
  - [2.3 Worked Example: Diagonalization](#23-worked-example-diagonalization)
  - [2.4 Remarks on Diagonalization](#24-remarks-on-diagonalization)
  - [2.5 Proof: Eigenvectors for Distinct Eigenvalues are LI](#25-proof-eigenvectors-for-distinct-eigenvalues-are-li)
  - [2.6 Powers of A (Markov Matrix Example)](#26-powers-of-a-markov-matrix-example)
  - [2.7 Similar Matrices](#27-similar-matrices)
  - [2.8 Matrix Powers and Fibonacci Numbers](#28-matrix-powers-and-fibonacci-numbers)
  - [2.9 Nondiagonalizable Matrices and Multiplicity](#29-nondiagonalizable-matrices-and-multiplicity)
- [3. Symmetric Positive Definite Matrices (6.3)](Concepts-Part2.md#3-symmetric-positive-definite-matrices-63)
  - [3.1 Symmetric Matrices: Key Properties](Concepts-Part2.md#31-symmetric-matrices-key-properties)
  - [3.2 Spectral Theorem](Concepts-Part2.md#32-spectral-theorem)
  - [3.3 Proof: Symmetric Matrices Have Orthonormal Eigenbasis](Concepts-Part2.md#33-proof-symmetric-matrices-have-orthonormal-eigenbasis)
  - [3.4 Positive Definite Matrices: Definitions](Concepts-Part2.md#34-positive-definite-matrices-definitions)
  - [3.5 Properties of Positive Definite Matrices](Concepts-Part2.md#35-properties-of-positive-definite-matrices)
  - [3.6 How to Check if a Matrix is Positive Definite](Concepts-Part2.md#36-how-to-check-if-a-matrix-is-positive-definite)
  - [3.7 Worked Examples: Positive Definite and Semidefinite](Concepts-Part2.md#37-worked-examples-positive-definite-and-semidefinite)
  - [3.8 The Ellipse and Quadratic Forms](Concepts-Part2.md#38-the-ellipse-and-quadratic-forms)
  - [3.9 Positive Definite Matrices and Minimum Problems](Concepts-Part2.md#39-positive-definite-matrices-and-minimum-problems)
  - [3.10 Positive Semidefinite Matrices](Concepts-Part2.md#310-positive-semidefinite-matrices)
  - [3.11 Congruent Matrices](Concepts-Part2.md#311-congruent-matrices)
  - [3.12 Optimization and Machine Learning](Concepts-Part2.md#312-optimization-and-machine-learning)
- [4. Solving Linear Differential Equations (6.5)](Concepts-Part3.md#4-solving-linear-differential-equations-65)
  - [4.1 Key Facts](Concepts-Part3.md#41-key-facts)
  - [4.2 Scalar ODE Review](Concepts-Part3.md#42-scalar-ode-review)
  - [4.3 Solution of du/dt = Au](Concepts-Part3.md#43-solution-of-dudt--au)
  - [4.4 General n x n Solution Procedure](Concepts-Part3.md#44-general-n-x-n-solution-procedure)
  - [4.5 Exponential of a Matrix](Concepts-Part3.md#45-exponential-of-a-matrix)
  - [4.6 Second Order Equations](Concepts-Part3.md#46-second-order-equations)
  - [4.7 Stability of 2 by 2 Matrices](Concepts-Part3.md#47-stability-of-2-by-2-matrices)
  - [4.8 Worked Examples](Concepts-Part3.md#48-worked-examples)
- [Summary](Concepts-Part3.md#summary)

---

<br>

## 1. Introduction to Eigenvalues (6.1)

### 1.1 Definition of Eigenvalues and Eigenvectors

When $`A`$ acts on $`\mathbf{x}`$, it only stretches or compresses the vector $`\mathbf{x}`$ by $`\lambda`$, without changing its direction.

```math
A\mathbf{x} = \lambda \mathbf{x}
```

- $`\mathbf{x}`$ is a **non-zero vector**, known as an **eigenvector**.
- $`\lambda`$ is the **eigenvalue** corresponding to $`\mathbf{x}`$.

**Formal statement:**

1. If $`A\mathbf{x} = \lambda\mathbf{x}`$, then $`\mathbf{x} \neq \mathbf{0}`$ is an eigenvector of $`A`$, and the number $`\lambda`$ is the eigenvalue.

2. $`A^n \mathbf{x} = \lambda^n \mathbf{x}`$ for every $`n`$; $`(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}`$; and $`A^{-1}\mathbf{x} = \frac{1}{\lambda}\mathbf{x} = \lambda^{-1}\mathbf{x}`$ if $`\lambda \neq 0`$.

**Proof** (for $`A^{-1}`$):

```math
A\mathbf{x} = \lambda\mathbf{x} \implies \mathbf{x} = A^{-1}(A\mathbf{x}) = A^{-1}(\lambda\mathbf{x}) = \lambda A^{-1}\mathbf{x} \quad \square
```

3. $`(A - \lambda I)\mathbf{x} = \mathbf{0} \implies \det(A - \lambda I) = 0`$. This equation produces $`n`$ $`\lambda`$'s.

4. Let $`A \in \mathbb{R}^{n \times n}`$:

```math
\det(A) = \lambda_1 \lambda_2 \cdots \lambda_n; \quad \text{trace}(A) = \lambda_1 + \lambda_2 + \cdots + \lambda_n
```

5. **Projection matrix** $`P`$ has $`\lambda = 1`$ or $`\lambda = 0`$. A square matrix $`P`$ is called a "projection matrix" if $`P^2 = P`$.

### 1.2 Powers of A and Eigenvalues

What happens if we multiply $`A`$ to the relation?

```math
A(A\mathbf{x}) = A(\lambda\mathbf{x}) \implies A^2\mathbf{x} = \lambda A\mathbf{x} = \lambda^2 \mathbf{x}
```

If we keep multiplying $`\mathbf{x}`$ by $`A`$:

```math
A^3\mathbf{x} = \lambda^3\mathbf{x}, \quad A^k\mathbf{x} = \lambda^k\mathbf{x}, \quad \ldots, \quad A^{100}\mathbf{x} = \lambda^{100}\mathbf{x}
```

```math
\boxed{A^k \mathbf{x} = \lambda^k \mathbf{x}}
```

**Behavior of $`A^k`$:**

- If $`|\lambda_i| < 1`$ for $`i = 1, 2, \ldots, n`$, then $`A^k`$ will eventually approach **zero**.
- If any $`|\lambda_i| > 1`$, then $`A^k`$ will eventually **grow**.
- If $`\lambda = 1`$, then the system state does not grow or decay over time, but rather stays **constant**.

When $`A\mathbf{x} = \mathbf{x}`$, then $`\mathbf{x}`$ is a **fixed point** or **equilibrium point** of a system, i.e., the system remains steady in $`\mathbf{x}`$ direction. The system reaches a **steady state**: $`A^k \mathbf{x} = \mathbf{x}`$.

### 1.3 Properties of Eigenvalues

An $`A \in \mathbb{R}^{n \times n}`$ matrix has $`n`$ eigenvalues. We can find the eigenvalues by solving the **characteristic polynomial**:

```math
\det(A - \lambda I) = 0
```

**Properties of a matrix significantly influence its eigenvalues:**

**(1)** The **trace** of a matrix $`A`$ is equal to the **sum** of its eigenvalues.

e.g., $`A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}`$, $`\text{trace}(A) = 1 + 4 = 5`$

$`\det(A - \lambda I) = 0 \implies \lambda_1 = 0, \lambda_2 = 5`$, and $`\lambda_1 + \lambda_2 = 5`$.

**(2)** The **determinant** of $`A`$ is the **product** of its eigenvalues.

e.g., $`\det(A) = 1 \cdot 4 - 2 \cdot 2 = 0`$, and $`\lambda_1 \cdot \lambda_2 = 0 \cdot 5 = 0`$.

**(3)** A **symmetric** matrix (i.e., $`A = A^T`$) has only **real eigenvalues**.

e.g., $`A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0`$, so $`\lambda = \pm 1`$.

**(4)** A **symmetric positive definite (SPD)** matrix has **real and positive** eigenvalues. SPD matrix is crucial in optimization because it guarantees that a quadratic form $`f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A\mathbf{x}`$ has a **unique minimum**. (See Section 6.3)

**(5)** For a **diagonal** matrix, the eigenvalues are the **diagonal elements**.

e.g., $`A = \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} a - \lambda & 0 \\ 0 & b - \lambda \end{vmatrix} = (a - \lambda)(b - \lambda) = 0`$, so $`\lambda_1 = a, \lambda_2 = b`$.

For a **triangular** matrix, the eigenvalues are the **diagonal elements**.

e.g., $`A = \begin{pmatrix} a & b \\ 0 & c \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} a - \lambda & b \\ 0 & c - \lambda \end{vmatrix} = (a - \lambda)(c - \lambda) = 0`$, so $`\lambda_1 = a, \lambda_2 = c`$.

**(6)** A **skew-symmetric** matrix (i.e., $`A = -A^T`$) has **purely imaginary** eigenvalues or **zero** eigenvalues.

e.g., $`A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ -1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0`$, so $`\lambda = \pm i`$.

### 1.4 The Equation for Eigenvalues

From $`A\mathbf{x} = \lambda\mathbf{x}`$:

```math
A\mathbf{x} - \lambda\mathbf{x} = (A - \lambda I)\mathbf{x} = \mathbf{0}
```

The eigenvectors make up the **nullspace** of $`A - \lambda I`$.

When we know an eigenvalue $`\lambda`$, we find an eigenvector by solving $`(A - \lambda I)\mathbf{x} = \mathbf{0}`$.

If $`\mathbf{x}`$ is a nonzero vector, then $`A - \lambda I`$ is **singular**.

```math
\boxed{\det(A - \lambda I) = 0}
```

This characteristic polynomial ($`\det(A - \lambda I) = 0`$) involves only $`\lambda`$, which appears all along the main diagonal of $`A - \lambda I`$. The determinant includes $`(-\lambda)^n`$.

- The equation has $`n`$ solutions $`\lambda_1`$ to $`\lambda_n`$.
- $`A_{n \times n}`$ has $`n`$ eigenvalues.

**Proof (Determinant = Product of Eigenvalues):**

```math
A\mathbf{x} = \lambda\mathbf{x} \implies \det(A - \lambda I) = 0
```

The characteristic polynomial for $`\lambda`$:

```math
\det(\lambda I - A) = 0 \iff \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0 = 0
```

```math
\iff (\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n) = 0
```

Take $`\lambda = 0`$:

```math
\det(0 \cdot I - A) = \det(-A) = (-1)^n \det(A)
```

```math
(-\lambda_1)(-\lambda_2)\cdots(-\lambda_n) = (-1)^n \lambda_1\lambda_2\cdots\lambda_n
```

```math
\therefore \det(A) = \lambda_1\lambda_2\cdots\lambda_n \quad \square
```

Also, from expanding $`(\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n)`$:

```math
= \lambda^n - (\lambda_1 + \lambda_2 + \cdots + \lambda_n)\lambda^{n-1} + \cdots + \det(A) = 0
```

The coefficient of $`\lambda^{n-1}`$ equals $`\text{trace}(A) = a_{11} + a_{22} + \cdots + a_{nn}`$ (skip proof).

**Derivative of a determinant of a matrix:**

i) $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$, $`\det(A) = ad - bc`$.

Let $`a = a(x), b = b(x), c = c(x), d = d(x)`$.

```math
\frac{d}{dx}\det(A) = \frac{d}{dx}(a(x)d(x) - b(x)c(x)) = a'd + ad' - b'c - bc'
```

```math
= \underbrace{a'd - b'c} + \underbrace{ad' - bc'}
```

```math
\begin{vmatrix} a & b \\ c & d \end{vmatrix}' = \begin{vmatrix} a' & b' \\ c & d \end{vmatrix} + \begin{vmatrix} a & b \\ c' & d' \end{vmatrix}
```

To differentiate a determinant, we differentiate **one row (or column) at a time**, keeping others unchanged.

For $`n \times n`$ matrix: writing rows as $`\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_n`$:

```math
\begin{vmatrix} a_{11} & \cdots & a_{1n} \\ \vdots & & \vdots \\ a_{n1} & \cdots & a_{nn} \end{vmatrix}' = \sum_{i=1}^{n} \begin{vmatrix} \mathbf{r}_1 \\ \vdots \\ \mathbf{r}_i' \\ \vdots \\ \mathbf{r}_n \end{vmatrix}
```

### 1.5 Determinant and Trace

**Observation 1:** Elimination does **not** preserve eigenvalues.

```math
A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}
```

$`\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 2 & 6 - \lambda \end{vmatrix} = \lambda^2 - 7\lambda = 0`$, so $`\lambda_1 = 7, \lambda_2 = 0`$.

$`\det(A_0 - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 - \lambda = 0`$, so $`\lambda_1 = 1, \lambda_2 = 0`$. (Different!)

**Observation 2:** The product $`\lambda_1 \cdot \lambda_2`$ and the sum $`\lambda_1 + \lambda_2`$ can be found from the matrix.

For $`A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}`$:

- $`\lambda_1 \lambda_2 = 7 \cdot 0 = 0 = \det(A) = 6 - 6 = 0`$
- $`\lambda_1 + \lambda_2 = 7 + 0 = 7 = \text{trace}(A) = 1 + 6 = 7`$

The product $`\lambda_1 \lambda_2 \cdots \lambda_n`$ of $`n`$ eigenvalues equals the **determinant** of $`A`$.

The sum $`\lambda_1 + \lambda_2 + \cdots + \lambda_n`$ equals the sum of the $`n`$ diagonal entries = **trace** of $`A`$.

### 1.6 Worked Examples

**Example 1:** Markov Matrix

```math
A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix}
```

```math
\det(A - \lambda I) = \begin{vmatrix} 0.8 - \lambda & 0.3 \\ 0.2 & 0.7 - \lambda \end{vmatrix} = 0.56 - 1.5\lambda + \lambda^2 - 0.06 = \lambda^2 - \frac{3}{2}\lambda + \frac{1}{2} = (\lambda - \frac{1}{2})(\lambda - 1) = 0
```

```math
\therefore \lambda_1 = 1, \quad \lambda_2 = \frac{1}{2}
```

This tells us that $`A - \lambda_1 I = A - I`$ and $`A - \lambda_2 I = A - \frac{1}{2}I`$ are **NOT invertible**.

The eigenvectors $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$ satisfy $`(A - I)\mathbf{x}_1 = \mathbf{0}`$ and $`(A - \frac{1}{2}I)\mathbf{x}_2 = \mathbf{0}`$.

That is, $`\mathbf{x}_1 \in \mathcal{N}(A - I)`$ and $`\mathbf{x}_2 \in \mathcal{N}(A - \frac{1}{2}I)`$.

i) $`(A - I)\mathbf{x}_1 = \begin{pmatrix} -0.2 & 0.3 \\ 0.2 & -0.3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`-2x_1 + 3x_2 = 0`$. Choose $`x_1 = 3`$, then $`x_2 = 2`$: $`\mathbf{x}_1 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}`$

ii) $`(A - \frac{1}{2}I)\mathbf{x}_2 = \begin{pmatrix} 0.3 & 0.3 \\ 0.2 & 0.2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_1 + x_2 = 0`$, so $`x_1 = -x_2`$. Choose $`x_1 = 1`$, then $`x_2 = -1`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

Multiply $`A`$ to $`\mathbf{x}_1`$: $`A\mathbf{x}_1 = \mathbf{x}_1`$, $`A^2\mathbf{x}_1 = \mathbf{x}_1`$, ..., $`A^{100}\mathbf{x}_1 = \mathbf{x}_1`$.

Multiply $`A`$ to $`\mathbf{x}_2`$: $`A\mathbf{x}_2 = \frac{1}{2}\mathbf{x}_2`$, $`A^2\mathbf{x}_2 = (\frac{1}{2})^2\mathbf{x}_2`$, ..., $`A^{100}\mathbf{x}_2 = (\frac{1}{2})^{100}\mathbf{x}_2`$.

Both $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$ stay in their own directions.

An eigenvector $`\mathbf{x}`$ of $`A`$ is also an eigenvector of every $`A^n`$ because of $`A^n\mathbf{x} = \lambda^n\mathbf{x}`$.

Eigenvectors $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$ span $`\mathbb{R}^2`$.

Any vector $`\mathbf{x}`$ is a linear combination of $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$: $`\mathbf{x} = c\mathbf{x}_1 + d\mathbf{x}_2 = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix}`$.

```math
A\mathbf{x} = cA\mathbf{x}_1 + dA\mathbf{x}_2 = c\mathbf{x}_1 + d\tfrac{1}{2}\mathbf{x}_2
```

```math
A^2\mathbf{x} = c\mathbf{x}_1 + d(\tfrac{1}{2})^2\mathbf{x}_2
```

```math
A^{100}\mathbf{x} = c\,\underbrace{\mathbf{x}_1}_{\text{steady state}} + d(\tfrac{1}{2})^{100}\underbrace{\mathbf{x}_2}_{\text{decaying mode}} \approx c\,\mathbf{x}_1
```

Let $`\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$:

```math
\begin{pmatrix} 1 \\ 0 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix}
```

```math
\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 & 1 \\ 2 & -3 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 0.2 \\ 0.4 \end{pmatrix}
```

```math
= 0.2\begin{pmatrix} 3 \\ 2 \end{pmatrix} + 0.4\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix} + \begin{pmatrix} 0.4 \\ -0.4 \end{pmatrix}
```

```math
A^{100}\begin{pmatrix} 1 \\ 0 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}
```

Similarly for $`\mathbf{x} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$: $`A^{100}\begin{pmatrix} 0 \\ 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}`$.

```math
A^{100}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_1)
```

The higher the power of $`A`$, the more closely its columns approach the **steady state**.

**Example 2:** Projection Matrix

```math
P = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}
```

has eigenvalues $`\lambda = 1`$ and $`\lambda = 0`$.

```math
\det(P - \lambda I) = \begin{vmatrix} \frac{1}{2} - \lambda & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} - \lambda \end{vmatrix} = (\frac{1}{2})^2 \begin{vmatrix} 1 - 2\lambda & 1 \\ 1 & 1 - 2\lambda \end{vmatrix} = \frac{1}{4}[(1 - 2\lambda)^2 - 1] = \frac{1}{4}(4\lambda^2 - 4\lambda) = \lambda(\lambda - 1) = 0
```

```math
\therefore \lambda_1 = 1, \; \lambda_2 = 0
```

i) $`\lambda_1 = 1`$: $`(P - I)\mathbf{x}_1 = \begin{pmatrix} -\frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & -\frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = x_2`$, choose $`x_1 = 1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`\lambda_2 = 0`$: $`P\mathbf{x}_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 + x_2 = 0`$, $`x_1 = -x_2`$, choose $`x_1 = 1`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) $`P\mathbf{x}_1 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \mathbf{x}_1`$

$`P\mathbf{x}_2 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}`$

$`\therefore \mathbf{x}_1 \in C(P)`$: the column space projects onto itself.

$`\mathbf{x}_2 \in \mathcal{N}(P)`$

iv) Let $`\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}`$

$`P\mathbf{w} = P\begin{pmatrix} 1 \\ -1 \end{pmatrix} + P\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \mathbf{0} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}`$

The only eigenvalues of a projection matrix are **0 and 1**.

**Example 3:** Exchange Matrix

```math
E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
```

has eigenvalues 1 and $`-1`$.

```math
\det(E - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = (\lambda - 1)(\lambda + 1) = 0
```

```math
\therefore \lambda_1 = 1, \; \lambda_2 = -1
```

i) $`\lambda_1 = 1`$: $`(E - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`-x_1 + x_2 = 0`$, $`x_1 = x_2`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`\lambda_2 = -1`$: $`(E + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = -x_2`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) The eigenvectors $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$ are the same as for $`P`$.

```math
E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} - \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 2P - I
```

When a matrix is shifted by $`I`$, each $`\lambda`$ is shifted by 1.

$`\det(2P - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & 1 - \lambda \end{vmatrix} = (1 - \lambda)^2 - 1 = \lambda^2 - 2\lambda = \lambda(\lambda - 2) = 0`$

$`\therefore \lambda_1 = 2, \lambda_2 = 0`$ for $`2P`$, meaning for $`E = 2P - I`$: $`\lambda_1 = 1, \lambda_2 = -1`$.

**Example 4:** Singular Matrix

```math
A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}
```

Find $`\lambda`$'s and the corresponding $`\mathbf{x}`$'s.

```math
\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 2 \\ 2 & 4 - \lambda \end{vmatrix} = (1 - \lambda)(4 - \lambda) - 4 = \lambda^2 - 5\lambda = \lambda(\lambda - 5) = 0
```

```math
\therefore \lambda_1 = 5 \text{ and } \lambda_2 = 0
```

i) $`\lambda_1 = 5`$: $`(A - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`2x_1 - x_2 = 0 \implies x_2 = 2x_1`$. Pick $`x_1 = 1 \implies x_2 = 2`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$

ii) $`\lambda_2 = 0`$: $`(A - 0I)\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_1 + 2x_2 = 0 \implies x_1 = -2x_2`$. Pick $`x_2 = 1 \implies x_1 = -2`$: $`\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}`$

**Remark 1:** $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$ are in the nullspace of $`(A - \lambda I)`$.

In this example, $`\lambda_2 = 0`$ is an eigenvalue because $`A`$ is singular. If $`A`$ is invertible, then "zero" is NOT an eigenvalue: $`A\mathbf{x} = \mathbf{0} \implies \mathbf{x} = \mathbf{0}`$.

**Remark 2:** For $`A \in \mathbb{R}^{2 \times 2}`$, when $`A - \lambda I`$ is singular, both rows are multiples of a vector $`(a, b)`$:

```math
\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

The eigenvector is any multiple of $`(b, -a)`$:

```math
\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} b \\ -a \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

**Remark 3:** This example has two distinct eigenvalues, which span $`\mathbb{R}^2`$. If a $`2 \times 2`$ matrix has only one eigenvalue, then it cannot span $`\mathbb{R}^2`$.

e.g., $`A = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}`$, $`\det(A - \lambda I) = (3 - \lambda)^2 = 0`$, so $`\lambda = 3`$. Geometric multiplicity = 1, algebraic multiplicity = 2. The matrix $`A`$ is **defective** and is **not diagonalizable**.

### 1.7 Imaginary Eigenvalues

**Example 5:** The 90-degree rotation

```math
R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
```

has no real eigenvectors.

```math
\det(R - \lambda I) = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0 \implies \lambda = \pm i
```

$`\lambda_1 + \lambda_2 = 0 = \text{trace}(R)`$, $`\lambda_1\lambda_2 = 1 = \det(R)`$.

$`R\mathbf{x} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -y \\ x \end{pmatrix}`$

After a rotation, no real vector $`R\mathbf{x}`$ stays in the same direction as $`\mathbf{x}`$, meaning that there is no real $`\lambda`$ such that $`R\mathbf{x} = \lambda\mathbf{x}`$.

i) $`\lambda_1 = i`$: $`(R - iI)\mathbf{x} = \begin{pmatrix} -i & -1 \\ 1 & -i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`-ix_1 = x_2`$, choose $`x_1 = 1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ -i \end{pmatrix}`$

ii) $`\lambda_2 = -i`$: $`(R + iI)\mathbf{x} = \begin{pmatrix} i & -1 \\ 1 & i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`ix_1 = x_2`$, choose $`x_1 = 1`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ i \end{pmatrix}`$

The complex eigenvectors $`\mathbf{x}_1`$ and $`\mathbf{x}_2`$ keep their direction as they are rotated in complex space.

### 1.8 Rotation Matrix Eigenvalues

A rotation matrix has $`\lambda = e^{i\theta}`$ and $`e^{-i\theta}`$.

```math
R\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}, \quad R\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} -\sin\theta \\ \cos\theta \end{pmatrix}
```

```math
\therefore R = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
```

```math
\det(R - \lambda I) = \begin{vmatrix} \cos\theta - \lambda & -\sin\theta \\ \sin\theta & \cos\theta - \lambda \end{vmatrix} = (\cos\theta - \lambda)^2 + \sin^2\theta
```

```math
= \lambda^2 - 2\cos\theta\,\lambda + \cos^2\theta + \sin^2\theta = \lambda^2 - 2\cos\theta\,\lambda + 1
```

```math
\therefore \lambda = \cos\theta \pm \sqrt{\cos^2\theta - 1} = \cos\theta \pm \sqrt{-\sin^2\theta} = \cos\theta \pm i\sin\theta
```

```math
\lambda_1 = \cos\theta + i\sin\theta = e^{i\theta}, \quad \lambda_2 = \cos\theta - i\sin\theta = e^{-i\theta}
```

**Two properties of the rotation matrix $`R`$:**

1. $`R`$ is an **orthogonal** matrix: $`R^T R = I`$, $`|\lambda| = 1`$.
2. $`|\lambda| = 1`$: Since $`R`$ is orthogonal ($`R^T = R^{-1}`$), all eigenvalues satisfy $`|\lambda| = 1`$. The eigenvalues are $`\lambda = e^{\pm i\theta}`$.

### 1.9 Eigenvalues of AB and A+B

Consider $`A\mathbf{x} = \lambda\mathbf{x}`$, $`B\mathbf{x} = \beta\mathbf{x}`$.

That is, $`\lambda`$ and $`\beta`$ are eigenvalues of $`A`$ and $`B`$.

```math
AB\mathbf{x} = A(\beta\mathbf{x}) = \beta A\mathbf{x} = \beta\lambda\mathbf{x}
```

**This is NOT true.** Why? In general, $`\mathbf{x}`$ is NOT the eigenvector of both $`A`$ and $`B`$.

Similarly, $`(A + B)\mathbf{x} \neq (\lambda + \beta)\mathbf{x}`$.

**Remark:** $`A`$ and $`B`$ share the same $`n`$ independent eigenvectors **if and only if** $`AB = BA`$.

---

<br>

## 2. Diagonalizing a Matrix (6.2)

### 2.1 Key Facts of Diagonalization

**(1)** The columns of $`AX = X\Lambda`$ are $`A\mathbf{x}_k = \lambda_k\mathbf{x}_k`$. The eigenvalue matrix $`\Lambda`$ is diagonal.

**(2)** $`n`$ independent eigenvectors in $`X`$ diagonalize $`A`$:

```math
A = X\Lambda X^{-1}
```

```math
AA = X\Lambda X^{-1} X\Lambda X^{-1} = X\Lambda^2 X^{-1}
```

```math
\boxed{A^k = X\Lambda^k X^{-1}}
```

**(3)** Solve $`\mathbf{u}_{k+1} = A\mathbf{u}_k`$ by $`\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0`$.

**(4)** No equal eigenvalues $`\implies`$ eigenvector $`X`$ is invertible $`\implies`$ $`A`$ can be diagonalized.

Repeated eigenvalues $`\implies`$ $`A`$ might have too few independent eigenvectors $`\implies`$ $`X^{-1}`$ fails.

**(5)** Every matrix $`C = B^{-1}AB`$ has the same eigenvalues as $`A`$. These $`C`$'s are **similar** to $`A`$.

### 2.2 Diagonalization Procedure

When $`\mathbf{x}`$ is an eigenvector, $`A\mathbf{x} = \lambda\mathbf{x}`$. Applying $`A`$ to $`\mathbf{x}`$ is just a multiplication by $`\lambda`$ — **very efficient**.

When $`A`$ is diagonalizable, $`A^{100}\mathbf{x} = X\Lambda^{100}X^{-1}\mathbf{x}`$ — **very efficient** as well.

**Diagonalization:** Suppose $`A_{n \times n}`$ has LI eigenvectors $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$.

Let $`X = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)`$.

```math
AX = (A\mathbf{x}_1 \; A\mathbf{x}_2 \; \cdots \; A\mathbf{x}_n) = (\lambda_1\mathbf{x}_1 \; \lambda_2\mathbf{x}_2 \; \cdots \; \lambda_n\mathbf{x}_n) = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix}
```

```math
AX = X\Lambda
```

```math
AXX^{-1} = X\Lambda X^{-1} \implies A = X\Lambda X^{-1}
```

```math
X^{-1}AX = X^{-1}X\Lambda = \Lambda
```

```math
\therefore \Lambda = X^{-1}AX
```

### 2.3 Worked Example: Diagonalization

```math
A = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}
```

$`\det(A - \lambda I) = \begin{vmatrix} 2 - \lambda & 4 \\ 0 & 6 - \lambda \end{vmatrix} = (6 - \lambda)(2 - \lambda) = 0`$, so $`\lambda_1 = 6, \lambda_2 = 2`$.

i) $`(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`-x_1 + x_2 = 0 \implies \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`(A - \lambda_2 I)\mathbf{x}_2 = \begin{pmatrix} 0 & 4 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_2 = 0`$, $`x_1 = 1`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

iii) $`A\mathbf{x}_1 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = 6\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`A\mathbf{x}_2 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = 2\begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

```math
A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix}
```

```math
\begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 6 & \\ & 2 \end{pmatrix}
```

Multiply: $`X^{-1} = -\begin{pmatrix} 0 & -1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & -1 \end{pmatrix}`$

$`X^{-1}AX = X^{-1}X\Lambda = \Lambda = \begin{pmatrix} 6 & \\ & 2 \end{pmatrix}`$

$`AXX^{-1} = X\Lambda X^{-1}`$

**Watch:** $`A^2 = (X\Lambda X^{-1})(X\Lambda X^{-1}) = X\Lambda^2 X^{-1}`$

$`\implies X^{-1}A^2 X = \Lambda^2`$

$`A^2`$ has the same eigenvectors in $`X`$, and the squared eigenvalues $`36, 4`$ in $`\Lambda^2`$.

### 2.4 Remarks on Diagonalization

**Remark:** The matrix $`X`$ has an inverse because the eigenvectors are **LI**.

**Remark:** Suppose eigenvalues are $`n`$ different numbers, which implies that the $`n`$ eigenvectors are LI. $`X^{-1}`$ exists. Any matrix that has no repeated eigenvalues can be diagonalized.

**Remark:** We can multiply eigenvectors by any nonzero constants: $`\alpha A\mathbf{x} = \alpha\lambda\mathbf{x} = A(\alpha\mathbf{x})`$.

**Remark:** The eigenvalues in $`\Lambda`$ come in the same order as the eigenvectors in $`X`$:

```math
A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix} \implies A(\mathbf{x}_2 \; \mathbf{x}_1) = (\mathbf{x}_2 \; \mathbf{x}_1)\begin{pmatrix} \lambda_2 & \\ & \lambda_1 \end{pmatrix}
```

**Remark:** Some matrices have too few eigenvectors. Those matrices **cannot be diagonalized**.

e.g., $`A = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}`$

$`\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & -1 \\ 1 & -1 - \lambda \end{vmatrix} = -(1 - \lambda)(1 + \lambda) + 1 = -(1 - \lambda^2) + 1 = \lambda^2 = 0`$

$`\therefore \lambda_1 = \lambda_2 = 0`$ (repetition of $`\lambda`$).

$`(A - 0I)\mathbf{x} = A\mathbf{x} = \mathbf{0}`$: $`\begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = x_2`$, $`x_1 = 1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

Only **one** eigenvector. $`A`$ **cannot be diagonalized**.

**Remark:** Invertibility is concerned with **nonzero** eigenvalues. If $`\lambda_i = 0`$, then $`\det(A) = \lambda_1\lambda_2\cdots\lambda_n = 0`$, meaning $`A`$ is **singular**.

### 2.5 Proof: Eigenvectors for Distinct Eigenvalues are LI

**Statement:** Let $`A`$ be an $`n \times n`$ matrix with $`n`$ distinct eigenvalues. Then the corresponding eigenvectors are **LI**.

**Proof:** Suppose the eigenvectors $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$ are **LD** (linearly dependent).

Let $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p`$ are LI, and $`\mathbf{x}_{p+1}, \mathbf{x}_{p+2}, \ldots, \mathbf{x}_n`$ are LD.

That is, there exist constants, **not all zero**, such that:

```math
\mathbf{x}_{p+1} = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_p\mathbf{x}_p
```

Multiply $`A`$ to the linear combination:

```math
A\mathbf{x}_{p+1} = \lambda_{p+1}\mathbf{x}_{p+1} = c_1 A\mathbf{x}_1 + c_2 A\mathbf{x}_2 + \cdots + c_p A\mathbf{x}_p = c_1\lambda_1\mathbf{x}_1 + c_2\lambda_2\mathbf{x}_2 + \cdots + c_p\lambda_p\mathbf{x}_p \quad \text{--- (1)}
```

Multiply $`\lambda_{p+1}`$ to the linear combination:

```math
\lambda_{p+1}\mathbf{x}_{p+1} = c_1\lambda_{p+1}\mathbf{x}_1 + c_2\lambda_{p+1}\mathbf{x}_2 + \cdots + c_p\lambda_{p+1}\mathbf{x}_p \quad \text{--- (2)}
```

Subtract (2) from (1):

```math
\mathbf{0} = c_1(\lambda_1 - \lambda_{p+1})\mathbf{x}_1 + \cdots + c_p(\lambda_p - \lambda_{p+1})\mathbf{x}_p
```

Since $`\lambda`$ are distinct and $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p`$ are LI:

```math
c_1 = c_2 = \cdots = c_p = 0
```

This implies that $`\mathbf{x}_{p+1} = \mathbf{0}`$. This contradicts our assumption.

Therefore $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$ are **LI**. $`\square`$

### 2.6 Powers of A (Markov Matrix Example)

```math
A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix} \quad \text{(Markov matrix)}
```

$`\lambda_1 = 1, \lambda_2 = 0.5 \implies \Lambda = \begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}`$

$`\mathbf{x}_1 = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}, \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \implies X = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}`$

```math
A = X\Lambda X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}
```

```math
A^2 = X\Lambda^2 X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.25 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}
```

```math
A^k = X\Lambda^k X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1^k & \\ & (0.5)^k \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}
```

As $`k \to \infty`$, $`(0.5)^k \to 0`$:

```math
A^{\infty} = X\Lambda^{\infty}X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix} = \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix}
```

**Q.** When does $`A^k \to`$ zero matrix?

**A.** All $`|\lambda_i| < 1`$.

### 2.7 Similar Matrices

Suppose $`\Lambda`$ is fixed. As we change the eigenvector matrix $`X`$, we get different matrices with the **same eigenvalues** $`\Lambda`$.

```math
A_1 = X_1 \Lambda X_1^{-1}, \quad A_2 = X_2 \Lambda X_2^{-1}, \quad \ldots
```

These are **similar matrices**: $`\Lambda, A_1, A_2, \ldots`$

All $`A_1, A_2, \ldots`$ are similar to $`C`$. They all share the eigenvalues of $`C`$.

Extend this idea to non-diagonalizable matrices. Choose a constant matrix $`C`$ and an invertible matrix $`B`$. We construct $`A = BCB^{-1}`$. $`A`$ and $`C`$ are **similar**.

$`A`$ and $`C`$ have the same $`n`$ eigenvalues.

**Statement:** If $`C\mathbf{x} = \lambda\mathbf{x}`$, then $`BCB^{-1}`$ has the same eigenvalue $`\lambda`$ with new eigenvector $`B\mathbf{x}`$.

**Proof:**

```math
(BCB^{-1})(B\mathbf{x}) = BCI\mathbf{x} = BC\mathbf{x} = B\lambda\mathbf{x} = \lambda(B\mathbf{x}) \quad \square
```

e.g., $`A_1 = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}`$ and $`A_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}`$ are similar.

$`\det(A_1 - \lambda I) = \lambda(1 - \lambda) = 0`$, so $`\lambda = 0, 1`$.

$`\det(A_2 - \lambda I) = \lambda^2 - \lambda + \frac{1}{4} - \frac{1}{4} = \lambda(\lambda - 1) = 0`$, so $`\lambda = 0, 1`$.

### 2.8 Matrix Powers and Fibonacci Numbers

Fibonacci number is the sum of the two previous numbers:

```math
0, 1, 1, 2, 3, 5, 8, 13, \ldots
```

```math
a_k + a_{k+1} = a_{k+2}, \quad k \geq 0
```

Introduce $`a_{k+1} = a_{k+1}`$:

```math
\begin{cases} a_k + a_{k+1} = a_{k+2} \\ a_{k+1} = a_{k+1} \end{cases}
```

Let $`\mathbf{u}_k = \begin{pmatrix} a_{k+1} \\ a_k \end{pmatrix}`$:

```math
\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\mathbf{u}_k = \mathbf{u}_{k+1}
```

```math
\therefore \boxed{\mathbf{u}_{k+1} = A\mathbf{u}_k}
```

with $`\mathbf{u}_0 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$:

$`\mathbf{u}_1 = A\mathbf{u}_0 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`\mathbf{u}_2 = A^2\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{u}_3 = A^3\mathbf{u}_0 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}`$, ...

$`\mathbf{u}_k = A^k\mathbf{u}_0`$

$`A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & -\lambda \end{vmatrix} = -\lambda(1 - \lambda) - 1 = \lambda^2 - \lambda - 1 = 0`$

```math
\therefore \lambda = \frac{1 \pm \sqrt{1 + 4}}{2} = \frac{1}{2} \pm \frac{\sqrt{5}}{2}
```

Two distinct eigenvalues $`\implies`$ two LI eigenvectors $`\implies`$ $`X^{-1}`$ exists $`\implies`$ $`A = X\Lambda X^{-1}`$.

```math
\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0
```

i) Write $`\mathbf{u}_0`$ as a linear combination $`X\mathbf{c}`$:

```math
\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = X\mathbf{c} \implies \mathbf{c} = X^{-1}\mathbf{u}_0
```

ii) Multiply $`\Lambda^k`$ to $`\mathbf{c}`$:

```math
\begin{pmatrix} \lambda_1^k & \\ & \lambda_2^k \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix}
```

iii) Multiply $`X`$ to $`\Lambda^k\mathbf{c}`$:

```math
\mathbf{u}_k = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix} = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2
```

**Generalize to $`A \in \mathbb{R}^{n \times n}`$:**

```math
\mathbf{u}_k = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2 + \cdots + c_n\lambda_n^k\mathbf{x}_n
```

Solution for $`\mathbf{u}_{k+1} = A\mathbf{u}_k`$.

### 2.9 Nondiagonalizable Matrices and Multiplicity

Suppose $`\lambda`$ is an eigenvalue of $`A`$.

1. **Eigenvectors (geometric):** $`A\mathbf{x} = \lambda\mathbf{x}`$, nonzero $`\mathbf{x}`$.
2. **Eigenvalues (algebraic):** $`\det(A - \lambda I) = 0`$.

$`\lambda`$ may be a simple eigenvalue, or a **multiple** eigenvalue (e.g., $`\lambda^2 = 0 \to \lambda = 0`$).

**Multiplicity:**

1. **Geometric multiplicity (GM):** count the independent vectors for $`\lambda`$ = $`\dim \mathcal{N}(A - \lambda I)`$.
2. **Algebraic multiplicity (AM):** count the repetition of $`\lambda`$ among the eigenvalues, i.e., the $`n`$ roots of $`\det(A - \lambda I) = 0`$.

e.g., $`A`$ has $`\lambda = 4, 4, 4`$: AM = 3, GM = 1, 2, or 3.

e.g., $`A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}`$, $`|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 = 0`$, so $`\lambda = 0, 0`$.

AM = 2, GM = 1 (1 eigenvector).

$`A\mathbf{x} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, so $`x_2 = 0`$: $`\mathbf{x} = \begin{pmatrix} c \\ 0 \end{pmatrix}`$.

When **GM < AM**, $`A`$ is **NOT diagonalizable**.

---

<br>
