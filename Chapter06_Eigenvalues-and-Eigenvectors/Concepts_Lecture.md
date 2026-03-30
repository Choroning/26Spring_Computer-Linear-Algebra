# Chapter 6 Lecture -- Eigenvalues and Eigenvectors

> **Last Updated:** 2026-03-30

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
- [4. Solving Linear Differential Equations (6.5)](#4-solving-linear-differential-equations-65)
  - [4.1 Key Facts](#41-key-facts)
  - [4.2 Scalar ODE Review](#42-scalar-ode-review)
  - [4.3 Solution of du/dt = Au](#43-solution-of-dudt--au)
  - [4.4 General n x n Solution Procedure](#44-general-n-x-n-solution-procedure)
  - [4.5 Exponential of a Matrix](#45-exponential-of-a-matrix)
  - [4.6 Second Order Equations](#46-second-order-equations)
  - [4.7 Stability of 2 by 2 Matrices](#47-stability-of-2-by-2-matrices)
  - [4.8 Worked Examples](#48-worked-examples)
- [Summary](#summary)

---

<br>

## 1. Introduction to Eigenvalues (6.1)

### 1.1 Definition of Eigenvalues and Eigenvectors

When $A$ acts on $\mathbf{x}$, it only stretches or compresses the vector $\mathbf{x}$ by $\lambda$, without changing its direction.

$$A\mathbf{x} = \lambda \mathbf{x}$$

- $\mathbf{x}$ is a **non-zero vector**, known as an **eigenvector**.
- $\lambda$ is the **eigenvalue** corresponding to $\mathbf{x}$.

**Formal statement:**

1. If $A\mathbf{x} = \lambda\mathbf{x}$, then $\mathbf{x} \neq \mathbf{0}$ is an eigenvector of $A$, and the number $\lambda$ is the eigenvalue.

2. $A^n \mathbf{x} = \lambda^n \mathbf{x}$ for every $n$; $(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}$; and $A^{-1}\mathbf{x} = \frac{1}{\lambda}\mathbf{x} = \lambda^{-1}\mathbf{x}$ if $\lambda \neq 0$.

**Proof** (for $A^{-1}$):

$$A\mathbf{x} = \lambda\mathbf{x} \implies \mathbf{x} = A^{-1}(A\mathbf{x}) = A^{-1}(\lambda\mathbf{x}) = \lambda A^{-1}\mathbf{x} \quad \square$$

3. $(A - \lambda I)\mathbf{x} = \mathbf{0} \implies \det(A - \lambda I) = 0$. This equation produces $n$ $\lambda$'s.

4. Let $A \in \mathbb{R}^{n \times n}$:

$$\det(A) = \lambda_1 \lambda_2 \cdots \lambda_n; \quad \text{trace}(A) = \lambda_1 + \lambda_2 + \cdots + \lambda_n$$

5. **Projection matrix** $P$ has $\lambda = 1$ or $\lambda = 0$. A square matrix $P$ is called a "projection matrix" if $P^2 = P$.

### 1.2 Powers of A and Eigenvalues

What happens if we multiply $A$ to the relation?

$$A(A\mathbf{x}) = A(\lambda\mathbf{x}) \implies A^2\mathbf{x} = \lambda A\mathbf{x} = \lambda^2 \mathbf{x}$$

If we keep multiplying $\mathbf{x}$ by $A$:

$$A^3\mathbf{x} = \lambda^3\mathbf{x}, \quad A^k\mathbf{x} = \lambda^k\mathbf{x}, \quad \ldots, \quad A^{100}\mathbf{x} = \lambda^{100}\mathbf{x}$$

$$\boxed{A^k \mathbf{x} = \lambda^k \mathbf{x}}$$

**Behavior of $A^k$:**

- If $|\lambda_i| < 1$ for $i = 1, 2, \ldots, n$, then $A^k$ will eventually approach **zero**.
- If any $|\lambda_i| > 1$, then $A^k$ will eventually **grow**.
- If $\lambda = 1$, then the system state does not grow or decay over time, but rather stays **constant**.

When $A\mathbf{x} = \mathbf{x}$, then $\mathbf{x}$ is a **fixed point** or **equilibrium point** of a system, i.e., the system remains steady in $\mathbf{x}$ direction. The system reaches a **steady state**: $A^k \mathbf{x} = \mathbf{x}$.

### 1.3 Properties of Eigenvalues

An $A \in \mathbb{R}^{n \times n}$ matrix has $n$ eigenvalues. We can find the eigenvalues by solving the **characteristic polynomial**:

$$\det(A - \lambda I) = 0$$

**Properties of a matrix significantly influence its eigenvalues:**

**(1)** The **trace** of a matrix $A$ is equal to the **sum** of its eigenvalues.

e.g., $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$, $\text{trace}(A) = 1 + 4 = 5$

$\det(A - \lambda I) = 0 \implies \lambda_1 = 0, \lambda_2 = 5$, and $\lambda_1 + \lambda_2 = 5$.

**(2)** The **determinant** of $A$ is the **product** of its eigenvalues.

e.g., $\det(A) = 1 \cdot 4 - 2 \cdot 2 = 0$, and $\lambda_1 \cdot \lambda_2 = 0 \cdot 5 = 0$.

**(3)** A **symmetric** matrix (i.e., $A = A^T$) has only **real eigenvalues**.

e.g., $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0$, so $\lambda = \pm 1$.

**(4)** A **symmetric positive definite (SPD)** matrix has **real and positive** eigenvalues. SPD matrix is crucial in optimization because it guarantees that a quadratic form $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A\mathbf{x}$ has a **unique minimum**. (See Section 6.3)

**(5)** For a **diagonal** matrix, the eigenvalues are the **diagonal elements**.

e.g., $A = \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} a - \lambda & 0 \\ 0 & b - \lambda \end{vmatrix} = (a - \lambda)(b - \lambda) = 0$, so $\lambda_1 = a, \lambda_2 = b$.

For a **triangular** matrix, the eigenvalues are the **diagonal elements**.

e.g., $A = \begin{pmatrix} a & b \\ 0 & c \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} a - \lambda & b \\ 0 & c - \lambda \end{vmatrix} = (a - \lambda)(c - \lambda) = 0$, so $\lambda_1 = a, \lambda_2 = c$.

**(6)** A **skew-symmetric** matrix (i.e., $A = -A^T$) has **purely imaginary** eigenvalues or **zero** eigenvalues.

e.g., $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ -1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0$, so $\lambda = \pm i$.

### 1.4 The Equation for Eigenvalues

From $A\mathbf{x} = \lambda\mathbf{x}$:

$$A\mathbf{x} - \lambda\mathbf{x} = (A - \lambda I)\mathbf{x} = \mathbf{0}$$

The eigenvectors make up the **nullspace** of $A - \lambda I$.

When we know an eigenvalue $\lambda$, we find an eigenvector by solving $(A - \lambda I)\mathbf{x} = \mathbf{0}$.

If $\mathbf{x}$ is a nonzero vector, then $A - \lambda I$ is **singular**.

$$\boxed{\det(A - \lambda I) = 0}$$

This characteristic polynomial ($\det(A - \lambda I) = 0$) involves only $\lambda$, which appears all along the main diagonal of $A - \lambda I$. The determinant includes $(-\lambda)^n$.

- The equation has $n$ solutions $\lambda_1$ to $\lambda_n$.
- $A_{n \times n}$ has $n$ eigenvalues.

**Proof (Determinant = Product of Eigenvalues):**

$$A\mathbf{x} = \lambda\mathbf{x} \implies \det(A - \lambda I) = 0$$

The characteristic polynomial for $\lambda$:

$$\det(\lambda I - A) = 0 \iff \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0 = 0$$

$$\iff (\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n) = 0$$

Take $\lambda = 0$:

$$\det(0 \cdot I - A) = \det(-A) = (-1)^n \det(A)$$

$$(-\lambda_1)(-\lambda_2)\cdots(-\lambda_n) = (-1)^n \lambda_1\lambda_2\cdots\lambda_n$$

$$\therefore \det(A) = \lambda_1\lambda_2\cdots\lambda_n \quad \square$$

Also, from expanding $(\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n)$:

$$= \lambda^n - (\lambda_1 + \lambda_2 + \cdots + \lambda_n)\lambda^{n-1} + \cdots + \det(A) = 0$$

The coefficient of $\lambda^{n-1}$ equals $\text{trace}(A) = a_{11} + a_{22} + \cdots + a_{nn}$ (skip proof).

**Derivative of a determinant of a matrix:**

i) $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, $\det(A) = ad - bc$.

Let $a = a(x), b = b(x), c = c(x), d = d(x)$.

$$\frac{d}{dx}\det(A) = \frac{d}{dx}(a(x)d(x) - b(x)c(x)) = a'd + ad' - b'c - bc'$$

$$= \underbrace{a'd - b'c} + \underbrace{ad' - bc'}$$

$$\begin{vmatrix} a & b \\ c & d \end{vmatrix}' = \begin{vmatrix} a' & b' \\ c & d \end{vmatrix} + \begin{vmatrix} a & b \\ c' & d' \end{vmatrix}$$

To differentiate a determinant, we differentiate **one row (or column) at a time**, keeping others unchanged.

For $n \times n$ matrix: writing rows as $\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_n$:

$$\begin{vmatrix} a_{11} & \cdots & a_{1n} \\ \vdots & & \vdots \\ a_{n1} & \cdots & a_{nn} \end{vmatrix}' = \sum_{i=1}^{n} \begin{vmatrix} \mathbf{r}_1 \\ \vdots \\ \mathbf{r}_i' \\ \vdots \\ \mathbf{r}_n \end{vmatrix}$$

### 1.5 Determinant and Trace

**Observation 1:** Elimination does **not** preserve eigenvalues.

$$A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}$$

$\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 2 & 6 - \lambda \end{vmatrix} = \lambda^2 - 7\lambda = 0$, so $\lambda_1 = 7, \lambda_2 = 0$.

$\det(A_0 - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 - \lambda = 0$, so $\lambda_1 = 1, \lambda_2 = 0$. (Different!)

**Observation 2:** The product $\lambda_1 \cdot \lambda_2$ and the sum $\lambda_1 + \lambda_2$ can be found from the matrix.

For $A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}$:

- $\lambda_1 \lambda_2 = 7 \cdot 0 = 0 = \det(A) = 6 - 6 = 0$
- $\lambda_1 + \lambda_2 = 7 + 0 = 7 = \text{trace}(A) = 1 + 6 = 7$

The product $\lambda_1 \lambda_2 \cdots \lambda_n$ of $n$ eigenvalues equals the **determinant** of $A$.

The sum $\lambda_1 + \lambda_2 + \cdots + \lambda_n$ equals the sum of the $n$ diagonal entries = **trace** of $A$.

### 1.6 Worked Examples

**Example 1:** Markov Matrix

$$A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix}$$

$$\det(A - \lambda I) = \begin{vmatrix} 0.8 - \lambda & 0.3 \\ 0.2 & 0.7 - \lambda \end{vmatrix} = 0.56 - 1.5\lambda + \lambda^2 - 0.06 = \lambda^2 - \frac{3}{2}\lambda + \frac{1}{2} = (\lambda - \frac{1}{2})(\lambda - 1) = 0$$

$$\therefore \lambda_1 = 1, \quad \lambda_2 = \frac{1}{2}$$

This tells us that $A - \lambda_1 I = A - I$ and $A - \lambda_2 I = A - \frac{1}{2}I$ are **NOT invertible**.

The eigenvectors $\mathbf{x}_1$ and $\mathbf{x}_2$ satisfy $(A - I)\mathbf{x}_1 = \mathbf{0}$ and $(A - \frac{1}{2}I)\mathbf{x}_2 = \mathbf{0}$.

That is, $\mathbf{x}_1 \in \mathcal{N}(A - I)$ and $\mathbf{x}_2 \in \mathcal{N}(A - \frac{1}{2}I)$.

i) $(A - I)\mathbf{x}_1 = \begin{pmatrix} -0.2 & 0.3 \\ 0.2 & -0.3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$-2x_1 + 3x_2 = 0$. Choose $x_1 = 3$, then $x_2 = 2$: $\mathbf{x}_1 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}$

ii) $(A - \frac{1}{2}I)\mathbf{x}_2 = \begin{pmatrix} 0.3 & 0.3 \\ 0.2 & 0.2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_1 + x_2 = 0$, so $x_1 = -x_2$. Choose $x_1 = 1$, then $x_2 = -1$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

Multiply $A$ to $\mathbf{x}_1$: $A\mathbf{x}_1 = \mathbf{x}_1$, $A^2\mathbf{x}_1 = \mathbf{x}_1$, ..., $A^{100}\mathbf{x}_1 = \mathbf{x}_1$.

Multiply $A$ to $\mathbf{x}_2$: $A\mathbf{x}_2 = \frac{1}{2}\mathbf{x}_2$, $A^2\mathbf{x}_2 = (\frac{1}{2})^2\mathbf{x}_2$, ..., $A^{100}\mathbf{x}_2 = (\frac{1}{2})^{100}\mathbf{x}_2$.

Both $\mathbf{x}_1$ and $\mathbf{x}_2$ stay in their own directions.

An eigenvector $\mathbf{x}$ of $A$ is also an eigenvector of every $A^n$ because of $A^n\mathbf{x} = \lambda^n\mathbf{x}$.

Eigenvectors $\mathbf{x}_1$ and $\mathbf{x}_2$ span $\mathbb{R}^2$.

Any vector $\mathbf{x}$ is a linear combination of $\mathbf{x}_1$ and $\mathbf{x}_2$: $\mathbf{x} = c\mathbf{x}_1 + d\mathbf{x}_2 = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix}$.

$$A\mathbf{x} = cA\mathbf{x}_1 + dA\mathbf{x}_2 = c\mathbf{x}_1 + d\tfrac{1}{2}\mathbf{x}_2$$

$$A^2\mathbf{x} = c\mathbf{x}_1 + d(\tfrac{1}{2})^2\mathbf{x}_2$$

$$A^{100}\mathbf{x} = c\,\underbrace{\mathbf{x}_1}_{\text{steady state}} + d(\tfrac{1}{2})^{100}\underbrace{\mathbf{x}_2}_{\text{decaying mode}} \approx c\,\mathbf{x}_1$$

Let $\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$:

$$\begin{pmatrix} 1 \\ 0 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix}$$

$$\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 & 1 \\ 2 & -3 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 0.2 \\ 0.4 \end{pmatrix}$$

$$= 0.2\begin{pmatrix} 3 \\ 2 \end{pmatrix} + 0.4\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix} + \begin{pmatrix} 0.4 \\ -0.4 \end{pmatrix}$$

$$A^{100}\begin{pmatrix} 1 \\ 0 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}$$

Similarly for $\mathbf{x} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$: $A^{100}\begin{pmatrix} 0 \\ 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}$.

$$A^{100}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_1)$$

The higher the power of $A$, the more closely its columns approach the **steady state**.

**Example 2:** Projection Matrix

$$P = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$$

has eigenvalues $\lambda = 1$ and $\lambda = 0$.

$$\det(P - \lambda I) = \begin{vmatrix} \frac{1}{2} - \lambda & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} - \lambda \end{vmatrix} = (\frac{1}{2})^2 \begin{vmatrix} 1 - 2\lambda & 1 \\ 1 & 1 - 2\lambda \end{vmatrix} = \frac{1}{4}[(1 - 2\lambda)^2 - 1] = \frac{1}{4}(4\lambda^2 - 4\lambda) = \lambda(\lambda - 1) = 0$$

$$\therefore \lambda_1 = 1, \; \lambda_2 = 0$$

i) $\lambda_1 = 1$: $(P - I)\mathbf{x}_1 = \begin{pmatrix} -\frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & -\frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = x_2$, choose $x_1 = 1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $\lambda_2 = 0$: $P\mathbf{x}_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 + x_2 = 0$, $x_1 = -x_2$, choose $x_1 = 1$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) $P\mathbf{x}_1 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \mathbf{x}_1$

$P\mathbf{x}_2 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}$

$\therefore \mathbf{x}_1 \in C(P)$: the column space projects onto itself.

$\mathbf{x}_2 \in \mathcal{N}(P)$

iv) Let $\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$

$P\mathbf{w} = P\begin{pmatrix} 1 \\ -1 \end{pmatrix} + P\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \mathbf{0} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}$

The only eigenvalues of a projection matrix are **0 and 1**.

**Example 3:** Exchange Matrix

$$E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$$

has eigenvalues 1 and $-1$.

$$\det(E - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = (\lambda - 1)(\lambda + 1) = 0$$

$$\therefore \lambda_1 = 1, \; \lambda_2 = -1$$

i) $\lambda_1 = 1$: $(E - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $-x_1 + x_2 = 0$, $x_1 = x_2$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $\lambda_2 = -1$: $(E + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = -x_2$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) The eigenvectors $\mathbf{x}_1$ and $\mathbf{x}_2$ are the same as for $P$.

$$E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} - \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 2P - I$$

When a matrix is shifted by $I$, each $\lambda$ is shifted by 1.

$\det(2P - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & 1 - \lambda \end{vmatrix} = (1 - \lambda)^2 - 1 = \lambda^2 - 2\lambda = \lambda(\lambda - 2) = 0$

$\therefore \lambda_1 = 2, \lambda_2 = 0$ for $2P$, meaning for $E = 2P - I$: $\lambda_1 = 1, \lambda_2 = -1$.

**Example 4:** Singular Matrix

$$A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$$

Find $\lambda$'s and the corresponding $\mathbf{x}$'s.

$$\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 2 \\ 2 & 4 - \lambda \end{vmatrix} = (1 - \lambda)(4 - \lambda) - 4 = \lambda^2 - 5\lambda = \lambda(\lambda - 5) = 0$$

$$\therefore \lambda_1 = 5 \text{ and } \lambda_2 = 0$$

i) $\lambda_1 = 5$: $(A - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$2x_1 - x_2 = 0 \implies x_2 = 2x_1$. Pick $x_1 = 1 \implies x_2 = 2$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$

ii) $\lambda_2 = 0$: $(A - 0I)\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_1 + 2x_2 = 0 \implies x_1 = -2x_2$. Pick $x_2 = 1 \implies x_1 = -2$: $\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}$

**Remark 1:** $\mathbf{x}_1$ and $\mathbf{x}_2$ are in the nullspace of $(A - \lambda I)$.

In this example, $\lambda_2 = 0$ is an eigenvalue because $A$ is singular. If $A$ is invertible, then "zero" is NOT an eigenvalue: $A\mathbf{x} = \mathbf{0} \implies \mathbf{x} = \mathbf{0}$.

**Remark 2:** For $A \in \mathbb{R}^{2 \times 2}$, when $A - \lambda I$ is singular, both rows are multiples of a vector $(a, b)$:

$$\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

The eigenvector is any multiple of $(b, -a)$:

$$\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} b \\ -a \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

**Remark 3:** This example has two distinct eigenvalues, which span $\mathbb{R}^2$. If a $2 \times 2$ matrix has only one eigenvalue, then it cannot span $\mathbb{R}^2$.

e.g., $A = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$, $\det(A - \lambda I) = (3 - \lambda)^2 = 0$, so $\lambda = 3$. Geometric multiplicity = 1, algebraic multiplicity = 2. The matrix $A$ is **defective** and is **not diagonalizable**.

### 1.7 Imaginary Eigenvalues

**Example 5:** The 90-degree rotation

$$R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$$

has no real eigenvectors.

$$\det(R - \lambda I) = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0 \implies \lambda = \pm i$$

$\lambda_1 + \lambda_2 = 0 = \text{trace}(R)$, $\lambda_1\lambda_2 = 1 = \det(R)$.

$R\mathbf{x} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -y \\ x \end{pmatrix}$

After a rotation, no real vector $R\mathbf{x}$ stays in the same direction as $\mathbf{x}$, meaning that there is no real $\lambda$ such that $R\mathbf{x} = \lambda\mathbf{x}$.

i) $\lambda_1 = i$: $(R - iI)\mathbf{x} = \begin{pmatrix} -i & -1 \\ 1 & -i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$-ix_1 = x_2$, choose $x_1 = 1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ -i \end{pmatrix}$

ii) $\lambda_2 = -i$: $(R + iI)\mathbf{x} = \begin{pmatrix} i & -1 \\ 1 & i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$ix_1 = x_2$, choose $x_1 = 1$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ i \end{pmatrix}$

The complex eigenvectors $\mathbf{x}_1$ and $\mathbf{x}_2$ keep their direction as they are rotated in complex space.

### 1.8 Rotation Matrix Eigenvalues

A rotation matrix has $\lambda = e^{i\theta}$ and $e^{-i\theta}$.

$$R\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}, \quad R\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} -\sin\theta \\ \cos\theta \end{pmatrix}$$

$$\therefore R = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

$$\det(R - \lambda I) = \begin{vmatrix} \cos\theta - \lambda & -\sin\theta \\ \sin\theta & \cos\theta - \lambda \end{vmatrix} = (\cos\theta - \lambda)^2 + \sin^2\theta$$

$$= \lambda^2 - 2\cos\theta\,\lambda + \cos^2\theta + \sin^2\theta = \lambda^2 - 2\cos\theta\,\lambda + 1$$

$$\therefore \lambda = \cos\theta \pm \sqrt{\cos^2\theta - 1} = \cos\theta \pm \sqrt{-\sin^2\theta} = \cos\theta \pm i\sin\theta$$

$$\lambda_1 = \cos\theta + i\sin\theta = e^{i\theta}, \quad \lambda_2 = \cos\theta - i\sin\theta = e^{-i\theta}$$

**Two properties of the rotation matrix $R$:**

1. $R$ is an **orthogonal** matrix: $R^T R = I$, $|\lambda| = 1$.
2. $R$ is a **skew-symmetric** matrix: $R = -R^T$, $\lambda$ is pure imaginary.

### 1.9 Eigenvalues of AB and A+B

Consider $A\mathbf{x} = \lambda\mathbf{x}$, $B\mathbf{x} = \beta\mathbf{x}$.

That is, $\lambda$ and $\beta$ are eigenvalues of $A$ and $B$.

$$AB\mathbf{x} = A(\beta\mathbf{x}) = \beta A\mathbf{x} = \beta\lambda\mathbf{x}$$

**This is NOT true.** Why? In general, $\mathbf{x}$ is NOT the eigenvector of both $A$ and $B$.

Similarly, $(A + B)\mathbf{x} \neq (\lambda + \beta)\mathbf{x}$.

**Remark:** $A$ and $B$ share the same $n$ independent eigenvectors **if and only if** $AB = BA$.

---

<br>

## 2. Diagonalizing a Matrix (6.2)

### 2.1 Key Facts of Diagonalization

**(1)** The columns of $AX = X\Lambda$ are $A\mathbf{x}_k = \lambda_k\mathbf{x}_k$. The eigenvalue matrix $\Lambda$ is diagonal.

**(2)** $n$ independent eigenvectors in $X$ diagonalize $A$:

$$A = X\Lambda X^{-1}$$

$$AA = X\Lambda X^{-1} X\Lambda X^{-1} = X\Lambda^2 X^{-1}$$

$$\boxed{A^k = X\Lambda^k X^{-1}}$$

**(3)** Solve $\mathbf{u}_{k+1} = A\mathbf{u}_k$ by $\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0$.

**(4)** No equal eigenvalues $\implies$ eigenvector $X$ is invertible $\implies$ $A$ can be diagonalized.

Repeated eigenvalues $\implies$ $A$ might have too few independent eigenvectors $\implies$ $X^{-1}$ fails.

**(5)** Every matrix $C = B^{-1}AB$ has the same eigenvalues as $A$. These $C$'s are **similar** to $A$.

### 2.2 Diagonalization Procedure

When $\mathbf{x}$ is an eigenvector, $A\mathbf{x} = \lambda\mathbf{x}$. Applying $A$ to $\mathbf{x}$ is just a multiplication by $\lambda$ -- **very efficient**.

When $A$ is diagonalizable, $A^{100}\mathbf{x} = X\Lambda^{100}X^{-1}\mathbf{x}$ -- **very efficient** as well.

**Diagonalization:** Suppose $A_{n \times n}$ has LI eigenvectors $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$.

Let $X = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)$.

$$AX = (A\mathbf{x}_1 \; A\mathbf{x}_2 \; \cdots \; A\mathbf{x}_n) = (\lambda_1\mathbf{x}_1 \; \lambda_2\mathbf{x}_2 \; \cdots \; \lambda_n\mathbf{x}_n) = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix}$$

$$AX = X\Lambda$$

$$AXX^{-1} = X\Lambda X^{-1} \implies A = X\Lambda X^{-1}$$

$$X^{-1}AX = X^{-1}X\Lambda = \Lambda$$

$$\therefore \Lambda = X^{-1}AX$$

### 2.3 Worked Example: Diagonalization

$$A = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}$$

$\det(A - \lambda I) = \begin{vmatrix} 2 - \lambda & 4 \\ 0 & 6 - \lambda \end{vmatrix} = (6 - \lambda)(2 - \lambda) = 0$, so $\lambda_1 = 6, \lambda_2 = 2$.

i) $(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$-x_1 + x_2 = 0 \implies \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $(A - \lambda_2 I)\mathbf{x}_2 = \begin{pmatrix} 0 & 4 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_2 = 0$, $x_1 = 1$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

iii) $A\mathbf{x}_1 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = 6\begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $A\mathbf{x}_2 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = 2\begin{pmatrix} 1 \\ 0 \end{pmatrix}$

$$A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix}$$

$$\begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 6 & \\ & 2 \end{pmatrix}$$

Multiply: $X^{-1} = -\begin{pmatrix} 0 & -1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & -1 \end{pmatrix}$

$X^{-1}AX = X^{-1}X\Lambda = \Lambda = \begin{pmatrix} 6 & \\ & 2 \end{pmatrix}$

$AXX^{-1} = X\Lambda X^{-1}$

**Watch:** $A^2 = (X\Lambda X^{-1})(X\Lambda X^{-1}) = X\Lambda^2 X^{-1}$

$\implies X^{-1}A^2 X = \Lambda^2$

$A^2$ has the same eigenvectors in $X$, and the squared eigenvalues $36, 4$ in $\Lambda^2$.

### 2.4 Remarks on Diagonalization

**Remark:** The matrix $X$ has an inverse because the eigenvectors are **LI**.

**Remark:** Suppose eigenvalues are $n$ different numbers, which implies that the $n$ eigenvectors are LI. $X^{-1}$ exists. Any matrix that has no repeated eigenvalues can be diagonalized.

**Remark:** We can multiply eigenvectors by any nonzero constants: $\alpha A\mathbf{x} = \alpha\lambda\mathbf{x} = A(\alpha\mathbf{x})$.

**Remark:** The eigenvalues in $\Lambda$ come in the same order as the eigenvectors in $X$:

$$A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix} \implies A(\mathbf{x}_2 \; \mathbf{x}_1) = (\mathbf{x}_2 \; \mathbf{x}_1)\begin{pmatrix} \lambda_2 & \\ & \lambda_1 \end{pmatrix}$$

**Remark:** Some matrices have too few eigenvectors. Those matrices **cannot be diagonalized**.

e.g., $A = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}$

$\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & -1 \\ 1 & -1 - \lambda \end{vmatrix} = -(1 - \lambda)(1 + \lambda) + 1 = -(1 - \lambda^2) + 1 = \lambda^2 = 0$

$\therefore \lambda_1 = \lambda_2 = 0$ (repetition of $\lambda$).

$(A - 0I)\mathbf{x} = A\mathbf{x} = \mathbf{0}$: $\begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = x_2$, $x_1 = 1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

Only **one** eigenvector. $A$ **cannot be diagonalized**.

**Remark:** Invertibility is concerned with **nonzero** eigenvalues. If $\lambda_i = 0$, then $\det(A) = \lambda_1\lambda_2\cdots\lambda_n = 0$, meaning $A$ is **singular**.

### 2.5 Proof: Eigenvectors for Distinct Eigenvalues are LI

**Statement:** Let $A$ be an $n \times n$ matrix with $n$ distinct eigenvalues. Then the corresponding eigenvectors are **LI**.

**Proof:** Suppose the eigenvectors $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$ are **LD** (linearly dependent).

Let $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p$ are LI, and $\mathbf{x}_{p+1}, \mathbf{x}_{p+2}, \ldots, \mathbf{x}_n$ are LD.

That is, there exist constants, **not all zero**, such that:

$$\mathbf{x}_{p+1} = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_p\mathbf{x}_p$$

Multiply $A$ to the linear combination:

$$A\mathbf{x}_{p+1} = \lambda_{p+1}\mathbf{x}_{p+1} = c_1 A\mathbf{x}_1 + c_2 A\mathbf{x}_2 + \cdots + c_p A\mathbf{x}_p = c_1\lambda_1\mathbf{x}_1 + c_2\lambda_2\mathbf{x}_2 + \cdots + c_p\lambda_p\mathbf{x}_p \quad \text{--- (1)}$$

Multiply $\lambda_{p+1}$ to the linear combination:

$$\lambda_{p+1}\mathbf{x}_{p+1} = c_1\lambda_{p+1}\mathbf{x}_1 + c_2\lambda_{p+1}\mathbf{x}_2 + \cdots + c_p\lambda_{p+1}\mathbf{x}_p \quad \text{--- (2)}$$

Subtract (2) from (1):

$$\mathbf{0} = c_1(\lambda_1 - \lambda_{p+1})\mathbf{x}_1 + \cdots + c_p(\lambda_p - \lambda_{p+1})\mathbf{x}_p$$

Since $\lambda$ are distinct and $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p$ are LI:

$$c_1 = c_2 = \cdots = c_p = 0$$

This implies that $\mathbf{x}_{p+1} = \mathbf{0}$. This contradicts our assumption.

Therefore $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$ are **LI**. $\square$

### 2.6 Powers of A (Markov Matrix Example)

$$A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix} \quad \text{(Markov matrix)}$$

$\lambda_1 = 1, \lambda_2 = 0.5 \implies \Lambda = \begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}$

$\mathbf{x}_1 = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}, \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \implies X = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}$

$$A = X\Lambda X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}$$

$$A^2 = X\Lambda^2 X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.25 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}$$

$$A^k = X\Lambda^k X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1^k & \\ & (0.5)^k \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}$$

As $k \to \infty$, $(0.5)^k \to 0$:

$$A^{\infty} = X\Lambda^{\infty}X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix} = \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix}$$

**Q.** When does $A^k \to$ zero matrix?

**A.** All $|\lambda_i| < 1$.

### 2.7 Similar Matrices

Suppose $\Lambda$ is fixed. As we change the eigenvector matrix $X$, we get different matrices with the **same eigenvalues** $\Lambda$.

$$A_1 = X_1 \Lambda X_1^{-1}, \quad A_2 = X_2 \Lambda X_2^{-1}, \quad \ldots$$

These are **similar matrices**: $\Lambda, A_1, A_2, \ldots$

All $A_1, A_2, \ldots$ are similar to $C$. They all share the eigenvalues of $C$.

Extend this idea to non-diagonalizable matrices. Choose a constant matrix $C$ and an invertible matrix $B$. We construct $A = BCB^{-1}$. $A$ and $C$ are **similar**.

$A$ and $C$ have the same $n$ eigenvalues.

**Statement:** If $C\mathbf{x} = \lambda\mathbf{x}$, then $BCB^{-1}$ has the same eigenvalue $\lambda$ with new eigenvector $B\mathbf{x}$.

**Proof:**

$$(BCB^{-1})(B\mathbf{x}) = BCI\mathbf{x} = BC\mathbf{x} = B\lambda\mathbf{x} = \lambda(B\mathbf{x}) \quad \square$$

e.g., $A_1 = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$ and $A_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}$ are similar.

$\det(A_1 - \lambda I) = \lambda(1 - \lambda) = 0$, so $\lambda = 0, 1$.

$\det(A_2 - \lambda I) = \lambda^2 - \lambda + \frac{1}{4} - \frac{1}{4} = \lambda(\lambda - 1) = 0$, so $\lambda = 0, 1$.

### 2.8 Matrix Powers and Fibonacci Numbers

Fibonacci number is the sum of the two previous numbers:

$$0, 1, 1, 2, 3, 5, 8, 13, \ldots$$

$$a_k + a_{k+1} = a_{k+2}, \quad k \geq 0$$

Introduce $a_{k+1} = a_{k+1}$:

$$\begin{cases} a_k + a_{k+1} = a_{k+2} \\ a_{k+1} = a_{k+1} \end{cases}$$

Let $\mathbf{u}_k = \begin{pmatrix} a_{k+1} \\ a_k \end{pmatrix}$:

$$\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\mathbf{u}_k = \mathbf{u}_{k+1}$$

$$\therefore \boxed{\mathbf{u}_{k+1} = A\mathbf{u}_k}$$

with $\mathbf{u}_0 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$:

$\mathbf{u}_1 = A\mathbf{u}_0 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $\mathbf{u}_2 = A^2\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$, $\mathbf{u}_3 = A^3\mathbf{u}_0 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}$, ...

$\mathbf{u}_k = A^k\mathbf{u}_0$

$A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & -\lambda \end{vmatrix} = -\lambda(1 - \lambda) - 1 = \lambda^2 - \lambda - 1 = 0$

$$\therefore \lambda = \frac{1 \pm \sqrt{1 + 4}}{2} = \frac{1}{2} \pm \frac{\sqrt{5}}{2}$$

Two distinct eigenvalues $\implies$ two LI eigenvectors $\implies$ $X^{-1}$ exists $\implies$ $A = X\Lambda X^{-1}$.

$$\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0$$

i) Write $\mathbf{u}_0$ as a linear combination $X\mathbf{c}$:

$$\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = X\mathbf{c} \implies \mathbf{c} = X^{-1}\mathbf{u}_0$$

ii) Multiply $\Lambda^k$ to $\mathbf{c}$:

$$\begin{pmatrix} \lambda_1^k & \\ & \lambda_2^k \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix}$$

iii) Multiply $X$ to $\Lambda^k\mathbf{c}$:

$$\mathbf{u}_k = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix} = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2$$

**Generalize to $A \in \mathbb{R}^{n \times n}$:**

$$\mathbf{u}_k = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2 + \cdots + c_n\lambda_n^k\mathbf{x}_n$$

Solution for $\mathbf{u}_{k+1} = A\mathbf{u}_k$.

### 2.9 Nondiagonalizable Matrices and Multiplicity

Suppose $\lambda$ is an eigenvalue of $A$.

1. **Eigenvectors (geometric):** $A\mathbf{x} = \lambda\mathbf{x}$, nonzero $\mathbf{x}$.
2. **Eigenvalues (algebraic):** $\det(A - \lambda I) = 0$.

$\lambda$ may be a simple eigenvalue, or a **multiple** eigenvalue (e.g., $\lambda^2 = 0 \to \lambda = 0$).

**Multiplicity:**

1. **Geometric multiplicity (GM):** count the independent vectors for $\lambda$ = $\dim \mathcal{N}(A - \lambda I)$.
2. **Algebraic multiplicity (AM):** count the repetition of $\lambda$ among the eigenvalues, i.e., the $n$ roots of $\det(A - \lambda I) = 0$.

e.g., $A$ has $\lambda = 4, 4, 4$: AM = 3, GM = 1, 2, or 3.

e.g., $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 = 0$, so $\lambda = 0, 0$.

AM = 2, GM = 1 (1 eigenvector).

$A\mathbf{x} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, so $x_2 = 0$: $\mathbf{x} = \begin{pmatrix} c \\ 0 \end{pmatrix}$.

When **GM < AM**, $A$ is **NOT diagonalizable**.

---

<br>

## 3. Symmetric Positive Definite Matrices (6.3)

### 3.1 Symmetric Matrices: Key Properties

**(1)** A symmetric matrix $S$ has $n$ **real** eigenvalues $\lambda_i$, and $n$ **orthonormal** eigenvectors $\mathbf{q}_i$.

**(2)** $S$ is diagonalized by an orthogonal eigenvector matrix $Q$:

$$S = Q\Lambda Q^{-1} = Q\Lambda Q^T$$

**(3a)** **Positive definite:** all $\lambda_i > 0$ and all pivots $> 0$ and all upper-left determinants $> 0$.

**(3b)** The energy test is $\mathbf{x}^T S\mathbf{x} > 0$ for all $\mathbf{x} \neq \mathbf{0}$. Then $S = A^T A$ with independent columns in $A$.

**(4)** **Positive semidefinite** allows $\lambda = 0$, pivot = 0, determinant = 0, energy $\mathbf{x}^T S\mathbf{x} = 0$, any $A$.

**Symmetric matrices ($S = S^T$) are special:**

1. All $n$ eigenvalues $\lambda$ are **real** numbers.
2. The $n$ eigenvectors $\mathbf{q}$ can be chosen **orthogonal**.

e.g., $I = I^T = \begin{pmatrix} 1 & & \\ & 1 & \\ & & \ddots \end{pmatrix}$, $\lambda = 1, 1, \ldots, 1$. $I\mathbf{x} = 1 \cdot \mathbf{x}$ -- every nonzero vector $\mathbf{x}$ is an eigenvector.

We can choose them to be **orthogonal**, and we can rescale them to be **unit vectors** $\implies$ **orthonormal**.

### 3.2 Spectral Theorem

Let $Q = (\mathbf{q}_1 \; \mathbf{q}_2 \; \cdots \; \mathbf{q}_n)$ where $\|\mathbf{q}_i\| = 1$, and $\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 1 & i = j \\ 0 & i \neq j \end{cases}$

$$Q^T Q = \begin{pmatrix} - \mathbf{q}_1^T - \\ - \mathbf{q}_2^T - \\ \vdots \\ - \mathbf{q}_n^T - \end{pmatrix}\begin{pmatrix} | & | & & | \\ \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_n \\ | & | & & | \end{pmatrix} = I$$

For a square matrix $Q$: $Q^T Q = I \implies Q^T = Q^{-1}$.

Recall that $AX = X\Lambda$. For a symmetric matrix $S$, we use orthonormal $Q$ instead of $X$:

$$SQ = Q\Lambda$$

$$\boxed{SQQ^T = Q\Lambda Q^{-1} = Q\Lambda Q^T}$$

**Spectral Theorem:** Every real symmetric matrix $S$ has the form

$$S = Q\Lambda Q^T$$

Proof of symmetry: $S^T = (Q\Lambda Q^T)^T = Q\Lambda^T Q^T = Q\Lambda Q^T = S$.

If $S$ has orthonormal eigenvectors, then $S$ is symmetric.

### 3.3 Proof: Symmetric Matrices Have Orthonormal Eigenbasis

**Statement:** If $S = S^T$, then $S$ has an orthonormal eigenbasis.

**Proof:** Consider two eigenvectors $\mathbf{u}, \mathbf{v}$.

i) $\mathbf{v} \cdot (S\mathbf{u}) = \mathbf{v}^T S\mathbf{u} = \mathbf{v}^T S^T\mathbf{u}$ (by symmetry) $= (S\mathbf{v})^T\mathbf{u} = (S\mathbf{v}) \cdot \mathbf{u}$

ii) Let $\alpha, \beta$ be the corresponding eigenvalues: $S\mathbf{u} = \alpha\mathbf{u}$, $S\mathbf{v} = \beta\mathbf{v}$.

From $\mathbf{v} \cdot (S\mathbf{u}) = (S\mathbf{v}) \cdot \mathbf{u}$, we have:

$$\mathbf{v} \cdot (\alpha\mathbf{u}) = (\beta\mathbf{v}) \cdot \mathbf{u}$$

$$(\alpha - \beta)\mathbf{v} \cdot \mathbf{u} = 0$$

When $\alpha \neq \beta$: $\mathbf{v} \cdot \mathbf{u} = 0$.

$\therefore \mathbf{v} \perp \mathbf{u}$ -- **orthogonal**.

Rescale $\mathbf{u}, \mathbf{v}$ to $\|\mathbf{u}\| = \|\mathbf{v}\| = 1$. Then $\mathbf{u}, \mathbf{v}$ are **orthonormal** eigenvectors. $\square$

### 3.4 Positive Definite Matrices: Definitions

**Def.** $n \times n$ symmetric real matrix $S$ is **positive definite** if

$$\mathbf{x}^T S\mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}, \; \mathbf{x} \in \mathbb{R}^n$$

(equivalently, $\forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$)

**Def.** $n \times n$ symmetric real matrix $S$ is **positive semidefinite** (= non-negative definite) if

$$\mathbf{x}^T S\mathbf{x} \geq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n$$

**Def.** $n \times n$ symmetric real matrix $S$ is **negative definite** if

$$\mathbf{x}^T S\mathbf{x} < 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$$

**Def.** $n \times n$ symmetric real matrix $S$ is **negative semidefinite** if

$$\mathbf{x}^T S\mathbf{x} \leq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n$$

### 3.5 Properties of Positive Definite Matrices

**Property 1.** Positive definite matrix $S$ has **all positive eigenvalues**.

**Proof:** $S$ is symmetric $\implies S = Q\Lambda Q^T$, $S\mathbf{x} = Q\Lambda Q^T\mathbf{x}$.

$$\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\Lambda Q^T\mathbf{x} = (Q^T\mathbf{x})^T \Lambda (Q^T\mathbf{x})$$

Let $\mathbf{y} = Q^T\mathbf{x}$:

$$= \mathbf{y}^T \Lambda \mathbf{y} = \mathbf{y}^T\begin{pmatrix} \lambda_1 y_1 \\ \lambda_2 y_2 \\ \vdots \\ \lambda_n y_n \end{pmatrix} = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2$$

Since $S$ is positive definite, $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$.

Choose $\mathbf{x} = Q\mathbf{e}_i$ where $\mathbf{e}_i$ is the $i$-th column of $I$. Then $\mathbf{y} = Q^T\mathbf{x} = Q^T Q\mathbf{e}_i = \mathbf{e}_i$, which leads to $\mathbf{x}^T S\mathbf{x} = \lambda_i > 0$.

Therefore all $\lambda_1, \lambda_2, \ldots, \lambda_n > 0$. $\square$

**Remark:** All positive eigenvalues does NOT imply the matrix is positive definite.

e.g., $A = \begin{pmatrix} 1 & -3 \\ 0 & 1 \end{pmatrix} \to \lambda_1 = \lambda_2 = 1 > 0$, but $\mathbf{x}^T A\mathbf{x} = (x_1, x_2)\begin{pmatrix} x_1 - 3x_2 \\ x_2 \end{pmatrix} = x_1^2 - 3x_1 x_2 + x_2^2 < 0$ when $x_1 = x_2 = 1$.

(The matrix must be **symmetric** for the equivalence to hold.)

**Positive energy is closely connected to positive eigenvalues ($\lambda > 0$):**

If $S\mathbf{x} = \lambda\mathbf{x}$, then $\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T\lambda\mathbf{x} = \lambda\mathbf{x}^T\mathbf{x} = \lambda\|\mathbf{x}\|^2$.

So $\lambda > 0 \implies \mathbf{x}^T S\mathbf{x} > 0$ for $\mathbf{x} \neq \mathbf{0}$.

**Statement:** If $\mathbf{x}^T S\mathbf{x} > 0$ for the eigenvectors of $S$, then $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \neq \mathbf{0}$.

**Proof:** Let $\mathbf{x} = Q\mathbf{c} = c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n$ where $\mathbf{q}_i$ is the $i$-th eigenvector of $S$ and an orthonormal vector.

$$\mathbf{x}^T S\mathbf{x} = (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T S(c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)$$

$$= (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T(c_1\lambda_1\mathbf{q}_1 + c_2\lambda_2\mathbf{q}_2 + \cdots + c_n\lambda_n\mathbf{q}_n)$$

$$= c_1^2\lambda_1\mathbf{q}_1^T\mathbf{q}_1 + c_2^2\lambda_2\mathbf{q}_2^T\mathbf{q}_2 + \cdots + c_n^2\lambda_n\mathbf{q}_n^T\mathbf{q}_n > 0$$

if $\lambda_1, \lambda_2, \ldots, \lambda_n > 0$. $\square$

**Property 2.** If $S$ is positive definite, then it is **invertible**, $\det(S) > 0$, and $S^{-1}$ is positive definite.

**Proof:**

i) $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$

$\implies S\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$

$\implies \mathcal{N}(S) = \{\mathbf{0}\}$

$\implies S$ has full rank $\iff \dim C(S) = n$

$\therefore S$ is **invertible**.

ii) Since $S$ is positive definite, $S$ has all positive eigenvalues.

Therefore $0 < \lambda_1\lambda_2\cdots\lambda_n = \det(S)$.

iii) $(S^{-1})^T = (S^T)^{-1} = S^{-1}$, so $S^{-1}$ is **symmetric**.

$\mathbf{x}^T S^{-1}\mathbf{x} = \mathbf{x}^T S^{-1}SS^{-1}\mathbf{x} = (S^{-T}\mathbf{x})^T S(S^{-1}\mathbf{x}) = (S^{-1}\mathbf{x})^T S(S^{-1}\mathbf{x})$

Let $\mathbf{z} = S^{-1}\mathbf{x}$: $= \mathbf{z}^T S\mathbf{z} > 0$. $\square$

**Property 3.** If $S$ is positive definite, then $S = A^T A$ with independent columns of $A$.

**Proof:** Symmetric matrix $S = Q\Lambda Q^T$. Since $S$ is positive definite, all eigenvalues are positive. The diagonal matrix $\Lambda$ has a square root $\sqrt{\Lambda}$:

$$\Lambda = \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix} = \begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix}\begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix} = \sqrt{\Lambda}\sqrt{\Lambda}$$

Let $A = Q\sqrt{\Lambda}Q^T$. Then:

$$A^T A = (Q\sqrt{\Lambda}Q^T)(Q\sqrt{\Lambda}Q^T) = Q\sqrt{\Lambda}\underbrace{Q^T Q}_{I}\sqrt{\Lambda}Q^T = Q\sqrt{\Lambda}\sqrt{\Lambda}Q^T = Q\Lambda Q^T = S$$

**Remark:** Energy $= \mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2$.

$\|A\mathbf{x}\| > 0$ if $A\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$ $\iff$ $A$ has full rank.

### 3.6 How to Check if a Matrix is Positive Definite

**(1) Positive eigenvalues:** A positive definite matrix $S$ has all positive eigenvalues.

e.g., $S = \mathbf{u}\mathbf{v}^T$, rank 1 matrix is NOT positive definite.

$A = \begin{pmatrix} a \\ b \end{pmatrix}(1 \; k) = \begin{pmatrix} a & ka \\ b & kb \end{pmatrix}$

Rank 1 matrix $\to$ 1 LI vector $\to$ $(n-1)$ dependent vectors.

$A\mathbf{u} = \mathbf{u}\mathbf{v}^T\mathbf{u} = (\mathbf{v}^T\mathbf{u})\mathbf{u}$, so $\lambda = \mathbf{v}^T\mathbf{u} = (1, k)\begin{pmatrix} a \\ b \end{pmatrix} = a + kb$.

Let $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_{n-1}$ be basis vectors of $\mathcal{N}(A)$. $A\mathbf{w}_i = \mathbf{0} = 0 \cdot \mathbf{w}_i$, so $\lambda = 0$ for $i = 1, 2, \ldots, n-1$.

$\therefore (n-1)$ zero eigenvalues.

**(2) Positive energy test:** $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$.

**(3) Positive energy test for $S = A^T A$:**

$$\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2$$

e.g., $S = A^T A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \end{pmatrix}^T\begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \end{pmatrix} = \begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix} \cdot \ldots$

Note: $\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2$, so $\text{rank}(A) = 2 < 3$, $\dim C(A^T) = 2$, $\dim \mathcal{N}(A) = 1$.

Since $\dim \mathcal{N}(A) \neq 0$, $\mathbf{x}^T S\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0$. $S$ is NOT positive definite. $S$ is **positive semidefinite**.

**(4) Determinant test:** Check if $\det(S) > 0$ -- but more precisely, check all **upper-left determinants**.

$$S = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & -1 & 2 & -1 \\ & & -1 & 2 \end{pmatrix}$$

$D_1 = |2| = 2 > 0$

$D_2 = \begin{vmatrix} 2 & -1 \\ -1 & 2 \end{vmatrix} = 3 > 0$

$D_3 = \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{vmatrix} = 8 - 4 = 4 > 0$

$D_4 = |S| = 2D_3 + \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & -1 & -1 \end{vmatrix} = 8 - 4 + 1 = 5 > 0$

All upper-left determinants are positive, so $S$ is positive definite.

**(5) Pivot test:** Check if the **pivots** are positive.

For $S = \begin{pmatrix} 2 & -1 \\ -1 & 2 & -1 \\ & -1 & 2 \end{pmatrix}$:

After elimination: pivots are $2, \frac{3}{2}, \frac{4}{3}$ -- all positive.

$S = LU$, $S = LDL^T$, $S = A^T A$.

Leading determinants are closely related to pivots: $D_2/D_1$, $D_3/D_1$, etc.

**For a $2 \times 2$ SPD matrix** $S = \begin{pmatrix} a & b \\ b & d \end{pmatrix}$:

- Determinant test: $D_1 = a > 0$, $D_2 = ad - b^2 > 0$.
- Pivot test: $d_1 = a > 0$, $d_2 = \frac{ad - b^2}{a} > 0$.
- Eigenvalues: $\lambda_1 > 0$, $\lambda_2 > 0$.
- Energy: $(x \; y)\begin{pmatrix} a & b \\ b & d \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = ax^2 + bxy + byx + dy^2 = ax^2 + 2bxy + dy^2 > 0$.

### 3.7 Worked Examples: Positive Definite and Semidefinite

**Example:** $S = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix} \implies \lambda_1 = 2 > 0, \lambda_2 = 6 > 0$

$S\mathbf{x} = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix}$

$\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix} = 2x_1^2 + 6x_2^2 > 0$. $\therefore S$ is positive definite.

**Example:** $S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}$

$|S - \lambda I| = \begin{vmatrix} 5 - \lambda & 4 \\ 4 & 5 - \lambda \end{vmatrix} = \lambda^2 - 10\lambda + 25 - 16 = (\lambda - 9)(\lambda - 1) = 0$

$\therefore \lambda_1 = 9 > 0$, $\lambda_2 = 1 > 0$.

i) $(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = x_2$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 + x_2 = 0$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) $Q = (\mathbf{x}_1 \; \mathbf{x}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$, $Q^{-1} = Q = Q^T$

iv) $\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T(Q\Lambda Q^T)\mathbf{x}$, let $\mathbf{y} = Q^T\mathbf{x}$:

$= \lambda_1 y_1^2 + \lambda_2 y_2^2 = 9y_1^2 + y_2^2 > 0$

$\therefore S$ is positive definite.

**Alternative approach (energy test):**

$\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = (x_1, x_2)\begin{pmatrix} 5x_1 + 4x_2 \\ 4x_1 + 5x_2 \end{pmatrix}$

$= 5x_1^2 + 4x_1x_2 + 4x_1x_2 + 5x_2^2 = x_1^2 + 4x_1^2 + 8x_1x_2 + 4x_2^2 + x_1^2$

$= x_1^2 + \underbrace{4(x_1^2 + 2x_1x_2 + x_2^2)}_{4(x_1 + x_2)^2} > 0$

**Example:** $S = \begin{pmatrix} 4 & 5 \\ 5 & 4 \end{pmatrix}$

$|S - \lambda I| = \lambda^2 - 8\lambda + 16 - 25 = \lambda^2 - 8\lambda - 9 = (\lambda - 9)(\lambda + 1) = 0$

$\therefore \lambda_1 = 9$ and $\lambda_2 = -1$.

Since $\lambda_2 < 0$, $S$ is **NOT** positive definite.

### 3.8 The Ellipse and Quadratic Forms

**The Ellipse** $ax^2 + 2bxy + cy^2 = 1$.

**Example 1:** Consider an ellipse $5x^2 + 8xy + 5y^2 = 1$.

$$(x \; y)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 1 \implies \mathbf{x}^T S\mathbf{x} = 1$$

$S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}$, $\lambda^2 - 10\lambda + 9 = (\lambda - 9)(\lambda - 1) = 0$, so $\lambda = 9, 1$.

i) $(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) $\mathbf{q}_1 = \frac{1}{\sqrt{2}}\mathbf{x}_1$, $\mathbf{q}_2 = \frac{1}{\sqrt{2}}\mathbf{x}_2$

$Q = (\mathbf{q}_1 \; \mathbf{q}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$, $Q^{-1} = Q^T = Q$

$S = Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T$

iv) $\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T\mathbf{x}$

Let $\mathbf{z} = Q^T\mathbf{x}$: $= \mathbf{z}^T\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}\mathbf{z} = (z_1, z_2)\begin{pmatrix} 9z_1 \\ z_2 \end{pmatrix} = 9z_1^2 + z_2^2$

$\mathbf{z} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix} x_1 + x_2 \\ x_1 - x_2 \end{pmatrix}$

$\therefore z_1 = \frac{x_1 + x_2}{\sqrt{2}}, \; z_2 = \frac{x_1 - x_2}{\sqrt{2}}$

$\mathbf{x}^T S\mathbf{x} = 9z_1^2 + z_2^2 = 1$

When $z_2 = 0$: $z_1^2 = \frac{1}{9}$, so $z_1 = \pm \frac{1}{3}$ (semi-axis along $\mathbf{q}_1$).

When $z_1 = 0$: $z_2^2 = 1$, so $z_2 = \pm 1$ (semi-axis along $\mathbf{q}_2$).

**Example 2:** Positive semidefinite

$$T = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}\begin{pmatrix} 3 & 1 \end{pmatrix} = A^T A$$

Determinants: $D_1 = 9 > 0$, $D_2 = 9 - 9 = 0 \to$ No inverse.

$\lambda^2 - \text{trace}(T)\lambda + \det(T) = \lambda^2 - 10\lambda = 0$, so $\lambda_1 = 10, \lambda_2 = 0$.

$(T - 10I)\mathbf{x}_1 = \begin{pmatrix} -1 & 3 \\ 3 & -9 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$

$(T - 0I)\mathbf{x}_2 = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -3 \end{pmatrix}$

Energy: $ax^2 + 2bxy + dy^2 = 9x^2 + 6xy + y^2 = (3x + y)^2 \geq 0$.

$E(x,y) = 1$ is a band: $3x + y = \pm 1$. The graph of energy $E(x,y)$ is a valley. Axes of the band along eigenvectors of $T$.

### 3.9 Positive Definite Matrices and Minimum Problems

Consider energy $E = \mathbf{x}^T S\mathbf{x} = \begin{pmatrix} x \\ y \end{pmatrix}^T\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 5x^2 + 8xy + 5y^2 > 0$

(Bowl shape, **convex**)

$E(x,y) = 0$ when $x = y = 0$. This connects to **minimum problems**.

The matrix of second derivatives is positive definite at all points.

1st derivatives: $\frac{\partial E}{\partial x} = 10x + 8y$, $\frac{\partial E}{\partial y} = 8x + 10y$.

2nd derivatives: $\frac{\partial^2 E}{\partial x^2} = 10 > 0$, $\frac{\partial^2 E}{\partial x \partial y} = 8 > 0$, $\frac{\partial^2 E}{\partial y^2} = 10 > 0$.

$$\nabla E = \left(\frac{\partial E}{\partial x}, \frac{\partial E}{\partial y}\right) \quad \xrightarrow{(x,y) = (0,0)} \quad \nabla E = (0, 0)$$

$$J(\nabla E^T) = \begin{pmatrix} \frac{\partial^2 E}{\partial x^2} & \frac{\partial^2 E}{\partial y \partial x} \\ \frac{\partial^2 E}{\partial x \partial y} & \frac{\partial^2 E}{\partial y^2} \end{pmatrix} = H \quad \text{(Hessian matrix)}$$

$H_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}$

If we define the energy by $\frac{1}{2}\mathbf{x}^T S\mathbf{x}$, then $H = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix} = S$.

e.g., $E = x^2 + y^2$ has minimum at $x = y = 0$.

What about $f = e^{x^2 + y^2}$?

$\nabla f = \left(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}\right)$ where $\frac{\partial f}{\partial x} = e^{x^2+y^2} \cdot 2x$, $\frac{\partial f}{\partial y} = e^{x^2+y^2} \cdot 2y$.

$J(\nabla f)^T = \begin{pmatrix} (2 + 4x^2)e^{x^2+y^2} & 4xy \, e^{x^2+y^2} \\ 4xy \, e^{x^2+y^2} & (2 + 4y^2)e^{x^2+y^2} \end{pmatrix} = S$

$S$ is positive definite, so $f$ is **strictly convex**.

**Example 1 (Positive Definite):**

$S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix}$

Determinants: $D_1 = 9 > 0$, $D_2 = 27 - 9 = 18 > 0$.

Pivots: $\begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} \xrightarrow{R_2 - \frac{1}{3}R_1} \begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix}$. Pivots: 9, 2 (both positive).

Energy: $E(x,y) = 9x^2 + 6xy + 3y^2 = (3x + y)^2 + 2y^2 > 0$. (Strictly convex function.)

Trace and determinants: $\text{trace}(S) = 12$, $\det(S) = 18$.

$\lambda^2 - 12\lambda + 18 = 0 \implies \lambda = 6 \pm \sqrt{36 - 18} = 6 \pm 3\sqrt{2}$

i) $\lambda_1 = 6 + 3\sqrt{2}$: $(S - \lambda_1 I)\mathbf{x} = \begin{pmatrix} 3 - 3\sqrt{2} & 3 \\ 3 & 3 - 3\sqrt{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_2 = (-1 + \sqrt{2})x_1$, choose $x_1 = 1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ -1 + \sqrt{2} \end{pmatrix}$

ii) $\lambda_2 = 6 - 3\sqrt{2}$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 - \sqrt{2} \end{pmatrix}$

**Decomposition:**

$S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} = LU = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix} = LDL^T = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & \frac{1}{3} \\ 0 & 1 \end{pmatrix}$

$= (L\sqrt{D})(\sqrt{D}L^T) = A^T A$ where $A = \begin{pmatrix} 3 & 1 \\ 0 & \sqrt{2} \end{pmatrix}$

### 3.10 Positive Semidefinite Matrices

$\mathbf{x}^T S\mathbf{x} \geq 0$

$$S = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$$

$\det(S) = 4 - 4 = 0$ -- singular matrix.

$\text{trace}(S) = 1 + 4 = 5$. $\lambda^2 - 5\lambda + 0 = \lambda(\lambda - 5) = 0$, so $\lambda_1 = 5, \lambda_2 = 0$.

i) $(S - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_2 = 2x_1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$

ii) $S\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 + 2x_2 = 0$: $\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}$

$E_{21}S = U$: $S = E_{21}^{-1}U = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = LDL^T$

$= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$

$= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}$

$S$ is decomposed into $A^T A$ with **dependent** columns in $A$.

$\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0$

### 3.11 Congruent Matrices

If $S$ is positive semidefinite, so is every matrix $A^T SA$:

$$\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} \geq 0 \quad \forall \mathbf{x}$$

where $\mathbf{y} = A\mathbf{x}$. $\square$

Suppose $\mathbf{x}^T S\mathbf{x} > 0$ and $A\mathbf{x} \neq \mathbf{0}$ for all $\mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$. Then $A^T SA$ is **positive definite**.

$$\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} > 0 \quad \forall \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\} \quad \square$$

The matrix $A^T SA$ is called **"congruent"** to $S$.

If $S^T = S$ has $P$ positive eigenvalues, $N$ negative eigenvalues, and $Z$ zero eigenvalues, then the same is true for $A^T SA$ provided $A$ is invertible.

**Proof of skew-symmetric eigenvalue property:**

$A = -\bar{A}^T$. $A\mathbf{x} = \lambda\mathbf{x}$, $\overline{A\mathbf{x}} = \bar{\lambda}\bar{\mathbf{x}}$.

$\bar{\mathbf{x}}^T A\mathbf{x} = \bar{\mathbf{x}}^T(-\bar{A}^T)\mathbf{x} = -(\bar{A}\bar{\mathbf{x}})^T\mathbf{x}$

$= \bar{\mathbf{x}}^T\lambda\mathbf{x}$ and $= -(\bar{\lambda}\bar{\mathbf{x}})^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}$

$\lambda\bar{\mathbf{x}}^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}$

$(\lambda + \bar{\lambda})\bar{\mathbf{x}}^T\mathbf{x} = 0$

$\therefore \lambda + \bar{\lambda} = (r + is) + (r - is) = 2r = 0$

Real part of $\lambda$ is **zero**. Antisymmetric matrix has zero or **pure imaginary** eigenvalues.

### 3.12 Optimization and Machine Learning

**Gradient descent** to minimize $f(x)$.

i) $f(x) = x^2 + 4x + 4 \implies f'(x) = 2x + 4$

ii) $x_0 = 10$, $\eta = 0.1$ (learning rate), stop criterion $|f'(x)| < 0.01$

iii) $x_{k+1} = x_k - \eta f'(x_k)$ -- **steepest direction**.

Iterate the approximation until $x_k \to x^*$ (the minimum point).

For $f(\mathbf{x})$: $\mathbf{x}_{k+1} = \mathbf{x}_k - \eta \nabla f(\mathbf{x}_k)$

If $\mathbf{x}_k \to \mathbf{x}^*$, then $\nabla f$ at $\mathbf{x}^*$ is zero $\implies \mathbf{x}_{k+1} \approx \mathbf{x}_k$.

If $J(\nabla f)^T = S$ is positive definite, then $f$ is **convex** -- easy to find $\mathbf{x}^*$.

---

<br>

## 4. Solving Linear Differential Equations (6.5)

### 4.1 Key Facts

**(1)** If $A\mathbf{x} = \lambda\mathbf{x}$, then $\mathbf{u}(t) = e^{\lambda t}\mathbf{x}$ will solve $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$. Each $\lambda$ and $\mathbf{x}$ give a solution $e^{\lambda t}\mathbf{x}$.

**(2)** If $A = X\Lambda X^{-1}$, then

$$\mathbf{u}(t) = e^{At}\mathbf{u}(0) = Xe^{\Lambda t}X^{-1}\mathbf{u}(0) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n$$

**(3)** Matrix exponential: $e^{At} = I + At + \cdots + \frac{(At)^n}{n!} + \cdots = Xe^{\Lambda t}X^{-1}$ if $A = X\Lambda X^{-1}$.

**(4)** $A$ is **stable** and $\mathbf{u}(t) \to \mathbf{0}$ and $e^{At} \to 0$ when all eigenvalues of $A$ have **real part $< 0$**.

**(5)** $\frac{d^2\mathbf{u}}{dt^2} + B\frac{d\mathbf{u}}{dt} + C\mathbf{u} = 0$ means $\frac{d}{dt}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix} = \begin{pmatrix} 0 & I \\ -C & -B \end{pmatrix}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix}$ where $\mathbf{v} = \frac{d\mathbf{u}}{dt}$.

### 4.2 Scalar ODE Review

Consider ordinary differential equations:

$$\frac{du}{dt} = u \implies u(t) = Ce^t$$

$$\boxed{\frac{du}{dt} = \lambda u} \implies \boxed{u(t) = Ce^{\lambda t}}$$

Check: $\frac{du}{dt} = \lambda Ce^{\lambda t} = \lambda u$.

At $t = 0$: $u(0) = Ce^{\lambda \cdot 0} = C \cdot 1 = C$. The initial value is $C$.

$\therefore u(t) = u(0)e^{\lambda t} = u_0 e^{\lambda t}$

**Behavior:**

- $\lambda > 0$: **unstable** (growing)
- $\lambda = 0$: **steady state** (constant)
- $\lambda < 0$: **stable** (decaying)

**What if $\lambda \in \mathbb{C}$?** e.g., $\lambda = -1 + 2i$

$e^{\lambda t} = e^{-t + 2it} = e^{-t}e^{2it}$

$|e^{\lambda t}| = |e^{-t}||e^{2it}| = e^{-t} \cdot \underbrace{1}_{|e^{i\theta}| = 1}$

**Observation:** Stability depends on the **real part** of $\lambda$.

- $\text{Re}(\lambda) > 0$: unstable
- $\text{Re}(\lambda) = 0$: steady (oscillation)
- $\text{Re}(\lambda) < 0$: stable

### 4.3 Solution of du/dt = Au

**Example 1:** Consider a coupled ODE:

$$\frac{dy}{dt} = z, \quad \frac{dz}{dt} = y$$

$$\implies \frac{d}{dt}\begin{pmatrix} y \\ z \end{pmatrix} = \begin{pmatrix} z \\ y \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y \\ z \end{pmatrix}$$

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u}$$

Given an initial condition $\mathbf{u}_0 = \begin{pmatrix} y(0) \\ z(0) \end{pmatrix}$, what is $\mathbf{u}(t) = ?$

Let $\lambda$ be an eigenvalue of $A$ and $\mathbf{x}$ be the corresponding eigenvector.

Choose $\mathbf{u} = e^{\lambda t}\mathbf{x}$ and plug into the coupled ODE:

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u} = e^{\lambda t}A\mathbf{x} = e^{\lambda t}\lambda\mathbf{x} = \lambda e^{\lambda t}\mathbf{x}$$

Since $\mathbf{u} = e^{\lambda t}\mathbf{x}$ satisfies the coupled ODE, $e^{\lambda t}\mathbf{x}$ is a solution to $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$.

$$\implies \mathbf{u} = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2$$

What are the eigenvalues of $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$?

$|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0$, so $\lambda = \pm 1$.

i) $\lambda_1 = 1$: $(A - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, so $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $\lambda_2 = -1$: $(A + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, so $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) The complete solution: $\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2$

iv) The constants $c_1$ and $c_2$ can be determined by the initial condition $\mathbf{u}_0 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$:

$$\mathbf{u}(t=0) = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$$

Write $\mathbf{u}_0$ as a combination of the eigenvectors of $A$:

$$= (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$$

$$\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = -\frac{1}{2}\begin{pmatrix} -1 & -1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$$

### 4.4 General n x n Solution Procedure

Generalize the idea to $n \times n$ matrix $A$:

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u}$$

$A$ has $n$ eigenvalues and $n$ eigenvectors ($\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$).

1. Write $\mathbf{u}_0$ as $\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_n\mathbf{x}_n$.
2. Multiply eigenvector $\mathbf{x}_i$ by $e^{\lambda_i t}$.
3. The solution to $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$ is:

$$\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}$$

**Example 2:**

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 1 \\ 0 & 0 & 3 \end{pmatrix}\mathbf{u}, \quad \mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}$$

$\lambda_1 = 1, \lambda_2 = 2, \lambda_3 = 3$ (upper triangular -- eigenvalues on diagonal).

i) $(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$

ii) $(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$

iii) $(A - 3I)\mathbf{x}_3 = \begin{pmatrix} -2 & 1 & 1 \\ 0 & -1 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_3 = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$

iv) $\mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2 \; \mathbf{x}_3)\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix}$

$$\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}^{-1}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}$$

$|X| = 1$, $X^{-1} = \frac{1}{|X|}C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}$

$$\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix}$$

v) $\therefore \mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + c_3 e^{\lambda_3 t}\mathbf{x}_3$

$$= 2e^t\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + 3e^{2t}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} + 4e^{3t}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$$

### 4.5 Exponential of a Matrix

Recall: $e^x = 1 + x + \frac{1}{2}x^2 + \frac{1}{6}x^3 + \cdots + \frac{x^n}{n!} + \cdots = \sum_{i=0}^{\infty}\frac{x^i}{i!}$

$\frac{1}{1-x} = 1 + x + x^2 + x^3 + \cdots + x^n + \cdots = \sum_{i=0}^{\infty}x^i$ (geometric series, $|x| < 1$)

**Matrix exponential:**

$$e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots + \frac{1}{n!}(At)^n + \cdots$$

$(I - At)^{-1} = I + At + (At)^2 + (At)^3 + \cdots + (At)^n + \cdots$

Let $\mathbf{u}(t) = e^{At}\mathbf{u}(0) = e^{At}\mathbf{u}_0$.

Check if $\mathbf{u}$ is the solution to $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$:

$$\frac{d\mathbf{u}}{dt} = \frac{d}{dt}(e^{At}\mathbf{u}_0) = \frac{de^{At}}{dt}\mathbf{u}_0 = \left(A + \frac{1}{2} \cdot 2 \cdot (At)A + \frac{3}{6}(At)^2 A + \cdots + \frac{n}{n!}(At)^{n-1}A + \cdots\right)\mathbf{u}_0$$

$$= A\left(I + At + \frac{1}{2}(At)^2 + \cdots + \frac{1}{(n-1)!}(At)^{n-1} + \cdots\right)\mathbf{u}_0 = Ae^{At}\mathbf{u}_0 = A\mathbf{u}$$

**If $A$ is diagonalizable,** $A = X\Lambda X^{-1}$, then:

$$e^{At} = Xe^{\Lambda t}X^{-1}$$

where

$$e^{\Lambda t} = I + \Lambda t + \frac{1}{2}(\Lambda t)^2 + \frac{1}{6}(\Lambda t)^3 + \cdots = \begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}$$

**Proof:**

$$e^{At} = I + At + \frac{1}{2}(At)^2 + \cdots$$

$$= I + X\Lambda X^{-1}t + \frac{1}{2}t^2(X\Lambda X^{-1})(X\Lambda X^{-1}) + \frac{1}{6}t^3(X\Lambda X^{-1})(X\Lambda X^{-1})(X\Lambda X^{-1}) + \cdots$$

$$= I + t(X\Lambda X^{-1}) + \frac{1}{2}t^2(X\Lambda^2 X^{-1}) + \frac{1}{6}t^3(X\Lambda^3 X^{-1}) + \cdots$$

$$= (XIX^{-1}) + (Xt\Lambda X^{-1}) + (X\frac{t^2}{2}\Lambda^2 X^{-1}) + (X\frac{t^3}{6}\Lambda^3 X^{-1}) + \cdots$$

$$= X\left(I + t\Lambda + \frac{t^2}{2}\Lambda^2 + \frac{t^3}{6}\Lambda^3 + \cdots\right)X^{-1}$$

$$\therefore e^{At} = Xe^{\Lambda t}X^{-1}$$

**Observation:** $A = X\Lambda X^{-1}$ and $e^{At} = Xe^{\Lambda t}X^{-1}$ $\implies$ $e^{At}$ and $A$ share the **same eigenvectors**.

**The solution:**

$$\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\mathbf{u}_0$$

Since $\mathbf{u}_0 = X\mathbf{c}$ (linear combination of eigenvectors):

$$= Xe^{\Lambda t}\underbrace{X^{-1}X}_{I}\mathbf{c} = Xe^{\Lambda t}\mathbf{c} = X\begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{pmatrix} = X\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}$$

$$= (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}$$

$$\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}$$

**Property:** $e^{A(s+t)} = e^{As} \cdot e^{At}$

**Proof:** (Detailed computation using double series and binomial theorem)

$$e^{As} \cdot e^{At} = \left(\sum_{j=0}^{\infty}\frac{(As)^j}{j!}\right)\left(\sum_{k=0}^{\infty}\frac{(At)^k}{k!}\right) = \sum_{j=0}^{\infty}\sum_{k=0}^{\infty}\frac{A^{j+k}s^j t^k}{j!\,k!}$$

Let $n = j + k$, $k = n - j$:

$$= \sum_{n=0}^{\infty}\frac{A^n}{n!}\sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^j t^{n-j} = \sum_{n=0}^{\infty}\frac{A^n}{n!}(s+t)^n = \sum_{n=0}^{\infty}\frac{(A(s+t))^n}{n!} = e^{A(s+t)} \quad \square$$

Using the binomial theorem: $(s + t)^n = \sum_{j=0}^{n}\binom{n}{j}s^{n-j}t^j = \sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^{n-j}t^j$.

From $e^{A(s+t)} = e^{As} \cdot e^{At}$, take $s = -t$:

$$e^{A(-t+t)} = e^0 = I = e^{-At} \cdot e^{At}$$

This implies that $e^{-At}$ is the **inverse** of $e^{At}$.

**Example 4:** Let $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$, what is $e^{At} = ?$

$|A - \lambda I| = \lambda^2 + 1 = 0$, so $\lambda = \pm i$ (antisymmetric matrix).

$A^2 = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}$

$A^3 = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$

$A^4 = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$

$$e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots$$

$$= \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}t + \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\frac{t^2}{2} + \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\frac{t^3}{6} + \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\frac{t^4}{24} + \cdots$$

$$= \begin{pmatrix} 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots & t - \frac{t^3}{6} + \cdots \\ -t + \frac{t^3}{6} - \cdots & 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots \end{pmatrix}$$

$$e^{At} = \begin{pmatrix} \cos(t) & \sin(t) \\ -\sin(t) & \cos(t) \end{pmatrix} \quad \text{(antisymmetric matrix gives rotation)}$$

**Example 5:** Solve $\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}\mathbf{u}$ with $\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$ at $t = 0$.

Upper triangular matrix: $\lambda_1 = 1, \lambda_2 = 2$.

i) $(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

ii) $(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

iii) $X = (\mathbf{x}_1 \; \mathbf{x}_2) = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \to X^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$

$A = X\Lambda X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$

$e^{At} = Xe^{\Lambda t}X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} e^t & \\ & e^{2t} \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$

iv) $\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\begin{pmatrix} 2 \\ 1 \end{pmatrix}$

$= Xe^{\Lambda t}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 2 \\ 1 \end{pmatrix} = Xe^{\Lambda t}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$

$= e^t\begin{pmatrix} 1 \\ 0 \end{pmatrix} + e^{2t}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} e^t + e^{2t} \\ e^{2t} \end{pmatrix}$

### 4.6 Second Order Equations

Consider $my'' + by' + ky = 0$ (mass-spring-damper system).

Take $y = e^{\lambda t}$: $y' = \lambda e^{\lambda t} = \lambda y$, $y'' = (y')' = \lambda^2 y$.

The spring equation becomes: $m\lambda^2 y + b\lambda y + ky = (m\lambda^2 + b\lambda + k)e^{\lambda t} = 0$

$$\implies m\lambda^2 + b\lambda + k = 0$$

$$\iff \lambda^2 + \frac{b}{m}\lambda + \frac{k}{m} = 0 \quad \text{---(*)}$$

$\lambda_1 + \lambda_2 = -\frac{b}{m}$, $\lambda_1\lambda_2 = \frac{k}{m}$

$\lambda = \frac{1}{2}\left(-\frac{b}{m} \pm \sqrt{\frac{b^2}{m^2} - \frac{4k}{m}}\right)$

Let $y_1 = e^{\lambda_1 t}$, $y_2 = e^{\lambda_2 t}$. Complete solution becomes:

$$y(t) = c_1 y_1 + c_2 y_2 \quad \text{if } \lambda_1 \neq \lambda_2 \quad \text{---(**)}$$

**Convert 2nd order ODE to 1st order ODEs:**

$my'' + by' + ky = 0 \iff y'' + \tilde{b}y' + \tilde{k}y = 0$ where $\tilde{b} = \frac{b}{m}$, $\tilde{k} = \frac{k}{m}$.

Let $\mathbf{y} = \begin{pmatrix} y' \\ y \end{pmatrix}$:

$y'' = -\tilde{b}y' - \tilde{k}y = (-\tilde{b} \; -\tilde{k})\begin{pmatrix} y' \\ y \end{pmatrix}$

$y' = y' = (1 \; 0)\begin{pmatrix} y' \\ y \end{pmatrix}$

$$\implies \frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -\tilde{b} & -\tilde{k} \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix} \implies \frac{d\mathbf{y}}{dt} = A\mathbf{y}$$

i) Eigenvalues of $A$:

$|A - \lambda I| = \begin{vmatrix} -\tilde{b} - \lambda & -\tilde{k} \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + \tilde{b}\lambda + \tilde{k} = 0$

The equation for the $\lambda$'s is the same as (*).

ii) Let $\lambda_1, \lambda_2$ be the eigenvalues of $A$:

$(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -\tilde{b} - \lambda_1 & -\tilde{k} \\ 1 & -\lambda_1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 - \lambda_1 x_2 = 0$: $\mathbf{x}_1 = \begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix}$

$(A - \lambda_2 I)\mathbf{x}_2$: $x_1 - \lambda_2 x_2 = 0$: $\mathbf{x}_2 = \begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}$

iii) $\mathbf{y}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2$

$$\begin{pmatrix} y'(t) \\ y(t) \end{pmatrix} = c_1 e^{\lambda_1 t}\begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix} + c_2 e^{\lambda_2 t}\begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}$$

$$\therefore y(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t} \quad \text{like (**)}$$

When $m = 1, b = 0, k = 1$: $y'' + y = 0 \iff y'' = -y$

$$\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -y \\ y' \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}$$

$|A - \lambda I| = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0$, so $\lambda = \pm i$ -- **oscillation**.

### 4.7 Stability of 2 by 2 Matrices

For the solution to $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$, there is a fundamental question:

Does the solution approach $\mathbf{u} = \mathbf{0}$ as $t \to \infty$?

Since the complete solution $\mathbf{u}(t)$ is built from pure solutions $e^{\lambda t}\mathbf{x}$, the stability depends on the eigenvalues of $A$.

$\lambda = r + is$

$e^{\lambda t} = e^{rt}e^{ist}$

$|e^{\lambda t}| = |e^{rt}|$

The **real part** of $\lambda$ controls the growth ($r > 0$) or the decay ($r < 0$).

For any $2 \times 2$ matrix $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$:

Negative $\lambda_1$ and $\lambda_2$ require:

i) $\lambda_1 + \lambda_2 = a + d < 0$ (trace negative)

ii) $\lambda_1 \lambda_2 > 0$ (determinant positive)

### 4.8 Worked Examples

**Example 6:** Solve $y'' + 4y' + 3y = 0$

Substitute $y$ with $e^{\lambda t}$: $(\lambda^2 + 4\lambda + 3)e^{\lambda t} = 0$

$\implies \lambda^2 + 4\lambda + 3 = 0 \implies (\lambda + 3)(\lambda + 1) = 0$

$\therefore \lambda = -1, -3$

$\implies y(t) = c_1 e^{-t} + c_2 e^{-3t}$ -- decaying solution $\to$ **stable** solution.

Introduce $\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}$: $y'' = -4y' - 3y$, $y' = y'$

$$\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -4 & -3 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad \frac{d\mathbf{u}}{dt} = A\mathbf{u}$$

$|A - \lambda I| = \begin{vmatrix} -4 - \lambda & -3 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 4\lambda + 3 = 0$, so $\lambda = -1, -3$.

Corresponding eigenvectors:

$\begin{pmatrix} -4+1 & -3 \\ 1 & +1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}$

$\begin{pmatrix} -4+3 & -3 \\ 1 & +3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} -3 \\ 1 \end{pmatrix}$

The complete solution is:

$$\begin{pmatrix} y' \\ y \end{pmatrix} = \mathbf{u}(t) = c_1 e^{-t}\begin{pmatrix} -1 \\ 1 \end{pmatrix} + c_2 e^{-3t}\begin{pmatrix} -3 \\ 1 \end{pmatrix}$$

This leads to $y(t) = c_1 e^{-t} + c_2 e^{-3t}$.

**Example 7:** $y'' - 2y' + y = 0$

Substitute $y$ with $e^{\lambda t}$: $(\lambda^2 - 2\lambda + 1) = (\lambda - 1)^2 = 0$, so $\lambda = 1, 1$ (repeated root).

$\implies y(t) = c_1 e^t + ?$

Since the root is repeated, the second solution is $c_2 \cdot t e^t$.

$$y(t) = c_1 e^t + c_2 t e^t$$

Introduce $\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}$: $y'' = 2y' - y$, $y' = y'$

$$\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad A = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}$$

$|A - \lambda I| = \begin{vmatrix} 2 - \lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 2\lambda + 1 = 0$, so $\lambda = 1, 1$.

$(A - I)\mathbf{x}_1 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

Since the number of eigenvectors is smaller than 2, diagonalization is **NOT possible**.

So, we compute $e^{At}$ from the definition of a series:

$$e^{At} = I + (At) + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots$$

Write $A = I + (A - I)$:

$$e^{At} = e^{It + (A-I)t} = e^{It} \cdot e^{(A-I)t}$$

Note that $(A - I)^2 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$, so $(A - I)^k = 0$ for $k \geq 2$.

$$e^{(A-I)t} = I + (A-I)t = I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t$$

$$e^{At} = e^t\left(I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t\right) = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}$$

$$\mathbf{u}(t) = e^{At}\mathbf{u}_0$$

$$\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}\begin{pmatrix} y_0' \\ y_0 \end{pmatrix}$$

$$\therefore y(t) = (e^t - te^t)y_0 + te^t y_0'$$

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Eigenvalue equation | $A\mathbf{x} = \lambda\mathbf{x}$, $\mathbf{x} \neq \mathbf{0}$ |
| Characteristic polynomial | $\det(A - \lambda I) = 0$ gives $n$ eigenvalues |
| Trace and determinant | $\text{trace}(A) = \sum \lambda_i$, $\det(A) = \prod \lambda_i$ |
| Powers of eigenvalues | $A^k\mathbf{x} = \lambda^k\mathbf{x}$ |
| Inverse eigenvalue | $A^{-1}\mathbf{x} = \lambda^{-1}\mathbf{x}$ (if $\lambda \neq 0$) |
| Shift rule | $(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}$ |
| Symmetric matrix | Real eigenvalues, orthogonal eigenvectors |
| Skew-symmetric matrix | Pure imaginary or zero eigenvalues |
| Rotation matrix | $\lambda = e^{\pm i\theta}$ |
| Projection matrix | $\lambda = 0$ or $1$ only |
| Diagonalization | $A = X\Lambda X^{-1}$ when $n$ LI eigenvectors exist |
| Matrix powers via diag. | $A^k = X\Lambda^k X^{-1}$ |
| Distinct eigenvalues | Eigenvectors are LI $\implies$ diagonalizable |
| Similar matrices | $C = B^{-1}AB$ has same eigenvalues as $A$ |
| GM vs AM | GM $\leq$ AM; if GM $<$ AM, not diagonalizable |
| Spectral theorem | $S = Q\Lambda Q^T$ for symmetric $S$ |
| Positive definite | $\mathbf{x}^T S\mathbf{x} > 0$ for all $\mathbf{x} \neq \mathbf{0}$; all $\lambda_i > 0$; all pivots $> 0$; all leading determinants $> 0$ |
| Positive semidefinite | $\mathbf{x}^T S\mathbf{x} \geq 0$; allows $\lambda = 0$ |
| Energy test for $A^T A$ | $\mathbf{x}^T A^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0$ |
| PD $\implies$ invertible | $\det(S) > 0$, $S^{-1}$ is also PD |
| PD $\implies S = A^T A$ | With independent columns in $A$ |
| Congruent matrices | $A^T SA$ preserves definiteness type |
| Ellipse equation | $\mathbf{x}^T S\mathbf{x} = 1$; axes along eigenvectors, lengths $1/\sqrt{\lambda}$ |
| Hessian and convexity | $H$ positive definite $\implies$ $f$ is strictly convex |
| Gradient descent | $\mathbf{x}_{k+1} = \mathbf{x}_k - \eta\nabla f(\mathbf{x}_k)$ |
| ODE solution | $\frac{d\mathbf{u}}{dt} = A\mathbf{u} \implies \mathbf{u}(t) = \sum c_i e^{\lambda_i t}\mathbf{x}_i$ |
| Matrix exponential | $e^{At} = I + At + \frac{(At)^2}{2!} + \cdots = Xe^{\Lambda t}X^{-1}$ |
| Stability (continuous) | Stable iff all $\text{Re}(\lambda_i) < 0$ |
| Stability (discrete) | $A^k \to 0$ iff all $|\lambda_i| < 1$ |
| $e^{A(s+t)} = e^{As}e^{At}$ | Matrix exponential additive property |
| $(e^{At})^{-1} = e^{-At}$ | Inverse of matrix exponential |
| 2nd order ODE | $my'' + by' + ky = 0 \iff$ 1st order system $\frac{d\mathbf{y}}{dt} = A\mathbf{y}$ |
| Repeated eigenvalue ODE | Non-diagonalizable case: use $e^{At} = e^{It} \cdot e^{(A-I)t}$ |

---
