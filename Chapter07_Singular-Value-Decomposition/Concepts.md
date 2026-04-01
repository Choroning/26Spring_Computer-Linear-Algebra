# Chapter 7 Lecture — Singular Value Decomposition

> **Last Updated:** 2026-03-31
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 7

> **Prerequisites**: [Linear Algebra] Eigenvalues, orthogonality (Ch 1-6).
>
> **Learning Objectives**:
> 1. Derive and compute the Singular Value Decomposition (SVD)
> 2. Interpret SVD geometrically in terms of the four fundamental subspaces
> 3. Apply SVD to low-rank approximation and data analysis

---

<br>

## Table of Contents

- [1. SVD (Singular Value Decomposition)](#1-svd-singular-value-decomposition)
  - [1.1 Motivation: Two Sets of Orthonormal Vectors](#11-motivation-two-sets-of-orthonormal-vectors)
  - [1.2 The Four Fundamental Subspaces in SVD](#12-the-four-fundamental-subspaces-in-svd)
  - [1.3 From Eigendecomposition to SVD](#13-from-eigendecomposition-to-svd)
  - [1.4 The SVD Equation](#14-the-svd-equation)
  - [1.5 Example: Computing SVD](#15-example-computing-svd)
  - [1.6 The Geometry of SVD](#16-the-geometry-of-svd)
  - [1.7 Full Size Form of SVD](#17-full-size-form-of-svd)
  - [1.8 Proof of the SVD](#18-proof-of-the-svd)
  - [1.9 Example: Computing SVD via the Proof Method](#19-example-computing-svd-via-the-proof-method)
  - [1.10 AB and BA: Equal Nonzero Eigenvalues](#110-ab-and-ba-equal-nonzero-eigenvalues)
- [2. Image Processing by Linear Algebra](#2-image-processing-by-linear-algebra)
  - [2.1 Images as Matrices](#21-images-as-matrices)
  - [2.2 Image Compression via Correlated Pixels](#22-image-compression-via-correlated-pixels)
  - [2.3 SVD and Rank-1 Decomposition for Compression](#23-svd-and-rank-1-decomposition-for-compression)
- [3. Principal Component Analysis (PCA)](#3-principal-component-analysis-pca)
  - [3.1 PCA by SVD](#31-pca-by-svd)
  - [3.2 Low-Rank Approximation](#32-low-rank-approximation)
  - [3.3 Norm of a Matrix](#33-norm-of-a-matrix)
  - [3.4 Norm Properties and Inequalities](#34-norm-properties-and-inequalities)
  - [3.5 PCA Example (P 7.3.1)](#35-pca-example-p-731)
  - [3.6 Perpendicular Least Squares](#36-perpendicular-least-squares)
- [Summary](#summary)

---

<br>

## 1. SVD (Singular Value Decomposition)

### 1.1 Motivation: Two Sets of Orthonormal Vectors

The SVD involves finding **two sets of orthonormal vectors**:

**Input vectors** (for the domain, $\mathbb{R}^n$):

$$\{\underset{\sim}{v}_1, \underset{\sim}{v}_2, \dots, \underset{\sim}{v}_r, \dots, \underset{\sim}{v}_n\}$$

- $\underset{\sim}{v}_1, \underset{\sim}{v}_2, \dots, \underset{\sim}{v}_r$ form a **basis for $C(A^T)$** (row space)
- $\underset{\sim}{v}_{r+1}, \dots, \underset{\sim}{v}_n$ form a **basis for $N(A)$** (null space)

**Output vectors** (for the codomain, $\mathbb{R}^m$):

$$\{\underset{\sim}{u}_1, \underset{\sim}{u}_2, \dots, \underset{\sim}{u}_r, \dots, \underset{\sim}{u}_m\}$$

- $\underset{\sim}{u}_1, \underset{\sim}{u}_2, \dots, \underset{\sim}{u}_r$ form a **basis for $C(A)$** (column space)
- $\underset{\sim}{u}_{r+1}, \dots, \underset{\sim}{u}_m$ form a **basis for $N(A^T)$** (left null space)

### 1.2 The Four Fundamental Subspaces in SVD

The dimensions of the four fundamental subspaces:

- $\dim C(A^T) = r$
- $\dim C(A) = r$
- $\dim N(A) = n - r$
- $\dim N(A^T) = m - r$

The SVD maps the input space to the output space:

$$A\underset{\sim}{v}_i = \sigma_i \underset{\sim}{u}_i$$

The matrix $A$ maps each right singular vector $\underset{\sim}{v}_i$ in the row space to a scaled left singular vector $\sigma_i \underset{\sim}{u}_i$ in the column space. The null space vectors are mapped to zero.

Writing this out for each basis vector:

$$\sigma_1 \underset{\sim}{u}_1 = A\underset{\sim}{v}_1$$

$$\sigma_2 \underset{\sim}{u}_2 = A\underset{\sim}{v}_2$$

Combining into matrix form:

$$\begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 \end{pmatrix} \begin{pmatrix} \sigma_1 & \\ & \sigma_2 \end{pmatrix} = A \begin{pmatrix} \underset{\sim}{v}_1 & \underset{\sim}{v}_2 \end{pmatrix}$$

### 1.3 From Eigendecomposition to SVD

Using only $r$ basis vectors:

$$\begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 & \cdots & \underset{\sim}{u}_r \end{pmatrix} \begin{pmatrix} \sigma_1 & & \\ & \sigma_2 & \\ & & \ddots \\ & & & \sigma_r \end{pmatrix} = A \begin{pmatrix} \underset{\sim}{v}_1 & \underset{\sim}{v}_2 & \cdots & \underset{\sim}{v}_r \end{pmatrix}$$

$$U \Sigma = A V$$

where $\sigma_1 \geq \sigma_2 \geq \sigma_3 \geq \cdots \geq \sigma_r > 0$ are the **singular values**.

The dimensions are:

$$\underset{m \times n}{A} \underset{n \times r}{V} = \underset{m \times r}{U} \underset{r \times r}{\Sigma}$$

**Recall** the eigendecomposition: for a matrix with eigenvalues $\lambda_1, \lambda_2, \dots, \lambda_n$ and eigenvectors $\underset{\sim}{x}_1, \underset{\sim}{x}_2, \dots, \underset{\sim}{x}_n$:

$$A(\underset{\sim}{x}_1 \; \underset{\sim}{x}_2 \; \cdots \; \underset{\sim}{x}_n) = (\lambda_1 \underset{\sim}{x}_1 \; \lambda_2 \underset{\sim}{x}_2 \; \cdots \; \lambda_n \underset{\sim}{x}_n) = (\underset{\sim}{x}_1 \; \underset{\sim}{x}_2 \; \cdots \; \underset{\sim}{x}_n) \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix}$$

$$AX = X\Lambda$$

If $A$ is a **symmetric matrix** $S$, then $X = Q$ (orthogonal matrix):

$$SQ = Q\Lambda$$

$$S = Q\Lambda Q^T$$

**Goal:** We want to go beyond symmetric matrices to **all matrices**.

$$S\underset{\sim}{x} = \lambda\underset{\sim}{x} \quad \Longrightarrow \quad A\underset{\sim}{v} = \sigma\underset{\sim}{u}$$

### 1.4 The SVD Equation

From $AV = U\Sigma$, since $V$ is orthogonal ($VV^T = I$):

$$AV = U\Sigma$$

$$AVV^T = U\Sigma V^T$$

$$\boxed{A = U\Sigma V^T}$$

**Every matrix $A$ is diagonalized by two sets of orthonormal vectors.**

This can also be written as a sum of rank-1 matrices:

$$A = \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T + \cdots + \sigma_r \underset{\sim}{u}_r \underset{\sim}{v}_r^T$$

### 1.5 Example: Computing SVD

**Example.** Let $A = \begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix}$.

Then $\text{rank}(A) = 2 = \text{rank}(A^T)$.

**Step 1:** Take orthonormal vectors in $C(A^T)$:

$$\underset{\sim}{v}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, \quad \underset{\sim}{v}_2 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}$$

Check: $\underset{\sim}{v}_1 \cdot \underset{\sim}{v}_2 = 1 - 1 = 0 \implies \underset{\sim}{v}_1 \perp \underset{\sim}{v}_2$.

**Step 2:** Compute $A\underset{\sim}{v}_1$ and $A\underset{\sim}{v}_2$:

$$A\underset{\sim}{v}_1 = \begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 9 \\ 3 \end{pmatrix} = 3\begin{pmatrix} 3 \\ 1 \end{pmatrix} = 3\sqrt{10} \cdot \frac{1}{\sqrt{10}}\begin{pmatrix} 3 \\ 1 \end{pmatrix} \quad \leftarrow \underset{\sim}{u}_1$$

$$A\underset{\sim}{v}_2 = \begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix}\begin{pmatrix} -1 \\ 1 \end{pmatrix} = \begin{pmatrix} -1 \\ 3 \end{pmatrix} = \frac{\sqrt{10}}{\sqrt{10}}\begin{pmatrix} -1 \\ 3 \end{pmatrix} \quad \leftarrow \underset{\sim}{u}_2$$

**Verify orthogonality:**

$$(A\underset{\sim}{v}_1) \cdot (A\underset{\sim}{v}_2) = 9 \cdot (-1) + 3 \cdot 3 = 0 \quad \checkmark$$

**Step 3:** Write in matrix form:

$$A\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 9 & -1 \\ 3 & 3 \end{pmatrix}$$

**Step 4:** Scale $\underset{\sim}{v}_1, \underset{\sim}{v}_2$ by dividing each by $\sqrt{2}$:

$$A \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix} 9 & -1 \\ 3 & 3 \end{pmatrix}$$

So the normalized vectors are:

$$\underset{\sim}{v}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \quad \underset{\sim}{v}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix} -1 \\ 1 \end{pmatrix}$$

This gives:

$$A \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} 3 & -1 \\ 1 & 3 \end{pmatrix} \begin{pmatrix} 3\sqrt{5} & \\ & \sqrt{5} \end{pmatrix}$$

$$AV = U\Sigma$$

where $\sigma_1 = 3\sqrt{5}$, $\sigma_2 = \sqrt{5}$.

**Verify $VV^T = I$:**

$$(\underset{\sim}{v}_1 \; \underset{\sim}{v}_2)\begin{pmatrix} \underset{\sim}{v}_1^T \\ \underset{\sim}{v}_2^T \end{pmatrix} = \underset{\sim}{v}_1\underset{\sim}{v}_1^T + \underset{\sim}{v}_2\underset{\sim}{v}_2^T$$

$$= \frac{1}{2}\begin{pmatrix} 1 \\ 1 \end{pmatrix}(1\;1) + \frac{1}{2}\begin{pmatrix} -1 \\ 1 \end{pmatrix}(-1\;1)$$

$$= \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} + \frac{1}{2}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I \quad \checkmark$$

**Final SVD:**

$$\boxed{A = U\Sigma V^T}$$

$$= \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 \end{pmatrix}\begin{pmatrix} \sigma_1 & \\ & \sigma_2 \end{pmatrix}\begin{pmatrix} \underset{\sim}{v}_1^T \\ \underset{\sim}{v}_2^T \end{pmatrix}$$

$$= \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T$$

### 1.6 The Geometry of SVD

For $A = U\Sigma V^T$:

$$\begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} 3 & -1 \\ 1 & 3 \end{pmatrix}\begin{pmatrix} 3\sqrt{5} & \\ & \sqrt{5} \end{pmatrix}\frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$$

Normalizing the singular values:

$$\frac{1}{\sqrt{5}}\begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} 3 & -1 \\ 1 & 3 \end{pmatrix}\begin{pmatrix} 3 & \\ & 1 \end{pmatrix}\frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$$

The geometric interpretation of $A = U\Sigma V^T$ acting on a vector:

1. **Rotate by $V^T$** — first rotation by angle $-\phi$ (aligns input to standard basis)
2. **Stretch by $\Sigma$** — scale along each axis by the singular values
3. **Rotate by $U$** — second rotation by angle $+\theta$ (rotates to output orientation)

**Geometric picture:**

- A unit circle in the input space is first rotated by $V^T$ to align with the coordinate axes.
- Then $\Sigma$ stretches it into an ellipse (semi-axes of lengths $\sigma_1$ and $\sigma_2$).
- Finally, $U$ rotates the ellipse to its final orientation.

In terms of rotation matrices:

$$A = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \sigma_1 & \\ & \sigma_2 \end{pmatrix}\begin{pmatrix} \cos(-\phi) & -\sin(-\phi) \\ \sin(-\phi) & \cos(-\phi) \end{pmatrix}$$

### 1.7 Full Size Form of SVD

The full size form includes **basis vectors for the nullspaces** of $A$ and $A^T$.

$$V = \begin{pmatrix} \underset{\sim}{v}_1 & \underset{\sim}{v}_2 & \cdots & \underset{\sim}{v}_r & \underset{\sim}{v}_{r+1} & \cdots & \underset{\sim}{v}_n \end{pmatrix}_{n \times n}$$

$$U = \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 & \cdots & \underset{\sim}{u}_r & \underset{\sim}{u}_{r+1} & \cdots & \underset{\sim}{u}_m \end{pmatrix}_{m \times m}$$

The full form equation:

$$\underset{m \times n}{A} \underset{n \times n}{V} = \underset{m \times m}{U} \underset{m \times n}{\Sigma}$$

Since $V$ is square orthogonal: $V^T = V^{-1}$. Similarly $U^T = U^{-1}$.

The null space vectors are mapped to zero:

$$A\underset{\sim}{v}_{r+1} = \underset{\sim}{0}, \quad A\underset{\sim}{v}_{r+2} = \underset{\sim}{0}, \quad \dots, \quad A\underset{\sim}{v}_n = \underset{\sim}{0}$$

So $AV$ has the form:

$$AV = A\begin{pmatrix} \underset{\sim}{v}_1 & \underset{\sim}{v}_2 & \cdots & \underset{\sim}{v}_r & \underset{\sim}{v}_{r+1} & \cdots & \underset{\sim}{v}_n \end{pmatrix} = \begin{pmatrix} \sigma_1 \underset{\sim}{u}_1 & \sigma_2 \underset{\sim}{u}_2 & \cdots & \sigma_r \underset{\sim}{u}_r & \underset{\sim}{0} & \cdots & \underset{\sim}{0} \end{pmatrix}_{m \times n}$$

$$= \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 & \cdots & \underset{\sim}{u}_r & \underset{\sim}{u}_{r+1} & \cdots & \underset{\sim}{u}_m \end{pmatrix}_{m \times m} \begin{pmatrix} \sigma_1 & & & 0 & \cdots & 0 \\ & \sigma_2 & & 0 & \cdots & 0 \\ & & \ddots & \vdots & & \vdots \\ & & & \sigma_r & 0 & \cdots & 0 \\ 0 & 0 & \cdots & 0 & 0 & \cdots & 0 \\ \vdots & \vdots & & \vdots & \vdots & & \vdots \\ 0 & 0 & \cdots & 0 & 0 & \cdots & 0 \end{pmatrix}_{m \times n}$$

**Two cases for the shape of $\Sigma$:**

**i) $m < n$ (short and wide):** The $\Sigma$ matrix is $m \times n$ with $\sigma_1, \dots, \sigma_r$ on the diagonal and zeros elsewhere, forming a wide rectangular matrix.

**ii) $m > n$ (thin and tall):** The $\Sigma$ matrix is $m \times n$ with $\sigma_1, \dots, \sigma_r$ on the diagonal and zeros elsewhere, forming a tall rectangular matrix with extra zero rows at the bottom.

**The full form has a lot of zeros**, which contribute nothing to the matrix multiplication. Therefore, we can take only the first $r$ vectors to get the **reduced form**:

$$\underset{m \times n}{A} \underset{n \times r}{V_r} = \underset{m \times r}{U_r} \underset{r \times r}{\Sigma_r}$$

where $V_r$ contains the basis for $C(A^T)$ and $U_r$ contains the basis for $C(A)$.

Note that:

$$V_r^T V_r = I_r \quad \text{and} \quad U_r^T U_r = I_r$$

**But:**

$$V_r V_r^T \neq I \quad \text{and} \quad U_r U_r^T \neq I$$

The reduced SVD:

$$A = U_r \Sigma_r V_r^T = \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \cdots + \sigma_r \underset{\sim}{u}_r \underset{\sim}{v}_r^T$$

### 1.8 Proof of the SVD

**How do we find $U$, $V$ (singular vectors)?**

Given $A = U\Sigma V^T$ (where $U$ is the **left** singular vectors and $V$ is the **right** singular vectors):

**Computing $A^T A$:**

$$A^T A = (U\Sigma V^T)^T (U\Sigma V^T) = V\Sigma^T U^T U \Sigma V^T = V\Sigma^2 V^T$$

(since $U^T U = I$)

**Computing $AA^T$:**

$$AA^T = (U\Sigma V^T)(U\Sigma V^T)^T = U\Sigma V^T V \Sigma^T U^T = U\Sigma^2 U^T$$

(since $V^T V = I$)

Since $A^T A$ is symmetric:

$$A^T A = Q\Lambda Q^T = V\Sigma^2 V^T$$

Since $AA^T$ is symmetric:

$$AA^T = Q\Lambda Q^T = U\Sigma^2 U^T$$

Therefore $\sigma_1^2, \sigma_2^2, \dots, \sigma_r^2$ are **nonzero eigenvalues of both $A^T A$ and $AA^T$**.

- **Choose $V$ from $A^T A$** (eigenvectors of $A^T A$)
- **Choose $U$ from $AA^T$** (eigenvectors of $AA^T$)

**Steps of the proof:**

**i)** Choose orthonormal eigenvectors $\underset{\sim}{v}_1, \underset{\sim}{v}_2, \dots, \underset{\sim}{v}_r$:

$$A^T A \underset{\sim}{v}_k = \lambda_k \underset{\sim}{v}_k = \sigma_k^2 \underset{\sim}{v}_k$$

$$\therefore \sigma_k = \sqrt{\lambda_k}$$

**ii)** Right singular vector $\underset{\sim}{v}_k$ is connected to left singular vector $\underset{\sim}{u}_k$:

$$A\underset{\sim}{v}_k = \sigma_k \underset{\sim}{u}_k$$

$$\therefore \underset{\sim}{u}_k = \frac{1}{\sigma_k} A\underset{\sim}{v}_k \quad \forall k = 1, 2, \dots, r$$

**iii)** Sanity check — verify $\underset{\sim}{u}_k$ is an eigenvector of $AA^T$:

$$AA^T \underset{\sim}{u}_k = AA^T \cdot \frac{1}{\sigma_k} A\underset{\sim}{v}_k = \frac{1}{\sigma_k} AA^T A \underset{\sim}{v}_k = \frac{1}{\sigma_k} A \sigma_k^2 \underset{\sim}{v}_k = \sigma_k A\underset{\sim}{v}_k = \sigma_k^2 \underset{\sim}{u}_k \quad \square$$

**iv)** Check if $\underset{\sim}{u}_k$ is orthonormal:

$$\underset{\sim}{u}_j^T \underset{\sim}{u}_k = \left(\frac{1}{\sigma_j} A\underset{\sim}{v}_j\right)^T \left(\frac{1}{\sigma_k} A\underset{\sim}{v}_k\right) = \frac{1}{\sigma_j \sigma_k} \underset{\sim}{v}_j^T A^T A \underset{\sim}{v}_k = \frac{\sigma_k}{\sigma_j} \underset{\sim}{v}_j^T \underset{\sim}{v}_k = \begin{cases} 1 & j = k \\ 0 & j \neq k \end{cases}$$

### 1.9 Example: Computing SVD via the Proof Method

**EX1.** Find $U, \Sigma, V$ for $A = \begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix}$.

$\text{rank}(A) = 2 \implies$ two singular values $\sigma_1, \sigma_2$.

**i) Compute $A^T A$ and find eigenvalues:**

$$A^T A = \begin{pmatrix} 5 & 0 \\ 4 & 3 \end{pmatrix}\begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix} = \begin{pmatrix} 25 & 20 \\ 20 & 25 \end{pmatrix}$$

Characteristic equation:

$$\lambda^2 - 50\lambda + 225 = 0$$

$$(\lambda - 45)(\lambda - 5) = 0 \implies \lambda_1 = 45, \; \lambda_2 = 5$$

$$\therefore \sigma_1 = \sqrt{45} = 3\sqrt{5}, \quad \sigma_2 = \sqrt{5}$$

**Find eigenvectors of $A^T A$:**

For $\lambda_1 = 45$:

$$(A^T A - \lambda_1 I)\underset{\sim}{x}_1 = \begin{pmatrix} -20 & 20 \\ 20 & -20 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \underset{\sim}{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$$

For $\lambda_2 = 5$:

$$(A^T A - \lambda_2 I)\underset{\sim}{x}_2 = \begin{pmatrix} 20 & 20 \\ 20 & 20 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \underset{\sim}{x}_2 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}$$

Normalize to get $V$:

$$V = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}$$

**ii) Compute $U$ from $A\underset{\sim}{v}_k = \sigma_k \underset{\sim}{u}_k$:**

$$\underset{\sim}{u}_1 = \frac{1}{\sigma_1} A\underset{\sim}{v}_1 = \frac{1}{3\sqrt{5}}\begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix}\frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{3\sqrt{10}}\begin{pmatrix} 9 \\ 3 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} 3 \\ 1 \end{pmatrix}$$

$$\underset{\sim}{u}_2 = \frac{1}{\sigma_2} A\underset{\sim}{v}_2 = \frac{1}{\sqrt{5}}\begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix}\frac{1}{\sqrt{2}}\begin{pmatrix} -1 \\ 1 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} -1 \\ 3 \end{pmatrix}$$

$$U = \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} 3 & -1 \\ 1 & 3 \end{pmatrix}$$

### 1.10 AB and BA: Equal Nonzero Eigenvalues

For $A_{m \times n}$ and $B_{n \times m}$:

$$AB_{m \times m} \neq BA_{n \times n}$$

**But $AB$ and $BA$ have the same nonzero eigenvalues.**

**Proof:** Let $\lambda$ be an eigenvalue of $AB$, with $\underset{\sim}{x}$ the corresponding eigenvector:

$$AB\underset{\sim}{x} = \lambda\underset{\sim}{x}$$

$$B(AB\underset{\sim}{x}) = B(\lambda\underset{\sim}{x}) = \lambda(B\underset{\sim}{x})$$

$$BA(B\underset{\sim}{x}) = \lambda(B\underset{\sim}{x})$$

Therefore $\lambda$ is an eigenvalue of $BA$, and $B\underset{\sim}{x}$ is the corresponding eigenvector.

Both $AB$ and $BA$ have the same $\lambda$.

---

<br>

## 2. Image Processing by Linear Algebra

### 2.1 Images as Matrices

An image is a **large matrix of grayscale values**.

- Each pixel has a value, e.g., **8-bit** (values from $0$ to $255$).
- The matrix entries represent the intensity at each pixel location.

### 2.2 Image Compression via Correlated Pixels

When **nearby pixels are correlated**, the image can be **compressed**.

**Example: French flag.**

The French flag image (blue, white, red vertical stripes) can be represented as a matrix:

$$\begin{pmatrix} b & b & w & w & r & r \\ b & b & w & w & r & r \\ b & b & w & w & r & r \\ b & b & w & w & r & r \\ b & b & w & w & r & r \\ b & b & w & w & r & r \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \\ 1 \\ 1 \end{pmatrix} \begin{pmatrix} b & b & w & w & r & r \end{pmatrix}$$

This is a **rank 1 matrix**.

Storage reduction: Replace $N^2$ entries with $2N$ entries (one column vector + one row vector).

> **Q.** Korean flag? (Much higher rank — more complex structure, harder to compress.)

### 2.3 SVD and Rank-1 Decomposition for Compression

**Singular Values with Diagonals:**

$$A = U\Sigma V^T$$

$$= \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 & \cdots & \underset{\sim}{u}_r \end{pmatrix}\begin{pmatrix} \sigma_1 & & \\ & \sigma_2 & \\ & & \ddots \\ & & & \sigma_r \end{pmatrix}\begin{pmatrix} \underset{\sim}{v}_1^T \\ \underset{\sim}{v}_2^T \\ \vdots \\ \underset{\sim}{v}_r^T \end{pmatrix}$$

$$= \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 & \cdots & \underset{\sim}{u}_r \end{pmatrix}\begin{pmatrix} \sigma_1 \underset{\sim}{v}_1^T \\ \sigma_2 \underset{\sim}{v}_2^T \\ \vdots \\ \sigma_r \underset{\sim}{v}_r^T \end{pmatrix}$$

$$= \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T + \cdots + \sigma_r \underset{\sim}{u}_r \underset{\sim}{v}_r^T$$

This is a **sum of rank-1 matrices**.

**Compression principle:** In compression, the **small $\sigma$'s can be discarded** with **no serious loss in image quality**.

$$A = U\Sigma V^T$$

$$= \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 & \cdots & \underset{\sim}{u}_r \end{pmatrix}\begin{pmatrix} \sigma_1 & & \\ & \sigma_2 & \\ & & \ddots \\ & & & \sigma_r \end{pmatrix}\begin{pmatrix} \underset{\sim}{v}_1^T \\ \underset{\sim}{v}_2^T \\ \vdots \\ \underset{\sim}{v}_r^T \end{pmatrix}$$

$$\approx \begin{pmatrix} \underset{\sim}{u}_1 & \underset{\sim}{u}_2 \end{pmatrix}\begin{pmatrix} \sigma_1 & \\ & \sigma_2 \end{pmatrix}\begin{pmatrix} \underset{\sim}{v}_1^T \\ \underset{\sim}{v}_2^T \end{pmatrix} = \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T$$

> **Example:** Image compression using SVD — search "timbaumann SVD" for an interactive demonstration.

---

<br>

## 3. Principal Component Analysis (PCA)

### 3.1 PCA by SVD

$$A = U\Sigma V^T$$

**PCA** uses the **largest $\sigma$'s** connected to the first $\underset{\sim}{u}$'s and $\underset{\sim}{v}$'s to understand the information in a matrix of data.

### 3.2 Low-Rank Approximation

We extract the most important part $A_k$:

$$A_k = \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T + \cdots + \sigma_k \underset{\sim}{u}_k \underset{\sim}{v}_k^T$$

$$\text{rank}(A_k) = k$$

$$A \approx A_k$$

- $A_k$ captures the **largest variance** in the data.
- $A_k$ is the **closest rank-$k$ matrix** to $A$.
- If $B$ has rank $k$, then:

$$\|A - A_k\| \leq \|A - B\|$$

### 3.3 Norm of a Matrix

**Three types of matrix norms:**

**i) Spectral norm (2-norm):**

$$\|A\|_2 = \max_{\underset{\sim}{x}} \frac{\|A\underset{\sim}{x}\|}{\|\underset{\sim}{x}\|} = \sigma_1$$

**ii) Frobenius norm:**

$$\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2 + \cdots + \sigma_r^2}$$

**iii) Nuclear norm (trace norm):**

$$\|A\|_N = \sigma_1 + \sigma_2 + \cdots + \sigma_r$$

**Example with the identity matrix** $I \in \mathbb{R}^{n \times n}$:

$$\|I\|_2 = \max \frac{\|I\underset{\sim}{x}\|}{\|\underset{\sim}{x}\|} = \max \frac{\|\underset{\sim}{x}\|}{\|\underset{\sim}{x}\|} = 1$$

$$\|I\|_F = \sqrt{1 + 1 + \cdots + 1} = \sqrt{n}$$

$$\|I\|_N = 1 + 1 + \cdots + 1 = n$$

**Example with SVD:**

For $A = U\Sigma V^T$:

$$\|A\|_2 = \|U\|_2 \|\Sigma\|_2 \|V^T\|_2$$

Since $U, V$ are orthonormal, the magnitude of $\underset{\sim}{x}$ is **NOT changed**:

$$\|U\|_2 = \max \frac{\|U\underset{\sim}{x}\|}{\|\underset{\sim}{x}\|} = 1$$

$$\|V^T\|_2 = 1$$

Therefore:

$$\|A\|_2 = \|\Sigma\|_2 = \sigma_1$$

### 3.4 Norm Properties and Inequalities

**Norm of a Matrix — extending from vector norms:**

Length of a vector $\underset{\sim}{v}$:

$$\|\underset{\sim}{v}\|^2 = v_1^2 + v_2^2 + \cdots + v_n^2$$

Extend the idea to a matrix:

$$\|A\|_F^2 = a_{11}^2 + a_{12}^2 + \cdots + a_{1n}^2 + \cdots + a_{mn}^2 \quad \text{(Frobenius Norm)}$$

**Basic properties:**

- $\|\underset{\sim}{v}\| \geq 0$ and $\|c\underset{\sim}{v}\| = |c| \|\underset{\sim}{v}\|$
- $\|A\|_F \geq 0$ and $\|cA\|_F = |c| \|A\|_F$

**Schwarz inequality:**

$$|\underset{\sim}{v}^T \underset{\sim}{w}| \leq \|\underset{\sim}{v}\| \|\underset{\sim}{w}\|$$

$$\|AB\|_F \leq \|A\|_F \|B\|_F$$

**Triangular inequality:**

$$\|\underset{\sim}{v} + \underset{\sim}{w}\| \leq \|\underset{\sim}{v}\| + \|\underset{\sim}{w}\|$$

$$\|A + B\|_F \leq \|A\|_F + \|B\|_F$$

### 3.5 PCA Example (P 7.3.1)

$A_0$ holds **2 measurements** of **5 samples**:

$$A_0 = \begin{pmatrix} 5 & 4 & 3 & 2 & 1 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}$$

**i) Find the average of each row and subtract it to produce the centered $A$.**

$$\text{row 1 avg} = \frac{1}{5}(5 + 4 + 3 + 2 + 1) = 3$$

$$\text{row 2 avg} = \frac{1}{5}(-1 + 1 + 0 + 1 + (-1)) = 0$$

$$A = A_0 - \begin{pmatrix} 3 \\ 0 \end{pmatrix}(1\;1\;1\;1\;1) = \begin{pmatrix} 2 & 1 & 0 & -1 & -2 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}$$

**ii) Compute the sample covariance matrix.**

$$S = \frac{AA^T}{n - 1}$$

$$AA^T = \begin{pmatrix} 2 & 1 & 0 & -1 & -2 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}\begin{pmatrix} 2 & -1 \\ 1 & 1 \\ 0 & 0 \\ -1 & 1 \\ -2 & -1 \end{pmatrix} = \begin{pmatrix} 10 & 0 \\ 0 & 4 \end{pmatrix}$$

$$S = \frac{1}{4}AA^T = \frac{1}{4}\begin{pmatrix} 10 & 0 \\ 0 & 4 \end{pmatrix} =: S$$

**iii) Find the eigenvalues of $S$.**

$$\lambda_1 = \frac{5}{2}, \quad \lambda_2 = 1$$

$$(S - \lambda I)\underset{\sim}{x}_1 = \begin{pmatrix} \frac{5}{2} - \lambda & 0 \\ 0 & 1 - \lambda \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

For $\lambda_1 = \frac{5}{2}$:

$$\begin{pmatrix} 0 & 0 \\ 0 & -\frac{3}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$$\therefore \underset{\sim}{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$$

**iv) What line through the origin is closest to the 5 samples in columns of $A$?**

The centered data matrix:

$$A = \begin{pmatrix} 2 & 1 & 0 & -1 & -2 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}$$

Plotting the 5 data points $(2, -1)$, $(1, 1)$, $(0, 0)$, $(-1, 1)$, $(-2, -1)$ in the $xy$-plane:

**The $x$-axis is closer to the 5 points.**

The first singular vector $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$ represents the direction in the data that accounts for the **greatest variability**.

### 3.6 Perpendicular Least Squares

Recall that the **least square solution** to $A\underset{\sim}{x} = \underset{\sim}{b}$ minimizes $\|A\underset{\sim}{\hat{x}} - \underset{\sim}{b}\|^2$, where the errors are defined by:

$$\underset{\sim}{e} = A\underset{\sim}{\hat{x}} - \underset{\sim}{b}$$

An error is the **vertical distance** from each data point to the fitted line.

In contrast, **perpendicular least squares** measures the **perpendicular distances** from data points to the line.

The sum of squared **perpendicular distances** from the data points to the $\underset{\sim}{u}_1$ line is a **minimum**.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| SVD Decomposition | $A = U\Sigma V^T$ where $U, V$ are orthonormal and $\Sigma$ is diagonal with singular values $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r > 0$ |
| Two Sets of Vectors | Input vectors $\{\underset{\sim}{v}_i\}$: basis for $C(A^T)$ and $N(A)$; Output vectors $\{\underset{\sim}{u}_i\}$: basis for $C(A)$ and $N(A^T)$ |
| SVD vs Eigendecomposition | Eigendecomposition $S = Q\Lambda Q^T$ works for symmetric matrices; SVD $A = U\Sigma V^T$ works for **all** matrices |
| Fundamental Relationship | $A\underset{\sim}{v}_i = \sigma_i \underset{\sim}{u}_i$ maps right singular vectors to scaled left singular vectors |
| Geometry of SVD | $A = U\Sigma V^T$: rotate by $V^T$, stretch by $\Sigma$, rotate by $U$ — transforms unit circle to ellipse |
| Full vs Reduced SVD | Full form: $A_{m \times n} V_{n \times n} = U_{m \times m} \Sigma_{m \times n}$; Reduced form keeps only $r$ vectors: $A V_r = U_r \Sigma_r$ |
| Rank-1 Decomposition | $A = \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T + \cdots + \sigma_r \underset{\sim}{u}_r \underset{\sim}{v}_r^T$ (sum of rank-1 matrices) |
| Finding Singular Vectors | $V$ from eigenvectors of $A^T A$; $U$ from eigenvectors of $AA^T$; $\sigma_k = \sqrt{\lambda_k}$; $\underset{\sim}{u}_k = \frac{1}{\sigma_k} A\underset{\sim}{v}_k$ |
| $A^T A$ and $AA^T$ | $A^T A = V\Sigma^2 V^T$; $AA^T = U\Sigma^2 U^T$; both share nonzero eigenvalues $\sigma_1^2, \dots, \sigma_r^2$ |
| AB and BA Eigenvalues | $AB$ and $BA$ have the same nonzero eigenvalues (different sizes but same $\lambda$'s) |
| Image as Matrix | An image is a large matrix of grayscale values (e.g., 8-bit: 0--255) |
| Image Compression | Discard small singular values $\sigma_i$ for compression with minimal quality loss; rank-1 images (e.g., flags) reduce $N^2 \to 2N$ storage |
| PCA | Uses the largest $\sigma$'s to extract the most important information from data |
| Low-Rank Approximation | $A_k$ (keeping top $k$ singular values) is the closest rank-$k$ matrix to $A$: $\|A - A_k\| \leq \|A - B\|$ for any rank-$k$ matrix $B$ |
| Spectral Norm | $\|A\|_2 = \sigma_1$ (largest singular value) |
| Frobenius Norm | $\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2 + \cdots + \sigma_r^2}$ (root sum of squares of all entries) |
| Nuclear Norm | $\|A\|_N = \sigma_1 + \sigma_2 + \cdots + \sigma_r$ (trace norm, sum of singular values) |
| Sample Covariance | $S = \frac{AA^T}{n-1}$ where $A$ is the centered data matrix |
| First Singular Vector | Represents the direction in the data accounting for the greatest variability |
| Perpendicular Least Squares | PCA minimizes perpendicular distances (not vertical) from data points to the $\underset{\sim}{u}_1$ line |

---
