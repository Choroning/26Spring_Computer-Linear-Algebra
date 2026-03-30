# Chapter 7 Lecture — Singular Value Decomposition

> **Last Updated:** 2026-03-30

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

$$(A\underset{\sim}{v}_1) \cdot (A\underset{\sim}{v}_2) = 9 \cdot 1 - 3 \cdot 3 = 0 \quad \checkmark$$

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

$$A = \frac{1}{\sqrt{5}}\begin{pmatrix} 5 & 4 \\ 0 & 3 \end{pmatrix} = \frac{1}{\sqrt{10}}\begin{pmatrix} 3 & -1 \\ 1 & 3 \end{pmatrix}\begin{pmatrix} 3 & \\ & 1 \end{pmatrix}\frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$$

The geometric interpretation of $A = U\Sigma V^T$ acting on a vector:

1. **Rotate by $V^T$** -- first rotation by angle $-\phi$ (aligns input to standard basis)
2. **Stretch by $\Sigma$** -- scale along each axis by the singular values
3. **Rotate by $U$** -- second rotation by angle $+\theta$ (rotates to output orientation)

**Geometric picture:**

- A unit circle in the input space is first rotated by $V^T$ to align with the coordinate axes.
- Then $\Sigma$ stretches it into an ellipse (semi-axes of lengths $\sigma_1$ and $\sigma_2$).
- Finally, $U$ rotates the ellipse to its final orientation.

In terms of rotation matrices:

$$A = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \sigma_1 & \\ & \sigma_2 \end{pmatrix}\begin{pmatrix} \cos(-\phi) & -\sin(-\phi) \\ \sin(-\phi) & \cos(-\phi) \end{pmatrix}$$

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

> **Q.** Korean flag? (Much higher rank -- more complex structure, harder to compress.)

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

> **Example:** Image compression using SVD -- search "timbaumann SVD" for an interactive demonstration.

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

**Norm of a Matrix -- extending from vector norms:**

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

$$\text{row 2 avg} = \frac{1}{5}(1 + 1 + 0 + 1 - 1) = 0$$

$$A = A_0 - \begin{pmatrix} 3 \\ 0 \end{pmatrix}(1\;1\;1\;1\;1) = \begin{pmatrix} 2 & 1 & 0 & -1 & -2 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}$$

**ii) Compute the sample covariance matrix.**

$$S = \frac{AA^T}{n - 1}$$

$$AA^T = \begin{pmatrix} 2 & 1 & 0 & -1 & -2 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}\begin{pmatrix} 2 & -1 \\ 1 & 1 \\ 0 & 0 \\ -1 & 1 \\ -2 & -1 \end{pmatrix} = \begin{pmatrix} 10 & 0 \\ 0 & 4 \end{pmatrix}$$

$$S = \frac{1}{4}AA^T = \frac{1}{4}\begin{pmatrix} 10 & 0 \\ 0 & 4 \end{pmatrix} =: S$$

**iii) Find the eigenvalues of $S$.**

$$\lambda_1 = \frac{5}{2}, \quad \lambda_2 = 1$$

$$(S - \lambda I)\underset{\sim}{x}_1 = \begin{pmatrix} 10 - 4\lambda & 0 \\ 0 & 1 - \lambda \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

For $\lambda_1 = \frac{5}{2}$:

$$\begin{pmatrix} 0 & 0 \\ 0 & -\frac{3}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$$\therefore \underset{\sim}{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$$

**iv) What line through the origin is closest to the 5 samples in columns of $A$?**

The centered data matrix:

$$A = \begin{pmatrix} 2 & 1 & 0 & -1 & -2 \\ -1 & 1 & 0 & 1 & -1 \end{pmatrix}$$

Plotting the 5 data points $(2, -1)$, $(1, 1)$, $(0, 0)$, $(-1, 1)$, $(-2, -1)$ in the $xy$-plane:

**The $x$-axis is closer to the 5 points.**

The first singular vector $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$ represents the direction in the data that accounts for the **greatest variability**.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| SVD Decomposition | $A = U\Sigma V^T$ where $U, V$ are orthonormal and $\Sigma$ is diagonal with singular values $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r > 0$ |
| Two Sets of Vectors | Input vectors $\{\underset{\sim}{v}_i\}$: basis for $C(A^T)$ and $N(A)$; Output vectors $\{\underset{\sim}{u}_i\}$: basis for $C(A)$ and $N(A^T)$ |
| SVD vs Eigendecomposition | Eigendecomposition $S = Q\Lambda Q^T$ works for symmetric matrices; SVD $A = U\Sigma V^T$ works for **all** matrices |
| Fundamental Relationship | $A\underset{\sim}{v}_i = \sigma_i \underset{\sim}{u}_i$ maps right singular vectors to scaled left singular vectors |
| Geometry of SVD | $A = U\Sigma V^T$: rotate by $V^T$, stretch by $\Sigma$, rotate by $U$ -- transforms unit circle to ellipse |
| Rank-1 Decomposition | $A = \sigma_1 \underset{\sim}{u}_1 \underset{\sim}{v}_1^T + \sigma_2 \underset{\sim}{u}_2 \underset{\sim}{v}_2^T + \cdots + \sigma_r \underset{\sim}{u}_r \underset{\sim}{v}_r^T$ (sum of rank-1 matrices) |
| Image as Matrix | An image is a large matrix of grayscale values (e.g., 8-bit: 0--255) |
| Image Compression | Discard small singular values $\sigma_i$ for compression with minimal quality loss; rank-1 images (e.g., flags) reduce $N^2 \to 2N$ storage |
| PCA | Uses the largest $\sigma$'s to extract the most important information from data |
| Low-Rank Approximation | $A_k$ (keeping top $k$ singular values) is the closest rank-$k$ matrix to $A$: $\|A - A_k\| \leq \|A - B\|$ for any rank-$k$ matrix $B$ |
| Spectral Norm | $\|A\|_2 = \sigma_1$ (largest singular value) |
| Frobenius Norm | $\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2 + \cdots + \sigma_r^2}$ (root sum of squares of all entries) |
| Nuclear Norm | $\|A\|_N = \sigma_1 + \sigma_2 + \cdots + \sigma_r$ (trace norm, sum of singular values) |
| Sample Covariance | $S = \frac{AA^T}{n-1}$ where $A$ is the centered data matrix |
| First Singular Vector | Represents the direction in the data accounting for the greatest variability |

---
