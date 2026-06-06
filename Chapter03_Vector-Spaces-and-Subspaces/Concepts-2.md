# Chapter 3 Lecture — Vector Spaces and Subspaces — Part 2/2

> **📑 This document is split into 2 parts.** GitHub renders only a limited number of math expressions per file, so it is split by section so that all math displays correctly.
>
> [Part 1](Concepts.md) · **Part 2**

---

<br>

## Table of Contents

- [5. Independence, Basis, and Dimension (3.4)](#5-independence-basis-and-dimension-34)
  - [5.1 Independent Vectors](#51-independent-vectors)
  - [5.2 Linear Independence via the Nullspace](#52-linear-independence-via-the-nullspace)
  - [5.3 Vectors that Span a Subspace](#53-vectors-that-span-a-subspace)
  - [5.4 Basis for a Vector Space](#54-basis-for-a-vector-space)
  - [5.5 Dimension of a Vector Space](#55-dimension-of-a-vector-space)
  - [5.6 Bases for Matrix Spaces and Function Spaces](#56-bases-for-matrix-spaces-and-function-spaces)
- [6. Dimensions of the Four Subspaces (3.5)](#6-dimensions-of-the-four-subspaces-35)
  - [6.1 Dimension Summary](#61-dimension-summary)
  - [6.2 Orthogonality of the Subspaces](#62-orthogonality-of-the-subspaces)
  - [6.3 The Four Subspaces for R_0](#63-the-four-subspaces-for-r_0)
  - [6.4 Relationship Between A and R_0](#64-relationship-between-a-and-r_0)
  - [6.5 The Fundamental Theorem of Linear Algebra](#65-the-fundamental-theorem-of-linear-algebra)
- [Summary](#summary)

---

<br>

## 5. Independence, Basis, and Dimension (3.4)

### 5.1 Independent Vectors

**(1) Independent vectors:** The only zero combination

```math
c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \cdots + c_n \mathbf{v}_n = \mathbf{0}
```

has all $`c_1 = c_2 = \cdots = c_n = 0`$.

This implies that if at least one of the scalars is nonzero, let's say $`c_1 \neq 0`$, then we can write:

```math
\mathbf{v}_1 = -\frac{c_2}{c_1}\mathbf{v}_2 - \frac{c_3}{c_1}\mathbf{v}_3 - \cdots - \frac{c_n}{c_1}\mathbf{v}_n
```

$`\mathbf{v}_1`$ is a linear combination of other vectors.

**(2) The vectors $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$ span the space $`\mathbb{S}`$ if $`\mathbb{S}`$ = all combinations of $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$.**

e.g., $`\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$ $`\implies`$ $`\begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}`$. $`\hat{i}, \hat{j}`$ span the space $`\mathbb{R}^2`$.

**(3)** The vectors $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$ are a **basis** for $`\mathbb{S}`$ if they are **linearly independent** and they **span** $`\mathbb{S}`$.

Every vector in the space is a **unique** combination of the basis vectors.

e.g., $`\mathbb{R}^2 \ni \begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}`$

**(4)** The **dimension** of a vector space $`\mathbb{S}`$ is the number $`n`$ of vectors in every basis for $`\mathbb{S}`$.

Consider $`A \in \mathbb{R}^{m \times n}`$. There are $`n`$ columns, of which $`r`$ are independent, meaning the remaining $`n - r`$ columns are dependent. The dimension of $`C(A)`$ is $`r`$, which is the rank of $`A`$.

Four essential ideas in this section are:
1. Independent Vectors
2. Spanning a Space
3. Basis for a Space
4. Dimension of a Space

### 5.2 Linear Independence via the Nullspace

**Def.** The columns of $`A`$ are **linearly independent** when the only solution to $`A\mathbf{x} = \mathbf{0}`$ is $`\mathbf{x} = \mathbf{0}`$.

```math
\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n = \mathbf{0}
```

iff $`x_1 = x_2 = \cdots = x_n = 0`$.

**Geometric interpretation:**

e.g., $`\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3 \in \mathbb{R}^3`$ are not in the same plane $`\implies`$ those vectors are independent. No combination of $`\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3`$ that equals $`\mathbf{0}`$ exists (except the trivial one).

e.g., $`\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3`$ are in the same plane in $`\mathbb{R}^3`$. They are dependent. For instance, $`\mathbf{w}_3 = \mathbf{w}_1 + \mathbf{w}_2`$ $`\iff`$ $`1 \cdot \mathbf{w}_1 + 1 \cdot \mathbf{w}_2 - 1 \cdot \mathbf{w}_3 = \mathbf{0}`$. The combination gives $`\mathbf{0}`$, but there are nonzero coefficients, meaning they are dependent.

**Quick checks for dependence:**

- Q. Are $`\begin{pmatrix} 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 10^{-5} \end{pmatrix}`$ dependent? **No.**
- Q. Are $`\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} -1 \\ -1 \end{pmatrix}`$ dependent? **Yes.**
- Q. Are $`\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$ dependent? **Yes.**
- Q. In $`\mathbb{R}^2`$, any three vectors are dependent. **True.**

**Ex 1.** Let $`A = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix}`$. Does $`A`$ have dependent columns?

Check the nullspace of $`A`$, $`N(A)`$.

```math
A\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`\xrightarrow{R_2 - 2R_1, \; R_3 - R_1}`$ $`\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -1 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$, $`\text{rank}(A) = 2 < 3`$.

Take $`x_3 = 1`$: $`x_1 + 3 = 0 \Rightarrow x_1 = -3`$, $`x_2 - 1 = 0 \Rightarrow x_2 = 1`$.

The nonzero coefficients result in the zero vector:

```math
-3\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + 1\begin{pmatrix} 3 \\ 5 \\ 3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

If the columns of $`A`$ are independent, then $`\text{rank}(A) = r = n`$ and the nullspace of $`A`$ has only the zero vector: $`N(A) = \lbrace\mathbf{0}\rbrace`$.

**Any set of $`n`$ vectors in $`\mathbb{R}^m`$ must be linearly dependent if $`n > m`$.**

e.g., Suppose you have 7 columns with 5 components ($`5 \times 7`$ matrix). 7 column vectors are from $`\mathbb{R}^5`$. There cannot be more than 5 pivots in 5 rows. $`A\mathbf{x} = \mathbf{0}`$ has at least $`2 \;(= 7 - 5)`$ free variables. That is, it has nonzero solutions.

**Note:** The columns might be dependent or independent if $`n \leq m`$. Elimination will reveal the $`r`$ pivots.

### 5.3 Vectors that Span a Subspace

Let $`A`$ be a matrix. $`C(A)`$ is the column space which consists of all combinations of $`A\mathbf{x}`$:

```math
\begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{v}_1 + x_2 \mathbf{v}_2 + \cdots + x_n \mathbf{v}_n
```

**e.g.,** $`A = \begin{pmatrix} 1 & 4 \\ 2 & 7 \\ 3 & 5 \end{pmatrix}`$, $`A^T = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 7 & 5 \end{pmatrix}`$

$`C(A)`$ is the plane in $`\mathbb{R}^3`$. $`C(A^T)`$, row space of $`A`$, is $`\mathbb{R}^2`$.

For $`A \in \mathbb{R}^{m \times n}`$: the row vectors are in $`\mathbb{R}^n`$, the column vectors are in $`\mathbb{R}^m`$.

### 5.4 Basis for a Vector Space

**Def.** A **basis** for a vector space is a sequence of vectors with two properties:
> i) The basis vectors are **linearly independent**.
> ii) They **span** the space.

**e.g.,** $`\begin{pmatrix} x \\ y \end{pmatrix} = x\begin{pmatrix} 1 \\ 0 \end{pmatrix} + y\begin{pmatrix} 0 \\ 1 \end{pmatrix}`$ where $`\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$.

i) $`\hat{i}`$ and $`\hat{j}`$ are LI.
ii) $`\hat{i}, \hat{j}`$ span $`\mathbb{R}^2`$.

The combination $`\mathbf{x} = x\hat{i} + y\hat{j}`$ is unique because $`\hat{i}, \hat{j}`$ are LI.

**There is one and only one way to write $`\mathbf{v}`$ as a combination of the basis vectors.**

**Proof.** Suppose $`\mathbf{v} = a_1 \mathbf{v}_1 + \cdots + a_n \mathbf{v}_n`$ and $`\mathbf{v} = b_1 \mathbf{v}_1 + \cdots + b_n \mathbf{v}_n`$.

By subtraction, we have the zero vector:

```math
\mathbf{0} = (a_1 - b_1)\mathbf{v}_1 + \cdots + (a_n - b_n)\mathbf{v}_n
```

From the LI of $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$: $`a_1 - b_1 = 0`$, $`a_2 - b_2 = 0`$, $`\ldots`$, $`a_n - b_n = 0`$.

Hence $`a_i = b_i`$ for $`i = 1, 2, \ldots, n`$. $`\square`$

---

**Ex 3.** The columns of $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}`$ produce the "standard basis" for $`\mathbb{R}^2`$.

Let $`\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$.

i) $`\hat{i}, \hat{j}`$ are LI.
ii) $`\hat{i}, \hat{j}`$ span $`\mathbb{R}^2`$.

$`\therefore \hat{i}, \hat{j}`$ are basis vectors in $`\mathbb{R}^2`$.

---

**Ex 4.** The columns of every invertible $`n`$ by $`n`$ matrix give a basis for $`\mathbb{R}^n`$.

**Invertible matrix example:**

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = R
```

3 LI column vectors. $`C(A) = \mathbb{R}^3`$. 3 nonzero pivots. $`\text{rank}(A) = 3`$. $`N(A) = \lbrace\mathbf{0}\rbrace`$.

**Note that the basis is NOT unique.**

The vectors $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$ are a basis for $`\mathbb{R}^n`$ when they are the columns of an $`n`$ by $`n`$ invertible matrix. $`\mathbb{R}^n`$ has infinitely many different bases.

---

**Singular matrix example:**

```math
B = \begin{pmatrix} 1 & 0 & 1 \\ 1 & 1 & 2 \\ 1 & 1 & 2 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 1 & 2 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 0 & 0 \end{pmatrix} = R_0
```

$`C(A) \neq \mathbb{R}^3`$. $`\text{rank}(A) = 2 < 3`$.

When columns are dependent, we only keep the **pivot columns**.

e.g., 1st, 2nd columns in $`B`$.

- Every set of independent vectors can be **extended** to a basis. (e.g., 1st, 2nd columns in $`B`$ span a plane in $`\mathbb{R}^3`$.)
- Every spanning set of vectors can be **reduced** to a basis.

---

**Ex 5.** $`A = \begin{pmatrix} 2 & 4 \\ 3 & 6 \end{pmatrix}`$

$`\xrightarrow{R_2 - \frac{3}{2}R_1}`$ $`\begin{pmatrix} 2 & 4 \\ 0 & 0 \end{pmatrix}`$ $`\xrightarrow{R_1 / 2}`$ $`\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0`$

$`\text{rank}(A) = 1 < 2`$. One pivot column, one pivot row.

---

**Ex 6.** $`R_0 = \begin{pmatrix} 1 & 2 & 0 & 3 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix}`$

1st, 3rd columns are pivot columns. $`\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}`$ are basis vectors for $`C(R_0)`$.

$`C(R_0)`$ is the $`xy`$ plane in $`\mathbb{R}^3`$.

Also, 2nd and 3rd column vectors are a basis of $`C(R_0)`$.

---

**All bases for a vector space contain the same number of vectors.** The number of vectors in any and every basis is the "**dimension**" of the space.

### 5.5 Dimension of a Vector Space

**Dimension of a Vector Space:**

If $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m`$ and $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n`$ are both bases for the same vector space, then $`m = n`$.

**Proof.** Let $`n > m`$. Since $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m`$ are a basis, each $`\mathbf{w}_i`$ for $`i = 1, 2, \ldots, n`$ must be a combination of the $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m`$:

```math
\mathbf{w}_1 = a_{11}\mathbf{v}_1 + a_{21}\mathbf{v}_2 + \cdots + a_{m1}\mathbf{v}_m
```
```math
\mathbf{w}_2 = a_{12}\mathbf{v}_1 + a_{22}\mathbf{v}_2 + \cdots + a_{m2}\mathbf{v}_m
```
```math
\vdots
```
```math
\mathbf{w}_n = a_{1n}\mathbf{v}_1 + a_{2n}\mathbf{v}_2 + \cdots + a_{mn}\mathbf{v}_m
```

This leads to:

```math
W = (\mathbf{w}_1 \quad \mathbf{w}_2 \quad \cdots \quad \mathbf{w}_n) = (\mathbf{v}_1 \quad \mathbf{v}_2 \quad \cdots \quad \mathbf{v}_m) \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} = VA
```

$`A`$ is $`m`$ by $`n`$ matrix (short and wide). Because $`n > m`$, $`A\mathbf{x} = \mathbf{0}`$ has nonzero solutions.

From $`A\mathbf{x} = \mathbf{0}`$: $`VA\mathbf{x} = V\mathbf{0} = \mathbf{0}`$.

That is $`W\mathbf{x} = \mathbf{0}`$, i.e., $`x_1\mathbf{w}_1 + x_2\mathbf{w}_2 + \cdots + x_n\mathbf{w}_n = \mathbf{0}`$.

Since $`\mathbf{x}`$ is a nonzero vector, $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n`$ are not LI. $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n`$ could NOT be a basis. This contradicts $`\mathbf{w}_i`$ for $`i = 1, 2, \ldots, n`$ is a basis. $`\square`$

---

**Def.** The **dimension** of a space is the number of vectors in every basis.

**e.g.,** The line through $`\mathbf{u} = \begin{pmatrix} 1 \\ 5 \\ 2 \end{pmatrix}`$ has 1 dimension.

Perpendicular to that line is the plane: $`\mathbf{u} \cdot \mathbf{x} = 0 \iff x + 5y + 2z = 0`$.

The plane null space of the matrix $`A = \begin{pmatrix} 1 & 5 & 2 \end{pmatrix}`$:

$`(1 \quad 5 \quad 2)\begin{pmatrix} x \\ y \\ z \end{pmatrix} = 0`$

$`n = 3`$, $`\text{rank}(A) = r = 1`$, $`n - r = 2`$ free variables.

$`n - r`$ special solutions give a basis for the nullspace: dimension $`n - r`$.

To find the special solutions:
- Take $`y = 1, z = 0 \implies x = -5`$
- Take $`y = 0, z = 1 \implies x = -2`$

$`\begin{pmatrix} -5 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} -2 \\ 0 \\ 1 \end{pmatrix}`$ are a basis (2 dimension).

**Summary of dimensions:**
- The row space has dimension $`r`$
- The column space has dimension $`r`$
- The nullspace has dimension $`n - r`$
- $`N(A^T)`$ has dimension $`m - r`$

### 5.6 Bases for Matrix Spaces and Function Spaces

**Bases for Matrix Spaces and Function Spaces**

The words "independence," "basis," "dimension" are not limited to column vectors.

**Matrix Spaces:**

The vector space $`\mathbb{M}`$ contains all $`2 \times 2`$ matrices. Its dimension is 4.

e.g., $`A_1, A_2, A_3, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

i) They are LI: $`x_1 A_1 + x_2 A_2 + x_3 A_3 + x_4 A_4 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}`$ iff $`x_1 = x_2 = x_3 = x_4 = 0`$.

ii) They span $`\mathbb{M}`$: $`\begin{pmatrix} a & b \\ c & d \end{pmatrix} = aA_1 + bA_2 + cA_3 + dA_4`$. The linear combinations of $`A_1, A_2, A_3, A_4`$ can produce any matrix in $`\mathbb{M}`$.

**Subspaces of matrix spaces:**

- $`A_1, A_2, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$ are a basis for the subspace $`\mathcal{U}`$ of upper triangular matrices: $`\mathcal{U} \ni \begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

- $`A_1, A_4`$ are a basis for the subspace $`\mathbb{D}`$ of diagonal matrices: $`\mathbb{D} \ni \begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

- $`A_1, A_4, A_2 + A_3`$ are a basis for symmetric matrices: $`\mathbb{S} \ni \begin{pmatrix} a & b \\ b & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

**Dimensions of matrix subspaces (for $`n \times n`$ matrices):**

| Space | Dimension |
|:------|:----------|
| Whole $`n \times n`$ matrix space | $`n^2`$ |
| Diagonal matrices | $`n`$ |
| Upper triangular matrices | $`\frac{1}{2}n^2 + \frac{1}{2}n`$ |
| Symmetric matrices | $`\frac{1}{2}n^2 + \frac{1}{2}n`$ |

For upper triangular matrices: Total number of zeros $`= \sum_{i=1}^{n}(i-1) = \sum_{i=1}^{n} i - n = \frac{n(n+1)}{2} - n`$. Total number of nonzeros $`= n^2 - \frac{n(n+1)}{2} + n = \frac{n^2}{2} + \frac{n}{2}`$.

The number of entries (dimension) for upper triangular and symmetric matrices is the same: $`\frac{1}{2}n^2 + \frac{1}{2}n`$.

---

**Function Spaces:**

```math
\frac{d^2y}{dx^2} = 0, \quad \frac{d^2y}{dx^2} = -y, \quad \frac{d^2y}{dx^2} = y
```

These involve the 2nd derivative. In calculus, we find the solution $`y(x)`$:

| ODE | Solution | Basis | Dimension |
|:----|:---------|:------|:----------|
| $`y'' = 0`$ | $`y = cx + d`$ | $`\lbrace1, x\rbrace`$ | 2 |
| $`y'' = -y`$ | $`y = c\sin x + d\cos x`$ | $`\lbrace\sin x, \cos x\rbrace`$ | 2 |
| $`y'' = y`$ | $`y = ce^x + de^{-x}`$ | $`\lbracee^x, e^{-x}\rbrace`$ | 2 |

The basis vectors are in the **nullspace** of the 2nd derivative.

$`y'' = 2`$ has the particular solution $`y_p = x^2`$: $`\frac{dy_p}{dx} = 2x`$, $`\frac{d^2 y_p}{dx^2} = 2`$.

Therefore, the general (complete) solution becomes:

```math
y(x) = y_p(x) + y_n(x) = x^2 + cx + d
```

---

<br>

## 6. Dimensions of the Four Subspaces (3.5)

### 6.1 Dimension Summary

1. The column space $`C(A)`$ and the row space $`C(A^T)`$ both have dimension $`r`$ (= rank of $`A`$).
2. The nullspace of $`A`$, $`N(A)`$, has dimension $`n - r`$.
3. The left nullspace of $`A`$, $`N(A^T)`$, has dimension $`m - r`$.
4. Elimination from $`A`$ to $`R_0`$ **changes** $`C(A)`$ and $`N(A^T)`$, but their dimensions don't change.

### 6.2 Orthogonality of the Subspaces

We are going to connect "**rank**" and "**dimension**":
- **Rank** of a matrix counts independent columns.
- **Dimension** of a subspace is the number of vectors in a basis.

The rank of $`A \in \mathbb{R}^{m \times n}`$ reveals the dimension of all four fundamental subspaces:

1. **Row space**, $`C(A^T)`$, is a subspace of $`\mathbb{R}^n`$, dimension $`r`$. (Each row vector $`\in \mathbb{R}^n`$.)
2. **Column space**, $`C(A)`$, is a subspace of $`\mathbb{R}^m`$, dimension $`r`$. (Each column vector $`\in \mathbb{R}^m`$.)
3. **Nullspace**, $`N(A)`$, is a subspace of $`\mathbb{R}^n`$, dimension $`n - r`$. ($`A\mathbf{x} = \mathbf{0}`$, $`\mathbf{x} \in \mathbb{R}^n`$.)
4. **Left nullspace**, $`N(A^T)`$, is a subspace of $`\mathbb{R}^m`$, dimension $`m - r`$. ($`A^T\mathbf{y} = \mathbf{0}`$, $`\mathbf{y} \in \mathbb{R}^m`$.)

**Orthogonal pairs:**

```math
C(A) \subset \mathbb{R}^m, \quad C(A^T) \subset \mathbb{R}^n
```
```math
N(A^T) \subset \mathbb{R}^m, \quad N(A) \subset \mathbb{R}^n
```

- $`C(A)`$ and $`C(A^T)`$ have $`r`$ dimensions.
- Column space of $`A`$ and row space of $`A`$ have the same dimension $`r`$.
- $`N(A)`$ has dimension $`n - r`$.
- $`N(A^T)`$ has dimension $`m - r`$.

**$`N(A)`$ is perpendicular to $`C(A^T)`$ (row space of $`A`$):**

e.g., $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$

The row vectors $`\lbrace(1, 2), (3, 4)\rbrace`$ span the row space of $`A`$. The solutions to $`\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$ are in the nullspace of $`A`$, $`N(A)`$.

$`N(A) \ni \mathbf{x}`$ is perpendicular to any vector in the row space of $`A`$ in the sense of inner product.

That is, take $`\mathbf{y} \in C(A^T)`$, $`\mathbf{x} \in N(A)`$: $`\mathbf{y} \cdot \mathbf{x} = 0`$.

**$`N(A^T)`$ is perpendicular to $`C(A)`$ (column space of $`A`$):**

e.g., $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$

The column vectors $`\lbrace\begin{pmatrix} 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ 4 \end{pmatrix}\rbrace`$ span the column space of $`A`$.

$`N(A^T) \ni \mathbf{y}`$ implies that $`A^T\mathbf{y} = \mathbf{0}`$:

$`\begin{pmatrix} 1 & 3 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`\Rightarrow C(A) \ni \mathbf{x} \perp \mathbf{y} \in N(A^T)`$

i.e., $`\alpha\begin{pmatrix} 1 \\ 3 \end{pmatrix} + \beta\begin{pmatrix} 2 \\ 4 \end{pmatrix} \perp \mathbf{y}`$

### 6.3 The Four Subspaces for R_0

**Suppose $`R_0 = \text{rref}(A)`$.** The four dimensions are the same for $`R_0`$ and $`A`$.

**Example:** Consider a $`3 \times 5`$ matrix $`R_0`$:

```math
R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}
```

Pivot rows: 1 and 2. Pivot columns: 1 and 4. $`\text{rank}(R_0) = r = 2`$.

**Row space:** Spanned by basis vectors $`\lbrace(1, 3, 5, 0, 7), \;(0, 0, 0, 1, 2)\rbrace`$. $`\dim C(R_0^T) = 2 = r`$.

**Column space:** The 1st and 4th column vectors form $`\left\lbrace\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\rbrace`$, a basis for $`C(R_0)`$. $`\dim C(R_0) = 2 = r`$.

**Nullspace:** $`R_0 \mathbf{x} = \mathbf{0}`$:

```math
\begin{pmatrix} x_1 + 3x_2 + 5x_3 + 7x_5 \\ x_4 + 2x_5 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`n - r = 5 - 2 = 3`$, 3 free variables.

i) Take $`x_2 = 1, x_3 = 0, x_5 = 0`$: $`x_1 + 3 = 0 \Rightarrow x_1 = -3`$, $`x_4 = 0`$. $`\mathbf{x} = \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}`$

ii) Take $`x_2 = 0, x_3 = 1, x_5 = 0`$: $`x_1 + 5 = 0 \Rightarrow x_1 = -5`$, $`x_4 = 0`$. $`\mathbf{x} = \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}`$

iii) Take $`x_2 = 0, x_3 = 0, x_5 = 1`$: $`x_1 + 7 = 0 \Rightarrow x_1 = -7`$, $`x_4 + 2 = 0 \Rightarrow x_4 = -2`$. $`\mathbf{x} = \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}`$

The 3 special solutions form a basis $`\left\lbrace\begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}\right\rbrace`$ for $`N(R_0)`$.

$`\dim N(R_0) = 3 = 5 - 2 = n - r`$.

**Left nullspace:** $`R_0^T \mathbf{y} = \mathbf{0}`$:

```math
\begin{pmatrix} 1 & 0 & 0 \\ 3 & 0 & 0 \\ 5 & 0 & 0 \\ 0 & 1 & 0 \\ 7 & 2 & 0 \end{pmatrix} \begin{pmatrix} y_1 \\ y_2 \\ y_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix} \iff \begin{pmatrix} y_1 \\ 3y_1 \\ 5y_1 \\ y_2 \\ 7y_1 + 2y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}
```

$`m - r = 3 - 2 = 1`$, 1 free variable.

Take $`y_3 = 1`$: $`\Rightarrow y_1 = y_2 = 0`$. $`\therefore \mathbf{y} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$

$`\left\lbrace\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\right\rbrace`$ forms a basis for $`N(R_0^T)`$.

$`\dim N(R_0^T) = 1 = m - r = 3 - 2`$.

$`R_0^T \mathbf{y} = \mathbf{0} \iff \mathbf{y}^T R_0 = \mathbf{0}`$. $`\mathbf{y}^T`$ is a row vector to the 'left' of $`R_0`$.

**Orthogonality summary:**

```math
C(A) \subset \mathbb{R}^m \perp N(A^T) \subset \mathbb{R}^m
```
```math
C(A^T) \subset \mathbb{R}^n \perp N(A) \subset \mathbb{R}^n
```

In $`\mathbb{R}^n`$: the row space and the null space have $`r`$ and $`n - r`$ dimensions.

In $`\mathbb{R}^m`$: the column space and the left nullspace have $`r`$ and $`m - r`$ dimensions.

### 6.4 Relationship Between A and R_0

**The four subspace dimensions for $`A`$ are the same as for $`R_0`$.**

```math
A = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 1 & 3 & 5 & 1 & 9 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}
```

$`A`$ has the same row space as $`R_0`$, but its column space differs from that of $`R_0`$.

A basis for $`C(A)`$: $`\left\lbrace\begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}\right\rbrace`$ (pivot columns of $`A`$, **not** $`R_0`$)

A basis for $`C(R_0)`$: $`\left\lbrace\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\rbrace`$

(Different planes in $`\mathbb{R}^3`$, but same dimension.)

**Key relationships:**

1. $`A`$ has the **same row space** as $`R_0`$, $`R`$: $`C(A^T) = C(R_0^T) = C(R^T)`$, dimension $`r`$.

2. The column space of $`A`$ has dimension $`r`$: $`\dim C(A) = \dim C(A^T)`$. $`C(A) \neq C(R_0)`$, but $`\dim C(A) = \dim C(R_0) = r`$.

3. $`A`$ has the **same nullspace** as $`R_0`$: $`A\mathbf{x} = \mathbf{0} \iff R_0 \mathbf{x} = \mathbf{0}`$. **Elimination does NOT change the solution.** The special solutions form a basis. $`n - r`$ free variables $`\implies`$ $`\dim N(A) = n - r`$.

```math
\dim C(A) + \dim N(A) = r + (n - r) = n
```

4. The left nullspace of $`A`$, $`N(A^T)`$: $`\dim C(A^T) + \dim N(A^T) = r + (m - r) = m`$.

### 6.5 The Fundamental Theorem of Linear Algebra

```math
\boxed{\textbf{Fundamental Theorem of Linear Algebra}}
```

> i) $`\dim C(A) = \dim C(A^T) = r`$
>
> ii) $`\dim N(A) = n - r`$, $`\quad \dim N(A^T) = m - r`$

**The Four Subspaces for $`A`$:**

```math
C(A^T) \subset \mathbb{R}^n \quad \perp \quad N(A) \subset \mathbb{R}^n
```
```math
C(A) \subset \mathbb{R}^m \quad \perp \quad N(A^T) \subset \mathbb{R}^m
```

In $`\mathbb{R}^n`$: $`C(A^T)`$ (dim $`r`$) and $`N(A)`$ (dim $`n - r`$) are orthogonal complements.

In $`\mathbb{R}^m`$: $`C(A)`$ (dim $`r`$) and $`N(A^T)`$ (dim $`m - r`$) are orthogonal complements.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Vector Space | A set closed under addition and scalar multiplication, satisfying 8 axioms |
| Field | A set ($`\mathbb{R}`$, $`\mathbb{C}`$) where $`+, -, \times, \div`$ are defined |
| Subspace | A subset of a vector space that is itself a vector space (must contain $`\mathbf{0}`$) |
| Column Space $`C(A)`$ | All linear combinations of columns of $`A`$; $`A\mathbf{x} = \mathbf{b}`$ solvable iff $`\mathbf{b} \in C(A)`$ |
| Row Space $`C(A^T)`$ | Column space of $`A^T`$; spanned by rows of $`A`$ |
| Nullspace $`N(A)`$ | All solutions to $`A\mathbf{x} = \mathbf{0}`$; subspace of $`\mathbb{R}^n`$ |
| Left Nullspace $`N(A^T)`$ | All $`\mathbf{y}`$ with $`A^T\mathbf{y} = \mathbf{0}`$; subspace of $`\mathbb{R}^m`$ |
| RREF $`R_0 = \text{rref}(A)`$ | Reduced row echelon form; contains $`I`$ in pivot columns, $`F`$ in free columns |
| $`A = CR`$ | $`C`$ = independent columns of $`A`$; $`R = (I \; F)`$ = reduced row echelon form (no zero rows) |
| Special Solutions | Columns of $`\begin{pmatrix} -F \\ I \end{pmatrix}`$; form a basis for $`N(A)`$ |
| Complete Solution | $`\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n`$ (particular + nullspace) |
| Particular Solution $`\mathbf{x}_p`$ | Set all free variables to 0, solve $`R_0\mathbf{x} = \mathbf{d}`$ |
| Full Column Rank ($`r = n`$) | $`N(A) = \lbrace\mathbf{0}\rbrace`$; at most 1 solution to $`A\mathbf{x} = \mathbf{b}`$; $`R_0 = \begin{pmatrix} I \\ 0 \end{pmatrix}`$ |
| Full Row Rank ($`r = m`$) | $`A\mathbf{x} = \mathbf{b}`$ always solvable; $`C(A) = \mathbb{R}^m`$; $`R = (I \; F)`$ |
| Linear Independence | $`c_1\mathbf{v}_1 + \cdots + c_n\mathbf{v}_n = \mathbf{0}`$ only when all $`c_i = 0`$ |
| Spanning | Vectors span $`\mathbb{S}`$ if $`\mathbb{S}`$ = all combinations of those vectors |
| Basis | Linearly independent vectors that span the space; representation is unique |
| Dimension | Number of vectors in any basis; invariant across all bases |
| $`\dim C(A) = \dim C(A^T) = r`$ | Column and row space share the same dimension (rank) |
| $`\dim N(A) = n - r`$ | Nullspace dimension equals number of free variables |
| $`\dim N(A^T) = m - r`$ | Left nullspace dimension |
| $`r + (n - r) = n`$ | Row space + nullspace fill $`\mathbb{R}^n`$ |
| $`r + (m - r) = m`$ | Column space + left nullspace fill $`\mathbb{R}^m`$ |
| Orthogonality | $`N(A) \perp C(A^T)`$ in $`\mathbb{R}^n`$; $`N(A^T) \perp C(A)`$ in $`\mathbb{R}^m`$ |
| Matrix Space dim | $`n \times n`$: $`n^2`$; diagonal: $`n`$; upper triangular: $`\frac{n^2+n}{2}`$; symmetric: $`\frac{n^2+n}{2}`$ |
| Function Space | Solutions to $`y'' = 0, -y, y`$ form 2-dim spaces with bases $`\lbrace1,x\rbrace`$, $`\lbrace\sin x, \cos x\rbrace`$, $`\lbracee^x, e^{-x}\rbrace`$ |
| Fundamental Theorem | $`\dim C(A) = \dim C(A^T) = r`$; $`\dim N(A) = n-r`$; $`\dim N(A^T) = m-r`$ |

---
