# Chapter 5 Lecture — Determinants

> **Last Updated:** 2026-03-30

---

<br>

## Table of Contents

- [1. Determinants (Overview)](#1-determinants-overview)
- [2. 3x3 Determinants and Cofactors (N5-1)](#2-3x3-determinants-and-cofactors-n5-1)
  - [2.1 Definition of Determinant](#21-definition-of-determinant)
  - [2.2 2x2 Determinants](#22-2x2-determinants)
  - [2.3 Row Exchange Reverses Sign](#23-row-exchange-reverses-sign)
  - [2.4 Linearity of Determinant in Each Row](#24-linearity-of-determinant-in-each-row)
  - [2.5 3x3 Determinants](#25-3x3-determinants)
  - [2.6 The Six Terms of a 3x3 Determinant (Permutation Approach)](#26-the-six-terms-of-a-3x3-determinant-permutation-approach)
  - [2.7 Cofactors and a Formula for A Inverse](#27-cofactors-and-a-formula-for-a-inverse)
  - [2.8 Cofactor Formula Along Row i](#28-cofactor-formula-along-row-i)
  - [2.9 Formula for A Inverse via Cofactors](#29-formula-for-a-inverse-via-cofactors)
  - [2.10 Example: 3x3 Inverse via Cofactors](#210-example-3x3-inverse-via-cofactors)
  - [2.11 4x4 Determinants](#211-4x4-determinants)
  - [2.12 Scalar Multiplication and Determinant](#212-scalar-multiplication-and-determinant)
- [3. Computing and Using Determinants (N5-2)](#3-computing-and-using-determinants-n5-2)
  - [3.1 Useful Facts](#31-useful-facts)
  - [3.2 Determinant of Triangular and Diagonal Matrices](#32-determinant-of-triangular-and-diagonal-matrices)
  - [3.3 det(A^T) = det(A)](#33-detat--deta)
  - [3.4 Elimination Matrices and Determinant](#34-elimination-matrices-and-determinant)
  - [3.5 Scaling Matrices](#35-scaling-matrices)
  - [3.6 Elementary Matrices](#36-elementary-matrices)
  - [3.7 Product Rule: det(AB) = det(A) det(B)](#37-product-rule-detab--deta-detb)
  - [3.8 Orthogonal Matrices](#38-orthogonal-matrices)
  - [3.9 Determinant via Pivots](#39-determinant-via-pivots)
  - [3.10 Determinant of Inverse](#310-determinant-of-inverse)
  - [3.11 det(A+B) is NOT det(A) + det(B)](#311-detab-is-not-deta--detb)
  - [3.12 Linearity in a Single Row](#312-linearity-in-a-single-row)
  - [3.13 Proof of det(A) = det(A^T)](#313-proof-of-deta--detat)
  - [3.14 Cramer's Rule](#314-cramers-rule)
- [4. Areas and Volumes by Determinants (N5-3)](#4-areas-and-volumes-by-determinants-n5-3)
  - [4.1 Parallelogram in 2D](#41-parallelogram-in-2d)
  - [4.2 Area of a Parallelogram](#42-area-of-a-parallelogram)
  - [4.3 Tilted Box in 3D](#43-tilted-box-in-3d)
  - [4.4 Area of a Triangle](#44-area-of-a-triangle)
  - [4.5 Cross Product](#45-cross-product)
  - [4.6 Geometric Proof of Parallelogram Area](#46-geometric-proof-of-parallelogram-area)
  - [4.7 Volume of a Box](#47-volume-of-a-box)
  - [4.8 Examples: Unit Cube and Axis-Aligned Box](#48-examples-unit-cube-and-axis-aligned-box)
- [Summary](#summary)

---

<br>

## 1. Determinants (Overview)

### Chapter Structure

- **5.1** &mdash; $3 \times 3$ determinants and cofactors
- **5.2** &mdash; Computing and using determinants
- **5.3** &mdash; Areas and volumes by determinants (connect geometry to linear algebra)

### Key Idea

Let $A$ be a square matrix.

- If $\det(A) = 0$, then $A$ is **NOT invertible**
  - $\Longleftrightarrow$ $A$ is **singular**
  - $\Longleftrightarrow$ $A$ has **dependent columns**

### Useful Formulas (Preview)

1. $PA = LU$

$$\det(PA) = \det(P)\det(A) = \pm\det(A)$$

$$\det(LU) = \underbrace{\det(L)}_{=1}\det(U) = \text{product of pivots in } U$$

2. Inverse formula via cofactors:

$$(A^{-1})_{ij} = \frac{1}{\det(A)}\bigl(\text{cofactor of } A_{ij}\bigr)^T$$

---

<br>

## 2. 3x3 Determinants and Cofactors (N5-1)

### 2.1 Definition of Determinant

A **determinant** is a scalar value associated with a square matrix, denoted as $\det(A)$ or $|A|$.

### 2.2 2x2 Determinants

**Formula:** The determinant of $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$ is:

$$\det(A) = ad - bc$$

**Examples:**

$$\det\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 1, \qquad \det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc$$

$$\det\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = -1, \qquad \det\begin{pmatrix} c & d \\ a & b \end{pmatrix} = bc - ad$$

$$\det\begin{pmatrix} b & a \\ d & c \end{pmatrix} = bc - ad$$

**Observation:** Determinants reverse sign when two rows are exchanged (also when two columns are exchanged).

**Singular case:**

$$\det\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = 0$$

$$\det\begin{pmatrix} a & b \\ 2a & 2b \end{pmatrix} = 2ab - 2ab = 0$$

If $ax + by = e$ has slope $-a/b$ and $2ax + 2by = f$ has slope $-a/b$ (same slope), the lines are parallel.

$\det A = 0$ means that the columns of $A$ are **NOT** linearly independent, i.e., $N(A) \neq \{\mathbf{0}\}$.

**Singular matrix example:**

$$\text{The singular matrix } \begin{pmatrix} a & 2a \\ c & 2c \end{pmatrix} \text{ has } \det = 0$$

### 2.3 Row Exchange Reverses Sign

**Property:** Row (column) exchange reverses the sign of the determinant.

$$PA = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \begin{pmatrix} c & d \\ a & b \end{pmatrix}$$

$$\det(PA) = bc - ad = -\det(A)$$

### 2.4 Linearity of Determinant in Each Row

The determinant of $\begin{pmatrix} xa + ye & xb + yf \\ c & d \end{pmatrix}$ is:

$$d(xa + ye) - c(xb + yf) = x(ad - bc) + y(de - cf)$$

This shows linearity in the first row.

**Remark:** $3 \times 3$ determinants have $3! = 6$ terms, which are separated into 3 terms (cofactors). Dividing those cofactors by $\det(A)$ produces $A^{-1}$.

### 2.5 3x3 Determinants

For $3 \times 3$ matrices, the determinant has $3! = 6$ terms, with one entry from each row (and each column), producing $3 \times 2 = 6$ terms.

**Permutation matrices and their determinants:**

Starting from $I_{3 \times 3}$:

$$\det\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = 1 \xrightarrow{C_1 \leftrightarrow C_2} \det\begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix} = -1$$

$$\xrightarrow{R_2 \leftrightarrow R_3} \det\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} = 1 \xrightarrow{C_2 \leftrightarrow C_3} \det\begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix} = -1$$

$$\xrightarrow{R_1 \leftrightarrow R_2} \det\begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = 1 \xrightarrow{R_1 \leftrightarrow R_2} \det\begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix} = -1$$

**Key property:** Row (column) exchange multiplies $\det(A)$ by $-1$.

### 2.6 The Six Terms of a 3x3 Determinant (Permutation Approach)

For $A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}$, each of the 6 permutation matrices picks out one term:

| Permutation Matrix | $\det(P)$ | Selected Entries | Term |
|:---|:---:|:---|:---:|
| $\begin{pmatrix} 1&0&0\\0&1&0\\0&0&1 \end{pmatrix}$ | $+1$ | $a, q, z$ | $+aqz$ |
| $\begin{pmatrix} 0&1&0\\1&0&0\\0&0&1 \end{pmatrix}$ | $-1$ | $b, p, z$ | $-bpz$ |
| $\begin{pmatrix} 0&1&0\\0&0&1\\1&0&0 \end{pmatrix}$ | $+1$ | $b, r, x$ | $+brx$ |
| $\begin{pmatrix} 0&0&1\\0&1&0\\1&0&0 \end{pmatrix}$ | $-1$ | $c, q, x$ | $-cqx$ |
| $\begin{pmatrix} 0&0&1\\1&0&0\\0&1&0 \end{pmatrix}$ | $+1$ | $c, p, y$ | $+cpy$ |
| $\begin{pmatrix} 1&0&0\\0&0&1\\0&1&0 \end{pmatrix}$ | $-1$ | $a, r, y$ | $-ary$ |

**Diagonal trick (Sarrus' Rule):** Write the matrix with the first two columns repeated to the right, then take products along diagonals (down = positive, up = negative).

**Combining the 6 terms:**

$$\det A = aqz - bpz + brx - cqx + cpy - ary$$

**Collecting by first-row entries** $a, b, c$:

$$\det A = a(qz - ry) - b(pz - rx) + c(py - qx)$$

**Submatrix determinants used in the formula:**

$$A_1 = (2), \quad \det(A_1) = 2$$

$$A_2 = \begin{pmatrix} 2 & 7 \\ -1 & 2 \end{pmatrix}, \quad \det(A_2) = 4 - (-7) = 3$$

$$A_3 = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}, \quad \det(A_3) = 8 - 2 - 2 = 4$$

### 2.7 Cofactors and a Formula for A Inverse

Recall that for a $3 \times 3$ matrix $A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}$:

$$\det(A) = a(qz - ry) - b(pz - rx) + c(py - qx)$$

We observe that $\det(A)$ is computed by using $2 \times 2$ determinants:

$$\det\begin{pmatrix} q & r \\ y & z \end{pmatrix}, \quad \det\begin{pmatrix} p & r \\ x & z \end{pmatrix}, \quad \det\begin{pmatrix} p & q \\ x & y \end{pmatrix}$$

So the determinant can be written as:

$$\det(A) = a \det\begin{pmatrix} q & r \\ y & z \end{pmatrix} + b\left(-\det\begin{pmatrix} p & r \\ x & z \end{pmatrix}\right) + c \det\begin{pmatrix} p & q \\ x & y \end{pmatrix}$$

Here, $a, b, c$ are **factors** (entries in row 1), and the $2 \times 2$ determinants are the **cofactors** for $a, b, c$:

| Entry | Cofactor |
|:------|:---------|
| $A_{11} = a$ | $C_{11} = +\det\begin{pmatrix} q & r \\ y & z \end{pmatrix}$ |
| $A_{12} = b$ | $C_{12} = -\det\begin{pmatrix} p & r \\ x & z \end{pmatrix}$ |
| $A_{13} = c$ | $C_{13} = +\det\begin{pmatrix} p & q \\ x & y \end{pmatrix}$ |

### 2.8 Cofactor Formula Along Row i

For the $(i,j)$ cofactor $C_{ij}$: **remove row $i$ and column $j$ from $A$**.

$$C_{ij} = (-1)^{i+j} \det\bigl(\text{remaining } (n{-}1) \times (n{-}1) \text{ matrix}\bigr)$$

**The cofactor formula along row $i$:**

$$\det A = A_{i1}C_{i1} + A_{i2}C_{i2} + \cdots + A_{in}C_{in}$$

$C_{ij}$ is a collection of all the terms in $\det(A)$ that are multiplied by $A_{ij}$.

**Example: Cofactors of a $2 \times 2$ matrix**

Given $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, the cofactors are:

$$A_{11} = a \to C_{11} = d, \quad A_{12} = b \to C_{12} = -c$$

$$A_{21} = c \to C_{21} = -b, \quad A_{22} = d \to C_{22} = a$$

The **cofactor matrix** of $A$ is:

$$C = \begin{pmatrix} d & -c \\ -b & a \end{pmatrix}$$

### 2.9 Formula for A Inverse via Cofactors

Verify for $2 \times 2$:

$$AC^T = \begin{pmatrix} a & b \\ c & d \end{pmatrix}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix} = \begin{pmatrix} ad - bc & 0 \\ 0 & ad - bc \end{pmatrix} = (ad - bc)\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = \det(A) \cdot I$$

Therefore:

$$A \cdot \frac{C^T}{\det(A)} = I$$

$$\boxed{A^{-1} = \frac{1}{\det(A)} C^T}$$

**Important remarks:**

- For $A^{-1}$ to exist, $\det(A)$ **cannot be zero**.
- Every entry in $A^{-1}$ is a ratio of two determinants:

$$(A^{-1})_{ij} = \frac{(C^T)_{ij}}{\det(A)} = \frac{C_{ji}}{\det(A)}$$

- $C^T = \text{adj}(A)$ is the **adjugate matrix** of $A$.

**General formula:**

$$AC^T = \begin{pmatrix} \det(A) & & 0 \\ & \det(A) & \\ 0 & & \ddots & \det(A) \end{pmatrix} = \det(A) \cdot I$$

### 2.10 Example: 3x3 Inverse via Cofactors

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}, \quad \det(A) = 1$$

**Using row reduction** $(A | I) \to (I | A^{-1})$:

$$\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 1 & 1 & | & 0 & 1 & 0 \\ 0 & 0 & 1 & | & 0 & 0 & 1 \end{pmatrix} \xrightarrow{R_1 - R_2,\; R_1 - R_3} \begin{pmatrix} 1 & 0 & 0 & | & 1 & -1 & 0 \\ 0 & 1 & 0 & | & 0 & 1 & -1 \\ 0 & 0 & 1 & | & 0 & 0 & 1 \end{pmatrix}$$

$$A^{-1} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}$$

**Using cofactor formula:**

$$A^{-1} = \frac{1}{\det(A)} C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T$$

Computing all cofactors:

$$C_{11} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1, \quad C_{12} = -\det\begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix} = 0, \quad C_{13} = \det\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = 0$$

$$C_{21} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = -1, \quad C_{22} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1, \quad C_{23} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = 0$$

$$C_{31} = \det\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = 0, \quad C_{32} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = -1, \quad C_{33} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1$$

### 2.11 4x4 Determinants

$A_4 \in \mathbb{R}^{4 \times 4}$ has $4! = 4 \cdot 3 \cdot 2 \cdot 1 = 24$ terms.

**Example:**

$$A_4 = \begin{pmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}$$

Using cofactor expansion along row 1:

$$\det(A_4) = 2 \det\begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix} + (-1)(-1)\det\begin{pmatrix} -1 & -1 & 0 \\ 0 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}$$

$$= 2 \cdot 4 + 1 \cdot (-3) = 5$$

- The value $2$ comes from $(A_4)_{11}$
- The value $-1$ comes from $(A_4)_{12}$

**Sign rule:** Multiply by $(-1)^{i+j}$: the minus sign appears when $i + j$ is an odd number for entry $A_{ij}$.

**For $2 \times 2$:**

$$\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \det\begin{pmatrix} a & \\ & d \end{pmatrix} + \det\begin{pmatrix} & b \\ c & \end{pmatrix} = ad - bc$$

### 2.12 Scalar Multiplication and Determinant

**$2 \times 2$ case:**

$$A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}, \quad \det(A) = ad - bc$$

$$2A = \begin{pmatrix} 2a & 2b \\ 2c & 2d \end{pmatrix}, \quad \det(2A) = (2a)(2d) - (2b)(2c) = 2^2(ad - bc) = 4\det(A)$$

$$3A = \begin{pmatrix} 3a & 3b \\ 3c & 3d \end{pmatrix}, \quad \det(3A) = 3^2(ad - bc) = 9\det(A)$$

$$\alpha A = \begin{pmatrix} \alpha a & \alpha b \\ \alpha c & \alpha d \end{pmatrix}, \quad \det(\alpha A) = \alpha^2 \det(A)$$

**$3 \times 3$ case:**

$$A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}, \quad \det(A) = a(qz - ry) - b(pz - rx) + c(py - qx)$$

$$\det(2A) = 2a(2q \cdot 2z - 2r \cdot 2y) - 2b(2p \cdot 2z - 2r \cdot 2x) + 2c(2p \cdot 2y - 2q \cdot 2x) = 2^3 \det(A)$$

$$\det(\alpha A) = \alpha^3 \det(A)$$

**General formula for $A \in \mathbb{R}^{n \times n}$:**

$$\boxed{\det(\alpha A) = \alpha^n \det(A)}$$

---

<br>

## 3. Computing and Using Determinants (N5-2)

### 3.1 Useful Facts

1. $\det(A^T) = \det(A)$
2. $\det(AB) = \det(A)\det(B)$
3. $|\det(Q)| = 1$ for orthogonal matrices $Q$
4. Elimination matrices have $\det(E) = 1$, so $\det(EA) = \det(E)\det(A) = \det(A)$

   e.g., $E = \begin{pmatrix} 1 & & \\ a & 1 & \\ b & c & 1 \end{pmatrix}$, $|E| = 1$

5. **Cramer's Rule** finds $\mathbf{x} = A^{-1}\mathbf{b}$ from ratios of determinants (a slow way)
6. $\det(A) = \pm$ product of the pivots in $A = LU$
7. **The big formula** for $\det(A)$ has $n!$ terms from $n!$ permutations (very slow if $n > 3$)

### What the determinant tells us

1. An invertible matrix has $\det(A) \neq 0$
2. A singular matrix has $\det(A) = 0$
3. Let $\lambda$ be eigenvalues, $\mathbf{x}$ eigenvectors, such that $A\mathbf{x} = \lambda \mathbf{x}$

$$\Rightarrow (A - \lambda I)\mathbf{x} = \mathbf{0}$$

Since $\mathbf{x}$ is not a zero vector, $A - \lambda I$ is singular, so:

$$\det(A - \lambda I) = 0$$

This gives us an equation for $\lambda$.

### 3.2 Determinant of Triangular and Diagonal Matrices

**Upper triangular:**

$$\det\begin{pmatrix} a & b & c \\ 0 & q & r \\ 0 & 0 & z \end{pmatrix} = aqz$$

**Diagonal:**

$$\det\begin{pmatrix} a & & \\ & q & \\ & & z \end{pmatrix} = aqz$$

We just multiply the diagonal entries to find the determinant.

### 3.3 det(A^T) = det(A)

If we transpose $A$, the determinant formula gives the same result.

**Example:** $A \in \mathbb{R}^{3 \times 3}$

$$A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix} \Rightarrow \det(A) = aqz + brx + cpy - cqx - ary - bpz$$

$$= a(qz - ry) - b(pz - rx) + c(py - qx)$$

$$A^T = \begin{pmatrix} a & p & x \\ b & q & y \\ c & r & z \end{pmatrix} \Rightarrow \det(A^T) = aqz + pyc - pbz + xby - xqc - ayr$$

$$= a(qz - yr) - b(pz - xy) + c(py - xq)$$

Both expressions are the same, confirming $\det(A^T) = \det(A)$.

### 3.4 Elimination Matrices and Determinant

**Example:**

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix}$$

Apply $R_2 - 2R_1$:

$$\begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 4 & 0 \end{pmatrix} = \underbrace{\begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}}_{E_{21}} \underbrace{\begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix}}_{A}$$

Computing determinants:

$$\det(A) = 8 - 20 = -12$$

$$\det(E_{21}A) = -12$$

$$\det(E_{21}) = 1$$

We observe that $\det(E_{21}A) = \det(E_{21})\det(A) = \det(A) = -12$.

**Full elimination to $U$:**

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 4 & 0 \end{pmatrix} \xrightarrow{R_3 - 2R_2} \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 0 & -6 \end{pmatrix} = U$$

$$E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -2 & 1 \end{pmatrix}$$

$$\det(E_{21}) = 1, \quad \det(E_{32}) = 1, \quad \det(U) = -12$$

$$\det(U) = \det(E_{32}E_{21}A) = \det(E_{32})\det(E_{21}A) = \det(E_{32})\det(E_{21})\det(A) = \det(A)$$

**Back-substitution to find $U^{-1}$:**

Starting from $(U|I)$:

$$\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 2 & 3 & | & 0 & 1 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}$$

Step 1: $R_2' = R_2/2$ using $S_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1/2 & 0 \\ 0 & 0 & 1 \end{pmatrix}$

$$\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}$$

Step 2: $R_1 - R_2$ using $E_{12} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$

$$\begin{pmatrix} 1 & 0 & -1/2 & | & 1 & -1/2 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}$$

Step 3: $R_3' = -R_3/6$ using $S_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -1/6 \end{pmatrix}$

$$\begin{pmatrix} 1 & 0 & -1/2 & | & 1 & -1/2 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & 1 & | & 0 & 0 & -1/6 \end{pmatrix}$$

Step 4: $R_1' = R_1 + \frac{1}{2}R_3$ and $R_2' = R_2 - \frac{3}{2}R_3$ using $E_{13}$ and $E_{23}$:

$$E_{13} = \begin{pmatrix} 1 & 0 & 1/2 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{23} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -3/2 \\ 0 & 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 1 & 0 & 0 & | & 1 & -1/2 & -1/12 \\ 0 & 1 & 0 & | & 0 & 1/2 & +1/4 \\ 0 & 0 & 1 & | & 0 & 0 & -1/6 \end{pmatrix} = (I | U^{-1})$$

Therefore: $E_{23}E_{13}S_3E_{12}S_2 \cdot U = I$

$$U = S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1} \cdot I$$

where $S_2^{-1} = \begin{pmatrix} 1 & \\ & 2 \\ & & 1 \end{pmatrix}$, $S_3^{-1} = \begin{pmatrix} 1 & \\ & 1 \\ & & -6 \end{pmatrix}$, $E_{12}^{-1} = \begin{pmatrix} 1 & 1 \\ & 1 \\ & & 1 \end{pmatrix}$, $E_{13}^{-1} = \begin{pmatrix} 1 & & -1/2 \\ & 1 & \\ & & 1 \end{pmatrix}$, $E_{23}^{-1} = \begin{pmatrix} 1 & & \\ & 1 & 3/4 \\ & & 1 \end{pmatrix}$

**Computing $\det(U)$:**

$$\det(U) = \det(S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1}) = 2 \cdot (-6) \det(E_{32}E_{21} I E_{12}^{-1} I E_{13}^{-1}E_{23}^{-1}) = -12$$

**Key result:**

$$\det(A) = \det(E_{32}^{-1}E_{21}^{-1} U) = \det(E_{32}^{-1}E_{21}^{-1} \cdot S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1} I)$$

$$= 2 \cdot (-6) \det(\cdots) = -12$$

### 3.5 Scaling Matrices

$$I = \begin{pmatrix} 1 & \\ & 1 \end{pmatrix}, \quad \det(I) = 1$$

**Scaling in $x$-direction:**

$$S_x = \begin{pmatrix} 2 & \\ & 1 \end{pmatrix}, \quad \det(S_x) = 2 = 2\det(I)$$

$$\begin{pmatrix} 2 & \\ & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x \\ y \end{pmatrix}$$

The unit square gets stretched to width 2 in the $x$-direction.

**Scaling in $y$-direction:**

$$S_y = \begin{pmatrix} 1 & \\ & 2 \end{pmatrix}, \quad \det(S_y) = 2 = 2\det(I)$$

$$\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x \\ 2y \end{pmatrix}$$

The unit square gets stretched to height 2 in the $y$-direction.

**General scaling matrix:**

$$S = \begin{pmatrix} a & & \\ & b & \\ & & c \end{pmatrix} \Rightarrow \det(S) = abc = abc \cdot \det(I)$$

### 3.6 Elementary Matrices

When $A$ has $A^{-1}$, $A$ can be decomposed into **elementary matrices** (elimination matrices, permutation matrices, and scaling matrices).

**Definition:** An **elementary matrix** is a square matrix obtained from performing a single elementary row operation on an identity matrix. Examples include elimination, permutation, and scaling matrices.

### 3.7 Product Rule: det(AB) = det(A) det(B)

$$\det(AB) = \det(A)\det(B)$$

**Proof.** Let $A, B \in \mathbb{R}^{n \times n}$.

**Case i)** $A$ or $B$ is NOT invertible.

- Suppose $B$ is NOT invertible. Then there exists a nontrivial vector $\mathbf{x}$ such that $B\mathbf{x} = \mathbf{0}$.
- Multiply $A$ on both sides: $AB\mathbf{x} = \mathbf{0}$.
- This means $AB$ is NOT invertible.
- Therefore: $\det(AB) = 0$ and $\det(A)\det(B) = 0$.

**Case ii)** $A, B$ are both invertible.

$$\text{rref}(A) = I \Rightarrow A = E_p E_{p-1} \cdots E_2 E_1 I$$

$$AB = E_p E_{p-1} \cdots E_2 E_1 B$$

$$\det(AB) = \det(E_p)\det(E_{p-1} \cdots E_2 E_1 B)$$

$$= \det(E_p)\det(E_{p-1})\det(E_{p-2} \cdots E_2 E_1 B)$$

$$= \cdots$$

$$= \det(E_p)\det(E_{p-1}) \cdots \det(E_2)\det(E_1)\det(B)$$

$$= \underbrace{\det(E_p E_{p-1} \cdots E_2 E_1)}_{\det(A)} \det(B)$$

$$= \det(A)\det(B) \qquad \square$$

### 3.8 Orthogonal Matrices

Orthogonal matrices $Q$ have determinant $\pm 1$.

$$Q^T Q = I$$

$$\det(Q^T Q) = \det(I) = 1$$

$$\det(Q^T)\det(Q) = 1$$

$$(\det Q)^2 = 1$$

$$\therefore \det(Q) = \pm 1$$

### 3.9 Determinant via Pivots

Invertible matrices have:

$$\det(A) = \pm (\text{product of pivots})$$

**Derivation:**

$$\det(A) = \det(LU) = \underbrace{\det(L)}_{= 1}\det(U) = \begin{vmatrix} * & * & * \\ & * & * \\ & & * \end{vmatrix} = \text{product of pivots}$$

With permutation:

$$\det(A) = \det(PLU) = \underbrace{\det(P)}_{\pm 1} \underbrace{\det(L)}_{1} \det(U) = \pm \text{product of pivots}$$

### 3.10 Determinant of Inverse

If $A$ is invertible, then:

$$AA^{-1} = I$$

$$\det(AA^{-1}) = \det(I) = 1$$

$$\det(A) \cdot \det(A^{-1}) = 1$$

$$\boxed{\det(A^{-1}) = \frac{1}{\det(A)}}$$

Also note:

$$\det(AB) = \det(A)\det(B) = \det(B)\det(A) = \det(BA)$$

### 3.11 det(A+B) is NOT det(A) + det(B)

$$\det(A + B) \stackrel{?}{=} \det(A) + \det(B) \quad \text{(FALSE in general!)}$$

**Counterexample:**

$$A = \begin{pmatrix} 5 & -6 \\ 0 & -12 \end{pmatrix}, \quad B = \begin{pmatrix} -3 & 0 \\ 1 & 9 \end{pmatrix}$$

$$A + B = \begin{pmatrix} 2 & -6 \\ 1 & -3 \end{pmatrix}$$

$$\det(A) = -60, \quad \det(B) = -27, \quad \det(A+B) = -6 + 6 = 0$$

$$0 \neq -60 + (-27) = -87$$

### 3.12 Linearity in a Single Row

The determinant **is** linear in a single row (while the other rows are fixed).

**$2 \times 2$ case — linearity in row 1:**

$$A = \begin{pmatrix} a_1 + a_2 & b_1 + b_2 \\ c & d \end{pmatrix}$$

$$\det(A) = (a_1 + a_2)d - c(b_1 + b_2) = (a_1 d - cb_1) + (a_2 d - cb_2)$$

$$= \begin{vmatrix} a_1 & b_1 \\ c & d \end{vmatrix} + \begin{vmatrix} a_2 & b_2 \\ c & d \end{vmatrix}$$

**$2 \times 2$ case — linearity in column 2 (row 2 entries split):**

$$A = \begin{pmatrix} a & b \\ c_1 + c_2 & d_1 + d_2 \end{pmatrix}$$

$$\det(A) = a(d_1 + d_2) - b(c_1 + c_2) = (ad_1 - bc_1) + (ad_2 - bc_2)$$

$$= \begin{vmatrix} a & b \\ c_1 & d_1 \end{vmatrix} + \begin{vmatrix} a & b \\ c_2 & d_2 \end{vmatrix}$$

**General $n \times n$ case:** If row $i$ of $A$ has entries $a_{ij} = \alpha_{ij} + \beta_{ij}$, then:

$$\det(A) = a_{i1}C_{i1} + \cdots + a_{in}C_{in}$$

$$= (\alpha_{i1} + \beta_{i1})C_{i1} + \cdots + (\alpha_{in} + \beta_{in})C_{in}$$

$$= (\alpha_{i1}C_{i1} + \alpha_{i2}C_{i2} + \cdots + \alpha_{in}C_{in}) + (\beta_{i1}C_{i1} + \beta_{i2}C_{i2} + \cdots + \beta_{in}C_{in})$$

$$= \det(A_\alpha) + \det(A_\beta)$$

where $A_\alpha$ has row $i$ replaced by $(\alpha_{i1}, \ldots, \alpha_{in})$ and $A_\beta$ has row $i$ replaced by $(\beta_{i1}, \ldots, \beta_{in})$.

### 3.13 Proof of det(A) = det(A^T)

$$\det(A) = \det(A^T)$$

**Proof.** Let $R_0 = \text{rref}(A)$ ("reduced row-echelon form").

Denote $A = E_p E_{p-1} \cdots E_2 E_1 R$ where $E_i$ are elementary matrices.

Transpose: $A^T = R^T E_1^T E_2^T \cdots E_{p-1}^T E_p^T$

$$\det(A^T) = \det(R^T)\det(E_1^T) \cdots \det(E_p^T)$$

**Considering $R$:**
- If $A$ is invertible, $R = I$, so $\det(R) = 1 = \det(R^T)$
- Otherwise, $R$ has rows of zeros, so $\det(R) = 0 = \det(R^T)$
- Therefore in all cases: $\det(R) = \det(R^T)$

Since $\det(E_i^T) = \det(E_i)$ for any elementary matrix:

$$\det(A^T) = \det(R)\det(E_1) \cdots \det(E_p) = \det(A) \qquad \square$$

### 3.14 Cramer's Rule

**Cramer's Rule to solve $A\mathbf{x} = \mathbf{b}$:**

Let $A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \mathbf{a}_3)$ with columns $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$.

**Finding $x_1$:** Replace the 1st column of $I$ with $\mathbf{x}$:

$$M_1 = \begin{pmatrix} x_1 & 0 & 0 \\ x_2 & 1 & 0 \\ x_3 & 0 & 1 \end{pmatrix}$$

Apply $A$ to $M_1$:

$$AM_1 = \bigl(A\mathbf{x} \;\; A\begin{pmatrix}0\\1\\0\end{pmatrix} \;\; A\begin{pmatrix}0\\0\\1\end{pmatrix}\bigr) = (\mathbf{b} \;\; \mathbf{a}_2 \;\; \mathbf{a}_3) = B_1$$

Take determinants:

$$\det(AM_1) = \det(B_1)$$

$$\det(A)\det(M_1) = \det(B_1)$$

$$\det(A) \cdot x_1 = \det(B_1)$$

$$\therefore x_1 = \frac{\det(B_1)}{\det(A)}$$

**Finding $x_2$:** Introduce $M_2$:

$$AM_2 = A\begin{pmatrix} 1 & x_1 & 0 \\ 0 & x_2 & 0 \\ 0 & x_3 & 1 \end{pmatrix} = (\mathbf{a}_1 \;\; \mathbf{b} \;\; \mathbf{a}_3) = B_2$$

$$\therefore x_2 = \frac{\det(B_2)}{\det(A)}$$

**Finding $x_3$:** Introduce $M_3$:

$$AM_3 = A\begin{pmatrix} 1 & 0 & x_1 \\ 0 & 1 & x_2 \\ 0 & 0 & x_3 \end{pmatrix} = (\mathbf{a}_1 \;\; \mathbf{a}_2 \;\; \mathbf{b}) = B_3$$

$$\therefore x_3 = \frac{\det(B_3)}{\det(A)}$$

**General Cramer's Rule:** If $\det(A) \neq 0$, then $A\mathbf{x} = \mathbf{b}$ is solved by:

$$\boxed{x_j = \frac{\det(B_j)}{\det(A)}}$$

where $B_j$ has the $j$-th column of $A$ replaced by the vector $\mathbf{b}$.

**Example 1:**

$$3x_1 + 4x_2 = 2, \quad 5x_1 + 6x_2 = 4$$

$$A = \begin{pmatrix} 3 & 4 \\ 5 & 6 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}$$

$$\det(A) = 18 - 20 = -2$$

$$B_1 = \begin{pmatrix} 2 & 4 \\ 4 & 6 \end{pmatrix}, \quad \det(B_1) = 12 - 16 = -4$$

$$B_2 = \begin{pmatrix} 3 & 2 \\ 5 & 4 \end{pmatrix}, \quad \det(B_2) = 12 - 10 = 2$$

$$x_1 = \frac{-4}{-2} = 2, \quad x_2 = \frac{2}{-2} = -1$$

**Check:**

$$\begin{pmatrix} 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 6 - 4 \\ 10 - 6 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix} \checkmark$$

**Example 2: Deriving $A^{-1}$ via Cramer's Rule**

$$A = \begin{pmatrix} a & b \\ c & d \end{pmatrix} \in \mathbb{R}^{2 \times 2}, \quad A^{-1} = (\mathbf{x} \;\; \mathbf{y})$$

$$AA^{-1} = A(\mathbf{x} \;\; \mathbf{y}) = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$$

**i)** $A\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

$$A\begin{pmatrix} x_1 & 0 \\ x_2 & 1 \end{pmatrix} = \bigl(\mathbf{b} \;\; A\begin{pmatrix}0\\1\end{pmatrix}\bigr) = \begin{pmatrix} 1 & b \\ 0 & d \end{pmatrix} = B_1$$

$$x_1 = \frac{\det(B_1)}{\det(A)} = \frac{d}{ad - bc}$$

$$A\begin{pmatrix} 1 & x_1 \\ 0 & x_2 \end{pmatrix} = \bigl(A\begin{pmatrix}1\\0\end{pmatrix} \;\; \mathbf{b}\bigr) = \begin{pmatrix} a & 1 \\ c & 0 \end{pmatrix} = B_2$$

$$x_2 = \frac{\det(B_2)}{\det(A)} = \frac{-c}{ad - bc}$$

This recovers the standard inverse formula.

---

<br>

## 4. Areas and Volumes by Determinants (N5-3)

### 4.1 Parallelogram in 2D

A parallelogram in 2D starts from $(0, 0)$ with sides:

$$\mathbf{e}_1 = (a, b), \quad \mathbf{e}_2 = (c, d)$$

The four vertices are $(0,0)$, $(a,b)$, $(c,d)$, and $(a+c, b+d)$.

### 4.2 Area of a Parallelogram

$$\text{Area of the parallelogram} = |ad - bc|$$

$$= \left|\begin{array}{cc} a & b \\ c & d \end{array}\right| = \left|\det\begin{pmatrix} a & b \\ c & d \end{pmatrix}\right| = \left|\det\begin{pmatrix} a & c \\ b & d \end{pmatrix}\right| = |\det(\mathbf{e}_1 \;\; \mathbf{e}_2)|$$

### 4.3 Tilted Box in 3D

A tilted box in 3D starts with three edges $\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3$ out from $(0, 0, 0)$.

$$\text{Volume of a tilted box} = \left|\det\begin{pmatrix} \mathbf{e}_1 & \mathbf{e}_2 & \mathbf{e}_3 \end{pmatrix}\right|$$

### 4.4 Area of a Triangle

$$\text{Area of a triangle} = \triangle = \frac{1}{2} b \cdot h$$

where $b = \|\mathbf{a}\|$ is the base and $h = \|\mathbf{b}\| \sin\theta$ is the height.

$$\triangle = \frac{1}{2}\|\mathbf{a}\| \|\mathbf{b}\| \sin\theta$$

**Method i)** Using the cross product:

$$\triangle = \frac{1}{2}\|\mathbf{a} \times \mathbf{b}\|$$

### 4.5 Cross Product

In 3D, the **cross product** is:

$$\mathbf{a} \times \mathbf{b} = \|\mathbf{a}\| \|\mathbf{b}\| \sin\theta \; \hat{n}$$

where $\hat{n}$ is the unit normal vector perpendicular to both $\mathbf{a}$ and $\mathbf{b}$ (i.e., $\hat{n} \cdot \mathbf{a} = 0$ and $\hat{n} \cdot \mathbf{b} = 0$).

**For vectors in 2D (embedded in 3D with $z = 0$):**

$$\mathbf{a} = \begin{pmatrix} x_1 \\ y_1 \\ 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} x_2 \\ y_2 \\ 0 \end{pmatrix}, \quad \hat{n} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$$

Unit vectors: $\hat{i} = \begin{pmatrix}1\\0\\0\end{pmatrix}$, $\hat{j} = \begin{pmatrix}0\\1\\0\end{pmatrix}$, $\hat{k} = \begin{pmatrix}0\\0\\1\end{pmatrix}$

$$\mathbf{a} \times \mathbf{b} = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ x_1 & y_1 & 0 \\ x_2 & y_2 & 0 \end{vmatrix} = \hat{k}\begin{vmatrix} x_1 & y_1 \\ x_2 & y_2 \end{vmatrix} = \hat{k}(x_1 y_2 - x_2 y_1)$$

$$\|\mathbf{a} \times \mathbf{b}\| = |x_1 y_2 - x_2 y_1|$$

**For edge vectors:**

Let $\mathbf{e}_1 = \begin{pmatrix} a \\ b \\ 0 \end{pmatrix}$, $\mathbf{e}_2 = \begin{pmatrix} c \\ d \\ 0 \end{pmatrix}$:

$$\mathbf{e}_1 \times \mathbf{e}_2 = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ a & b & 0 \\ c & d & 0 \end{vmatrix}$$

$$\boxed{\|\mathbf{e}_1 \times \mathbf{e}_2\| = |ad - bc| = 2\triangle = \text{area of parallelogram}}$$

### 4.6 Geometric Proof of Parallelogram Area

Consider the parallelogram with vertices $(0,0)$, $(a,b)$, $(c,d)$, $(a+c, b+d)$.

Enclose it in a rectangle of dimensions $(a+c) \times (b+d)$:

$$\text{Rectangle area} = (a+c)(b+d) = ab + ad + bc + cd$$

Subtract the corner regions:
- Two yellow triangles: $2 \times \frac{1}{2}ab = ab$
- Two green rectangles: $2 \times bc = 2bc$
- Two pink triangles: $2 \times \frac{1}{2}cd = cd$

$$\text{Parallelogram area} = (ab + ad + bc + cd) - ab - 2bc - cd = ad - bc$$

Taking absolute value: $\text{Area} = |ad - bc|$

### 4.7 Volume of a Box

**Volume = base area $\times$ height**

$$\|\mathbf{e}_1 \times \mathbf{e}_2\| = \text{base area (parallelogram)}$$

$$h = \|\mathbf{e}_3 \cdot \hat{n}\| = \text{height (projection onto normal)}$$

$$\text{Volume} = h \cdot \|\mathbf{e}_1 \times \mathbf{e}_2\| = \|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\|$$

**Computing $\mathbf{e}_1 \times \mathbf{e}_2$:**

Let $\mathbf{e}_1 = \begin{pmatrix} a \\ b \\ c \end{pmatrix}$, $\mathbf{e}_2 = \begin{pmatrix} p \\ q \\ r \end{pmatrix}$, $\mathbf{e}_3 = \begin{pmatrix} x \\ y \\ z \end{pmatrix}$

$$\mathbf{e}_1 \times \mathbf{e}_2 = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ a & b & c \\ p & q & r \end{vmatrix} = \hat{i}\begin{vmatrix}b&c\\q&r\end{vmatrix} - \hat{j}\begin{vmatrix}a&c\\p&r\end{vmatrix} + \hat{k}\begin{vmatrix}a&b\\p&q\end{vmatrix} = \begin{pmatrix} br - cq \\ cp - ar \\ aq - bp \end{pmatrix}$$

**Triple scalar product:**

$$\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2) = x(br - cq) + y(cp - ar) + z(aq - bp)$$

$$= \begin{vmatrix} a & b & c \\ p & q & r \\ x & y & z \end{vmatrix}$$

Therefore:

$$\boxed{\text{Volume} = \left|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\right| = \left|\det\begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}\right|}$$

### 4.8 Examples: Unit Cube and Axis-Aligned Box

**Unit cube:** edges along the standard basis vectors.

$$\text{Volume} = \det\begin{pmatrix} 1 & & \\ & 1 & \\ & & 1 \end{pmatrix} = 1$$

**Axis-aligned box** with side lengths $a, b, c$:

$$\text{Volume} = \det\begin{pmatrix} a & & \\ & b & \\ & & c \end{pmatrix} = abc$$

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Determinant definition | Scalar value associated with a square matrix; $\det(A)$ or $\|A\|$ |
| $2 \times 2$ determinant | $\det\begin{pmatrix}a&b\\c&d\end{pmatrix} = ad - bc$ |
| $3 \times 3$ determinant | 6 terms from $3! = 6$ permutations; $a(qz-ry) - b(pz-rx) + c(py-qx)$ |
| Row/column exchange | Reverses the sign: $\det \to -\det$ |
| Singular matrix | $\det(A) = 0$ iff $A$ is not invertible iff columns are dependent |
| Cofactor $C_{ij}$ | $(-1)^{i+j}\det(\text{matrix with row } i, \text{col } j \text{ removed})$ |
| Cofactor expansion | $\det A = \sum_j A_{ij}C_{ij}$ along any row $i$ |
| Inverse via cofactors | $A^{-1} = \frac{1}{\det(A)}C^T$ where $C^T = \text{adj}(A)$ |
| Scalar multiplication | $\det(\alpha A) = \alpha^n \det(A)$ for $A \in \mathbb{R}^{n \times n}$ |
| Transpose | $\det(A^T) = \det(A)$ |
| Product rule | $\det(AB) = \det(A)\det(B)$ |
| Elimination matrices | $\det(E) = 1$, so $\det(EA) = \det(A)$ |
| Orthogonal matrices | $\det(Q) = \pm 1$ |
| Triangular/diagonal | Determinant = product of diagonal entries |
| Pivots formula | $\det(A) = \pm(\text{product of pivots})$ via $PA = LU$ |
| Inverse determinant | $\det(A^{-1}) = 1/\det(A)$ |
| $\det(A+B)$ | $\neq \det(A) + \det(B)$ in general |
| Row linearity | Determinant is linear in each row separately |
| Elementary matrices | $A$ decomposes into elimination, permutation, scaling matrices |
| Cramer's Rule | $x_j = \det(B_j)/\det(A)$ where $B_j$ has column $j$ replaced by $\mathbf{b}$ |
| Parallelogram area (2D) | $\|ad - bc\| = \|\det(\mathbf{e}_1\;\mathbf{e}_2)\|$ |
| Triangle area | $\frac{1}{2}\|\mathbf{a} \times \mathbf{b}\|$ |
| Cross product | $\mathbf{a} \times \mathbf{b} = \|\mathbf{a}\|\|\mathbf{b}\|\sin\theta\;\hat{n}$; computed via $3 \times 3$ determinant |
| Box volume (3D) | $\|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\| = \|\det(\mathbf{e}_1\;\mathbf{e}_2\;\mathbf{e}_3)\|$ (triple scalar product) |

---
