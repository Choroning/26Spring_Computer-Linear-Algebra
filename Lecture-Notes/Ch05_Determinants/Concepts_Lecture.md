# Chapter 5 — Determinants

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-18
>
> Introduction to Linear Algebra, Strang, Ch 5

> **Learning Objectives:**
> 1. Compute determinants using cofactor expansion and properties
> 2. Apply the three key properties of determinants
> 3. Use Cramer's Rule to solve linear systems

---

<br>

## Table of Contents

- [1. Determinants (Overview)](#1-determinants-overview)
  - [1.1 Chapter Structure](#11-chapter-structure)
  - [1.2 Key Idea](#12-key-idea)
  - [1.3 Useful Formulas (Preview)](#13-useful-formulas-preview)
- [2. 3x3 Determinants and Cofactors](#2-3x3-determinants-and-cofactors)
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
- [3. Computing and Using Determinants](#3-computing-and-using-determinants)
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
  - [3.11 det(A+B) is not det(A) + det(B)](#311-detab-is-not-deta--detb)
  - [3.12 Linearity in a Single Row](#312-linearity-in-a-single-row)
  - [3.13 Proof of det(A) = det(A^T)](#313-proof-of-deta--detat)
  - [3.14 Cramer's Rule](#314-cramers-rule)
- [4. Areas and Volumes by Determinants](#4-areas-and-volumes-by-determinants)
  - [4.1 Parallelogram in 2D](#41-parallelogram-in-2d)
  - [4.2 Area of a Parallelogram](#42-area-of-a-parallelogram)
  - [4.3 Tilted Box in 3D](#43-tilted-box-in-3d)
  - [4.4 Area of a Triangle](#44-area-of-a-triangle)
  - [4.5 Cross Product](#45-cross-product)
  - [4.6 Geometric Proof of Parallelogram Area](#46-geometric-proof-of-parallelogram-area)
  - [4.7 Volume of a Box](#47-volume-of-a-box)
  - [4.8 Examples: Unit Cube and Axis-Aligned Box](#48-examples-unit-cube-and-axis-aligned-box)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Determinants (Overview)

### 1.1 Chapter Structure

- Textbook Section 5.1, $`3 \times 3`$ determinants and cofactors, is covered in Section 2.
- Textbook Section 5.2, computing and using determinants, is covered in Section 3.
- Textbook Section 5.3, areas and volumes by determinants, is covered in Section 4 and connects geometry to linear algebra.

### 1.2 Key Idea

Let $`A`$ be a square matrix.

- If $`\det(A) = 0`$, then $`A`$ is **not invertible**, and the converse also holds.
  - $`\Longleftrightarrow`$ $`A`$ is **singular**.
  - $`\Longleftrightarrow`$ $`A`$ has **dependent columns**.

### 1.3 Useful Formulas (Preview)

1. $`PA = LU`$

```math
\det(PA) = \det(P)\det(A) = \pm\det(A)
```

```math
\det(LU) = \underbrace{\det(L)}_{=1}\det(U) = \text{product of pivots in } U
```

2. Inverse formula via cofactors:

```math
A^{-1} = \frac{1}{\det(A)}\,C^T, \qquad (A^{-1})_{ij} = \frac{C_{ji}}{\det(A)} \quad (C_{ji} \text{ is the cofactor of } A_{ji})
```

---

<br>

## 2. 3x3 Determinants and Cofactors

### 2.1 Definition of Determinant

A **determinant** is a scalar value associated with a square matrix, denoted as $`\det(A)`$ or $`|A|`$.

### 2.2 2x2 Determinants

**Formula:** The determinant of $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$ is:

```math
\det(A) = ad - bc
```

**Examples:**

```math
\det\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 1, \qquad \det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc
```

```math
\det\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = -1, \qquad \det\begin{pmatrix} c & d \\ a & b \end{pmatrix} = bc - ad
```

```math
\det\begin{pmatrix} b & a \\ d & c \end{pmatrix} = bc - ad
```

**Observation:** Determinants reverse sign when two rows are exchanged (also when two columns are exchanged).

**Singular case:**

```math
\det\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = 0
```

```math
\det\begin{pmatrix} a & b \\ 2a & 2b \end{pmatrix} = 2ab - 2ab = 0
```

If $`b \neq 0`$, the lines $`ax + by = e`$ and $`2ax + 2by = f`$ both have slope $`-a/b`$ (if $`b = 0`$, both are vertical). So the two lines are parallel ($`f \neq 2e`$) or coincide ($`f = 2e`$).

$`\det A = 0`$ means that the columns of $`A`$ are **not** linearly independent, i.e., $`N(A) \neq \lbrace\mathbf{0}\rbrace`$.

**Singular matrix example:**

```math
\text{The singular matrix } \begin{pmatrix} a & 2a \\ c & 2c \end{pmatrix} \text{ has } \det = 0
```

### 2.3 Row Exchange Reverses Sign

**Property:** Row (column) exchange reverses the sign of the determinant.

```math
PA = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \begin{pmatrix} c & d \\ a & b \end{pmatrix}
```

```math
\det(PA) = bc - ad = -\det(A)
```

### 2.4 Linearity of Determinant in Each Row

The determinant of $`\begin{pmatrix} xa + ye & xb + yf \\ c & d \end{pmatrix}`$ is:

```math
d(xa + ye) - c(xb + yf) = x(ad - bc) + y(de - cf)
```

This shows linearity in the first row.

> **Note:** The $`3! = 6`$ terms of a $`3 \times 3`$ determinant fall into 3 groups by row 1 entry, and each group is a row 1 entry times its cofactor (Section 2.7). Transposing the matrix of all cofactors and dividing by $`\det(A)`$ produces $`A^{-1}`$ (Section 2.9).

> **Key Point:** The three key properties of the determinant are $`\det(I) = 1`$, sign reversal under a row exchange (Section 2.3), and linearity in each row (this section and Section 3.12). All other properties of the determinant follow from these three.

### 2.5 3x3 Determinants

For $`3 \times 3`$ matrices, the determinant has $`3! = 6`$ terms. Each term is a product with one entry from each row (and each column); there are 3 choices in row 1, 2 in row 2, and 1 in row 3, producing $`3 \times 2 \times 1 = 6`$ terms.

**Permutation matrices and their determinants:**

Starting from $`I_{3 \times 3}`$:

```math
\det\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = 1 \xrightarrow{C_1 \leftrightarrow C_2} \det\begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix} = -1
```

```math
\xrightarrow{R_2 \leftrightarrow R_3} \det\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} = 1 \xrightarrow{C_2 \leftrightarrow C_3} \det\begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix} = -1
```

```math
\xrightarrow{R_2 \leftrightarrow R_3} \det\begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = 1 \xrightarrow{R_1 \leftrightarrow R_2} \det\begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix} = -1
```

**Key property:** Row (column) exchange multiplies $`\det(A)`$ by $`-1`$.

### 2.6 The Six Terms of a 3x3 Determinant (Permutation Approach)

For $`A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}`$, each of the 6 permutation matrices picks out one term:

| Permutation Matrix | $`\det(P)`$ | Selected Entries | Term |
|:---|:---:|:---|:---:|
| $`\begin{pmatrix} 1&0&0\\0&1&0\\0&0&1 \end{pmatrix}`$ | $`+1`$ | $`a, q, z`$ | $`+aqz`$ |
| $`\begin{pmatrix} 0&1&0\\1&0&0\\0&0&1 \end{pmatrix}`$ | $`-1`$ | $`b, p, z`$ | $`-bpz`$ |
| $`\begin{pmatrix} 0&1&0\\0&0&1\\1&0&0 \end{pmatrix}`$ | $`+1`$ | $`b, r, x`$ | $`+brx`$ |
| $`\begin{pmatrix} 0&0&1\\0&1&0\\1&0&0 \end{pmatrix}`$ | $`-1`$ | $`c, q, x`$ | $`-cqx`$ |
| $`\begin{pmatrix} 0&0&1\\1&0&0\\0&1&0 \end{pmatrix}`$ | $`+1`$ | $`c, p, y`$ | $`+cpy`$ |
| $`\begin{pmatrix} 1&0&0\\0&0&1\\0&1&0 \end{pmatrix}`$ | $`-1`$ | $`a, r, y`$ | $`-ary`$ |

*Table 1. Entries picked by the six permutation matrices and the resulting determinant terms*

**Diagonal trick (Sarrus' Rule):** Write the matrix with the first two columns repeated to the right, then take products along diagonals (down = positive, up = negative). This trick works only for $`3 \times 3`$ matrices and does not extend to $`4 \times 4`$ or larger.

**Combining the 6 terms:**

```math
\det A = aqz - bpz + brx - cqx + cpy - ary
```

**Collecting by first-row entries $`a, b, c`$:**

```math
\det A = a(qz - ry) - b(pz - rx) + c(py - qx)
```

**Submatrix determinants (used in the $`4 \times 4`$ example of Section 2.11):**

```math
A_1 = (2), \quad \det(A_1) = 2
```

```math
A_2 = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}, \quad \det(A_2) = 4 - 1 = 3
```

```math
A_3 = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}, \quad \det(A_3) = 8 - 2 - 2 = 4
```

### 2.7 Cofactors and a Formula for A Inverse

Recall that for a $`3 \times 3`$ matrix $`A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}`$:

```math
\det(A) = a(qz - ry) - b(pz - rx) + c(py - qx)
```

We observe that $`\det(A)`$ is computed by using $`2 \times 2`$ determinants:

```math
\det\begin{pmatrix} q & r \\ y & z \end{pmatrix}, \quad \det\begin{pmatrix} p & r \\ x & z \end{pmatrix}, \quad \det\begin{pmatrix} p & q \\ x & y \end{pmatrix}
```

So the determinant can be written as:

```math
\det(A) = a \det\begin{pmatrix} q & r \\ y & z \end{pmatrix} + b\left(-\det\begin{pmatrix} p & r \\ x & z \end{pmatrix}\right) + c \det\begin{pmatrix} p & q \\ x & y \end{pmatrix}
```

Here, $`a, b, c`$ are **factors** (entries in row 1), and the signed $`2 \times 2`$ determinants are the **cofactors** for $`a, b, c`$:

| Entry | Cofactor |
|:------|:---------|
| $`A_{11} = a`$ | $`C_{11} = +\det\begin{pmatrix} q & r \\ y & z \end{pmatrix}`$ |
| $`A_{12} = b`$ | $`C_{12} = -\det\begin{pmatrix} p & r \\ x & z \end{pmatrix}`$ |
| $`A_{13} = c`$ | $`C_{13} = +\det\begin{pmatrix} p & q \\ x & y \end{pmatrix}`$ |

*Table 2. Row 1 entries and their cofactors*

### 2.8 Cofactor Formula Along Row i

To find the $`(i,j)`$ cofactor $`C_{ij}`$, **remove row $`i`$ and column $`j`$ from $`A`$** and multiply the determinant of the remaining matrix by the sign $`(-1)^{i+j}`$:

```math
C_{ij} = (-1)^{i+j} \det\bigl(\text{remaining } (n{-}1) \times (n{-}1) \text{ matrix}\bigr)
```

**The cofactor formula along row $`i`$:**

```math
\det A = A_{i1}C_{i1} + A_{i2}C_{i2} + \cdots + A_{in}C_{in}
```

$`C_{ij}`$ is the coefficient left after collecting all the terms in $`\det(A)`$ that contain $`A_{ij}`$ and factoring out $`A_{ij}`$.

**Example: Cofactors of a $`2 \times 2`$ matrix**

Given $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$, the cofactors are:

```math
A_{11} = a \to C_{11} = d, \quad A_{12} = b \to C_{12} = -c
```

```math
A_{21} = c \to C_{21} = -b, \quad A_{22} = d \to C_{22} = a
```

The **cofactor matrix** of $`A`$ is:

```math
C = \begin{pmatrix} d & -c \\ -b & a \end{pmatrix}
```

### 2.9 Formula for A Inverse via Cofactors

Verify for $`2 \times 2`$:

```math
AC^T = \begin{pmatrix} a & b \\ c & d \end{pmatrix}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix} = \begin{pmatrix} ad - bc & 0 \\ 0 & ad - bc \end{pmatrix} = (ad - bc)\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = \det(A) \cdot I
```

Therefore:

```math
A \cdot \frac{C^T}{\det(A)} = I
```

```math
\boxed{A^{-1} = \frac{1}{\det(A)} C^T}
```

**Important remarks:**

- For $`A^{-1}`$ to exist, $`\det(A)`$ **cannot be zero**.
- Every entry in $`A^{-1}`$ is a ratio of two determinants:

```math
(A^{-1})_{ij} = \frac{(C^T)_{ij}}{\det(A)} = \frac{C_{ji}}{\det(A)}
```

- $`C^T = \text{adj}(A)`$ is the **adjugate matrix** of $`A`$.

**General formula:**

```math
AC^T = \begin{pmatrix} \det(A) & & \\ & \ddots & \\ & & \det(A) \end{pmatrix} = \det(A) \cdot I
```

### 2.10 Example: 3x3 Inverse via Cofactors

```math
A = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}, \quad \det(A) = 1
```

**Using row reduction** $`(A | I) \to (I | A^{-1})`$:

```math
\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 1 & 1 & | & 0 & 1 & 0 \\ 0 & 0 & 1 & | & 0 & 0 & 1 \end{pmatrix} \xrightarrow{R_1 - R_2,\; R_2 - R_3} \begin{pmatrix} 1 & 0 & 0 & | & 1 & -1 & 0 \\ 0 & 1 & 0 & | & 0 & 1 & -1 \\ 0 & 0 & 1 & | & 0 & 0 & 1 \end{pmatrix}
```

```math
A^{-1} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}
```

**Using cofactor formula:**

```math
A^{-1} = \frac{1}{\det(A)} C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T
```

Computing all cofactors:

```math
C_{11} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1, \quad C_{12} = -\det\begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix} = 0, \quad C_{13} = \det\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = 0
```

```math
C_{21} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = -1, \quad C_{22} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1, \quad C_{23} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = 0
```

```math
C_{31} = \det\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = 0, \quad C_{32} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = -1, \quad C_{33} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1
```

### 2.11 4x4 Determinants

The determinant of $`A_4 \in \mathbb{R}^{4 \times 4}`$ has $`4! = 4 \cdot 3 \cdot 2 \cdot 1 = 24`$ terms.

**Example:**

```math
A_4 = \begin{pmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}
```

Using cofactor expansion along row 1:

```math
\det(A_4) = 2 \det\begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix} + (-1)(-1)\det\begin{pmatrix} -1 & -1 & 0 \\ 0 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}
```

```math
= 2 \cdot 4 + 1 \cdot (-3) = 5
```

- The value $`2`$ in the first term comes from $`(A_4)_{11}`$.
- In the second term, the first $`-1`$ comes from $`(A_4)_{12}`$ and the second $`-1`$ is the sign $`(-1)^{1+2}`$.

**Sign rule:** A cofactor is multiplied by $`(-1)^{i+j}`$, so the minus sign appears when $`i + j`$ is odd for entry $`A_{ij}`$.

**For $`2 \times 2`$:**

```math
\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \det\begin{pmatrix} a & \\ & d \end{pmatrix} + \det\begin{pmatrix} & b \\ c & \end{pmatrix} = ad - bc
```

### 2.12 Scalar Multiplication and Determinant

**$`2 \times 2`$ case:**

```math
A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}, \quad \det(A) = ad - bc
```

```math
2A = \begin{pmatrix} 2a & 2b \\ 2c & 2d \end{pmatrix}, \quad \det(2A) = (2a)(2d) - (2b)(2c) = 2^2(ad - bc) = 4\det(A)
```

```math
3A = \begin{pmatrix} 3a & 3b \\ 3c & 3d \end{pmatrix}, \quad \det(3A) = 3^2(ad - bc) = 9\det(A)
```

```math
\alpha A = \begin{pmatrix} \alpha a & \alpha b \\ \alpha c & \alpha d \end{pmatrix}, \quad \det(\alpha A) = \alpha^2 \det(A)
```

**$`3 \times 3`$ case:**

```math
A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}, \quad \det(A) = a(qz - ry) - b(pz - rx) + c(py - qx)
```

```math
\det(2A) = 2a(2q \cdot 2z - 2r \cdot 2y) - 2b(2p \cdot 2z - 2r \cdot 2x) + 2c(2p \cdot 2y - 2q \cdot 2x) = 2^3 \det(A)
```

```math
\det(\alpha A) = \alpha^3 \det(A)
```

**General formula for $`A \in \mathbb{R}^{n \times n}`$:**

```math
\boxed{\det(\alpha A) = \alpha^n \det(A)}
```

---

<br>

## 3. Computing and Using Determinants

### 3.1 Useful Facts

1. $`\det(A^T) = \det(A)`$.
2. $`\det(AB) = \det(A)\det(B)`$.
3. Orthogonal matrices $`Q`$ have $`|\det(Q)| = 1`$.
4. Elimination matrices have $`\det(E) = 1`$, so $`\det(EA) = \det(E)\det(A) = \det(A)`$.

   For example, $`E = \begin{pmatrix} 1 & & \\ a & 1 & \\ b & c & 1 \end{pmatrix}`$ has $`|E| = 1`$.

5. **Cramer's Rule** finds $`\mathbf{x} = A^{-1}\mathbf{b}`$ from ratios of determinants, which is a slow way.
6. $`\det(A) = \pm(\text{product of the pivots})`$, where the pivots come from $`U`$ in $`PA = LU`$.
7. **The big formula** for $`\det(A)`$ has $`n!`$ terms from $`n!`$ permutations and is very slow if $`n > 3`$.

#### What the Determinant Tells Us

1. An invertible matrix has $`\det(A) \neq 0`$.
2. A singular matrix has $`\det(A) = 0`$.
3. Let $`\lambda`$ be an eigenvalue and $`\mathbf{x}`$ an eigenvector, so that $`A\mathbf{x} = \lambda \mathbf{x}`$. Then:

```math
\Rightarrow (A - \lambda I)\mathbf{x} = \mathbf{0}
```

Since $`\mathbf{x}`$ is not a zero vector, $`A - \lambda I`$ is singular, so:

```math
\det(A - \lambda I) = 0
```

This gives us an equation for $`\lambda`$.

### 3.2 Determinant of Triangular and Diagonal Matrices

**Upper triangular:**

```math
\det\begin{pmatrix} a & b & c \\ 0 & q & r \\ 0 & 0 & z \end{pmatrix} = aqz
```

**Diagonal:**

```math
\det\begin{pmatrix} a & & \\ & q & \\ & & z \end{pmatrix} = aqz
```

We just multiply the diagonal entries to find the determinant.

### 3.3 det(A^T) = det(A)

If we transpose $`A`$, the determinant formula gives the same result.

**Example:** $`A \in \mathbb{R}^{3 \times 3}`$

```math
A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix} \Rightarrow \det(A) = aqz + brx + cpy - cqx - ary - bpz
```

```math
= a(qz - ry) - b(pz - rx) + c(py - qx)
```

```math
A^T = \begin{pmatrix} a & p & x \\ b & q & y \\ c & r & z \end{pmatrix} \Rightarrow \det(A^T) = aqz + pyc + xbr - xqc - ayr - pbz
```

```math
= a(qz - yr) - b(pz - xr) + c(py - xq)
```

Both expressions are the same, confirming $`\det(A^T) = \det(A)`$.

### 3.4 Elimination Matrices and Determinant

**Example:**

```math
A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix}
```

Apply $`R_2 - 2R_1`$:

```math
\begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 4 & 0 \end{pmatrix} = \underbrace{\begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}}_{E_{21}} \underbrace{\begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix}}_{A}
```

Computing determinants:

```math
\det(A) = 8 - 20 = -12
```

```math
\det(E_{21}A) = -12
```

```math
\det(E_{21}) = 1
```

We observe that $`\det(E_{21}A) = \det(E_{21})\det(A) = \det(A) = -12`$.

**Full elimination to $`U`$:**

```math
A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 4 & 0 \end{pmatrix} \xrightarrow{R_3 - 2R_2} \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 0 & -6 \end{pmatrix} = U
```

```math
E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -2 & 1 \end{pmatrix}
```

```math
\det(E_{21}) = 1, \quad \det(E_{32}) = 1, \quad \det(U) = -12
```

```math
\det(U) = \det(E_{32}E_{21}A) = \det(E_{32})\det(E_{21}A) = \det(E_{32})\det(E_{21})\det(A) = \det(A)
```

**Back-substitution to find $`U^{-1}`$:**

Starting from $`(U|I)`$:

```math
\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 2 & 3 & | & 0 & 1 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}
```

Step 1: $`R_2' = R_2/2`$ using $`S_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1/2 & 0 \\ 0 & 0 & 1 \end{pmatrix}`$

```math
\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}
```

Step 2: $`R_1 - R_2`$ using $`E_{12} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}`$

```math
\begin{pmatrix} 1 & 0 & -1/2 & | & 1 & -1/2 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}
```

Step 3: $`R_3' = -R_3/6`$ using $`S_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -1/6 \end{pmatrix}`$

```math
\begin{pmatrix} 1 & 0 & -1/2 & | & 1 & -1/2 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & 1 & | & 0 & 0 & -1/6 \end{pmatrix}
```

Step 4: $`R_1' = R_1 + \frac{1}{2}R_3`$ and $`R_2' = R_2 - \frac{3}{2}R_3`$ using $`E_{13}`$ and $`E_{23}`$:

```math
E_{13} = \begin{pmatrix} 1 & 0 & 1/2 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{23} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -3/2 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 1 & 0 & 0 & | & 1 & -1/2 & -1/12 \\ 0 & 1 & 0 & | & 0 & 1/2 & +1/4 \\ 0 & 0 & 1 & | & 0 & 0 & -1/6 \end{pmatrix} = (I | U^{-1})
```

Therefore: $`E_{23}E_{13}S_3E_{12}S_2 \cdot U = I`$

```math
U = S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1} \cdot I
```

where $`S_2^{-1} = \begin{pmatrix} 1 & \\ & 2 \\ & & 1 \end{pmatrix}`$, $`S_3^{-1} = \begin{pmatrix} 1 & \\ & 1 \\ & & -6 \end{pmatrix}`$, $`E_{12}^{-1} = \begin{pmatrix} 1 & 1 \\ & 1 \\ & & 1 \end{pmatrix}`$, $`E_{13}^{-1} = \begin{pmatrix} 1 & & -1/2 \\ & 1 & \\ & & 1 \end{pmatrix}`$, $`E_{23}^{-1} = \begin{pmatrix} 1 & & \\ & 1 & 3/2 \\ & & 1 \end{pmatrix}`$

**Computing $`\det(U)`$:**

```math
\det(U) = \det(S_2^{-1})\det(E_{12}^{-1})\det(S_3^{-1})\det(E_{13}^{-1})\det(E_{23}^{-1}) = 2 \cdot 1 \cdot (-6) \cdot 1 \cdot 1 = -12
```

**Key result:**

```math
\det(A) = \det(E_{21}^{-1}E_{32}^{-1} U) = \det(E_{21}^{-1}E_{32}^{-1} \cdot S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1} I)
```

```math
= 2 \cdot (-6) \det(\cdots) = -12
```

### 3.5 Scaling Matrices

```math
I = \begin{pmatrix} 1 & \\ & 1 \end{pmatrix}, \quad \det(I) = 1
```

**Scaling in $`x`$-direction:**

```math
S_x = \begin{pmatrix} 2 & \\ & 1 \end{pmatrix}, \quad \det(S_x) = 2 = 2\det(I)
```

```math
\begin{pmatrix} 2 & \\ & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x \\ y \end{pmatrix}
```

The unit square gets stretched to width 2 in the $`x`$-direction.

**Scaling in $`y`$-direction:**

```math
S_y = \begin{pmatrix} 1 & \\ & 2 \end{pmatrix}, \quad \det(S_y) = 2 = 2\det(I)
```

```math
\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x \\ 2y \end{pmatrix}
```

The unit square gets stretched to height 2 in the $`y`$-direction.

**General scaling matrix:**

```math
S = \begin{pmatrix} a & & \\ & b & \\ & & c \end{pmatrix} \Rightarrow \det(S) = abc = abc \cdot \det(I)
```

### 3.6 Elementary Matrices

When $`A`$ has $`A^{-1}`$, $`A`$ can be decomposed into **elementary matrices** (elimination matrices, permutation matrices, and scaling matrices).

> **Definition:** An **elementary matrix** is a square matrix obtained from performing a single elementary row operation on an identity matrix. Examples include elimination, permutation, and scaling matrices.

### 3.7 Product Rule: det(AB) = det(A) det(B)

```math
\det(AB) = \det(A)\det(B)
```

**Proof.** Let $`A, B \in \mathbb{R}^{n \times n}`$.

**Case i)** $`A`$ or $`B`$ is not invertible.

- Suppose $`B`$ is not invertible. Then there exists a nontrivial vector $`\mathbf{x}`$ such that $`B\mathbf{x} = \mathbf{0}`$.
- Multiplying both sides by $`A`$ gives $`AB\mathbf{x} = \mathbf{0}`$.
- This means $`AB`$ is not invertible.
- If $`B`$ is invertible and $`A`$ is not, there is a nonzero vector $`\mathbf{y}`$ with $`A\mathbf{y} = \mathbf{0}`$. Setting $`\mathbf{x} = B^{-1}\mathbf{y} \neq \mathbf{0}`$ gives $`AB\mathbf{x} = A\mathbf{y} = \mathbf{0}`$, so $`AB`$ is again not invertible.
- Therefore, in either case $`\det(AB) = 0`$ and $`\det(A)\det(B) = 0`$, so the two sides are equal.

**Case ii)** $`A, B`$ are both invertible.

```math
\text{rref}(A) = I \Rightarrow A = E_p E_{p-1} \cdots E_2 E_1 I
```

```math
AB = E_p E_{p-1} \cdots E_2 E_1 B
```

```math
\det(AB) = \det(E_p)\det(E_{p-1} \cdots E_2 E_1 B)
```

```math
= \det(E_p)\det(E_{p-1})\det(E_{p-2} \cdots E_2 E_1 B)
```

```math
= \cdots
```

```math
= \det(E_p)\det(E_{p-1}) \cdots \det(E_2)\det(E_1)\det(B)
```

```math
= \underbrace{\det(E_p E_{p-1} \cdots E_2 E_1)}_{\det(A)} \det(B)
```

```math
= \det(A)\det(B) \qquad \square
```

### 3.8 Orthogonal Matrices

Orthogonal matrices $`Q`$ have determinant $`\pm 1`$.

```math
Q^T Q = I
```

```math
\det(Q^T Q) = \det(I) = 1
```

```math
\det(Q^T)\det(Q) = 1
```

```math
(\det Q)^2 = 1
```

```math
\therefore \det(Q) = \pm 1
```

### 3.9 Determinant via Pivots

Invertible matrices have:

```math
\det(A) = \pm (\text{product of pivots})
```

**Derivation:**

```math
\det(A) = \det(LU) = \underbrace{\det(L)}_{= 1}\det(U) = \begin{vmatrix} * & * & * \\ & * & * \\ & & * \end{vmatrix} = \text{product of pivots}
```

With permutation:

```math
\det(A) = \det(PLU) = \underbrace{\det(P)}_{\pm 1} \underbrace{\det(L)}_{1} \det(U) = \pm \text{product of pivots}
```

### 3.10 Determinant of Inverse

If $`A`$ is invertible, then:

```math
AA^{-1} = I
```

```math
\det(AA^{-1}) = \det(I) = 1
```

```math
\det(A) \cdot \det(A^{-1}) = 1
```

```math
\boxed{\det(A^{-1}) = \frac{1}{\det(A)}}
```

Also note:

```math
\det(AB) = \det(A)\det(B) = \det(B)\det(A) = \det(BA)
```

### 3.11 det(A+B) is not det(A) + det(B)

```math
\det(A + B) \stackrel{?}{=} \det(A) + \det(B) \quad \text{(FALSE in general!)}
```

**Counterexample:**

```math
A = \begin{pmatrix} 5 & -6 \\ 0 & -12 \end{pmatrix}, \quad B = \begin{pmatrix} -3 & 0 \\ 1 & 9 \end{pmatrix}
```

```math
A + B = \begin{pmatrix} 2 & -6 \\ 1 & -3 \end{pmatrix}
```

```math
\det(A) = -60, \quad \det(B) = -27, \quad \det(A+B) = -6 + 6 = 0
```

```math
0 \neq -60 + (-27) = -87
```

### 3.12 Linearity in a Single Row

The determinant **is** linear in a single row (while the other rows are fixed).

**$`2 \times 2`$ case, linearity in row 1:**

```math
A = \begin{pmatrix} a_1 + a_2 & b_1 + b_2 \\ c & d \end{pmatrix}
```

```math
\det(A) = (a_1 + a_2)d - c(b_1 + b_2) = (a_1 d - cb_1) + (a_2 d - cb_2)
```

```math
= \begin{vmatrix} a_1 & b_1 \\ c & d \end{vmatrix} + \begin{vmatrix} a_2 & b_2 \\ c & d \end{vmatrix}
```

**$`2 \times 2`$ case, linearity in row 2 (row 2 entries split):**

```math
A = \begin{pmatrix} a & b \\ c_1 + c_2 & d_1 + d_2 \end{pmatrix}
```

```math
\det(A) = a(d_1 + d_2) - b(c_1 + c_2) = (ad_1 - bc_1) + (ad_2 - bc_2)
```

```math
= \begin{vmatrix} a & b \\ c_1 & d_1 \end{vmatrix} + \begin{vmatrix} a & b \\ c_2 & d_2 \end{vmatrix}
```

**General $`n \times n`$ case:** If row $`i`$ of $`A`$ has entries $`a_{ij} = \alpha_{ij} + \beta_{ij}`$, then:

```math
\det(A) = a_{i1}C_{i1} + \cdots + a_{in}C_{in}
```

```math
= (\alpha_{i1} + \beta_{i1})C_{i1} + \cdots + (\alpha_{in} + \beta_{in})C_{in}
```

```math
= (\alpha_{i1}C_{i1} + \alpha_{i2}C_{i2} + \cdots + \alpha_{in}C_{in}) + (\beta_{i1}C_{i1} + \beta_{i2}C_{i2} + \cdots + \beta_{in}C_{in})
```

```math
= \det(A_\alpha) + \det(A_\beta)
```

where $`A_\alpha`$ has row $`i`$ replaced by $`(\alpha_{i1}, \ldots, \alpha_{in})`$ and $`A_\beta`$ has row $`i`$ replaced by $`(\beta_{i1}, \ldots, \beta_{in})`$.

### 3.13 Proof of det(A) = det(A^T)

```math
\det(A) = \det(A^T)
```

**Proof.** Let $`R = \text{rref}(A)`$ be the **reduced row-echelon form** of $`A`$.

Then we can write $`A = E_p E_{p-1} \cdots E_2 E_1 R`$, where $`E_i`$ are elementary matrices.

Transposing gives $`A^T = R^T E_1^T E_2^T \cdots E_{p-1}^T E_p^T`$.

```math
\det(A^T) = \det(R^T)\det(E_1^T) \cdots \det(E_p^T)
```

**Considering $`R`$:**
- If $`A`$ is invertible, $`R = I`$, so $`\det(R) = 1 = \det(R^T)`$.
- Otherwise, $`R`$ has a row of zeros, so $`\det(R) = 0 = \det(R^T)`$.
- Therefore, in all cases $`\det(R) = \det(R^T)`$.

Since $`\det(E_i^T) = \det(E_i)`$ for any elementary matrix:

```math
\det(A^T) = \det(R)\det(E_1) \cdots \det(E_p) = \det(A) \qquad \square
```

### 3.14 Cramer's Rule

**Cramer's Rule to solve $`A\mathbf{x} = \mathbf{b}`$:**

Let $`A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \mathbf{a}_3)`$ with columns $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3`$.

**Finding $`x_1`$:** Replace the 1st column of $`I`$ with $`\mathbf{x}`$:

```math
M_1 = \begin{pmatrix} x_1 & 0 & 0 \\ x_2 & 1 & 0 \\ x_3 & 0 & 1 \end{pmatrix}
```

Apply $`A`$ to $`M_1`$:

```math
AM_1 = \bigl(A\mathbf{x} \;\; A\begin{pmatrix}0\\1\\0\end{pmatrix} \;\; A\begin{pmatrix}0\\0\\1\end{pmatrix}\bigr) = (\mathbf{b} \;\; \mathbf{a}_2 \;\; \mathbf{a}_3) = B_1
```

Take determinants:

```math
\det(AM_1) = \det(B_1)
```

```math
\det(A)\det(M_1) = \det(B_1)
```

```math
\det(A) \cdot x_1 = \det(B_1)
```

```math
\therefore x_1 = \frac{\det(B_1)}{\det(A)}
```

**Finding $`x_2`$:** Introduce $`M_2`$:

```math
AM_2 = A\begin{pmatrix} 1 & x_1 & 0 \\ 0 & x_2 & 0 \\ 0 & x_3 & 1 \end{pmatrix} = (\mathbf{a}_1 \;\; \mathbf{b} \;\; \mathbf{a}_3) = B_2
```

```math
\therefore x_2 = \frac{\det(B_2)}{\det(A)}
```

**Finding $`x_3`$:** Introduce $`M_3`$:

```math
AM_3 = A\begin{pmatrix} 1 & 0 & x_1 \\ 0 & 1 & x_2 \\ 0 & 0 & x_3 \end{pmatrix} = (\mathbf{a}_1 \;\; \mathbf{a}_2 \;\; \mathbf{b}) = B_3
```

```math
\therefore x_3 = \frac{\det(B_3)}{\det(A)}
```

**General Cramer's Rule:** If $`\det(A) \neq 0`$, then $`A\mathbf{x} = \mathbf{b}`$ is solved by:

```math
\boxed{x_j = \frac{\det(B_j)}{\det(A)}}
```

where $`B_j`$ has the $`j`$-th column of $`A`$ replaced by the vector $`\mathbf{b}`$.

**Example 1:**

```math
3x_1 + 4x_2 = 2, \quad 5x_1 + 6x_2 = 4
```

```math
A = \begin{pmatrix} 3 & 4 \\ 5 & 6 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}
```

```math
\det(A) = 18 - 20 = -2
```

```math
B_1 = \begin{pmatrix} 2 & 4 \\ 4 & 6 \end{pmatrix}, \quad \det(B_1) = 12 - 16 = -4
```

```math
B_2 = \begin{pmatrix} 3 & 2 \\ 5 & 4 \end{pmatrix}, \quad \det(B_2) = 12 - 10 = 2
```

```math
x_1 = \frac{-4}{-2} = 2, \quad x_2 = \frac{2}{-2} = -1
```

**Check:**

```math
\begin{pmatrix} 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 6 - 4 \\ 10 - 6 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix} \checkmark
```

**Example 2: Deriving $`A^{-1}`$ via Cramer's Rule**

```math
A = \begin{pmatrix} a & b \\ c & d \end{pmatrix} \in \mathbb{R}^{2 \times 2}, \quad A^{-1} = (\mathbf{x} \;\; \mathbf{y})
```

```math
AA^{-1} = A(\mathbf{x} \;\; \mathbf{y}) = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
```

**i)** $`A\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

```math
A\begin{pmatrix} x_1 & 0 \\ x_2 & 1 \end{pmatrix} = \bigl(\mathbf{b} \;\; A\begin{pmatrix}0\\1\end{pmatrix}\bigr) = \begin{pmatrix} 1 & b \\ 0 & d \end{pmatrix} = B_1
```

```math
x_1 = \frac{\det(B_1)}{\det(A)} = \frac{d}{ad - bc}
```

```math
A\begin{pmatrix} 1 & x_1 \\ 0 & x_2 \end{pmatrix} = \bigl(A\begin{pmatrix}1\\0\end{pmatrix} \;\; \mathbf{b}\bigr) = \begin{pmatrix} a & 1 \\ c & 0 \end{pmatrix} = B_2
```

```math
x_2 = \frac{\det(B_2)}{\det(A)} = \frac{-c}{ad - bc}
```

In the same way, **ii)** solving $`A\mathbf{y} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$ gives $`y_1 = \frac{-b}{ad - bc}`$ and $`y_2 = \frac{a}{ad - bc}`$. Putting the two columns together recovers the standard inverse formula $`A^{-1} = \frac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}`$.

---

<br>

## 4. Areas and Volumes by Determinants

### 4.1 Parallelogram in 2D

A parallelogram in 2D starts from $`(0, 0)`$ with sides:

```math
\mathbf{e}_1 = (a, b), \quad \mathbf{e}_2 = (c, d)
```

The four vertices are $`(0,0)`$, $`(a,b)`$, $`(c,d)`$, and $`(a+c, b+d)`$.

### 4.2 Area of a Parallelogram

```math
\text{Area of the parallelogram} = |ad - bc|
```

```math
= \left|\,\begin{vmatrix} a & b \\ c & d \end{vmatrix}\,\right| = \left|\det\begin{pmatrix} a & b \\ c & d \end{pmatrix}\right| = \left|\det\begin{pmatrix} a & c \\ b & d \end{pmatrix}\right| = |\det(\mathbf{e}_1 \;\; \mathbf{e}_2)|
```

### 4.3 Tilted Box in 3D

A tilted box in 3D starts with three edges $`\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3`$ out from $`(0, 0, 0)`$.

```math
\text{Volume of a tilted box} = \left|\det\begin{pmatrix} \mathbf{e}_1 & \mathbf{e}_2 & \mathbf{e}_3 \end{pmatrix}\right|
```

### 4.4 Area of a Triangle

```math
\text{Area of a triangle} = \triangle = \frac{1}{2} b \cdot h
```

where $`b = \|\mathbf{a}\|`$ is the base, $`h = \|\mathbf{b}\| \sin\theta`$ is the height, and $`\theta`$ is the angle between $`\mathbf{a}`$ and $`\mathbf{b}`$.

```math
\triangle = \frac{1}{2}\|\mathbf{a}\| \|\mathbf{b}\| \sin\theta
```

**Method i)** Using the cross product:

```math
\triangle = \frac{1}{2}\|\mathbf{a} \times \mathbf{b}\|
```

**Method ii)** For 2D vectors $`\mathbf{a} = (x_1, y_1)`$ and $`\mathbf{b} = (x_2, y_2)`$, use the determinant: $`\triangle = \frac{1}{2}|x_1 y_2 - x_2 y_1|`$ (Section 4.5).

### 4.5 Cross Product

In 3D, the **cross product** is:

```math
\mathbf{a} \times \mathbf{b} = \|\mathbf{a}\| \|\mathbf{b}\| \sin\theta \; \hat{n}
```

where $`\hat{n}`$ is the unit normal vector perpendicular to both $`\mathbf{a}`$ and $`\mathbf{b}`$ (i.e., $`\hat{n} \cdot \mathbf{a} = 0`$ and $`\hat{n} \cdot \mathbf{b} = 0`$), with its direction chosen by the right-hand rule.

**For vectors in 2D (embedded in 3D with $`z = 0`$):**

```math
\mathbf{a} = \begin{pmatrix} x_1 \\ y_1 \\ 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} x_2 \\ y_2 \\ 0 \end{pmatrix}, \quad \hat{n} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
```

Unit vectors: $`\hat{i} = \begin{pmatrix}1\\0\\0\end{pmatrix}`$, $`\hat{j} = \begin{pmatrix}0\\1\\0\end{pmatrix}`$, $`\hat{k} = \begin{pmatrix}0\\0\\1\end{pmatrix}`$

```math
\mathbf{a} \times \mathbf{b} = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ x_1 & y_1 & 0 \\ x_2 & y_2 & 0 \end{vmatrix} = \hat{k}\begin{vmatrix} x_1 & y_1 \\ x_2 & y_2 \end{vmatrix} = \hat{k}(x_1 y_2 - x_2 y_1)
```

```math
\|\mathbf{a} \times \mathbf{b}\| = |x_1 y_2 - x_2 y_1|
```

The choice $`\hat{n} = (0, 0, 1)`$ above holds when $`x_1 y_2 - x_2 y_1 > 0`$; if this value is negative, $`\mathbf{a} \times \mathbf{b}`$ points in the $`-\hat{k}`$ direction.

**For edge vectors:**

Let $`\mathbf{e}_1 = \begin{pmatrix} a \\ b \\ 0 \end{pmatrix}`$, $`\mathbf{e}_2 = \begin{pmatrix} c \\ d \\ 0 \end{pmatrix}`$:

```math
\mathbf{e}_1 \times \mathbf{e}_2 = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ a & b & 0 \\ c & d & 0 \end{vmatrix}
```

```math
\boxed{\|\mathbf{e}_1 \times \mathbf{e}_2\| = |ad - bc| = 2\triangle = \text{area of parallelogram}}
```

### 4.6 Geometric Proof of Parallelogram Area

Consider the parallelogram with vertices $`O = (0,0)`$, $`P = (a,b)`$, $`Q = (c,d)`$, $`R = (a+c, b+d)`$. The figure below shows the case $`a, b, c, d > 0`$ and $`ad > bc`$.

Enclosing it in a rectangle of dimensions $`(a+c) \times (b+d)`$ splits the part of the rectangle outside the parallelogram into six corner regions.

```text
+--------+----------------------....R
|  bc    |  ab/2       .........   /|
|        |    .........           / |
+--------Q....                   /  |
|       /                       /   |
|      /                       /    |
|cd/2 /      ad - bc          /     |
|    /                       / cd/2 |
|   /                       /       |
|  /                   ....P--------+
| /           .........    |        |
|/   .........   ab/2      |  bc    |
O....----------------------+--------+
```

*Figure 1. Parallelogram inside the enclosing rectangle and the six corner regions*

The rectangle area is:

```math
\text{Rectangle area} = (a+c)(b+d) = ab + ad + bc + cd
```

Subtract the corner regions:
- Two triangles at the bottom and top (legs $`a`$ and $`b`$): $`2 \times \frac{1}{2}ab = ab`$
- Two rectangles at the upper left and lower right (width $`c`$, height $`b`$): $`2 \times bc = 2bc`$
- Two triangles at the left and right (legs $`c`$ and $`d`$): $`2 \times \frac{1}{2}cd = cd`$

```math
\text{Parallelogram area} = (ab + ad + bc + cd) - ab - 2bc - cd = ad - bc
```

When $`ad < bc`$ (the two edges swap positions), the same computation gives $`bc - ad`$, so in general $`\text{Area} = |ad - bc|`$.

### 4.7 Volume of a Box

**Volume = base area $`\times`$ height**

```math
\|\mathbf{e}_1 \times \mathbf{e}_2\| = \text{base area (parallelogram)}
```

```math
h = |\mathbf{e}_3 \cdot \hat{n}| = \text{height (projection onto normal)}
```

```math
\text{Volume} = h \cdot \|\mathbf{e}_1 \times \mathbf{e}_2\| = |\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)|
```

**Computing $`\mathbf{e}_1 \times \mathbf{e}_2`$:**

Let $`\mathbf{e}_1 = \begin{pmatrix} a \\ b \\ c \end{pmatrix}`$, $`\mathbf{e}_2 = \begin{pmatrix} p \\ q \\ r \end{pmatrix}`$, $`\mathbf{e}_3 = \begin{pmatrix} x \\ y \\ z \end{pmatrix}`$

```math
\mathbf{e}_1 \times \mathbf{e}_2 = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ a & b & c \\ p & q & r \end{vmatrix} = \hat{i}\begin{vmatrix}b&c\\q&r\end{vmatrix} - \hat{j}\begin{vmatrix}a&c\\p&r\end{vmatrix} + \hat{k}\begin{vmatrix}a&b\\p&q\end{vmatrix} = \begin{pmatrix} br - cq \\ cp - ar \\ aq - bp \end{pmatrix}
```

**Triple scalar product:**

```math
\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2) = x(br - cq) + y(cp - ar) + z(aq - bp)
```

```math
= \begin{vmatrix} a & b & c \\ p & q & r \\ x & y & z \end{vmatrix}
```

Therefore:

```math
\boxed{\text{Volume} = \left|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\right| = \left|\det\begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}\right|}
```

### 4.8 Examples: Unit Cube and Axis-Aligned Box

**Unit cube:** the edges lie along the standard basis vectors.

```math
\text{Volume} = \det\begin{pmatrix} 1 & & \\ & 1 & \\ & & 1 \end{pmatrix} = 1
```

**Axis-aligned box** with side lengths $`a, b, c`$:

```math
\text{Volume} = \det\begin{pmatrix} a & & \\ & b & \\ & & c \end{pmatrix} = abc
```

---

<br>

## Summary

| Concept | Key Summary |
|:--------|:---------|
| Determinant definition | The determinant is a scalar value associated with a square matrix, written $`\det(A)`$ or $`\vert A \vert`$. |
| $`2 \times 2`$ determinant | The determinant of $`\begin{pmatrix}a&b\\c&d\end{pmatrix}`$ is $`ad - bc`$. |
| $`3 \times 3`$ determinant | It consists of 6 terms from the $`3! = 6`$ permutations and equals $`a(qz-ry) - b(pz-rx) + c(py-qx)`$. |
| Row/column exchange | Exchanging two rows or two columns reverses the sign of the determinant ($`\det \to -\det`$). |
| Singular matrix | $`\det(A) = 0`$ holds if and only if $`A`$ is not invertible, which holds if and only if its columns are dependent. |
| Cofactor $`C_{ij}`$ | It is $`(-1)^{i+j}`$ times the determinant of the matrix with row $`i`$ and column $`j`$ removed. |
| Cofactor expansion | The determinant expands along any row $`i`$ as $`\det A = \sum_j A_{ij}C_{ij}`$. |
| Inverse via cofactors | $`A^{-1} = \frac{1}{\det(A)}C^T`$, where $`C^T = \text{adj}(A)`$ is the adjugate matrix. |
| Scalar multiplication | For $`A \in \mathbb{R}^{n \times n}`$, $`\det(\alpha A) = \alpha^n \det(A)`$. |
| Transpose | Transposing does not change the determinant, since $`\det(A^T) = \det(A)`$. |
| Product rule | The determinant of a product is the product of the determinants, $`\det(AB) = \det(A)\det(B)`$. |
| Elimination matrices | Since $`\det(E) = 1`$, $`\det(EA) = \det(A)`$. |
| Orthogonal matrices | An orthogonal matrix $`Q`$ has $`\det(Q) = \pm 1`$. |
| Triangular/diagonal | The determinant of a triangular or diagonal matrix is the product of its diagonal entries. |
| Pivots formula | Using $`PA = LU`$, $`\det(A) = \pm(\text{product of pivots})`$. |
| Inverse determinant | For an invertible $`A`$, $`\det(A^{-1}) = 1/\det(A)`$. |
| $`\det(A+B)`$ | In general, $`\det(A+B)`$ is not equal to $`\det(A) + \det(B)`$. |
| Row linearity | The determinant is linear in each row separately while the other rows are fixed. |
| Elementary matrices | An invertible matrix $`A`$ decomposes into a product of elimination, permutation, and scaling matrices. |
| Cramer's Rule | If $`\det(A) \neq 0`$, then $`x_j = \det(B_j)/\det(A)`$, where $`B_j`$ has column $`j`$ of $`A`$ replaced by $`\mathbf{b}`$. |
| Parallelogram area (2D) | The area is $`\vert ad - bc \vert = \vert \det(\mathbf{e}_1\;\mathbf{e}_2) \vert`$. |
| Triangle area | The area is $`\frac{1}{2}\Vert\mathbf{a} \times \mathbf{b}\Vert`$. |
| Cross product | The cross product is $`\mathbf{a} \times \mathbf{b} = \Vert\mathbf{a}\Vert\Vert\mathbf{b}\Vert\sin\theta\;\hat{n}`$ and is computed with a $`3 \times 3`$ determinant. |
| Box volume (3D) | The volume is the absolute value of the triple scalar product, $`\vert\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\vert = \vert\det(\mathbf{e}_1\;\mathbf{e}_2\;\mathbf{e}_3)\vert`$. |

---

<br>

## Self-Check Questions

1. **Cofactor expansion:** Compute the determinant of $`A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 2 \end{pmatrix}`$ by cofactor expansion along row 1, and decide whether $`A`$ is invertible.

   <details>
   <summary>Answer</summary>

   The row 1 cofactors are $`C_{11} = +\det\begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix} = 5`$, $`C_{12} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix} = -2`$, and $`C_{13} = +\det\begin{pmatrix} 1 & 3 \\ 0 & 1 \end{pmatrix} = 1`$. So $`\det A = 2 \cdot 5 + 1 \cdot (-2) + 0 \cdot 1 = 8`$, and since $`\det A \neq 0`$, $`A`$ is invertible.

   </details>

2. **Properties of determinants:** Let $`A`$ be a $`3 \times 3`$ matrix with $`\det(A) = 2`$. Find $`\det(3A)`$, $`\det(A^{-1})`$, $`\det(A^T A)`$, and the determinant of the matrix obtained by exchanging two rows of $`A`$.

   <details>
   <summary>Answer</summary>

   $`\det(3A) = 3^3 \det(A) = 54`$. $`\det(A^{-1}) = 1/\det(A) = 1/2`$. $`\det(A^T A) = \det(A^T)\det(A) = 2 \cdot 2 = 4`$. Exchanging two rows reverses the sign, so that determinant is $`-2`$.

   </details>

3. **Cramer's Rule:** Solve the system $`2x_1 + x_2 = 5`$, $`x_1 + 3x_2 = 5`$ by Cramer's Rule.

   <details>
   <summary>Answer</summary>

   With $`A = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix}`$, $`\det A = 6 - 1 = 5`$. $`B_1 = \begin{pmatrix} 5 & 1 \\ 5 & 3 \end{pmatrix}`$ has determinant $`15 - 5 = 10`$, and $`B_2 = \begin{pmatrix} 2 & 5 \\ 1 & 5 \end{pmatrix}`$ has determinant $`10 - 5 = 5`$. So $`x_1 = 10/5 = 2`$ and $`x_2 = 5/5 = 1`$, which checks out: $`2 \cdot 2 + 1 = 5`$ and $`2 + 3 \cdot 1 = 5`$.

   </details>

4. **Determinant via pivots:** Elimination with one row exchange gives $`PA = LU`$, and the diagonal entries of $`U`$ are $`1, 2, -6`$. Find $`\det(A)`$.

   <details>
   <summary>Answer</summary>

   One row exchange gives $`\det(P) = -1`$, and $`\det(L) = 1`$, $`\det(U) = 1 \cdot 2 \cdot (-6) = -12`$. From $`\det(P)\det(A) = \det(L)\det(U)`$, $`-\det(A) = -12`$, so $`\det(A) = 12`$.

   </details>

5. **Area and volume:** Find the area of the parallelogram with edge vectors $`\mathbf{e}_1 = (3, 1)`$ and $`\mathbf{e}_2 = (1, 2)`$, and the area of the triangle formed by the two vectors. Also find the volume of the box with edges $`(1, 0, 0)`$, $`(1, 2, 0)`$, and $`(1, 1, 3)`$.

   <details>
   <summary>Answer</summary>

   The parallelogram area is $`|3 \cdot 2 - 1 \cdot 1| = 5`$, and the triangle area is half of it, $`5/2`$. The box volume is the absolute value of the determinant of the matrix $`\begin{pmatrix} 1 & 0 & 0 \\ 1 & 2 & 0 \\ 1 & 1 & 3 \end{pmatrix}`$ with the edges as rows; this matrix is lower triangular, so the volume is $`1 \cdot 2 \cdot 3 = 6`$.

   </details>
