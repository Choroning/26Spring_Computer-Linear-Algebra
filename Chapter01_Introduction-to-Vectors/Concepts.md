# Chapter 1 Lecture — Vectors and Matrices

> **Last Updated:** 2026-06-06
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 1

> **Prerequisites**: [Linear Algebra] Basic high school algebra (equations, graphs).
>
> **Learning Objectives**:
> 1. Define vectors and perform linear combinations
> 2. Compute dot products and understand geometric meaning of orthogonality
> 3. Perform basic matrix-vector multiplication

> **📑 This document is split into 2 parts.**
>
> **Part 1** · [Part 2](Concepts-2.md)

---

<br>

## Table of Contents

- [0. Course Overview](#0-course-overview)
- [1. Vectors and Matrices (Warmup)](#1-vectors-and-matrices-warmup)
- [1.1 Vectors and Linear Combinations](#11-vectors-and-linear-combinations)
  - [1.1.1 Linear Combinations in R^2](#111-linear-combinations-in-r2)
  - [1.1.2 Two Key Questions](#112-two-key-questions)
  - [1.1.3 Extending to Higher Dimensions](#113-extending-to-higher-dimensions)
  - [1.1.4 Solving Two Equations](#114-solving-two-equations)
  - [1.1.5 Can Elimination Fail?](#115-can-elimination-fail)
  - [1.1.6 Vectors in Three Dimensions](#116-vectors-in-three-dimensions)
- [1.2 Lengths and Angles from Dot Products](#12-lengths-and-angles-from-dot-products)
  - [1.2.1 Dot Product Definition](#121-dot-product-definition)
  - [1.2.2 Length of a Vector](#122-length-of-a-vector)
  - [1.2.3 Unit Vectors](#123-unit-vectors)
  - [1.2.4 Perpendicular Vectors](#124-perpendicular-vectors)
  - [1.2.5 The Angle Between Two Vectors](#125-the-angle-between-two-vectors)
  - [1.2.6 Schwarz Inequality](#126-schwarz-inequality)
  - [1.2.7 Triangle Inequality](#127-triangle-inequality)
  - [1.2.8 Planes in 3D](#128-planes-in-3d)
- [1.3 Matrices and Their Column Spaces](Concepts-2.md#13-matrices-and-their-column-spaces)
  - [1.3.1 Matrix-Vector Multiplication](Concepts-2.md#131-matrix-vector-multiplication)
  - [1.3.2 Column Space](Concepts-2.md#132-column-space)
  - [1.3.3 Independent, Dependent, and Column Space](Concepts-2.md#133-independent-dependent-and-column-space)
  - [1.3.4 Span](Concepts-2.md#134-span)
  - [1.3.5 Rank and Basis](Concepts-2.md#135-rank-and-basis)
  - [1.3.6 Matrices of Rank One](Concepts-2.md#136-matrices-of-rank-one)
- [1.4 Matrix Multiplication AB and CR](Concepts-2.md#14-matrix-multiplication-ab-and-cr)
  - [1.4.1 Rules for Matrix Multiplication](Concepts-2.md#141-rules-for-matrix-multiplication)
  - [1.4.2 Column Interpretation of AB](Concepts-2.md#142-column-interpretation-of-ab)
  - [1.4.3 Computational Cost](Concepts-2.md#143-computational-cost)
  - [1.4.4 Properties of Matrix Multiplication](Concepts-2.md#144-properties-of-matrix-multiplication)
  - [1.4.5 Rank One Matrices and A = CR](Concepts-2.md#145-rank-one-matrices-and-a--cr)
  - [1.4.6 Finding C and R](Concepts-2.md#146-finding-c-and-r)
  - [1.4.7 Columns of A times Rows of B (Outer Product)](Concepts-2.md#147-columns-of-a-times-rows-of-b-outer-product)
- [Summary](Concepts-2.md#summary)

---

<br>

## 0. Course Overview

**Course:** DCSS321 — Computer Linear Algebra (Shinhoo Kang)

**What you will learn in 321:**
- Vector spaces
- Linear Transformations
- Solve System of Linear Equations
- Eigenvalues and Eigenvectors

**Course Plan (Part 1):**

| Week | Course Plan | Key Concepts | Materials |
|:----:|:------------|:-------------|:----------|
| 1 | Introduction to Vectors | Vectors | Chap. 1 |
| 2 | Solving Linear Equations | $`Ax = b`$ | Chap. 2 |
| 3 | Solving Linear Equations | $`Ax = b`$ | Chap. 3 |
| 4 | Four Fundamental Subspaces | Column space | Chap. 3.1--3.2 |
| 5 | Four Fundamental Subspaces | Null space | Chap. 3.3--3.5 |
| 6 | Orthogonality | Orthogonality | Chap. 4.1--4.3 |
| 7 | Orthogonality | Orthogonality | Chap. 4.4--4.5 |
| 8 | **Midterm Exam** | | |

**Course Plan (Part 2):**

| Week | Course Plan | Key Concepts | Materials |
|:----:|:------------|:-------------|:----------|
| 9 | Determinants | Determinants | Chap. 5 |
| 10 | Eigenvalues and Eigenvectors | Eigen decomposition | Chap. 6.1--6.3 |
| 11 | Eigenvalues and Eigenvectors | Eigen decomposition | Chap. 6.4--6.5 |
| 12 | Singular Value Decomposition | SVD | Chap. 7 |
| 13 | Linear Transformation | Linear map | Chap. 8 |
| 14 | Linear Transformation | Linear map | Chap. 8 |
| 15 | Linear Algebra in Optimization | SGD | Chap. 9 |
| 16 | **Final Exam** | | |

---

<br>

## 1. Vectors and Matrices (Warmup)

**Chapter 1 covers:**
- 1.1 Vectors and Linear Combination
- 1.2 Lengths and Angles from Dot Products
- 1.3 Matrices and their Column Spaces
- 1.4 Matrix Multiplication $`AB`$ and $`CR`$

### Key Ideas

- Linear algebra is about **vectors** $`\mathbf{v}`$ and **matrices** $`A`$.
- We define operations: $`+, -, \cdot`$ (addition, subtraction, multiplication).
- Consider vectors and scalars:

```math
\mathbf{v} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

```math
\mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 7 \end{pmatrix} \in \mathbb{R}^2
```

- For scalars $`c, d \in \mathbb{R}`$:

```math
c\mathbf{v} + d\mathbf{w} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix} \in \mathbb{R}^2
```

> The linear combinations fill the $`xy`$ plane.

### Length of a Vector

The length of a vector $`\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} \in \mathbb{R}^2`$:

```math
\|\mathbf{v}\| = \sqrt{v_1^2 + v_2^2}
```

**Example:** $`\mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}`$, $`\|\mathbf{w}\| = \sqrt{1^2 + 3^2} = \sqrt{10}`$

### Dot Product

The dot product of $`\mathbf{v}`$ and $`\mathbf{w}`$ is:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}^T \begin{pmatrix} w_1 \\ w_2 \end{pmatrix} = v_1 w_1 + v_2 w_2
```

**Example:**

```math
\begin{pmatrix} 2 \\ 4 \end{pmatrix} \cdot \begin{pmatrix} 1 \\ 3 \end{pmatrix} = (2)(1) + (4)(3) = 2 + 12 = 14
```

### Matrix from Column Vectors

A matrix $`A`$ contains two columns:

```math
A = \begin{pmatrix} \mathbf{v} & \mathbf{w} \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}
```

### Matrix-Vector Product as Linear Combination

Multiply the matrix $`A`$ by a vector $`\begin{pmatrix} c \\ d \end{pmatrix}`$:

```math
A \begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix} \begin{pmatrix} c \\ d \end{pmatrix} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

This is a **linear combination of $`\mathbf{v}`$ and $`\mathbf{w}`$**.

Let $`\mathbf{x} = \begin{pmatrix} c \\ d \end{pmatrix}`$. All combinations $`A\mathbf{x}`$ produce the **column space** of the matrix $`A`$. (The column space is a plane!)

### Dependent Vectors

Let $`\mathbf{z} = \mathbf{v} + \mathbf{w}`$.

```math
B = \begin{pmatrix} \mathbf{v} & \mathbf{w} & \mathbf{z} \end{pmatrix} = \begin{pmatrix} 2 & 1 & 3 \\ 4 & 3 & 7 \end{pmatrix}
```

- The column space of $`B`$ is still the $`xy`$ plane.
- $`\mathbf{v}`$ and $`\mathbf{w}`$ are **independent** vectors.
- $`\mathbf{z}`$ is a **dependent** vector.

### Matrix Multiplication Preview

A matrix multiplication $`AB`$ can be interpreted as **$`A`$ times each column of $`B`$**.

---

<br>

## 1.1 Vectors and Linear Combinations

### 1.1.1 Linear Combinations in R^2

**(1)** $`2\mathbf{v} - 3\mathbf{w}`$ is a linear combination $`c\mathbf{v} + d\mathbf{w}`$ of the vectors $`\mathbf{v}`$ and $`\mathbf{w}`$.

**(2)** Let $`\mathbf{v} = \begin{pmatrix} 4 \\ 1 \end{pmatrix}`$ and $`\mathbf{w} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, then:

```math
2\mathbf{v} - 3\mathbf{w} = 2\begin{pmatrix} 4 \\ 1 \end{pmatrix} - 3\begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ -5 \end{pmatrix}
```

**(3)** All combinations $`c\begin{pmatrix} 4 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$ fill the $`xy`$ plane.

**(4)** The vectors $`c\begin{pmatrix} 4 \\ 1 \\ 0 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}`$ fill a **plane** in $`xyz`$ space.

$`\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$ is **NOT** on the plane.

### Building Linear Combinations

Linear combinations of vectors $`\mathbf{v}`$ and $`\mathbf{w}`$ are built from two basic operations:
1. **Scalar multiplication:** $`c\mathbf{v}`$, $`d\mathbf{w}`$
2. **Vector addition:** $`\mathbf{v} + \mathbf{w}`$

This leads to a **linear combination** $`c\mathbf{v} + d\mathbf{w}`$ of $`\mathbf{v}`$ and $`\mathbf{w}`$.

### 1.1.2 Two Key Questions

It opens up two questions:

**(1) Describe** all the combinations $`c\mathbf{v} + d\mathbf{w}`$.
- Result is either a **plane** or a **line**.

**(2) Find** the numbers $`c`$ and $`d`$ such that $`c\mathbf{v} + d\mathbf{w} = \mathbf{x}`$.

**Example:** Find $`c`$ and $`d`$ such that:

```math
c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 4 \\ 3 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}
```

### 1.1.3 Extending to Higher Dimensions

So far $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$.

Let's extend the dimension of the vectors and increase the number of vectors:

```math
\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \in \mathbb{R}^m
```

$`n`$ vectors in $`m`$-dimensional space.

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}
```

$`m`$ rows and $`n`$ columns: an $`m`$ by $`n`$ matrix.

This opens up the same questions in higher dimensions:

**(1) Describe** all the combinations:

```math
A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{v}_1 + x_2\mathbf{v}_2 + \cdots + x_n\mathbf{v}_n
```

of the columns of $`A`$.

**(2) Find** the numbers $`x_1`$ to $`x_n`$ such that:

```math
A\mathbf{x} = \mathbf{b}
```

### Linear Combinations $`c\mathbf{v} + d\mathbf{w}`$

Consider $`\mathbf{v} \in \mathbb{R}^2`$, which has two components: $`\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}`$.

Geometrically, the vectors $`\mathbf{v} + \mathbf{w} = \mathbf{w} + \mathbf{v}`$ (commutativity is shown by the parallelogram rule).

The vectors $`c\mathbf{v}`$ fill an infinitely long line in the $`xy`$ plane.

If $`\mathbf{w}`$ is not on that line, then the vectors $`d\mathbf{w}`$ fill the 2nd line.

**The linear combinations $`c\mathbf{v} + d\mathbf{w}`$ fill the plane.**

**Special cases:**
- $`1\mathbf{v} + 1\mathbf{w}`$ = sum of vectors
- $`1\mathbf{v} - 1\mathbf{w}`$ = difference of vectors
- $`0\mathbf{v} + 0\mathbf{w}`$ = zero vector
- $`c\mathbf{v} + 0\mathbf{w}`$ = vector $`c\mathbf{v}`$, in the direction of $`\mathbf{v}`$

### 1.1.4 Solving Two Equations

**Solve:**

```math
c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \quad \cdots (*)
```

This is equivalent to the system:

```math
\begin{cases} 2c + 2d = 8 \\ c - d = 2 \end{cases}
```

**Solution by elimination:**

Divide the first equation by 2 to get $`c + d = 4`$. Then add the two equations:

```math
c + d = 4
```

```math
c - d = 2
```

Adding: $`2c = 6`$, so $`c = 3`$.

Then $`3 + d = 4`$, so $`d = 1`$.

**Verification:** $`(*)`$ becomes:

```math
3\begin{pmatrix} 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \checkmark
```

**In matrix form:**

```math
\begin{pmatrix} 2 & 2 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix}
```

Let $`\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}`$, $`\mathbf{b} = \begin{pmatrix} 8 \\ 2 \end{pmatrix}`$.

```math
c\mathbf{v} + d\mathbf{w} = \mathbf{b}
```

```math
\Updownarrow
```

```math
cv_1 + dw_1 = b_1
```

```math
cv_2 + dw_2 = b_2
```

```math
\Updownarrow
```

```math
\begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
```

**What does this mean geometrically?**

The system $`xv_1 + yw_1 = b_1`$ and $`xv_2 + yw_2 = b_2`$ represents two lines. The solution $`(c, d)`$ is the intersection point:

- Two linear equations are met at the point $`(c, d)`$.
- $`\mathbf{v}`$ and $`\mathbf{w}`$ are **linearly independent**.
- $`A = \begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}`$ is **invertible**.

### 1.1.5 Can Elimination Fail?

**Yes**, when $`\mathbf{v} \| \mathbf{w}`$ (i.e., $`\frac{v_1}{v_2} = \frac{w_1}{w_2}`$).

The system becomes:

```math
xv_1 + yw_1 = b_1
```

```math
xv_2 + yw_2 = b_2
```

Two possibilities:
- **No solution**: All combinations of $`\mathbf{v}`$ and $`\mathbf{w}`$ lie on the same line. If $`\mathbf{b}`$ is NOT on the line, there is no solution. No combination of $`\mathbf{v}`$ and $`\mathbf{w}`$ equals $`\mathbf{b}`$.
- **Infinite solutions**: If $`\mathbf{b}`$ is on the line, there exist infinite solutions.

### 1.1.6 Vectors in Three Dimensions

Suppose $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^3`$ (three-dimensional space).

```math
\mathbf{v} = \begin{pmatrix} 2 \\ 3 \\ 1 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 4 \\ 1 \end{pmatrix}
```

```math
c\mathbf{v} + d\mathbf{w} = \begin{pmatrix} 2c + d \\ 3c + d \\ c \end{pmatrix}
```

> $`c\mathbf{v} + d\mathbf{w}`$ **DO NOT** fill the whole 3D space. Can at most fill a **2D plane**!

**Need three independent vectors for filling 3D space.**

**Example:** $`\hat{i} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}`$, $`\hat{k} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$

```math
c\hat{i} + d\hat{j} + e\hat{k} = \begin{pmatrix} c \\ d \\ e \end{pmatrix}
```

$`\hat{i}, \hat{j}, \hat{k}`$ correspond to the unit vectors along $`x, y, z`$ axes in 3D space.

```math
\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix} = v_1\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + v_2\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + v_3\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
```

**Vector form:** $`= v_1\hat{i} + v_2\hat{j} + v_3\hat{k}`$

**Matrix form:**

```math
= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}
```

> $`\mathbf{v} = I\mathbf{v}`$ (identity matrix)

### How Do We Know If It Is a Plane?

Suppose nonzero vectors $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^3`$ which are **independent** ($`\mathbf{w}`$ is not a multiple $`c\mathbf{v}`$).

Then their linear combinations fill a **plane** inside 3D space (a flat surface).

---

<br>

## 1.2 Lengths and Angles from Dot Products

### 1.2.1 Dot Product Definition

**(1)** Let $`\mathbf{v} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 4 \\ 6 \end{pmatrix}`$.

Dot product of $`\mathbf{v}`$ and $`\mathbf{w}`$ is:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}^T \begin{pmatrix} 4 \\ 6 \end{pmatrix} = 1 \cdot 4 + 2 \cdot 6 = 16
```

**General definition for $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$:**

```math
\mathbf{v} \cdot \mathbf{w} = v_1 w_1 + v_2 w_2
```

**Extension to $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^n`$:**

```math
\mathbf{v} \cdot \mathbf{w} = v_1 w_1 + v_2 w_2 + \cdots + v_n w_n = \sum_{i=1}^{n} v_i w_i
```

The dot product $`\mathbf{v} \cdot \mathbf{v}`$ tells us the **squared length**:

```math
\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + \cdots + v_n^2
```

of a vector $`\mathbf{v}`$.

**Example:** $`\mathbf{v} \in \mathbb{R}^2`$, $`\|\mathbf{v}\|^2 = v_1^2 + v_2^2`$ (Pythagoras formula).

**Example:** $`\mathbf{v} \in \mathbb{R}^3`$, $`\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + v_3^2 = (v_1^2 + v_2^2) + v_3^2`$ where $`(v_1^2 + v_2^2)`$ is in the $`xy`$ plane.

### 1.2.2 Length of a Vector

**(2)** The length squared of $`\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}`$ is:

```math
\|\mathbf{v}\|^2 = \mathbf{v} \cdot \mathbf{v} = 1^2 + 3^2 + 2^2 = 1 + 9 + 4 = 14
```

```math
\therefore \|\mathbf{v}\| = \sqrt{14}
```

### 1.2.3 Unit Vectors

- $`\mathbf{v}`$ is a **unit vector** when $`\|\mathbf{v}\| = 1`$.
- If $`\mathbf{v} \neq \mathbf{0}`$, then $`\frac{\mathbf{v}}{\|\mathbf{v}\|}`$ is a unit vector.

**Example (ex1):** $`\mathbf{u} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}`$ is a unit vector.

```math
\|\mathbf{u}\|^2 = \cos^2\theta + \sin^2\theta = 1
```

### 1.2.4 Perpendicular Vectors

**(3)** $`\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}`$ is perpendicular to $`\mathbf{w} = \begin{pmatrix} 4 \\ -4 \\ 4 \end{pmatrix}`$:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}^T \begin{pmatrix} 4 \\ -4 \\ 4 \end{pmatrix} = (1)(4) + (3)(-4) + (2)(4) = 4 - 12 + 8 = 0
```

Suppose the angle between $`\mathbf{v}`$ and $`\mathbf{w}`$ is $`90°`$:

```math
\cos\theta \Rightarrow \cos 90° = 0
```

```math
\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\| \|\mathbf{w}\| \cos\theta = 0
```

**Pythagorean Theorem for perpendicular vectors:**

```math
\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w})
```

```math
= \mathbf{v} \cdot \mathbf{v} + \mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{v} + \mathbf{w} \cdot \mathbf{w}
```

```math
= \|\mathbf{v}\|^2 + 2\mathbf{v} \cdot \mathbf{w} + \|\mathbf{w}\|^2
```

When $`\mathbf{v} \cdot \mathbf{w} = 0`$:

```math
\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 \quad \text{(Pythagorean theorem)}
```

Note that $`\|\mathbf{v} - \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2`$ also holds when $`\mathbf{v} \perp \mathbf{w}`$.

**Example (ex2):** $`\mathbf{v} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$ is at a $`45°`$ angle with the $`x`$-axis.

$`\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$ is at a $`-45°`$ angle with the $`x`$-axis.

```math
\mathbf{v} + \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} + \begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 2 \\ 0 \end{pmatrix}
```

```math
\mathbf{v} - \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} - \begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0 \\ 2 \end{pmatrix}
```

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}^T \begin{pmatrix} 1 \\ -1 \end{pmatrix} = 1 \cdot 1 + 1 \cdot (-1) = 0
```

```math
\|\mathbf{v}\| = \sqrt{2}, \quad \|\mathbf{w}\| = \sqrt{2}
```

```math
\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = 4
```

```math
\|\mathbf{v} - \mathbf{w}\|^2 = 4
```

**Example (ex3):** $`\mathbf{v} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} -1 \\ 2 \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}^T \begin{pmatrix} -1 \\ 2 \end{pmatrix} = -4 + 4 = 0
```

```math
\Rightarrow \mathbf{v} \perp \mathbf{w}
```

The weights times the distances $`v_1 w_1`$ and $`v_2 w_2`$ are balanced.

**Example (ex4):** $`\mathbf{v} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$ (unit vector), $`\mathbf{w} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \cos\theta
```

The angle between $`\mathbf{v}`$ and $`\mathbf{w}`$ has $`\cos\theta = \mathbf{v} \cdot \mathbf{w}`$, if $`\|\mathbf{v}\| = \|\mathbf{w}\| = 1`$.

### 1.2.5 The Angle Between Two Vectors

**(4)** The angle $`\theta = 45°`$ between $`\mathbf{v} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$ and $`\mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$ has:

```math
\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|} = \frac{1}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}
```

**(5)** All angles have $`|\cos\theta| \leq 1`$.

All vectors have $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$ and $`\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|`$.

**Example (ex5) — Application: Total Cost**

Price vector $`\mathbf{p} = \begin{pmatrix} p_1 \\ p_2 \\ p_3 \end{pmatrix}`$, quantity vector $`\mathbf{q} = \begin{pmatrix} q_1 \\ q_2 \\ q_3 \end{pmatrix}`$

$`p_i q_i \Rightarrow`$ Buying $`q_i`$ units at the price $`p_i`$.

```math
\mathbf{p} \cdot \mathbf{q} = p_1 q_1 + p_2 q_2 + p_3 q_3 \Rightarrow \text{total cost}
```

**The dot product $`\mathbf{v} \cdot \mathbf{w}`$ finds the angle between any two nonzero vectors $`\mathbf{v}`$ and $`\mathbf{w}`$.**

**Example (ex6):** $`\mathbf{v} = \begin{pmatrix} \cos\alpha \\ \sin\alpha \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} \cos\beta \\ \sin\beta \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \cos\alpha\cos\beta + \sin\alpha\sin\beta
```

By trigonometry: $`= \cos(\beta - \alpha)`$

The angle between the vectors is $`\theta = \beta - \alpha`$.

**The sign of $`\mathbf{v} \cdot \mathbf{w}`$ tells whether we are below or above a right angle:**

| Case | Angle | Dot Product |
|:-----|:------|:------------|
| i) | $`\theta < 90°`$ | $`\mathbf{v} \cdot \mathbf{w} > 0`$ |
| ii) | $`\theta = 90°`$ | $`\mathbf{v} \cdot \mathbf{w} = 0`$ |
| iii) | $`\theta > 90°`$ | $`\mathbf{v} \cdot \mathbf{w} < 0`$ |

We know $`|\cos\theta| \leq 1`$. So, for nonzero $`\mathbf{v}`$ and $`\mathbf{w}`$, we can measure the angle by:

```math
\frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \cos\theta
```

### 1.2.6 Schwarz Inequality

```math
\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\|\|\mathbf{w}\|\cos\theta
```

```math
|\mathbf{v} \cdot \mathbf{w}| = \|\mathbf{v}\|\|\mathbf{w}\||\cos\theta| \leq \|\mathbf{v}\|\|\mathbf{w}\|
```

> **Schwarz Inequality (Cauchy--Schwarz--Bunyakovsky):**
> ```math
> |\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|
> ```

**Example (ex7):** Find $`\cos\theta`$ for $`\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = 2 \cdot 1 + 1 \cdot 2 = 4
```

```math
\|\mathbf{v}\| = \sqrt{5}, \quad \|\mathbf{w}\| = \sqrt{5}
```

```math
\cos\theta = \frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \frac{4}{5}
```

By Schwarz inequality: $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$, i.e., $`4 < 5`$. $`\checkmark`$

**Example (ex8):** $`\mathbf{v} = \begin{pmatrix} a \\ b \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} b \\ a \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = ab + ba = 2ab
```

```math
\|\mathbf{v}\| = \|\mathbf{w}\| = \sqrt{a^2 + b^2}
```

The Schwarz inequality $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$ is:

```math
|2ab| \leq a^2 + b^2
```

```math
\Leftrightarrow 0 \leq a^2 + b^2 - 2|ab|
```

- If $`ab \geq 0`$: $`a^2 + b^2 - 2ab = (a - b)^2 \geq 0`$ $`\checkmark`$
- If $`ab < 0`$: $`a^2 + b^2 + 2ab = (a + b)^2 \geq 0`$ $`\checkmark`$

### 1.2.7 Triangle Inequality

**Triangle inequality comes directly from Schwarz inequality.**

```math
\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w}) = \mathbf{v} \cdot \mathbf{v} + 2\mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{w}
```

```math
\leq \|\mathbf{v}\|^2 + 2|\mathbf{v} \cdot \mathbf{w}| + \|\mathbf{w}\|^2
```

```math
\leq \|\mathbf{v}\|^2 + 2\|\mathbf{v}\|\|\mathbf{w}\| + \|\mathbf{w}\|^2 \quad \text{(by Schwarz)}
```

```math
= (\|\mathbf{v}\| + \|\mathbf{w}\|)^2
```

> **Triangle Inequality:**
> ```math
> \|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
> ```

**Verification from ex7:** By triangle inequality:

```math
\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
```

```math
\sqrt{3^2 + 3^2} = \sqrt{18} = 3\sqrt{2} < 2\sqrt{5}
```

```math
\sqrt{18} < \sqrt{20} \quad \checkmark
```

### 1.2.8 Planes in 3D

**Plane in 3D.** Suppose $`\mathbf{n}`$ is a unit vector $`\mathbf{n} = \begin{pmatrix} n_1 \\ n_2 \\ n_3 \end{pmatrix}`$.

Look at all vectors $`\mathbf{w} \in \mathbb{R}^3`$ such that $`\mathbf{w} \perp \mathbf{n}`$, i.e., $`\mathbf{w} \cdot \mathbf{n} = 0`$.

The vectors $`\mathbf{w}`$ with $`\mathbf{w} \cdot \mathbf{n} = 0`$ fill a **2D plane** in $`\mathbb{R}^3`$.

```math
\mathbf{w} \cdot \mathbf{n} = w_1 n_1 + w_2 n_2 + w_3 n_3 = 0
```

**Example:** $`\mathbf{n} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$

```math
\mathbf{w} \cdot \mathbf{n} = w_3 = 0
```

or $`z = 0`$, which represents the $`xy`$ plane.

---

<br>
