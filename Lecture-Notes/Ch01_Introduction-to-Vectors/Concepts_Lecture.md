# Chapter 1 — Vectors and Matrices

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-10-27
>
> Introduction to Linear Algebra, Strang, Ch 1

> **Learning Objectives:**
> 1. Define vectors and perform linear combinations
> 2. Compute dot products and understand the geometric meaning of orthogonality
> 3. Perform basic matrix-vector multiplication
> 4. Explain the meaning of the column space, rank, and basis
> 5. Compute the matrix product $`AB`$ from several viewpoints and find the factorization $`A = CR`$

---

<br>

## Table of Contents

- [1. Vectors and Matrices (Warmup)](#1-vectors-and-matrices-warmup)
  - [1.1 Key Ideas](#11-key-ideas)
  - [1.2 Length of a Vector](#12-length-of-a-vector)
  - [1.3 Dot Product](#13-dot-product)
  - [1.4 Matrix from Column Vectors](#14-matrix-from-column-vectors)
  - [1.5 Matrix-Vector Product as Linear Combination](#15-matrix-vector-product-as-linear-combination)
  - [1.6 Dependent Vectors](#16-dependent-vectors)
  - [1.7 Matrix Multiplication Preview](#17-matrix-multiplication-preview)
- [2. Vectors and Linear Combinations (1.1)](#2-vectors-and-linear-combinations-11)
  - [2.1 Linear Combinations in R^2](#21-linear-combinations-in-r2)
  - [2.2 Two Key Questions](#22-two-key-questions)
  - [2.3 Extending to Higher Dimensions](#23-extending-to-higher-dimensions)
  - [2.4 Solving Two Equations](#24-solving-two-equations)
  - [2.5 Can Elimination Fail?](#25-can-elimination-fail)
  - [2.6 Vectors in Three Dimensions](#26-vectors-in-three-dimensions)
- [3. Lengths and Angles from Dot Products (1.2)](#3-lengths-and-angles-from-dot-products-12)
  - [3.1 Dot Product Definition](#31-dot-product-definition)
  - [3.2 Length of a Vector](#32-length-of-a-vector)
  - [3.3 Unit Vectors](#33-unit-vectors)
  - [3.4 Perpendicular Vectors](#34-perpendicular-vectors)
  - [3.5 The Angle Between Two Vectors](#35-the-angle-between-two-vectors)
  - [3.6 Schwarz Inequality](#36-schwarz-inequality)
  - [3.7 Triangle Inequality](#37-triangle-inequality)
  - [3.8 Planes in 3D](#38-planes-in-3d)
- [4. Matrices and Their Column Spaces (1.3)](#4-matrices-and-their-column-spaces-13)
  - [4.1 Matrix-Vector Multiplication](#41-matrix-vector-multiplication)
  - [4.2 Column Space](#42-column-space)
  - [4.3 Independent, Dependent, and Column Space](#43-independent-dependent-and-column-space)
  - [4.4 Span](#44-span)
  - [4.5 Rank and Basis](#45-rank-and-basis)
  - [4.6 Matrices of Rank One](#46-matrices-of-rank-one)
- [5. Matrix Multiplication AB and CR (1.4)](#5-matrix-multiplication-ab-and-cr-14)
  - [5.1 Rules for Matrix Multiplication](#51-rules-for-matrix-multiplication)
  - [5.2 Column Interpretation of AB](#52-column-interpretation-of-ab)
  - [5.3 Computational Cost](#53-computational-cost)
  - [5.4 Properties of Matrix Multiplication](#54-properties-of-matrix-multiplication)
  - [5.5 Rank One Matrices and A = CR](#55-rank-one-matrices-and-a--cr)
  - [5.6 Finding C and R](#56-finding-c-and-r)
  - [5.7 Columns of A times Rows of B (Outer Product)](#57-columns-of-a-times-rows-of-b-outer-product)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Vectors and Matrices (Warmup)

**Chapter 1 of the textbook (Sections 2-5 of these notes, in order):**
- 1.1 Vectors and Linear Combination
- 1.2 Lengths and Angles from Dot Products
- 1.3 Matrices and their Column Spaces
- 1.4 Matrix Multiplication $`AB`$ and $`CR`$

### 1.1 Key Ideas

- Linear algebra is about **vectors** $`\mathbf{v}`$ and **matrices** $`A`$.
- We define operations: $`+, -, \cdot`$ (addition, subtraction, multiplication).
- Consider the following vectors and scalars:

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

> **Key Point:** Since $`\mathbf{v}`$ and $`\mathbf{w}`$ do not lie on one line, the linear combinations $`c\mathbf{v} + d\mathbf{w}`$ of the two vectors fill the whole $`xy`$ plane.

### 1.2 Length of a Vector

The length of a vector $`\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} \in \mathbb{R}^2`$:

```math
\|\mathbf{v}\| = \sqrt{v_1^2 + v_2^2}
```

**Example:** $`\mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}`$, $`\|\mathbf{w}\| = \sqrt{1^2 + 3^2} = \sqrt{10}`$

### 1.3 Dot Product

The dot product of $`\mathbf{v}`$ and $`\mathbf{w}`$ is:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}^T \begin{pmatrix} w_1 \\ w_2 \end{pmatrix} = v_1 w_1 + v_2 w_2
```

**Example:**

```math
\begin{pmatrix} 2 \\ 4 \end{pmatrix} \cdot \begin{pmatrix} 1 \\ 3 \end{pmatrix} = (2)(1) + (4)(3) = 2 + 12 = 14
```

### 1.4 Matrix from Column Vectors

A matrix $`A`$ contains two columns:

```math
A = \begin{pmatrix} \mathbf{v} & \mathbf{w} \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}
```

### 1.5 Matrix-Vector Product as Linear Combination

Multiply the matrix $`A`$ by a vector $`\begin{pmatrix} c \\ d \end{pmatrix}`$:

```math
A \begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix} \begin{pmatrix} c \\ d \end{pmatrix} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

This is a **linear combination of $`\mathbf{v}`$ and $`\mathbf{w}`$**.

Let $`\mathbf{x} = \begin{pmatrix} c \\ d \end{pmatrix}`$. All combinations $`A\mathbf{x}`$ form the **column space** of the matrix $`A`$. In this example, the column space is the whole $`xy`$ plane.

### 1.6 Dependent Vectors

Let $`\mathbf{z} = \mathbf{v} + \mathbf{w}`$.

```math
B = \begin{pmatrix} \mathbf{v} & \mathbf{w} & \mathbf{z} \end{pmatrix} = \begin{pmatrix} 2 & 1 & 3 \\ 4 & 3 & 7 \end{pmatrix}
```

- The column space of $`B`$ is still the $`xy`$ plane.
- $`\mathbf{v}`$ and $`\mathbf{w}`$ are **independent** vectors.
- $`\mathbf{z}`$ is a **dependent** vector, since it is a combination of $`\mathbf{v}`$ and $`\mathbf{w}`$.

### 1.7 Matrix Multiplication Preview

A matrix multiplication $`AB`$ can be interpreted as **$`A`$ times each column of $`B`$**.

---

<br>

## 2. Vectors and Linear Combinations (1.1)

### 2.1 Linear Combinations in R^2

**(1)** $`2\mathbf{v} - 3\mathbf{w}`$ is a linear combination $`c\mathbf{v} + d\mathbf{w}`$ of the vectors $`\mathbf{v}`$ and $`\mathbf{w}`$.

**(2)** Let $`\mathbf{v} = \begin{pmatrix} 4 \\ 1 \end{pmatrix}`$ and $`\mathbf{w} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, then:

```math
2\mathbf{v} - 3\mathbf{w} = 2\begin{pmatrix} 4 \\ 1 \end{pmatrix} - 3\begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}
```

**(3)** All combinations $`c\begin{pmatrix} 4 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$ fill the $`xy`$ plane.

**(4)** The vectors $`c\begin{pmatrix} 4 \\ 1 \\ 0 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}`$ fill a **plane** in $`xyz`$ space.

$`\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$ is **not** on the plane.

#### Building Linear Combinations

Linear combinations of vectors $`\mathbf{v}`$ and $`\mathbf{w}`$ are built from two basic operations:
1. **Scalar multiplication:** $`c\mathbf{v}`$, $`d\mathbf{w}`$
2. **Vector addition:** $`\mathbf{v} + \mathbf{w}`$

This leads to a **linear combination** $`c\mathbf{v} + d\mathbf{w}`$ of $`\mathbf{v}`$ and $`\mathbf{w}`$.

### 2.2 Two Key Questions

It opens up two questions:

**(1) Describe** all the combinations $`c\mathbf{v} + d\mathbf{w}`$.
- The result is either a **plane** or a **line**. If both vectors are zero, it is the single point at the origin.

**(2) Find** the numbers $`c`$ and $`d`$ such that $`c\mathbf{v} + d\mathbf{w} = \mathbf{x}`$.

**Example:** Find $`c`$ and $`d`$ such that the following holds (Self-Check Question 1 gives the solution):

```math
c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 4 \\ 3 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}
```

### 2.3 Extending to Higher Dimensions

So far $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$.

Let us now extend the dimension of the vectors and increase the number of vectors:

```math
\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \in \mathbb{R}^m
```

These are $`n`$ vectors in $`m`$-dimensional space.

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}
```

$`A`$ is an $`m`$ by $`n`$ matrix with $`m`$ rows and $`n`$ columns.

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

#### Linear Combinations $`c\mathbf{v} + d\mathbf{w}`$

Let $`\mathbf{v} \in \mathbb{R}^2`$. It has two components: $`\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}`$.

Geometrically, $`\mathbf{v} + \mathbf{w} = \mathbf{w} + \mathbf{v}`$, and the parallelogram rule shows this commutativity (Figure 1).

```text
                    v
         w ●─────────────────● v + w = w + v
          ╱                 ╱
       w ╱                 ╱ w
        ╱                 ╱
     0 ●─────────────────● v
                v
```

*Figure 1. Parallelogram showing that vector addition is commutative*

If $`\mathbf{v} \neq \mathbf{0}`$, the vectors $`c\mathbf{v}`$ fill an infinitely long line in the $`xy`$ plane.

If $`\mathbf{w}`$ is not on that line, then the vectors $`d\mathbf{w}`$ fill a second line.

**Then the linear combinations $`c\mathbf{v} + d\mathbf{w}`$ fill the plane.**

**Special cases:**
- $`1\mathbf{v} + 1\mathbf{w}`$ is the sum of the two vectors.
- $`1\mathbf{v} - 1\mathbf{w}`$ is the difference of the two vectors.
- $`0\mathbf{v} + 0\mathbf{w}`$ is the zero vector.
- $`c\mathbf{v} + 0\mathbf{w}`$ is the vector $`c\mathbf{v}`$ in the direction of $`\mathbf{v}`$.

### 2.4 Solving Two Equations

**Solve the following equation:**

```math
c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \quad \cdots (*)
```

This is equivalent to the system:

```math
\begin{cases} 2c + 2d = 8 \\ c - d = 2 \end{cases}
```

**Solution by elimination:**

Divide the first equation by 2 to get $`c + d = 4`$. Then write the two equations together:

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

- The two linear equations meet at the point $`(c, d)`$.
- $`\mathbf{v}`$ and $`\mathbf{w}`$ are **linearly independent**.
- $`A = \begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}`$ is **invertible**.

### 2.5 Can Elimination Fail?

**Yes**, when $`\mathbf{v} \| \mathbf{w}`$ (i.e., $`v_1 w_2 = v_2 w_1`$, which is $`\frac{v_1}{v_2} = \frac{w_1}{w_2}`$ when the components are nonzero).

The system becomes:

```math
xv_1 + yw_1 = b_1
```

```math
xv_2 + yw_2 = b_2
```

Two possibilities:
- **No solution**: All combinations of $`\mathbf{v}`$ and $`\mathbf{w}`$ lie on the same line. If $`\mathbf{b}`$ is not on the line, there is no solution. No combination of $`\mathbf{v}`$ and $`\mathbf{w}`$ equals $`\mathbf{b}`$.
- **Infinitely many solutions**: If $`\mathbf{b}`$ is on the line, there are infinitely many solutions.

### 2.6 Vectors in Three Dimensions

Suppose $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^3`$ (three-dimensional space).

```math
\mathbf{v} = \begin{pmatrix} 2 \\ 3 \\ 1 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 4 \\ 1 \end{pmatrix}
```

```math
c\mathbf{v} + d\mathbf{w} = \begin{pmatrix} 2c + d \\ 3c + d \\ c \end{pmatrix}
```

> **Key Point:** $`c\mathbf{v} + d\mathbf{w}`$ **does not** fill the whole 3D space. It can fill at most a **2D plane**.

**Filling 3D space requires three independent vectors.**

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

> **Definition:** The **identity matrix** $`I`$ is the square matrix with ones on the diagonal and zeros elsewhere. Every vector $`\mathbf{v}`$ satisfies $`\mathbf{v} = I\mathbf{v}`$.

#### How Do We Know If It Is a Plane?

Suppose nonzero vectors $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^3`$ are **independent**, that is, $`\mathbf{w}`$ is not a multiple $`c\mathbf{v}`$ of $`\mathbf{v}`$.

Then their linear combinations fill a **plane** inside 3D space (a flat surface).

---

<br>

## 3. Lengths and Angles from Dot Products (1.2)

### 3.1 Dot Product Definition

**(1)** Let $`\mathbf{v} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 4 \\ 6 \end{pmatrix}`$.

The dot product of $`\mathbf{v}`$ and $`\mathbf{w}`$ is:

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

**Example:** $`\mathbf{v} \in \mathbb{R}^3`$, $`\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + v_3^2 = (v_1^2 + v_2^2) + v_3^2`$ where $`(v_1^2 + v_2^2)`$ is the squared length of the part in the $`xy`$ plane.

### 3.2 Length of a Vector

**(2)** The length squared of $`\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}`$ is:

```math
\|\mathbf{v}\|^2 = \mathbf{v} \cdot \mathbf{v} = 1^2 + 3^2 + 2^2 = 1 + 9 + 4 = 14
```

```math
\therefore \|\mathbf{v}\| = \sqrt{14}
```

### 3.3 Unit Vectors

- $`\mathbf{v}`$ is a **unit vector** when $`\|\mathbf{v}\| = 1`$.
- If $`\mathbf{v} \neq \mathbf{0}`$, then $`\frac{\mathbf{v}}{\|\mathbf{v}\|}`$ is a unit vector.

**Example (ex1):** $`\mathbf{u} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}`$ is a unit vector.

```math
\|\mathbf{u}\|^2 = \cos^2\theta + \sin^2\theta = 1
```

### 3.4 Perpendicular Vectors

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

### 3.5 The Angle Between Two Vectors

**(4)** The angle $`\theta = 45°`$ between $`\mathbf{v} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$ and $`\mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$ has:

```math
\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|} = \frac{1}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}
```

**(5)** All angles have $`|\cos\theta| \leq 1`$.

All vectors have $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$ and $`\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|`$.

**Example (ex5): Application to Total Cost**

Price vector $`\mathbf{p} = \begin{pmatrix} p_1 \\ p_2 \\ p_3 \end{pmatrix}`$, quantity vector $`\mathbf{q} = \begin{pmatrix} q_1 \\ q_2 \\ q_3 \end{pmatrix}`$

$`p_i q_i`$ is the cost of buying $`q_i`$ units at the price $`p_i`$.

```math
\mathbf{p} \cdot \mathbf{q} = p_1 q_1 + p_2 q_2 + p_3 q_3 \Rightarrow \text{total cost}
```

**The dot product $`\mathbf{v} \cdot \mathbf{w}`$ finds the angle between any two nonzero vectors $`\mathbf{v}`$ and $`\mathbf{w}`$.**

**Example (ex6):** $`\mathbf{v} = \begin{pmatrix} \cos\alpha \\ \sin\alpha \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} \cos\beta \\ \sin\beta \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \cos\alpha\cos\beta + \sin\alpha\sin\beta
```

By the trigonometric addition formula, this equals $`\cos(\beta - \alpha)`$.

So, if $`0 \leq \beta - \alpha \leq \pi`$, the angle between the vectors is $`\theta = \beta - \alpha`$.

**For nonzero vectors, the sign of $`\mathbf{v} \cdot \mathbf{w}`$ tells whether the angle is below or above a right angle:**

| Case | Angle | Dot Product |
|:-----|:------|:------------|
| i) | $`\theta < 90°`$ | $`\mathbf{v} \cdot \mathbf{w} > 0`$ |
| ii) | $`\theta = 90°`$ | $`\mathbf{v} \cdot \mathbf{w} = 0`$ |
| iii) | $`\theta > 90°`$ | $`\mathbf{v} \cdot \mathbf{w} < 0`$ |

*Table 1. Sign of the dot product and the angle between two vectors*

We know $`|\cos\theta| \leq 1`$. So, for nonzero $`\mathbf{v}`$ and $`\mathbf{w}`$, we can measure the angle by:

```math
\frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \cos\theta
```

### 3.6 Schwarz Inequality

```math
\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\|\|\mathbf{w}\|\cos\theta
```

```math
|\mathbf{v} \cdot \mathbf{w}| = \|\mathbf{v}\|\|\mathbf{w}\||\cos\theta| \leq \|\mathbf{v}\|\|\mathbf{w}\|
```

> **Key Point:** **Schwarz inequality (Cauchy-Schwarz-Bunyakovsky).** For all vectors $`\mathbf{v}, \mathbf{w}`$, the following holds.
> ```math
> |\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|
> ```

**Example (ex7):** Find $`\cos\theta`$ for $`\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$.

```math
\mathbf{v} \cdot \mathbf{w} = 2 \cdot 1 + 1 \cdot 2 = 4
```

```math
\|\mathbf{v}\| = \sqrt{5}, \quad \|\mathbf{w}\| = \sqrt{5}
```

```math
\cos\theta = \frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \frac{4}{5}
```

The Schwarz inequality $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$ holds as $`4 < 5`$. $`\checkmark`$

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

### 3.7 Triangle Inequality

**The triangle inequality comes directly from the Schwarz inequality.**

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

> **Key Point:** **Triangle inequality.** For all vectors $`\mathbf{v}, \mathbf{w}`$, the following holds.
> ```math
> \|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
> ```

**Verification with the vectors of ex7:** Since $`\mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 3 \end{pmatrix}`$, the triangle inequality holds as follows:

```math
\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
```

```math
\sqrt{3^2 + 3^2} = \sqrt{18} = 3\sqrt{2} < 2\sqrt{5}
```

```math
\sqrt{18} < \sqrt{20} \quad \checkmark
```

### 3.8 Planes in 3D

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

that is, $`z = 0`$, which represents the $`xy`$ plane.

---

<br>

## 4. Matrices and Their Column Spaces (1.3)

### 4.1 Matrix-Vector Multiplication

**(1)** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}`$ is a **3 by 2 matrix** with 3 rows and 2 columns, and its rank is 2.

**(2)** The 3 components of $`A\mathbf{x}`$ are dot products of the 3 rows of $`A`$ with the vector $`\mathbf{x}`$.

**Example:** $`\mathbf{x} = \begin{pmatrix} 7 \\ 8 \end{pmatrix}`$

```math
A\mathbf{x} = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 7 \\ 8 \end{pmatrix} = \begin{pmatrix} 7 + 16 \\ 21 + 32 \\ 35 + 48 \end{pmatrix} = \begin{pmatrix} 23 \\ 53 \\ 83 \end{pmatrix}
```

**(3)** $`A\mathbf{x}`$ is a **combination of the columns** of $`A`$:

```math
\begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 7 \\ 8 \end{pmatrix} = 7\begin{pmatrix} 1 \\ 3 \\ 5 \end{pmatrix} + 8\begin{pmatrix} 2 \\ 4 \\ 6 \end{pmatrix}
```

**(4)** The column space of $`A`$ contains all combinations $`A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2`$ of columns.

**(5)** In a rank one matrix, all columns of $`A`$ are on one line.

#### Consider $`m`$ by $`n`$ Matrices

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}
```

If $`m = n`$, $`A`$ is a **square matrix**.

**Examples of square matrices ($`A \in \mathbb{R}^{3 \times 3}`$):**

```math
\text{Identity matrix: } \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \qquad \text{Diagonal matrix: } \begin{pmatrix} 2 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 0 & 5 \end{pmatrix}
```

```math
\text{Triangular matrix: } \begin{pmatrix} 2 & 1 & -3 \\ 0 & 4 & 7 \\ 0 & 0 & 5 \end{pmatrix} \qquad \text{Symmetric matrix: } \begin{pmatrix} 2 & 1 & -3 \\ 1 & 4 & 7 \\ -3 & 7 & 5 \end{pmatrix}
```

#### Interpreting Columns as Vectors

We can interpret the columns of $`A`$ as vectors:

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}
```

where $`\mathbf{a}_i \in \mathbb{R}^m`$, $`i = 1, 2, \ldots, n`$.

**Question:** How does an $`m`$ by $`n`$ matrix $`A`$ multiply an $`n`$ by 1 vector $`\mathbf{x}`$?

1. By **dot products** of $`\mathbf{x}`$ with the rows of $`A`$.
2. As a **linear combination** of the columns of $`A`$.

#### Two Viewpoints of $`A\mathbf{x}`$

**(1) Row viewpoint (dot products):**

```math
A\mathbf{x} = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}
```

Three rows give three dot products.

Note that each row of $`A`$ has the same number of components as the vector $`\mathbf{x}`$.

**(2) Column viewpoint (linear combinations):**

```math
A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n
```

This is a **combination of the columns of $`A`$**.

**Example:** $`A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}`$

```math
A\mathbf{x} = x_1\begin{pmatrix} -1 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix} + x_4\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}
```

### 4.2 Column Space

When $`A = (\mathbf{a}_1 \ \mathbf{a}_2)`$ with linearly independent $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$:

```math
A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2
```

shows **the plane of all combinations** (the span).

### 4.3 Independent, Dependent, and Column Space

**Example (ex1):**

```math
A_1 = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 4 & 0 \\ 3 & 5 & 6 \end{pmatrix}
```

Each column gives a new direction. Their combinations fill **3D space**.

```math
A_1\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = x_1\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix} + x_2\begin{pmatrix} 0 \\ 4 \\ 5 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 0 \\ 6 \end{pmatrix}
```

**Independent columns:** $`A_1\mathbf{x} = \mathbf{0}`$ holds only when $`\mathbf{x} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$.

**Example (ex2):**

```math
A_2 = \begin{pmatrix} 1 & 2 & 3 \\ 1 & 4 & 5 \\ 6 & 0 & 6 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

$`\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2`$. The combinations of $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3`$ **do not** fill 3D space. $`\mathbf{a}_3`$ lies on the plane of $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ in 3D.

**Example (ex3):**

```math
A_3 = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 6 & 8 \\ 5 & 15 & 20 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

$`\mathbf{a}_2 = 3\mathbf{a}_1`$ and $`\mathbf{a}_3 = 4\mathbf{a}_1`$.

All three columns of $`A_3`$ lie on the **same line** in 3D space.

#### The Column Space $`C(A)`$

The column space $`C(A)`$ contains all vectors $`A\mathbf{x}`$, that is, all combinations of the columns.

**Thinking about the column space of $`A`$:**

```math
A_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix} \quad \text{(Independent columns)}
```

```math
A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix} \quad \text{Is the 4th column dependent?}
```

For $`A_5`$, $`\mathbf{a}_4 = \mathbf{a}_1 - \mathbf{a}_2 + \mathbf{a}_3`$, so the 4th column is dependent.

Consider $`\mathbf{v} = A_4\mathbf{x}`$:

```math
\mathbf{v} = x_1\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + x_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
```

Solving $`A_4\mathbf{x} = \mathbf{v}`$, the equations are:

```math
v_4 = x_4, \quad v_3 = x_4 + x_3, \quad v_2 = x_4 + x_3 + x_2, \quad v_1 = x_4 + x_3 + x_2 + x_1
```

So $`x_4 = v_4`$, $`x_3 = v_3 - v_4`$, $`x_2 = v_2 - v_3`$, and $`x_1 = v_1 - v_2`$.

```math
\therefore \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \\ v_4 \end{pmatrix} = (v_1 - v_2)\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + (v_2 - v_3)\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + (v_3 - v_4)\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + v_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
```

This implies that **every $`\mathbf{v}`$ is in the column space**. The four equations $`A_4\mathbf{x} = \mathbf{v}`$ are solved for every $`\mathbf{v}`$.

### 4.4 Span

The **span** describes all the linear combinations of a set of vectors.

**The span of the columns of $`A`$ is the column space.**

**Example:** $`\mathbf{a}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}`$, $`\mathbf{a}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} \in \mathbb{R}^4`$

- $`c_1\mathbf{a}_1`$ corresponds to a line in 4D.
- $`c_2\mathbf{a}_2`$ corresponds to a line in 4D.
- The combinations $`c_1\mathbf{a}_1 + c_2\mathbf{a}_2`$ fill a 2D plane in 4D space.
- The plane is the span of columns $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$.

**Example:** $`A_4`$ has 4 independent columns. The column space of $`A_4`$ is all of $`\mathbb{R}^4`$.

**Example:** $`A_5`$ has one dependent column (only three independent columns).

```math
A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix}
```

The column space of $`A_5`$ is a **3D subspace** inside $`\mathbb{R}^4`$. The 4th column is in that subspace.

We can only solve $`A_5\mathbf{x} = \mathbf{v}`$ when $`\mathbf{v} \in C(A_5)`$.

### 4.5 Rank and Basis

Let $`A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}`$.

$`C(A)`$ is the column space of $`A`$, $`\mathbf{a}_i \in \mathbb{R}^m`$.

The column space $`C(A)`$ might fill all of $`\mathbb{R}^m`$ or might not.

**Example:** Take $`m = 3`$ and $`A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}`$

$`C(A)`$ is one of the following:
- The whole space $`\mathbb{R}^3`$, if $`A`$ has 3 independent columns (I.C.).
- A plane in $`\mathbb{R}^3`$, if $`A`$ has 2 I.C.
- A line in $`\mathbb{R}^3`$, if $`A`$ has 1 I.C.
- A single point $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$, if $`A`$ is a zero matrix.

**More examples ($`3 \times 3`$ matrices):**

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}: C(A) = \mathbb{R}^3
```

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}: C(A) = xy \text{ plane}
```

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}: C(A) = x \text{ axis}
```

**Key Questions:**

- **Q1:** How many columns of $`A`$ are independent? That number $`r`$ is the **rank** of $`A`$.
- **Q2:** Which are the first $`r`$ independent columns? They are a **basis** for the column space.
- **Q3:** What combinations of those $`r`$ basis vectors produce the remaining $`n - r`$ columns?
- **Q4:** Write any $`A`$ as an $`m \times r`$ column matrix $`C`$ times an $`r \times n`$ matrix $`R`$: $`A = CR`$.
- **Q5:** The $`r`$ rows of $`R`$ are a basis for the **row space** of $`A`$. The rows of $`R`$ do not come directly from $`A`$.

### 4.6 Matrices of Rank One

For a rank one matrix, **all column vectors lie along the same line**.

**Example:**

```math
A_6 = \begin{pmatrix} 1 & 3 & -2 \\ 4 & 12 & -8 \\ 2 & 6 & -4 \end{pmatrix}
```

$`\mathbf{a}_2 = 3\mathbf{a}_1`$, $`\mathbf{a}_3 = -2\mathbf{a}_1`$.

- $`A_6`$ has rank $`r = 1`$.
- $`C(A_6)`$ consists of all $`c_1\mathbf{a}_1`$, a line passing through the origin.
- All rows of $`A_6`$ are multiples of one row.

> **Key Point:** When the column space is a single line in $`\mathbb{R}^m`$, the row space is a single line in $`\mathbb{R}^n`$.

**Question: Why does this happen?**

If all columns are in the same direction, then all rows are in the same direction.

**Example:** $`A = \begin{pmatrix} a & ma \\ b & mb \end{pmatrix}`$

If $`a \neq 0`$, row 2 is $`\frac{b}{a}`$ times row 1. If the column rank is 1, then the row rank is 1.

**Example:** $`A = \begin{pmatrix} a & ma & pa \\ b & mb & pb \\ c & mc & pc \end{pmatrix}`$

If $`a \neq 0`$, row 2 = $`\frac{b}{a}`$ row 1 and row 3 = $`\frac{c}{a}`$ row 1.

If the column rank is 1, then the row rank is 1.

> **Key Point:** **Does row rank equal column rank for every matrix?** Yes. For every matrix, the row rank equals the column rank.

---

<br>

## 5. Matrix Multiplication AB and CR (1.4)

### 5.1 Rules for Matrix Multiplication

**(1)** To multiply $`AB`$, we need $`A \in \mathbb{R}^{m \times n}`$ and $`B \in \mathbb{R}^{n \times p}`$.

**The row length of $`A`$ must equal the column length of $`B`$.**

**(2) Dot product viewpoint:**

```math
(AB)_{ij} = (\text{row } i \text{ of } A) \cdot (\text{column } j \text{ of } B)
```

**(3) Column viewpoint:**

```math
AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p)
```

Column $`j`$ of $`AB`$ is $`A\mathbf{b}_j`$.

**(4) Non-commutativity:** $`AB \neq BA`$ in general.

**(5)** If the $`r`$ independent columns of $`A`$ are collected in $`C`$, then $`A = CR`$, where $`C \in \mathbb{R}^{m \times r}`$ and $`R \in \mathbb{R}^{r \times n}`$.

### 5.2 Column Interpretation of AB

Let $`B = \begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix}`$.

```math
AB = A\begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix} = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}
```

Each column of $`AB`$ is a **combination of the columns of $`A`$**.

**Example:** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$, $`B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}`$

```math
AB = (A\mathbf{b}_1 \ A\mathbf{b}_2) = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}
```

```math
A\mathbf{b}_1 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}
```

```math
A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

**Example:** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$, $`B = \begin{pmatrix} 5 & 6 \\ 7 & 8 \end{pmatrix}`$

```math
A\mathbf{b}_1 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 5 \\ 7 \end{pmatrix}
```

**(1) Dot product approach:**

```math
= \begin{pmatrix} 1 \cdot 5 + 2 \cdot 7 \\ 3 \cdot 5 + 4 \cdot 7 \end{pmatrix} = \begin{pmatrix} 19 \\ 43 \end{pmatrix}
```

**(2) Column combination approach:**

```math
= 5\begin{pmatrix} 1 \\ 3 \end{pmatrix} + 7\begin{pmatrix} 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 19 \\ 43 \end{pmatrix}
```

Both approaches use 4 multiplications.

```math
A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 6 \\ 8 \end{pmatrix} = \begin{pmatrix} 6 + 16 \\ 18 + 32 \end{pmatrix} = \begin{pmatrix} 22 \\ 50 \end{pmatrix}
```

```math
\therefore AB = \begin{pmatrix} 19 & 22 \\ 43 & 50 \end{pmatrix}
```

### 5.3 Computational Cost

Let $`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$.

```math
AB = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}
```

**(1) Dot product count:** $`A\mathbf{b}_1`$ takes $`m`$ dot products, and $`AB`$ takes $`mp`$ dot products. Each dot product has $`n`$ multiplications, so **$`mnp`$ multiplications** are needed in total.

**(2) Column combination count:** In $`A\mathbf{b}_1 = b_{11}\mathbf{a}_1 + b_{21}\mathbf{a}_2 + \cdots + b_{n1}\mathbf{a}_n`$, each term takes $`m`$ multiplications, so one column takes $`mn`$ multiplications. So $`AB`$ needs **$`mnp`$ multiplications**.

> **Key Point:** The cost of matrix multiplication is $`mnp`$.
>
> If $`A, B \in \mathbb{R}^{n \times n}`$, then the cost of $`AB`$ is $`n^3`$.

### 5.4 Properties of Matrix Multiplication

**Non-commutativity:**

```math
AB \neq BA \quad \text{in general.}
```

Matrix multiplication is **not commutative**.

**Associative law:**

```math
(AB)C = A(BC)
```

For matrix multiplication, the associative law **is** true.

**Distributive law:**

```math
A(B + C) = AB + AC
```

### 5.5 Rank One Matrices and A = CR

**Review of $`AB`$:**

**(1) Dot product:** $`(AB)_{ij} = (\text{row } i \text{ of } A) \cdot (\text{column } j \text{ of } B)`$. The product $`(m \times n)(n \times p)`$ needs $`mp`$ dot products.

**(2) Combine columns:** $`AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p) = (A\mathbf{b}_1 \ A\mathbf{b}_2 \ \cdots \ A\mathbf{b}_p)`$. Column $`j`$ of $`AB`$ is $`A\mathbf{b}_j`$.

**Rank one matrices and $`A = CR`$:**

All columns of a rank one matrix lie on the same line.

When all the columns of $`A`$ are in the same column direction, then all the rows of $`A`$ are in the same row direction.

**Example:** A matrix of rank $`r = 1`$.

```math
A = \begin{pmatrix} 1 & 2 & 10 & 100 \\ 3 & 6 & 30 & 300 \\ 2 & 4 & 20 & 200 \end{pmatrix}
```

It has one independent column and one independent row.

If the column space of $`A`$ is a line, the row space of $`A`$ is also a line.

We can decompose $`A`$ into:

```math
A = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}\begin{pmatrix} 1 & 2 & 10 & 100 \end{pmatrix} = CR
```

where $`C \in \mathbb{R}^{3 \times 1}`$ and $`R \in \mathbb{R}^{1 \times 4}`$.

### 5.6 Finding C and R

**$`C`$ contains the first $`r`$ independent columns of $`A`$.**

Given $`A`$, we look for independent columns **from left to right**:

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}
```

1. If $`\mathbf{a}_1 \neq \mathbf{0}`$, then put $`\mathbf{a}_1`$ into $`C`$.
2. If $`\mathbf{a}_2 \neq c_1\mathbf{a}_1`$ (not a multiple of $`\mathbf{a}_1`$), then put $`\mathbf{a}_2`$ into $`C`$.
3. If $`\mathbf{a}_3 \neq c_1\mathbf{a}_1 + c_2\mathbf{a}_2`$ (not a combination of $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$), then put $`\mathbf{a}_3`$ into $`C`$.
4. Continue until $`C`$ has $`r`$ columns.

**The number $`r`$ is the rank of $`A`$ and $`C`$.**

So $`C\mathbf{x} = \mathbf{0}`$ only when $`\mathbf{x} = \mathbf{0}`$. No nontrivial combination of the columns of $`C`$ gives the zero vector.

**Example:**

```math
A = \begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

1. Start with $`C = \begin{pmatrix} 2 \\ 4 \\ 1 \end{pmatrix}`$.
2. $`\mathbf{a}_2 = 3\mathbf{a}_1`$, so it is not put into $`C`$.
3. $`\mathbf{a}_3 \neq c_1\mathbf{a}_1`$, so $`C = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}`$.

$`\therefore`$ The rank is $`r = 2`$.

**Question:** What is $`R`$ for $`A = CR`$?

```math
\begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 3 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
3 \times 3 = (3 \times 2)(2 \times 3)
```

Rearranging the columns so that the independent columns come first:

```math
\begin{pmatrix} 2 & 4 & 6 \\ 4 & 8 & 12 \\ 1 & 5 & 3 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 0 \end{pmatrix}
```

**We observe that:**

1. $`C`$ contains a full set of $`r`$ independent columns of $`A`$.
2. After the rearrangement, $`R = (I \ F)`$ contains the identity matrix $`I \in \mathbb{R}^{r \times r}`$.
3. The dependent columns of $`A`$ are combinations of the independent columns in $`C`$.
4. Then $`A = CR = C(I \ F) = (C \ CF)`$.
5. $`C`$ has the same column space as $`A`$. $`R`$ has the same row space as $`A`$.

**Example (ex9):**

```math
\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 4 & 5 \\ 7 & 8 \end{pmatrix}\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}
```

This matrix has rank $`r = 2`$.

Column $`j`$ of $`A`$ is $`\mathbf{a}_j = C\mathbf{r}_j`$, where $`\mathbf{r}_j`$ is column $`j`$ of $`R`$.

$`(1 \ 2 \ 3) = (1 \ 2)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`(4 \ 5 \ 6) = (4 \ 5)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`(7 \ 8 \ 9) = (7 \ 8)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

Row $`i`$ of $`A`$ is row $`i`$ of $`C`$ times $`R`$, that is, a combination of the rows of $`R`$.

**Example:** $`A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 1 & 2 & 4 & 5 \end{pmatrix}_{2 \times 4}`$

```math
= \begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}\begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}
```

$`C \in \mathbb{R}^{2 \times 2}`$, $`R \in \mathbb{R}^{2 \times 4}`$. The rank is $`r = 2`$, and we recover all the columns of $`A`$ from $`C`$ by using $`R`$.

#### How to Find the Matrix R

We can use **elimination**. Chapter 2 covers elimination, and Chapter 3 covers finding $`R`$ by elimination in detail.

**Example:**

```math
A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}
```

**Step 1:** $`R2 - 2R1`$, $`R3 - 3R1`$:

```math
\begin{pmatrix} 1 & 3 & 4 \\ 0 & -2 & -6 \\ 0 & -2 & -6 \end{pmatrix}
```

**Step 2:** $`R3 - R2`$:

```math
\begin{pmatrix} 1 & 3 & 4 \\ 0 & -2 & -6 \\ 0 & 0 & 0 \end{pmatrix} \quad \text{row rank } r = 2
```

**Step 3:** $`R2 / (-2)`$:

```math
\begin{pmatrix} 1 & 3 & 4 \\ 0 & 1 & 3 \\ 0 & 0 & 0 \end{pmatrix}
```

**Step 4:** $`R1 - 3R2`$:

```math
\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \\ 0 & 0 & 0 \end{pmatrix} \quad \leftarrow R \text{ matrix}
```

Columns $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ are independent, and removing the zero row of the last matrix gives the $`R`$ of $`A = CR`$.

```math
A = \begin{pmatrix} 1 & 3 \\ 2 & 4 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \end{pmatrix} = CR
```

#### Key Properties of A = CR

- $`A = CR`$.
- The $`r`$ columns of $`C`$ are a **basis** for the column space of $`A`$.
- The $`r`$ rows of $`R`$ are a **basis** for the row space of $`A`$.
- Therefore both the column space and the row space have dimension $`r`$.

> **Key Point:** In $`A = CR`$, $`R`$ is the **reduced row echelon form** of $`A`$ with its zero rows removed.

**Example:**

```math
A = \begin{pmatrix} | & | & | \\ \mathbf{a}_1 & \mathbf{a}_2 & 3\mathbf{a}_1 + 4\mathbf{a}_2 \\ | & | & | \end{pmatrix}_{m \times 3}
```

```math
= \begin{pmatrix} | & | \\ \mathbf{a}_1 & \mathbf{a}_2 \\ | & | \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 4 \end{pmatrix}
```

```math
= \mathbf{a}_1(1 \ 0 \ 3) + \mathbf{a}_2(0 \ 1 \ 4) = (\mathbf{a}_1 \ \mathbf{a}_2 \ 3\mathbf{a}_1 + 4\mathbf{a}_2) = A
```

### 5.7 Columns of A times Rows of B (Outer Product)

```math
AB = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}_{m \times n} \begin{pmatrix} - & \mathbf{b}_1^* & - \\ - & \mathbf{b}_2^* & - \\ & \vdots & \\ - & \mathbf{b}_n^* & - \end{pmatrix}_{n \times p}
```

```math
= \mathbf{a}_1\mathbf{b}_1^* + \mathbf{a}_2\mathbf{b}_2^* + \cdots + \mathbf{a}_n\mathbf{b}_n^*
```

where $`\mathbf{b}_k^*`$ is row $`k`$ of $`B`$ and $`\mathbf{a}_k\mathbf{b}_k^* = \begin{pmatrix} | \\ \mathbf{a}_k \\ | \end{pmatrix}_{m \times 1}(- \ \mathbf{b}_k^* \ -)_{1 \times p}`$

This is a **column times row = rank 1 matrix** (when $`\mathbf{a}_k`$ and $`\mathbf{b}_k^*`$ are nonzero), with $`mp`$ entries.

```math
AB = \sum_{k=1}^{n} \mathbf{a}_k\mathbf{b}_k^*
```

So $`AB`$ is the sum of $`n`$ rank one matrices. Each needs $`mp`$ multiplications, so the total is $`mnp`$ multiplications.

> **Note:** For the matrix $`A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}`$, $`\text{rank}(A) = 2`$. Since $`A`$ is a $`3 \times 3`$ matrix whose rank is less than 3, the following holds.
>
> $`\Rightarrow`$ $`A`$ has no inverse matrix.
> $`\Rightarrow`$ The determinant of $`A`$ is zero (Chapter 5 covers determinants).

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Linear Combination | $`c\mathbf{v} + d\mathbf{w}`$ is built from scalar multiplication and vector addition. |
| Column Space $`C(A)`$ | It consists of all vectors $`A\mathbf{x}`$, that is, all combinations of the columns of $`A`$. |
| Span | It is the set of all linear combinations of a set of vectors. |
| Independent Vectors | The columns of $`A`$ are independent when $`A\mathbf{x} = \mathbf{0}`$ holds only for $`\mathbf{x} = \mathbf{0}`$. |
| Dependent Vector | It can be written as a combination of other vectors. |
| Rank | It is the number of independent columns, which equals the number of independent rows. |
| Basis | The first $`r`$ independent columns form a basis for $`C(A)`$. |
| Dot Product | $`\mathbf{v} \cdot \mathbf{w} = \sum v_i w_i`$ gives length and angle information. |
| Length (Norm) | The length is $`\Vert\mathbf{v}\Vert = \sqrt{\mathbf{v} \cdot \mathbf{v}}`$. |
| Unit Vector | A unit vector has $`\Vert\mathbf{v}\Vert = 1`$, and any nonzero vector is normalized as $`\mathbf{v}/\Vert\mathbf{v}\Vert`$. |
| Perpendicularity | $`\mathbf{v} \perp \mathbf{w} \Leftrightarrow \mathbf{v} \cdot \mathbf{w} = 0`$. |
| Angle Formula | For two nonzero vectors, $`\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\Vert\mathbf{v}\Vert\Vert\mathbf{w}\Vert}`$. |
| Schwarz Inequality | $`\vert\mathbf{v} \cdot \mathbf{w}\vert \leq \Vert\mathbf{v}\Vert\Vert\mathbf{w}\Vert`$. |
| Triangle Inequality | $`\Vert\mathbf{v} + \mathbf{w}\Vert \leq \Vert\mathbf{v}\Vert + \Vert\mathbf{w}\Vert`$. |
| Pythagorean Theorem | If $`\mathbf{v} \perp \mathbf{w}`$, then $`\Vert\mathbf{v} + \mathbf{w}\Vert^2 = \Vert\mathbf{v}\Vert^2 + \Vert\mathbf{w}\Vert^2`$. |
| Matrix Multiplication $`AB`$ | $`(AB)_{ij} = \text{row}_i(A) \cdot \text{col}_j(B)`$, and column $`j`$ of $`AB`$ is $`A\mathbf{b}_j`$. |
| Outer Product | $`AB = \sum \mathbf{a}_k \mathbf{b}_k^*`$ is a sum of rank-1 matrices. |
| Cost of $`AB`$ | It takes $`mnp`$ multiplications for $`A \in \mathbb{R}^{m \times n}`$ and $`B \in \mathbb{R}^{n \times p}`$. |
| Non-commutativity | $`AB \neq BA`$ in general. |
| Associativity | $`(AB)C = A(BC)`$. |
| Distributivity | $`A(B+C) = AB + AC`$. |
| $`A = CR`$ Factorization | $`C`$ holds the independent columns of $`A`$, and $`R`$ is the reduced row echelon form with its zero rows removed. |
| Rank One Matrix | All columns lie on one line, and $`A = \mathbf{c}\mathbf{r}^T`$. |
| Row Rank = Column Rank | This always holds for any matrix. |
| Identity Matrix | $`\mathbf{v} = I\mathbf{v}`$, and its columns are the standard basis vectors. |
| Plane via Normal | $`\mathbf{w} \cdot \mathbf{n} = 0`$ defines the plane through the origin perpendicular to $`\mathbf{n}`$. |
| Invertibility | For a square matrix $`A`$, $`A`$ is invertible $`\Leftrightarrow`$ its columns are independent $`\Leftrightarrow`$ $`A\mathbf{x} = \mathbf{b}`$ has a unique solution. |
| Elimination Failure | When $`\mathbf{v} \parallel \mathbf{w}`$, there is no solution or there are infinitely many solutions. |

---

<br>

## Self-Check Questions

1. **Coefficients of a linear combination:** Find $`c`$ and $`d`$ such that $`c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 4 \\ 3 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}`$ (the example in Section 2.2).

   <details>
   <summary>Answer</summary>

   The system is $`2c + 4d = 2`$ and $`c + 3d = -1`$. Dividing the first equation by 2 gives $`c + 2d = 1`$, and subtracting it from the second equation gives $`d = -2`$. So $`c = 1 - 2d = 5`$. Check: $`5\begin{pmatrix} 2 \\ 1 \end{pmatrix} - 2\begin{pmatrix} 4 \\ 3 \end{pmatrix} = \begin{pmatrix} 10 - 8 \\ 5 - 6 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}`$. The two vectors are not multiples of each other (they are independent), so the solution is unique.

   </details>

2. **Length, angle, and inequalities:** For $`\mathbf{v} = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}`$ and $`\mathbf{w} = \begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}`$, find $`\mathbf{v} \cdot \mathbf{w}`$, $`\|\mathbf{v}\|`$, $`\|\mathbf{w}\|`$, and the angle $`\theta`$ between the vectors, and check the Schwarz and triangle inequalities.

   <details>
   <summary>Answer</summary>

   $`\mathbf{v} \cdot \mathbf{w} = 1 \cdot 1 + 1 \cdot 0 + 0 \cdot 1 = 1`$ and $`\|\mathbf{v}\| = \|\mathbf{w}\| = \sqrt{2}`$. Since $`\cos\theta = \frac{1}{\sqrt{2}\sqrt{2}} = \frac{1}{2}`$, $`\theta = 60°`$. The Schwarz inequality holds as $`|\mathbf{v} \cdot \mathbf{w}| = 1 \leq \|\mathbf{v}\|\|\mathbf{w}\| = 2`$. Since $`\mathbf{v} + \mathbf{w} = \begin{pmatrix} 2 \\ 1 \\ 1 \end{pmatrix}`$, $`\|\mathbf{v} + \mathbf{w}\| = \sqrt{6} \approx 2.449`$, which is less than $`\|\mathbf{v}\| + \|\mathbf{w}\| = 2\sqrt{2} \approx 2.828`$, so the triangle inequality also holds.

   </details>

3. **Column space and A = CR:** For $`A = \begin{pmatrix} 1 & 2 & 1 \\ 2 & 4 & 3 \\ 3 & 6 & 4 \end{pmatrix}`$, find the independent columns from left to right to get $`C`$ and the rank $`r`$, and write $`R`$ with $`A = CR`$. What shape is $`C(A)`$ in $`\mathbb{R}^3`$?

   <details>
   <summary>Answer</summary>

   $`\mathbf{a}_1 = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix} \neq \mathbf{0}`$, so it goes into $`C`$. $`\mathbf{a}_2 = 2\mathbf{a}_1`$, so it does not. $`\mathbf{a}_3 = \begin{pmatrix} 1 \\ 3 \\ 4 \end{pmatrix}`$ is not a multiple of $`\mathbf{a}_1`$, so it goes into $`C`$. So $`C = \begin{pmatrix} 1 & 1 \\ 2 & 3 \\ 3 & 4 \end{pmatrix}`$, $`r = 2`$, and $`R = \begin{pmatrix} 1 & 2 & 0 \\ 0 & 0 & 1 \end{pmatrix}`$. Indeed, the three columns of $`CR`$ are $`\mathbf{a}_1`$, $`2\mathbf{a}_1`$, and $`\mathbf{a}_3`$. With two independent columns, $`C(A)`$ is a plane through the origin in $`\mathbb{R}^3`$.

   </details>

4. **Matrix multiplication and non-commutativity:** For $`A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}`$ and $`B = \begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix}`$, compute $`AB`$ with the column viewpoint (column $`j`$ of $`AB`$ is $`A\mathbf{b}_j`$) and compare it with $`BA`$. Also give the number of multiplications needed for $`AB`$.

   <details>
   <summary>Answer</summary>

   $`A\mathbf{b}_1 = 3\begin{pmatrix} 1 \\ 0 \end{pmatrix} + 1\begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 5 \\ 1 \end{pmatrix}`$ and $`A\mathbf{b}_2 = 1\begin{pmatrix} 1 \\ 0 \end{pmatrix} + 2\begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 5 \\ 2 \end{pmatrix}`$, so $`AB = \begin{pmatrix} 5 & 5 \\ 1 & 2 \end{pmatrix}`$. On the other hand, $`BA = \begin{pmatrix} 3 & 7 \\ 1 & 4 \end{pmatrix}`$, so $`AB \neq BA`$. Since $`m = n = p = 2`$, the number of multiplications is $`mnp = 8`$.

   </details>

5. **Rank one matrix:** Write $`A = \begin{pmatrix} 2 & 4 & 6 \\ 1 & 2 & 3 \end{pmatrix}`$ as a column times a row, $`A = \mathbf{c}\mathbf{r}^T`$, and describe in which spaces the column space and the row space are lines.

   <details>
   <summary>Answer</summary>

   All three columns are multiples of $`\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$ (1, 2, and 3 times), so $`A = \begin{pmatrix} 2 \\ 1 \end{pmatrix}\begin{pmatrix} 1 & 2 & 3 \end{pmatrix}`$ and the rank is 1. The column space is the line through the origin in $`\mathbb{R}^2`$ along $`\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, and the row space is the line through the origin in $`\mathbb{R}^3`$ along $`(1, 2, 3)`$.

   </details>
