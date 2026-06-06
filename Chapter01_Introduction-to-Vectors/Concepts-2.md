# Chapter 1 Lecture — Vectors and Matrices — Part 2/2

> **📑 This document is split into 2 parts.** GitHub renders only a limited number of math expressions per file, so it is split by section so that all math displays correctly.
>
> [Part 1](Concepts.md) · **Part 2**

---

<br>

## Table of Contents

- [1.3 Matrices and Their Column Spaces](#13-matrices-and-their-column-spaces)
  - [1.3.1 Matrix-Vector Multiplication](#131-matrix-vector-multiplication)
  - [1.3.2 Column Space](#132-column-space)
  - [1.3.3 Independent, Dependent, and Column Space](#133-independent-dependent-and-column-space)
  - [1.3.4 Span](#134-span)
  - [1.3.5 Rank and Basis](#135-rank-and-basis)
  - [1.3.6 Matrices of Rank One](#136-matrices-of-rank-one)
- [1.4 Matrix Multiplication AB and CR](#14-matrix-multiplication-ab-and-cr)
  - [1.4.1 Rules for Matrix Multiplication](#141-rules-for-matrix-multiplication)
  - [1.4.2 Column Interpretation of AB](#142-column-interpretation-of-ab)
  - [1.4.3 Computational Cost](#143-computational-cost)
  - [1.4.4 Properties of Matrix Multiplication](#144-properties-of-matrix-multiplication)
  - [1.4.5 Rank One Matrices and A = CR](#145-rank-one-matrices-and-a--cr)
  - [1.4.6 Finding C and R](#146-finding-c-and-r)
  - [1.4.7 Columns of A times Rows of B (Outer Product)](#147-columns-of-a-times-rows-of-b-outer-product)
- [Summary](#summary)

---

<br>

## 1.3 Matrices and Their Column Spaces

### 1.3.1 Matrix-Vector Multiplication

**(1)** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}`$ is a **3 by 2 matrix**: 3 rows and 2 columns, Rank 2.

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

**(5)** Rank one matrices: All columns of $`A`$ are on ONE line.

### Consider $`m`$ by $`n`$ Matrices

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}
```

If $`m = n`$: **square matrix**.

**Examples of square matrices ($`A \in \mathbb{R}^{3 \times 3}`$):**

```math
\text{Identity matrix: } \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \qquad \text{Diagonal matrix: } \begin{pmatrix} 2 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 0 & 5 \end{pmatrix}
```

```math
\text{Triangular matrix: } \begin{pmatrix} 2 & 1 & -3 \\ 0 & 4 & 7 \\ 0 & 0 & 5 \end{pmatrix} \qquad \text{Symmetric matrix: } \begin{pmatrix} 2 & 1 & -3 \\ 1 & 4 & 7 \\ -3 & 7 & 5 \end{pmatrix}
```

### Interpreting Columns as Vectors

We can interpret the columns of $`A`$ as vectors:

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}
```

where $`\mathbf{a}_i \in \mathbb{R}^m`$, $`i = 1, 2, \ldots, n`$.

**Q.** How does an $`m`$ by $`n`$ matrix $`A`$ multiply an $`n`$ by 1 vector $`\mathbf{x}`$?

1. **Dot products** of $`\mathbf{x}`$ with the rows of $`A`$
2. **Linear combinations** of the columns of $`A`$

### Two Viewpoints of $`A\mathbf{x}`$

**(1) Row viewpoint (dot products):**

```math
A\mathbf{x} = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}
```

3 rows $`\Rightarrow`$ 3 dot products.

Note that each row of $`A`$ has the same number of components as the vector $`\mathbf{x}`$.

**(2) Column viewpoint (linear combinations):**

```math
A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n
```

A **combination of the columns of $`A`$**.

**Example:** $`A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}`$

```math
A\mathbf{x} = x_1\begin{pmatrix} -1 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix} + x_4\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}
```

### 1.3.2 Column Space

When $`A = (\mathbf{a}_1 \ \mathbf{a}_2)`$ with linearly independent $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$:

```math
A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2
```

shows **the plane of all combinations** (the span).

### 1.3.3 Independent, Dependent, and Column Space

**Example (ex1):**

```math
A_1 = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 4 & 0 \\ 3 & 5 & 6 \end{pmatrix}
```

Each column gives a new direction. Their combinations fill **3D space**.

```math
A_1\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = x_1\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix} + x_2\begin{pmatrix} 0 \\ 4 \\ 5 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 0 \\ 6 \end{pmatrix}
```

**Independent columns:** $`A_1\mathbf{x} = \mathbf{0}`$ iff $`\mathbf{x} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$

**Example (ex2):**

```math
A_2 = \begin{pmatrix} 1 & 2 & 3 \\ 1 & 4 & 5 \\ 6 & 0 & 6 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

$`\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2`$. The combinations of $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3`$ **DO NOT** fill 3D space. $`\mathbf{a}_3`$ lies on the plane of $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ in 3D.

**Example (ex3):**

```math
A_3 = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 6 & 8 \\ 5 & 15 & 20 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

$`\mathbf{a}_2 = 3\mathbf{a}_1`$ and $`\mathbf{a}_3 = 4\mathbf{a}_1`$.

All three columns of $`A_3`$ lie on the **same line** in 3D space.

### The Column Space $`C(A)`$

The column space $`C(A)`$ contains all vectors $`A\mathbf{x}`$: All combinations of the columns.

**Thinking about the column space of $`A`$:**

```math
A_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix} \quad \text{(Independent columns)}
```

```math
A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix} \quad \text{Is the 4th column dependent?}
```

For $`A_5`$: $`\mathbf{a}_4 = \mathbf{a}_1 - \mathbf{a}_2 + \mathbf{a}_3`$.

Consider $`\mathbf{v} = A_4\mathbf{x}`$:

```math
\mathbf{v} = x_1\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + x_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
```

Solving $`A_4\mathbf{x} = \mathbf{v}`$: From the equations,

```math
v_4 = x_4, \quad v_3 = x_4 + x_3, \quad v_2 = x_4 + x_3 + x_2, \quad v_1 = x_4 + x_3 + x_2 + x_1
```

So: $`x_4 = v_4`$, $`x_3 = v_3 - v_4`$, $`x_2 = v_2 - v_3`$, $`x_1 = v_1 - v_2`$

```math
\therefore \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \\ v_4 \end{pmatrix} = (v_1 - v_2)\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + (v_2 - v_3)\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + (v_3 - v_4)\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + v_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
```

This implies that **every $`\mathbf{v}`$ is in the column space**. This solved the four equations $`A_4\mathbf{x} = \mathbf{v}`$.

### 1.3.4 Span

**SPAN** describes all the linear combinations of a set of vectors.

**The span of the columns of $`A`$ is the column space.**

**Example:** $`\mathbf{a}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}`$, $`\mathbf{a}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} \in \mathbb{R}^4`$

- $`c_1\mathbf{a}_1`$ corresponds to a line in 4D.
- $`c_2\mathbf{a}_2`$ corresponds to a line in 4D.
- $`c_1\mathbf{a}_1 + c_2\mathbf{a}_2`$ fill a 2D plane in 4D space.
- The plane is the span of columns $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$.

**Example:** $`A_4`$ has 4 independent columns. The column space of $`A_4`$ is all of $`\mathbb{R}^4`$.

**Example:** $`A_5`$ has one dependent column (only three independent columns).

```math
A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix}
```

The column space of $`A_5`$ is a **3D subspace** inside $`\mathbb{R}^4`$. The 4th column is in that subspace.

We can only solve $`A_5\mathbf{x} = \mathbf{v}`$ when $`\mathbf{v} \in C(A_5)`$.

### 1.3.5 Rank and Basis

Let $`A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}`$

$`C(A)`$ is the column space of $`A`$, $`\mathbf{a}_i \in \mathbb{R}^m`$.

The column space $`C(A)`$ might fill all of $`\mathbb{R}^m`$ or might not.

**Example:** Take $`m = 3`$. $`A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}`$

$`C(A)`$ is:
- The whole space $`\mathbb{R}^3`$ — if $`A`$ has 3 independent columns (I.C.)
- A plane in $`\mathbb{R}^3`$ — if $`A`$ has 2 I.C.
- A line in $`\mathbb{R}^3`$ — if $`A`$ has 1 I.C.
- A single point $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$ — if $`A`$ is a zero matrix

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

- **Q1:** How many columns of $`A`$ are independent? That number $`r`$ is the **"rank"** of $`A`$.
- **Q2:** Which are the first $`r`$ independent columns? They are a **"basis"** for the column space.
- **Q3:** What combinations of those $`r`$ basis vectors produce the remaining $`n - r`$ columns?
- **Q4:** Write any $`A`$ as an $`m \times r`$ column matrix $`C`$ times an $`r \times n`$ matrix $`R`$: $`A = CR`$.
- **Q5:** The $`r`$ rows of $`R`$ are a basis for the **row space** of $`A`$. The rows of $`R`$ DO NOT come directly from $`A`$.

### 1.3.6 Matrices of Rank One

For a rank one matrix, **all column vectors lie along the same line**.

**Example:**

```math
A_6 = \begin{pmatrix} 1 & 3 & -2 \\ 4 & 12 & -8 \\ 2 & 6 & -4 \end{pmatrix}
```

$`\mathbf{a}_2 = 3\mathbf{a}_1`$, $`\mathbf{a}_3 = -2\mathbf{a}_1`$.

- $`A_6`$ has rank $`r = 1`$.
- $`C(A_6) = c_1\mathbf{a}_1`$, a line passing through the origin.
- All rows of $`A_6`$ are multiples of one row.

> When the column space is a single line $`\in \mathbb{R}^m`$, the row space is a single line $`\in \mathbb{R}^n`$.

**Q: Why does this happen?**

If all columns are in the same direction, then all rows are in the same direction.

**Example:** $`A = \begin{pmatrix} a & ma \\ b & mb \end{pmatrix}`$

Row 2 is a multiple of row 1. If the column rank is 1, then the row rank is 1.

**Example:** $`A = \begin{pmatrix} a & ma & pa \\ b & mb & pb \\ c & mc & pc \end{pmatrix}`$

Row 2 = $`\frac{b}{a}`$ row 1, Row 3 = $`\frac{c}{a}`$ row 1.

If the column rank is 1, then the row rank is 1.

> **Q. Row rank equals column rank for every matrix.** Answer: **Yes.**

---

<br>

## 1.4 Matrix Multiplication AB and CR

### 1.4.1 Rules for Matrix Multiplication

**(1)** To multiply $`AB`$: $`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$.

**Row length for $`A`$ should be equal to column length for $`B`$.**

**(2) Dot product viewpoint:**

```math
(AB)_{ij} = (\text{row } i \text{ of } A) \cdot (\text{column } j \text{ of } B)
```

**(3) Column viewpoint:**

```math
AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p)
```

Column $`j`$ of $`AB`$ = $`A\mathbf{b}_j`$

**(4) Non-commutativity:** $`AB \neq BA`$ in general.

**(5)** If $`A`$ has $`r`$ independent columns in $`C`$: $`A = CR`$ where $`C \in \mathbb{R}^{m \times r}`$, $`R \in \mathbb{R}^{r \times n}`$.

### 1.4.2 Column Interpretation of AB

Let $`B = \begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix}`$

```math
AB = A\begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix} = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}
```

These are **combinations of $`A`$**.

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

Both approaches use "4" multiplications.

```math
A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 6 \\ 8 \end{pmatrix} = \begin{pmatrix} 6 + 16 \\ 18 + 32 \end{pmatrix} = \begin{pmatrix} 22 \\ 50 \end{pmatrix}
```

```math
\therefore AB = \begin{pmatrix} 19 & 22 \\ 43 & 50 \end{pmatrix}
```

### 1.4.3 Computational Cost

Let $`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$.

```math
AB = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}
```

**(1) Dot product count:** $`A\mathbf{b}_1 \Rightarrow`$ $`m`$ dot products. $`AB \Rightarrow`$ $`mp`$ dot products. Each dot has $`n`$ multiplications $`\Rightarrow`$ **$`mnp`$ multiplications** total.

**(2) Column combination count:** $`A\mathbf{b}_1 = b_{11}\mathbf{a}_1 + b_{12}\mathbf{a}_2 + \cdots + b_{1n}\mathbf{a}_n \Rightarrow`$ $`m`$ multiplications each, $`mn`$ multiplications for one column. $`AB \Rightarrow`$ **$`mnp`$ multiplications**.

> The cost for matrix multiplication is $`mnp`$.
>
> If $`A, B \in \mathbb{R}^{n \times n}`$, then the cost of $`AB`$ is $`n^3`$.

### 1.4.4 Properties of Matrix Multiplication

**Non-commutativity:**

```math
AB \neq BA \quad \text{in general.}
```

Matrix multiplication is **NOT commutative**.

**Associative law:**

```math
(AB)C = A(BC)
```

For matrix multiplication, the associative law **is** true.

**Distributive law:**

```math
A(B + C) = AB + AC
```

### 1.4.5 Rank One Matrices and A = CR

**Review of $`AB`$:**

**(1) Dot product:** $`(AB)_{ij} = (\text{row } i \text{ of } A) \cdot (\text{column } j \text{ of } B)`$. $`(m \times n)(n \times p) \Rightarrow mp`$ dot products.

**(2) Combine columns:** $`AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p) = (A\mathbf{b}_1 \ A\mathbf{b}_2 \ \cdots \ A\mathbf{b}_p)`$. Column $`j`$ of $`AB`$ = $`A\mathbf{b}_j`$.

**Rank one matrices and $`A = CR`$:**

All columns of a rank one matrix lie on the same line.

When all the columns of $`A`$ are in the same column direction, then all the rows of $`A`$ are in the same row direction.

**Example:** Rank $`r = 1`$.

```math
A = \begin{pmatrix} 1 & 2 & 10 & 100 \\ 3 & 6 & 30 & 300 \\ 2 & 4 & 20 & 200 \end{pmatrix}
```

One independent column, one independent row.

If the column space of $`A`$ is a line, the row space of $`A`$ is also a line.

We can decompose $`A`$ into:

```math
A = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}\begin{pmatrix} 1 & 2 & 10 & 100 \end{pmatrix} = CR
```

where $`C \in \mathbb{R}^{3 \times 1}`$ and $`R \in \mathbb{R}^{1 \times 4}`$.

### 1.4.6 Finding C and R

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

$`\Rightarrow C\mathbf{x} = \mathbf{0}`$ iff $`\mathbf{x} = \mathbf{0}`$. No combinations of columns gives the zero vector.

**Example:**

```math
A = \begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

1. $`C = \begin{pmatrix} 2 \\ 4 \\ 1 \end{pmatrix}`$
2. $`\mathbf{a}_2 = 3\mathbf{a}_1`$
3. $`\mathbf{a}_3 \neq c_1\mathbf{a}_1`$, so $`C = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}`$

$`\therefore`$ Rank $`r = 2`$.

**Q. What is $`R`$ for $`A = CR`$?**

```math
\begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 3 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
3 \times 3 = (3 \times 2)(2 \times 3)
```

Rearranging columns:

```math
\begin{pmatrix} 2 & 4 & 6 \\ 4 & 8 & 12 \\ 1 & 5 & 3 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 0 \end{pmatrix}
```

**We observe that:**

1. $`C`$ contains a full set of $`r`$ independent columns of $`A`$.
2. $`R = (I \ F)`$ contains the identity matrix $`I \in \mathbb{R}^{r \times r}`$.
3. The dependent columns of $`A`$ are combinations of the independent columns in $`C`$.
4. $`A = CR = C(I \ F) = (C \ CF)`$
5. $`C`$ has the same column space of $`A`$. $`R`$ has the same row space of $`A`$.

**Example (ex9):**

```math
\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 4 & 5 \\ 7 & 8 \end{pmatrix}\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}
```

Rank $`r = 2`$ matrix.

$`\mathbf{a}_j = C\mathbf{r}_j`$

$`(1 \ 2 \ 3) = (1 \ 2)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`(4 \ 5 \ 6) = (4 \ 5)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`(7 \ 8 \ 9) = (7 \ 8)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

Row $`i`$ of $`A`$ = row $`i`$ of $`C`$ times $`R`$: a combination of the rows of $`R`$.

**Example:** $`A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 1 & 2 & 4 & 5 \end{pmatrix}_{2 \times 4}`$

```math
= \begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}\begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}
```

$`C \in \mathbb{R}^{2 \times 2}`$, $`R \in \mathbb{R}^{2 \times 4}`$. The rank $`r = 2`$, we recover all the columns of $`A`$ from $`C`$ by using $`R`$.

### How to Find the Matrix R

We can use **"elimination"**, which will be covered in Chapter 3.

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

Columns of $`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ are independent.

```math
A = \begin{pmatrix} 1 & 3 \\ 2 & 4 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \end{pmatrix} = CR
```

### Key Properties of A = CR

- $`A = CR`$
- $`r`$ columns of $`C`$ are a **basis** for the column space of $`A`$.
- $`r`$ rows of $`R`$ are a **basis** for the row space of $`A`$.
- $`\Rightarrow`$ $`r`$ dimension.

> In $`A = CR`$, "$`R`$" is the **reduced row echelon form**.

**Example (ex5):**

```math
A = \begin{pmatrix} | & | & | \\ \mathbf{a}_1 & \mathbf{a}_2 & 3\mathbf{a}_1 + 4\mathbf{a}_2 \\ | & | & | \end{pmatrix}_{m \times 3}
```

```math
= \begin{pmatrix} | & | \\ \mathbf{a}_1 & \mathbf{a}_2 \\ | & | \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 4 \end{pmatrix}
```

```math
= \mathbf{a}_1(1 \ 0 \ 3) + \mathbf{a}_2(0 \ 1 \ 4) = (\mathbf{a}_1 \ \mathbf{a}_2 \ 3\mathbf{a}_1 + 4\mathbf{a}_2) = A
```

### 1.4.7 Columns of A times Rows of B (Outer Product)

```math
AB = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}_{m \times n} \begin{pmatrix} - & \mathbf{b}_1^* & - \\ - & \mathbf{b}_2^* & - \\ & \vdots & \\ - & \mathbf{b}_n^* & - \end{pmatrix}_{n \times p}
```

```math
= \mathbf{a}_1\mathbf{b}_1^* + \mathbf{a}_2\mathbf{b}_2^* + \cdots + \mathbf{a}_n\mathbf{b}_n^*
```

where $`\mathbf{a}_k\mathbf{b}_k^* = \begin{pmatrix} | \\ \mathbf{a}_k \\ | \end{pmatrix}_{m \times 1}(- \ \mathbf{b}_k^* \ -)_{1 \times p}`$

This is a **column times row = rank 1 matrix**, with $`mp`$ entries.

```math
AB = \sum_{k=1}^{n} \mathbf{a}_k\mathbf{b}_k^*
```

The sum of $`n`$ rank one matrices. Each has $`mp`$ multiplications $`\Rightarrow`$ total $`mnp`$ multiplications.

**Remark:** For the matrix $`A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}`$, $`\text{rank}(A) = 2`$.

$`\Rightarrow`$ $`A`$ has no inverse matrix.
$`\Rightarrow`$ The determinant of $`A`$ is zero.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| **Linear Combination** | $`c\mathbf{v} + d\mathbf{w}`$: scalar multiplication + vector addition |
| **Column Space $`C(A)`$** | All vectors $`A\mathbf{x}`$ = all combinations of columns of $`A`$ |
| **Span** | All linear combinations of a set of vectors |
| **Independent Vectors** | $`A\mathbf{x} = \mathbf{0}`$ only when $`\mathbf{x} = \mathbf{0}`$ |
| **Dependent Vector** | Can be written as a combination of other vectors |
| **Rank** | Number of independent columns (= number of independent rows) |
| **Basis** | The first $`r`$ independent columns form a basis for $`C(A)`$ |
| **Dot Product** | $`\mathbf{v} \cdot \mathbf{w} = \sum v_i w_i`$; gives length and angle info |
| **Length (Norm)** | $`\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}}`$ |
| **Unit Vector** | $`\|\mathbf{v}\| = 1`$; any nonzero vector normalized: $`\mathbf{v}/\|\mathbf{v}\|`$ |
| **Perpendicularity** | $`\mathbf{v} \perp \mathbf{w} \Leftrightarrow \mathbf{v} \cdot \mathbf{w} = 0`$ |
| **Angle Formula** | $`\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|}`$ |
| **Schwarz Inequality** | $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$ |
| **Triangle Inequality** | $`\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|`$ |
| **Pythagorean Theorem** | If $`\mathbf{v} \perp \mathbf{w}`$: $`\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2`$ |
| **Matrix Multiplication $`AB`$** | $`(AB)_{ij} = \text{row}_i(A) \cdot \text{col}_j(B)`$; column $`j`$ of $`AB = A\mathbf{b}_j`$ |
| **Outer Product** | $`AB = \sum \mathbf{a}_k \mathbf{b}_k^*`$ (sum of rank-1 matrices) |
| **Cost of $`AB`$** | $`mnp`$ multiplications for $`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$ |
| **Non-commutativity** | $`AB \neq BA`$ in general |
| **Associativity** | $`(AB)C = A(BC)`$ |
| **Distributivity** | $`A(B+C) = AB + AC`$ |
| **$`A = CR`$ Factorization** | $`C`$: independent columns of $`A`$; $`R`$: reduced row echelon form |
| **Rank One Matrix** | All columns on one line; $`A = \mathbf{c}\mathbf{r}^T`$ |
| **Row Rank = Column Rank** | Always true for any matrix |
| **Identity Matrix** | $`\mathbf{v} = I\mathbf{v}`$; standard basis vectors as columns |
| **Plane via Normal** | $`\mathbf{w} \cdot \mathbf{n} = 0`$ defines a plane perpendicular to $`\mathbf{n}`$ |
| **Invertibility** | $`A`$ invertible $`\Leftrightarrow`$ columns are independent $`\Leftrightarrow`$ unique solution to $`A\mathbf{x} = \mathbf{b}`$ |
| **Elimination Failure** | When $`\mathbf{v} \| \mathbf{w}`$: no solution or infinite solutions |

---
