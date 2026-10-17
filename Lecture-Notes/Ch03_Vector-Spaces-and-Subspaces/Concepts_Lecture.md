# Chapter 3 Lecture — Vector Spaces and Subspaces

> **Last Updated:** 2026-06-19
>
> Introduction to Linear Algebra, Strang (6th Ed.) - Ch 3

> **Prerequisites**: [Linear Algebra] Elimination and matrix operations (Ch 1-2).
>
> **Learning Objectives**:
> 1. Define vector spaces, subspaces, and verify closure properties
> 2. Compute column space, null space, row space, and left null space
> 3. Determine basis, dimension, and rank of a matrix

> **📑 This document is split into 2 parts.**
>
> **Part 1** · [Part 2](Concepts-Part2.md)

---

<br>

## Table of Contents

- [1. Overview: The Four Fundamental Subspaces](#1-overview-the-four-fundamental-subspaces)
- [2. Vector Spaces and Subspaces (3.1)](#2-vector-spaces-and-subspaces-31)
  - [2.1 Requirements for a Vector Space](#21-requirements-for-a-vector-space)
  - [2.2 The Space R^n](#22-the-space-rn)
  - [2.3 Definition of a Vector Space](#23-definition-of-a-vector-space)
  - [2.4 The Eight Axioms](#24-the-eight-axioms)
  - [2.5 Consequences of the Axioms](#25-consequences-of-the-axioms)
  - [2.6 Examples: Is It a Vector Space?](#26-examples-is-it-a-vector-space)
  - [2.7 Examples of Vector Spaces](#27-examples-of-vector-spaces)
  - [2.8 Generalized Vector Spaces](#28-generalized-vector-spaces)
  - [2.9 Subspaces of Vector Spaces](#29-subspaces-of-vector-spaces)
  - [2.10 Definition of a Subspace](#210-definition-of-a-subspace)
  - [2.11 Examples of Subspaces and Non-Subspaces](#211-examples-of-subspaces-and-non-subspaces)
  - [2.12 Column Space and Row Space](#212-column-space-and-row-space)
  - [2.13 Spanning](#213-spanning)
- [3. Computing the Nullspace by Elimination (3.2)](#3-computing-the-nullspace-by-elimination-32)
  - [3.1 Key Facts: A = CR](#31-key-facts-a--cr)
  - [3.2 Finding All Solutions to Ax = 0](#32-finding-all-solutions-to-ax--0)
  - [3.3 Reduced Row Echelon Form](#33-reduced-row-echelon-form)
  - [3.4 Special Solutions and the Nullspace Basis](#34-special-solutions-and-the-nullspace-basis)
  - [3.5 The Nullspace Matrix: Columns of (-F; I)](#35-the-nullspace-matrix-columns-of--fi)
  - [3.6 The Matrix Factorization A = CR and N(A)](#36-the-matrix-factorization-a--cr-and-na)
- [4. The Complete Solution to Ax = b (3.3)](#4-the-complete-solution-to-ax--b-33)
  - [4.1 Structure of the Complete Solution](#41-structure-of-the-complete-solution)
  - [4.2 Worked Example: Finding Particular Solutions](#42-worked-example-finding-particular-solutions)
  - [4.3 The Complete Solution Decomposition](#43-the-complete-solution-decomposition)
  - [4.4 Full Column Rank: r = n](#44-full-column-rank-r--n)
  - [4.5 Solvability Conditions](#45-solvability-conditions)
  - [4.6 Full Row Rank and the Complete Solution](#46-full-row-rank-and-the-complete-solution)
  - [4.7 Four Possibilities for Linear Equations](#47-four-possibilities-for-linear-equations)
- [5. Independence, Basis, and Dimension (3.4)](Concepts-Part2.md#5-independence-basis-and-dimension-34)
  - [5.1 Independent Vectors](Concepts-Part2.md#51-independent-vectors)
  - [5.2 Linear Independence via the Nullspace](Concepts-Part2.md#52-linear-independence-via-the-nullspace)
  - [5.3 Vectors that Span a Subspace](Concepts-Part2.md#53-vectors-that-span-a-subspace)
  - [5.4 Basis for a Vector Space](Concepts-Part2.md#54-basis-for-a-vector-space)
  - [5.5 Dimension of a Vector Space](Concepts-Part2.md#55-dimension-of-a-vector-space)
  - [5.6 Bases for Matrix Spaces and Function Spaces](Concepts-Part2.md#56-bases-for-matrix-spaces-and-function-spaces)
- [6. Dimensions of the Four Subspaces (3.5)](Concepts-Part2.md#6-dimensions-of-the-four-subspaces-35)
  - [6.1 Dimension Summary](Concepts-Part2.md#61-dimension-summary)
  - [6.2 Orthogonality of the Subspaces](Concepts-Part2.md#62-orthogonality-of-the-subspaces)
  - [6.3 The Four Subspaces for R_0](Concepts-Part2.md#63-the-four-subspaces-for-r_0)
  - [6.4 Relationship Between A and R_0](Concepts-Part2.md#64-relationship-between-a-and-r_0)
  - [6.5 The Fundamental Theorem of Linear Algebra](Concepts-Part2.md#65-the-fundamental-theorem-of-linear-algebra)
- [Summary](Concepts-Part2.md#summary)

---

<br>

## 1. Overview: The Four Fundamental Subspaces

Chapter 3 covers five major topics:

- **3.1** Vector Spaces and Subspaces — How do we define a vector space? Key operations are $`\mathbf{u} + \mathbf{v}`$ and $`c\mathbf{u}`$. There are 8 rules that vectors $`\mathbf{u}`$ and scalar $`c`$ must satisfy.
- **3.2** The Nullspace of $`A`$: Solving $`A\mathbf{x} = \mathbf{0}`$
- **3.3** The Complete Solution to $`A\mathbf{x} = \mathbf{b}`$
- **3.4** Independence, **Basis**, and **Dimension** — A set of vectors that describes the space. Let $`A \in \mathbb{R}^{n \times n}`$. $`A`$ has $`r`$ independent columns $`\Rightarrow`$ $`C(A)`$ has dimension $`r`$. $`(n - r)`$ special solutions to $`A\mathbf{x} = \mathbf{0}`$ are a basis for $`N(A)`$, the null space of $`A`$.
- **3.5** Dimensions of the Four Subspaces:

| Subspace | Dimension |
|:---------|:----------|
| Column space of $`A`$ | $`r`$ |
| Row space of $`A`$ | $`r`$ |
| Null space of $`A`$ | $`n - r`$ |
| Null space of $`A^T`$ (left nullspace) | $`m - r`$ |

---

<br>

## 2. Vector Spaces and Subspaces (3.1)

### 2.1 Requirements for a Vector Space

1. **Requirement:** All the linear combinations $`c\mathbf{u} + d\mathbf{w}`$ must stay in the vector space.
2. The **row space** of $`A`$ is "spanned" by the rows of $`A`$. The columns of $`A`$ span the column space $`C(A)`$.
3. **Matrices** $`M_1`$ to $`M_n`$ and **functions** $`f_1`$ to $`f_n`$ span **matrix spaces** and **function spaces**.

### 2.2 The Space R^n

The space $`\mathbb{R}^n`$ contains all column vectors $`\mathbf{v}`$ of length $`n`$.

```math
\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} \in \mathbb{R}^n \quad \text{where } x_1, x_2, \ldots, x_n \in \mathbb{R} \text{ (real numbers)}
```

Examples:
- $`x \in \mathbb{R}^1`$
- $`\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} \in \mathbb{R}^2`$
- $`\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} \in \mathbb{R}^3`$

**Note:** When $`x_1, x_2, \ldots, x_n \in \mathbb{C}`$ (complex numbers), the space becomes $`\mathbb{C}^n`$.

**Note:** If $`c, d \in \mathbb{R}`$ and $`\mathbf{u}, \mathbf{w} \in \mathbb{R}^n`$, then $`c\mathbf{u} + d\mathbf{w} \in \mathbb{R}^n`$. The linear combinations stay in the vector space $`\mathbb{R}^n`$.

### 2.3 Definition of a Vector Space

**Vector space** (= linear space) $`V`$ is a set whose elements (vectors) can be added together and multiplied by numbers.

Let $`\mathbf{u}, \mathbf{w} \in V`$ (vector space):
1. $`\mathbf{u} + \mathbf{w} \in V`$ — vector addition
2. $`\alpha \mathbf{u} \in V \quad \forall \alpha \in \mathbb{F}`$ — scalar multiplication

**A field** is a set on which addition, subtraction, multiplication, and division are defined.

Examples of fields: $`\mathbb{R}`$ (real numbers), $`\mathbb{C}`$ (complex numbers).

### 2.4 The Eight Axioms

Let $`\mathbf{u}, \mathbf{v}, \mathbf{w} \in V`$ and $`c, d \in \mathbb{F}`$. A vector space must satisfy the following eight axioms:

**(1) Associativity of vector addition:**

```math
\mathbf{u} + (\mathbf{v} + \mathbf{w}) = (\mathbf{u} + \mathbf{v}) + \mathbf{w}
```

**(2) Commutativity of vector addition:**

```math
\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}
```

**(3) Identity element of vector addition (zero vector):**

```math
\exists!\; \mathbf{0} \in V \text{ s.t. } \mathbf{v} + \mathbf{0} = \mathbf{v} \quad \forall \mathbf{v} \in V
```

**(4) Inverse element of vector addition:**

```math
\forall \mathbf{v} \in V, \; \exists!\; -\mathbf{v} \in V \text{ s.t. } \mathbf{v} + (-\mathbf{v}) = \mathbf{0}
```

($`-\mathbf{v}`$ is the additive inverse of $`\mathbf{v}`$)

**(5) Compatibility of scalar multiplication with field multiplication:**

```math
c(d\mathbf{v}) = (cd)\mathbf{v}
```

**(6) Identity element of scalar multiplication:**

```math
1 \cdot \mathbf{v} = \mathbf{v}
```

**(7) Distributivity of scalar multiplication w.r.t. vector addition:**

```math
c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}
```

**(8) Distributivity of scalar multiplication w.r.t. field addition:**

```math
(c + d)\mathbf{u} = c\mathbf{u} + d\mathbf{u}
```

### 2.5 Consequences of the Axioms

From the eight axioms, the following properties hold:

```math
0\mathbf{u} = \mathbf{0}
```

```math
c\mathbf{0} = \mathbf{0}
```

```math
(-1)\mathbf{u} = -\mathbf{u}
```

```math
c\mathbf{v} = \mathbf{0} \implies c = 0 \text{ or } \mathbf{v} = \mathbf{0}
```

### 2.6 Examples: Is It a Vector Space?

**Q.** Is the set of all positive vectors $`\mathbb{X}`$, where $`\mathbb{X} \ni \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}`$ with every $`v_i > 0`$, a vector space?

**A.** No. $`-\mathbf{v} \notin \mathbb{X}`$.

---

**Q.** Let $`\mathbb{X}`$ be the set of solutions to $`A\mathbf{x} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$. Is $`\mathbb{X}`$ a vector space?

**A.** No. If $`\mathbf{u}, \mathbf{w} \in \mathbb{X}`$, then $`A\mathbf{u} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$ and $`A\mathbf{w} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$. But $`\mathbf{u} + \mathbf{w} \notin \mathbb{X}`$ because $`A(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w} = \begin{pmatrix} 2 \\ \vdots \\ 2 \end{pmatrix} \neq \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$.

---

**Q.** Is a line in $`\mathbb{R}^n`$ a vector space?

A line is a collection of points. $`\mathbf{q}`$ represents a position on the line:

```math
\mathbf{q} = \mathbf{p} + t\mathbf{d}, \quad t \in \mathbb{R}
```

$`\mathbf{p}`$ and $`\mathbf{d}`$ are fixed. $`t`$ varies from $`-\infty`$ to $`\infty`$.

Let $`\mathbb{X}`$ be the set of all vectors $`\mathbf{q}`$: $`\mathbb{X} = \lbrace\mathbf{q} : \mathbf{p} + t\mathbf{d} \;\forall\; t \in \mathbb{R}\rbrace`$.

Take $`\mathbf{a}, \mathbf{b} \in \mathbb{X}`$. Does $`\mathbf{a} + \mathbf{b}`$ belong to $`\mathbb{X}`$?

```math
\mathbf{a} = \mathbf{p} + t_1 \mathbf{d}
```
```math
\mathbf{b} = \mathbf{p} + t_2 \mathbf{d}
```
```math
\mathbf{a} + \mathbf{b} = 2\mathbf{p} + (t_1 + t_2)\mathbf{d} \notin \mathbb{X} \quad \text{if } \mathbf{p} \neq \mathbf{0}
```

However:
- $`\mathbf{a} + \mathbf{b} = (t_1 + t_2)\mathbf{d} \in \mathbb{X}`$ if $`\mathbf{p} = \mathbf{0}`$
- $`c\mathbf{a} = ct_1\mathbf{d} \in \mathbb{X}`$ if $`\mathbf{p} = \mathbf{0}`$

**A line passing through the origin is a vector space.** The line through $`\mathbf{0}`$ in $`\mathbb{R}^n`$ is a **subspace** of $`\mathbb{R}^n`$ — a vector space inside another vector space.

### 2.7 Examples of Vector Spaces

- $`\mathbb{R}^n`$ is a vector space.
- $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$ is a vector space (the smallest vector space):
  - i) $`\mathbf{0} + \mathbf{0} = \mathbf{0} \in \mathbb{Z}`$
  - ii) $`c\mathbf{0} = \mathbf{0} \in \mathbb{Z}`$

$`\mathbb{Z}`$ is the **null space of an invertible matrix**. If the only solution to $`A\mathbf{x} = \mathbf{0}`$ is the zero vector $`\mathbf{x} = \mathbf{0}`$, then the columns of $`A`$ are linearly independent (LI) and the nullspace of $`A`$ is $`\mathbb{Z}`$.

### 2.8 Generalized Vector Spaces

**Note:** A vector space (= a linear space) is a set whose elements can be added together and multiplied by numbers.

i.e., $`\mathbf{u}, \mathbf{w} \in V \implies c\mathbf{u} + d\mathbf{w} \in V`$.

**Matrices and functions can be considered as vectors.**

**Vector space of matrices:**

$`A, B \in \mathbb{R}^{m \times n} \implies cA + dB \in \mathbb{R}^{m \times n}`$

The set of all matrices of a fixed size forms a vector space. (Check if they satisfy the eight rules.)

**Vector space of functions:**

Let $`\mathbb{F}`$ be a set of functions that take elements from $`\mathbb{R}`$ and map them to a real number:

```math
\mathbb{F} = \lbracef \mid f: \mathbb{R} \to \mathbb{R}\rbrace
```

Let $`f, g \in \mathbb{F}`$, $`c, d \in \mathbb{R}`$. Then $`cf + dg \in \mathbb{F}`$.

Define:
- $`(f + g)(x) = f(x) + g(x)`$
- $`(cf)(x) = c(f(x))`$

**Example:** $`\mathbb{F}`$ = the line of functions $`y = ce^x`$.

$`\mathbb{F} = \lbracef \mid f: \mathbb{R} \ni x \to ce^x \in \mathbb{R}, \;\forall c \in \mathbb{R}\rbrace`$

$`f(x) = e^x`$, $`g(x) = 2e^x`$.

$`(f + g)(x) = f(x) + g(x) = e^x + 2e^x = 3e^x \in \mathbb{F}`$

$`(cf)(x) = c(f(x)) = ce^x \in \mathbb{F}`$

**Remark:** A **set** is just a collection of elements without any additional structure. A **space** is a set along with additional structures defined on it. Example: a vector space is a set of vectors where vector addition and scalar multiplication is defined. A vector space should satisfy the 8 rules.

### 2.9 Subspaces of Vector Spaces

Consider a vector space $`\mathbb{R}^n`$, where $`\mathbf{v} \in \mathbb{R}^n`$ is a column vector with $`n`$ components.

There are important vector spaces inside $`\mathbb{R}^n`$. Those are **subspaces** of $`\mathbb{R}^n`$.

A plane is a vector space. If $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$, then $`c\mathbf{v} + d\mathbf{w} \in \mathbb{R}^2`$.

A plane through the origin $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$ in $`\mathbb{R}^3`$ is a vector space. This plane is not $`\mathbb{R}^2`$ because $`\mathbf{u}, \mathbf{v} \in \mathbb{R}^3`$. The plane is a vector space inside $`\mathbb{R}^3`$. This is a **subspace** of the full vector space $`\mathbb{R}^3`$.

### 2.10 Definition of a Subspace

**Def.** A **subspace** of a vector space is a set of vectors (including $`\mathbf{0}`$) that satisfies:

> i) $`\mathbf{u} + \mathbf{w}`$ is in the subspace
>
> ii) $`c\mathbf{u}`$ is in the subspace

if $`\mathbf{u}, \mathbf{w} \in`$ subspace and $`c \in \mathbb{R}`$ (or $`\mathbb{F}`$).

The set of vectors is "**closed**" under addition $`\mathbf{u} + \mathbf{w}`$ and multiplication $`c\mathbf{u}`$.

Conditions i) and ii) mean: **all linear combinations stay in the subspace.**

**Note:** Every subspace contains the zero vector. From ii), $`c = 0`$ implies that $`0\mathbf{u} = \mathbf{0}`$ is in the subspace.

### 2.11 Examples of Subspaces and Non-Subspaces

**Q:** Is the plane $`z = 5`$ in $`\mathbb{R}^3`$ a subspace? **No.** (It does not contain the origin.)

Lines through the origin are subspaces.

$`\mathbb{R}^3`$ is a subspace of itself.

The single vector $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$ is a subspace of $`\mathbb{R}^3`$.

---

**Ex 1.** $`\mathbb{R}^2`$ is a vector space. Is the 1st quadrant a subspace?

```math
\mathcal{U} = \left\lbrace\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0, y \geq 0 \right\rbrace
```

Take $`\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}`$, $`c = -1`$.

Then $`c\mathbf{u} = \begin{pmatrix} -2 \\ -3 \end{pmatrix} \notin \mathcal{U}`$.

This violates rule ii). $`\mathcal{U}`$ is **not** a subspace of $`\mathbb{R}^2`$.

---

**Ex 2.** $`\mathcal{U} = \left\lbrace\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0 \text{ and } y \geq 0, \text{ or } x \leq 0 \text{ and } y \leq 0 \right\rbrace`$

Is $`\mathcal{U}`$ a subspace?

Take $`\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}`$, $`\mathbf{w} = \begin{pmatrix} -3 \\ -2 \end{pmatrix} \in \mathcal{U}`$.

But the sum $`\mathbf{u} + \mathbf{w} = \begin{pmatrix} -1 \\ 1 \end{pmatrix} \notin \mathcal{U}`$.

$`\mathcal{U}`$ is **not** a subspace.

---

**Ex 3.** $`\mathbb{M}`$ is a vector space of $`2 \times 2`$ matrices $`\begin{pmatrix} a & b \\ c & d \end{pmatrix}`$.

$`\mathcal{U}`$ is a set of all upper triangular matrices $`\begin{pmatrix} a & b \\ 0 & d \end{pmatrix}`$.

$`\mathbb{D}`$ is a set of all diagonal matrices $`\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix}`$.

**Both $`\mathcal{U}`$ and $`\mathbb{D}`$ are subspaces of $`\mathbb{M}`$:**

- $`A, B \in \mathcal{U} \implies cA + dB \in \mathcal{U}`$
- $`A, B \in \mathbb{D} \implies cA + dB \in \mathbb{D}`$

Note that:

```math
\begin{pmatrix} a & b \\ c & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + c\begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
```

### 2.12 Column Space and Row Space

**Column Space of $`A`$:**

```math
A\mathbf{x} = \mathbf{b}
```

```math
\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}
```

```math
= x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n
```

The **column space of $`A`$** is a vector space made up of column vectors.

To solve $`A\mathbf{x} = \mathbf{b}`$ is to express $`\mathbf{b}`$ as a linear combination of the columns. The right side vector $`\mathbf{b}`$ has to be in the column space of $`A`$.

> **The equations $`A\mathbf{x} = \mathbf{b}`$ are solvable iff $`\mathbf{b}`$ is in the column space of $`A`$.**

**The Row Space of $`A`$:**

The rows of $`A`$ are the columns of $`A^T`$.

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} \implies A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}
```

$`A^T`$ has $`m`$ column vectors.

**The row space of $`A`$ is the column space of $`A^T`$.**

The equations $`A^T \mathbf{y} = \mathbf{c}`$ are solvable iff $`\mathbf{c}`$ is in the column space of $`A^T`$ ($`= \mathbf{c}`$ is in the row space of $`A`$).

**Example:** Consider a rank 1 matrix $`A = \mathbf{u}\mathbf{v}^T`$:

```math
A = \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix} \begin{pmatrix} v_1 & v_2 & \cdots & v_n \end{pmatrix} = \begin{pmatrix} v_1\mathbf{u} & v_2\mathbf{u} & \cdots & v_n\mathbf{u} \end{pmatrix}
```

$`C(A)`$ is the line of all column vectors of $`A`$: $`c\mathbf{u}`$.

$`A^T = \mathbf{v}\mathbf{u}^T = \begin{pmatrix} u_1\mathbf{v} & u_2\mathbf{v} & \cdots & u_m\mathbf{v} \end{pmatrix}`$

$`C(A^T)`$ is the line of all column vectors of $`A^T`$: $`c\mathbf{v}`$.

### 2.13 Spanning

**The Columns of $`A`$ Span the Vector Space $`C(A)`$.**

Let $`\mathbb{S}`$ be a set of vectors in $`\mathbb{R}^m`$. If $`\mathbb{S} = \left\lbrace \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix}, \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_m \end{pmatrix} \right\rbrace`$, then $`\mathbb{S}`$ is **not** a subspace of $`\mathbb{R}^m`$ because $`\mathbf{u} + \mathbf{v} \notin \mathbb{S}`$.

If we include all combinations of the vectors in $`\mathbb{S}`$, then we have a vector space $`V`$.

**The set $`\mathbb{S}`$ spans $`V`$.** $`V`$ is the smallest vector space containing $`\mathbb{S}`$.

Consider a matrix $`A \in \mathbb{R}^{m \times n}`$:
- $`n`$ columns span the column space $`C(A)`$
- $`m`$ columns of $`A^T`$ span the row space $`C(A^T)`$

---

<br>

## 3. Computing the Nullspace by Elimination (3.2)

### 3.1 Key Facts: A = CR

```math
A = CR
```

1. The nullspace $`N(A)`$ in $`\mathbb{R}^n`$ contains all solutions $`\mathbf{x}`$ to $`A\mathbf{x} = \mathbf{0}`$. This includes $`\mathbf{x} = \mathbf{0}`$.
2. Elimination from $`A`$ to $`R_0`$ to $`R`$ does **not** change the nullspace.
3. The reduced row echelon form $`R_0 = \text{rref}(A)`$ has $`I`$ in $`r`$ columns and $`F`$ in $`n - r`$ columns.
4. If column $`j`$ is dependent on previous columns, $`A\mathbf{x} = \mathbf{0}`$ has a "special solution" with $`x_j = 1`$.
5. The $`n - r`$ special solutions to $`A\mathbf{x} = \mathbf{0}`$ contain $`-F`$ and $`I`$.
6. Every short wide matrix with $`m < n`$ has nonzero solutions to $`A\mathbf{x} = \mathbf{0}`$ in its nullspace.

### 3.2 Finding All Solutions to Ax = 0

We would like to find all solutions to $`A\mathbf{x} = \mathbf{0}`$.

If $`A \in \mathbb{R}^{n \times n}`$ is invertible ($`\text{rank}(A) = n`$), then the only solution is $`\mathbf{x} = \mathbf{0}`$. The nullspace of $`A`$ only contains the zero vector: $`N(A) = \lbrace\mathbf{0}\rbrace`$.

In general, however, $`\text{rank}(A) = r`$. That is, $`A`$ has $`r`$ independent columns. The other $`n - r`$ dependent columns of $`A`$ are combinations of those independent columns. We will find $`n - r`$ vectors in $`N(A)`$ which are special solutions to $`A\mathbf{x} = \mathbf{0}`$.

In Chapter 2, an invertible matrix $`A`$ is reduced to an upper triangular matrix $`U`$. For $`A \in \mathbb{R}^{m \times n}`$, $`A\mathbf{x} = \mathbf{0}`$ is simplified to $`R\mathbf{x} = \mathbf{0}`$ (echelon form).

### 3.3 Reduced Row Echelon Form

In this section, we consider $`A \in \mathbb{R}^{m \times n}`$, and eliminate $`A`$ into reduced row echelon form: $`R_0 = \text{rref}(A)`$.

$`R_0`$ may have zero rows. Removing all zero rows of $`R_0`$ gives $`R`$.

**Ex 1:** $`R = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = (I \quad F)`$

$`\text{rank}(R) = 2`$, $`n = 4`$, $`n - r = 2`$ dependent columns.

```math
R\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_1 + 3x_3 + 5x_4 \\ x_2 + 4x_3 + 6x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}
```

How do we find the special solutions? Fix two variables and use $`R\mathbf{x} = \mathbf{0}`$:

Take $`x_3 = 1, x_4 = 0`$: $`\Rightarrow x_1 + 3 = 0, \; x_2 + 4 = 0`$. $`\therefore \mathbf{s}_1 = \begin{pmatrix} -3 \\ -4 \\ 1 \\ 0 \end{pmatrix}`$

Take $`x_3 = 0, x_4 = 1`$: $`\Rightarrow x_1 + 5 = 0, \; x_2 + 6 = 0`$. $`\therefore \mathbf{s}_2 = \begin{pmatrix} -5 \\ -6 \\ 0 \\ 1 \end{pmatrix}`$

### 3.4 Special Solutions and the Nullspace Basis

The two special solutions $`\mathbf{s}_1, \mathbf{s}_2`$ are in the nullspace of $`R`$:

```math
\mathbf{s}_1, \mathbf{s}_2 \in N(R)
```

```math
R\mathbf{s}_1 = \mathbf{0}, \quad R\mathbf{s}_2 = \mathbf{0}
```

```math
\Rightarrow R(c_1 \mathbf{s}_1 + c_2 \mathbf{s}_2) = \mathbf{0}
```

Linear combinations of $`\mathbf{s}_1`$ and $`\mathbf{s}_2`$ are in the nullspace of $`R`$. **The special solutions $`\mathbf{s}_1`$ and $`\mathbf{s}_2`$ are a basis for the nullspace.**

---

**Ex 2:** $`R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix}`$ (reduced row echelon form with a row of zeros)

```math
R_0 \mathbf{x} = \begin{pmatrix} x_1 + 7x_2 + 8x_4 \\ x_3 + 9x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`\text{rank}(R_0) = 2`$, $`n = 4`$, $`n - r = 2`$. Find two special solutions.

Take $`x_1 = 0, x_2 = 1`$:

$`7 + 8x_4 = 0 \Rightarrow x_4 = -7/8`$, $`x_3 = -9(-7/8) = 63/8`$

$`\Rightarrow \mathbf{x} = \begin{pmatrix} 0 \\ 1 \\ 63/8 \\ -7/8 \end{pmatrix} \Rightarrow \mathbf{s}_1 = \begin{pmatrix} 0 \\ 8 \\ 63 \\ -7 \end{pmatrix}`$

Take $`x_1 = 1, x_2 = 0`$:

$`1 + 8x_4 = 0 \Rightarrow x_4 = -1/8`$, $`x_3 = +9/8`$

$`\Rightarrow \mathbf{x} = \begin{pmatrix} 1 \\ 0 \\ 9/8 \\ -1/8 \end{pmatrix} \Rightarrow \mathbf{s}_2 = \begin{pmatrix} 8 \\ 0 \\ 9 \\ -1 \end{pmatrix}`$

We can also take $`x_2 = 1, x_4 = 0`$ and $`x_2 = 0, x_4 = 1`$ to find $`x_1`$ and $`x_3`$.

### 3.5 The Nullspace Matrix: Columns of (-F; I)

In this section, we consider $`A \in \mathbb{R}^{m \times n}`$, eliminate $`A`$ into reduced row echelon form $`R_0 = \text{rref}(A)`$.

$`R_0`$ may have zero rows. Removing all zero rows of $`R_0`$ gives $`R`$.

```math
R = (I \quad F)Q
```

where $`Q`$ is a permutation matrix (reordering columns so that pivot columns come first).

**Recall** $`R\mathbf{x} = \mathbf{0}`$, i.e., $`(I \quad F)\mathbf{x} = \mathbf{0}`$ (for Ex 1):

```math
\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

Setting free variables to $`\begin{pmatrix} 1 \\ 0 \end{pmatrix}`$ and $`\begin{pmatrix} 0 \\ 1 \end{pmatrix}`$:

```math
I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 3 \\ 4 \end{pmatrix}
```

```math
I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 5 \\ 6 \end{pmatrix}
```

Equivalently:

```math
\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} -3 & -5 \\ -4 & -6 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}
```

That is:

```math
\boxed{(I \quad F) \begin{pmatrix} -F \\ I \end{pmatrix} = O}
```

**The two special solutions to $`(I \quad F)\mathbf{x} = \mathbf{0}`$ are the columns of $`\begin{pmatrix} -F \\ I \end{pmatrix}`$.**

### 3.6 The Matrix Factorization A = CR and N(A)

**Elimination from $`A`$ to $`\text{rref}(A)`$: Reduced Row Echelon Form**

Apply elimination to reduce $`A`$ to $`R_0`$. Then "$`I`$" in $`R_0`$ locates the matrix $`C`$ of independent columns in $`A`$. Removing any zero rows from $`R_0`$ produces the row matrix $`R`$ in $`A = CR`$.

In Chapter 2, $`A`$ was square and invertible:
- $`A\mathbf{x} = \mathbf{b}`$ $`\xrightarrow{a) \text{ elimination}}`$ $`U\mathbf{x} = \mathbf{c}`$ $`\xrightarrow{b) \text{ back substitution}}`$ $`\mathbf{x} = U^{-1}\mathbf{c}`$

For any matrix $`A`$ of rank $`r`$:

```math
A = CR = C(I \quad F)
```

Elimination continues until we reach an $`I_{r \times r}`$ identity matrix.

---

**Ex 1:**

```math
A = \begin{pmatrix} 1 & 2 & 11 & 17 \\ 3 & 7 & 37 & 57 \end{pmatrix}
```

$`\xrightarrow{R_2 - 3R_1}`$ $`\begin{pmatrix} 1 & 2 & 11 & 17 \\ 0 & 1 & 4 & 6 \end{pmatrix}`$ $`\xrightarrow{R_1 - 2R_2}`$ $`\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = R`$

$`\text{rank}(A) = 2`$, $`n = 4`$, $`n - r = 2`$.

$`A = (W \quad H)`$ where $`W = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}`$ (independent columns) and $`H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix}`$ (dependent columns).

$`R = (I \quad F)`$. The elimination inverted $`W`$. This is the same as multiplying $`W^{-1}`$ to $`A`$:

```math
W^{-1}A = W^{-1}(W \quad H) = (I \quad W^{-1}H) = (I \quad F) = R
```

$`W^{-1}H = F \implies H = WF`$.

Note that $`W`$ consists of independent columns. $`F`$ tells how to combine the independent columns of $`A`$.

```math
H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix} = WF = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 3 & 5 \\ 4 & 6 \end{pmatrix}
```

---

**Ex 2:**

```math
A = \begin{pmatrix} 1 & 7 & 3 & 35 \\ 2 & 14 & 6 & 70 \\ 2 & 14 & 9 & 97 \end{pmatrix}
```

$`\xrightarrow{R_2 - 2R_1, \; R_3 - 2R_1}`$ $`\begin{pmatrix} 1 & 7 & 3 & 35 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 3 & 27 \end{pmatrix}`$ $`\xrightarrow{R_1 - R_3, \; R_3/3}`$ $`\begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 9 \end{pmatrix}`$ $`\xrightarrow{\text{swap } R_2 \text{ and } R_3}`$

```math
R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q
```

where $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}`$, $`F = \begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}`$, and $`Q`$ is a permutation matrix reordering columns.

Remove the row of zeros: $`R = (I \quad F)Q`$.

After column permutation $`PAQ`$:

$`PAQ = \begin{pmatrix} 1 & 3 & 7 & 35 \\ 2 & 9 & 14 & 97 \\ 2 & 6 & 14 & 70 \end{pmatrix}`$ with $`W = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}`$, $`H = \begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix}`$

$`W^{-1}H = F \implies H = WF`$:

```math
\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}
```

**Recall Ex 2:** The identity matrix appears in the 1st and 3rd columns of $`R_0`$, meaning that the 1st and 3rd columns of $`A`$ are independent. The 2nd and 4th columns of $`A`$ are linear combinations of the independent columns:

```math
\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}
```

**The Matrix Factorization $`A = CR`$ and $`N(A)`$:**

Apply elimination to reduce $`A`$ to $`R_0`$. Then "$`I`$" in $`R_0`$ locates the matrix $`C`$ of independent columns in $`A`$. Removing any zero rows from $`R_0`$ produces the row matrix $`R`$ in $`A = CR`$.

Let's find two special solutions $`\mathbf{s}_1, \mathbf{s}_2`$ for Ex 2:

```math
R = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \end{pmatrix}
```

```math
R\mathbf{s} = \begin{pmatrix} s_1 + 7s_2 + 8s_4 \\ s_3 + 9s_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

Take $`s_2 = 1, s_4 = 0`$: $`\Rightarrow s_3 = 0, s_1 = -7`$. $`\therefore \mathbf{s}_1 = \begin{pmatrix} -7 \\ 1 \\ 0 \\ 0 \end{pmatrix}`$

Take $`s_2 = 0, s_4 = 1`$: $`\Rightarrow s_1 = -8, s_3 = -9`$. $`\therefore \mathbf{s}_2 = \begin{pmatrix} -8 \\ 0 \\ -9 \\ 1 \end{pmatrix}`$

---

<br>

## 4. The Complete Solution to Ax = b (3.3)

### 4.1 Structure of the Complete Solution

1. **Complete Solution to $`A\mathbf{x} = \mathbf{b}`$:**

```math
\mathbf{x} = \text{one particular solution } \mathbf{x}_p + \text{any } \mathbf{x}_n \text{ in the null space}
```

2. **Elimination** on $`A\mathbf{x} = \mathbf{b}`$ leads to $`R_0 \mathbf{x} = \mathbf{d}`$. Solvable when zero rows of $`R_0`$ have zeros in $`\mathbf{d}`$.

```math
\begin{pmatrix} R \\ 0 \; 0 \; \cdots \; 0 \\ 0 \; 0 \; \cdots \; 0 \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ 0 \\ 0 \end{pmatrix}
```

3. When $`R_0 \mathbf{x} = \mathbf{d}`$ is solvable, one particular solution $`\mathbf{x}_p`$ has **all free variables equal to zero**.

4. $`A`$ has **full column rank** $`r = n`$ when its nullspace $`N(A) = \lbrace\text{zero vector}\rbrace`$. No free variables.

5. $`A`$ has **full row rank** $`r = m`$ when its column space $`C(A)`$ is $`\mathbb{R}^m`$: $`A\mathbf{x} = \mathbf{b}`$ is always solvable.

### 4.2 Worked Example: Finding Particular Solutions

In the previous section, we found solutions to $`A\mathbf{x} = \mathbf{0}`$. In this section, we find solutions to $`A\mathbf{x} = \mathbf{b}`$.

Row operations on the left side must act on the right side. Use the **augmented matrix** $`(A | \mathbf{b})`$.

```math
\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 1 & 3 & 1 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 7 \end{pmatrix} \quad \Rightarrow \quad (A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & 1 \\ 0 & 0 & 1 & 4 & | & 6 \\ 1 & 3 & 1 & 6 & | & 7 \end{pmatrix}
```

After elimination ($`R_3 - R_1 - R_2`$):

```math
\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} \quad \Rightarrow \quad (R_0 | \mathbf{d})
```

The last equation is $`0 = 0`$. Let's consider a general $`\mathbf{b}`$:

```math
(A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 1 & 3 & 1 & 6 & | & b_3 \end{pmatrix} \xrightarrow{R_3 - R_1 - R_2} \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 0 & 0 & 0 & 0 & | & b_3 - b_1 - b_2 \end{pmatrix} = (R_0 | \mathbf{d})
```

We can get $`0 = 0`$ in the third equation if $`b_3 - b_1 - b_2 = 0`$.

**One Particular Solution $`A\mathbf{x}_p = \mathbf{b}`$:**

```math
R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} = \mathbf{d}
```

$`\text{rank}(R_0) = 2`$, $`n = 4`$, $`n - r = 2`$ free variables.

**Take $`x_2 = 1, x_4 = 0`$:**

$`x_1 + 3 = 1 \Rightarrow x_1 = -2`$, $`x_3 = 6`$

$`\therefore \mathbf{x} = \begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix}`$

**Take $`x_2 = 0, x_4 = 1`$:**

$`x_1 + 2 = 1 \Rightarrow x_1 = -1`$, $`x_3 + 4 = 6 \Rightarrow x_3 = 2`$

$`\therefore \mathbf{x} = \begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix}`$

**Take $`x_2 = 0, x_4 = 2`$:**

$`x_1 + 2 \cdot 2 = 1 \Rightarrow x_1 = -3`$, $`x_3 + 4 \cdot 2 = 6 \Rightarrow x_3 = -2`$

$`\therefore \mathbf{x} = \begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix}`$

We can find infinite number of solutions due to the free variables.

### 4.3 The Complete Solution Decomposition

The solutions can be decomposed into the **particular solution** and the **special solutions** in the nullspace.

```math
R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`R_0 = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q`$ where $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}`$, $`F = \begin{pmatrix} 3 & 2 \\ 0 & 4 \end{pmatrix}`$

The nullspace vectors are the columns of $`Q^T \begin{pmatrix} -F \\ I \end{pmatrix}`$:

```math
\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} -3 & -2 \\ 0 & -4 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} -3 & -2 \\ 1 & 0 \\ 0 & -4 \\ 0 & 1 \end{pmatrix}
```

**Revisit the solutions to $`R_0 \mathbf{x} = \mathbf{d}`$:**

```math
\begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix}
```

```math
\begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}
```

```math
\begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -4 \\ 0 \\ -8 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + 2\begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}
```

**Q. How can we find the particular solution?**

Take $`x_2 = 0, x_4 = 0`$ (set all free variables to zero):

```math
R_0 \mathbf{x} = \begin{pmatrix} x_1 + 3x_2 + 2x_4 \\ x_3 + 4x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix}
```

$`\Rightarrow x_1 = 1, \; x_3 = 6`$

```math
\therefore \mathbf{x}_p = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix}
```

**The complete solution** $`\mathbf{x}_p + \mathbf{x}_n`$ to $`A\mathbf{x} = \mathbf{b}`$ becomes:

```math
\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + x_2 \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_4 \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}
```

where $`\mathbf{x}_p`$ is the particular solution and the other terms are nullspace vectors.

### 4.4 Full Column Rank: r = n

**Q. What happens if $`m = n = r`$ for $`\mathbf{x}_p, \mathbf{x}_n`$?**

**A.** $`\mathbf{x}_n = \mathbf{0}`$.

$`A\mathbf{x} = \mathbf{b} \iff A(\mathbf{x}_p + \mathbf{x}_n) = \mathbf{b} \implies A(\mathbf{x}_p + \mathbf{0}) = \mathbf{b} \implies \mathbf{x}_p = A^{-1}\mathbf{b}`$

### 4.5 Solvability Conditions

**Ex 1.** Find the condition on $`(b_1, b_2, b_3)`$ for $`A\mathbf{x} = \mathbf{b}`$ to be solvable if:

```math
A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ -2 & -3 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \\ b_3 \end{pmatrix}
```

```math
(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & | & b_1 \\ 1 & 2 & | & b_2 \\ -2 & -3 & | & b_3 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 + 2R_1} \begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & -1 & | & b_3 + 2b_1 \end{pmatrix} \xrightarrow{R_3 + R_2}
```

```math
\begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & | & 2b_1 - b_2 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} = (R_0, \mathbf{d})
```

If $`b_3 + b_2 + b_1 = 0`$, then $`A\mathbf{x} = \mathbf{b}`$ is solvable. i.e., $`b_3 + b_2 + b_1 = 0`$ is the condition to put $`\mathbf{b}`$ in the column space of $`A`$.

$`\text{rank}(A) = 2`$, $`n = 2`$, $`n - r = 0`$. No free variables. $`N(A) = \lbrace\mathbf{0}\rbrace`$, $`\mathbf{x}_n = \mathbf{0}`$.

```math
\mathbf{x}_p = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix}
```

```math
\therefore \mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

If $`b_3 + b_2 + b_1 \neq 0`$, then there is no solution to $`A\mathbf{x} = \mathbf{b}`$.

**For full rank $`r = n`$:**

```math
R_0 = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix}
```

1. The matrix $`A`$ has $`n`$ independent columns.
2. The null space of $`A`$ is $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$.
3. If $`A\mathbf{x} = \mathbf{b}`$ has a solution, it has only **one** solution.

```math
R_0 \mathbf{x} = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ O_{(m-n) \times 1} \end{pmatrix}
```

The bottom $`m - n`$ rows give $`m - n`$ conditions for $`\mathbf{b}`$ to be in the column space of $`A`$.

### 4.6 Full Row Rank and the Complete Solution

$`A \in \mathbb{R}^{m \times n}`$, $`m \leq n`$ (short and wide matrix).

A matrix has **full row rank** if $`r = m`$.

**Ex 2.** $`A\mathbf{x} = \mathbf{b}`$ has $`n = 3`$ unknowns but only $`m = 2`$ equations:

```math
x + y + z = 3, \quad x + 2y - z = 4
```

```math
(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 1 & 2 & -1 & | & 4 \end{pmatrix} \xrightarrow{R_2 - R_1} \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & 3 & | & 2 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} = (R | \mathbf{d})
```

$`\text{rank}(A) = 2`$, $`n = 3`$, $`n - r = 1`$. 1 free variable, 1 special solution.

**i) $`\mathbf{x}_p`$:**

```math
\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}
```

```math
x + 3z = 2, \quad y - 2z = 1
```

Take $`z = 0`$: $`\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{x}_p = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}`$ (particular solution).

**ii) $`\mathbf{x}_n`$:**

```math
\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

Take $`z = 1`$: $`x + 3 = 0 \Rightarrow x = -3`$, $`y - 2 = 0 \Rightarrow y = 2`$.

$`\therefore \mathbf{x}_n = \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}`$

**iii) Complete solution:**

```math
\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix} + \alpha \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}
```

This is a **line** through $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$ direction (for $`A\mathbf{x}_n = \mathbf{0}`$) shifted by $`\mathbf{x}_p`$ (for $`A\mathbf{x} = \mathbf{b}`$).

**Every matrix $`A`$ with full row rank ($`r = m`$) has the following properties:**
1. All rows have pivots, and $`R_0`$ has no zero rows. $`R_0 = R`$.
2. $`A\mathbf{x} = \mathbf{b}`$ has a solution for $`\forall \mathbf{b}`$.
3. The column space of $`A`$ is the whole space $`\mathbb{R}^m`$.

$`A\mathbf{x} = \mathbf{b} \to R\mathbf{x} = \mathbf{d}`$ where $`R = (I_{m \times m} \quad F_{m \times (n-m)})`$.

4. If $`m < n`$, $`A\mathbf{x} = \mathbf{b}`$ is **underdetermined** (many solutions).

With full row rank ($`r = m`$), $`m`$ rows are linearly independent. The columns of $`A^T`$ are LI. The nullspace of $`A^T`$ is $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$.

### 4.7 Four Possibilities for Linear Equations

Four possibilities for linear equations depend on the rank $`r`$:

| Case | Shape | $`R_0`$ form | Solutions to $`A\mathbf{x} = \mathbf{b}`$ |
|:-----|:------|:-----------|:----------------------------------------|
| $`r = m`$ and $`r = n`$ | Square, invertible | $`(I)`$ | 1 solution |
| $`r = m`$ and $`r < n`$ | Short, wide | $`(I \quad F)`$ | $`\infty`$ solutions |
| $`r < m`$ and $`r = n`$ | Tall, thin | $`\begin{pmatrix} I \\ 0 \end{pmatrix}`$ | 0 or 1 solution |
| $`r < m`$ and $`r < n`$ | Not full rank | $`\begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix}`$ | 0 or $`\infty`$ solutions |

---

<br>
