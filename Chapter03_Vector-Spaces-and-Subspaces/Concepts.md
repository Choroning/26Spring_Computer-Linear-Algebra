# Chapter 3 Lecture — Vector Spaces and Subspaces

> **Last Updated:** 2026-03-31

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

## 1. Overview: The Four Fundamental Subspaces

Chapter 3 covers five major topics:

- **3.1** Vector Spaces and Subspaces — How do we define a vector space? Key operations are $\mathbf{u} + \mathbf{v}$ and $c\mathbf{u}$. There are 8 rules that vectors $\mathbf{u}$ and scalar $c$ must satisfy.
- **3.2** The Nullspace of $A$: Solving $A\mathbf{x} = \mathbf{0}$
- **3.3** The Complete Solution to $A\mathbf{x} = \mathbf{b}$
- **3.4** Independence, **Basis**, and **Dimension** — A set of vectors that describes the space. Let $A \in \mathbb{R}^{n \times n}$. $A$ has $r$ independent columns $\Rightarrow$ $C(A)$ has dimension $r$. $(n - r)$ special solutions to $A\mathbf{x} = \mathbf{0}$ are a basis for $N(A)$, the null space of $A$.
- **3.5** Dimensions of the Four Subspaces:

| Subspace | Dimension |
|:---------|:----------|
| Column space of $A$ | $r$ |
| Row space of $A$ | $r$ |
| Null space of $A$ | $n - r$ |
| Null space of $A^T$ (left nullspace) | $m - r$ |

---

<br>

## 2. Vector Spaces and Subspaces (3.1)

### 2.1 Requirements for a Vector Space

1. **Requirement:** All the linear combinations $c\mathbf{u} + d\mathbf{w}$ must stay in the vector space.
2. The **row space** of $A$ is "spanned" by the rows of $A$. The columns of $A$ span the column space $C(A)$.
3. **Matrices** $M_1$ to $M_n$ and **functions** $f_1$ to $f_n$ span **matrix spaces** and **function spaces**.

### 2.2 The Space R^n

The space $\mathbb{R}^n$ contains all column vectors $\mathbf{v}$ of length $n$.

$$
\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} \in \mathbb{R}^n \quad \text{where } x_1, x_2, \ldots, x_n \in \mathbb{R} \text{ (real numbers)}
$$

Examples:
- $x \in \mathbb{R}^1$
- $\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} \in \mathbb{R}^2$
- $\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} \in \mathbb{R}^3$

**Note:** When $x_1, x_2, \ldots, x_n \in \mathbb{C}$ (complex numbers), the space becomes $\mathbb{C}^n$.

**Note:** If $c, d \in \mathbb{R}$ and $\mathbf{u}, \mathbf{w} \in \mathbb{R}^n$, then $c\mathbf{u} + d\mathbf{w} \in \mathbb{R}^n$. The linear combinations stay in the vector space $\mathbb{R}^n$.

### 2.3 Definition of a Vector Space

**Vector space** (= linear space) $V$ is a set whose elements (vectors) can be added together and multiplied by numbers.

Let $\mathbf{u}, \mathbf{w} \in V$ (vector space):
1. $\mathbf{u} + \mathbf{w} \in V$ — vector addition
2. $\alpha \mathbf{u} \in V \quad \forall \alpha \in \mathbb{F}$ — scalar multiplication

**A field** is a set on which addition, subtraction, multiplication, and division are defined.

Examples of fields: $\mathbb{R}$ (real numbers), $\mathbb{C}$ (complex numbers).

### 2.4 The Eight Axioms

Let $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ and $c, d \in \mathbb{F}$. A vector space must satisfy the following eight axioms:

**(1) Associativity of vector addition:**

$$\mathbf{u} + (\mathbf{v} + \mathbf{w}) = (\mathbf{u} + \mathbf{v}) + \mathbf{w}$$

**(2) Commutativity of vector addition:**

$$\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$$

**(3) Identity element of vector addition (zero vector):**

$$\exists!\; \mathbf{0} \in V \text{ s.t. } \mathbf{v} + \mathbf{0} = \mathbf{v} \quad \forall \mathbf{v} \in V$$

**(4) Inverse element of vector addition:**

$$\forall \mathbf{v} \in V, \; \exists!\; -\mathbf{v} \in V \text{ s.t. } \mathbf{v} + (-\mathbf{v}) = \mathbf{0}$$

($-\mathbf{v}$ is the additive inverse of $\mathbf{v}$)

**(5) Compatibility of scalar multiplication with field multiplication:**

$$c(d\mathbf{v}) = (cd)\mathbf{v}$$

**(6) Identity element of scalar multiplication:**

$$1 \cdot \mathbf{v} = \mathbf{v}$$

**(7) Distributivity of scalar multiplication w.r.t. vector addition:**

$$c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$$

**(8) Distributivity of scalar multiplication w.r.t. field addition:**

$$(c + d)\mathbf{u} = c\mathbf{u} + d\mathbf{u}$$

### 2.5 Consequences of the Axioms

From the eight axioms, the following properties hold:

$$0\mathbf{u} = \mathbf{0}$$

$$c\mathbf{0} = \mathbf{0}$$

$$(-1)\mathbf{u} = -\mathbf{u}$$

$$c\mathbf{v} = \mathbf{0} \implies c = 0 \text{ or } \mathbf{v} = \mathbf{0}$$

### 2.6 Examples: Is It a Vector Space?

**Q.** Is the set of all positive vectors $\mathbb{X}$, where $\mathbb{X} \ni \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}$ with every $v_i > 0$, a vector space?

**A.** No. $-\mathbf{v} \notin \mathbb{X}$.

---

**Q.** Let $\mathbb{X}$ be the set of solutions to $A\mathbf{x} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$. Is $\mathbb{X}$ a vector space?

**A.** No. If $\mathbf{u}, \mathbf{w} \in \mathbb{X}$, then $A\mathbf{u} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$ and $A\mathbf{w} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$. But $\mathbf{u} + \mathbf{w} \notin \mathbb{X}$ because $A(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w} = \begin{pmatrix} 2 \\ \vdots \\ 2 \end{pmatrix} \neq \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$.

---

**Q.** Is a line in $\mathbb{R}^n$ a vector space?

A line is a collection of points. $\mathbf{q}$ represents a position on the line:

$$\mathbf{q} = \mathbf{p} + t\mathbf{d}, \quad t \in \mathbb{R}$$

$\mathbf{p}$ and $\mathbf{d}$ are fixed. $t$ varies from $-\infty$ to $\infty$.

Let $\mathbb{X}$ be the set of all vectors $\mathbf{q}$: $\mathbb{X} = \{\mathbf{q} : \mathbf{p} + t\mathbf{d} \;\forall\; t \in \mathbb{R}\}$.

Take $\mathbf{a}, \mathbf{b} \in \mathbb{X}$. Does $\mathbf{a} + \mathbf{b}$ belong to $\mathbb{X}$?

$$\mathbf{a} = \mathbf{p} + t_1 \mathbf{d}$$
$$\mathbf{b} = \mathbf{p} + t_2 \mathbf{d}$$
$$\mathbf{a} + \mathbf{b} = 2\mathbf{p} + (t_1 + t_2)\mathbf{d} \notin \mathbb{X} \quad \text{if } \mathbf{p} \neq \mathbf{0}$$

However:
- $\mathbf{a} + \mathbf{b} = (t_1 + t_2)\mathbf{d} \in \mathbb{X}$ if $\mathbf{p} = \mathbf{0}$
- $c\mathbf{a} = ct_1\mathbf{d} \in \mathbb{X}$ if $\mathbf{p} = \mathbf{0}$

**A line passing through the origin is a vector space.** The line through $\mathbf{0}$ in $\mathbb{R}^n$ is a **subspace** of $\mathbb{R}^n$ — a vector space inside another vector space.

### 2.7 Examples of Vector Spaces

- $\mathbb{R}^n$ is a vector space.
- $\mathbb{Z} = \{\mathbf{0}\}$ is a vector space (the smallest vector space):
  - i) $\mathbf{0} + \mathbf{0} = \mathbf{0} \in \mathbb{Z}$
  - ii) $c\mathbf{0} = \mathbf{0} \in \mathbb{Z}$

$\mathbb{Z}$ is the **null space of an invertible matrix**. If the only solution to $A\mathbf{x} = \mathbf{0}$ is the zero vector $\mathbf{x} = \mathbf{0}$, then the columns of $A$ are linearly independent (LI) and the nullspace of $A$ is $\mathbb{Z}$.

### 2.8 Generalized Vector Spaces

**Note:** A vector space (= a linear space) is a set whose elements can be added together and multiplied by numbers.

i.e., $\mathbf{u}, \mathbf{w} \in V \implies c\mathbf{u} + d\mathbf{w} \in V$.

**Matrices and functions can be considered as vectors.**

**Vector space of matrices:**

$A, B \in \mathbb{R}^{m \times n} \implies cA + dB \in \mathbb{R}^{m \times n}$

The set of all matrices of a fixed size forms a vector space. (Check if they satisfy the eight rules.)

**Vector space of functions:**

Let $\mathbb{F}$ be a set of functions that take elements from $\mathbb{R}$ and map them to a real number:

$$\mathbb{F} = \{f \mid f: \mathbb{R} \to \mathbb{R}\}$$

Let $f, g \in \mathbb{F}$, $c, d \in \mathbb{R}$. Then $cf + dg \in \mathbb{F}$.

Define:
- $(f + g)(x) = f(x) + g(x)$
- $(cf)(x) = c(f(x))$

**Example:** $\mathbb{F}$ = the line of functions $y = ce^x$.

$\mathbb{F} = \{f \mid f: \mathbb{R} \ni x \to ce^x \in \mathbb{R}, \;\forall c \in \mathbb{R}\}$

$f(x) = e^x$, $g(x) = 2e^x$.

$(f + g)(x) = f(x) + g(x) = e^x + 2e^x = 3e^x \in \mathbb{F}$

$(cf)(x) = c(f(x)) = ce^x \in \mathbb{F}$

**Remark:** A **set** is just a collection of elements without any additional structure. A **space** is a set along with additional structures defined on it. Example: a vector space is a set of vectors where vector addition and scalar multiplication is defined. A vector space should satisfy the 8 rules.

### 2.9 Subspaces of Vector Spaces

Consider a vector space $\mathbb{R}^n$, where $\mathbf{v} \in \mathbb{R}^n$ is a column vector with $n$ components.

There are important vector spaces inside $\mathbb{R}^n$. Those are **subspaces** of $\mathbb{R}^n$.

A plane is a vector space. If $\mathbf{v}, \mathbf{w} \in \mathbb{R}^2$, then $c\mathbf{v} + d\mathbf{w} \in \mathbb{R}^2$.

A plane through the origin $\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$ in $\mathbb{R}^3$ is a vector space. This plane is not $\mathbb{R}^2$ because $\mathbf{u}, \mathbf{v} \in \mathbb{R}^3$. The plane is a vector space inside $\mathbb{R}^3$. This is a **subspace** of the full vector space $\mathbb{R}^3$.

### 2.10 Definition of a Subspace

**Def.** A **subspace** of a vector space is a set of vectors (including $\mathbf{0}$) that satisfies:

> i) $\mathbf{u} + \mathbf{w}$ is in the subspace
>
> ii) $c\mathbf{u}$ is in the subspace

if $\mathbf{u}, \mathbf{w} \in$ subspace and $c \in \mathbb{R}$ (or $\mathbb{F}$).

The set of vectors is "**closed**" under addition $\mathbf{u} + \mathbf{w}$ and multiplication $c\mathbf{u}$.

Conditions i) and ii) mean: **all linear combinations stay in the subspace.**

**Note:** Every subspace contains the zero vector. From ii), $c = 0$ implies that $0\mathbf{u} = \mathbf{0}$ is in the subspace.

### 2.11 Examples of Subspaces and Non-Subspaces

**Q:** Is the plane $z = 5$ in $\mathbb{R}^3$ a subspace? **No.** (It does not contain the origin.)

Lines through the origin are subspaces.

$\mathbb{R}^3$ is a subspace of itself.

The single vector $\mathbb{Z} = \{\mathbf{0}\}$ is a subspace of $\mathbb{R}^3$.

---

**Ex 1.** $\mathbb{R}^2$ is a vector space. Is the 1st quadrant a subspace?

$$\mathcal{U} = \left\{\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0, y \geq 0 \right\}$$

Take $\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}$, $c = -1$.

Then $c\mathbf{u} = \begin{pmatrix} -2 \\ -3 \end{pmatrix} \notin \mathcal{U}$.

This violates rule ii). $\mathcal{U}$ is **not** a subspace of $\mathbb{R}^2$.

---

**Ex 2.** $\mathcal{U} = \left\{\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0 \text{ and } y \geq 0, \text{ or } x \leq 0 \text{ and } y \leq 0 \right\}$

Is $\mathcal{U}$ a subspace?

Take $\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}$, $\mathbf{w} = \begin{pmatrix} -3 \\ -2 \end{pmatrix} \in \mathcal{U}$.

But the sum $\mathbf{u} + \mathbf{w} = \begin{pmatrix} -1 \\ 1 \end{pmatrix} \notin \mathcal{U}$.

$\mathcal{U}$ is **not** a subspace.

---

**Ex 3.** $\mathbb{M}$ is a vector space of $2 \times 2$ matrices $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$.

$\mathcal{U}$ is a set of all upper triangular matrices $\begin{pmatrix} a & b \\ 0 & d \end{pmatrix}$.

$\mathbb{D}$ is a set of all diagonal matrices $\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix}$.

**Both $\mathcal{U}$ and $\mathbb{D}$ are subspaces of $\mathbb{M}$:**

- $A, B \in \mathcal{U} \implies cA + dB \in \mathcal{U}$
- $A, B \in \mathbb{D} \implies cA + dB \in \mathbb{D}$

Note that:

$$\begin{pmatrix} a & b \\ c & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + c\begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$$

### 2.12 Column Space and Row Space

**Column Space of $A$:**

$$A\mathbf{x} = \mathbf{b}$$

$$\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}$$

$$= x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n$$

The **column space of $A$** is a vector space made up of column vectors.

To solve $A\mathbf{x} = \mathbf{b}$ is to express $\mathbf{b}$ as a linear combination of the columns. The right side vector $\mathbf{b}$ has to be in the column space of $A$.

> **The equations $A\mathbf{x} = \mathbf{b}$ are solvable iff $\mathbf{b}$ is in the column space of $A$.**

**The Row Space of $A$:**

The rows of $A$ are the columns of $A^T$.

$$A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} \implies A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}$$

$A^T$ has $m$ column vectors.

**The row space of $A$ is the column space of $A^T$.**

The equations $A^T \mathbf{y} = \mathbf{c}$ are solvable iff $\mathbf{c}$ is in the column space of $A^T$ ($= \mathbf{c}$ is in the row space of $A$).

**Example:** Consider a rank 1 matrix $A = \mathbf{u}\mathbf{v}^T$:

$$A = \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix} \begin{pmatrix} v_1 & v_2 & \cdots & v_n \end{pmatrix} = \begin{pmatrix} v_1\mathbf{u} & v_2\mathbf{u} & \cdots & v_n\mathbf{u} \end{pmatrix}$$

$C(A)$ is the line of all column vectors of $A$: $c\mathbf{u}$.

$A^T = \mathbf{v}\mathbf{u}^T = \begin{pmatrix} u_1\mathbf{v} & u_2\mathbf{v} & \cdots & u_m\mathbf{v} \end{pmatrix}$

$C(A^T)$ is the line of all column vectors of $A^T$: $c\mathbf{v}$.

### 2.13 Spanning

**The Columns of $A$ Span the Vector Space $C(A)$.**

Let $\mathbb{S}$ be a set of vectors in $\mathbb{R}^m$. If $\mathbb{S} = \left\{ \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix}, \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_m \end{pmatrix} \right\}$, then $\mathbb{S}$ is **not** a subspace of $\mathbb{R}^m$ because $\mathbf{u} + \mathbf{v} \notin \mathbb{S}$.

If we include all combinations of the vectors in $\mathbb{S}$, then we have a vector space $V$.

**The set $\mathbb{S}$ spans $V$.** $V$ is the smallest vector space containing $\mathbb{S}$.

Consider a matrix $A \in \mathbb{R}^{m \times n}$:
- $n$ columns span the column space $C(A)$
- $m$ columns of $A^T$ span the row space $C(A^T)$

---

<br>

## 3. Computing the Nullspace by Elimination (3.2)

### 3.1 Key Facts: A = CR

$$A = CR$$

1. The nullspace $N(A)$ in $\mathbb{R}^n$ contains all solutions $\mathbf{x}$ to $A\mathbf{x} = \mathbf{0}$. This includes $\mathbf{x} = \mathbf{0}$.
2. Elimination from $A$ to $R_0$ to $R$ does **not** change the nullspace.
3. The reduced row echelon form $R_0 = \text{rref}(A)$ has $I$ in $r$ columns and $F$ in $n - r$ columns.
4. If column $j$ is dependent on previous columns, $A\mathbf{x} = \mathbf{0}$ has a "special solution" with $x_j = 1$.
5. The $n - r$ special solutions to $A\mathbf{x} = \mathbf{0}$ contain $-F$ and $I$.
6. Every short wide matrix with $m < n$ has nonzero solutions to $A\mathbf{x} = \mathbf{0}$ in its nullspace.

### 3.2 Finding All Solutions to Ax = 0

We would like to find all solutions to $A\mathbf{x} = \mathbf{0}$.

If $A \in \mathbb{R}^{n \times n}$ is invertible ($\text{rank}(A) = n$), then the only solution is $\mathbf{x} = \mathbf{0}$. The nullspace of $A$ only contains the zero vector: $N(A) = \{\mathbf{0}\}$.

In general, however, $\text{rank}(A) = r$. That is, $A$ has $r$ independent columns. The other $n - r$ dependent columns of $A$ are combinations of those independent columns. We will find $n - r$ vectors in $N(A)$ which are special solutions to $A\mathbf{x} = \mathbf{0}$.

In Chapter 2, an invertible matrix $A$ is reduced to an upper triangular matrix $U$. For $A \in \mathbb{R}^{m \times n}$, $A\mathbf{x} = \mathbf{0}$ is simplified to $R\mathbf{x} = \mathbf{0}$ (echelon form).

### 3.3 Reduced Row Echelon Form

In this section, we consider $A \in \mathbb{R}^{m \times n}$, and eliminate $A$ into reduced row echelon form: $R_0 = \text{rref}(A)$.

$R_0$ may have zero rows. Removing all zero rows of $R_0$ gives $R$.

**Ex 1:** $R = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = (I \quad F)$

$\text{rank}(R) = 2$, $n = 4$, $n - r = 2$ dependent columns.

$$R\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_1 + 3x_3 + 5x_4 \\ x_2 + 4x_3 + 6x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}$$

How do we find the special solutions? Fix two variables and use $R\mathbf{x} = \mathbf{0}$:

Take $x_3 = 1, x_4 = 0$: $\Rightarrow x_1 + 3 = 0, \; x_2 + 4 = 0$. $\therefore \mathbf{s}_1 = \begin{pmatrix} -3 \\ -4 \\ 1 \\ 0 \end{pmatrix}$

Take $x_3 = 0, x_4 = 1$: $\Rightarrow x_1 + 5 = 0, \; x_2 + 6 = 0$. $\therefore \mathbf{s}_2 = \begin{pmatrix} -5 \\ -6 \\ 0 \\ 1 \end{pmatrix}$

### 3.4 Special Solutions and the Nullspace Basis

The two special solutions $\mathbf{s}_1, \mathbf{s}_2$ are in the nullspace of $R$:

$$\mathbf{s}_1, \mathbf{s}_2 \in N(R)$$

$$R\mathbf{s}_1 = \mathbf{0}, \quad R\mathbf{s}_2 = \mathbf{0}$$

$$\Rightarrow R(c_1 \mathbf{s}_1 + c_2 \mathbf{s}_2) = \mathbf{0}$$

Linear combinations of $\mathbf{s}_1$ and $\mathbf{s}_2$ are in the nullspace of $R$. **The special solutions $\mathbf{s}_1$ and $\mathbf{s}_2$ are a basis for the nullspace.**

---

**Ex 2:** $R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix}$ (reduced row echelon form with a row of zeros)

$$R_0 \mathbf{x} = \begin{pmatrix} x_1 + 7x_2 + 8x_4 \\ x_3 + 9x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$\text{rank}(R_0) = 2$, $n = 4$, $n - r = 2$. Find two special solutions.

Take $x_1 = 0, x_2 = 1$:

$7 + 8x_4 = 0 \Rightarrow x_4 = -7/8$, $x_3 = -9(-7/8) = 63/8$

$\Rightarrow \mathbf{x} = \begin{pmatrix} 0 \\ 1 \\ 63/8 \\ -7/8 \end{pmatrix} \Rightarrow \mathbf{s}_1 = \begin{pmatrix} 0 \\ 8 \\ 63 \\ -7 \end{pmatrix}$

Take $x_1 = 1, x_2 = 0$:

$1 + 8x_4 = 0 \Rightarrow x_4 = -1/8$, $x_3 = +9/8$

$\Rightarrow \mathbf{x} = \begin{pmatrix} 1 \\ 0 \\ 9/8 \\ -1/8 \end{pmatrix} \Rightarrow \mathbf{s}_2 = \begin{pmatrix} 8 \\ 0 \\ 9 \\ -1 \end{pmatrix}$

We can also take $x_2 = 1, x_4 = 0$ and $x_2 = 0, x_4 = 1$ to find $x_1$ and $x_3$.

### 3.5 The Nullspace Matrix: Columns of (-F; I)

In this section, we consider $A \in \mathbb{R}^{m \times n}$, eliminate $A$ into reduced row echelon form $R_0 = \text{rref}(A)$.

$R_0$ may have zero rows. Removing all zero rows of $R_0$ gives $R$.

$$R = (I \quad F)Q$$

where $Q$ is a permutation matrix (reordering columns so that pivot columns come first).

**Recall** $R\mathbf{x} = \mathbf{0}$, i.e., $(I \quad F)\mathbf{x} = \mathbf{0}$ (for Ex 1):

$$\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

Setting free variables to $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$ and $\begin{pmatrix} 0 \\ 1 \end{pmatrix}$:

$$I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 3 \\ 4 \end{pmatrix}$$

$$I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 5 \\ 6 \end{pmatrix}$$

Equivalently:

$$\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} -3 & -5 \\ -4 & -6 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$$

That is:

$$\boxed{(I \quad F) \begin{pmatrix} -F \\ I \end{pmatrix} = O}$$

**The two special solutions to $(I \quad F)\mathbf{x} = \mathbf{0}$ are the columns of $\begin{pmatrix} -F \\ I \end{pmatrix}$.**

### 3.6 The Matrix Factorization A = CR and N(A)

**Elimination from $A$ to $\text{rref}(A)$: Reduced Row Echelon Form**

Apply elimination to reduce $A$ to $R_0$. Then "$I$" in $R_0$ locates the matrix $C$ of independent columns in $A$. Removing any zero rows from $R_0$ produces the row matrix $R$ in $A = CR$.

In Chapter 2, $A$ was square and invertible:
- $A\mathbf{x} = \mathbf{b}$ $\xrightarrow{a) \text{ elimination}}$ $U\mathbf{x} = \mathbf{c}$ $\xrightarrow{b) \text{ back substitution}}$ $\mathbf{x} = U^{-1}\mathbf{c}$

For any matrix $A$ of rank $r$:

$$A = CR = C(I \quad F)$$

Elimination continues until we reach an $I_{r \times r}$ identity matrix.

---

**Ex 1:**

$$A = \begin{pmatrix} 1 & 2 & 11 & 17 \\ 3 & 7 & 37 & 57 \end{pmatrix}$$

$\xrightarrow{R_2 - 3R_1}$ $\begin{pmatrix} 1 & 2 & 11 & 17 \\ 0 & 1 & 4 & 6 \end{pmatrix}$ $\xrightarrow{R_1 - 2R_2}$ $\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = R$

$\text{rank}(A) = 2$, $n = 4$, $n - r = 2$.

$A = (W \quad H)$ where $W = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}$ (independent columns) and $H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix}$ (dependent columns).

$R = (I \quad F)$. The elimination inverted $W$. This is the same as multiplying $W^{-1}$ to $A$:

$$W^{-1}A = W^{-1}(W \quad H) = (I \quad W^{-1}H) = (I \quad F) = R$$

$W^{-1}H = F \implies H = WF$.

Note that $W$ consists of independent columns. $F$ tells how to combine the independent columns of $A$.

$$H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix} = WF = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 3 & 5 \\ 4 & 6 \end{pmatrix}$$

---

**Ex 2:**

$$A = \begin{pmatrix} 1 & 7 & 3 & 35 \\ 2 & 14 & 6 & 70 \\ 2 & 14 & 9 & 97 \end{pmatrix}$$

$\xrightarrow{R_2 - 2R_1, \; R_3 - 2R_1}$ $\begin{pmatrix} 1 & 7 & 3 & 35 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 3 & 27 \end{pmatrix}$ $\xrightarrow{R_1 - R_3, \; R_3/3}$ $\begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 9 \end{pmatrix}$ $\xrightarrow{\text{swap } R_2 \text{ and } R_3}$

$$R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q$$

where $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$, $F = \begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}$, and $Q$ is a permutation matrix reordering columns.

Remove the row of zeros: $R = (I \quad F)Q$.

After column permutation $PAQ$:

$PAQ = \begin{pmatrix} 1 & 3 & 7 & 35 \\ 2 & 9 & 14 & 97 \\ 2 & 6 & 14 & 70 \end{pmatrix}$ with $W = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}$, $H = \begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix}$

$W^{-1}H = F \implies H = WF$:

$$\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}$$

**Recall Ex 2:** The identity matrix appears in the 1st and 3rd columns of $R_0$, meaning that the 1st and 3rd columns of $A$ are independent. The 2nd and 4th columns of $A$ are linear combinations of the independent columns:

$$\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}$$

**The Matrix Factorization $A = CR$ and $N(A)$:**

Apply elimination to reduce $A$ to $R_0$. Then "$I$" in $R_0$ locates the matrix $C$ of independent columns in $A$. Removing any zero rows from $R_0$ produces the row matrix $R$ in $A = CR$.

Let's find two special solutions $\mathbf{s}_1, \mathbf{s}_2$ for Ex 2:

$$R = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \end{pmatrix}$$

$$R\mathbf{s} = \begin{pmatrix} s_1 + 7s_2 + 8s_4 \\ s_3 + 9s_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

Take $s_2 = 1, s_4 = 0$: $\Rightarrow s_3 = 0, s_1 = -7$. $\therefore \mathbf{s}_1 = \begin{pmatrix} -7 \\ 1 \\ 0 \\ 0 \end{pmatrix}$

Take $s_2 = 0, s_4 = 1$: $\Rightarrow s_1 = -8, s_3 = -9$. $\therefore \mathbf{s}_2 = \begin{pmatrix} -8 \\ 0 \\ -9 \\ 1 \end{pmatrix}$

---

<br>

## 4. The Complete Solution to Ax = b (3.3)

### 4.1 Structure of the Complete Solution

1. **Complete Solution to $A\mathbf{x} = \mathbf{b}$:**

$$\mathbf{x} = \text{one particular solution } \mathbf{x}_p + \text{any } \mathbf{x}_n \text{ in the null space}$$

2. **Elimination** on $A\mathbf{x} = \mathbf{b}$ leads to $R_0 \mathbf{x} = \mathbf{d}$. Solvable when zero rows of $R_0$ have zeros in $\mathbf{d}$.

$$\begin{pmatrix} R \\ 0 \; 0 \; \cdots \; 0 \\ 0 \; 0 \; \cdots \; 0 \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ 0 \\ 0 \end{pmatrix}$$

3. When $R_0 \mathbf{x} = \mathbf{d}$ is solvable, one particular solution $\mathbf{x}_p$ has **all free variables equal to zero**.

4. $A$ has **full column rank** $r = n$ when its nullspace $N(A) = \{\text{zero vector}\}$. No free variables.

5. $A$ has **full row rank** $r = m$ when its column space $C(A)$ is $\mathbb{R}^m$: $A\mathbf{x} = \mathbf{b}$ is always solvable.

### 4.2 Worked Example: Finding Particular Solutions

In the previous section, we found solutions to $A\mathbf{x} = \mathbf{0}$. In this section, we find solutions to $A\mathbf{x} = \mathbf{b}$.

Row operations on the left side must act on the right side. Use the **augmented matrix** $(A | \mathbf{b})$.

$$\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 1 & 3 & 1 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 7 \end{pmatrix} \quad \Rightarrow \quad (A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & 1 \\ 0 & 0 & 1 & 4 & | & 6 \\ 1 & 3 & 1 & 6 & | & 7 \end{pmatrix}$$

After elimination ($R_3 - R_1 - R_2$):

$$\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} \quad \Rightarrow \quad (R_0 | \mathbf{d})$$

The last equation is $0 = 0$. Let's consider a general $\mathbf{b}$:

$$(A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 1 & 3 & 1 & 6 & | & b_3 \end{pmatrix} \xrightarrow{R_3 - R_1 - R_2} \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 0 & 0 & 0 & 0 & | & b_3 - b_1 - b_2 \end{pmatrix} = (R_0 | \mathbf{d})$$

We can get $0 = 0$ in the third equation if $b_3 - b_1 - b_2 = 0$.

**One Particular Solution $A\mathbf{x}_p = \mathbf{b}$:**

$$R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} = \mathbf{d}$$

$\text{rank}(R_0) = 2$, $n = 4$, $n - r = 2$ free variables.

**Take $x_2 = 1, x_4 = 0$:**

$x_1 + 3 = 1 \Rightarrow x_1 = -2$, $x_3 = 6$

$\therefore \mathbf{x} = \begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix}$

**Take $x_2 = 0, x_4 = 1$:**

$x_1 + 2 = 1 \Rightarrow x_1 = -1$, $x_3 + 4 = 6 \Rightarrow x_3 = 2$

$\therefore \mathbf{x} = \begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix}$

**Take $x_2 = 0, x_4 = 2$:**

$x_1 + 2 \cdot 2 = 1 \Rightarrow x_1 = -3$, $x_3 + 4 \cdot 2 = 6 \Rightarrow x_3 = -2$

$\therefore \mathbf{x} = \begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix}$

We can find infinite number of solutions due to the free variables.

### 4.3 The Complete Solution Decomposition

The solutions can be decomposed into the **particular solution** and the **special solutions** in the nullspace.

$$R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$R_0 = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q$ where $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$, $F = \begin{pmatrix} 3 & 2 \\ 0 & 4 \end{pmatrix}$

The nullspace vectors are the columns of $Q^T \begin{pmatrix} -F \\ I \end{pmatrix}$:

$$\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} -3 & -2 \\ 0 & -4 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} -3 & -2 \\ 1 & 0 \\ 0 & -4 \\ 0 & 1 \end{pmatrix}$$

**Revisit the solutions to $R_0 \mathbf{x} = \mathbf{d}$:**

$$\begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix}$$

$$\begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}$$

$$\begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -4 \\ 0 \\ -8 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + 2\begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}$$

**Q. How can we find the particular solution?**

Take $x_2 = 0, x_4 = 0$ (set all free variables to zero):

$$R_0 \mathbf{x} = \begin{pmatrix} x_1 + 3x_2 + 2x_4 \\ x_3 + 4x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix}$$

$\Rightarrow x_1 = 1, \; x_3 = 6$

$$\therefore \mathbf{x}_p = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix}$$

**The complete solution** $\mathbf{x}_p + \mathbf{x}_n$ to $A\mathbf{x} = \mathbf{b}$ becomes:

$$\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + x_2 \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_4 \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}$$

where $\mathbf{x}_p$ is the particular solution and the other terms are nullspace vectors.

### 4.4 Full Column Rank: r = n

**Q. What happens if $m = n = r$ for $\mathbf{x}_p, \mathbf{x}_n$?**

**A.** $\mathbf{x}_n = \mathbf{0}$.

$A\mathbf{x} = \mathbf{b} \iff A(\mathbf{x}_p + \mathbf{x}_n) = \mathbf{b} \implies A(\mathbf{x}_p + \mathbf{0}) = \mathbf{b} \implies \mathbf{x}_p = A^{-1}\mathbf{b}$

### 4.5 Solvability Conditions

**Ex 1.** Find the condition on $(b_1, b_2, b_3)$ for $A\mathbf{x} = \mathbf{b}$ to be solvable if:

$$A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ -2 & -3 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \\ b_3 \end{pmatrix}$$

$$(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & | & b_1 \\ 1 & 2 & | & b_2 \\ -2 & -3 & | & b_3 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 + 2R_1} \begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & -1 & | & b_3 + 2b_1 \end{pmatrix} \xrightarrow{R_3 + R_2}$$

$$\begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & | & 2b_1 - b_2 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} = (R_0, \mathbf{d})$$

If $b_3 + b_2 + b_1 = 0$, then $A\mathbf{x} = \mathbf{b}$ is solvable. i.e., $b_3 + b_2 + b_1 = 0$ is the condition to put $\mathbf{b}$ in the column space of $A$.

$\text{rank}(A) = 2$, $n = 2$, $n - r = 0$. No free variables. $N(A) = \{\mathbf{0}\}$, $\mathbf{x}_n = \mathbf{0}$.

$$\mathbf{x}_p = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix}$$

$$\therefore \mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

If $b_3 + b_2 + b_1 \neq 0$, then there is no solution to $A\mathbf{x} = \mathbf{b}$.

**For full rank $r = n$:**

$$R_0 = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix}$$

1. The matrix $A$ has $n$ independent columns.
2. The null space of $A$ is $\mathbb{Z} = \{\mathbf{0}\}$.
3. If $A\mathbf{x} = \mathbf{b}$ has a solution, it has only **one** solution.

$$R_0 \mathbf{x} = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ O_{(m-n) \times 1} \end{pmatrix}$$

The bottom $m - n$ rows give $m - n$ conditions for $\mathbf{b}$ to be in the column space of $A$.

### 4.6 Full Row Rank and the Complete Solution

$A \in \mathbb{R}^{m \times n}$, $m \leq n$ (short and wide matrix).

A matrix has **full row rank** if $r = m$.

**Ex 2.** $A\mathbf{x} = \mathbf{b}$ has $n = 3$ unknowns but only $m = 2$ equations:

$$x + y + z = 3, \quad x + 2y - z = 4$$

$$(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 1 & 2 & -1 & | & 4 \end{pmatrix} \xrightarrow{R_2 - R_1} \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & 3 & | & 2 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} = (R | \mathbf{d})$$

$\text{rank}(A) = 2$, $n = 3$, $n - r = 1$. 1 free variable, 1 special solution.

**i) $\mathbf{x}_p$:**

$$\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$$

$$x + 3z = 2, \quad y - 2z = 1$$

Take $z = 0$: $\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$, $\mathbf{x}_p = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}$ (particular solution).

**ii) $\mathbf{x}_n$:**

$$\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

Take $z = 1$: $x + 3 = 0 \Rightarrow x = -3$, $y - 2 = 0 \Rightarrow y = 2$.

$\therefore \mathbf{x}_n = \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}$

**iii) Complete solution:**

$$\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix} + \alpha \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}$$

This is a **line** through $\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$ direction (for $A\mathbf{x}_n = \mathbf{0}$) shifted by $\mathbf{x}_p$ (for $A\mathbf{x} = \mathbf{b}$).

**Every matrix $A$ with full row rank ($r = m$) has the following properties:**
1. All rows have pivots, and $R_0$ has no zero rows. $R_0 = R$.
2. $A\mathbf{x} = \mathbf{b}$ has a solution for $\forall \mathbf{b}$.
3. The column space of $A$ is the whole space $\mathbb{R}^m$.

$A\mathbf{x} = \mathbf{b} \to R\mathbf{x} = \mathbf{d}$ where $R = (I_{m \times m} \quad F_{m \times (n-m)})$.

4. If $m < n$, $A\mathbf{x} = \mathbf{b}$ is **underdetermined** (many solutions).

With full row rank ($r = m$), $m$ rows are linearly independent. The columns of $A^T$ are LI. The nullspace of $A^T$ is $\mathbb{Z} = \{\mathbf{0}\}$.

### 4.7 Four Possibilities for Linear Equations

Four possibilities for linear equations depend on the rank $r$:

| Case | Shape | $R_0$ form | Solutions to $A\mathbf{x} = \mathbf{b}$ |
|:-----|:------|:-----------|:----------------------------------------|
| $r = m$ and $r = n$ | Square, invertible | $(I)$ | 1 solution |
| $r = m$ and $r < n$ | Short, wide | $(I \quad F)$ | $\infty$ solutions |
| $r < m$ and $r = n$ | Tall, thin | $\begin{pmatrix} I \\ 0 \end{pmatrix}$ | 0 or 1 solution |
| $r < m$ and $r < n$ | Not full rank | $\begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix}$ | 0 or $\infty$ solutions |

---

<br>

## 5. Independence, Basis, and Dimension (3.4)

### 5.1 Independent Vectors

**(1) Independent vectors:** The only zero combination

$$c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \cdots + c_n \mathbf{v}_n = \mathbf{0}$$

has all $c_1 = c_2 = \cdots = c_n = 0$.

This implies that if at least one of the scalars is nonzero, let's say $c_1 \neq 0$, then we can write:

$$\mathbf{v}_1 = -\frac{c_2}{c_1}\mathbf{v}_2 - \frac{c_3}{c_1}\mathbf{v}_3 - \cdots - \frac{c_n}{c_1}\mathbf{v}_n$$

$\mathbf{v}_1$ is a linear combination of other vectors.

**(2) The vectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$ span the space $\mathbb{S}$ if $\mathbb{S}$ = all combinations of $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$.**

e.g., $\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$ $\implies$ $\begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}$. $\hat{i}, \hat{j}$ span the space $\mathbb{R}^2$.

**(3)** The vectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$ are a **basis** for $\mathbb{S}$ if they are **linearly independent** and they **span** $\mathbb{S}$.

Every vector in the space is a **unique** combination of the basis vectors.

e.g., $\mathbb{R}^2 \ni \begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}$

**(4)** The **dimension** of a vector space $\mathbb{S}$ is the number $n$ of vectors in every basis for $\mathbb{S}$.

Consider $A \in \mathbb{R}^{m \times n}$. There are $n$ columns, of which $r$ are independent, meaning the remaining $n - r$ columns are dependent. The dimension of $C(A)$ is $r$, which is the rank of $A$.

Four essential ideas in this section are:
1. Independent Vectors
2. Spanning a Space
3. Basis for a Space
4. Dimension of a Space

### 5.2 Linear Independence via the Nullspace

**Def.** The columns of $A$ are **linearly independent** when the only solution to $A\mathbf{x} = \mathbf{0}$ is $\mathbf{x} = \mathbf{0}$.

$$\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n = \mathbf{0}$$

iff $x_1 = x_2 = \cdots = x_n = 0$.

**Geometric interpretation:**

e.g., $\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3 \in \mathbb{R}^3$ are not in the same plane $\implies$ those vectors are independent. No combination of $\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$ that equals $\mathbf{0}$ exists (except the trivial one).

e.g., $\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3$ are in the same plane in $\mathbb{R}^3$. They are dependent. For instance, $\mathbf{w}_3 = \mathbf{w}_1 + \mathbf{w}_2$ $\iff$ $1 \cdot \mathbf{w}_1 + 1 \cdot \mathbf{w}_2 - 1 \cdot \mathbf{w}_3 = \mathbf{0}$. The combination gives $\mathbf{0}$, but there are nonzero coefficients, meaning they are dependent.

**Quick checks for dependence:**

- Q. Are $\begin{pmatrix} 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 10^{-5} \end{pmatrix}$ dependent? **No.**
- Q. Are $\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} -1 \\ -1 \end{pmatrix}$ dependent? **Yes.**
- Q. Are $\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 0 \end{pmatrix}$ dependent? **Yes.**
- Q. In $\mathbb{R}^2$, any three vectors are dependent. **True.**

**Ex 1.** Let $A = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix}$. Does $A$ have dependent columns?

Check the nullspace of $A$, $N(A)$.

$$A\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$\xrightarrow{R_2 - 2R_1, \; R_3 - R_1}$ $\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -1 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$, $\text{rank}(A) = 2 < 3$.

Take $x_3 = 1$: $x_1 + 3 = 0 \Rightarrow x_1 = -3$, $x_2 - 1 = 0 \Rightarrow x_2 = 1$.

The nonzero coefficients result in the zero vector:

$$-3\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + 1\begin{pmatrix} 3 \\ 5 \\ 3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

If the columns of $A$ are independent, then $\text{rank}(A) = r = n$ and the nullspace of $A$ has only the zero vector: $N(A) = \{\mathbf{0}\}$.

**Any set of $n$ vectors in $\mathbb{R}^m$ must be linearly dependent if $n > m$.**

e.g., Suppose you have 7 columns with 5 components ($5 \times 7$ matrix). 7 column vectors are from $\mathbb{R}^5$. There cannot be more than 5 pivots in 5 rows. $A\mathbf{x} = \mathbf{0}$ has at least $2 \;(= 7 - 5)$ free variables. That is, it has nonzero solutions.

**Note:** The columns might be dependent or independent if $n \leq m$. Elimination will reveal the $r$ pivots.

### 5.3 Vectors that Span a Subspace

Let $A$ be a matrix. $C(A)$ is the column space which consists of all combinations of $A\mathbf{x}$:

$$\begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{v}_1 + x_2 \mathbf{v}_2 + \cdots + x_n \mathbf{v}_n$$

**e.g.,** $A = \begin{pmatrix} 1 & 4 \\ 2 & 7 \\ 3 & 5 \end{pmatrix}$, $A^T = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 7 & 5 \end{pmatrix}$

$C(A)$ is the plane in $\mathbb{R}^3$. $C(A^T)$, row space of $A$, is $\mathbb{R}^2$.

For $A \in \mathbb{R}^{m \times n}$: the row vectors are in $\mathbb{R}^n$, the column vectors are in $\mathbb{R}^m$.

### 5.4 Basis for a Vector Space

**Def.** A **basis** for a vector space is a sequence of vectors with two properties:
> i) The basis vectors are **linearly independent**.
> ii) They **span** the space.

**e.g.,** $\begin{pmatrix} x \\ y \end{pmatrix} = x\begin{pmatrix} 1 \\ 0 \end{pmatrix} + y\begin{pmatrix} 0 \\ 1 \end{pmatrix}$ where $\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$.

i) $\hat{i}$ and $\hat{j}$ are LI.
ii) $\hat{i}, \hat{j}$ span $\mathbb{R}^2$.

The combination $\mathbf{x} = x\hat{i} + y\hat{j}$ is unique because $\hat{i}, \hat{j}$ are LI.

**There is one and only one way to write $\mathbf{v}$ as a combination of the basis vectors.**

**Proof.** Suppose $\mathbf{v} = a_1 \mathbf{v}_1 + \cdots + a_n \mathbf{v}_n$ and $\mathbf{v} = b_1 \mathbf{v}_1 + \cdots + b_n \mathbf{v}_n$.

By subtraction, we have the zero vector:

$$\mathbf{0} = (a_1 - b_1)\mathbf{v}_1 + \cdots + (a_n - b_n)\mathbf{v}_n$$

From the LI of $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$: $a_1 - b_1 = 0$, $a_2 - b_2 = 0$, $\ldots$, $a_n - b_n = 0$.

Hence $a_i = b_i$ for $i = 1, 2, \ldots, n$. $\square$

---

**Ex 3.** The columns of $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$ produce the "standard basis" for $\mathbb{R}^2$.

Let $\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$.

i) $\hat{i}, \hat{j}$ are LI.
ii) $\hat{i}, \hat{j}$ span $\mathbb{R}^2$.

$\therefore \hat{i}, \hat{j}$ are basis vectors in $\mathbb{R}^2$.

---

**Ex 4.** The columns of every invertible $n$ by $n$ matrix give a basis for $\mathbb{R}^n$.

**Invertible matrix example:**

$$A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = R$$

3 LI column vectors. $C(A) = \mathbb{R}^3$. 3 nonzero pivots. $\text{rank}(A) = 3$. $N(A) = \{\mathbf{0}\}$.

**Note that the basis is NOT unique.**

The vectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$ are a basis for $\mathbb{R}^n$ when they are the columns of an $n$ by $n$ invertible matrix. $\mathbb{R}^n$ has infinitely many different bases.

---

**Singular matrix example:**

$$B = \begin{pmatrix} 1 & 0 & 1 \\ 1 & 1 & 2 \\ 1 & 1 & 2 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 1 & 2 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 0 & 0 \end{pmatrix} = R_0$$

$C(A) \neq \mathbb{R}^3$. $\text{rank}(A) = 2 < 3$.

When columns are dependent, we only keep the **pivot columns**.

e.g., 1st, 2nd columns in $B$.

- Every set of independent vectors can be **extended** to a basis. (e.g., 1st, 2nd columns in $B$ span a plane in $\mathbb{R}^3$.)
- Every spanning set of vectors can be **reduced** to a basis.

---

**Ex 5.** $A = \begin{pmatrix} 2 & 4 \\ 3 & 6 \end{pmatrix}$

$\xrightarrow{R_2 - \frac{3}{2}R_1}$ $\begin{pmatrix} 2 & 4 \\ 0 & 0 \end{pmatrix}$ $\xrightarrow{R_1 / 2}$ $\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0$

$\text{rank}(A) = 1 < 2$. One pivot column, one pivot row.

---

**Ex 6.** $R_0 = \begin{pmatrix} 1 & 2 & 0 & 3 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix}$

1st, 3rd columns are pivot columns. $\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$ are basis vectors for $C(R_0)$.

$C(R_0)$ is the $xy$ plane in $\mathbb{R}^3$.

Also, 2nd and 3rd column vectors are a basis of $C(R_0)$.

---

**All bases for a vector space contain the same number of vectors.** The number of vectors in any and every basis is the "**dimension**" of the space.

### 5.5 Dimension of a Vector Space

**Dimension of a Vector Space:**

If $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m$ and $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n$ are both bases for the same vector space, then $m = n$.

**Proof.** Let $n > m$. Since $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m$ are a basis, each $\mathbf{w}_i$ for $i = 1, 2, \ldots, n$ must be a combination of the $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m$:

$$\mathbf{w}_1 = a_{11}\mathbf{v}_1 + a_{21}\mathbf{v}_2 + \cdots + a_{m1}\mathbf{v}_m$$
$$\mathbf{w}_2 = a_{12}\mathbf{v}_1 + a_{22}\mathbf{v}_2 + \cdots + a_{m2}\mathbf{v}_m$$
$$\vdots$$
$$\mathbf{w}_n = a_{1n}\mathbf{v}_1 + a_{2n}\mathbf{v}_2 + \cdots + a_{mn}\mathbf{v}_m$$

This leads to:

$$W = (\mathbf{w}_1 \quad \mathbf{w}_2 \quad \cdots \quad \mathbf{w}_n) = (\mathbf{v}_1 \quad \mathbf{v}_2 \quad \cdots \quad \mathbf{v}_m) \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} = VA$$

$A$ is $m$ by $n$ matrix (short and wide). Because $n > m$, $A\mathbf{x} = \mathbf{0}$ has nonzero solutions.

From $A\mathbf{x} = \mathbf{0}$: $VA\mathbf{x} = V\mathbf{0} = \mathbf{0}$.

That is $W\mathbf{x} = \mathbf{0}$, i.e., $x_1\mathbf{w}_1 + x_2\mathbf{w}_2 + \cdots + x_n\mathbf{w}_n = \mathbf{0}$.

Since $\mathbf{x}$ is a nonzero vector, $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n$ are not LI. $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n$ could NOT be a basis. This contradicts $\mathbf{w}_i$ for $i = 1, 2, \ldots, n$ is a basis. $\square$

---

**Def.** The **dimension** of a space is the number of vectors in every basis.

**e.g.,** The line through $\mathbf{u} = \begin{pmatrix} 1 \\ 5 \\ 2 \end{pmatrix}$ has 1 dimension.

Perpendicular to that line is the plane: $\mathbf{u} \cdot \mathbf{x} = 0 \iff x + 5y + 2z = 0$.

The plane null space of the matrix $A = \begin{pmatrix} 1 & 5 & 2 \end{pmatrix}$:

$(1 \quad 5 \quad 2)\begin{pmatrix} x \\ y \\ z \end{pmatrix} = 0$

$n = 3$, $\text{rank}(A) = r = 1$, $n - r = 2$ free variables.

$n - r$ special solutions give a basis for the nullspace: dimension $n - r$.

To find the special solutions:
- Take $y = 1, z = 0 \implies x = -5$
- Take $y = 0, z = 1 \implies x = -2$

$\begin{pmatrix} -5 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} -2 \\ 0 \\ 1 \end{pmatrix}$ are a basis (2 dimension).

**Summary of dimensions:**
- The row space has dimension $r$
- The column space has dimension $r$
- The nullspace has dimension $n - r$
- $N(A^T)$ has dimension $m - r$

### 5.6 Bases for Matrix Spaces and Function Spaces

**Bases for Matrix Spaces and Function Spaces**

The words "independence," "basis," "dimension" are not limited to column vectors.

**Matrix Spaces:**

The vector space $\mathbb{M}$ contains all $2 \times 2$ matrices. Its dimension is 4.

e.g., $A_1, A_2, A_3, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

i) They are LI: $x_1 A_1 + x_2 A_2 + x_3 A_3 + x_4 A_4 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$ iff $x_1 = x_2 = x_3 = x_4 = 0$.

ii) They span $\mathbb{M}$: $\begin{pmatrix} a & b \\ c & d \end{pmatrix} = aA_1 + bA_2 + cA_3 + dA_4$. The linear combinations of $A_1, A_2, A_3, A_4$ can produce any matrix in $\mathbb{M}$.

**Subspaces of matrix spaces:**

- $A_1, A_2, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$ are a basis for the subspace $\mathcal{U}$ of upper triangular matrices: $\mathcal{U} \ni \begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

- $A_1, A_4$ are a basis for the subspace $\mathbb{D}$ of diagonal matrices: $\mathbb{D} \ni \begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

- $A_1, A_4, A_2 + A_3$ are a basis for symmetric matrices: $\mathbb{S} \ni \begin{pmatrix} a & b \\ b & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

**Dimensions of matrix subspaces (for $n \times n$ matrices):**

| Space | Dimension |
|:------|:----------|
| Whole $n \times n$ matrix space | $n^2$ |
| Diagonal matrices | $n$ |
| Upper triangular matrices | $\frac{1}{2}n^2 + \frac{1}{2}n$ |
| Symmetric matrices | $\frac{1}{2}n^2 + \frac{1}{2}n$ |

For upper triangular matrices: Total number of zeros $= \sum_{i=1}^{n}(i-1) = \sum_{i=1}^{n} i - n = \frac{n(n+1)}{2} - n$. Total number of nonzeros $= n^2 - \frac{n(n+1)}{2} + n = \frac{n^2}{2} + \frac{n}{2}$.

The number of entries (dimension) for upper triangular and symmetric matrices is the same: $\frac{1}{2}n^2 + \frac{1}{2}n$.

---

**Function Spaces:**

$$\frac{d^2y}{dx^2} = 0, \quad \frac{d^2y}{dx^2} = -y, \quad \frac{d^2y}{dx^2} = y$$

These involve the 2nd derivative. In calculus, we find the solution $y(x)$:

| ODE | Solution | Basis | Dimension |
|:----|:---------|:------|:----------|
| $y'' = 0$ | $y = cx + d$ | $\{1, x\}$ | 2 |
| $y'' = -y$ | $y = c\sin x + d\cos x$ | $\{\sin x, \cos x\}$ | 2 |
| $y'' = y$ | $y = ce^x + de^{-x}$ | $\{e^x, e^{-x}\}$ | 2 |

The basis vectors are in the **nullspace** of the 2nd derivative.

$y'' = 2$ has the particular solution $y_p = x^2$: $\frac{dy_p}{dx} = 2x$, $\frac{d^2 y_p}{dx^2} = 2$.

Therefore, the general (complete) solution becomes:

$$y(x) = y_p(x) + y_n(x) = x^2 + cx + d$$

---

<br>

## 6. Dimensions of the Four Subspaces (3.5)

### 6.1 Dimension Summary

1. The column space $C(A)$ and the row space $C(A^T)$ both have dimension $r$ (= rank of $A$).
2. The nullspace of $A$, $N(A)$, has dimension $n - r$.
3. The left nullspace of $A$, $N(A^T)$, has dimension $m - r$.
4. Elimination from $A$ to $R_0$ **changes** $C(A)$ and $N(A^T)$, but their dimensions don't change.

### 6.2 Orthogonality of the Subspaces

We are going to connect "**rank**" and "**dimension**":
- **Rank** of a matrix counts independent columns.
- **Dimension** of a subspace is the number of vectors in a basis.

The rank of $A \in \mathbb{R}^{m \times n}$ reveals the dimension of all four fundamental subspaces:

1. **Row space**, $C(A^T)$, is a subspace of $\mathbb{R}^n$, dimension $r$. (Each row vector $\in \mathbb{R}^n$.)
2. **Column space**, $C(A)$, is a subspace of $\mathbb{R}^m$, dimension $r$. (Each column vector $\in \mathbb{R}^m$.)
3. **Nullspace**, $N(A)$, is a subspace of $\mathbb{R}^n$, dimension $n - r$. ($A\mathbf{x} = \mathbf{0}$, $\mathbf{x} \in \mathbb{R}^n$.)
4. **Left nullspace**, $N(A^T)$, is a subspace of $\mathbb{R}^m$, dimension $m - r$. ($A^T\mathbf{y} = \mathbf{0}$, $\mathbf{y} \in \mathbb{R}^m$.)

**Orthogonal pairs:**

$$C(A) \subset \mathbb{R}^m, \quad C(A^T) \subset \mathbb{R}^n$$
$$N(A^T) \subset \mathbb{R}^m, \quad N(A) \subset \mathbb{R}^n$$

- $C(A)$ and $C(A^T)$ have $r$ dimensions.
- Column space of $A$ and row space of $A$ have the same dimension $r$.
- $N(A)$ has dimension $n - r$.
- $N(A^T)$ has dimension $m - r$.

**$N(A)$ is perpendicular to $C(A^T)$ (row space of $A$):**

e.g., $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$

The row vectors $\{(1, 2), (3, 4)\}$ span the row space of $A$. The solutions to $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$ are in the nullspace of $A$, $N(A)$.

$N(A) \ni \mathbf{x}$ is perpendicular to any vector in the row space of $A$ in the sense of inner product.

That is, take $\mathbf{y} \in C(A^T)$, $\mathbf{x} \in N(A)$: $\mathbf{y} \cdot \mathbf{x} = 0$.

**$N(A^T)$ is perpendicular to $C(A)$ (column space of $A$):**

e.g., $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$

The column vectors $\{\begin{pmatrix} 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ 4 \end{pmatrix}\}$ span the column space of $A$.

$N(A^T) \ni \mathbf{y}$ implies that $A^T\mathbf{y} = \mathbf{0}$:

$\begin{pmatrix} 1 & 3 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$\Rightarrow C(A) \ni \mathbf{x} \perp \mathbf{y} \in N(A^T)$

i.e., $\alpha\begin{pmatrix} 1 \\ 3 \end{pmatrix} + \beta\begin{pmatrix} 2 \\ 4 \end{pmatrix} \perp \mathbf{y}$

### 6.3 The Four Subspaces for R_0

**Suppose $R_0 = \text{rref}(A)$.** The four dimensions are the same for $R_0$ and $A$.

**Example:** Consider a $3 \times 5$ matrix $R_0$:

$$R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}$$

Pivot rows: 1 and 2. Pivot columns: 1 and 4. $\text{rank}(R_0) = r = 2$.

**Row space:** Spanned by basis vectors $\{(1, 3, 5, 0, 7), \;(0, 0, 0, 1, 2)\}$. $\dim C(R_0^T) = 2 = r$.

**Column space:** The 1st and 4th column vectors form $\left\{\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\}$, a basis for $C(R_0)$. $\dim C(R_0) = 2 = r$.

**Nullspace:** $R_0 \mathbf{x} = \mathbf{0}$:

$$\begin{pmatrix} x_1 + 3x_2 + 5x_3 + 7x_5 \\ x_4 + 2x_5 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$n - r = 5 - 2 = 3$, 3 free variables.

i) Take $x_2 = 1, x_3 = 0, x_5 = 0$: $x_1 + 3 = 0 \Rightarrow x_1 = -3$, $x_4 = 0$. $\mathbf{x} = \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}$

ii) Take $x_2 = 0, x_3 = 1, x_5 = 0$: $x_1 + 5 = 0 \Rightarrow x_1 = -5$, $x_4 = 0$. $\mathbf{x} = \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}$

iii) Take $x_2 = 0, x_3 = 0, x_5 = 1$: $x_1 + 7 = 0 \Rightarrow x_1 = -7$, $x_4 + 2 = 0 \Rightarrow x_4 = -2$. $\mathbf{x} = \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}$

The 3 special solutions form a basis $\left\{\begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}\right\}$ for $N(R_0)$.

$\dim N(R_0) = 3 = 5 - 2 = n - r$.

**Left nullspace:** $R_0^T \mathbf{y} = \mathbf{0}$:

$$\begin{pmatrix} 1 & 0 & 0 \\ 3 & 0 & 0 \\ 5 & 0 & 0 \\ 0 & 1 & 0 \\ 7 & 2 & 0 \end{pmatrix} \begin{pmatrix} y_1 \\ y_2 \\ y_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix} \iff \begin{pmatrix} y_1 \\ 3y_1 \\ 5y_1 \\ y_2 \\ 7y_1 + 2y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}$$

$m - r = 3 - 2 = 1$, 1 free variable.

Take $y_3 = 1$: $\Rightarrow y_1 = y_2 = 0$. $\therefore \mathbf{y} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$

$\left\{\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\right\}$ forms a basis for $N(R_0^T)$.

$\dim N(R_0^T) = 1 = m - r = 3 - 2$.

$R_0^T \mathbf{y} = \mathbf{0} \iff \mathbf{y}^T R_0 = \mathbf{0}$. $\mathbf{y}^T$ is a row vector to the 'left' of $R_0$.

**Orthogonality summary:**

$$C(A) \subset \mathbb{R}^m \perp N(A^T) \subset \mathbb{R}^m$$
$$C(A^T) \subset \mathbb{R}^n \perp N(A) \subset \mathbb{R}^n$$

In $\mathbb{R}^n$: the row space and the null space have $r$ and $n - r$ dimensions.

In $\mathbb{R}^m$: the column space and the left nullspace have $r$ and $m - r$ dimensions.

### 6.4 Relationship Between A and R_0

**The four subspace dimensions for $A$ are the same as for $R_0$.**

$$A = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 1 & 3 & 5 & 1 & 9 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}$$

$A$ has the same row space as $R_0$, but its column space differs from that of $R_0$.

A basis for $C(A)$: $\left\{\begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}\right\}$ (pivot columns of $A$, **not** $R_0$)

A basis for $C(R_0)$: $\left\{\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\}$

(Different planes in $\mathbb{R}^3$, but same dimension.)

**Key relationships:**

1. $A$ has the **same row space** as $R_0$, $R$: $C(A^T) = C(R_0^T) = C(R^T)$, dimension $r$.

2. The column space of $A$ has dimension $r$: $\dim C(A) = \dim C(A^T)$. $C(A) \neq C(R_0)$, but $\dim C(A) = \dim C(R_0) = r$.

3. $A$ has the **same nullspace** as $R_0$: $A\mathbf{x} = \mathbf{0} \iff R_0 \mathbf{x} = \mathbf{0}$. **Elimination does NOT change the solution.** The special solutions form a basis. $n - r$ free variables $\implies$ $\dim N(A) = n - r$.

$$\dim C(A) + \dim N(A) = r + (n - r) = n$$

4. The left nullspace of $A$, $N(A^T)$: $\dim C(A^T) + \dim N(A^T) = r + (m - r) = m$.

### 6.5 The Fundamental Theorem of Linear Algebra

$$\boxed{\textbf{Fundamental Theorem of Linear Algebra}}$$

> i) $\dim C(A) = \dim C(A^T) = r$
>
> ii) $\dim N(A) = n - r$, $\quad \dim N(A^T) = m - r$

**The Four Subspaces for $A$:**

$$C(A^T) \subset \mathbb{R}^n \quad \perp \quad N(A) \subset \mathbb{R}^n$$
$$C(A) \subset \mathbb{R}^m \quad \perp \quad N(A^T) \subset \mathbb{R}^m$$

In $\mathbb{R}^n$: $C(A^T)$ (dim $r$) and $N(A)$ (dim $n - r$) are orthogonal complements.

In $\mathbb{R}^m$: $C(A)$ (dim $r$) and $N(A^T)$ (dim $m - r$) are orthogonal complements.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Vector Space | A set closed under addition and scalar multiplication, satisfying 8 axioms |
| Field | A set ($\mathbb{R}$, $\mathbb{C}$) where $+, -, \times, \div$ are defined |
| Subspace | A subset of a vector space that is itself a vector space (must contain $\mathbf{0}$) |
| Column Space $C(A)$ | All linear combinations of columns of $A$; $A\mathbf{x} = \mathbf{b}$ solvable iff $\mathbf{b} \in C(A)$ |
| Row Space $C(A^T)$ | Column space of $A^T$; spanned by rows of $A$ |
| Nullspace $N(A)$ | All solutions to $A\mathbf{x} = \mathbf{0}$; subspace of $\mathbb{R}^n$ |
| Left Nullspace $N(A^T)$ | All $\mathbf{y}$ with $A^T\mathbf{y} = \mathbf{0}$; subspace of $\mathbb{R}^m$ |
| RREF $R_0 = \text{rref}(A)$ | Reduced row echelon form; contains $I$ in pivot columns, $F$ in free columns |
| $A = CR$ | $C$ = independent columns of $A$; $R = (I \; F)$ = reduced row echelon form (no zero rows) |
| Special Solutions | Columns of $\begin{pmatrix} -F \\ I \end{pmatrix}$; form a basis for $N(A)$ |
| Complete Solution | $\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n$ (particular + nullspace) |
| Particular Solution $\mathbf{x}_p$ | Set all free variables to 0, solve $R_0\mathbf{x} = \mathbf{d}$ |
| Full Column Rank ($r = n$) | $N(A) = \{\mathbf{0}\}$; at most 1 solution to $A\mathbf{x} = \mathbf{b}$; $R_0 = \begin{pmatrix} I \\ 0 \end{pmatrix}$ |
| Full Row Rank ($r = m$) | $A\mathbf{x} = \mathbf{b}$ always solvable; $C(A) = \mathbb{R}^m$; $R = (I \; F)$ |
| Linear Independence | $c_1\mathbf{v}_1 + \cdots + c_n\mathbf{v}_n = \mathbf{0}$ only when all $c_i = 0$ |
| Spanning | Vectors span $\mathbb{S}$ if $\mathbb{S}$ = all combinations of those vectors |
| Basis | Linearly independent vectors that span the space; representation is unique |
| Dimension | Number of vectors in any basis; invariant across all bases |
| $\dim C(A) = \dim C(A^T) = r$ | Column and row space share the same dimension (rank) |
| $\dim N(A) = n - r$ | Nullspace dimension equals number of free variables |
| $\dim N(A^T) = m - r$ | Left nullspace dimension |
| $r + (n - r) = n$ | Row space + nullspace fill $\mathbb{R}^n$ |
| $r + (m - r) = m$ | Column space + left nullspace fill $\mathbb{R}^m$ |
| Orthogonality | $N(A) \perp C(A^T)$ in $\mathbb{R}^n$; $N(A^T) \perp C(A)$ in $\mathbb{R}^m$ |
| Matrix Space dim | $n \times n$: $n^2$; diagonal: $n$; upper triangular: $\frac{n^2+n}{2}$; symmetric: $\frac{n^2+n}{2}$ |
| Function Space | Solutions to $y'' = 0, -y, y$ form 2-dim spaces with bases $\{1,x\}$, $\{\sin x, \cos x\}$, $\{e^x, e^{-x}\}$ |
| Fundamental Theorem | $\dim C(A) = \dim C(A^T) = r$; $\dim N(A) = n-r$; $\dim N(A^T) = m-r$ |

---
