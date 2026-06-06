# Chapter 4 Lecture — Orthogonality

> **Last Updated:** 2026-06-06
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 4

> **Prerequisites**: [Linear Algebra] Vector spaces, basis, and dimension (Ch 1-3).
>
> **Learning Objectives**:
> 1. Define orthogonality of vectors and subspaces
> 2. Compute projections onto lines and subspaces
> 3. Apply Gram-Schmidt process to construct orthogonal bases
> 4. Solve least squares problems using projections

> **📑 This document is split into 2 parts.**
>
> **Part 1** · [Part 2](Concepts-Part2.md)

---

<br>

## Table of Contents

- [1. Orthogonality of Vectors and Subspaces (4.1)](#1-orthogonality-of-vectors-and-subspaces-41)
- [2. Projections onto Lines and Subspaces (4.2)](#2-projections-onto-lines-and-subspaces-42)
- [3. Least Square Approximations (4.3)](Concepts-Part2.md#3-least-square-approximations-43)
- [4. Orthogonal Bases and Gram-Schmidt (4.4)](Concepts-Part2.md#4-orthogonal-bases-and-gram-schmidt-44)
- [5. The Pseudoinverse of a Matrix (4.5)](Concepts-Part2.md#5-the-pseudoinverse-of-a-matrix-45)
- [Summary](Concepts-Part2.md#summary)

---

<br>

## 1. Orthogonality of Vectors and Subspaces (4.1)

### 1.1 Orthogonal Vectors

Orthogonal vectors satisfy:

```math
\mathbf{v}^T \mathbf{w} = 0
```

and the **Pythagorean theorem** in vector form:

```math
\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2
```

**Proof:**

```math
\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w})
```

```math
= \mathbf{v} \cdot \mathbf{v} + \mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{v} + \mathbf{w} \cdot \mathbf{w}
```

```math
= \|\mathbf{v}\|^2 + 2\mathbf{v} \cdot \mathbf{w} + \|\mathbf{w}\|^2
```

When $`\mathbf{v}^T \mathbf{w} = 0`$, the cross term $`2\mathbf{v} \cdot \mathbf{w} = 0`$, giving:

```math
\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2
```

This is analogous to $`a^2 + b^2 = c^2`$.

**Recall from Chapter 1:** The dot product connects to the angle between $`\mathbf{v}`$ and $`\mathbf{w}`$:

```math
\mathbf{v}^T \mathbf{w} = \|\mathbf{v}\| \, \|\mathbf{w}\| \cos\theta
```

When $`\theta = 90°`$, $`\mathbf{v}^T \mathbf{w} = 0`$.

The vectors $`\mathbf{v}`$, $`\mathbf{w}`$, $`\mathbf{v} + \mathbf{w}`$ produce a right triangle. From the Pythagorean theorem:

```math
\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2
```

### 1.2 The Fundamental Subspaces are Orthogonal

**(1)** The nullspace of $`A`$, $`\mathcal{N}(A)`$, contains all vectors orthogonal to the row space $`C(A^T)`$.

```math
A\mathbf{x} = \mathbf{0}
```

```math
\begin{pmatrix} \text{--- row}_1 \text{ of } A \text{ ---} \\ \text{--- row}_2 \text{ of } A \text{ ---} \\ \vdots \\ \text{--- row}_m \text{ of } A \text{ ---} \end{pmatrix} \begin{pmatrix} \\ \mathbf{x} \\ {} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ \vdots \\ 0 \end{pmatrix}
```

$`\mathbf{x}`$ is orthogonal to each row of $`A`$. Every row has zero dot product with $`\mathbf{x}`$.

**(2)** The nullspace of $`A^T`$, $`\mathcal{N}(A^T)`$, contains all vectors orthogonal to the column space $`C(A)`$.

```math
A^T \mathbf{y} = \mathbf{0} \quad \Longleftrightarrow \quad \mathbf{y}^T A = \mathbf{0}
```

```math
\mathbf{y}^T \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} = \mathbf{0}
```

$`\Rightarrow \mathbf{y}^T \mathbf{a}_1 = 0, \quad \mathbf{y}^T \mathbf{a}_2 = 0, \quad \ldots, \quad \mathbf{y}^T \mathbf{a}_n = 0`$

$`\Rightarrow \mathbf{y}^T (c_1 \mathbf{a}_1 + c_2 \mathbf{a}_2 + \cdots + c_n \mathbf{a}_n) = 0`$

$`\mathbf{y}`$ is orthogonal to each column of $`A`$.

### 1.3 Projection Preview

If $`\mathbf{b}`$ is outside the column space of $`A`$, find the closest point $`\mathbf{p}`$ that is inside.

```math
A\mathbf{x} \neq \mathbf{b}
```

The error $`\mathbf{e}`$ is perpendicular to $`C(A)`$.

The least squares equation:

```math
A^T A \hat{\mathbf{x}} = A^T \mathbf{b}
```

produces the closest $`\mathbf{p} = A\hat{\mathbf{x}}`$. This gives the **best** solution $`\hat{\mathbf{x}}`$ when $`A\mathbf{x} = \mathbf{b}`$ is unsolvable:

```math
\min \|A\hat{\mathbf{x}} - \mathbf{b}\|^2
```

This is **least squares**. When $`A^T A = I`$, the problem becomes easy. Constructing such a matrix $`A`$ is the topic of Section 4.4.

### 1.4 Orthogonal Subspaces

**Definition:** Subspaces $`V`$ and $`W`$ are **orthogonal** when:

```math
\mathbf{v}^T \mathbf{w} = 0 \quad \forall \, \mathbf{v} \in V, \, \mathbf{w} \in W
```

Consider a line (1-dimensional subspace) passing through a vertical plane (2D subspace). Every vector on the line is perpendicular to every vector in the plane.

**(3)** The row space of $`A`$ is orthogonal to the nullspace of $`A`$. The column space of $`A`$ is orthogonal to the nullspace of $`A^T`$ (= left nullspace of $`A`$).

**Proof that $`C(A^T) \perp \mathcal{N}(A)`$:**

$`\mathbf{x}`$ in the nullspace of $`A`$ is orthogonal to $`A^T\mathbf{y}`$ in the row space of $`A`$:

```math
\mathbf{x} \cdot (A^T \mathbf{y}) = (A^T \mathbf{y}) \cdot \mathbf{x} = (A^T \mathbf{y})^T \mathbf{x} = \mathbf{y}^T A \mathbf{x} = \mathbf{y}^T (A\mathbf{x}) = \mathbf{y}^T \mathbf{0} = 0
```

Also, from $`A\mathbf{x} = \mathbf{0}`$:

```math
(\text{row}_1) \cdot \mathbf{x} = 0, \quad (\text{row}_2) \cdot \mathbf{x} = 0, \quad \ldots, \quad (\text{row}_m) \cdot \mathbf{x} = 0
```

```math
\Rightarrow [c_1(\text{row}_1) + c_2(\text{row}_2) + \cdots + c_m(\text{row}_m)] \cdot \mathbf{x} = 0
```

Therefore $`C(A^T)`$ is perpendicular to $`\mathcal{N}(A)`$.

**Proof that $`C(A) \perp \mathcal{N}(A^T)`$:**

Similarly, $`\mathbf{y}`$ in the nullspace of $`A^T`$ is orthogonal to $`A\mathbf{x}`$ in the column space of $`A`$:

```math
\mathbf{y} \cdot (A\mathbf{x}) = (A\mathbf{x}) \cdot \mathbf{y} = (A\mathbf{x})^T \mathbf{y} = \mathbf{x}^T A^T \mathbf{y} = \mathbf{x}^T (A^T \mathbf{y}) = \mathbf{x}^T \mathbf{0} = 0
```

### 1.5 Orthogonal Complements and Dimension

**(4)** The dimensions add to:

```math
r + (n - r) = n \quad \text{and} \quad r + (m - r) = m
```

These are **orthogonal complements**.

**Example:**

```math
A = \begin{pmatrix} 1 & -2 & 1 \\ 1 & 0 & -1 \end{pmatrix}
```

Row reduction: $`R_2 - R_1`$, then $`R_2/2`$, then $`R_1 + 2R_2`$:

```math
R = \begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & -1 \end{pmatrix}
```

$`\text{rank}(A) = \text{rank}(R) = 2`$, $`n - r = 3 - 2 = 1`$ (1 free variable).

```math
\mathcal{N}(A) = \left\lbrace c\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} \right\rbrace, \quad C(A) = \mathbb{R}^2
```

For $`A^T`$:

```math
A^T = \begin{pmatrix} 1 & 1 \\ -2 & 0 \\ 1 & -1 \end{pmatrix}
```

Row reduction: $`R_2 + 2R_1`$, $`R_3 - R_1`$, then $`R_3 + R_2`$, $`R_2/2`$, then $`R_1 - R_2`$:

```math
R_0 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}
```

$`\text{rank}(A^T) = \text{rank}(R_0) = 2`$, $`m - r = 2 - 2 = 0`$ (no free variable).

```math
\mathcal{N}(A^T) = \left\lbrace \begin{pmatrix} 0 \\ 0 \end{pmatrix} \right\rbrace, \quad C(A^T) \text{ is a plane in } \mathbb{R}^3
```

Verification:

```math
\dim C(A^T) + \dim \mathcal{N}(A) = 2 + 1 = 3 = n
```

```math
\dim C(A) + \dim \mathcal{N}(A^T) = 2 + 0 = 2 = m
```

**Important restriction:** If $`V`$ and $`W`$ are orthogonal subspaces in $`\mathbb{R}^n`$, then:

```math
\dim V + \dim W \leq n
```

Two orthogonal subspaces that account for the whole space have a special name: **"orthogonal complements"**. The orthogonal complement $`V^\perp`$ of $`V`$ contains all vectors orthogonal to $`V`$.

**Orthogonal complement pairs:**

| Pair | Dimensions | Space |
|:-----|:-----------|:------|
| Row space and nullspace | $`r + (n-r) = n`$ | $`\mathbb{R}^n`$ |
| Column space and left nullspace | $`r + (m-r) = m`$ | $`\mathbb{R}^m`$ |

```math
\mathcal{N}(A) \text{ is the orthogonal complement of row space } C(A^T) \text{ in } \mathbb{R}^n
```

```math
\mathcal{N}(A^T) \text{ is the orthogonal complement of column space } C(A) \text{ in } \mathbb{R}^m
```

Every $`\mathbf{x}`$ can be split into a row space component $`\mathbf{x}_r`$ and a null space component $`\mathbf{x}_n`$:

```math
\mathbb{R}^n \ni \mathbf{x} = \mathbf{x}_r + \mathbf{x}_n
```

```math
\mathbb{R}^m \ni \mathbf{y} = \mathbf{y}_{\text{col}} + \mathbf{y}_{\text{left null}}
```

### 1.6 The Big Picture of Four Subspaces

Multiply $`A`$ to $`\mathbf{x}`$:

```math
A\mathbf{x} = A(\mathbf{x}_r + \mathbf{x}_n) = A\mathbf{x}_r + A\mathbf{x}_n = \mathbf{b}
```

since $`A\mathbf{x}_n = \mathbf{0}`$.

- $`A\mathbf{x}_r = \mathbf{b}`$ is in the column space of $`A`$.
- $`A\mathbf{x}_n = \mathbf{0}`$.

The complete solution to $`A\mathbf{x} = \mathbf{b}`$ is:

```math
\mathbf{x} = \mathbf{x}_r + \mathbf{x}_n
```

where $`\mathbf{x}_r`$ is the **unique** row space component.

### 1.7 Minimum Norm Solution

```math
\|\mathbf{x}\|^2 = \mathbf{x} \cdot \mathbf{x} = (\mathbf{x}_r + \mathbf{x}_n) \cdot (\mathbf{x}_r + \mathbf{x}_n)
```

```math
= \mathbf{x}_r \cdot \mathbf{x}_r + \mathbf{x}_r \cdot \mathbf{x}_n + \mathbf{x}_n \cdot \mathbf{x}_r + \mathbf{x}_n \cdot \mathbf{x}_n
```

```math
= \|\mathbf{x}_r\|^2 + \|\mathbf{x}_n\|^2
```

(since the cross terms are zero because $`\mathbf{x}_r \perp \mathbf{x}_n`$).

Therefore the **minimum norm solution** to $`A\mathbf{x} = \mathbf{b}`$ is $`\mathbf{x} = \mathbf{x}_r`$ and $`\mathbf{x}_n = \mathbf{0}`$.

**Uniqueness:** Every vector $`\mathbf{b} \in C(A)`$ comes from exactly **one** vector $`\mathbf{x}_r`$ in the row space.

**Proof:** Let $`A\mathbf{x}_r = A\mathbf{x}_r' = \mathbf{b}`$. Then $`A\mathbf{x}_r - A\mathbf{x}_r' = A(\mathbf{x}_r - \mathbf{x}_r') = \mathbf{0}`$. So $`\mathbf{x}_r - \mathbf{x}_r' \in \mathcal{N}(A)`$. Both $`\mathbf{x}_r`$ and $`\mathbf{x}_r'`$ come from $`C(A^T)`$, so $`\mathbf{x}_r - \mathbf{x}_r' \in C(A^T)`$. Because $`\mathbf{0} \in C(A^T)`$ and $`\mathbf{0} \in \mathcal{N}(A)`$, we have $`\mathbf{x}_r - \mathbf{x}_r' = \mathbf{0}`$, i.e., $`\mathbf{x}_r = \mathbf{x}_r'`$. $`\square`$

### 1.8 Example: Invertible Submatrix

**Example 2:** Every matrix of rank $`r`$ has an $`r`$ by $`r`$ invertible submatrix.

```math
A = \begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 1 & 2 & 4 & 5 & 6 \\ 1 & 2 & 4 & 5 & 6 \end{pmatrix}
```

Row reduction: $`R_2 - R_1`$, $`R_3 - R_1`$:

```math
\begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 & 1 \end{pmatrix}
```

$`R_3 - R_2`$:

```math
\begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}
```

$`R_1 - 3R_2`$:

```math
\begin{pmatrix} 1 & 2 & 0 & 1 & 2 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}, \quad \text{rank}(A) = 2
```

$`A`$ contains the $`2 \times 2`$ invertible submatrix $`\begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}`$ (from pivot columns 1 and 3).

### 1.9 Combining Bases from Subspaces

A **basis** contains linearly independent vectors that span the space.

**Standard basis** of $`\mathbb{R}^n`$ is $`\lbrace\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\rbrace`$ where:

```math
\mathbb{R}^n \ni \mathbf{e}_i = \begin{pmatrix} 0 \\ \vdots \\ 0 \\ 1 \\ 0 \\ \vdots \\ 0 \end{pmatrix} \leftarrow i\text{th row}
```

is the vector in $`\mathbb{R}^n`$ which has a one in the $`i`$th entry and zeros elsewhere. That is the $`i`$th column of $`I \in \mathbb{R}^{n \times n}`$.

The dimension of $`\mathbb{R}^n`$ is $`n`$ because the number of basis vectors of $`\mathbb{R}^n`$ is $`n`$, e.g., $`\lbrace\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\rbrace`$.

### 1.10 Two Properties in $`\mathbb{R}^n`$

**i)** Suppose $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is LI. Then $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is a basis for $`\mathbb{R}^n`$.

**ii)** Suppose $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m\rbrace`$ spans $`\mathbb{R}^n`$. Then $`m \geq n`$. If $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ spans $`\mathbb{R}^n`$, then $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is LI.

```math
\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace \text{ is LI} \iff \lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace \text{ spans } \mathbb{R}^n
```

Therefore $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is a basis for $`\mathbb{R}^n`$.

**Proof:**

**i)** $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is LI. Need to show the set spans $`\mathbb{R}^n`$.

Define $`A = (\mathbf{u}_1 \; \mathbf{u}_2 \; \cdots \; \mathbf{u}_n) \in \mathbb{R}^{n \times n}`$. Our goal is to find a unique $`\mathbf{x}`$ such that $`A\mathbf{x} = \mathbf{v}`$. Since the square matrix $`A`$ has full rank, it has an inverse: $`A^{-1}A\mathbf{x} = A^{-1}\mathbf{v}`$, so $`\mathbf{x} = A^{-1}\mathbf{v}`$. Therefore $`\mathbf{v} = A\mathbf{x}`$ is a linear combination of $`\mathbf{u}_i`$.

**ii)** Since $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m\rbrace`$ spans $`\mathbb{R}^n`$: $`\mathbf{v} \in \mathbb{R}^n = c_1\mathbf{u}_1 + c_2\mathbf{u}_2 + \cdots + c_m\mathbf{u}_m`$. Let $`A = (\mathbf{u}_1 \; \mathbf{u}_2 \; \cdots \; \mathbf{u}_m)`$ which is $`n \times m`$. The reduced echelon form of $`A`$ reveals $`r \leq m`$ independent column vectors. $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_r\rbrace`$ would be a basis for $`\mathbb{R}^n`$, but this is contrary to $`\dim \mathbb{R}^n = n`$. Therefore $`m \geq n`$.

Finally suppose $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is NOT LI. Then $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_r\rbrace`$ would form a basis for $`\mathbb{R}^n`$ with $`r < n`$, which contradicts $`\dim \mathbb{R}^n = n`$. Therefore $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$ is LI. $`\square`$

### 1.11 Properties of $`A \in \mathbb{R}^{n \times n}`$

**i)** If the $`n`$ columns of $`A`$ are LI, they span $`\mathbb{R}^n`$.

**ii)** If the $`n`$ columns span $`\mathbb{R}^n`$, they are LI (see the proof above).

For a square full rank matrix $`A \in \mathbb{R}^{n \times n}`$:
- $`A\mathbf{x} = \mathbf{b}`$ is solvable (existence)
- $`\mathbf{x} = A^{-1}\mathbf{b}`$ is unique (uniqueness)

**iii)** $`AB = I`$ for full rank $`A, B \in \mathbb{R}^{n \times n}`$. Then $`BA = I`$.

**Proof:** $`AB = I`$. Then $`B(AB) = BI = B`$ (associative law). $`(BA)B = IB`$ (distributive law). $`(BA - I)B = 0`$, therefore $`BA = I`$.

### 1.12 Example 3: Four Subspaces

```math
A = \begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}
```

$`R_2 - 3R_1`$:

```math
\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0, \quad \text{rank}(A) = \text{rank}(R_0) = 1 = r
```

$`n = 2`$, $`n - r = 2 - 1 = 1`$ (1 free variable).

```math
C(A^T) = \text{span}\lbrace(1, 2)\rbrace, \quad C(A) = \text{span}\left\lbrace\begin{pmatrix} 1 \\ 3 \end{pmatrix}\right\rbrace
```

$`A\mathbf{x} = \mathbf{0}`$: choose $`x_2 = 1`$, $`x_1 + 2 = 0`$, $`x_1 = -2`$.

```math
\mathcal{N}(A) = \text{span}\left\lbrace\begin{pmatrix} -2 \\ 1 \end{pmatrix}\right\rbrace
```

$`A^T\mathbf{y} = \mathbf{0}`$: $`\begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$. After reduction: $`\begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}`$. Take $`y_2 = 1`$, $`y_1 = -3`$.

```math
\mathcal{N}(A^T) = \text{span}\left\lbrace\begin{pmatrix} -3 \\ 1 \end{pmatrix}\right\rbrace
```

Suppose $`\mathbf{b} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}`$. The solution $`\mathbf{x}`$ to $`A\mathbf{x} = \mathbf{b}`$ is:

```math
\mathbf{x} = \mathbf{x}_p + c\,\mathbf{x}_n
```

```math
\begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}
```

Take $`x_2 = 0`$, $`x_1 = 10`$. So $`\mathbf{x}_p = \begin{pmatrix} 10 \\ 0 \end{pmatrix}`$.

```math
\mathbf{x} = \begin{pmatrix} 10 \\ 0 \end{pmatrix} + c\begin{pmatrix} -2 \\ 1 \end{pmatrix}
```

When $`c = 1`$: $`\begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}\begin{pmatrix} 8 \\ 1 \end{pmatrix} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}`$. Verified.

**Note:** $`\mathbf{x}_p`$ is **not** orthogonal to $`\mathbf{x}_n`$: $`(10 \quad 0)\begin{pmatrix} -2 \\ 1 \end{pmatrix} = -20 \neq 0`$.

We can further decompose $`\mathbf{x}_p`$ into $`\mathbf{x}_r`$ and $`\mathbf{x}_n`$:

```math
\begin{pmatrix} 10 \\ 0 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix} - 4\begin{pmatrix} -2 \\ 1 \end{pmatrix}
```

where $`\mathbf{x}_r = \begin{pmatrix} 2 \\ 4 \end{pmatrix}`$ and the null component is $`-4\begin{pmatrix} -2 \\ 1 \end{pmatrix}`$.

Then $`\mathbf{x} = \mathbf{x}_r + (c - 4)\mathbf{x}_n = \mathbf{x}_r + c'\mathbf{x}_n`$.

Check: $`(2 \quad 4)\begin{pmatrix} -2 \\ 1 \end{pmatrix} = -4 + 4 = 0`$, so $`\mathbf{x}_r \perp \mathbf{x}_n`$.

---

<br>

## 2. Projections onto Lines and Subspaces (4.2)

### 2.1 Key Facts Summary

**(1)** The projection of $`\mathbf{b}`$ onto the line through $`\mathbf{a}`$ is the closest point to $`\mathbf{b}`$:

```math
\mathbf{p} = \mathbf{a}\frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}
```

**(2)** The error $`\mathbf{e} = \mathbf{b} - \mathbf{p}`$ is perpendicular to $`\mathbf{a}`$. Right triangle $`\mathbf{b}`$, $`\mathbf{p}`$, $`\mathbf{e}`$ has:

```math
\|\mathbf{p}\|^2 + \|\mathbf{e}\|^2 = \|\mathbf{b}\|^2
```

**(3)** The projection of $`\mathbf{b}`$ onto a subspace $`S`$ is the closest vector $`\mathbf{p}`$ in $`S`$; $`\mathbf{b} - \mathbf{p}`$ is orthogonal to $`S`$.

**(4)** $`A^T A`$ is invertible (and symmetric) when $`A`$ has independent columns: $`\mathcal{N}(A^T A) = \mathcal{N}(A)`$.

**(5)** Then the projection of $`\mathbf{b}`$ onto the column space $`C(A)`$ is:

```math
\mathbf{p} = A(A^T A)^{-1} A^T \mathbf{b}
```

**(6)** The projection matrix onto $`C(A)`$ is:

```math
P = A(A^T A)^{-1} A^T
```

It has $`\mathbf{p} = P\mathbf{b}`$ and $`P^2 = P = P^T`$.

### 2.2 Projection onto a Line

A line goes through the origin in the direction of $`\mathbf{a}`$. We project $`\mathbf{b}`$ onto the line. The line from $`\mathbf{b}`$ to $`\mathbf{p}`$ is perpendicular to the vector $`\mathbf{a}`$:

```math
\mathbf{e} = \mathbf{b} - \mathbf{p} \perp \mathbf{a}
```

The projection $`\mathbf{p}`$ is a multiple of $`\mathbf{a}`$: $`\mathbf{p} = \alpha \mathbf{a}`$.

**Derivation:**

```math
\mathbf{e} = \mathbf{b} - \mathbf{p} = \mathbf{b} - \alpha\mathbf{a}
```

```math
\mathbf{e} \cdot \mathbf{a} = (\mathbf{b} - \alpha\mathbf{a}) \cdot \mathbf{a} = \mathbf{b} \cdot \mathbf{a} - \alpha\,\mathbf{a} \cdot \mathbf{a} = \mathbf{a}^T\mathbf{b} - \alpha\,\mathbf{a}^T\mathbf{a} = 0
```

```math
\therefore \alpha = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}
```

```math
\mathbf{p} = \alpha\mathbf{a} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\,\mathbf{a}
```

**Special cases:**
- If $`\mathbf{b} = \mathbf{a}`$, then $`\alpha = 1`$, $`\mathbf{p} = \mathbf{a}`$.
- If $`\mathbf{b} \perp \mathbf{a}`$ (i.e., $`\mathbf{a} \cdot \mathbf{b} = 0`$), then $`\alpha = 0`$, $`\mathbf{p} = \mathbf{0}`$.

**Example 1:** Project $`\mathbf{b} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$ onto $`\mathbf{a} = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}`$, to find $`\mathbf{p} = \alpha\mathbf{a}`$.

```math
\mathbf{a}^T\mathbf{a} = (1\;2\;2)\begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix} = 1 + 4 + 4 = 9
```

```math
\mathbf{a}^T\mathbf{b} = (1\;2\;2)\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = 1 + 2 + 2 = 5
```

```math
\alpha = \frac{5}{9}
```

```math
\mathbf{p} = \frac{5}{9}\mathbf{a} = \begin{pmatrix} 5/9 \\ 10/9 \\ 10/9 \end{pmatrix}, \quad \mathbf{e} = \mathbf{b} - \mathbf{p} = \begin{pmatrix} 4/9 \\ -1/9 \\ -1/9 \end{pmatrix}
```

$`\mathbf{b}`$ is split into two parts: $`\mathbf{b} = \mathbf{p} + \mathbf{e}`$, where $`\mathbf{p} \perp \mathbf{e}`$.

```math
\|\mathbf{p}\| = \|\mathbf{b}\|\cos\theta, \quad \|\mathbf{e}\| = \|\mathbf{b}\|\sin\theta
```

### 2.3 Projection Matrix P (onto a Line)

```math
\mathbf{p} = \alpha\mathbf{a} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\,\mathbf{a} = \mathbf{a}\frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}} = \left(\frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}\right)\mathbf{b}
```

```math
P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}
```

Note: $`\mathbf{a}\mathbf{a}^T`$ is a **rank 1** matrix (column times row). We are projecting onto a one-dimensional subspace, the line through $`\mathbf{a}`$, which is $`C(P)`$.

**Example:** Find the projection matrix $`P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}`$ onto a line through $`\mathbf{a} = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}`$.

```math
\mathbf{a}^T\mathbf{a} = 9
```

```math
\mathbf{a}\mathbf{a}^T = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}(1\;2\;2) = \begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}
```

```math
P = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}
```

What happens if $`\mathbf{a} = \begin{pmatrix} 2 \\ 4 \\ 4 \end{pmatrix}`$? Then $`\mathbf{a}^T\mathbf{a} = 4(1\;2\;2)\begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix} = 36`$, and $`\mathbf{a}\mathbf{a}^T = 4\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}`$. So $`P = \frac{1}{4 \cdot 9} \cdot 4\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}`$ — **same!**

**Verification of $`P^2 = P`$:**

```math
P^2 = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}\frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = \frac{1}{81}\begin{pmatrix} 9 & 18 & 18 \\ 18 & 36 & 36 \\ 18 & 36 & 36 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = P
```

**Trace:** $`\text{diag}(P) \cdot \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \frac{1}{9}(1\;4\;4)\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \frac{1}{9} \cdot 9 = 1`$.

**Complementary projection:** When $`P`$ projects onto one subspace, $`I - P`$ projects onto the perpendicular subspace (orthogonal complement). $`I - P`$ projects onto the plane perpendicular to $`\mathbf{a}`$.

### 2.4 Projection in $`\mathbb{R}^3`$ Example

Consider $`\mathbf{b} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix} \in \mathbb{R}^3`$.

- $`\mathbf{p}_1 = P_1\mathbf{b} = \begin{pmatrix} 0 \\ 0 \\ 4 \end{pmatrix}`$ is the projection of $`\mathbf{b}`$ onto the $`z`$-axis.
- $`\mathbf{p}_2 = P_2\mathbf{b} = \begin{pmatrix} 2 \\ 3 \\ 0 \end{pmatrix}`$ is the projection of $`\mathbf{b}`$ onto the $`xy`$-plane.

```math
P_1 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad P_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
```

**Observations:**
- $`P_1 + P_2 = I_{3 \times 3}`$
- $`P_1 P_2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}`$ (zero matrix)
- $`P_1, P_2`$ are perpendicular; $`xy`$-plane and $`z`$-axis are orthogonal subspaces
- The line and the plane are orthogonal complements: $`\dim(\text{line}) + \dim(\text{plane}) = 1 + 2 = 3`$
- Every vector $`\mathbf{b}`$ is the sum of its parts in the two subspaces:

```math
\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x \\ y \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \\ z \end{pmatrix} = \mathbf{p}_2 + \mathbf{p}_1
```

### 2.5 Projection onto a Subspace

Every subspace of $`\mathbb{R}^m`$ has its own $`m`$ by $`m`$ projection matrix $`P`$. The projection matrix $`P`$ produces the part: $`\mathbf{p} = P\mathbf{b}`$.

A subspace is constructed by a basis. For example:

```math
A_1 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}, \quad A_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}, \quad A_3 = \begin{pmatrix} 2 & 3 \\ 2 & 3 \\ 0 & 0 \end{pmatrix}
```

$`C(A_1)`$ is the $`z`$-axis, $`C(A_2)`$ is the $`xy`$-plane, $`C(A_3)`$ is the $`xy`$-plane.

Start with LI $`n`$-vectors $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3, \ldots, \mathbf{a}_n \in \mathbb{R}^m`$. Find the combination:

```math
\mathbf{p} = \alpha_1\mathbf{a}_1 + \alpha_2\mathbf{a}_2 + \cdots + \alpha_n\mathbf{a}_n
```

closest to a given vector $`\mathbf{b}`$. We are projecting each $`\mathbf{b}`$ in $`\mathbb{R}^m`$ onto the $`n`$-dimensional subspace spanned by the $`\mathbf{a}`$'s:

```math
A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)
```

$`C(A)`$ is a subspace of $`\mathbb{R}^m`$. $`A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n \in C(A)`$.

We are looking for the particular combination $`\mathbf{p} = A\hat{\boldsymbol{\alpha}}`$ that is closest to $`\mathbf{b}`$. $`\hat{\boldsymbol{\alpha}}`$ is the best vector in $`C(A)`$.

When $`n = 1`$: $`\alpha = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}`$.

For $`n > 1`$: $`\hat{\boldsymbol{\alpha}} = \begin{pmatrix} \alpha_1 \\ \alpha_2 \\ \vdots \\ \alpha_n \end{pmatrix}`$ is to be found.

### 2.6 Derivation for $`n = 2`$

Let $`n = 2`$, $`A = (\mathbf{a}_1 \; \mathbf{a}_2)`$.

The subspace $`S`$ is spanned by $`\mathbf{a}_1, \mathbf{a}_2`$:

```math
\mathbf{p} = \alpha_1\mathbf{a}_1 + \alpha_2\mathbf{a}_2 = A\hat{\boldsymbol{\alpha}} \in C(A)
```

The error vector $`\mathbf{e} = \mathbf{b} - \mathbf{p}`$ is perpendicular to the subspace $`S`$:

```math
\mathbf{a}_1 \cdot \mathbf{e} = 0, \quad \mathbf{a}_2 \cdot \mathbf{e} = 0
```

```math
\mathbf{a}_1^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \mathbf{a}_2^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0
```

```math
\begin{pmatrix} \text{---} \; \mathbf{a}_1^T \; \text{---} \\ \text{---} \; \mathbf{a}_2^T \; \text{---} \end{pmatrix}(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

```math
A^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \mathbf{0}
```

```math
A^T\mathbf{b} = A^T A\hat{\boldsymbol{\alpha}}
```

If $`A^T A`$ is invertible:

```math
\hat{\boldsymbol{\alpha}} = (A^T A)^{-1} A^T \mathbf{b}
```

(When $`n = 1`$: $`\alpha = (\mathbf{a}^T\mathbf{a})^{-1}\mathbf{a}^T\mathbf{b}`$.)

```math
\mathbf{p} = A\hat{\boldsymbol{\alpha}} = A(A^T A)^{-1} A^T \mathbf{b}
```

```math
P = A(A^T A)^{-1} A^T
```

### 2.7 Extension to $`n`$-dimensional Subspace

We can easily extend this to $`n`$-dimensional subspace, where we have $`n`$ equations:

```math
\mathbf{a}_1^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \mathbf{a}_2^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \ldots, \quad \mathbf{a}_n^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0
```

```math
A^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \mathbf{0}
```

**Remarks:**
- $`\mathbf{e} = \mathbf{b} - A\hat{\boldsymbol{\alpha}}`$ is in $`\mathcal{N}(A^T)`$, which is perpendicular to $`C(A)`$.
- Left nullspace of $`A`$ contains the error vector.
- $`\mathbf{b}`$ is split into the projection $`\mathbf{p}`$ and the error $`\mathbf{e}`$.

### 2.8 Example 2: Projection onto a Subspace

```math
A = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix}
```

Find $`\hat{\boldsymbol{\alpha}}`$, $`\mathbf{p}`$, and $`P`$.

```math
A^T A = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix} = \begin{pmatrix} 3 & 3 \\ 3 & 5 \end{pmatrix}
```

```math
A^T\mathbf{b} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{pmatrix}\begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} 6 \\ 0 \end{pmatrix}
```

```math
\hat{\boldsymbol{\alpha}} = (A^T A)^{-1} A^T\mathbf{b} = \frac{1}{15 - 9}\begin{pmatrix} 5 & -3 \\ -3 & 3 \end{pmatrix}\begin{pmatrix} 6 \\ 0 \end{pmatrix} = \begin{pmatrix} 5 & -3 \\ -3 & 3 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 5 \\ -3 \end{pmatrix}
```

```math
\mathbf{p} = A\hat{\boldsymbol{\alpha}} = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} 5 \\ -3 \end{pmatrix} = 5\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} - 3\begin{pmatrix} 0 \\ 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 5 \\ 2 \\ -1 \end{pmatrix}
```

```math
P = A(A^T A)^{-1} A^T = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}\frac{1}{6}\begin{pmatrix} 5 & -3 \\ -3 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 5 & 2 & -1 \\ 2 & 2 & 2 \\ -1 & 2 & 5 \end{pmatrix}
```

### 2.9 $`A^T A`$ is Invertible iff $`A`$ has LI Columns

**Note:** $`(A^T A)^{-1} \neq A^{-1}(A^T)^{-1}`$ because $`A^{-1}`$ does NOT exist (when $`A`$ is not square).

**Theorem:** $`A^T A`$ is invertible if and only if $`A`$ has linearly independent columns.

**Proof:**

Let $`A \in \mathbb{R}^{m \times n}`$.

($`\Rightarrow`$) $`A\mathbf{x} = \mathbf{0}`$ implies $`\mathbf{x} \in \mathcal{N}(A)`$. Multiplying by $`A^T`$ gives $`A^T A\mathbf{x} = \mathbf{0}`$, meaning $`\mathbf{x} \in \mathcal{N}(A^T A)`$. That is, $`\mathcal{N}(A) \ni \mathbf{x} \longrightarrow \mathbf{x} \in \mathcal{N}(A^T A)`$.

($`\Leftarrow`$) From $`A^T A\mathbf{x} = \mathbf{0}`$, multiply by $`\mathbf{x}^T`$:

```math
\mathbf{x}^T(A^T A\mathbf{x}) = \mathbf{x}^T\mathbf{0}
```

```math
(\mathbf{x}^T A^T)(A\mathbf{x}) = 0
```

```math
(A\mathbf{x})^T(A\mathbf{x}) = 0 \iff \|A\mathbf{x}\|^2 = 0
```

```math
\therefore A\mathbf{x} = \mathbf{0} \longrightarrow \mathbf{x} \in \mathcal{N}(A)
```

Therefore $`\mathcal{N}(A) = \lbrace\mathbf{0}\rbrace \iff \mathcal{N}(A^T A) = \lbrace\mathbf{0}\rbrace`$.

$`\iff`$ $`A`$ is invertible $`\iff`$ $`A^T A`$ is invertible. $`\square`$

**Size considerations:** $`A \in \mathbb{R}^{m \times n}`$, $`A^T \in \mathbb{R}^{n \times m}`$, $`A^T A \in \mathbb{R}^{n \times n}`$.

$`A^T A`$ is **symmetric**.

**Counter-example:** $`A = \begin{pmatrix} 1 & 2 \\ 1 & 2 \\ 0 & 0 \end{pmatrix}`$ has a dependent column.

```math
A^T A = \begin{pmatrix} 1 & 1 & 0 \\ 2 & 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 1 & 2 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \end{pmatrix}
```

$`\det(A^T A) = 2 \cdot 8 - 4 \cdot 4 = 0`$. $`A^T A`$ is NOT invertible, $`A^T A`$ is singular.

### 2.10 Worked Example (4.2A)

**Problem:** Project the vector $`\mathbf{b} = (3, 4, 4)`$ onto the line through $`\mathbf{a} = (2, 2, 1)`$ and then onto the plane that also contains $`\mathbf{a}^* = (1, 0, 0)`$. Check that the first error vector $`\mathbf{b} - \mathbf{p}`$ is perpendicular to $`\mathbf{a}`$, and the second error vector $`\mathbf{e}^* = \mathbf{b} - \mathbf{p}^*`$ is also perpendicular to $`\mathbf{a}`$ and $`\mathbf{a}^*`$. Find the $`3 \times 3`$ projection matrix $`P`$ onto that plane of $`\mathbf{a}`$ and $`\mathbf{a}^*`$. Find a vector whose projection onto the plane is $`\mathbf{p} = 0`$. Note $`P^2 = P = P^T`$.

**Solution:**

**Onto a line:**

```math
\mathbf{p} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\mathbf{a} = \frac{18}{9}(2, 2, 1) = (4, 4, 2) = 2\mathbf{a}
```

The error vector $`\mathbf{e} = \mathbf{b} - \mathbf{p} = (-1, 0, 2)`$ is perpendicular to $`\mathbf{a} = (2, 2, 1)`$. Check: $`\mathbf{e}^T\mathbf{a} = -2 + 0 + 2 = 0`$. So $`\mathbf{p}`$ is correct.

**Onto the plane:** The plane of $`\mathbf{a} = (2, 2, 1)`$ and $`\mathbf{a}^* = (1, 0, 0)`$ is the column space of $`A = [\mathbf{a} \;\; \mathbf{a}^*]`$:

```math
A = \begin{pmatrix} 2 & 1 \\ 2 & 0 \\ 1 & 0 \end{pmatrix}, \quad A^T A = \begin{pmatrix} 9 & 2 \\ 2 & 1 \end{pmatrix}, \quad (A^T A)^{-1} = \frac{1}{5}\begin{pmatrix} 1 & -2 \\ -2 & 9 \end{pmatrix}
```

```math
P = A(A^T A)^{-1}A^T = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0.8 & 0.4 \\ 0 & 0.4 & 0.2 \end{pmatrix}
```

Now $`\mathbf{p}^* = P\mathbf{b} = (3, 4.8, 2.4)`$. The error $`\mathbf{e}^* = \mathbf{b} - \mathbf{p}^* = (0, -0.8, 1.6)`$ is perpendicular to $`\mathbf{a}`$ and $`\mathbf{a}^*`$.

**Verification:**

```math
(\mathbf{e}^*)^T\mathbf{a} = \begin{pmatrix} 0 \\ -0.8 \\ 1.6 \end{pmatrix}^T\begin{pmatrix} 2 \\ 2 \\ 1 \end{pmatrix} = -1.6 + 1.6 = 0
```

```math
(\mathbf{e}^*)^T\mathbf{a}^* = \begin{pmatrix} 0 \\ -0.8 \\ 1.6 \end{pmatrix}^T\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} = 0 + 0 + 0 = 0
```

Since $`A^T\mathbf{e}^* = \mathbf{0}`$, we have $`A(A^T A)^{-1}A^T\mathbf{e}^* = \mathbf{0}`$, i.e., $`P\mathbf{e}^* = \mathbf{0}`$, so $`\mathbf{e}^* \in \mathcal{N}(P)`$.

---

<br>
