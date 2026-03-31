# Chapter 4 Lecture — Orthogonality

> **Last Updated:** 2026-03-31

---

<br>

## Table of Contents

- [1. Orthogonality of Vectors and Subspaces (4.1)](#1-orthogonality-of-vectors-and-subspaces-41)
- [2. Projections onto Lines and Subspaces (4.2)](#2-projections-onto-lines-and-subspaces-42)
- [3. Least Square Approximations (4.3)](#3-least-square-approximations-43)
- [4. Orthogonal Bases and Gram-Schmidt (4.4)](#4-orthogonal-bases-and-gram-schmidt-44)
- [5. The Pseudoinverse of a Matrix (4.5)](#5-the-pseudoinverse-of-a-matrix-45)
- [Summary](#summary)

---

<br>

## 1. Orthogonality of Vectors and Subspaces (4.1)

### 1.1 Orthogonal Vectors

Orthogonal vectors satisfy:

$$\mathbf{v}^T \mathbf{w} = 0$$

and the **Pythagorean theorem** in vector form:

$$\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2$$

**Proof:**

$$\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w})$$

$$= \mathbf{v} \cdot \mathbf{v} + \mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{v} + \mathbf{w} \cdot \mathbf{w}$$

$$= \|\mathbf{v}\|^2 + 2\mathbf{v} \cdot \mathbf{w} + \|\mathbf{w}\|^2$$

When $\mathbf{v}^T \mathbf{w} = 0$, the cross term $2\mathbf{v} \cdot \mathbf{w} = 0$, giving:

$$\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2$$

This is analogous to $a^2 + b^2 = c^2$.

**Recall from Chapter 1:** The dot product connects to the angle between $\mathbf{v}$ and $\mathbf{w}$:

$$\mathbf{v}^T \mathbf{w} = \|\mathbf{v}\| \, \|\mathbf{w}\| \cos\theta$$

When $\theta = 90°$, $\mathbf{v}^T \mathbf{w} = 0$.

The vectors $\mathbf{v}$, $\mathbf{w}$, $\mathbf{v} + \mathbf{w}$ produce a right triangle. From the Pythagorean theorem:

$$\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2$$

### 1.2 The Fundamental Subspaces are Orthogonal

**(1)** The nullspace of $A$, $\mathcal{N}(A)$, contains all vectors orthogonal to the row space $C(A^T)$.

$$A\mathbf{x} = \mathbf{0}$$

$$\begin{pmatrix} \text{--- row}_1 \text{ of } A \text{ ---} \\ \text{--- row}_2 \text{ of } A \text{ ---} \\ \vdots \\ \text{--- row}_m \text{ of } A \text{ ---} \end{pmatrix} \begin{pmatrix} \\ \mathbf{x} \\ \phantom{x} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ \vdots \\ 0 \end{pmatrix}$$

$\mathbf{x}$ is orthogonal to each row of $A$. Every row has zero dot product with $\mathbf{x}$.

**(2)** The nullspace of $A^T$, $\mathcal{N}(A^T)$, contains all vectors orthogonal to the column space $C(A)$.

$$A^T \mathbf{y} = \mathbf{0} \quad \Longleftrightarrow \quad \mathbf{y}^T A = \mathbf{0}$$

$$\mathbf{y}^T \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} = \mathbf{0}$$

$\Rightarrow \mathbf{y}^T \mathbf{a}_1 = 0, \quad \mathbf{y}^T \mathbf{a}_2 = 0, \quad \ldots, \quad \mathbf{y}^T \mathbf{a}_n = 0$

$\Rightarrow \mathbf{y}^T (c_1 \mathbf{a}_1 + c_2 \mathbf{a}_2 + \cdots + c_n \mathbf{a}_n) = 0$

$\mathbf{y}$ is orthogonal to each column of $A$.

### 1.3 Projection Preview

If $\mathbf{b}$ is outside the column space of $A$, find the closest point $\mathbf{p}$ that is inside.

$$A\mathbf{x} \neq \mathbf{b}$$

The error $\mathbf{e}$ is perpendicular to $C(A)$.

The least squares equation:

$$A^T A \hat{\mathbf{x}} = A^T \mathbf{b}$$

produces the closest $\mathbf{p} = A\hat{\mathbf{x}}$. This gives the **best** solution $\hat{\mathbf{x}}$ when $A\mathbf{x} = \mathbf{b}$ is unsolvable:

$$\min \|A\hat{\mathbf{x}} - \mathbf{b}\|^2$$

This is **least squares**. When $A^T A = I$, the problem becomes easy. Constructing such a matrix $A$ is the topic of Section 4.4.

### 1.4 Orthogonal Subspaces

**Definition:** Subspaces $V$ and $W$ are **orthogonal** when:

$$\mathbf{v}^T \mathbf{w} = 0 \quad \forall \, \mathbf{v} \in V, \, \mathbf{w} \in W$$

Consider a line (1-dimensional subspace) passing through a vertical plane (2D subspace). Every vector on the line is perpendicular to every vector in the plane.

**(3)** The row space of $A$ is orthogonal to the nullspace of $A$. The column space of $A$ is orthogonal to the nullspace of $A^T$ (= left nullspace of $A$).

**Proof that $C(A^T) \perp \mathcal{N}(A)$:**

$\mathbf{x}$ in the nullspace of $A$ is orthogonal to $A^T\mathbf{y}$ in the row space of $A$:

$$\mathbf{x} \cdot (A^T \mathbf{y}) = (A^T \mathbf{y}) \cdot \mathbf{x} = (A^T \mathbf{y})^T \mathbf{x} = \mathbf{y}^T A \mathbf{x} = \mathbf{y}^T (A\mathbf{x}) = \mathbf{y}^T \mathbf{0} = 0$$

Also, from $A\mathbf{x} = \mathbf{0}$:

$$(\text{row}_1) \cdot \mathbf{x} = 0, \quad (\text{row}_2) \cdot \mathbf{x} = 0, \quad \ldots, \quad (\text{row}_m) \cdot \mathbf{x} = 0$$

$$\Rightarrow [c_1(\text{row}_1) + c_2(\text{row}_2) + \cdots + c_m(\text{row}_m)] \cdot \mathbf{x} = 0$$

Therefore $C(A^T)$ is perpendicular to $\mathcal{N}(A)$.

**Proof that $C(A) \perp \mathcal{N}(A^T)$:**

Similarly, $\mathbf{y}$ in the nullspace of $A^T$ is orthogonal to $A\mathbf{x}$ in the column space of $A$:

$$\mathbf{y} \cdot (A\mathbf{x}) = (A\mathbf{x}) \cdot \mathbf{y} = (A\mathbf{x})^T \mathbf{y} = \mathbf{x}^T A^T \mathbf{y} = \mathbf{x}^T (A^T \mathbf{y}) = \mathbf{x}^T \mathbf{0} = 0$$

### 1.5 Orthogonal Complements and Dimension

**(4)** The dimensions add to:

$$r + (n - r) = n \quad \text{and} \quad r + (m - r) = m$$

These are **orthogonal complements**.

**Example:**

$$A = \begin{pmatrix} 1 & -2 & 1 \\ 1 & 0 & -1 \end{pmatrix}$$

Row reduction: $R_2 - R_1$, then $R_2/2$, then $R_1 + 2R_2$:

$$R = \begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & -1 \end{pmatrix}$$

$\text{rank}(A) = \text{rank}(R) = 2$, $n - r = 3 - 2 = 1$ (1 free variable).

$$\mathcal{N}(A) = \left\{ c\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} \right\}, \quad C(A) = \mathbb{R}^2$$

For $A^T$:

$$A^T = \begin{pmatrix} 1 & 1 \\ -2 & 0 \\ 1 & -1 \end{pmatrix}$$

Row reduction: $R_2 + 2R_1$, $R_3 - R_1$, then $R_3 + R_2$, $R_2/2$, then $R_1 - R_2$:

$$R_0 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$$

$\text{rank}(A^T) = \text{rank}(R_0) = 2$, $m - r = 2 - 2 = 0$ (no free variable).

$$\mathcal{N}(A^T) = \left\{ \begin{pmatrix} 0 \\ 0 \end{pmatrix} \right\}, \quad C(A^T) \text{ is a plane in } \mathbb{R}^3$$

Verification:

$$\dim C(A^T) + \dim \mathcal{N}(A) = 2 + 1 = 3 = n$$

$$\dim C(A) + \dim \mathcal{N}(A^T) = 2 + 0 = 2 = m$$

**Important restriction:** If $V$ and $W$ are orthogonal subspaces in $\mathbb{R}^n$, then:

$$\dim V + \dim W \leq n$$

Two orthogonal subspaces that account for the whole space have a special name: **"orthogonal complements"**. The orthogonal complement $V^\perp$ of $V$ contains all vectors orthogonal to $V$.

**Orthogonal complement pairs:**

| Pair | Dimensions | Space |
|:-----|:-----------|:------|
| Row space and nullspace | $r + (n-r) = n$ | $\mathbb{R}^n$ |
| Column space and left nullspace | $r + (m-r) = m$ | $\mathbb{R}^m$ |

$$\mathcal{N}(A) \text{ is the orthogonal complement of row space } C(A^T) \text{ in } \mathbb{R}^n$$

$$\mathcal{N}(A^T) \text{ is the orthogonal complement of column space } C(A) \text{ in } \mathbb{R}^m$$

Every $\mathbf{x}$ can be split into a row space component $\mathbf{x}_r$ and a null space component $\mathbf{x}_n$:

$$\mathbb{R}^n \ni \mathbf{x} = \mathbf{x}_r + \mathbf{x}_n$$

$$\mathbb{R}^m \ni \mathbf{y} = \mathbf{y}_{\text{col}} + \mathbf{y}_{\text{left null}}$$

### 1.6 The Big Picture of Four Subspaces

Multiply $A$ to $\mathbf{x}$:

$$A\mathbf{x} = A(\mathbf{x}_r + \mathbf{x}_n) = A\mathbf{x}_r + A\mathbf{x}_n = \mathbf{b}$$

since $A\mathbf{x}_n = \mathbf{0}$.

- $A\mathbf{x}_r = \mathbf{b}$ is in the column space of $A$.
- $A\mathbf{x}_n = \mathbf{0}$.

The complete solution to $A\mathbf{x} = \mathbf{b}$ is:

$$\mathbf{x} = \mathbf{x}_r + \mathbf{x}_n$$

where $\mathbf{x}_r$ is the **unique** row space component.

### 1.7 Minimum Norm Solution

$$\|\mathbf{x}\|^2 = \mathbf{x} \cdot \mathbf{x} = (\mathbf{x}_r + \mathbf{x}_n) \cdot (\mathbf{x}_r + \mathbf{x}_n)$$

$$= \mathbf{x}_r \cdot \mathbf{x}_r + \mathbf{x}_r \cdot \mathbf{x}_n + \mathbf{x}_n \cdot \mathbf{x}_r + \mathbf{x}_n \cdot \mathbf{x}_n$$

$$= \|\mathbf{x}_r\|^2 + \|\mathbf{x}_n\|^2$$

(since the cross terms are zero because $\mathbf{x}_r \perp \mathbf{x}_n$).

Therefore the **minimum norm solution** to $A\mathbf{x} = \mathbf{b}$ is $\mathbf{x} = \mathbf{x}_r$ and $\mathbf{x}_n = \mathbf{0}$.

**Uniqueness:** Every vector $\mathbf{b} \in C(A)$ comes from exactly **one** vector $\mathbf{x}_r$ in the row space.

**Proof:** Let $A\mathbf{x}_r = A\mathbf{x}_r' = \mathbf{b}$. Then $A\mathbf{x}_r - A\mathbf{x}_r' = A(\mathbf{x}_r - \mathbf{x}_r') = \mathbf{0}$. So $\mathbf{x}_r - \mathbf{x}_r' \in \mathcal{N}(A)$. Both $\mathbf{x}_r$ and $\mathbf{x}_r'$ come from $C(A^T)$, so $\mathbf{x}_r - \mathbf{x}_r' \in C(A^T)$. Because $\mathbf{0} \in C(A^T)$ and $\mathbf{0} \in \mathcal{N}(A)$, we have $\mathbf{x}_r - \mathbf{x}_r' = \mathbf{0}$, i.e., $\mathbf{x}_r = \mathbf{x}_r'$. $\square$

### 1.8 Example: Invertible Submatrix

**Example 2:** Every matrix of rank $r$ has an $r$ by $r$ invertible submatrix.

$$A = \begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 1 & 2 & 4 & 5 & 6 \\ 1 & 2 & 4 & 5 & 6 \end{pmatrix}$$

Row reduction: $R_2 - R_1$, $R_3 - R_1$:

$$\begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 & 1 \end{pmatrix}$$

$R_3 - R_2$:

$$\begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}$$

$R_1 - 3R_2$:

$$\begin{pmatrix} 1 & 2 & 0 & 1 & 2 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}, \quad \text{rank}(A) = 2$$

$A$ contains the $2 \times 2$ invertible submatrix $\begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}$ (from pivot columns 1 and 3).

### 1.9 Combining Bases from Subspaces

A **basis** contains linearly independent vectors that span the space.

**Standard basis** of $\mathbb{R}^n$ is $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$ where:

$$\mathbb{R}^n \ni \mathbf{e}_i = \begin{pmatrix} 0 \\ \vdots \\ 0 \\ 1 \\ 0 \\ \vdots \\ 0 \end{pmatrix} \leftarrow i\text{th row}$$

is the vector in $\mathbb{R}^n$ which has a one in the $i$th entry and zeros elsewhere. That is the $i$th column of $I \in \mathbb{R}^{n \times n}$.

The dimension of $\mathbb{R}^n$ is $n$ because the number of basis vectors of $\mathbb{R}^n$ is $n$, e.g., $\{\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\}$.

### 1.10 Two Properties in $\mathbb{R}^n$

**i)** Suppose $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is LI. Then $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is a basis for $\mathbb{R}^n$.

**ii)** Suppose $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m\}$ spans $\mathbb{R}^n$. Then $m \geq n$. If $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ spans $\mathbb{R}^n$, then $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is LI.

$$\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\} \text{ is LI} \iff \{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\} \text{ spans } \mathbb{R}^n$$

Therefore $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is a basis for $\mathbb{R}^n$.

**Proof:**

**i)** $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is LI. Need to show the set spans $\mathbb{R}^n$.

Define $A = (\mathbf{u}_1 \; \mathbf{u}_2 \; \cdots \; \mathbf{u}_n) \in \mathbb{R}^{n \times n}$. Our goal is to find a unique $\mathbf{x}$ such that $A\mathbf{x} = \mathbf{v}$. Since the square matrix $A$ has full rank, it has an inverse: $A^{-1}A\mathbf{x} = A^{-1}\mathbf{v}$, so $\mathbf{x} = A^{-1}\mathbf{v}$. Therefore $\mathbf{v} = A\mathbf{x}$ is a linear combination of $\mathbf{u}_i$.

**ii)** Since $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m\}$ spans $\mathbb{R}^n$: $\mathbf{v} \in \mathbb{R}^n = c_1\mathbf{u}_1 + c_2\mathbf{u}_2 + \cdots + c_m\mathbf{u}_m$. Let $A = (\mathbf{u}_1 \; \mathbf{u}_2 \; \cdots \; \mathbf{u}_m)$ which is $n \times m$. The reduced echelon form of $A$ reveals $r \leq m$ independent column vectors. $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_r\}$ would be a basis for $\mathbb{R}^n$, but this is contrary to $\dim \mathbb{R}^n = n$. Therefore $m \geq n$.

Finally suppose $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is NOT LI. Then $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_r\}$ would form a basis for $\mathbb{R}^n$ with $r < n$, which contradicts $\dim \mathbb{R}^n = n$. Therefore $\{\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\}$ is LI. $\square$

### 1.11 Properties of $A \in \mathbb{R}^{n \times n}$

**i)** If the $n$ columns of $A$ are LI, they span $\mathbb{R}^n$.

**ii)** If the $n$ columns span $\mathbb{R}^n$, they are LI (see the proof above).

For a square full rank matrix $A \in \mathbb{R}^{n \times n}$:
- $A\mathbf{x} = \mathbf{b}$ is solvable (existence)
- $\mathbf{x} = A^{-1}\mathbf{b}$ is unique (uniqueness)

**iii)** $AB = I$ for full rank $A, B \in \mathbb{R}^{n \times n}$. Then $BA = I$.

**Proof:** $AB = I$. Then $B(AB) = BI = B$ (associative law). $(BA)B = IB$ (distributive law). $(BA - I)B = 0$, therefore $BA = I$.

### 1.12 Example 3: Four Subspaces

$$A = \begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}$$

$R_2 - 3R_1$:

$$\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0, \quad \text{rank}(A) = \text{rank}(R_0) = 1 = r$$

$n = 2$, $n - r = 2 - 1 = 1$ (1 free variable).

$$C(A^T) = \text{span}\{(1, 2)\}, \quad C(A) = \text{span}\left\{\begin{pmatrix} 1 \\ 3 \end{pmatrix}\right\}$$

$A\mathbf{x} = \mathbf{0}$: choose $x_2 = 1$, $x_1 + 2 = 0$, $x_1 = -2$.

$$\mathcal{N}(A) = \text{span}\left\{\begin{pmatrix} -2 \\ 1 \end{pmatrix}\right\}$$

$A^T\mathbf{y} = \mathbf{0}$: $\begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$. After reduction: $\begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}$. Take $y_2 = 1$, $y_1 = -3$.

$$\mathcal{N}(A^T) = \text{span}\left\{\begin{pmatrix} -3 \\ 1 \end{pmatrix}\right\}$$

Suppose $\mathbf{b} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}$. The solution $\mathbf{x}$ to $A\mathbf{x} = \mathbf{b}$ is:

$$\mathbf{x} = \mathbf{x}_p + c\,\mathbf{x}_n$$

$$\begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}$$

Take $x_2 = 0$, $x_1 = 10$. So $\mathbf{x}_p = \begin{pmatrix} 10 \\ 0 \end{pmatrix}$.

$$\mathbf{x} = \begin{pmatrix} 10 \\ 0 \end{pmatrix} + c\begin{pmatrix} -2 \\ 1 \end{pmatrix}$$

When $c = 1$: $\begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}\begin{pmatrix} 8 \\ 1 \end{pmatrix} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}$. Verified.

**Note:** $\mathbf{x}_p$ is **not** orthogonal to $\mathbf{x}_n$: $(10 \quad 0)\begin{pmatrix} -2 \\ 1 \end{pmatrix} = -20 \neq 0$.

We can further decompose $\mathbf{x}_p$ into $\mathbf{x}_r$ and $\mathbf{x}_n$:

$$\begin{pmatrix} 10 \\ 0 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix} - 4\begin{pmatrix} -2 \\ 1 \end{pmatrix}$$

where $\mathbf{x}_r = \begin{pmatrix} 2 \\ 4 \end{pmatrix}$ and the null component is $-4\begin{pmatrix} -2 \\ 1 \end{pmatrix}$.

Then $\mathbf{x} = \mathbf{x}_r + (c - 4)\mathbf{x}_n = \mathbf{x}_r + c'\mathbf{x}_n$.

Check: $(2 \quad 4)\begin{pmatrix} -2 \\ 1 \end{pmatrix} = -4 + 4 = 0$, so $\mathbf{x}_r \perp \mathbf{x}_n$.

---

<br>

## 2. Projections onto Lines and Subspaces (4.2)

### 2.1 Key Facts Summary

**(1)** The projection of $\mathbf{b}$ onto the line through $\mathbf{a}$ is the closest point to $\mathbf{b}$:

$$\mathbf{p} = \mathbf{a}\frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}$$

**(2)** The error $\mathbf{e} = \mathbf{b} - \mathbf{p}$ is perpendicular to $\mathbf{a}$. Right triangle $\mathbf{b}$, $\mathbf{p}$, $\mathbf{e}$ has:

$$\|\mathbf{p}\|^2 + \|\mathbf{e}\|^2 = \|\mathbf{b}\|^2$$

**(3)** The projection of $\mathbf{b}$ onto a subspace $S$ is the closest vector $\mathbf{p}$ in $S$; $\mathbf{b} - \mathbf{p}$ is orthogonal to $S$.

**(4)** $A^T A$ is invertible (and symmetric) when $A$ has independent columns: $\mathcal{N}(A^T A) = \mathcal{N}(A)$.

**(5)** Then the projection of $\mathbf{b}$ onto the column space $C(A)$ is:

$$\mathbf{p} = A(A^T A)^{-1} A^T \mathbf{b}$$

**(6)** The projection matrix onto $C(A)$ is:

$$P = A(A^T A)^{-1} A^T$$

It has $\mathbf{p} = P\mathbf{b}$ and $P^2 = P = P^T$.

### 2.2 Projection onto a Line

A line goes through the origin in the direction of $\mathbf{a}$. We project $\mathbf{b}$ onto the line. The line from $\mathbf{b}$ to $\mathbf{p}$ is perpendicular to the vector $\mathbf{a}$:

$$\mathbf{e} = \mathbf{b} - \mathbf{p} \perp \mathbf{a}$$

The projection $\mathbf{p}$ is a multiple of $\mathbf{a}$: $\mathbf{p} = \alpha \mathbf{a}$.

**Derivation:**

$$\mathbf{e} = \mathbf{b} - \mathbf{p} = \mathbf{b} - \alpha\mathbf{a}$$

$$\mathbf{e} \cdot \mathbf{a} = (\mathbf{b} - \alpha\mathbf{a}) \cdot \mathbf{a} = \mathbf{b} \cdot \mathbf{a} - \alpha\,\mathbf{a} \cdot \mathbf{a} = \mathbf{a}^T\mathbf{b} - \alpha\,\mathbf{a}^T\mathbf{a} = 0$$

$$\therefore \alpha = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}$$

$$\mathbf{p} = \alpha\mathbf{a} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\,\mathbf{a}$$

**Special cases:**
- If $\mathbf{b} = \mathbf{a}$, then $\alpha = 1$, $\mathbf{p} = \mathbf{a}$.
- If $\mathbf{b} \perp \mathbf{a}$ (i.e., $\mathbf{a} \cdot \mathbf{b} = 0$), then $\alpha = 0$, $\mathbf{p} = \mathbf{0}$.

**Example 1:** Project $\mathbf{b} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$ onto $\mathbf{a} = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}$, to find $\mathbf{p} = \alpha\mathbf{a}$.

$$\mathbf{a}^T\mathbf{a} = (1\;2\;2)\begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix} = 1 + 4 + 4 = 9$$

$$\mathbf{a}^T\mathbf{b} = (1\;2\;2)\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = 1 + 2 + 2 = 5$$

$$\alpha = \frac{5}{9}$$

$$\mathbf{p} = \frac{5}{9}\mathbf{a} = \begin{pmatrix} 5/9 \\ 10/9 \\ 10/9 \end{pmatrix}, \quad \mathbf{e} = \mathbf{b} - \mathbf{p} = \begin{pmatrix} 4/9 \\ -1/9 \\ -1/9 \end{pmatrix}$$

$\mathbf{b}$ is split into two parts: $\mathbf{b} = \mathbf{p} + \mathbf{e}$, where $\mathbf{p} \perp \mathbf{e}$.

$$\|\mathbf{p}\| = \|\mathbf{b}\|\cos\theta, \quad \|\mathbf{e}\| = \|\mathbf{b}\|\sin\theta$$

### 2.3 Projection Matrix P (onto a Line)

$$\mathbf{p} = \alpha\mathbf{a} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\,\mathbf{a} = \mathbf{a}\frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}} = \left(\frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}\right)\mathbf{b}$$

$$P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}$$

Note: $\mathbf{a}\mathbf{a}^T$ is a **rank 1** matrix (column times row). We are projecting onto a one-dimensional subspace, the line through $\mathbf{a}$, which is $C(P)$.

**Example:** Find the projection matrix $P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}$ onto a line through $\mathbf{a} = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}$.

$$\mathbf{a}^T\mathbf{a} = 9$$

$$\mathbf{a}\mathbf{a}^T = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}(1\;2\;2) = \begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}$$

$$P = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}$$

What happens if $\mathbf{a} = \begin{pmatrix} 2 \\ 4 \\ 4 \end{pmatrix}$? Then $\mathbf{a}^T\mathbf{a} = 4(1\;2\;2)\begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix} = 36$, and $\mathbf{a}\mathbf{a}^T = 4\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}$. So $P = \frac{1}{4 \cdot 9} \cdot 4\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}$ — **same!**

**Verification of $P^2 = P$:**

$$P^2 = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}\frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = \frac{1}{81}\begin{pmatrix} 9 & 18 & 18 \\ 18 & 36 & 36 \\ 18 & 36 & 36 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = P$$

**Trace:** $\text{diag}(P) \cdot \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \frac{1}{9}(1\;4\;4)\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \frac{1}{9} \cdot 9 = 1$.

**Complementary projection:** When $P$ projects onto one subspace, $I - P$ projects onto the perpendicular subspace (orthogonal complement). $I - P$ projects onto the plane perpendicular to $\mathbf{a}$.

### 2.4 Projection in $\mathbb{R}^3$ Example

Consider $\mathbf{b} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix} \in \mathbb{R}^3$.

- $\mathbf{p}_1 = P_1\mathbf{b} = \begin{pmatrix} 0 \\ 0 \\ 4 \end{pmatrix}$ is the projection of $\mathbf{b}$ onto the $z$-axis.
- $\mathbf{p}_2 = P_2\mathbf{b} = \begin{pmatrix} 2 \\ 3 \\ 0 \end{pmatrix}$ is the projection of $\mathbf{b}$ onto the $xy$-plane.

$$P_1 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad P_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

**Observations:**
- $P_1 + P_2 = I_{3 \times 3}$
- $P_1 P_2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}$ (zero matrix)
- $P_1, P_2$ are perpendicular; $xy$-plane and $z$-axis are orthogonal subspaces
- The line and the plane are orthogonal complements: $\dim(\text{line}) + \dim(\text{plane}) = 1 + 2 = 3$
- Every vector $\mathbf{b}$ is the sum of its parts in the two subspaces:

$$\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x \\ y \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \\ z \end{pmatrix} = \mathbf{p}_2 + \mathbf{p}_1$$

### 2.5 Projection onto a Subspace

Every subspace of $\mathbb{R}^m$ has its own $m$ by $m$ projection matrix $P$. The projection matrix $P$ produces the part: $\mathbf{p} = P\mathbf{b}$.

A subspace is constructed by a basis. For example:

$$A_1 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}, \quad A_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}, \quad A_3 = \begin{pmatrix} 2 & 3 \\ 2 & 3 \\ 0 & 0 \end{pmatrix}$$

$C(A_1)$ is the $z$-axis, $C(A_2)$ is the $xy$-plane, $C(A_3)$ is the $xy$-plane.

Start with LI $n$-vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3, \ldots, \mathbf{a}_n \in \mathbb{R}^m$. Find the combination:

$$\mathbf{p} = \alpha_1\mathbf{a}_1 + \alpha_2\mathbf{a}_2 + \cdots + \alpha_n\mathbf{a}_n$$

closest to a given vector $\mathbf{b}$. We are projecting each $\mathbf{b}$ in $\mathbb{R}^m$ onto the $n$-dimensional subspace spanned by the $\mathbf{a}$'s:

$$A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)$$

$C(A)$ is a subspace of $\mathbb{R}^m$. $A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n \in C(A)$.

We are looking for the particular combination $\mathbf{p} = A\hat{\boldsymbol{\alpha}}$ that is closest to $\mathbf{b}$. $\hat{\boldsymbol{\alpha}}$ is the best vector in $C(A)$.

When $n = 1$: $\alpha = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}$.

For $n > 1$: $\hat{\boldsymbol{\alpha}} = \begin{pmatrix} \alpha_1 \\ \alpha_2 \\ \vdots \\ \alpha_n \end{pmatrix}$ is to be found.

### 2.6 Derivation for $n = 2$

Let $n = 2$, $A = (\mathbf{a}_1 \; \mathbf{a}_2)$.

The subspace $S$ is spanned by $\mathbf{a}_1, \mathbf{a}_2$:

$$\mathbf{p} = \alpha_1\mathbf{a}_1 + \alpha_2\mathbf{a}_2 = A\hat{\boldsymbol{\alpha}} \in C(A)$$

The error vector $\mathbf{e} = \mathbf{b} - \mathbf{p}$ is perpendicular to the subspace $S$:

$$\mathbf{a}_1 \cdot \mathbf{e} = 0, \quad \mathbf{a}_2 \cdot \mathbf{e} = 0$$

$$\mathbf{a}_1^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \mathbf{a}_2^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0$$

$$\begin{pmatrix} \text{---} \; \mathbf{a}_1^T \; \text{---} \\ \text{---} \; \mathbf{a}_2^T \; \text{---} \end{pmatrix}(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$$A^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \mathbf{0}$$

$$A^T\mathbf{b} = A^T A\hat{\boldsymbol{\alpha}}$$

If $A^T A$ is invertible:

$$\hat{\boldsymbol{\alpha}} = (A^T A)^{-1} A^T \mathbf{b}$$

(When $n = 1$: $\alpha = (\mathbf{a}^T\mathbf{a})^{-1}\mathbf{a}^T\mathbf{b}$.)

$$\mathbf{p} = A\hat{\boldsymbol{\alpha}} = A(A^T A)^{-1} A^T \mathbf{b}$$

$$P = A(A^T A)^{-1} A^T$$

### 2.7 Extension to $n$-dimensional Subspace

We can easily extend this to $n$-dimensional subspace, where we have $n$ equations:

$$\mathbf{a}_1^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \mathbf{a}_2^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \ldots, \quad \mathbf{a}_n^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0$$

$$A^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \mathbf{0}$$

**Remarks:**
- $\mathbf{e} = \mathbf{b} - A\hat{\boldsymbol{\alpha}}$ is in $\mathcal{N}(A^T)$, which is perpendicular to $C(A)$.
- Left nullspace of $A$ contains the error vector.
- $\mathbf{b}$ is split into the projection $\mathbf{p}$ and the error $\mathbf{e}$.

### 2.8 Example 2: Projection onto a Subspace

$$A = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix}$$

Find $\hat{\boldsymbol{\alpha}}$, $\mathbf{p}$, and $P$.

$$A^T A = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix} = \begin{pmatrix} 3 & 3 \\ 3 & 5 \end{pmatrix}$$

$$A^T\mathbf{b} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{pmatrix}\begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} 6 \\ 0 \end{pmatrix}$$

$$\hat{\boldsymbol{\alpha}} = (A^T A)^{-1} A^T\mathbf{b} = \frac{1}{15 - 9}\begin{pmatrix} 5 & -3 \\ -3 & 3 \end{pmatrix}\begin{pmatrix} 6 \\ 0 \end{pmatrix} = \begin{pmatrix} 5 & -3 \\ -3 & 3 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 5 \\ -3 \end{pmatrix}$$

$$\mathbf{p} = A\hat{\boldsymbol{\alpha}} = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} 5 \\ -3 \end{pmatrix} = 5\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} - 3\begin{pmatrix} 0 \\ 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 5 \\ 2 \\ -1 \end{pmatrix}$$

$$P = A(A^T A)^{-1} A^T = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}\frac{1}{6}\begin{pmatrix} 5 & -3 \\ -3 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 2 \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 5 & 2 & -1 \\ 2 & 2 & 2 \\ -1 & 2 & 5 \end{pmatrix}$$

### 2.9 $A^T A$ is Invertible iff $A$ has LI Columns

**Note:** $(A^T A)^{-1} \neq A^{-1}(A^T)^{-1}$ because $A^{-1}$ does NOT exist (when $A$ is not square).

**Theorem:** $A^T A$ is invertible if and only if $A$ has linearly independent columns.

**Proof:**

Let $A \in \mathbb{R}^{m \times n}$.

($\Rightarrow$) $A\mathbf{x} = \mathbf{0}$ implies $\mathbf{x} \in \mathcal{N}(A)$. Multiplying by $A^T$ gives $A^T A\mathbf{x} = \mathbf{0}$, meaning $\mathbf{x} \in \mathcal{N}(A^T A)$. That is, $\mathcal{N}(A) \ni \mathbf{x} \longrightarrow \mathbf{x} \in \mathcal{N}(A^T A)$.

($\Leftarrow$) From $A^T A\mathbf{x} = \mathbf{0}$, multiply by $\mathbf{x}^T$:

$$\mathbf{x}^T(A^T A\mathbf{x}) = \mathbf{x}^T\mathbf{0}$$

$$(\mathbf{x}^T A^T)(A\mathbf{x}) = 0$$

$$(A\mathbf{x})^T(A\mathbf{x}) = 0 \iff \|A\mathbf{x}\|^2 = 0$$

$$\therefore A\mathbf{x} = \mathbf{0} \longrightarrow \mathbf{x} \in \mathcal{N}(A)$$

Therefore $\mathcal{N}(A) = \{\mathbf{0}\} \iff \mathcal{N}(A^T A) = \{\mathbf{0}\}$.

$\iff$ $A$ is invertible $\iff$ $A^T A$ is invertible. $\square$

**Size considerations:** $A \in \mathbb{R}^{m \times n}$, $A^T \in \mathbb{R}^{n \times m}$, $A^T A \in \mathbb{R}^{n \times n}$.

$A^T A$ is **symmetric**.

**Counter-example:** $A = \begin{pmatrix} 1 & 2 \\ 1 & 2 \\ 0 & 0 \end{pmatrix}$ has a dependent column.

$$A^T A = \begin{pmatrix} 1 & 1 & 0 \\ 2 & 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 1 & 2 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \end{pmatrix}$$

$\det(A^T A) = 2 \cdot 8 - 4 \cdot 4 = 0$. $A^T A$ is NOT invertible, $A^T A$ is singular.

### 2.10 Worked Example (4.2A)

**Problem:** Project the vector $\mathbf{b} = (3, 4, 4)$ onto the line through $\mathbf{a} = (2, 2, 1)$ and then onto the plane that also contains $\mathbf{a}^* = (1, 0, 0)$. Check that the first error vector $\mathbf{b} - \mathbf{p}$ is perpendicular to $\mathbf{a}$, and the second error vector $\mathbf{e}^* = \mathbf{b} - \mathbf{p}^*$ is also perpendicular to $\mathbf{a}$ and $\mathbf{a}^*$. Find the $3 \times 3$ projection matrix $P$ onto that plane of $\mathbf{a}$ and $\mathbf{a}^*$. Find a vector whose projection onto the plane is $\mathbf{p} = 0$. Note $P^2 = P = P^T$.

**Solution:**

**Onto a line:**

$$\mathbf{p} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\mathbf{a} = \frac{18}{9}(2, 2, 1) = (4, 4, 2) = 2\mathbf{a}$$

The error vector $\mathbf{e} = \mathbf{b} - \mathbf{p} = (-1, 0, 2)$ is perpendicular to $\mathbf{a} = (2, 2, 1)$. Check: $\mathbf{e}^T\mathbf{a} = -2 + 0 + 2 = 0$. So $\mathbf{p}$ is correct.

**Onto the plane:** The plane of $\mathbf{a} = (2, 2, 1)$ and $\mathbf{a}^* = (1, 0, 0)$ is the column space of $A = [\mathbf{a} \;\; \mathbf{a}^*]$:

$$A = \begin{pmatrix} 2 & 1 \\ 2 & 0 \\ 1 & 0 \end{pmatrix}, \quad A^T A = \begin{pmatrix} 9 & 2 \\ 2 & 1 \end{pmatrix}, \quad (A^T A)^{-1} = \frac{1}{5}\begin{pmatrix} 1 & -2 \\ -2 & 9 \end{pmatrix}$$

$$P = A(A^T A)^{-1}A^T = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0.8 & 0.4 \\ 0 & 0.4 & 0.2 \end{pmatrix}$$

Now $\mathbf{p}^* = P\mathbf{b} = (3, 4.8, 2.4)$. The error $\mathbf{e}^* = \mathbf{b} - \mathbf{p}^* = (0, -0.8, 1.6)$ is perpendicular to $\mathbf{a}$ and $\mathbf{a}^*$.

**Verification:**

$$(\mathbf{e}^*)^T\mathbf{a} = \begin{pmatrix} 0 \\ -0.8 \\ 1.6 \end{pmatrix}^T\begin{pmatrix} 2 \\ 2 \\ 1 \end{pmatrix} = -1.6 + 1.6 = 0$$

$$(\mathbf{e}^*)^T\mathbf{a}^* = \begin{pmatrix} 0 \\ -0.8 \\ 1.6 \end{pmatrix}^T\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} = 0 + 0 + 0 = 0$$

Since $A^T\mathbf{e}^* = \mathbf{0}$, we have $A(A^T A)^{-1}A^T\mathbf{e}^* = \mathbf{0}$, i.e., $P\mathbf{e}^* = \mathbf{0}$, so $\mathbf{e}^* \in \mathcal{N}(P)$.

---

<br>

## 3. Least Square Approximations (4.3)

### 3.1 Key Facts Summary

**(1)** Solving $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$ gives the projection $\mathbf{p} = A\hat{\boldsymbol{\alpha}}$ of $\mathbf{b}$ onto the column space of $A$.

**(2)** When $A\mathbf{x} = \mathbf{b}$ has no solution, $\hat{\boldsymbol{\alpha}}$ is the "least-squares solution": $\|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2 = \text{minimum}$.

**(3)** Setting derivatives of $E = \|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2$ to zero ($\frac{\partial E}{\partial \alpha_i} = 0$) also produces $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$.

**(4)** To fit points $(t_1, b_1), (t_2, b_2), \ldots, (t_m, b_m)$ by a straight line, $A$ has columns $(1, 1, \ldots, 1)$ and $(t_1, t_2, \ldots, t_m)$.

**(5)** In that case $A^T A$ is the $2 \times 2$ matrix $\begin{pmatrix} m & \sum t_i \\ \sum t_i & \sum t_i^2 \end{pmatrix}$ and $A^T\mathbf{b}$ is $\begin{pmatrix} \sum b_i \\ \sum t_i b_i \end{pmatrix}$.

### 3.2 Overdetermined Systems

It often happens that $A\mathbf{x} = \mathbf{b}$ has no solution. The matrix $A$ has more rows than columns ($m > n$). There are more equations than unknowns. The $n$ columns span a small part of $m$-dimensional space.

$\mathbf{b}$ can be outside of $C(A)$. Even in this case, we can find the approximation to $A\mathbf{x} = \mathbf{b}$ by solving $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$.

$\hat{\boldsymbol{\alpha}}$ is a **least square solution** because it minimizes:

$$E = \|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2$$

### 3.3 Example 1: Fitting a Line

Find the closest line to the points $\begin{pmatrix} 0 \\ 6 \end{pmatrix}$, $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$, and $\begin{pmatrix} 2 \\ 0 \end{pmatrix}$.

$$y = ax + b$$

$$y_1 = ax_1 + b \Rightarrow 6 = a \cdot 0 + b$$

$$y_2 = ax_2 + b \Rightarrow 0 = a \cdot 1 + b$$

$$y_3 = ax_3 + b \Rightarrow 0 = a \cdot 2 + b$$

$$\begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 1 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}$$

Unfortunately, $\mathbf{b} \notin C(A)$. $A\mathbf{x} = \mathbf{b}$ is NOT solvable.

Instead, we find the approximation $\hat{\boldsymbol{\alpha}}$ to $\hat{\mathbf{x}}$ by solving $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$. In fact, $\hat{\boldsymbol{\alpha}}$ is the minimizer of $E$.

### 3.4 Minimizing the Error: Calculus Approach

$$E = \|\mathbf{y} - \hat{\mathbf{y}}\|^2 = (y_1 - \hat{y}_1)^2 + (y_2 - \hat{y}_2)^2 + (y_3 - \hat{y}_3)^2 = \sum_{i=1}^{3}(y_i - \hat{y}_i)^2$$

where $\hat{y} = ax + b$.

$$\frac{\partial E}{\partial a} = 2\sum_{i=1}^{3}(y_i - \hat{y}_i)x_i = 0$$

$$\frac{\partial E}{\partial b} = 2\sum_{i=1}^{3}(y_i - \hat{y}_i) \cdot 1 = 0$$

From $\sum(y_i - ax_i - b)x_i = 0$: $\sum y_i x_i = a\sum x_i^2 + b\sum x_i$

From $\sum(y_i - \hat{y}_i) = 0$: $\sum y_i = a\sum x_i + b\sum 1$

$$\begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix} = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}$$

Compute:
- $\sum x_i = 0 + 1 + 2 = 3$
- $\sum x_i^2 = 0 + 1 + 4 = 5$
- $\sum y_i = 6 + 0 + 0 = 6$
- $\sum x_i y_i = 6 \cdot 0 + 0 \cdot 1 + 0 \cdot 2 = 0$

$$\begin{pmatrix} 0 \\ 6 \end{pmatrix} = \begin{pmatrix} 5 & 3 \\ 3 & 3 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}$$

$$\begin{pmatrix} a \\ b \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 3 & -3 \\ -3 & 5 \end{pmatrix}\begin{pmatrix} 0 \\ 6 \end{pmatrix} = \begin{pmatrix} -3 \\ 5 \end{pmatrix}$$

So the best line is $\hat{y} = -3x + 5$.

### 3.5 General Normal Equations for Line Fitting

$$A^T A = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}, \quad A^T\mathbf{b} = \begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix}$$

The normal equation $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$:

$$\begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix}$$

### 3.6 Minimizing the Error: Three Approaches

**How do we make the error $\mathbf{e} = \mathbf{b} - A\hat{\boldsymbol{\alpha}}$ as small as possible?**

**(1) Geometry:** We look for the closest point to $\mathbf{b}$. $\|\mathbf{e}\|$ becomes the minimum when $\mathbf{e} \perp \mathbf{a}$, i.e., $\mathbf{e} \perp C(A)$.

**(2) Algebra:** $A\mathbf{x} = \mathbf{b} = \mathbf{p} + \mathbf{e}$ is NOT solvable. $A\hat{\boldsymbol{\alpha}} = \mathbf{p}$ is solvable. $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{p}$, so $\hat{\boldsymbol{\alpha}} = (A^T A)^{-1}A^T\mathbf{p}$.

**(3) Square error for any $\mathbf{x}$:**

$$\|A\mathbf{x} - \mathbf{b}\|^2 = \|A\mathbf{x} - \mathbf{p} - \mathbf{e}\|^2$$

$$= (A\mathbf{x} - \mathbf{p} - \mathbf{e})^T(A\mathbf{x} - \mathbf{p} - \mathbf{e})$$

Since $A\mathbf{x} - \mathbf{p} \in C(A)$ and $\mathbf{e} \in \mathcal{N}(A^T)$, the cross terms vanish:

$$= \|A\mathbf{x} - \mathbf{p}\|^2 + \|\mathbf{e}\|^2$$

Choose $\mathbf{x} = \hat{\boldsymbol{\alpha}}$ so that $A\hat{\boldsymbol{\alpha}} - \mathbf{p} = \mathbf{0}$. Then the squared length $A\mathbf{x} - \mathbf{b}$ is minimized:

$$\|A\mathbf{x} - \mathbf{b}\|^2 = \|\mathbf{e}\|^2$$

The least squares solution $\hat{\boldsymbol{\alpha}}$ makes $E = \|A\mathbf{x} - \mathbf{b}\|^2$ as small as possible.

**(4) Calculus:**

$$E = \|A\mathbf{x} - \mathbf{b}\|^2 = (A\mathbf{x} - \mathbf{b})^T(A\mathbf{x} - \mathbf{b})$$

$$= (A\mathbf{x})^T A\mathbf{x} - (A\mathbf{x})^T\mathbf{b} - \mathbf{b}^T(A\mathbf{x}) + \mathbf{b}^T\mathbf{b}$$

$$= \mathbf{x}^T A^T A\mathbf{x} - 2\mathbf{x}^T A^T\mathbf{b} + \mathbf{b}^T\mathbf{b}$$

Using index notation: $E = x_i(A^T A)_{ij}x_j - 2x_i(A^T\mathbf{b})_i + b_ib_i$.

$$\frac{\partial E}{\partial x_i} = 2(A^T A)_{ij}x_j - 2(A^T\mathbf{b})_i = 0$$

This gives $A^T A\mathbf{x} = A^T\mathbf{b}$.

The partial derivatives of $E = \|A\mathbf{x} - \mathbf{b}\|^2$ are zero when $A^T A\mathbf{x} = A^T\mathbf{b}$.

### 3.7 The Big Picture for Least Squares

Split $\mathbf{b}$ into $\mathbf{p}$ and $\mathbf{e}$:

- $\hat{\boldsymbol{\alpha}} \in C(A^T)$ (row space), $A\hat{\boldsymbol{\alpha}} = \mathbf{p}$ is solvable
- $\mathbf{p} = P\mathbf{b}$ is the nearest $\mathbf{b}$ in $C(A)$
- $\mathbf{e}$ is the minimum error in $\mathcal{N}(A^T)$
- When $\mathbf{b} \notin C(A)$, $A\mathbf{x} = \mathbf{b}$ is NOT solvable

### 3.8 Example 2: Another Line Fit

Points $(-2, 1)$, $(0, 2)$, $(2, 4)$ are given. Find the straight line that minimizes the least squares error.

$$\hat{y} = ax + b$$

$$1 = a(-2) + b, \quad 2 = a(0) + b, \quad 4 = a(2) + b$$

$$\begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 1 & -2 \\ 1 & 0 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} b \\ a \end{pmatrix}$$

$$A^T A = \begin{pmatrix} 1 & 1 & 1 \\ -2 & 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & -2 \\ 1 & 0 \\ 1 & 2 \end{pmatrix} = \begin{pmatrix} 3 & 0 \\ 0 & 8 \end{pmatrix}$$

$$A^T\mathbf{b} = \begin{pmatrix} 1 & 1 & 1 \\ -2 & 0 & 2 \end{pmatrix}\begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 7 \\ 6 \end{pmatrix}$$

$$\begin{pmatrix} 3 & 0 \\ 0 & 8 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} 7 \\ 6 \end{pmatrix}$$

$$\alpha_1 = 7/3, \quad \alpha_2 = 3/4$$

Note: $\mathbf{a}_2 \perp \mathbf{a}_1$ (orthogonal column vectors), which makes $A^T A$ diagonal.

### 3.9 Dependent Columns in $A$: What is $\hat{\boldsymbol{\alpha}}$?

Which $\hat{\boldsymbol{\alpha}}$ is best if $A$ has dependent columns?

$$A\mathbf{x} = \mathbf{b}: \quad \begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 3 \\ 3 \end{pmatrix}, \quad \mathbf{b} \notin C(A)$$

$$\mathbf{p} = \begin{pmatrix} 2 \\ 2 \\ 2 \end{pmatrix} \in C(A)$$

$A\hat{\boldsymbol{\alpha}} = \mathbf{p}$ has **many solutions**: $\begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \\ 2 \end{pmatrix}$.

The problem is that $A$ has a dependent column and $\mathcal{N}(A) \ni \begin{pmatrix} 1 \\ -1 \end{pmatrix}$.

How can we find the best solution? We will learn "pseudoinverse" of $A$ in Section 4.5.

### 3.10 Fitting by a Parabola

Try to fit heights $b_1, b_2, \ldots, b_m$ at times $t_1, t_2, \ldots, t_m$ by a parabola:

$$\hat{b} = \alpha_1 + \alpha_2 t + \alpha_3 t^2$$

with $m > 3$:

$$b_1 = \alpha_1 + \alpha_2 t_1 + \alpha_3 t_1^2, \quad b_2 = \alpha_1 + \alpha_2 t_2 + \alpha_3 t_2^2, \quad \ldots, \quad b_m = \alpha_1 + \alpha_2 t_m + \alpha_3 t_m^2$$

$$\begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix} = \begin{pmatrix} 1 & t_1 & t_1^2 \\ 1 & t_2 & t_2^2 \\ \vdots & \vdots & \vdots \\ 1 & t_m & t_m^2 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \\ \alpha_3 \end{pmatrix}$$

We can find $\hat{\boldsymbol{\alpha}}$ by using least squares: $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$.

### 3.11 Example 3: Parabola Fit

$\hat{b}(t) = \alpha_1 + \alpha_2 t + \alpha_3 t^2$. Three heights $b_1, b_2, b_3$ are measured when $t = 0, 1, 2$, where $b_1 = 6$, $b_2 = 0$, $b_3 = 0$.

$$A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 1 \\ 1 & 2 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix}$$

$A\hat{\boldsymbol{\alpha}} = \mathbf{b}$. $\text{rank}(A) = 3$, so $\hat{\boldsymbol{\alpha}} = A^{-1}\mathbf{b}$.

Row reduction of $(A|\mathbf{b})$:

$$\begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 1 & 1 & 1 & | & 0 \\ 1 & 2 & 4 & | & 0 \end{pmatrix} \xrightarrow{R_2 - R_1, \, R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 2 & 4 & | & -6 \end{pmatrix}$$

$$\xrightarrow{R_3/2} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 1 & 2 & | & -3 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 0 & 1 & | & 3 \end{pmatrix}$$

$$\xrightarrow{R_2 - R_3} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 0 & | & -9 \\ 0 & 0 & 1 & | & 3 \end{pmatrix} = (I | A^{-1}\mathbf{b})$$

$$\therefore \hat{b}(t) = 6 - 9t + 3t^2$$

---

<br>

## 4. Orthogonal Bases and Gram-Schmidt (4.4)

### 4.1 Key Facts Summary

**(1)** The columns $\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n$ are **orthonormal** if:

$$\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 0 & \text{when } i \neq j \text{ (orthogonal vectors)} \\ 1 & \text{when } i = j \text{ (unit vectors: } \|\mathbf{q}_i\| = 1) \end{cases}$$

Then $Q^T Q = I_{n \times n}$.

**(2)** If $Q$ is also square, then $QQ^T = I$ and $Q^T = Q^{-1}$. Now $Q$ is an **"orthogonal matrix"**.

**(3)** The least squares solution to $Q\mathbf{x} = \mathbf{b}$ is $\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}$. Projection of $\mathbf{b}$:

$$\mathbf{p} = QQ^T\mathbf{b} = P\mathbf{b} = \mathbf{q}_1 d_1 + \mathbf{q}_2 d_2 + \cdots + \mathbf{q}_n d_n$$

**(4)** The Gram-Schmidt process takes independent $\mathbf{a}_i$ to orthogonal $\mathbf{q}_i$. Start with $\mathbf{q}_1 = \mathbf{r}_1/\|\mathbf{r}_1\|$.

**(5)** $\mathbf{q}_i$ is $(\mathbf{a}_i - \text{its projection } \mathbf{p}_i) / \|\mathbf{a}_i - \mathbf{p}_i\|$, where projection $\mathbf{p}_i = (\mathbf{a}_i^T\mathbf{q}_1)\mathbf{q}_1 + \cdots + (\mathbf{a}_i^T\mathbf{q}_{i-1})\mathbf{q}_{i-1}$.

**(6)** Each $\mathbf{a}_i$ will be a combination of $\mathbf{q}_1$ to $\mathbf{q}_n$. Then $A = QR$: orthogonal $Q$ and triangular $R$.

### 4.2 Goal: Orthogonal Columns

**i)** Orthogonal columns in $A$ are good.

$$A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)$$

$$\mathbf{a}_1^T\mathbf{a}_2 = 0, \quad \ldots, \quad \mathbf{a}_1^T\mathbf{a}_n = 0$$

$$A^T A = \begin{pmatrix} \mathbf{a}_1^T\mathbf{a}_1 & 0 & \cdots & 0 \\ 0 & \mathbf{a}_2^T\mathbf{a}_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \mathbf{a}_n^T\mathbf{a}_n \end{pmatrix}$$

This is a **diagonal matrix**.

**ii)** Construct orthogonal vectors $\mathbf{q}_i$ through Gram-Schmidt process.

### 4.3 Orthonormal Vectors and Matrices

**Definition:** The $n$ vectors are **orthonormal** if:

$$\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 0 & \text{when } i \neq j \\ 1 & \text{when } i = j \end{cases}$$

A matrix $Q$ with orthonormal columns has $Q^T Q = I$. Typically $m > n$.

Note: $Q$ is NOT required to be square. $QQ^T \neq I$ in general.

**Example:** $Q \in \mathbb{R}^{3 \times 2}$ orthonormal matrix:

$$Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix}$$

$$Q^T Q = \frac{1}{2}\begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = I_{2 \times 2}$$

$$QQ^T = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{pmatrix} \neq I_{3 \times 3}$$

When $Q$ is square, $Q^T Q = I \Longrightarrow Q^T = Q^{-1}$ (due to the definition of inverse matrix).

### 4.4 Examples of Orthogonal Matrices

**Example 1: Rotation.** $Q$ rotates every vector in the plane by the angle $\theta$.

$$Q = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}, \quad Q^T = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}$$

$$Q^T Q = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$$

$\mathbf{a}_1$ and $\mathbf{a}_2$ are an orthonormal basis for the plane $\mathbb{R}^2$. The standard basis vectors $\mathbf{i}, \mathbf{j}$ are rotated through $\theta$. $Q^{-1}$ rotates vectors back through $-\theta$.

$$Q^{-1} = \begin{pmatrix} \cos(-\theta) & -\sin(-\theta) \\ \sin(-\theta) & \cos(-\theta) \end{pmatrix}$$

**Example 2: Permutation.**

$$\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} y \\ z \\ x \end{pmatrix}, \quad \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} y \\ x \end{pmatrix}$$

All columns are unit vectors. They are also orthogonal. The inverse of a permutation matrix is its transpose: $Q^{-1} = Q^T$. **Every permutation matrix is an orthogonal matrix.**

**Example 3: Reflection.**

Given a unit normal vector $\mathbf{n}$ ($\|\mathbf{n}\| = 1$), decompose $\mathbf{u}$:

$$\mathbf{u} = \mathbf{u}_{\text{normal}} + \mathbf{u}_{\text{tangential}} = \mathbf{u}_\perp + \mathbf{u}_\parallel$$

$$\mathbf{u}_\perp = (\mathbf{u} \cdot \mathbf{n})\mathbf{n} = \mathbf{n}(\mathbf{n}^T\mathbf{u}) = (\mathbf{n}\mathbf{n}^T)\mathbf{u}$$

$$\mathbf{u}_\parallel = \mathbf{u} - \mathbf{u}_\perp = I\mathbf{u} - (\mathbf{n}\mathbf{n}^T)\mathbf{u} = (I - \mathbf{n}\mathbf{n}^T)\mathbf{u}$$

The reflection $\mathbf{v}$:

$$\mathbf{v} = \mathbf{u}_\parallel - \mathbf{u}_\perp = (I - \mathbf{n}\mathbf{n}^T)\mathbf{u} - (\mathbf{n}\mathbf{n}^T)\mathbf{u} = (I - 2\mathbf{n}\mathbf{n}^T)\mathbf{u}$$

$$Q = I - 2\mathbf{n}\mathbf{n}^T$$

**Verification that $Q$ is orthogonal:**

$$Q^T Q = (I - 2\mathbf{n}\mathbf{n}^T)^T(I - 2\mathbf{n}\mathbf{n}^T) = (I - 2\mathbf{n}\mathbf{n}^T)(I - 2\mathbf{n}\mathbf{n}^T)$$

$$= I - 4\mathbf{n}\mathbf{n}^T + 4(\mathbf{n}\mathbf{n}^T)(\mathbf{n}\mathbf{n}^T) = I - 4\mathbf{n}\mathbf{n}^T + 4\mathbf{n}\mathbf{n}^T = I$$

(since $\mathbf{n}^T\mathbf{n} = \|\mathbf{n}\|^2 = 1$).

**Key property:** Rotation, permutation, reflection matrices **preserve the length and the angle** of every vector:

$$\|Q\mathbf{x}\|^2 = (Q\mathbf{x})^T(Q\mathbf{x}) = \mathbf{x}^T Q^T Q\mathbf{x} = \mathbf{x}^T I\mathbf{x} = \|\mathbf{x}\|^2$$

If $Q$ has orthonormal columns ($Q^T Q = I$), it leaves lengths unchanged.

### 4.5 Projections Using Orthonormal Bases

$Q$ replaces $A$.

Recall the projection of $\mathbf{b}$ onto $C(A)$ is $\mathbf{p} = A(A^T A)^{-1}A^T\mathbf{b}$.

When $A = Q$ and $A^T A = Q^T Q = I$:

$$\mathbf{p} = QQ^T\mathbf{b}$$

$$\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}$$

$$\mathbf{p} = Q\hat{\boldsymbol{\alpha}} = \mathbf{q}_1 d_1 + \mathbf{q}_2 d_2 + \cdots + \mathbf{q}_n d_n$$

$$= \mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \cdots + \mathbf{q}_n(\mathbf{q}_n^T\mathbf{b})$$

**Remark:** When $Q \in \mathbb{R}^{n \times n}$ (square), $Q^T = Q^{-1}$, so $\hat{\boldsymbol{\alpha}} = Q^{-1}\mathbf{b}$ and $\mathbf{p} = QQ^{-1}\mathbf{b} = \mathbf{b}$, $P = I$.

$$\mathbf{b} = \mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \cdots + \mathbf{q}_n(\mathbf{q}_n^T\mathbf{b}) = QQ^T\mathbf{b}$$

$QQ^T = I$ is the foundation of **Fourier Series**, where they break $\mathbf{b}$ into perpendicular pieces. Then the linear combination of basis vectors puts $\mathbf{b}$ back together.

### 4.6 Example 4: Square Orthogonal Matrix

$$Q = \frac{1}{3}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix}, \quad Q^T = Q$$

$$Q^T Q = \frac{1}{9}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 9 & 0 & 0 \\ 0 & 9 & 0 \\ 0 & 0 & 9 \end{pmatrix} = I = QQ^T$$

Let $\mathbf{b} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$:

$$QQ^T\mathbf{b} = I\mathbf{b} = \mathbf{b}$$

$$\mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \mathbf{q}_3(\mathbf{q}_3^T\mathbf{b})$$

$$= \frac{1}{3}\begin{pmatrix} -1 \\ 2 \\ 2 \end{pmatrix} \cdot 2 + \frac{1}{3}\begin{pmatrix} 2 \\ -1 \\ 2 \end{pmatrix} \cdot 2 + \frac{1}{3}\begin{pmatrix} 2 \\ 2 \\ -1 \end{pmatrix} \cdot (-1) = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$$

### 4.7 The Gram-Schmidt Process

Given $\mathbf{a}, \mathbf{b}, \mathbf{c}$ vectors, create three orthogonal vectors $\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3$, which produce three orthonormal vectors:

$$\mathbf{q}_1 = \frac{\mathbf{r}_1}{\|\mathbf{r}_1\|}, \quad \mathbf{q}_2 = \frac{\mathbf{r}_2}{\|\mathbf{r}_2\|}, \quad \mathbf{q}_3 = \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|}$$

**i)** $\mathbf{r}_1 = \mathbf{a}$, so $\mathbf{q}_1 = \mathbf{r}_1/\|\mathbf{r}_1\|$.

All vectors $\mathbf{a}$, $\mathbf{r}_1$, $\mathbf{q}_1$ are on a line.

**ii)** $\mathbf{r}_2$ must be perpendicular to $\mathbf{r}_1$.

Split $\mathbf{b}$ into $\mathbf{b}_\parallel$ and $\mathbf{b}_\perp$ to $\mathbf{q}_1$:

$$\mathbf{b}_\parallel = (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1 = \left(\mathbf{b} \cdot \frac{\mathbf{r}_1}{\|\mathbf{r}_1\|}\right)\frac{\mathbf{r}_1}{\|\mathbf{r}_1\|} = \frac{1}{\|\mathbf{r}_1\|^2}(\mathbf{r}_1^T\mathbf{b})\,\mathbf{r}_1 = \left(\frac{\mathbf{r}_1^T\mathbf{b}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1$$

$$\mathbf{b}_\perp = \mathbf{b} - \mathbf{b}_\parallel = \mathbf{b} - \left(\frac{\mathbf{r}_1^T\mathbf{b}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1 = \mathbf{r}_2$$

$$\mathbf{q}_2 = \frac{\mathbf{r}_2}{\|\mathbf{r}_2\|}$$

**iii)** $\mathbf{r}_3$ must be perpendicular to $\mathbf{r}_1, \mathbf{r}_2$.

Subtract from every new vector its projection in the direction already set:

$$\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2$$

$$= \mathbf{c} - \left(\frac{\mathbf{r}_1^T\mathbf{c}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1 - \left(\frac{\mathbf{r}_2^T\mathbf{c}}{\mathbf{r}_2^T\mathbf{r}_2}\right)\mathbf{r}_2$$

$$\mathbf{q}_3 = \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|}$$

All vectors $\mathbf{a}, \mathbf{b}, \mathbf{c}, \mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3, \mathbf{q}_1, \mathbf{q}_2, \mathbf{q}_3$ are in one subspace ($\mathbb{R}^3$).

### 4.8 Gram-Schmidt Example

$$\mathbf{a} = \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix}, \quad \mathbf{c} = \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix}$$

**i)** $\mathbf{r}_1 = \mathbf{a} = \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}$, $\|\mathbf{r}_1\| = \sqrt{2}$:

$$\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}$$

**ii)** $\mathbf{r}_2 = \mathbf{b} - (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1$:

$$= \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix} - \frac{1}{\sqrt{2}} \cdot 2 \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} = \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix} - \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}$$

$$\mathbf{q}_2 = \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}$$

**iii)** $\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2$:

$$= \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix} - \frac{6}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} - \frac{-6}{\sqrt{6}} \cdot \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}$$

$$= \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix} - 3\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + \begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$$

$$\mathbf{q}_3 = \frac{1}{\sqrt{3}}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$$

### 4.9 Factorization $A = QR$

$$(\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = (\mathbf{q}_1 \;\; \mathbf{q}_2 \;\; \mathbf{q}_3)\,R$$

The vectors $\mathbf{a}, \mathbf{b}, \mathbf{c}$ are combinations of $\mathbf{q}_1, \mathbf{q}_2, \mathbf{q}_3$:

$$\mathbf{a} = (\mathbf{q}_1 \cdot \mathbf{a})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{a})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{a})\mathbf{q}_3$$

$$\mathbf{b} = (\mathbf{q}_1 \cdot \mathbf{b})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{b})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{b})\mathbf{q}_3$$

$$\mathbf{c} = (\mathbf{q}_1 \cdot \mathbf{c})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{c})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{c})\mathbf{q}_3$$

Note: Later $\mathbf{q}$'s are orthogonal to earlier $\mathbf{a}$'s, so many entries in $R$ are zero (below diagonal):

$$(\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = (\mathbf{q}_1 \;\; \mathbf{q}_2 \;\; \mathbf{q}_3)\begin{pmatrix} \mathbf{q}_1^T\mathbf{a} & \mathbf{q}_1^T\mathbf{b} & \mathbf{q}_1^T\mathbf{c} \\ 0 & \mathbf{q}_2^T\mathbf{b} & \mathbf{q}_2^T\mathbf{c} \\ 0 & 0 & \mathbf{q}_3^T\mathbf{c} \end{pmatrix}$$

$$A = QR$$

$Q^T A = Q^T Q R = IR = R$.

From LI $\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n$, Gram-Schmidt conducts orthonormal vectors $\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n$. The matrices with them satisfy $A = QR$. Then $R = Q^T A$ is **upper triangular** because later $\mathbf{q}$'s are orthogonal to earlier $\mathbf{a}$'s.

**From the example:**

$$A = (\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = \begin{pmatrix} 1 & 2 & 3 \\ -1 & 0 & -3 \\ 0 & -2 & 3 \end{pmatrix}$$

$$= \frac{1}{\sqrt{6}}\begin{pmatrix} \sqrt{3} & 1 & \sqrt{2} \\ -\sqrt{3} & 1 & \sqrt{2} \\ 0 & -2 & \sqrt{2} \end{pmatrix} \begin{pmatrix} \sqrt{2} & \sqrt{2} & \sqrt{18} \\ 0 & \sqrt{6} & \sqrt{6} \\ 0 & 0 & \sqrt{3} \end{pmatrix}$$

$$Q \qquad \qquad R$$

The **lengths** of $\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3$ are the diagonal of $R$, which is positive.

Any $m$ by $n$ matrix $A$ with independent columns can be factored into $A = QR$. $Q \in \mathbb{R}^{m \times n}$ has orthogonal columns, $R \in \mathbb{R}^{n \times n}$ is upper triangular with positive diagonal.

### 4.10 Application: Least Squares via QR

$$A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$$

Substituting $A = QR$: $A^T A = (QR)^T QR = R^T Q^T QR = R^T R$.

$$R^T R\hat{\boldsymbol{\alpha}} = R^T Q^T\mathbf{b}$$

Since $R$ is invertible:

$$R\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b} \quad \text{(very fast)}$$

$$\hat{\boldsymbol{\alpha}} = R^{-1}Q^T\mathbf{b}$$

---

<br>

## 5. The Pseudoinverse of a Matrix (4.5)

### 5.1 Two-Sided and One-Sided Inverses

**(1)** **Two-sided inverse:**

$$A^{-1}A = AA^{-1} = I$$

**(2)** **One-sided inverse:**

$$A^+ A = I \quad (A^+ \text{ is the left inverse of } A)$$

$$AA^+ = I \quad (A^+ \text{ is the right inverse of } A)$$

Every matrix $A$ has a pseudoinverse $A^+$.

### 5.2 Three Cases with Identity-like Matrices

**(1)** $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}_{2 \times 2}$, $m = n = 2$, $\text{rank}(I) = 2 = r$.

$\dim C(I) = 2$, $\dim \mathcal{N}(I^T) = 0$, $\dim C(I^T) = 2$, $\dim \mathcal{N}(I) = 0$.

All four subspaces are either full or trivial.

**(2)** $I_L = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}_{2 \times 3}$, $m = 2$, $n = 3$, $\text{rank}(I_L) = 2 = r = m$.

Every row is LI, but $r = m < n$. Nullspace has nontrivial elements.

$\dim C(I_L) = 2$, $\dim \mathcal{N}(I_L^T) = 0$, $\dim C(I_L^T) = 2$, $\dim \mathcal{N}(I_L) = 1$.

$I_L I_R = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I_{2 \times 2}$

$I_L$ = left inverse of $I_R$, $I_R$ = right inverse of $I_L$.

**(3)** $I_R = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}_{3 \times 2}$, $m = 3$, $n = 2$, $\text{rank}(I_R) = 2 = r = n$.

Every column is LI, but $r < m$.

$\dim C(I_R) = 2$, $\dim \mathcal{N}(I_R^T) = 1$, $\dim C(I_R^T) = 2$, $\dim \mathcal{N}(I_R) = 0$.

### 5.3 Left Inverse and Right Inverse

| | $A^+A = I_{n \times n}$ (Left Inverse) | $AA^+ = I_{m \times m}$ (Right Inverse) |
|:---|:---|:---|
| **Condition** | Full column rank: $r = n < m$ | Full row rank: $r = m < n$ |
| **Interpretation** | # unknowns < # equations | # equations < # unknowns |
| **Solution to $A\mathbf{x} = \mathbf{b}$** | 0 or 1 solution | Infinitely many solutions |
| **Nullspace** | $\mathcal{N}(A) = \{\mathbf{0}\}$ | $\mathcal{N}(A^T) = \{\mathbf{0}\}$ |
| **$A^T A$ or $AA^T$** | $A^T A$ is $n \times n$ and invertible | $AA^T$ is $m \times m$ and invertible |
| **Formula** | $A^+ = (A^T A)^{-1}A^T$ | $A^+ = A^T(AA^T)^{-1}$ |

**Left inverse case:** $A^+A = I$ describes the matrices in least squares. $\hat{\boldsymbol{\alpha}} = A^+\mathbf{b}$ is the solution to $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$, where $\mathbf{b}$ may not be in $C(A)$ and $\hat{\boldsymbol{\alpha}}$ is in $C(A^T)$.

**Right inverse case:** $AA^+ = I$. $\mathbf{x}^+ = A^+\mathbf{b}$ is the least length solution to $A\mathbf{x} = \mathbf{b}$ where $\mathbf{b}$ is in $C(A)$ and $\mathbf{x}^+$ is in $C(A^T)$.

### 5.4 The Pseudoinverse $A^+$ of a General Matrix $A_{m \times n}$

**Step 1:** Every vector $\mathbf{b} \in \mathbb{R}^m$ has two perpendicular parts: $\mathbf{p}$ and $\mathbf{z}$, such that $\mathbf{b} = \mathbf{p} + \mathbf{z}$.

**Step 2:** $\mathbf{p} \in C(A)$ and $\mathbf{z} \in \mathcal{N}(A^T)$ ($\Leftrightarrow A^T\mathbf{z} = \mathbf{0}$, $A^+\mathbf{z} = \mathbf{0}$).

**Step 3:** $C(A^T) \ni \mathbf{x}^+ \longrightarrow A\mathbf{x}^+ = \mathbf{p}$. Invert this part: $\mathbf{x}^+ = A^+\mathbf{p}$.

**Step 4:** $A^+\mathbf{b} = A^+(\mathbf{p} + \mathbf{z}) = A^+\mathbf{p} + A^+\mathbf{z} = \mathbf{x}^+$.

### 5.5 Big Picture for $A^+$

- $A^+$ shares the same four subspaces as $A^T$.
- $A^+$ inverts $A$ when it can: from column space to row space.

$$C(A^T) \xrightarrow{A} C(A) \xrightarrow{A^+} C(A^T)$$

- $A^+A$ brings $\mathbf{x}^+ \in C(A^T)$ back to the same $\mathbf{x}^+$.
- $A^+A$ is an $n \times n$ projection matrix $P_{\text{row}}$ onto the row space of $A$, $C(A^T)$.

$$P_{\text{row}} = A^+A = (A^+A)^2 = (A^+A)^T$$

- $AA^+$ is an $m \times m$ projection matrix $P_{\text{col}}$ onto the column space of $A$, $C(A)$.

$$P_{\text{col}} = AA^+ = (AA^+)^2 = (AA^+)^T$$

### 5.6 Example 1: Finding $A^+$

$$A = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ b_2 \end{pmatrix} = \mathbf{p} + \mathbf{z}$$

$r = 1$.

**i)** $\mathbf{p} \in C(A)$: $A\mathbf{x} = \mathbf{p}$: $\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ 0 \end{pmatrix}$

$2x_1 = b_1$, $x_1 = b_1/2$.

$$\mathbf{x} = \begin{pmatrix} b_1/2 \\ x_2 \end{pmatrix} = \underbrace{\begin{pmatrix} b_1/2 \\ 0 \end{pmatrix}}_{\mathbf{x}^+} + \underbrace{\begin{pmatrix} 0 \\ x_2 \end{pmatrix}}_{\mathbf{x}_n}$$

**ii)** $\mathbf{z} \in \mathcal{N}(A^T)$: $A^T\mathbf{z} = \mathbf{0}$: $\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} z_1 \\ z_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, so $z_1 = 0$, $\mathbf{z} = \begin{pmatrix} 0 \\ z_2 \end{pmatrix}$.

**iii)** $A^+\mathbf{b} = A^+(\mathbf{p} + \mathbf{z}) = A^+\mathbf{p} + A^+\mathbf{z} = \mathbf{x}^+ = \begin{pmatrix} b_1/2 \\ 0 \end{pmatrix}$

$$= \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$$

$$A^+ = \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}$$

**iv)** Verification:

$$AA^+ = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = P_{\text{col}}$$

$$A^+A = \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = P_{\text{row}}$$

### 5.7 Pseudoinverse of a Diagonal Matrix

$$D = \begin{pmatrix} 2 & & \\ & 3 & \\ & & 0 \end{pmatrix}, \quad D^+ = \begin{pmatrix} 1/2 & & \\ & 1/3 & \\ & & 0 \end{pmatrix}$$

Invert the nonzero diagonal entries; leave zeros as zeros.

### 5.8 Pseudoinverse via Factorization

Let $U, V$ be invertible matrices. Let $A = UDV^T$. Then:

$$A^+ = (UDV^T)^+ = V^{-T}D^+U^{-1}$$

**Example (full column rank):** Let $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $\text{rank}(A) = 1 = r = n < m$. Full column rank, $\exists$ left inverse.

$A^+A = 1$. $A^+ = (A^T A)^{-1}A^T = ((1\;1)\begin{pmatrix} 1 \\ 1 \end{pmatrix})^{-1}(1\;1) = \frac{1}{2}(1\;1)$.

**Example (full row rank):** Let $A = (1\;1)$, $\text{rank}(A) = 1 = m < n$. Full row rank, $\exists$ right inverse.

$AA^+ = 1$. $A^+ = A^T(AA^T)^{-1} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}((1\;1)\begin{pmatrix} 1 \\ 1 \end{pmatrix})^{-1} = \frac{1}{2}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$.

### 5.9 Important Action of $A$: Row Space to Column Space

$$C(A^T) \xrightarrow{A} C(A)$$

$$C(A) \xrightarrow{A^+} C(A^T)$$

**i)** $\mathbf{x}_1, \mathbf{x}_2 \in C(A^T)$, $\mathbf{x}_1 \neq \mathbf{x}_2 \Rightarrow A\mathbf{x}_1 \neq A\mathbf{x}_2$.

**Proof:** Suppose $A\mathbf{x}_1 = \mathbf{b}$ and $A\mathbf{x}_2 = \mathbf{b}$. Subtract: $A(\mathbf{x}_1 - \mathbf{x}_2) = \mathbf{0}$. By definition of nullspace, $\mathbf{x}_1 - \mathbf{x}_2 \in \mathcal{N}(A)$. But $\mathbf{x}_1 - \mathbf{x}_2 \in C(A^T)$. Since $\mathcal{N}(A)$ is orthogonal to $C(A^T)$, $\mathbf{x}_1 - \mathbf{x}_2 = \mathbf{0}$. This contradicts $\mathbf{x}_1 \neq \mathbf{x}_2$. Therefore $A\mathbf{x}_1 \neq A\mathbf{x}_2$. $\square$

**ii)** $\forall\,\mathbf{b} \in C(A)$, $\exists!\,\mathbf{x}^+ \in C(A^T)$ such that $A^+\mathbf{b} = \mathbf{x}^+$.

### 5.10 Pseudoinverse $A^+ = R^+C^+$ of $A = CR$

Given $A = CR$ where $C$ has $r$ independent columns (full column rank) and $R$ has $r$ independent rows (full row rank):

- $C$ full column rank: $\exists (C^T C)^{-1}$, $\exists$ left inverse of $C$: $C^+C = I$, $C^+ = (C^T C)^{-1}C^T$.
- $R$ full row rank: $\exists (RR^T)^{-1}$, $\exists$ right inverse of $R$: $RR^+ = I$, $R^+ = R^T(RR^T)^{-1}$.

$$A^+ = R^+C^+ = R^T(RR^T)^{-1}(C^T C)^{-1}C^T$$

$$= R^T(C^T A R^T)^{-1}C^T$$

(since $C^T C R R^T = C^T(CR)R^T = C^T A R^T$).

### 5.11 Example: Pseudoinverse of a $5 \times 4$ Matrix

$$A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ -1 & 0 & 1 & 0 \\ 0 & -1 & 1 & 0 \\ -1 & 0 & 0 & 1 \\ 0 & -1 & 0 & 1 \end{pmatrix}_{5 \times 4}$$

Row reduction steps: $R_2 - R_1$, $R_4 - R_1$, then $R_3 - R_2$, $R_4 - R_2$, $R_5 - R_2$, swap $R_3$ and $R_5$, $R_4 - R_3$, then $R_1 + R_2$, $R_2 + R_3$, $\times(-1)$:

$$R_0 = \begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}$$

Columns 1, 2, 3 are LI. $\text{rank}(A) = 3$.

$$R = \begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{pmatrix}$$

$$A = CR = \begin{pmatrix} -1 & 1 & 0 \\ -1 & 0 & 1 \\ 0 & -1 & 1 \\ -1 & 0 & 0 \\ 0 & -1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{pmatrix}$$

$A^+ = (CR)^+ = R^+C^+ = R^T(C^T AR^T)^{-1}C^T$.

It turns out $C^T AR^T = \begin{pmatrix} 4 & 0 & 0 \\ 0 & 4 & 0 \\ -1 & -1 & 2 \end{pmatrix}$.

$$A^+ = \frac{1}{8}\begin{pmatrix} -2 & -2 & 0 & -2 & 0 \\ 2 & 0 & -2 & 0 & -2 \\ 0 & 3 & 3 & -1 & -1 \\ 0 & -1 & -1 & 3 & 3 \end{pmatrix}$$

The mapping through $A = CR$:

$$C(A^T) \ni \mathbf{x} \xrightarrow{R} R\mathbf{x} \xrightarrow{C} CR\mathbf{x} = A\mathbf{x} \in C(A)$$

$R$ has full row rank: $\mathcal{N}(R^T) = \{\mathbf{0}\}$. $C$ has full column rank: $\mathcal{N}(C) = \{\mathbf{0}\}$.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Orthogonal vectors | $\mathbf{v}^T\mathbf{w} = 0$ implies $\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2$ |
| Fundamental subspace orthogonality | $C(A^T) \perp \mathcal{N}(A)$ in $\mathbb{R}^n$; $C(A) \perp \mathcal{N}(A^T)$ in $\mathbb{R}^m$ |
| Orthogonal complement | $\dim V + \dim V^\perp = n$; every $\mathbf{x} = \mathbf{x}_r + \mathbf{x}_n$ |
| Projection onto a line | $\mathbf{p} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\mathbf{a}$; projection matrix $P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}$ |
| Projection onto a subspace | $\mathbf{p} = A(A^T A)^{-1}A^T\mathbf{b}$; $P = A(A^T A)^{-1}A^T$; $P^2 = P = P^T$ |
| $A^T A$ invertibility | $A^T A$ is invertible $\iff$ $A$ has LI columns $\iff$ $\mathcal{N}(A) = \{\mathbf{0}\}$ |
| Least squares solution | $\hat{\boldsymbol{\alpha}} = (A^T A)^{-1}A^T\mathbf{b}$ minimizes $E = \|A\mathbf{x} - \mathbf{b}\|^2$ |
| Normal equations | $A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}$ (geometry, algebra, or calculus derivation) |
| Line fitting | $A = \begin{pmatrix} x_1 & 1 \\ \vdots & \vdots \\ x_m & 1 \end{pmatrix}$; $A^T A = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & m \end{pmatrix}$ |
| Parabola fitting | $A$ has columns $(1, 1, \ldots)$, $(t_1, t_2, \ldots)$, $(t_1^2, t_2^2, \ldots)$ |
| Orthonormal vectors | $\mathbf{q}_i^T\mathbf{q}_j = \delta_{ij}$; $Q^T Q = I$ |
| Orthogonal matrix (square $Q$) | $Q^T Q = QQ^T = I$; $Q^T = Q^{-1}$; preserves lengths and angles |
| Examples of orthogonal matrices | Rotation, permutation, reflection ($Q = I - 2\mathbf{n}\mathbf{n}^T$) |
| Projection with orthonormal basis | $\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}$; $\mathbf{p} = QQ^T\mathbf{b}$; if $Q$ square, $P = I$ |
| Gram-Schmidt process | $\mathbf{r}_1 = \mathbf{a}$; $\mathbf{r}_2 = \mathbf{b} - (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1$; $\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2$; normalize each |
| QR Factorization | $A = QR$; $Q$ orthogonal columns, $R$ upper triangular with positive diagonal |
| Least squares via QR | $R\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}$ (much faster than normal equations) |
| Pseudoinverse $A^+$ | Inverts $A$ from $C(A) \to C(A^T)$; $A^+\mathbf{z} = \mathbf{0}$ for $\mathbf{z} \in \mathcal{N}(A^T)$ |
| Left inverse ($r = n < m$) | $A^+ = (A^T A)^{-1}A^T$; $A^+A = I$; least squares solution |
| Right inverse ($r = m < n$) | $A^+ = A^T(AA^T)^{-1}$; $AA^+ = I$; minimum length solution |
| $A^+A$ and $AA^+$ | $A^+A = P_{\text{row}}$ (projection onto row space); $AA^+ = P_{\text{col}}$ (projection onto column space) |
| Pseudoinverse via $A = CR$ | $A^+ = R^+C^+ = R^T(C^T AR^T)^{-1}C^T$ |

---
