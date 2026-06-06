# Chapter 4 Lecture — Orthogonality — Part 2/2

> **📑 This document is split into 2 parts.**
>
> [Part 1](Concepts.md) · **Part 2**

---

<br>

## Table of Contents

- [3. Least Square Approximations (4.3)](#3-least-square-approximations-43)
- [4. Orthogonal Bases and Gram-Schmidt (4.4)](#4-orthogonal-bases-and-gram-schmidt-44)
- [5. The Pseudoinverse of a Matrix (4.5)](#5-the-pseudoinverse-of-a-matrix-45)
- [Summary](#summary)

---

<br>

## 3. Least Square Approximations (4.3)

### 3.1 Key Facts Summary

**(1)** Solving $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$ gives the projection $`\mathbf{p} = A\hat{\boldsymbol{\alpha}}`$ of $`\mathbf{b}`$ onto the column space of $`A`$.

**(2)** When $`A\mathbf{x} = \mathbf{b}`$ has no solution, $`\hat{\boldsymbol{\alpha}}`$ is the "least-squares solution": $`\|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2 = \text{minimum}`$.

**(3)** Setting derivatives of $`E = \|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2`$ to zero ($`\frac{\partial E}{\partial \alpha_i} = 0`$) also produces $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$.

**(4)** To fit points $`(t_1, b_1), (t_2, b_2), \ldots, (t_m, b_m)`$ by a straight line, $`A`$ has columns $`(1, 1, \ldots, 1)`$ and $`(t_1, t_2, \ldots, t_m)`$.

**(5)** In that case $`A^T A`$ is the $`2 \times 2`$ matrix $`\begin{pmatrix} m & \sum t_i \\ \sum t_i & \sum t_i^2 \end{pmatrix}`$ and $`A^T\mathbf{b}`$ is $`\begin{pmatrix} \sum b_i \\ \sum t_i b_i \end{pmatrix}`$.

### 3.2 Overdetermined Systems

It often happens that $`A\mathbf{x} = \mathbf{b}`$ has no solution. The matrix $`A`$ has more rows than columns ($`m > n`$). There are more equations than unknowns. The $`n`$ columns span a small part of $`m`$-dimensional space.

$`\mathbf{b}`$ can be outside of $`C(A)`$. Even in this case, we can find the approximation to $`A\mathbf{x} = \mathbf{b}`$ by solving $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$.

$`\hat{\boldsymbol{\alpha}}`$ is a **least square solution** because it minimizes:

```math
E = \|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2
```

### 3.3 Example 1: Fitting a Line

Find the closest line to the points $`\begin{pmatrix} 0 \\ 6 \end{pmatrix}`$, $`\begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, and $`\begin{pmatrix} 2 \\ 0 \end{pmatrix}`$.

```math
y = ax + b
```

```math
y_1 = ax_1 + b \Rightarrow 6 = a \cdot 0 + b
```

```math
y_2 = ax_2 + b \Rightarrow 0 = a \cdot 1 + b
```

```math
y_3 = ax_3 + b \Rightarrow 0 = a \cdot 2 + b
```

```math
\begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 1 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}
```

Unfortunately, $`\mathbf{b} \notin C(A)`$. $`A\mathbf{x} = \mathbf{b}`$ is NOT solvable.

Instead, we find the approximation $`\hat{\boldsymbol{\alpha}}`$ to $`\hat{\mathbf{x}}`$ by solving $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$. In fact, $`\hat{\boldsymbol{\alpha}}`$ is the minimizer of $`E`$.

### 3.4 Minimizing the Error: Calculus Approach

```math
E = \|\mathbf{y} - \hat{\mathbf{y}}\|^2 = (y_1 - \hat{y}_1)^2 + (y_2 - \hat{y}_2)^2 + (y_3 - \hat{y}_3)^2 = \sum_{i=1}^{3}(y_i - \hat{y}_i)^2
```

where $`\hat{y} = ax + b`$.

```math
\frac{\partial E}{\partial a} = 2\sum_{i=1}^{3}(y_i - \hat{y}_i)x_i = 0
```

```math
\frac{\partial E}{\partial b} = 2\sum_{i=1}^{3}(y_i - \hat{y}_i) \cdot 1 = 0
```

From $`\sum(y_i - ax_i - b)x_i = 0`$: $`\sum y_i x_i = a\sum x_i^2 + b\sum x_i`$

From $`\sum(y_i - \hat{y}_i) = 0`$: $`\sum y_i = a\sum x_i + b\sum 1`$

```math
\begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix} = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}
```

Compute:
- $`\sum x_i = 0 + 1 + 2 = 3`$
- $`\sum x_i^2 = 0 + 1 + 4 = 5`$
- $`\sum y_i = 6 + 0 + 0 = 6`$
- $`\sum x_i y_i = 6 \cdot 0 + 0 \cdot 1 + 0 \cdot 2 = 0`$

```math
\begin{pmatrix} 0 \\ 6 \end{pmatrix} = \begin{pmatrix} 5 & 3 \\ 3 & 3 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}
```

```math
\begin{pmatrix} a \\ b \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 3 & -3 \\ -3 & 5 \end{pmatrix}\begin{pmatrix} 0 \\ 6 \end{pmatrix} = \begin{pmatrix} -3 \\ 5 \end{pmatrix}
```

So the best line is $`\hat{y} = -3x + 5`$.

### 3.5 General Normal Equations for Line Fitting

```math
A^T A = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}, \quad A^T\mathbf{b} = \begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix}
```

The normal equation $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$:

```math
\begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix}
```

### 3.6 Minimizing the Error: Three Approaches

**How do we make the error $`\mathbf{e} = \mathbf{b} - A\hat{\boldsymbol{\alpha}}`$ as small as possible?**

**(1) Geometry:** We look for the closest point to $`\mathbf{b}`$. $`\|\mathbf{e}\|`$ becomes the minimum when $`\mathbf{e} \perp \mathbf{a}`$, i.e., $`\mathbf{e} \perp C(A)`$.

**(2) Algebra:** $`A\mathbf{x} = \mathbf{b} = \mathbf{p} + \mathbf{e}`$ is NOT solvable. $`A\hat{\boldsymbol{\alpha}} = \mathbf{p}`$ is solvable. $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{p}`$, so $`\hat{\boldsymbol{\alpha}} = (A^T A)^{-1}A^T\mathbf{p}`$.

**(3) Square error for any $`\mathbf{x}`$:**

```math
\|A\mathbf{x} - \mathbf{b}\|^2 = \|A\mathbf{x} - \mathbf{p} - \mathbf{e}\|^2
```

```math
= (A\mathbf{x} - \mathbf{p} - \mathbf{e})^T(A\mathbf{x} - \mathbf{p} - \mathbf{e})
```

Since $`A\mathbf{x} - \mathbf{p} \in C(A)`$ and $`\mathbf{e} \in \mathcal{N}(A^T)`$, the cross terms vanish:

```math
= \|A\mathbf{x} - \mathbf{p}\|^2 + \|\mathbf{e}\|^2
```

Choose $`\mathbf{x} = \hat{\boldsymbol{\alpha}}`$ so that $`A\hat{\boldsymbol{\alpha}} - \mathbf{p} = \mathbf{0}`$. Then the squared length $`A\mathbf{x} - \mathbf{b}`$ is minimized:

```math
\|A\mathbf{x} - \mathbf{b}\|^2 = \|\mathbf{e}\|^2
```

The least squares solution $`\hat{\boldsymbol{\alpha}}`$ makes $`E = \|A\mathbf{x} - \mathbf{b}\|^2`$ as small as possible.

**(4) Calculus:**

```math
E = \|A\mathbf{x} - \mathbf{b}\|^2 = (A\mathbf{x} - \mathbf{b})^T(A\mathbf{x} - \mathbf{b})
```

```math
= (A\mathbf{x})^T A\mathbf{x} - (A\mathbf{x})^T\mathbf{b} - \mathbf{b}^T(A\mathbf{x}) + \mathbf{b}^T\mathbf{b}
```

```math
= \mathbf{x}^T A^T A\mathbf{x} - 2\mathbf{x}^T A^T\mathbf{b} + \mathbf{b}^T\mathbf{b}
```

Using index notation: $`E = x_i(A^T A)_{ij}x_j - 2x_i(A^T\mathbf{b})_i + b_ib_i`$.

```math
\frac{\partial E}{\partial x_i} = 2(A^T A)_{ij}x_j - 2(A^T\mathbf{b})_i = 0
```

This gives $`A^T A\mathbf{x} = A^T\mathbf{b}`$.

The partial derivatives of $`E = \|A\mathbf{x} - \mathbf{b}\|^2`$ are zero when $`A^T A\mathbf{x} = A^T\mathbf{b}`$.

### 3.7 The Big Picture for Least Squares

Split $`\mathbf{b}`$ into $`\mathbf{p}`$ and $`\mathbf{e}`$:

- $`\hat{\boldsymbol{\alpha}} \in C(A^T)`$ (row space), $`A\hat{\boldsymbol{\alpha}} = \mathbf{p}`$ is solvable
- $`\mathbf{p} = P\mathbf{b}`$ is the nearest $`\mathbf{b}`$ in $`C(A)`$
- $`\mathbf{e}`$ is the minimum error in $`\mathcal{N}(A^T)`$
- When $`\mathbf{b} \notin C(A)`$, $`A\mathbf{x} = \mathbf{b}`$ is NOT solvable

### 3.8 Example 2: Another Line Fit

Points $`(-2, 1)`$, $`(0, 2)`$, $`(2, 4)`$ are given. Find the straight line that minimizes the least squares error.

```math
\hat{y} = ax + b
```

```math
1 = a(-2) + b, \quad 2 = a(0) + b, \quad 4 = a(2) + b
```

```math
\begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 1 & -2 \\ 1 & 0 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} b \\ a \end{pmatrix}
```

```math
A^T A = \begin{pmatrix} 1 & 1 & 1 \\ -2 & 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & -2 \\ 1 & 0 \\ 1 & 2 \end{pmatrix} = \begin{pmatrix} 3 & 0 \\ 0 & 8 \end{pmatrix}
```

```math
A^T\mathbf{b} = \begin{pmatrix} 1 & 1 & 1 \\ -2 & 0 & 2 \end{pmatrix}\begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 7 \\ 6 \end{pmatrix}
```

```math
\begin{pmatrix} 3 & 0 \\ 0 & 8 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} 7 \\ 6 \end{pmatrix}
```

```math
\alpha_1 = 7/3, \quad \alpha_2 = 3/4
```

Note: $`\mathbf{a}_2 \perp \mathbf{a}_1`$ (orthogonal column vectors), which makes $`A^T A`$ diagonal.

### 3.9 Dependent Columns in $`A`$: What is $`\hat{\boldsymbol{\alpha}}`$?

Which $`\hat{\boldsymbol{\alpha}}`$ is best if $`A`$ has dependent columns?

```math
A\mathbf{x} = \mathbf{b}: \quad \begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 3 \\ 3 \end{pmatrix}, \quad \mathbf{b} \notin C(A)
```

```math
\mathbf{p} = \begin{pmatrix} 2 \\ 2 \\ 2 \end{pmatrix} \in C(A)
```

$`A\hat{\boldsymbol{\alpha}} = \mathbf{p}`$ has **many solutions**: $`\begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \\ 2 \end{pmatrix}`$.

The problem is that $`A`$ has a dependent column and $`\mathcal{N}(A) \ni \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$.

How can we find the best solution? We will learn "pseudoinverse" of $`A`$ in Section 4.5.

### 3.10 Fitting by a Parabola

Try to fit heights $`b_1, b_2, \ldots, b_m`$ at times $`t_1, t_2, \ldots, t_m`$ by a parabola:

```math
\hat{b} = \alpha_1 + \alpha_2 t + \alpha_3 t^2
```

with $`m > 3`$:

```math
b_1 = \alpha_1 + \alpha_2 t_1 + \alpha_3 t_1^2, \quad b_2 = \alpha_1 + \alpha_2 t_2 + \alpha_3 t_2^2, \quad \ldots, \quad b_m = \alpha_1 + \alpha_2 t_m + \alpha_3 t_m^2
```

```math
\begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix} = \begin{pmatrix} 1 & t_1 & t_1^2 \\ 1 & t_2 & t_2^2 \\ \vdots & \vdots & \vdots \\ 1 & t_m & t_m^2 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \\ \alpha_3 \end{pmatrix}
```

We can find $`\hat{\boldsymbol{\alpha}}`$ by using least squares: $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$.

### 3.11 Example 3: Parabola Fit

$`\hat{b}(t) = \alpha_1 + \alpha_2 t + \alpha_3 t^2`$. Three heights $`b_1, b_2, b_3`$ are measured when $`t = 0, 1, 2`$, where $`b_1 = 6`$, $`b_2 = 0`$, $`b_3 = 0`$.

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 1 \\ 1 & 2 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix}
```

$`A\hat{\boldsymbol{\alpha}} = \mathbf{b}`$. $`\text{rank}(A) = 3`$, so $`\hat{\boldsymbol{\alpha}} = A^{-1}\mathbf{b}`$.

Row reduction of $`(A|\mathbf{b})`$:

```math
\begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 1 & 1 & 1 & | & 0 \\ 1 & 2 & 4 & | & 0 \end{pmatrix} \xrightarrow{R_2 - R_1, \, R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 2 & 4 & | & -6 \end{pmatrix}
```

```math
\xrightarrow{R_3/2} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 1 & 2 & | & -3 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 0 & 1 & | & 3 \end{pmatrix}
```

```math
\xrightarrow{R_2 - R_3} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 0 & | & -9 \\ 0 & 0 & 1 & | & 3 \end{pmatrix} = (I | A^{-1}\mathbf{b})
```

```math
\therefore \hat{b}(t) = 6 - 9t + 3t^2
```

---

<br>

## 4. Orthogonal Bases and Gram-Schmidt (4.4)

### 4.1 Key Facts Summary

**(1)** The columns $`\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n`$ are **orthonormal** if:

```math
\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 0 & \text{when } i \neq j \text{ (orthogonal vectors)} \\ 1 & \text{when } i = j \text{ (unit vectors: } \|\mathbf{q}_i\| = 1) \end{cases}
```

Then $`Q^T Q = I_{n \times n}`$.

**(2)** If $`Q`$ is also square, then $`QQ^T = I`$ and $`Q^T = Q^{-1}`$. Now $`Q`$ is an **"orthogonal matrix"**.

**(3)** The least squares solution to $`Q\mathbf{x} = \mathbf{b}`$ is $`\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}`$. Projection of $`\mathbf{b}`$:

```math
\mathbf{p} = QQ^T\mathbf{b} = P\mathbf{b} = \mathbf{q}_1 d_1 + \mathbf{q}_2 d_2 + \cdots + \mathbf{q}_n d_n
```

**(4)** The Gram-Schmidt process takes independent $`\mathbf{a}_i`$ to orthogonal $`\mathbf{q}_i`$. Start with $`\mathbf{q}_1 = \mathbf{r}_1/\|\mathbf{r}_1\|`$.

**(5)** $`\mathbf{q}_i`$ is $`(\mathbf{a}_i - \text{its projection } \mathbf{p}_i) / \|\mathbf{a}_i - \mathbf{p}_i\|`$, where projection $`\mathbf{p}_i = (\mathbf{a}_i^T\mathbf{q}_1)\mathbf{q}_1 + \cdots + (\mathbf{a}_i^T\mathbf{q}_{i-1})\mathbf{q}_{i-1}`$.

**(6)** Each $`\mathbf{a}_i`$ will be a combination of $`\mathbf{q}_1`$ to $`\mathbf{q}_n`$. Then $`A = QR`$: orthogonal $`Q`$ and triangular $`R`$.

### 4.2 Goal: Orthogonal Columns

**i)** Orthogonal columns in $`A`$ are good.

```math
A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)
```

```math
\mathbf{a}_1^T\mathbf{a}_2 = 0, \quad \ldots, \quad \mathbf{a}_1^T\mathbf{a}_n = 0
```

```math
A^T A = \begin{pmatrix} \mathbf{a}_1^T\mathbf{a}_1 & 0 & \cdots & 0 \\ 0 & \mathbf{a}_2^T\mathbf{a}_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \mathbf{a}_n^T\mathbf{a}_n \end{pmatrix}
```

This is a **diagonal matrix**.

**ii)** Construct orthogonal vectors $`\mathbf{q}_i`$ through Gram-Schmidt process.

### 4.3 Orthonormal Vectors and Matrices

**Definition:** The $`n`$ vectors are **orthonormal** if:

```math
\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 0 & \text{when } i \neq j \\ 1 & \text{when } i = j \end{cases}
```

A matrix $`Q`$ with orthonormal columns has $`Q^T Q = I`$. Typically $`m > n`$.

Note: $`Q`$ is NOT required to be square. $`QQ^T \neq I`$ in general.

**Example:** $`Q \in \mathbb{R}^{3 \times 2}`$ orthonormal matrix:

```math
Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix}
```

```math
Q^T Q = \frac{1}{2}\begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = I_{2 \times 2}
```

```math
QQ^T = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{pmatrix} \neq I_{3 \times 3}
```

When $`Q`$ is square, $`Q^T Q = I \Longrightarrow Q^T = Q^{-1}`$ (due to the definition of inverse matrix).

### 4.4 Examples of Orthogonal Matrices

**Example 1: Rotation.** $`Q`$ rotates every vector in the plane by the angle $`\theta`$.

```math
Q = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}, \quad Q^T = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}
```

```math
Q^T Q = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
```

$`\mathbf{a}_1`$ and $`\mathbf{a}_2`$ are an orthonormal basis for the plane $`\mathbb{R}^2`$. The standard basis vectors $`\mathbf{i}, \mathbf{j}`$ are rotated through $`\theta`$. $`Q^{-1}`$ rotates vectors back through $`-\theta`$.

```math
Q^{-1} = \begin{pmatrix} \cos(-\theta) & -\sin(-\theta) \\ \sin(-\theta) & \cos(-\theta) \end{pmatrix}
```

**Example 2: Permutation.**

```math
\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} y \\ z \\ x \end{pmatrix}, \quad \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} y \\ x \end{pmatrix}
```

All columns are unit vectors. They are also orthogonal. The inverse of a permutation matrix is its transpose: $`Q^{-1} = Q^T`$. **Every permutation matrix is an orthogonal matrix.**

**Example 3: Reflection.**

Given a unit normal vector $`\mathbf{n}`$ ($`\|\mathbf{n}\| = 1`$), decompose $`\mathbf{u}`$:

```math
\mathbf{u} = \mathbf{u}_{\text{normal}} + \mathbf{u}_{\text{tangential}} = \mathbf{u}_\perp + \mathbf{u}_\parallel
```

```math
\mathbf{u}_\perp = (\mathbf{u} \cdot \mathbf{n})\mathbf{n} = \mathbf{n}(\mathbf{n}^T\mathbf{u}) = (\mathbf{n}\mathbf{n}^T)\mathbf{u}
```

```math
\mathbf{u}_\parallel = \mathbf{u} - \mathbf{u}_\perp = I\mathbf{u} - (\mathbf{n}\mathbf{n}^T)\mathbf{u} = (I - \mathbf{n}\mathbf{n}^T)\mathbf{u}
```

The reflection $`\mathbf{v}`$:

```math
\mathbf{v} = \mathbf{u}_\parallel - \mathbf{u}_\perp = (I - \mathbf{n}\mathbf{n}^T)\mathbf{u} - (\mathbf{n}\mathbf{n}^T)\mathbf{u} = (I - 2\mathbf{n}\mathbf{n}^T)\mathbf{u}
```

```math
Q = I - 2\mathbf{n}\mathbf{n}^T
```

**Verification that $`Q`$ is orthogonal:**

```math
Q^T Q = (I - 2\mathbf{n}\mathbf{n}^T)^T(I - 2\mathbf{n}\mathbf{n}^T) = (I - 2\mathbf{n}\mathbf{n}^T)(I - 2\mathbf{n}\mathbf{n}^T)
```

```math
= I - 4\mathbf{n}\mathbf{n}^T + 4(\mathbf{n}\mathbf{n}^T)(\mathbf{n}\mathbf{n}^T) = I - 4\mathbf{n}\mathbf{n}^T + 4\mathbf{n}\mathbf{n}^T = I
```

(since $`\mathbf{n}^T\mathbf{n} = \|\mathbf{n}\|^2 = 1`$).

**Key property:** Rotation, permutation, reflection matrices **preserve the length and the angle** of every vector:

```math
\|Q\mathbf{x}\|^2 = (Q\mathbf{x})^T(Q\mathbf{x}) = \mathbf{x}^T Q^T Q\mathbf{x} = \mathbf{x}^T I\mathbf{x} = \|\mathbf{x}\|^2
```

If $`Q`$ has orthonormal columns ($`Q^T Q = I`$), it leaves lengths unchanged.

### 4.5 Projections Using Orthonormal Bases

$`Q`$ replaces $`A`$.

Recall the projection of $`\mathbf{b}`$ onto $`C(A)`$ is $`\mathbf{p} = A(A^T A)^{-1}A^T\mathbf{b}`$.

When $`A = Q`$ and $`A^T A = Q^T Q = I`$:

```math
\mathbf{p} = QQ^T\mathbf{b}
```

```math
\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}
```

```math
\mathbf{p} = Q\hat{\boldsymbol{\alpha}} = \mathbf{q}_1 d_1 + \mathbf{q}_2 d_2 + \cdots + \mathbf{q}_n d_n
```

```math
= \mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \cdots + \mathbf{q}_n(\mathbf{q}_n^T\mathbf{b})
```

**Remark:** When $`Q \in \mathbb{R}^{n \times n}`$ (square), $`Q^T = Q^{-1}`$, so $`\hat{\boldsymbol{\alpha}} = Q^{-1}\mathbf{b}`$ and $`\mathbf{p} = QQ^{-1}\mathbf{b} = \mathbf{b}`$, $`P = I`$.

```math
\mathbf{b} = \mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \cdots + \mathbf{q}_n(\mathbf{q}_n^T\mathbf{b}) = QQ^T\mathbf{b}
```

$`QQ^T = I`$ is the foundation of **Fourier Series**, where they break $`\mathbf{b}`$ into perpendicular pieces. Then the linear combination of basis vectors puts $`\mathbf{b}`$ back together.

### 4.6 Example 4: Square Orthogonal Matrix

```math
Q = \frac{1}{3}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix}, \quad Q^T = Q
```

```math
Q^T Q = \frac{1}{9}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 9 & 0 & 0 \\ 0 & 9 & 0 \\ 0 & 0 & 9 \end{pmatrix} = I = QQ^T
```

Let $`\mathbf{b} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$:

```math
QQ^T\mathbf{b} = I\mathbf{b} = \mathbf{b}
```

```math
\mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \mathbf{q}_3(\mathbf{q}_3^T\mathbf{b})
```

```math
= \frac{1}{3}\begin{pmatrix} -1 \\ 2 \\ 2 \end{pmatrix} \cdot 2 + \frac{1}{3}\begin{pmatrix} 2 \\ -1 \\ 2 \end{pmatrix} \cdot 2 + \frac{1}{3}\begin{pmatrix} 2 \\ 2 \\ -1 \end{pmatrix} \cdot (-1) = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
```

### 4.7 The Gram-Schmidt Process

Given $`\mathbf{a}, \mathbf{b}, \mathbf{c}`$ vectors, create three orthogonal vectors $`\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3`$, which produce three orthonormal vectors:

```math
\mathbf{q}_1 = \frac{\mathbf{r}_1}{\|\mathbf{r}_1\|}, \quad \mathbf{q}_2 = \frac{\mathbf{r}_2}{\|\mathbf{r}_2\|}, \quad \mathbf{q}_3 = \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|}
```

**i)** $`\mathbf{r}_1 = \mathbf{a}`$, so $`\mathbf{q}_1 = \mathbf{r}_1/\|\mathbf{r}_1\|`$.

All vectors $`\mathbf{a}`$, $`\mathbf{r}_1`$, $`\mathbf{q}_1`$ are on a line.

**ii)** $`\mathbf{r}_2`$ must be perpendicular to $`\mathbf{r}_1`$.

Split $`\mathbf{b}`$ into $`\mathbf{b}_\parallel`$ and $`\mathbf{b}_\perp`$ to $`\mathbf{q}_1`$:

```math
\mathbf{b}_\parallel = (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1 = \left(\mathbf{b} \cdot \frac{\mathbf{r}_1}{\|\mathbf{r}_1\|}\right)\frac{\mathbf{r}_1}{\|\mathbf{r}_1\|} = \frac{1}{\|\mathbf{r}_1\|^2}(\mathbf{r}_1^T\mathbf{b})\,\mathbf{r}_1 = \left(\frac{\mathbf{r}_1^T\mathbf{b}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1
```

```math
\mathbf{b}_\perp = \mathbf{b} - \mathbf{b}_\parallel = \mathbf{b} - \left(\frac{\mathbf{r}_1^T\mathbf{b}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1 = \mathbf{r}_2
```

```math
\mathbf{q}_2 = \frac{\mathbf{r}_2}{\|\mathbf{r}_2\|}
```

**iii)** $`\mathbf{r}_3`$ must be perpendicular to $`\mathbf{r}_1, \mathbf{r}_2`$.

Subtract from every new vector its projection in the direction already set:

```math
\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2
```

```math
= \mathbf{c} - \left(\frac{\mathbf{r}_1^T\mathbf{c}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1 - \left(\frac{\mathbf{r}_2^T\mathbf{c}}{\mathbf{r}_2^T\mathbf{r}_2}\right)\mathbf{r}_2
```

```math
\mathbf{q}_3 = \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|}
```

All vectors $`\mathbf{a}, \mathbf{b}, \mathbf{c}, \mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3, \mathbf{q}_1, \mathbf{q}_2, \mathbf{q}_3`$ are in one subspace ($`\mathbb{R}^3`$).

### 4.8 Gram-Schmidt Example

```math
\mathbf{a} = \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix}, \quad \mathbf{c} = \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix}
```

**i)** $`\mathbf{r}_1 = \mathbf{a} = \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}`$, $`\|\mathbf{r}_1\| = \sqrt{2}`$:

```math
\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}
```

**ii)** $`\mathbf{r}_2 = \mathbf{b} - (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1`$:

```math
= \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix} - \frac{1}{\sqrt{2}} \cdot 2 \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} = \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix} - \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}
```

```math
\mathbf{q}_2 = \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}
```

**iii)** $`\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2`$:

```math
= \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix} - \frac{6}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} - \frac{-6}{\sqrt{6}} \cdot \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}
```

```math
= \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix} - 3\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + \begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}
```

```math
\mathbf{q}_3 = \frac{1}{\sqrt{3}}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}
```

### 4.9 Factorization $`A = QR`$

```math
(\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = (\mathbf{q}_1 \;\; \mathbf{q}_2 \;\; \mathbf{q}_3)\,R
```

The vectors $`\mathbf{a}, \mathbf{b}, \mathbf{c}`$ are combinations of $`\mathbf{q}_1, \mathbf{q}_2, \mathbf{q}_3`$:

```math
\mathbf{a} = (\mathbf{q}_1 \cdot \mathbf{a})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{a})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{a})\mathbf{q}_3
```

```math
\mathbf{b} = (\mathbf{q}_1 \cdot \mathbf{b})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{b})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{b})\mathbf{q}_3
```

```math
\mathbf{c} = (\mathbf{q}_1 \cdot \mathbf{c})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{c})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{c})\mathbf{q}_3
```

Note: Later $`\mathbf{q}`$'s are orthogonal to earlier $`\mathbf{a}`$'s, so many entries in $`R`$ are zero (below diagonal):

```math
(\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = (\mathbf{q}_1 \;\; \mathbf{q}_2 \;\; \mathbf{q}_3)\begin{pmatrix} \mathbf{q}_1^T\mathbf{a} & \mathbf{q}_1^T\mathbf{b} & \mathbf{q}_1^T\mathbf{c} \\ 0 & \mathbf{q}_2^T\mathbf{b} & \mathbf{q}_2^T\mathbf{c} \\ 0 & 0 & \mathbf{q}_3^T\mathbf{c} \end{pmatrix}
```

```math
A = QR
```

$`Q^T A = Q^T Q R = IR = R`$.

From LI $`\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n`$, Gram-Schmidt conducts orthonormal vectors $`\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n`$. The matrices with them satisfy $`A = QR`$. Then $`R = Q^T A`$ is **upper triangular** because later $`\mathbf{q}`$'s are orthogonal to earlier $`\mathbf{a}`$'s.

**From the example:**

```math
A = (\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = \begin{pmatrix} 1 & 2 & 3 \\ -1 & 0 & -3 \\ 0 & -2 & 3 \end{pmatrix}
```

```math
= \frac{1}{\sqrt{6}}\begin{pmatrix} \sqrt{3} & 1 & \sqrt{2} \\ -\sqrt{3} & 1 & \sqrt{2} \\ 0 & -2 & \sqrt{2} \end{pmatrix} \begin{pmatrix} \sqrt{2} & \sqrt{2} & \sqrt{18} \\ 0 & \sqrt{6} & \sqrt{6} \\ 0 & 0 & \sqrt{3} \end{pmatrix}
```

```math
Q \qquad \qquad R
```

The **lengths** of $`\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3`$ are the diagonal of $`R`$, which is positive.

Any $`m`$ by $`n`$ matrix $`A`$ with independent columns can be factored into $`A = QR`$. $`Q \in \mathbb{R}^{m \times n}`$ has orthogonal columns, $`R \in \mathbb{R}^{n \times n}`$ is upper triangular with positive diagonal.

### 4.10 Application: Least Squares via QR

```math
A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}
```

Substituting $`A = QR`$: $`A^T A = (QR)^T QR = R^T Q^T QR = R^T R`$.

```math
R^T R\hat{\boldsymbol{\alpha}} = R^T Q^T\mathbf{b}
```

Since $`R`$ is invertible:

```math
R\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b} \quad \text{(very fast)}
```

```math
\hat{\boldsymbol{\alpha}} = R^{-1}Q^T\mathbf{b}
```

---

<br>

## 5. The Pseudoinverse of a Matrix (4.5)

### 5.1 Two-Sided and One-Sided Inverses

**(1)** **Two-sided inverse:**

```math
A^{-1}A = AA^{-1} = I
```

**(2)** **One-sided inverse:**

```math
A^+ A = I \quad (A^+ \text{ is the left inverse of } A)
```

```math
AA^+ = I \quad (A^+ \text{ is the right inverse of } A)
```

Every matrix $`A`$ has a pseudoinverse $`A^+`$.

### 5.2 Three Cases with Identity-like Matrices

**(1)** $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}_{2 \times 2}`$, $`m = n = 2`$, $`\text{rank}(I) = 2 = r`$.

$`\dim C(I) = 2`$, $`\dim \mathcal{N}(I^T) = 0`$, $`\dim C(I^T) = 2`$, $`\dim \mathcal{N}(I) = 0`$.

All four subspaces are either full or trivial.

**(2)** $`I_L = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}_{2 \times 3}`$, $`m = 2`$, $`n = 3`$, $`\text{rank}(I_L) = 2 = r = m`$.

Every row is LI, but $`r = m < n`$. Nullspace has nontrivial elements.

$`\dim C(I_L) = 2`$, $`\dim \mathcal{N}(I_L^T) = 0`$, $`\dim C(I_L^T) = 2`$, $`\dim \mathcal{N}(I_L) = 1`$.

$`I_L I_R = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I_{2 \times 2}`$

$`I_L`$ = left inverse of $`I_R`$, $`I_R`$ = right inverse of $`I_L`$.

**(3)** $`I_R = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}_{3 \times 2}`$, $`m = 3`$, $`n = 2`$, $`\text{rank}(I_R) = 2 = r = n`$.

Every column is LI, but $`r < m`$.

$`\dim C(I_R) = 2`$, $`\dim \mathcal{N}(I_R^T) = 1`$, $`\dim C(I_R^T) = 2`$, $`\dim \mathcal{N}(I_R) = 0`$.

### 5.3 Left Inverse and Right Inverse

| | $`A^+A = I_{n \times n}`$ (Left Inverse) | $`AA^+ = I_{m \times m}`$ (Right Inverse) |
|:---|:---|:---|
| **Condition** | Full column rank: $`r = n < m`$ | Full row rank: $`r = m < n`$ |
| **Interpretation** | # unknowns < # equations | # equations < # unknowns |
| **Solution to $`A\mathbf{x} = \mathbf{b}`$** | 0 or 1 solution | Infinitely many solutions |
| **Nullspace** | $`\mathcal{N}(A) = \lbrace\mathbf{0}\rbrace`$ | $`\mathcal{N}(A^T) = \lbrace\mathbf{0}\rbrace`$ |
| **$`A^T A`$ or $`AA^T`$** | $`A^T A`$ is $`n \times n`$ and invertible | $`AA^T`$ is $`m \times m`$ and invertible |
| **Formula** | $`A^+ = (A^T A)^{-1}A^T`$ | $`A^+ = A^T(AA^T)^{-1}`$ |

**Left inverse case:** $`A^+A = I`$ describes the matrices in least squares. $`\hat{\boldsymbol{\alpha}} = A^+\mathbf{b}`$ is the solution to $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$, where $`\mathbf{b}`$ may not be in $`C(A)`$ and $`\hat{\boldsymbol{\alpha}}`$ is in $`C(A^T)`$.

**Right inverse case:** $`AA^+ = I`$. $`\mathbf{x}^+ = A^+\mathbf{b}`$ is the least length solution to $`A\mathbf{x} = \mathbf{b}`$ where $`\mathbf{b}`$ is in $`C(A)`$ and $`\mathbf{x}^+`$ is in $`C(A^T)`$.

### 5.4 The Pseudoinverse $`A^+`$ of a General Matrix $`A_{m \times n}`$

**Step 1:** Every vector $`\mathbf{b} \in \mathbb{R}^m`$ has two perpendicular parts: $`\mathbf{p}`$ and $`\mathbf{z}`$, such that $`\mathbf{b} = \mathbf{p} + \mathbf{z}`$.

**Step 2:** $`\mathbf{p} \in C(A)`$ and $`\mathbf{z} \in \mathcal{N}(A^T)`$ ($`\Leftrightarrow A^T\mathbf{z} = \mathbf{0}`$, $`A^+\mathbf{z} = \mathbf{0}`$).

**Step 3:** $`C(A^T) \ni \mathbf{x}^+ \longrightarrow A\mathbf{x}^+ = \mathbf{p}`$. Invert this part: $`\mathbf{x}^+ = A^+\mathbf{p}`$.

**Step 4:** $`A^+\mathbf{b} = A^+(\mathbf{p} + \mathbf{z}) = A^+\mathbf{p} + A^+\mathbf{z} = \mathbf{x}^+`$.

### 5.5 Big Picture for $`A^+`$

- $`A^+`$ shares the same four subspaces as $`A^T`$.
- $`A^+`$ inverts $`A`$ when it can: from column space to row space.

```math
C(A^T) \xrightarrow{A} C(A) \xrightarrow{A^+} C(A^T)
```

- $`A^+A`$ brings $`\mathbf{x}^+ \in C(A^T)`$ back to the same $`\mathbf{x}^+`$.
- $`A^+A`$ is an $`n \times n`$ projection matrix $`P_{\text{row}}`$ onto the row space of $`A`$, $`C(A^T)`$.

```math
P_{\text{row}} = A^+A = (A^+A)^2 = (A^+A)^T
```

- $`AA^+`$ is an $`m \times m`$ projection matrix $`P_{\text{col}}`$ onto the column space of $`A`$, $`C(A)`$.

```math
P_{\text{col}} = AA^+ = (AA^+)^2 = (AA^+)^T
```

### 5.6 Example 1: Finding $`A^+`$

```math
A = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ b_2 \end{pmatrix} = \mathbf{p} + \mathbf{z}
```

$`r = 1`$.

**i)** $`\mathbf{p} \in C(A)`$: $`A\mathbf{x} = \mathbf{p}`$: $`\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ 0 \end{pmatrix}`$

$`2x_1 = b_1`$, $`x_1 = b_1/2`$.

```math
\mathbf{x} = \begin{pmatrix} b_1/2 \\ x_2 \end{pmatrix} = \underbrace{\begin{pmatrix} b_1/2 \\ 0 \end{pmatrix}}_{\mathbf{x}^+} + \underbrace{\begin{pmatrix} 0 \\ x_2 \end{pmatrix}}_{\mathbf{x}_n}
```

**ii)** $`\mathbf{z} \in \mathcal{N}(A^T)`$: $`A^T\mathbf{z} = \mathbf{0}`$: $`\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} z_1 \\ z_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, so $`z_1 = 0`$, $`\mathbf{z} = \begin{pmatrix} 0 \\ z_2 \end{pmatrix}`$.

**iii)** $`A^+\mathbf{b} = A^+(\mathbf{p} + \mathbf{z}) = A^+\mathbf{p} + A^+\mathbf{z} = \mathbf{x}^+ = \begin{pmatrix} b_1/2 \\ 0 \end{pmatrix}`$

```math
= \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
```

```math
A^+ = \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}
```

**iv)** Verification:

```math
AA^+ = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = P_{\text{col}}
```

```math
A^+A = \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = P_{\text{row}}
```

### 5.7 Pseudoinverse of a Diagonal Matrix

```math
D = \begin{pmatrix} 2 & & \\ & 3 & \\ & & 0 \end{pmatrix}, \quad D^+ = \begin{pmatrix} 1/2 & & \\ & 1/3 & \\ & & 0 \end{pmatrix}
```

Invert the nonzero diagonal entries; leave zeros as zeros.

### 5.8 Pseudoinverse via Factorization

Let $`U, V`$ be invertible matrices. Let $`A = UDV^T`$. Then:

```math
A^+ = (UDV^T)^+ = V^{-T}D^+U^{-1}
```

**Example (full column rank):** Let $`A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`\text{rank}(A) = 1 = r = n < m`$. Full column rank, $`\exists`$ left inverse.

$`A^+A = 1`$. $`A^+ = (A^T A)^{-1}A^T = ((1\;1)\begin{pmatrix} 1 \\ 1 \end{pmatrix})^{-1}(1\;1) = \frac{1}{2}(1\;1)`$.

**Example (full row rank):** Let $`A = (1\;1)`$, $`\text{rank}(A) = 1 = m < n`$. Full row rank, $`\exists`$ right inverse.

$`AA^+ = 1`$. $`A^+ = A^T(AA^T)^{-1} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}((1\;1)\begin{pmatrix} 1 \\ 1 \end{pmatrix})^{-1} = \frac{1}{2}\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$.

### 5.9 Important Action of $`A`$: Row Space to Column Space

```math
C(A^T) \xrightarrow{A} C(A)
```

```math
C(A) \xrightarrow{A^+} C(A^T)
```

**i)** $`\mathbf{x}_1, \mathbf{x}_2 \in C(A^T)`$, $`\mathbf{x}_1 \neq \mathbf{x}_2 \Rightarrow A\mathbf{x}_1 \neq A\mathbf{x}_2`$.

**Proof:** Suppose $`A\mathbf{x}_1 = \mathbf{b}`$ and $`A\mathbf{x}_2 = \mathbf{b}`$. Subtract: $`A(\mathbf{x}_1 - \mathbf{x}_2) = \mathbf{0}`$. By definition of nullspace, $`\mathbf{x}_1 - \mathbf{x}_2 \in \mathcal{N}(A)`$. But $`\mathbf{x}_1 - \mathbf{x}_2 \in C(A^T)`$. Since $`\mathcal{N}(A)`$ is orthogonal to $`C(A^T)`$, $`\mathbf{x}_1 - \mathbf{x}_2 = \mathbf{0}`$. This contradicts $`\mathbf{x}_1 \neq \mathbf{x}_2`$. Therefore $`A\mathbf{x}_1 \neq A\mathbf{x}_2`$. $`\square`$

**ii)** $`\forall\,\mathbf{b} \in C(A)`$, $`\exists!\,\mathbf{x}^+ \in C(A^T)`$ such that $`A^+\mathbf{b} = \mathbf{x}^+`$.

### 5.10 Pseudoinverse $`A^+ = R^+C^+`$ of $`A = CR`$

Given $`A = CR`$ where $`C`$ has $`r`$ independent columns (full column rank) and $`R`$ has $`r`$ independent rows (full row rank):

- $`C`$ full column rank: $`\exists (C^T C)^{-1}`$, $`\exists`$ left inverse of $`C`$: $`C^+C = I`$, $`C^+ = (C^T C)^{-1}C^T`$.
- $`R`$ full row rank: $`\exists (RR^T)^{-1}`$, $`\exists`$ right inverse of $`R`$: $`RR^+ = I`$, $`R^+ = R^T(RR^T)^{-1}`$.

```math
A^+ = R^+C^+ = R^T(RR^T)^{-1}(C^T C)^{-1}C^T
```

```math
= R^T(C^T A R^T)^{-1}C^T
```

(since $`C^T C R R^T = C^T(CR)R^T = C^T A R^T`$).

### 5.11 Example: Pseudoinverse of a $`5 \times 4`$ Matrix

```math
A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ -1 & 0 & 1 & 0 \\ 0 & -1 & 1 & 0 \\ -1 & 0 & 0 & 1 \\ 0 & -1 & 0 & 1 \end{pmatrix}_{5 \times 4}
```

Row reduction steps: $`R_2 - R_1`$, $`R_4 - R_1`$, then $`R_3 - R_2`$, $`R_4 - R_2`$, $`R_5 - R_2`$, swap $`R_3`$ and $`R_5`$, $`R_4 - R_3`$, then $`R_1 + R_2`$, $`R_2 + R_3`$, $`\times(-1)`$:

```math
R_0 = \begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}
```

Columns 1, 2, 3 are LI. $`\text{rank}(A) = 3`$.

```math
R = \begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{pmatrix}
```

```math
A = CR = \begin{pmatrix} -1 & 1 & 0 \\ -1 & 0 & 1 \\ 0 & -1 & 1 \\ -1 & 0 & 0 \\ 0 & -1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{pmatrix}
```

$`A^+ = (CR)^+ = R^+C^+ = R^T(C^T AR^T)^{-1}C^T`$.

It turns out $`C^T AR^T = \begin{pmatrix} 4 & 0 & 0 \\ 0 & 4 & 0 \\ -1 & -1 & 2 \end{pmatrix}`$.

```math
A^+ = \frac{1}{8}\begin{pmatrix} -2 & -2 & 0 & -2 & 0 \\ 2 & 0 & -2 & 0 & -2 \\ 0 & 3 & 3 & -1 & -1 \\ 0 & -1 & -1 & 3 & 3 \end{pmatrix}
```

The mapping through $`A = CR`$:

```math
C(A^T) \ni \mathbf{x} \xrightarrow{R} R\mathbf{x} \xrightarrow{C} CR\mathbf{x} = A\mathbf{x} \in C(A)
```

$`R`$ has full row rank: $`\mathcal{N}(R^T) = \lbrace\mathbf{0}\rbrace`$. $`C`$ has full column rank: $`\mathcal{N}(C) = \lbrace\mathbf{0}\rbrace`$.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Orthogonal vectors | $`\mathbf{v}^T\mathbf{w} = 0`$ implies $`\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2`$ |
| Fundamental subspace orthogonality | $`C(A^T) \perp \mathcal{N}(A)`$ in $`\mathbb{R}^n`$; $`C(A) \perp \mathcal{N}(A^T)`$ in $`\mathbb{R}^m`$ |
| Orthogonal complement | $`\dim V + \dim V^\perp = n`$; every $`\mathbf{x} = \mathbf{x}_r + \mathbf{x}_n`$ |
| Projection onto a line | $`\mathbf{p} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\mathbf{a}`$; projection matrix $`P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}`$ |
| Projection onto a subspace | $`\mathbf{p} = A(A^T A)^{-1}A^T\mathbf{b}`$; $`P = A(A^T A)^{-1}A^T`$; $`P^2 = P = P^T`$ |
| $`A^T A`$ invertibility | $`A^T A`$ is invertible $`\iff`$ $`A`$ has LI columns $`\iff`$ $`\mathcal{N}(A) = \lbrace\mathbf{0}\rbrace`$ |
| Least squares solution | $`\hat{\boldsymbol{\alpha}} = (A^T A)^{-1}A^T\mathbf{b}`$ minimizes $`E = \|A\mathbf{x} - \mathbf{b}\|^2`$ |
| Normal equations | $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$ (geometry, algebra, or calculus derivation) |
| Line fitting | $`A = \begin{pmatrix} x_1 & 1 \\ \vdots & \vdots \\ x_m & 1 \end{pmatrix}`$; $`A^T A = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & m \end{pmatrix}`$ |
| Parabola fitting | $`A`$ has columns $`(1, 1, \ldots)`$, $`(t_1, t_2, \ldots)`$, $`(t_1^2, t_2^2, \ldots)`$ |
| Orthonormal vectors | $`\mathbf{q}_i^T\mathbf{q}_j = \delta_{ij}`$; $`Q^T Q = I`$ |
| Orthogonal matrix (square $`Q`$) | $`Q^T Q = QQ^T = I`$; $`Q^T = Q^{-1}`$; preserves lengths and angles |
| Examples of orthogonal matrices | Rotation, permutation, reflection ($`Q = I - 2\mathbf{n}\mathbf{n}^T`$) |
| Projection with orthonormal basis | $`\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}`$; $`\mathbf{p} = QQ^T\mathbf{b}`$; if $`Q`$ square, $`P = I`$ |
| Gram-Schmidt process | $`\mathbf{r}_1 = \mathbf{a}`$; $`\mathbf{r}_2 = \mathbf{b} - (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1`$; $`\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2`$; normalize each |
| QR Factorization | $`A = QR`$; $`Q`$ orthogonal columns, $`R`$ upper triangular with positive diagonal |
| Least squares via QR | $`R\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}`$ (much faster than normal equations) |
| Pseudoinverse $`A^+`$ | Inverts $`A`$ from $`C(A) \to C(A^T)`$; $`A^+\mathbf{z} = \mathbf{0}`$ for $`\mathbf{z} \in \mathcal{N}(A^T)`$ |
| Left inverse ($`r = n < m`$) | $`A^+ = (A^T A)^{-1}A^T`$; $`A^+A = I`$; least squares solution |
| Right inverse ($`r = m < n`$) | $`A^+ = A^T(AA^T)^{-1}`$; $`AA^+ = I`$; minimum length solution |
| $`A^+A`$ and $`AA^+`$ | $`A^+A = P_{\text{row}}`$ (projection onto row space); $`AA^+ = P_{\text{col}}`$ (projection onto column space) |
| Pseudoinverse via $`A = CR`$ | $`A^+ = R^+C^+ = R^T(C^T AR^T)^{-1}C^T`$ |

---
