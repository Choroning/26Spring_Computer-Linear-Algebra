# Chapter 8 Lecture — Linear Transformations

> **Last Updated:** 2026-03-31

---

<br>

## Table of Contents

- [1. Idea of a Linear Transformation (8.1)](#1-idea-of-a-linear-transformation-81)
  - [1.1 Definition and Linearity Condition](#11-definition-and-linearity-condition)
  - [1.2 Formal Definition of a Linear Transformation](#12-formal-definition-of-a-linear-transformation)
  - [1.3 Non-Linear Example: Shift Transformation and Affine Mappings](#13-non-linear-example-shift-transformation-and-affine-mappings)
  - [1.4 Examples of Linear and Non-Linear Transformations](#14-examples-of-linear-and-non-linear-transformations)
  - [1.5 Linear Transformation Determined by Basis](#15-linear-transformation-determined-by-basis)
  - [1.6 Geometric Interpretation: Lines to Lines](#16-geometric-interpretation-lines-to-lines)
  - [1.7 Linear Transformations in Calculus (Derivative)](#17-linear-transformations-in-calculus-derivative)
  - [1.8 Example 5: Integration is Also Linear](#18-example-5-integration-is-also-linear)
  - [1.9 Example 6: Projection onto z=1 Plane (Non-Linear)](#19-example-6-projection-onto-z1-plane-non-linear)
  - [1.10 Example 7: Invertible Matrix and Range/Kernel](#110-example-7-invertible-matrix-and-rangekernel)
- [2. The Matrix of a Linear Transformation (8.2)](#2-the-matrix-of-a-linear-transformation-82)
  - [2.1 Key Ideas Summary](#21-key-ideas-summary)
  - [2.2 Choice of Bases and the Standard Matrix](#22-choice-of-bases-and-the-standard-matrix)
  - [2.3 Example 1: Standard Basis in R^2 to R^3](#23-example-1-standard-basis-in-r2-to-r3)
  - [2.4 Construction of the Matrix for T](#24-construction-of-the-matrix-for-t)
  - [2.5 Change of Basis: Matrix B](#25-change-of-basis-matrix-b)
  - [2.6 Example 3: Derivative Matrix for Polynomials](#26-example-3-derivative-matrix-for-polynomials)
  - [2.7 Example 4: Integral Matrix (Pseudoinverse of Derivative)](#27-example-4-integral-matrix-pseudoinverse-of-derivative)
  - [2.8 Matrix Product AB Matches Transformation TS](#28-matrix-product-ab-matches-transformation-ts)
  - [2.9 Example 5: Rotation Composition](#29-example-5-rotation-composition)
  - [2.10 Example 6: Inverse Rotation](#210-example-6-inverse-rotation)
- [3. The Search for a Good Basis (8.3)](#3-the-search-for-a-good-basis-83)
  - [3.1 Key Ideas Summary](#31-key-ideas-summary)
  - [3.2 The Change of Basis Formula](#32-the-change-of-basis-formula)
  - [3.3 Best Basis 1: Eigenvectors (Diagonalization)](#33-best-basis-1-eigenvectors-diagonalization)
  - [3.4 Best Basis 2: Singular Vectors (SVD)](#34-best-basis-2-singular-vectors-svd)
  - [3.5 Best Basis 3: Generalized Eigenvectors (Jordan Form)](#35-best-basis-3-generalized-eigenvectors-jordan-form)
  - [3.6 Jordan Form: Structure and Definition](#36-jordan-form-structure-and-definition)
  - [3.7 Example: 2x2 Jordan Form](#37-example-2x2-jordan-form)
  - [3.8 Example: 3x3 Jordan Form](#38-example-3x3-jordan-form)
  - [3.9 Basis for Function Space](#39-basis-for-function-space)
  - [3.10 Orthogonal Bases for Function Space](#310-orthogonal-bases-for-function-space)
  - [3.11 Constructing the Legendre Basis via Gram-Schmidt](#311-constructing-the-legendre-basis-via-gram-schmidt)
- [Summary](#summary)

---

<br>

## 1. Idea of a Linear Transformation (8.1)

### 1.1 Definition and Linearity Condition

**Chapter overview:**

- **8.1** Idea of a Linear Transformation
- **8.2** The Matrix of a Linear Transformation
- **8.3** The Search for a Good Basis

A **linear transformation** $T$ takes vectors $\mathbf{u}$ to vectors $T(\mathbf{u})$:

$$T: \mathbf{u} \longrightarrow T(\mathbf{u})$$

**Linearity requires:**

$$T(c\mathbf{u} + d\mathbf{w}) = c\,T(\mathbf{u}) + d\,T(\mathbf{w})$$

**Note:** $T(\mathbf{0}) = \mathbf{0}$.

So $T(\mathbf{u}) = \mathbf{u} + \mathbf{u}_0$ is **NOT** linear (because $T(\mathbf{0}) = \mathbf{u}_0 \neq \mathbf{0}$).

The input vector $\mathbf{u}$ and output $T(\mathbf{u})$ can be in $\mathbb{R}^n$ or matrix space or function space.

If $A$ is $m$ by $n$, then $T(\mathbf{u}) = A\mathbf{u}$ is linear, from the input space $\mathbb{R}^n$ to the output space $\mathbb{R}^m$.

The **derivative** $T(f) = \dfrac{df}{dx}$ is linear.

The **integral** $T^+(f) = \displaystyle\int_0^x f(t)\,dt$ is its pseudoinverse.

The **product** $ST$ of two linear transformations is still linear:

$$(ST)(\mathbf{u}) = S(T(\mathbf{u}))$$

---

### 1.2 Formal Definition of a Linear Transformation

$A: \mathbf{u} \longrightarrow A\mathbf{u}$

$A$ transforms $\mathbf{u}$ into another vector $A\mathbf{u}$.

Similar to $f: x \longrightarrow f(x)$, a transformation $T$ takes a vector $\mathbf{u}$ and maps it to another vector $T(\mathbf{u})$:

$$T: \mathbf{u} \longrightarrow T(\mathbf{u})$$

Start with $A$:

$$A: \mathbf{u} \longrightarrow A\mathbf{u}$$
$$A: \mathbf{w} \longrightarrow A\mathbf{w}$$

What happens to $\mathbf{y} = \mathbf{u} + \mathbf{w}$?

$$A\mathbf{y} = A(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w}$$

Matrix multiplication $T(\mathbf{u}) = A\mathbf{u}$ gives a linear transformation.

A transformation $T$ assigns an output $T(\mathbf{u})$ to each input vector $\mathbf{u}$ in $V$.

**The transformation is linear if it meets these requirements for all $\mathbf{u}$ and $\mathbf{w} \in V$:**

**(a)** $T(\mathbf{u} + \mathbf{w}) = T(\mathbf{u}) + T(\mathbf{w})$

**(b)** $T(c\mathbf{u}) = c\,T(\mathbf{u}) \quad \forall\, c \in \mathbb{C}$

We say $T$ is **additive** and **homogeneous**.

$$T(\mathbf{0}) = T(0\mathbf{u}) = 0\,T(\mathbf{u}) = \mathbf{0}$$

We can combine (a) and (b) into:

$$T(c\mathbf{u} + d\mathbf{w}) = c\,T(\mathbf{u}) + d\,T(\mathbf{w})$$

**Generalize:**

$$T(c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_n\mathbf{x}_n) = c_1\,T(\mathbf{x}_1) + c_2\,T(\mathbf{x}_2) + \cdots + c_n\,T(\mathbf{x}_n)$$

---

### 1.3 Non-Linear Example: Shift Transformation and Affine Mappings

**Example.** Define $T$ that takes a vector $\mathbf{u} \in V$ and adds $\mathbf{u}_0 \in V$:

$$T: \mathbf{u} \longrightarrow \mathbf{u} + \mathbf{u}_0$$
$$T: \mathbf{w} \longrightarrow \mathbf{w} + \mathbf{u}_0$$

Check additivity:

$$T(\mathbf{u} + \mathbf{w}) \stackrel{?}{=} T(\mathbf{u}) + T(\mathbf{w})$$

$$\mathbf{u} + \mathbf{w} + \mathbf{u}_0 \neq \mathbf{u} + \mathbf{u}_0 + \mathbf{w} + \mathbf{u}_0$$

Therefore $T$ is **NOT** linear.

**Linear plus shift transformation:**

$$T(\mathbf{u}) = A\mathbf{u} + \mathbf{u}_0$$

is called **"affine"**.

In computer graphics, affine mapping is used:
- Start with shape $\mathbf{x}$
- Rotate: $\mathbf{y} = R\mathbf{x}$ (linear)
- Translate: $\mathbf{z} = \mathbf{y} + \begin{pmatrix} 0 \\ u_0 \end{pmatrix}$ (shift)

---

### 1.4 Examples of Linear and Non-Linear Transformations

**Example 1 (Linear — Dot product):**

$$\mathbf{a} = \begin{pmatrix} 1 \\ 3 \\ 4 \end{pmatrix}, \quad \mathbf{u} = \begin{pmatrix} u_1 \\ u_2 \\ u_3 \end{pmatrix}$$

$$T(\mathbf{u}) = \mathbf{a} \cdot \mathbf{u} = \mathbf{a}^T\mathbf{u} = (\mathbf{a}^T)\mathbf{u} = u_1 + 3u_2 + 4u_3$$

**Dot products are linear.**

---

**Example 2 (Non-Linear — Length/Norm):**

$T(\mathbf{u}) = \|\mathbf{u}\|$ is **NOT** linear.

**Check (i) — Additivity:**

$$T(\mathbf{u} + \mathbf{w}) \stackrel{?}{=} T(\mathbf{u}) + T(\mathbf{w})$$

$$\|\mathbf{u} + \mathbf{w}\| \neq \|\mathbf{u}\| + \|\mathbf{w}\|$$

**Check (ii) — Homogeneity:**

$$T(c\mathbf{u}) \stackrel{?}{=} c\,T(\mathbf{u})$$

$$T(-\mathbf{u}) \stackrel{?}{=} -T(\mathbf{u})$$

$$\|-\mathbf{u}\| = \|\mathbf{u}\| \neq -\|\mathbf{u}\|$$

Fails homogeneity (norm is always non-negative).

---

**Example 3 (Linear — Rotation):**

$T$ is the rotation matrix that rotates every vector by $30°$.

$$T: \mathbf{u} \longrightarrow \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\mathbf{u} \quad \text{with } \theta = 30°$$

where the rotation matrix is $R$.

**Check (i):**

$$T(\mathbf{u} + \mathbf{w}) = R(\mathbf{u} + \mathbf{w}) = R\mathbf{u} + R\mathbf{w} = T(\mathbf{u}) + T(\mathbf{w})$$

**Check (ii):**

$$T(c\mathbf{u}) = R(c\mathbf{u}) = cR\mathbf{u} = c\,T(\mathbf{u})$$

**Recall the generalization:**

$$T(c_1\mathbf{u}_1 + c_2\mathbf{u}_2 + \cdots + c_n\mathbf{u}_n) = c_1\,T(\mathbf{u}_1) + c_2\,T(\mathbf{u}_2) + \cdots + c_n\,T(\mathbf{u}_n)$$

---

### 1.5 Linear Transformation Determined by Basis

Suppose you know $T(\mathbf{u})$ for all vectors $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n$ in a basis.

Let $\mathbf{y} = c_1\mathbf{u}_1 + \cdots + c_n\mathbf{u}_n$.

Then you know $T(\mathbf{y})$ for all $\mathbf{y}$ in the space.

> A linear transformation is completely determined by its action on a basis.

---

### 1.6 Geometric Interpretation: Lines to Lines

Linear transformations map:
- **Lines to lines**
- **Triangles to triangles**
- **Equally spaced points go to equally spaced points**

---

### 1.7 Linear Transformations in Calculus (Derivative)

**Example 4.** $T(u) = \dfrac{d}{dt}(u)$

**Check (i):**

$$T(cu + dw) = \frac{d}{dt}(cu + dw) = c\frac{du}{dt} + d\frac{dw}{dt} = c\,T(u) + d\,T(w) \quad \text{Linear}$$

where $c, d \in \mathbb{R}$.

**Nullspace of $T$:**

$$T(u) = \frac{d}{dt}u = 0$$

This only happens when $u$ is a constant $c$.

$$dc \in \mathcal{N}(T), \quad d \in \mathbb{R}$$

The nullspace is a **line** in function space.

**Column space of $T$:**

Let $u = a + bt + ct^2$, then:

$$T(u) = b + 2ct, \quad \text{a linear function.}$$

**Dimensions:**

$$\dim C(T) + \dim \mathcal{N}(T) = 2 + 1 = 3$$

**Matrix for the derivative $T = \dfrac{d}{dt}$:**

View $u = (1 \;\; t \;\; t^2) \begin{pmatrix} a \\ b \\ c \end{pmatrix}$.

Let $\mathbf{u}_1 = 1$, $\mathbf{u}_2 = t$, $\mathbf{u}_3 = t^2$.

$$\frac{d}{dt}\mathbf{u}_1 = \frac{d}{dt}(1) = 0$$

$$\frac{d}{dt}\mathbf{u}_2 = \frac{d}{dt}(t) = 1$$

$$\frac{d}{dt}\mathbf{u}_3 = \frac{d}{dt}(t^2) = 2t$$

$$T: \begin{pmatrix} a \\ b \\ c \end{pmatrix} \longrightarrow \begin{pmatrix} b \\ 2c \end{pmatrix}$$

$$a + bx + cx^2 \longrightarrow b + 2cx$$

Linear transformation $\dfrac{dy}{dx}$ is connected to $A\mathbf{u}$:

$$\begin{pmatrix} b \\ 2c \end{pmatrix} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} a \\ b \\ c \end{pmatrix}$$

The matrix is $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$.

---

### 1.8 Example 5: Integration is Also Linear

Integration $T^+$ is also linear:

$$\int_0^x (D + Ex)\,dx = Dx + \frac{1}{2}Ex^2$$

$$T^+: \begin{pmatrix} D \\ E \end{pmatrix} \longrightarrow \begin{pmatrix} 0 \\ D \\ \frac{1}{2}E \end{pmatrix}$$

The input $\mathbf{u} = D + Ex$ and the output is $Dx + \frac{1}{2}Ex^2$.

$$\begin{pmatrix} 0 \\ D \\ \frac{1}{2}E \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix}\begin{pmatrix} D \\ E \end{pmatrix}$$

The integral matrix is $A^+$.

**Products $A^+A$ and $AA^+$:**

$$A^+A = \begin{pmatrix} 0 & 0 \\ 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix}\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$AA^+ = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} 0 & 0 \\ 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$$

---

### 1.9 Example 6: Projection onto z=1 Plane (Non-Linear)

Project a vector $\mathbf{u} \in \mathbb{R}^3$ onto the horizontal plane $z = 1$:

$$T: \begin{pmatrix} x \\ y \\ z \end{pmatrix} \longrightarrow \begin{pmatrix} x \\ y \\ 1 \end{pmatrix}$$

$T(\mathbf{u})$ is **NOT** linear because:

$$T(\mathbf{0}) = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \longrightarrow \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} \neq \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

---

### 1.10 Example 7: Invertible Matrix and Range/Kernel

Suppose $A$ is invertible:

$$T(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w} = T(\mathbf{u}) + T(\mathbf{w})$$

The inverse transformation: $T^{-1}(\mathbf{u}) = A^{-1}(\mathbf{u})$

$$T^{-1}(T(\mathbf{u})) = A^{-1}(A\mathbf{u}) = \mathbf{u}$$

If $T(\mathbf{u}) = A\mathbf{u}$, $S(\mathbf{u}) = B\mathbf{u}$, then $T \circ S(\mathbf{u})$ corresponds to $AB\mathbf{u}$.

**Q:** Are all linear transformations from $V = \mathbb{R}^n$ to $W = \mathbb{R}^m$ produced by matrices? **Yes.**

$A\mathbf{u}$ is the **column space**.

The **null space** of $A$ contains all the input for which $A\mathbf{u} = \mathbf{0}$.

$$\text{Range of } T = T(\mathbf{u}) \quad \longleftrightarrow \quad \text{column space } A\mathbf{u}$$

$$\text{Kernel of } T = \text{all inputs for which } T(\mathbf{u}) = \mathbf{0} \quad \longleftrightarrow \quad \text{nullspace of } A$$

---

<br>

## 2. The Matrix of a Linear Transformation (8.2)

### 2.1 Key Ideas Summary

1. **Linearity** tells us all $T(\mathbf{u})$ if we know $T(\mathbf{u}_1), \ldots, T(\mathbf{u}_n)$ for an input basis $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n$.

2. **Column $j$** in the matrix for $T$ comes from applying $T$ to the input basis vector $\mathbf{u}_j$.

3. Write $T(\mathbf{u}_j) = a_{1j}\,\mathbf{w}_1 + \cdots + a_{mj}\,\mathbf{w}_m = \displaystyle\sum_{i=1}^{m} a_{ij}\,\mathbf{w}_i$ in the output basis of $\mathbf{w}$'s. Those $a_{ij}$ go into column $j$.

4. The matrix $T(\mathbf{x}) = A\mathbf{x}$ is $A$, if the input and the output bases are the columns of $I_{n \times n}$ and $I_{m \times m}$.

5. When the bases change to $\mathbf{v}$'s and $\mathbf{w}$'s, the matrix for the same $T$ changes from $A$ to $W^{-1}AV$.

6. **Best bases:** $V = W$ = eigenvectors and $V, W$ = singular vectors: $A \to \Lambda$ and $\Sigma$.

---

### 2.2 Choice of Bases and the Standard Matrix

A transformation $T$ maps from input space $V$ to output space $W$:

$$V \xrightarrow{T} W$$

Choose $V = \mathbb{R}^n$, $W = \mathbb{R}^m$. Then the matrix $A$ for this transformation will be $m$ by $n$.

The choice of bases in $V$ and $W$ will decide $A$.

The standard basis vectors for $\mathbb{R}^n$ and $\mathbb{R}^m$ are the columns of $I$.

This choice leads to a **standard matrix**: $T(\mathbf{u}) = A\mathbf{u}$.

We can choose other bases for $V$ and $W$.
$\Rightarrow$ The same transformation $T$ is represented by other matrices.

**Q.** How should we choose the bases that give the **best matrix** for $T$?

---

### 2.3 Example 1: Standard Basis in R^2 to R^3

$T: \mathbb{R}^2 \to \mathbb{R}^3$

$$\mathbf{u}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix} \longrightarrow T(\mathbf{u}_1) = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix}$$

$$\mathbf{u}_2 = \begin{pmatrix} 0 \\ 1 \end{pmatrix} \longrightarrow T(\mathbf{u}_2) = \begin{pmatrix} 5 \\ 5 \\ 5 \end{pmatrix}$$

Express this linear transformation as matrix $A$:

$$A = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) \end{pmatrix} = \begin{pmatrix} 2 & 5 \\ 3 & 5 \\ 4 & 5 \end{pmatrix}$$

The outputs go into the **columns of $A$**. Then $T(\mathbf{u}) = A\mathbf{u}$.

**Verification:**

$$T(\mathbf{u}_1 + \mathbf{u}_2) = A(\mathbf{u}_1 + \mathbf{u}_2) = A\mathbf{u}_1 + A\mathbf{u}_2 = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix} + \begin{pmatrix} 5 \\ 5 \\ 5 \end{pmatrix} = A\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 7 \\ 8 \\ 9 \end{pmatrix}$$

**General formula:**

$$T(c_1\mathbf{u}_1 + c_2\mathbf{u}_2) = c_1\,T(\mathbf{u}_1) + c_2\,T(\mathbf{u}_2) = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = A\mathbf{c}$$

---

### 2.4 Construction of the Matrix for T

Construct a matrix for **any** linear transformation.

$$V \xrightarrow{T} W, \quad \dim V = n, \quad \dim W = m$$

Choose a basis $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n$ for $V$.

Choose a basis $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_m$ for $W$.

$A$ will be $m$ by $n$.

$T(\mathbf{u}_1)$ is a combination of the output basis for $W$:

$$T(\mathbf{u}_1) = d_1\mathbf{w}_1 + d_2\mathbf{w}_2 + \cdots + d_m\mathbf{w}_m = a_{11}\mathbf{w}_1 + a_{21}\mathbf{w}_2 + \cdots + a_{m1}\mathbf{w}_m = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{11} \\ a_{21} \\ \vdots \\ a_{m1} \end{pmatrix}$$

$$T(\mathbf{u}_2) = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{12} \\ a_{22} \\ \vdots \\ a_{m2} \end{pmatrix}$$

$$\vdots$$

$$T(\mathbf{u}_n) = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{1n} \\ a_{2n} \\ \vdots \\ a_{mn} \end{pmatrix}$$

Combining:

$$\begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) & \cdots & T(\mathbf{u}_n) \end{pmatrix} = A$$

$$= (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

> **The $j$th column of $A$ is found by applying $T$ to the $j$th basis vector $\mathbf{u}_j$, which is a linear combination of output basis vectors.**

Suppose we know outputs $T(\mathbf{u})$ for the input basis vectors $\mathbf{u}_1$ to $\mathbf{u}_n$:

$$T(\mathbf{u}_1),\; T(\mathbf{u}_2),\; \ldots,\; T(\mathbf{u}_n)$$

$$A = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) & \cdots & T(\mathbf{u}_n) \end{pmatrix}$$

$$A\mathbf{c} = c_1\,T(\mathbf{u}_1) + c_2\,T(\mathbf{u}_2) + \cdots + c_n\,T(\mathbf{u}_n)$$

---

### 2.5 Change of Basis: Matrix B

**Example 2.** $V = W = \mathbb{R}^2$

Suppose $T(\mathbf{u}) = \mathbf{u}$ (Identity transformation).

**Input basis** ($V$): Standard basis $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$.

**Output basis** ($W$): Rotated basis $\mathbf{w}_1 = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}$, $\mathbf{w}_2 = \begin{pmatrix} -\sin\theta \\ \cos\theta \end{pmatrix}$.

Express $\mathbf{u}$ in both bases:

$$\mathbf{u} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 = (\mathbf{v}_1 \;\; \mathbf{v}_2)\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = V\mathbf{c}$$

$$\mathbf{u} = d_1\mathbf{w}_1 + d_2\mathbf{w}_2 = (\mathbf{w}_1 \;\; \mathbf{w}_2)\begin{pmatrix} d_1 \\ d_2 \end{pmatrix} = W\mathbf{d}$$

Since both represent the same vector:

$$\boxed{V\mathbf{c} = W\mathbf{d}}$$

$$\mathbf{c} = V^{-1}W\mathbf{d} \qquad \mathbf{d} = W^{-1}V\mathbf{c}$$

$$B = W^{-1}V$$

**In this example:**

$$\mathbf{d} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}^{-1}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\mathbf{c}$$

$$= \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\mathbf{c}$$

$$= \begin{pmatrix} \cos(-\theta) & -\sin(-\theta) \\ \sin(-\theta) & \cos(-\theta) \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\mathbf{c}$$

When $V = I$, the coefficients $\mathbf{d}$ in $W$ become:

$$\mathbf{d} = W^{-1}\begin{pmatrix} x \\ y \end{pmatrix}$$

---

### 2.6 Example 3: Derivative Matrix for Polynomials

$T(u) = \dfrac{du}{dx}$

Input space: $V = \text{span}\{1, x, x^2, x^3\}$ with basis $\mathbf{u}_1 = 1,\; \mathbf{u}_2 = x,\; \mathbf{u}_3 = x^2,\; \mathbf{u}_4 = x^3$.

Output space: $W = \text{span}\{1, x, x^2\}$ with basis $\mathbf{w}_1 = 1,\; \mathbf{w}_2 = x,\; \mathbf{w}_3 = x^2$.

Apply $T$ to each basis vector:

$$T(\mathbf{u}_1) = T(1) = 0$$
$$T(\mathbf{u}_2) = T(x) = 1 = \mathbf{w}_1$$
$$T(\mathbf{u}_3) = T(x^2) = 2x = 2\mathbf{w}_2$$
$$T(\mathbf{u}_4) = T(x^3) = 3x^2 = 3\mathbf{w}_3$$

$$A = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) & T(\mathbf{u}_3) & T(\mathbf{u}_4) \end{pmatrix}$$

$$= (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \mathbf{w}_3)\begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$$

The **derivative matrix** is:

$$D = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$$

For input $u = c_1 + c_2 x + c_3 x^2 + c_4 x^3$ with coefficient vector $\mathbf{c} = (c_1, c_2, c_3, c_4)^T$:

$$T(u) = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \mathbf{w}_3) \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \\ c_3 \\ c_4 \end{pmatrix} = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \mathbf{w}_3)\begin{pmatrix} c_2 \\ 2c_3 \\ 3c_4 \end{pmatrix}$$

$\mathbf{d} = D\mathbf{c}$ produces the coefficients in the $T(u)$ combination.

---

### 2.7 Example 4: Integral Matrix (Pseudoinverse of Derivative)

Integral of $d_1 + d_2 x + d_3 x^2$:

$$(1 \;\; x \;\; x^2)\begin{pmatrix} d_1 \\ d_2 \\ d_3 \end{pmatrix}$$

is $d_1 x + \frac{1}{2}d_2 x^2 + \frac{1}{3}d_3 x^3$:

$$(1 \;\; x \;\; x^2 \;\; x^3)\begin{pmatrix} 0 \\ d_1 \\ \frac{1}{2}d_2 \\ \frac{1}{3}d_3 \end{pmatrix}$$

$$= (1 \;\; x \;\; x^2 \;\; x^3)\begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix}\begin{pmatrix} d_1 \\ d_2 \\ d_3 \end{pmatrix}$$

The **integral matrix** is:

$$D^+ = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix}$$

**Products $DD^+$ and $D^+D$:**

$$DD^+ = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}\begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = I$$

$$D^+D = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix}\begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix} = \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

The derivative $T$ has a **kernel** (constant function): $T(1) = 0$.

Its matrix $D$ has a **nullspace**.

---

### 2.8 Matrix Product AB Matches Transformation TS

Consider $TS$:

$$(TS)(\mathbf{u}) := T \circ S(\mathbf{u}) = T(S(\mathbf{u}))$$

$$(AB)(\mathbf{x}) := A(B\mathbf{x})$$

Matrix multiplication gives the correct matrix $AB$ to represent $TS$.

$$U \xrightarrow{S} V \xrightarrow{T} W$$

$$\dim U = p, \quad \dim V = n, \quad \dim W = m$$

$$B: n \times p, \quad A: m \times n$$

The linear transformation $TS$:
- Start with any vector $\mathbf{u}$ in $U$
- Goes to $S(\mathbf{u})$ in $V$
- Then to $T(S(\mathbf{u}))$ in $W$

The matrix $AB$:
- Starts with any $\mathbf{x} \in \mathbb{R}^p$
- Goes to $B\mathbf{x} \in \mathbb{R}^n$
- Then to $AB\mathbf{x} \in \mathbb{R}^m$

The matrix $AB$ correctly represents $TS$:

$$TS: U \to V \to W$$
$$AB: (m \times n)(n \times p) = (m \times p)$$

---

### 2.9 Example 5: Rotation Composition

$S$ rotates the plane by $\theta$.

$T$ also rotates by $\theta$.

Then $TS$ rotates by $2\theta$.

This transformation $T^2$ corresponds to the rotation matrix $A^2$ through $2\theta$:

$$T = S, \quad A = B, \quad T^2 = \text{rotation by } 2\theta$$

$$A^2 = \begin{pmatrix} \cos 2\theta & -\sin 2\theta \\ \sin 2\theta & \cos 2\theta \end{pmatrix}$$

**Verification by matrix multiplication:**

$$A^2 = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

$$= \begin{pmatrix} \cos^2\theta - \sin^2\theta & -2\sin\theta\cos\theta \\ 2\sin\theta\cos\theta & \cos^2\theta - \sin^2\theta \end{pmatrix}$$

Using double-angle formulas:
- $\cos^2\theta - \sin^2\theta = \cos 2\theta$
- $2\sin\theta\cos\theta = \sin 2\theta$

---

### 2.10 Example 6: Inverse Rotation

$S$ rotates by the angle $\theta$.

$T$ rotates by $-\theta$.

$$TS = T \circ S = I \longrightarrow AB = I$$

$$AB\mathbf{x} = I\mathbf{x} = \mathbf{x}$$

$$AB = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} \cos^2\theta + \sin^2\theta & 0 \\ 0 & \cos^2\theta + \sin^2\theta \end{pmatrix} = I$$

---

<br>

## 3. The Search for a Good Basis (8.3)

### 3.1 Key Ideas Summary

1. With a new input basis $B_{\text{in}}$ and output basis $B_{\text{out}}$, every matrix $A$ becomes $B_{\text{out}}^{-1} A B_{\text{in}}$.

2. $B_{\text{in}} = B_{\text{out}}$ = generalized eigenvectors of $A$ produces the **Jordan form** $J = B^{-1}AB$.

3. The **Fourier matrix** $F = B_{\text{in}} = B_{\text{out}}$ diagonalizes every **circulant matrix**.

4. Sines, cosines, $e^{ikx}$, **Legendre** and **Chebyshev**: those are great bases for function space.

---

### 3.2 The Change of Basis Formula

Input basis vectors form $B_{\text{in}} = (\mathbf{b}_1 \;\; \mathbf{b}_2 \;\; \cdots \;\; \mathbf{b}_n)$.

Output basis vectors form $B_{\text{out}} = (\mathbf{b}_1' \;\; \mathbf{b}_2' \;\; \cdots \;\; \mathbf{b}_m')$.

Always, $B_{\text{in}}$ and $B_{\text{out}}$ are **invertible**.

The new matrix representation is:

$$B_{\text{out}}^{-1} \; A \; B_{\text{in}}$$

where $B_{\text{out}}$ is $m \times m$, $A$ is $m \times n$, $B_{\text{in}}$ is $n \times n$.

**Example:** If $B_{\text{in}} = I_{n \times n}$ and $B_{\text{out}} = I_{m \times m}$, the matrix stays as $A$.

**When $B = B_{\text{in}} = B_{\text{out}}$:**

$$B^{-1}AB \quad \text{is similar to } A$$

Similar matrices have the **same eigenvalues**.

---

### 3.3 Best Basis 1: Eigenvectors (Diagonalization)

$B_{\text{in}} = B_{\text{out}}$ = eigenvector matrix $X$.

$$A = X\Lambda X^{-1} \quad \Longleftrightarrow \quad X^{-1}AX = \Lambda$$

$A$ is a square matrix with $n$ independent eigenvectors. $A$ must be **diagonalizable**.

$\Lambda$ is obtained when $B_{\text{in}} = B_{\text{out}} = X$ (eigenvectors).

---

### 3.4 Best Basis 2: Singular Vectors (SVD)

$B_{\text{in}} = V$, $B_{\text{out}} = U$: singular vectors of $A$.

$$A = U\Sigma V^T$$

$$U^{-1}AV = \Sigma \quad \text{(singular values)}$$

$U, V$ are **orthonormal eigenvectors** of $A^TA$ and $AA^T$.

---

### 3.5 Best Basis 3: Generalized Eigenvectors (Jordan Form)

$B_{\text{in}} = B_{\text{out}}$ = generalized eigenvectors of $A$.

$$B^{-1}AB = \text{Jordan form } J$$

$A$ is a square matrix, but it may only have "$s$ independent eigenvectors":

- If $s = n$: $B = X$ (eigenvectors), $J = \Lambda$
- If $s < n$: there are "$n - s$" dependent columns $\to$ construct "$n - s$" additional **generalized eigenvectors**

Properties of the Jordan form:
- (i) There are $s$ square blocks along the diagonal of $J$
- (ii) Each block has one eigenvalue $\lambda$, one eigenvector, and 1's above the diagonal

---

### 3.6 Jordan Form: Structure and Definition

When all eigenvectors are independent ($s = n$):

$$J = \Lambda = \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots & \\ & & & \lambda_n \end{pmatrix} \quad \text{(n 1x1 blocks)}$$

**Example 1:**

$$J = \begin{pmatrix} \boxed{2} & & & \\ & \boxed{2} & & \\ & & \boxed{3} & 1 \\ & & & \boxed{3} \end{pmatrix}$$

This has **2** 1x1 blocks and **1** 2x2 block.

$B^{-1}AB = J$ is **nearly diagonal**.

**Jordan Form — General Statement:**

For every $A$, we want to choose $B$ such that $B^{-1}AB$ is nearly diagonal as possible.

- When $A$ is diagonalizable (i.e., $n$ independent eigenvectors), then $B = X$ (eigenvectors): $X^{-1}AX = \Lambda$.
- When $A$ is not diagonalizable (i.e., $s < n$ independent eigenvectors):

$$B^{-1}AB = J$$

The Jordan form of $A$ is:

$$J = \begin{pmatrix} J_1 & & \\ & J_2 & \\ & & \ddots & \\ & & & J_s \end{pmatrix}$$

where each **Jordan block** $J_i$ is:

$$J_i = \begin{pmatrix} \lambda_i & 1 & 0 & \cdots & 0 \\ 0 & \lambda_i & 1 & \cdots & 0 \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & \cdots & 0 & \lambda_i & 1 \\ 0 & \cdots & 0 & 0 & \lambda_i \end{pmatrix}$$

> **The best basis $B$ gives $B^{-1}AB = J$.**

---

### 3.7 Example: 2x2 Jordan Form

$$A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, \quad \lambda^2 - 2\lambda + 1 = (\lambda - 1)^2 = 0$$

**(i)** $\lambda_1 = 1$

$$(A - \lambda I)\mathbf{x} = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$x_2 = 0$, free variable $x_1 \Rightarrow \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

**(ii)** Can we find the generalized eigenvector?

$$(A - \lambda_1 I)\mathbf{x}_2 = \mathbf{x}_1$$

(Note: $(A - \lambda_1 I)^2\mathbf{x}_2 = (A - \lambda_1 I)\mathbf{x}_1 = \mathbf{0}$)

$$(A - \lambda I)\mathbf{x}_2 = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$$

$$2x_2 = 1 \quad \Rightarrow \quad x_2 = \frac{1}{2}$$

Free variable $x_1 \to x_1 = 0$

$$\therefore \mathbf{x}_2 = \begin{pmatrix} 0 \\ \frac{1}{2} \end{pmatrix} \perp \mathbf{x}_1$$

**(iii)** Form the basis matrix:

$$B = (\mathbf{x}_1 \;\; \mathbf{x}_2) = \begin{pmatrix} 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix}, \quad B^{-1} = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$$

$$B^{-1}AB = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & \frac{1}{2} \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = J$$

This is **1 Jordan block** (2x2).

---

### 3.8 Example: 3x3 Jordan Form

$$A = \begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}, \quad \lambda^3 = 0 \;\Rightarrow\; \lambda = 0$$

**(i)** Find eigenvector: $(A - \lambda I)\mathbf{x} = \mathbf{0}$

$$(A - 0 \cdot I)\mathbf{x} = \begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$x_3 = 0$, $x_2 + 2x_3 = 0 \Rightarrow x_2 = 0$, free variable $x_1$

$$\Rightarrow \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$$

**(ii)** Find first generalized eigenvector: $(A - \lambda I)\mathbf{x}_2 = \mathbf{x}_1$

$$\begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$$

$x_3 = 0$, $x_2 + 2x_3 = 1 \Rightarrow x_2 = 1$, $x_1$ is free $\to x_1 = 0$

$$\Rightarrow \mathbf{x}_2 = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$$

**(iii)** Find second generalized eigenvector: $(A - \lambda I)\mathbf{x}_3 = \mathbf{x}_2$

$$\begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$$

$x_3 = 1$, $x_2 + 2x_3 = 0 \Rightarrow x_2 = -2$, choose $x_1 = 0$

$$\Rightarrow \mathbf{x}_3 = \begin{pmatrix} 0 \\ -2 \\ 1 \end{pmatrix}$$

**(iv)** Jordan chains:

$$(A - \lambda I)\mathbf{x}_1 = \mathbf{0}$$
$$(A - \lambda I)^2\mathbf{x}_2 = (A - \lambda I)\mathbf{x}_1 = \mathbf{0}$$
$$(A - \lambda I)^3\mathbf{x}_3 = (A - \lambda I)^2\mathbf{x}_2 = \mathbf{0}$$

Form the basis:

$$B = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -2 \\ 0 & 0 & 1 \end{pmatrix}, \quad B^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}$$

$$B^{-1}AB = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -2 \\ 0 & 0 & 1 \end{pmatrix}$$

$$= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix} = J$$

This is **1 Jordan block** (3x3) with eigenvalue $\lambda = 0$.

---

### 3.9 Basis for Function Space

Consider the standard polynomial basis: $1, x, x^2, x^3, \ldots$

$$x^{10} \approx \text{span}\{1, x, x^2, \ldots, x^9\}$$

This is an **ill-conditioned basis**.

How can we check if a basis is good or not?

$$B = \begin{pmatrix} \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_n \end{pmatrix} \quad \text{(basis vectors)}$$

$$B^TB = I \quad \text{when the basis is orthonormal. This is the best.}$$

**Inner product of two vectors:** $\mathbf{b}_i^T\mathbf{b}_j$

**Inner product of two functions:** $\displaystyle\int_0^1 x^i \cdot x^j \, dx$

$$\int_0^1 x^i \cdot x^j \, dx = \frac{1}{i + j + 1} \longrightarrow 0 \quad \text{as } i, j \to \infty$$

**Inner product of $f(x)$ and $g(x)$:**

$$(f, g) := \int f(x)\,g(x)\,dx$$

**Complex inner product:**

$$(f, g) := \int \overline{f(x)}\,g(x)\,dx$$

**Weighted inner product:**

$$(f, g)_w := \int w(x)\,\overline{f(x)}\,g(x)\,dx$$

---

### 3.10 Orthogonal Bases for Function Space

**Fourier basis:**

$$1, \;\sin x, \;\cos x, \;\sin 2x, \;\ldots$$

**Legendre basis:**

$$1, \;x, \;x^2 - \frac{1}{3}, \;x^3 - \frac{3}{5}x, \;\ldots$$

**Chebyshev basis:**

$$1, \;x, \;2x^2 - 1, \;4x^3 - 3x, \;\ldots$$

From the powers of $x$: $1, x, x^2, \ldots$, how can we construct the Legendre basis?

---

### 3.11 Constructing the Legendre Basis via Gram-Schmidt

Compute inner products on $[-1, 1]$:

$$(1, 1) = \int_{-1}^{1} 1\,dx = 2$$

$$(1, x) = \int_{-1}^{1} 1 \cdot x\,dx = 0 \quad \text{(orthogonal)}$$

$$(1, x^2) = \int_{-1}^{1} 1 \cdot x^2\,dx = 2\int_0^1 x^2\,dx = \frac{2}{3}$$

$x^2$ has a nonzero component of $1$.

$$(x, x^2) = \int_{-1}^{1} x^3\,dx = 0 \quad \text{(orthogonal)}$$

By **Gram-Schmidt process**, construct the orthogonal basis:

$$b_1 = 1$$
$$b_2 = x \quad \perp\; b_1$$
$$b_3 = x^2 - \frac{1}{3} \quad \perp\; b_1,\; b_2$$
$$\vdots$$

By using the **Gram-Schmidt process**, we can obtain the **Legendre basis**.

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Linear Transformation | $T(c\mathbf{u} + d\mathbf{w}) = cT(\mathbf{u}) + dT(\mathbf{w})$; must satisfy additivity and homogeneity |
| $T(\mathbf{0}) = \mathbf{0}$ | Every linear transformation maps zero to zero |
| Non-linear examples | Norm $\|\mathbf{u}\|$, shift $\mathbf{u} + \mathbf{u}_0$ are NOT linear |
| Affine transformation | $T(\mathbf{u}) = A\mathbf{u} + \mathbf{u}_0$ (linear + shift); used in computer graphics |
| Dot product | $T(\mathbf{u}) = \mathbf{a}^T\mathbf{u}$ is linear |
| Derivative | $T(u) = du/dt$ is a linear transformation on function space |
| Integration | $T^+(f) = \int_0^x f(t)\,dt$ is also linear; its matrix $A^+$ satisfies $AA^+ = I$ |
| Projection onto $z=1$ | NOT linear because $T(\mathbf{0}) \neq \mathbf{0}$ |
| Range and Kernel | Range of $T$ = column space of $A$; Kernel of $T$ = nullspace of $A$ |
| Determined by basis | A linear transformation is fully determined by its values on basis vectors |
| Matrix of $T$ | Column $j$ = $T(\mathbf{u}_j)$ expressed in output basis; $A = (T(\mathbf{u}_1) \;\cdots\; T(\mathbf{u}_n))$ |
| Standard matrix | Uses columns of $I$ as input/output bases: $T(\mathbf{u}) = A\mathbf{u}$ |
| Change of basis | New matrix = $B_{\text{out}}^{-1} A B_{\text{in}}$; identity transformation gives $B = W^{-1}V$ |
| Derivative matrix $D$ | For $\{1, x, x^2, x^3\} \to \{1, x, x^2\}$: $D = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$ |
| Integral matrix $D^+$ | Pseudoinverse of $D$; $DD^+ = I$ but $D^+D \neq I$ (nullspace of $D$) |
| Composition $TS$ | Matrix product $AB$ represents $TS$; $(m \times n)(n \times p) = (m \times p)$ |
| Rotation composition | Rotation by $\theta$ twice $= A^2 =$ rotation by $2\theta$ |
| Best basis: Eigenvectors | $B = X$ (eigenvectors) $\Rightarrow$ $X^{-1}AX = \Lambda$ (diagonal) |
| Best basis: Singular vectors | $B_{\text{in}} = V$, $B_{\text{out}} = U$ $\Rightarrow$ $U^{-1}AV = \Sigma$ |
| Best basis: Generalized eigenvectors | $B^{-1}AB = J$ (Jordan form); used when $A$ has fewer than $n$ independent eigenvectors |
| Jordan block $J_i$ | $\lambda_i$ on diagonal, 1's on superdiagonal; $s$ blocks for $s$ independent eigenvectors |
| Jordan chains | $\mathbf{x}_1$: eigenvector; $(A - \lambda I)\mathbf{x}_2 = \mathbf{x}_1$; $(A - \lambda I)\mathbf{x}_3 = \mathbf{x}_2$; etc. |
| Fourier matrix | Diagonalizes every circulant matrix |
| Function space bases | Fourier ($1, \sin x, \cos x, \ldots$), Legendre ($1, x, x^2 - 1/3, \ldots$), Chebyshev ($1, x, 2x^2-1, \ldots$) |
| Gram-Schmidt for functions | Inner product $(f,g) = \int f\,g\,dx$; produces Legendre polynomials from $1, x, x^2, \ldots$ |
| Orthonormal basis | $B^TB = I$; best conditioned basis for computation |

---
