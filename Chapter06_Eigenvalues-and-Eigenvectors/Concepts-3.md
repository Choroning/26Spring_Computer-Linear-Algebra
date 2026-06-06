# Chapter 6 Lecture — Eigenvalues and Eigenvectors — Part 3/3

> **📑 This document is split into 3 parts.** GitHub renders only a limited number of math expressions per file, so it is split by section so that all math displays correctly.
>
> [Part 1](Concepts.md) · [Part 2](Concepts-2.md) · **Part 3**

---

<br>

## Table of Contents

- [4. Solving Linear Differential Equations (6.5)](#4-solving-linear-differential-equations-65)
  - [4.1 Key Facts](#41-key-facts)
  - [4.2 Scalar ODE Review](#42-scalar-ode-review)
  - [4.3 Solution of du/dt = Au](#43-solution-of-dudt--au)
  - [4.4 General n x n Solution Procedure](#44-general-n-x-n-solution-procedure)
  - [4.5 Exponential of a Matrix](#45-exponential-of-a-matrix)
  - [4.6 Second Order Equations](#46-second-order-equations)
  - [4.7 Stability of 2 by 2 Matrices](#47-stability-of-2-by-2-matrices)
  - [4.8 Worked Examples](#48-worked-examples)
- [Summary](#summary)

---

<br>

## 4. Solving Linear Differential Equations (6.5)

### 4.1 Key Facts

**(1)** If $`A\mathbf{x} = \lambda\mathbf{x}`$, then $`\mathbf{u}(t) = e^{\lambda t}\mathbf{x}`$ will solve $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$. Each $`\lambda`$ and $`\mathbf{x}`$ give a solution $`e^{\lambda t}\mathbf{x}`$.

**(2)** If $`A = X\Lambda X^{-1}`$, then

```math
\mathbf{u}(t) = e^{At}\mathbf{u}(0) = Xe^{\Lambda t}X^{-1}\mathbf{u}(0) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n
```

**(3)** Matrix exponential: $`e^{At} = I + At + \cdots + \frac{(At)^n}{n!} + \cdots = Xe^{\Lambda t}X^{-1}`$ if $`A = X\Lambda X^{-1}`$.

**(4)** $`A`$ is **stable** and $`\mathbf{u}(t) \to \mathbf{0}`$ and $`e^{At} \to 0`$ when all eigenvalues of $`A`$ have **real part $`< 0`$**.

**(5)** $`\frac{d^2\mathbf{u}}{dt^2} + B\frac{d\mathbf{u}}{dt} + C\mathbf{u} = 0`$ means $`\frac{d}{dt}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix} = \begin{pmatrix} 0 & I \\ -C & -B \end{pmatrix}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix}`$ where $`\mathbf{v} = \frac{d\mathbf{u}}{dt}`$.

### 4.2 Scalar ODE Review

Consider ordinary differential equations:

```math
\frac{du}{dt} = u \implies u(t) = Ce^t
```

```math
\boxed{\frac{du}{dt} = \lambda u} \implies \boxed{u(t) = Ce^{\lambda t}}
```

Check: $`\frac{du}{dt} = \lambda Ce^{\lambda t} = \lambda u`$.

At $`t = 0`$: $`u(0) = Ce^{\lambda \cdot 0} = C \cdot 1 = C`$. The initial value is $`C`$.

$`\therefore u(t) = u(0)e^{\lambda t} = u_0 e^{\lambda t}`$

**Behavior:**

- $`\lambda > 0`$: **unstable** (growing)
- $`\lambda = 0`$: **steady state** (constant)
- $`\lambda < 0`$: **stable** (decaying)

**What if $`\lambda \in \mathbb{C}`$?** e.g., $`\lambda = -1 + 2i`$

$`e^{\lambda t} = e^{-t + 2it} = e^{-t}e^{2it}`$

$`|e^{\lambda t}| = |e^{-t}||e^{2it}| = e^{-t} \cdot \underbrace{1}_{|e^{i\theta}| = 1}`$

**Observation:** Stability depends on the **real part** of $`\lambda`$.

- $`\text{Re}(\lambda) > 0`$: unstable
- $`\text{Re}(\lambda) = 0`$: steady (oscillation)
- $`\text{Re}(\lambda) < 0`$: stable

### 4.3 Solution of du/dt = Au

**Example 1:** Consider a coupled ODE:

```math
\frac{dy}{dt} = z, \quad \frac{dz}{dt} = y
```

```math
\implies \frac{d}{dt}\begin{pmatrix} y \\ z \end{pmatrix} = \begin{pmatrix} z \\ y \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y \\ z \end{pmatrix}
```

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u}
```

Given an initial condition $`\mathbf{u}_0 = \begin{pmatrix} y(0) \\ z(0) \end{pmatrix}`$, what is $`\mathbf{u}(t) = ?`$

Let $`\lambda`$ be an eigenvalue of $`A`$ and $`\mathbf{x}`$ be the corresponding eigenvector.

Choose $`\mathbf{u} = e^{\lambda t}\mathbf{x}`$ and plug into the coupled ODE:

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u} = e^{\lambda t}A\mathbf{x} = e^{\lambda t}\lambda\mathbf{x} = \lambda e^{\lambda t}\mathbf{x}
```

Since $`\mathbf{u} = e^{\lambda t}\mathbf{x}`$ satisfies the coupled ODE, $`e^{\lambda t}\mathbf{x}`$ is a solution to $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$.

```math
\implies \mathbf{u} = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2
```

What are the eigenvalues of $`A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}`$?

$`|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0`$, so $`\lambda = \pm 1`$.

i) $`\lambda_1 = 1`$: $`(A - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, so $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`\lambda_2 = -1`$: $`(A + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, so $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) The complete solution: $`\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2`$

iv) The constants $`c_1`$ and $`c_2`$ can be determined by the initial condition $`\mathbf{u}_0 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}`$:

```math
\mathbf{u}(t=0) = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}
```

Write $`\mathbf{u}_0`$ as a combination of the eigenvectors of $`A`$:

```math
= (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}
```

```math
\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = -\frac{1}{2}\begin{pmatrix} -1 & -1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}
```

### 4.4 General n x n Solution Procedure

Generalize the idea to $`n \times n`$ matrix $`A`$:

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u}
```

$`A`$ has $`n`$ eigenvalues and $`n`$ eigenvectors ($`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$).

1. Write $`\mathbf{u}_0`$ as $`\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_n\mathbf{x}_n`$.
2. Multiply eigenvector $`\mathbf{x}_i`$ by $`e^{\lambda_i t}`$.
3. The solution to $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$ is:

```math
\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}
```

**Example 2:**

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 1 \\ 0 & 0 & 3 \end{pmatrix}\mathbf{u}, \quad \mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}
```

$`\lambda_1 = 1, \lambda_2 = 2, \lambda_3 = 3`$ (upper triangular — eigenvalues on diagonal).

i) $`(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}`$

ii) $`(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}`$

iii) $`(A - 3I)\mathbf{x}_3 = \begin{pmatrix} -2 & 1 & 1 \\ 0 & -1 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_3 = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$

iv) $`\mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2 \; \mathbf{x}_3)\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix}`$

```math
\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}^{-1}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}
```

$`|X| = 1`$, $`X^{-1} = \frac{1}{|X|}C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}`$

```math
\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix}
```

v) $`\therefore \mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + c_3 e^{\lambda_3 t}\mathbf{x}_3`$

```math
= 2e^t\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + 3e^{2t}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} + 4e^{3t}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}
```

### 4.5 Exponential of a Matrix

Recall: $`e^x = 1 + x + \frac{1}{2}x^2 + \frac{1}{6}x^3 + \cdots + \frac{x^n}{n!} + \cdots = \sum_{i=0}^{\infty}\frac{x^i}{i!}`$

$`\frac{1}{1-x} = 1 + x + x^2 + x^3 + \cdots + x^n + \cdots = \sum_{i=0}^{\infty}x^i`$ (geometric series, $`|x| < 1`$)

**Matrix exponential:**

```math
e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots + \frac{1}{n!}(At)^n + \cdots
```

$`(I - At)^{-1} = I + At + (At)^2 + (At)^3 + \cdots + (At)^n + \cdots`$

Let $`\mathbf{u}(t) = e^{At}\mathbf{u}(0) = e^{At}\mathbf{u}_0`$.

Check if $`\mathbf{u}`$ is the solution to $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$:

```math
\frac{d\mathbf{u}}{dt} = \frac{d}{dt}(e^{At}\mathbf{u}_0) = \frac{de^{At}}{dt}\mathbf{u}_0 = \left(A + \frac{1}{2} \cdot 2 \cdot (At)A + \frac{3}{6}(At)^2 A + \cdots + \frac{n}{n!}(At)^{n-1}A + \cdots\right)\mathbf{u}_0
```

```math
= A\left(I + At + \frac{1}{2}(At)^2 + \cdots + \frac{1}{(n-1)!}(At)^{n-1} + \cdots\right)\mathbf{u}_0 = Ae^{At}\mathbf{u}_0 = A\mathbf{u}
```

**If $`A`$ is diagonalizable,** $`A = X\Lambda X^{-1}`$, then:

```math
e^{At} = Xe^{\Lambda t}X^{-1}
```

where

```math
e^{\Lambda t} = I + \Lambda t + \frac{1}{2}(\Lambda t)^2 + \frac{1}{6}(\Lambda t)^3 + \cdots = \begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}
```

**Proof:**

```math
e^{At} = I + At + \frac{1}{2}(At)^2 + \cdots
```

```math
= I + X\Lambda X^{-1}t + \frac{1}{2}t^2(X\Lambda X^{-1})(X\Lambda X^{-1}) + \frac{1}{6}t^3(X\Lambda X^{-1})(X\Lambda X^{-1})(X\Lambda X^{-1}) + \cdots
```

```math
= I + t(X\Lambda X^{-1}) + \frac{1}{2}t^2(X\Lambda^2 X^{-1}) + \frac{1}{6}t^3(X\Lambda^3 X^{-1}) + \cdots
```

```math
= (XIX^{-1}) + (Xt\Lambda X^{-1}) + (X\frac{t^2}{2}\Lambda^2 X^{-1}) + (X\frac{t^3}{6}\Lambda^3 X^{-1}) + \cdots
```

```math
= X\left(I + t\Lambda + \frac{t^2}{2}\Lambda^2 + \frac{t^3}{6}\Lambda^3 + \cdots\right)X^{-1}
```

```math
\therefore e^{At} = Xe^{\Lambda t}X^{-1}
```

**Observation:** $`A = X\Lambda X^{-1}`$ and $`e^{At} = Xe^{\Lambda t}X^{-1}`$ $`\implies`$ $`e^{At}`$ and $`A`$ share the **same eigenvectors**.

**The solution:**

```math
\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\mathbf{u}_0
```

Since $`\mathbf{u}_0 = X\mathbf{c}`$ (linear combination of eigenvectors):

```math
= Xe^{\Lambda t}\underbrace{X^{-1}X}_{I}\mathbf{c} = Xe^{\Lambda t}\mathbf{c} = X\begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{pmatrix} = X\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}
```

```math
= (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}
```

```math
\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}
```

**Property:** $`e^{A(s+t)} = e^{As} \cdot e^{At}`$

**Proof:** (Detailed computation using double series and binomial theorem)

```math
e^{As} \cdot e^{At} = \left(\sum_{j=0}^{\infty}\frac{(As)^j}{j!}\right)\left(\sum_{k=0}^{\infty}\frac{(At)^k}{k!}\right) = \sum_{j=0}^{\infty}\sum_{k=0}^{\infty}\frac{A^{j+k}s^j t^k}{j!\,k!}
```

Let $`n = j + k`$, $`k = n - j`$:

```math
= \sum_{n=0}^{\infty}\frac{A^n}{n!}\sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^j t^{n-j} = \sum_{n=0}^{\infty}\frac{A^n}{n!}(s+t)^n = \sum_{n=0}^{\infty}\frac{(A(s+t))^n}{n!} = e^{A(s+t)} \quad \square
```

Using the binomial theorem: $`(s + t)^n = \sum_{j=0}^{n}\binom{n}{j}s^{n-j}t^j = \sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^{n-j}t^j`$.

From $`e^{A(s+t)} = e^{As} \cdot e^{At}`$, take $`s = -t`$:

```math
e^{A(-t+t)} = e^0 = I = e^{-At} \cdot e^{At}
```

This implies that $`e^{-At}`$ is the **inverse** of $`e^{At}`$.

**Example 4:** Let $`A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}`$, what is $`e^{At} = ?`$

$`|A - \lambda I| = \lambda^2 + 1 = 0`$, so $`\lambda = \pm i`$ (antisymmetric matrix).

$`A^2 = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}`$

$`A^3 = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}`$

$`A^4 = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I`$

```math
e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots
```

```math
= \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}t + \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\frac{t^2}{2} + \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\frac{t^3}{6} + \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\frac{t^4}{24} + \cdots
```

```math
= \begin{pmatrix} 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots & t - \frac{t^3}{6} + \cdots \\ -t + \frac{t^3}{6} - \cdots & 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots \end{pmatrix}
```

```math
e^{At} = \begin{pmatrix} \cos(t) & \sin(t) \\ -\sin(t) & \cos(t) \end{pmatrix} \quad \text{(antisymmetric matrix gives rotation)}
```

**Example 5:** Solve $`\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}\mathbf{u}`$ with $`\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$ at $`t = 0`$.

Upper triangular matrix: $`\lambda_1 = 1, \lambda_2 = 2`$.

i) $`(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

ii) $`(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

iii) $`X = (\mathbf{x}_1 \; \mathbf{x}_2) = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \to X^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}`$

$`A = X\Lambda X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}`$

$`e^{At} = Xe^{\Lambda t}X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} e^t & \\ & e^{2t} \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}`$

iv) $`\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$

$`= Xe^{\Lambda t}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 2 \\ 1 \end{pmatrix} = Xe^{\Lambda t}\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

$`= e^t\begin{pmatrix} 1 \\ 0 \end{pmatrix} + e^{2t}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} e^t + e^{2t} \\ e^{2t} \end{pmatrix}`$

### 4.6 Second Order Equations

Consider $`my'' + by' + ky = 0`$ (mass-spring-damper system).

Take $`y = e^{\lambda t}`$: $`y' = \lambda e^{\lambda t} = \lambda y`$, $`y'' = (y')' = \lambda^2 y`$.

The spring equation becomes: $`m\lambda^2 y + b\lambda y + ky = (m\lambda^2 + b\lambda + k)e^{\lambda t} = 0`$

```math
\implies m\lambda^2 + b\lambda + k = 0
```

```math
\iff \lambda^2 + \frac{b}{m}\lambda + \frac{k}{m} = 0 \quad \text{---(*)}
```

$`\lambda_1 + \lambda_2 = -\frac{b}{m}`$, $`\lambda_1\lambda_2 = \frac{k}{m}`$

$`\lambda = \frac{1}{2}\left(-\frac{b}{m} \pm \sqrt{\frac{b^2}{m^2} - \frac{4k}{m}}\right)`$

Let $`y_1 = e^{\lambda_1 t}`$, $`y_2 = e^{\lambda_2 t}`$. Complete solution becomes:

```math
y(t) = c_1 y_1 + c_2 y_2 \quad \text{if } \lambda_1 \neq \lambda_2 \quad \text{---(**)}
```

**Convert 2nd order ODE to 1st order ODEs:**

$`my'' + by' + ky = 0 \iff y'' + \tilde{b}y' + \tilde{k}y = 0`$ where $`\tilde{b} = \frac{b}{m}`$, $`\tilde{k} = \frac{k}{m}`$.

Let $`\mathbf{y} = \begin{pmatrix} y' \\ y \end{pmatrix}`$:

$`y'' = -\tilde{b}y' - \tilde{k}y = (-\tilde{b} \; -\tilde{k})\begin{pmatrix} y' \\ y \end{pmatrix}`$

$`y' = y' = (1 \; 0)\begin{pmatrix} y' \\ y \end{pmatrix}`$

```math
\implies \frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -\tilde{b} & -\tilde{k} \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix} \implies \frac{d\mathbf{y}}{dt} = A\mathbf{y}
```

i) Eigenvalues of $`A`$:

$`|A - \lambda I| = \begin{vmatrix} -\tilde{b} - \lambda & -\tilde{k} \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + \tilde{b}\lambda + \tilde{k} = 0`$

The equation for the $`\lambda`$'s is the same as (*).

ii) Let $`\lambda_1, \lambda_2`$ be the eigenvalues of $`A`$:

$`(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -\tilde{b} - \lambda_1 & -\tilde{k} \\ 1 & -\lambda_1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 - \lambda_1 x_2 = 0`$: $`\mathbf{x}_1 = \begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix}`$

$`(A - \lambda_2 I)\mathbf{x}_2`$: $`x_1 - \lambda_2 x_2 = 0`$: $`\mathbf{x}_2 = \begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}`$

iii) $`\mathbf{y}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2`$

```math
\begin{pmatrix} y'(t) \\ y(t) \end{pmatrix} = c_1 e^{\lambda_1 t}\begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix} + c_2 e^{\lambda_2 t}\begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}
```

```math
\therefore y(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t} \quad \text{like (**)}
```

When $`m = 1, b = 0, k = 1`$: $`y'' + y = 0 \iff y'' = -y`$

```math
\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -y \\ y' \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}
```

$`|A - \lambda I| = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0`$, so $`\lambda = \pm i`$ — **oscillation**.

### 4.7 Stability of 2 by 2 Matrices

For the solution to $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$, there is a fundamental question:

Does the solution approach $`\mathbf{u} = \mathbf{0}`$ as $`t \to \infty`$?

Since the complete solution $`\mathbf{u}(t)`$ is built from pure solutions $`e^{\lambda t}\mathbf{x}`$, the stability depends on the eigenvalues of $`A`$.

$`\lambda = r + is`$

$`e^{\lambda t} = e^{rt}e^{ist}`$

$`|e^{\lambda t}| = |e^{rt}|`$

The **real part** of $`\lambda`$ controls the growth ($`r > 0`$) or the decay ($`r < 0`$).

For any $`2 \times 2`$ matrix $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$:

Negative $`\lambda_1`$ and $`\lambda_2`$ require:

i) $`\lambda_1 + \lambda_2 = a + d < 0`$ (trace negative)

ii) $`\lambda_1 \lambda_2 > 0`$ (determinant positive)

### 4.8 Worked Examples

**Example 6:** Solve $`y'' + 4y' + 3y = 0`$

Substitute $`y`$ with $`e^{\lambda t}`$: $`(\lambda^2 + 4\lambda + 3)e^{\lambda t} = 0`$

$`\implies \lambda^2 + 4\lambda + 3 = 0 \implies (\lambda + 3)(\lambda + 1) = 0`$

$`\therefore \lambda = -1, -3`$

$`\implies y(t) = c_1 e^{-t} + c_2 e^{-3t}`$ — decaying solution $`\to`$ **stable** solution.

Introduce $`\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}`$: $`y'' = -4y' - 3y`$, $`y' = y'`$

```math
\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -4 & -3 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad \frac{d\mathbf{u}}{dt} = A\mathbf{u}
```

$`|A - \lambda I| = \begin{vmatrix} -4 - \lambda & -3 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 4\lambda + 3 = 0`$, so $`\lambda = -1, -3`$.

Corresponding eigenvectors:

$`\begin{pmatrix} -4+1 & -3 \\ 1 & +1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}`$

$`\begin{pmatrix} -4+3 & -3 \\ 1 & +3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} -3 \\ 1 \end{pmatrix}`$

The complete solution is:

```math
\begin{pmatrix} y' \\ y \end{pmatrix} = \mathbf{u}(t) = c_1 e^{-t}\begin{pmatrix} -1 \\ 1 \end{pmatrix} + c_2 e^{-3t}\begin{pmatrix} -3 \\ 1 \end{pmatrix}
```

This leads to $`y(t) = c_1 e^{-t} + c_2 e^{-3t}`$.

**Example 7:** $`y'' - 2y' + y = 0`$

Substitute $`y`$ with $`e^{\lambda t}`$: $`(\lambda^2 - 2\lambda + 1) = (\lambda - 1)^2 = 0`$, so $`\lambda = 1, 1`$ (repeated root).

$`\implies y(t) = c_1 e^t + ?`$

Since the root is repeated, the second solution is $`c_2 \cdot t e^t`$.

```math
y(t) = c_1 e^t + c_2 t e^t
```

Introduce $`\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}`$: $`y'' = 2y' - y`$, $`y' = y'`$

```math
\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad A = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}
```

$`|A - \lambda I| = \begin{vmatrix} 2 - \lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 2\lambda + 1 = 0`$, so $`\lambda = 1, 1`$.

$`(A - I)\mathbf{x}_1 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

Since the number of eigenvectors is smaller than 2, diagonalization is **NOT possible**.

So, we compute $`e^{At}`$ from the definition of a series:

```math
e^{At} = I + (At) + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots
```

Write $`A = I + (A - I)`$:

```math
e^{At} = e^{It + (A-I)t} = e^{It} \cdot e^{(A-I)t}
```

Note that $`(A - I)^2 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}`$, so $`(A - I)^k = 0`$ for $`k \geq 2`$.

```math
e^{(A-I)t} = I + (A-I)t = I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t
```

```math
e^{At} = e^t\left(I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t\right) = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}
```

```math
\mathbf{u}(t) = e^{At}\mathbf{u}_0
```

```math
\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}\begin{pmatrix} y_0' \\ y_0 \end{pmatrix}
```

```math
\therefore y(t) = (e^t - te^t)y_0 + te^t y_0'
```

---

<br>

## Summary

| Concept | Key Idea |
|:--------|:---------|
| Eigenvalue equation | $`A\mathbf{x} = \lambda\mathbf{x}`$, $`\mathbf{x} \neq \mathbf{0}`$ |
| Characteristic polynomial | $`\det(A - \lambda I) = 0`$ gives $`n`$ eigenvalues |
| Trace and determinant | $`\text{trace}(A) = \sum \lambda_i`$, $`\det(A) = \prod \lambda_i`$ |
| Powers of eigenvalues | $`A^k\mathbf{x} = \lambda^k\mathbf{x}`$ |
| Inverse eigenvalue | $`A^{-1}\mathbf{x} = \lambda^{-1}\mathbf{x}`$ (if $`\lambda \neq 0`$) |
| Shift rule | $`(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}`$ |
| Symmetric matrix | Real eigenvalues, orthogonal eigenvectors |
| Skew-symmetric matrix | Pure imaginary or zero eigenvalues |
| Rotation matrix | $`\lambda = e^{\pm i\theta}`$ |
| Projection matrix | $`\lambda = 0`$ or $`1`$ only |
| Diagonalization | $`A = X\Lambda X^{-1}`$ when $`n`$ LI eigenvectors exist |
| Matrix powers via diag. | $`A^k = X\Lambda^k X^{-1}`$ |
| Distinct eigenvalues | Eigenvectors are LI $`\implies`$ diagonalizable |
| Similar matrices | $`C = B^{-1}AB`$ has same eigenvalues as $`A`$ |
| GM vs AM | GM $`\leq`$ AM; if GM $`<`$ AM, not diagonalizable |
| Spectral theorem | $`S = Q\Lambda Q^T`$ for symmetric $`S`$ |
| Positive definite | $`\mathbf{x}^T S\mathbf{x} > 0`$ for all $`\mathbf{x} \neq \mathbf{0}`$; all $`\lambda_i > 0`$; all pivots $`> 0`$; all leading determinants $`> 0`$ |
| Positive semidefinite | $`\mathbf{x}^T S\mathbf{x} \geq 0`$; allows $`\lambda = 0`$ |
| Energy test for $`A^T A`$ | $`\mathbf{x}^T A^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0`$ |
| PD $`\implies`$ invertible | $`\det(S) > 0`$, $`S^{-1}`$ is also PD |
| PD $`\implies S = A^T A`$ | With independent columns in $`A`$ |
| Congruent matrices | $`A^T SA`$ preserves definiteness type |
| Ellipse equation | $`\mathbf{x}^T S\mathbf{x} = 1`$; axes along eigenvectors, lengths $`1/\sqrt{\lambda}`$ |
| Hessian and convexity | $`H`$ positive definite $`\implies`$ $`f`$ is strictly convex |
| Gradient descent | $`\mathbf{x}_{k+1} = \mathbf{x}_k - \eta\nabla f(\mathbf{x}_k)`$ |
| ODE solution | $`\frac{d\mathbf{u}}{dt} = A\mathbf{u} \implies \mathbf{u}(t) = \sum c_i e^{\lambda_i t}\mathbf{x}_i`$ |
| Matrix exponential | $`e^{At} = I + At + \frac{(At)^2}{2!} + \cdots = Xe^{\Lambda t}X^{-1}`$ |
| Stability (continuous) | Stable iff all $`\text{Re}(\lambda_i) < 0`$ |
| Stability (discrete) | $`A^k \to 0`$ iff all $`|\lambda_i| < 1`$ |
| $`e^{A(s+t)} = e^{As}e^{At}`$ | Matrix exponential additive property |
| $`(e^{At})^{-1} = e^{-At}`$ | Inverse of matrix exponential |
| 2nd order ODE | $`my'' + by' + ky = 0 \iff`$ 1st order system $`\frac{d\mathbf{y}}{dt} = A\mathbf{y}`$ |
| Repeated eigenvalue ODE | Non-diagonalizable case: use $`e^{At} = e^{It} \cdot e^{(A-I)t}`$ |

---
