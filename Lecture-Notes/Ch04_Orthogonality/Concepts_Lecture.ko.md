# 제4장 강의 — 직교성

> **최종 수정일:** 2026-06-19
>
> Introduction to Linear Algebra, Strang (6th Ed.) - Ch 4

> **선수 지식**: [선형대수학] 벡터 공간, 기저, 차원 (제1-3장).
>
> **학습 목표**:
> 1. 벡터와 부분공간의 직교성을 정의할 수 있다
> 2. 직선과 부분공간으로의 사영을 계산할 수 있다
> 3. 그람-슈미트 과정을 적용하여 직교 기저를 구성할 수 있다
> 4. 사영을 이용하여 최소제곱 문제를 풀 수 있다

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.**
>
> **Part 1** · [Part 2](Concepts-Part2.ko.md)

---

<br>

## 목차

- [1. 벡터와 부분공간의 직교성 (4.1)](#1-벡터와-부분공간의-직교성-41)
- [2. 직선과 부분공간으로의 사영 (4.2)](#2-직선과-부분공간으로의-사영-42)
- [3. 최소제곱 근사 (4.3)](Concepts-Part2.ko.md#3-최소제곱-근사-43)
- [4. 직교 기저와 그람-슈미트 (4.4)](Concepts-Part2.ko.md#4-직교-기저와-그람-슈미트-44)
- [5. 행렬의 유사역행렬 (4.5)](Concepts-Part2.ko.md#5-행렬의-유사역행렬-45)
- [요약](Concepts-Part2.ko.md#요약)

---

<br>

## 1. 벡터와 부분공간의 직교성 (4.1)

### 1.1 직교 벡터 (Orthogonal Vectors)

직교 벡터는 다음을 만족한다:

```math
\mathbf{v}^T \mathbf{w} = 0
```

그리고 벡터 형태의 **피타고라스 정리** (Pythagorean theorem):

```math
\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2
```

**증명:**

```math
\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w})
```

```math
= \mathbf{v} \cdot \mathbf{v} + \mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{v} + \mathbf{w} \cdot \mathbf{w}
```

```math
= \|\mathbf{v}\|^2 + 2\mathbf{v} \cdot \mathbf{w} + \|\mathbf{w}\|^2
```

$`\mathbf{v}^T \mathbf{w} = 0`$일 때 교차항 $`2\mathbf{v} \cdot \mathbf{w} = 0`$이 되어 다음이 성립한다:

```math
\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2
```

이는 $`a^2 + b^2 = c^2`$와 동일한 관계이다.

**제1장 복습:** 내적(dot product)은 $`\mathbf{v}`$와 $`\mathbf{w}`$ 사이의 각도와 연결된다:

```math
\mathbf{v}^T \mathbf{w} = \|\mathbf{v}\| \, \|\mathbf{w}\| \cos\theta
```

$`\theta = 90°`$일 때, $`\mathbf{v}^T \mathbf{w} = 0`$이다.

벡터 $`\mathbf{v}`$, $`\mathbf{w}`$, $`\mathbf{v} + \mathbf{w}`$는 직각삼각형을 이룬다. 피타고라스 정리에 의해:

```math
\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2
```

### 1.2 기본 부분공간은 직교한다 (The Fundamental Subspaces are Orthogonal)

**(1)** $`A`$의 영공간(nullspace) $`\mathcal{N}(A)`$는 행 공간(row space) $`C(A^T)`$에 직교하는 모든 벡터를 포함한다.

```math
A\mathbf{x} = \mathbf{0}
```

```math
\begin{pmatrix} \text{--- } A\text{의 행}_1 \text{ ---} \\ \text{--- } A\text{의 행}_2 \text{ ---} \\ \vdots \\ \text{--- } A\text{의 행}_m \text{ ---} \end{pmatrix} \begin{pmatrix} \\ \mathbf{x} \\ {} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ \vdots \\ 0 \end{pmatrix}
```

$`\mathbf{x}`$는 $`A`$의 각 행에 직교한다. 모든 행은 $`\mathbf{x}`$와의 내적이 0이다.

**(2)** $`A^T`$의 영공간 $`\mathcal{N}(A^T)`$는 열 공간(column space) $`C(A)`$에 직교하는 모든 벡터를 포함한다.

```math
A^T \mathbf{y} = \mathbf{0} \quad \Longleftrightarrow \quad \mathbf{y}^T A = \mathbf{0}
```

```math
\mathbf{y}^T \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} = \mathbf{0}
```

$`\Rightarrow \mathbf{y}^T \mathbf{a}_1 = 0, \quad \mathbf{y}^T \mathbf{a}_2 = 0, \quad \ldots, \quad \mathbf{y}^T \mathbf{a}_n = 0`$

$`\Rightarrow \mathbf{y}^T (c_1 \mathbf{a}_1 + c_2 \mathbf{a}_2 + \cdots + c_n \mathbf{a}_n) = 0`$

$`\mathbf{y}`$는 $`A`$의 각 열에 직교한다.

### 1.3 사영 미리보기 (Projection Preview)

$`\mathbf{b}`$가 $`A`$의 열 공간 밖에 있을 때, 열 공간 안에서 가장 가까운 점 $`\mathbf{p}`$를 찾는다.

```math
A\mathbf{x} \neq \mathbf{b}
```

오차(error) $`\mathbf{e}`$는 $`C(A)`$에 수직이다.

최소제곱 방정식 (least squares equation):

```math
A^T A \hat{\mathbf{x}} = A^T \mathbf{b}
```

은 가장 가까운 $`\mathbf{p} = A\hat{\mathbf{x}}`$를 만든다. $`A\mathbf{x} = \mathbf{b}`$가 풀 수 없을 때 **최적** 해 $`\hat{\mathbf{x}}`$를 제공한다:

```math
\min \|A\hat{\mathbf{x}} - \mathbf{b}\|^2
```

이것이 **최소제곱** (least squares)이다. $`A^T A = I`$이면 문제가 쉬워진다. 그러한 행렬 $`A`$를 구성하는 것이 4.4절의 주제이다.

### 1.4 직교 부분공간 (Orthogonal Subspaces)

**정의:** 부분공간 $`V`$와 $`W`$가 다음을 만족할 때 **직교** 한다:

```math
\mathbf{v}^T \mathbf{w} = 0 \quad \forall \, \mathbf{v} \in V, \, \mathbf{w} \in W
```

직선(1차원 부분공간)이 수직 평면(2차원 부분공간)을 통과하는 것을 생각하자. 직선 위의 모든 벡터는 평면 위의 모든 벡터에 수직이다.

**(3)** $`A`$의 행 공간은 $`A`$의 영공간에 직교한다. $`A`$의 열 공간은 $`A^T`$의 영공간(= $`A`$의 좌영공간)에 직교한다.

**$`C(A^T) \perp \mathcal{N}(A)`$ 증명:**

$`A`$의 영공간에 있는 $`\mathbf{x}`$는 $`A`$의 행 공간에 있는 $`A^T\mathbf{y}`$에 직교한다:

```math
\mathbf{x} \cdot (A^T \mathbf{y}) = (A^T \mathbf{y}) \cdot \mathbf{x} = (A^T \mathbf{y})^T \mathbf{x} = \mathbf{y}^T A \mathbf{x} = \mathbf{y}^T (A\mathbf{x}) = \mathbf{y}^T \mathbf{0} = 0
```

또한, $`A\mathbf{x} = \mathbf{0}`$으로부터:

```math
(\text{행}_1) \cdot \mathbf{x} = 0, \quad (\text{행}_2) \cdot \mathbf{x} = 0, \quad \ldots, \quad (\text{행}_m) \cdot \mathbf{x} = 0
```

```math
\Rightarrow [c_1(\text{행}_1) + c_2(\text{행}_2) + \cdots + c_m(\text{행}_m)] \cdot \mathbf{x} = 0
```

따라서 $`C(A^T)`$는 $`\mathcal{N}(A)`$에 수직이다.

**$`C(A) \perp \mathcal{N}(A^T)`$ 증명:**

마찬가지로, $`A^T`$의 영공간에 있는 $`\mathbf{y}`$는 $`A`$의 열 공간에 있는 $`A\mathbf{x}`$에 직교한다:

```math
\mathbf{y} \cdot (A\mathbf{x}) = (A\mathbf{x}) \cdot \mathbf{y} = (A\mathbf{x})^T \mathbf{y} = \mathbf{x}^T A^T \mathbf{y} = \mathbf{x}^T (A^T \mathbf{y}) = \mathbf{x}^T \mathbf{0} = 0
```

### 1.5 직교 여공간과 차원 (Orthogonal Complements and Dimension)

**(4)** 차원의 합:

```math
r + (n - r) = n \quad \text{그리고} \quad r + (m - r) = m
```

이것들이 **직교 여공간** (orthogonal complements)이다.

**예제:**

```math
A = \begin{pmatrix} 1 & -2 & 1 \\ 1 & 0 & -1 \end{pmatrix}
```

행 축소: $`R_2 - R_1`$, 그 다음 $`R_2/2`$, 그 다음 $`R_1 + 2R_2`$:

```math
R = \begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & -1 \end{pmatrix}
```

$`\text{rank}(A) = \text{rank}(R) = 2`$, $`n - r = 3 - 2 = 1`$ (자유 변수 1개).

```math
\mathcal{N}(A) = \left\lbrace c\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} \right\rbrace, \quad C(A) = \mathbb{R}^2
```

$`A^T`$에 대해:

```math
A^T = \begin{pmatrix} 1 & 1 \\ -2 & 0 \\ 1 & -1 \end{pmatrix}
```

행 축소: $`R_2 + 2R_1`$, $`R_3 - R_1`$, 그 다음 $`R_3 + R_2`$, $`R_2/2`$, 그 다음 $`R_1 - R_2`$:

```math
R_0 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}
```

$`\text{rank}(A^T) = \text{rank}(R_0) = 2`$, $`m - r = 2 - 2 = 0`$ (자유 변수 없음).

```math
\mathcal{N}(A^T) = \left\lbrace \begin{pmatrix} 0 \\ 0 \end{pmatrix} \right\rbrace, \quad C(A^T) \text{는 } \mathbb{R}^3\text{의 평면}
```

검증:

```math
\dim C(A^T) + \dim \mathcal{N}(A) = 2 + 1 = 3 = n
```

```math
\dim C(A) + \dim \mathcal{N}(A^T) = 2 + 0 = 2 = m
```

**중요한 제약:** $`V`$와 $`W`$가 $`\mathbb{R}^n`$의 직교 부분공간이면:

```math
\dim V + \dim W \leq n
```

전체 공간을 설명하는 두 직교 부분공간은 특별한 이름을 갖는다: **"직교 여공간"** (orthogonal complements). $`V`$의 직교 여공간 $`V^\perp`$는 $`V`$에 직교하는 모든 벡터를 포함한다.

**직교 여공간 쌍:**

| 쌍 | 차원 | 공간 |
|:-----|:-----------|:------|
| 행 공간과 영공간 | $`r + (n-r) = n`$ | $`\mathbb{R}^n`$ |
| 열 공간과 좌영공간 | $`r + (m-r) = m`$ | $`\mathbb{R}^m`$ |

```math
\mathcal{N}(A) \text{는 } \mathbb{R}^n\text{에서 행 공간 } C(A^T)\text{의 직교 여공간이다}
```

```math
\mathcal{N}(A^T) \text{는 } \mathbb{R}^m\text{에서 열 공간 } C(A)\text{의 직교 여공간이다}
```

모든 $`\mathbf{x}`$는 행 공간 성분 $`\mathbf{x}_r`$과 영공간 성분 $`\mathbf{x}_n`$으로 분리할 수 있다:

```math
\mathbb{R}^n \ni \mathbf{x} = \mathbf{x}_r + \mathbf{x}_n
```

```math
\mathbb{R}^m \ni \mathbf{y} = \mathbf{y}_{\text{col}} + \mathbf{y}_{\text{left null}}
```

### 1.6 네 부분공간의 큰 그림 (The Big Picture of Four Subspaces)

$`A`$를 $`\mathbf{x}`$에 곱하면:

```math
A\mathbf{x} = A(\mathbf{x}_r + \mathbf{x}_n) = A\mathbf{x}_r + A\mathbf{x}_n = \mathbf{b}
```

$`A\mathbf{x}_n = \mathbf{0}`$이므로.

- $`A\mathbf{x}_r = \mathbf{b}`$는 $`A`$의 열 공간에 있다.
- $`A\mathbf{x}_n = \mathbf{0}`$.

$`A\mathbf{x} = \mathbf{b}`$의 완전해(complete solution)는:

```math
\mathbf{x} = \mathbf{x}_r + \mathbf{x}_n
```

여기서 $`\mathbf{x}_r`$은 **유일한** 행 공간 성분이다.

### 1.7 최소 노름 해 (Minimum Norm Solution)

```math
\|\mathbf{x}\|^2 = \mathbf{x} \cdot \mathbf{x} = (\mathbf{x}_r + \mathbf{x}_n) \cdot (\mathbf{x}_r + \mathbf{x}_n)
```

```math
= \mathbf{x}_r \cdot \mathbf{x}_r + \mathbf{x}_r \cdot \mathbf{x}_n + \mathbf{x}_n \cdot \mathbf{x}_r + \mathbf{x}_n \cdot \mathbf{x}_n
```

```math
= \|\mathbf{x}_r\|^2 + \|\mathbf{x}_n\|^2
```

($`\mathbf{x}_r \perp \mathbf{x}_n`$이므로 교차항이 0이다.)

따라서 $`A\mathbf{x} = \mathbf{b}`$의 **최소 노름 해** 는 $`\mathbf{x} = \mathbf{x}_r`$이고 $`\mathbf{x}_n = \mathbf{0}`$이다.

**유일성:** 모든 벡터 $`\mathbf{b} \in C(A)`$는 행 공간의 정확히 **하나의** 벡터 $`\mathbf{x}_r`$에서 온다.

**증명:** $`A\mathbf{x}_r = A\mathbf{x}_r' = \mathbf{b}`$라 하자. 그러면 $`A\mathbf{x}_r - A\mathbf{x}_r' = A(\mathbf{x}_r - \mathbf{x}_r') = \mathbf{0}`$. 따라서 $`\mathbf{x}_r - \mathbf{x}_r' \in \mathcal{N}(A)`$. $`\mathbf{x}_r`$과 $`\mathbf{x}_r'`$ 모두 $`C(A^T)`$에서 오므로, $`\mathbf{x}_r - \mathbf{x}_r' \in C(A^T)`$. $`\mathbf{0} \in C(A^T)`$이고 $`\mathbf{0} \in \mathcal{N}(A)`$이므로, $`\mathbf{x}_r - \mathbf{x}_r' = \mathbf{0}`$, 즉, $`\mathbf{x}_r = \mathbf{x}_r'`$. $`\square`$

### 1.8 예제: 가역 부분행렬 (Invertible Submatrix)

**예제 2:** 랭크 $`r`$인 모든 행렬은 $`r \times r`$ 가역 부분행렬을 갖는다.

```math
A = \begin{pmatrix} 1 & 2 & 3 & 4 & 5 \\ 1 & 2 & 4 & 5 & 6 \\ 1 & 2 & 4 & 5 & 6 \end{pmatrix}
```

행 축소: $`R_2 - R_1`$, $`R_3 - R_1`$:

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

$`A`$는 피벗 열 1과 3에서 $`2 \times 2`$ 가역 부분행렬 $`\begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}`$를 포함한다.

### 1.9 부분공간의 기저 결합 (Combining Bases from Subspaces)

**기저** (basis)는 공간을 생성하는 선형독립 벡터들을 포함한다.

$`\mathbb{R}^n`$의 **표준 기저** (standard basis)는 $`\lbrace\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\rbrace`$이며:

```math
\mathbb{R}^n \ni \mathbf{e}_i = \begin{pmatrix} 0 \\ \vdots \\ 0 \\ 1 \\ 0 \\ \vdots \\ 0 \end{pmatrix} \leftarrow i\text{번째 행}
```

이는 $`i`$번째 성분이 1이고 나머지가 0인 $`\mathbb{R}^n`$의 벡터이다. 즉 $`I \in \mathbb{R}^{n \times n}`$의 $`i`$번째 열이다.

$`\mathbb{R}^n`$의 차원은 $`n`$이다. 이는 $`\mathbb{R}^n`$의 기저 벡터 수가 $`n`$개이기 때문이다, 예를 들어 $`\lbrace\mathbf{e}_1, \mathbf{e}_2, \ldots, \mathbf{e}_n\rbrace`$.

### 1.10 $`\mathbb{R}^n`$에서의 두 성질 (Two Properties in $`\mathbb{R}^n`$)

**i)** $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$이 선형독립(LI)이라고 가정하자. 그러면 $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$은 $`\mathbb{R}^n`$의 기저이다.

**ii)** $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m\rbrace`$이 $`\mathbb{R}^n`$을 생성한다고 가정하자. 그러면 $`m \geq n`$. 만약 $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$이 $`\mathbb{R}^n`$을 생성하면, $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$은 LI이다.

```math
\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace \text{이 LI} \iff \lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace \text{이 } \mathbb{R}^n\text{을 생성}
```

따라서 $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$은 $`\mathbb{R}^n`$의 기저이다.

**증명:**

**i)** $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$은 LI이다. 이 집합이 $`\mathbb{R}^n`$을 생성함을 보여야 한다.

$`A = (\mathbf{u}_1 \; \mathbf{u}_2 \; \cdots \; \mathbf{u}_n) \in \mathbb{R}^{n \times n}`$으로 정의하자. 목표는 $`A\mathbf{x} = \mathbf{v}`$를 만족하는 유일한 $`\mathbf{x}`$를 찾는 것이다. 정사각 행렬 $`A`$가 풀 랭크이므로, 역행렬이 존재한다: $`A^{-1}A\mathbf{x} = A^{-1}\mathbf{v}`$, 따라서 $`\mathbf{x} = A^{-1}\mathbf{v}`$. 그러므로 $`\mathbf{v} = A\mathbf{x}`$는 $`\mathbf{u}_i`$의 선형결합이다.

**ii)** $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m\rbrace`$이 $`\mathbb{R}^n`$을 생성하므로: $`\mathbf{v} \in \mathbb{R}^n = c_1\mathbf{u}_1 + c_2\mathbf{u}_2 + \cdots + c_m\mathbf{u}_m`$. $`A = (\mathbf{u}_1 \; \mathbf{u}_2 \; \cdots \; \mathbf{u}_m)`$이라 하면 이는 $`n \times m`$ 행렬이다. $`A`$의 기약 행 사다리꼴은 $`r \leq m`$개의 독립 열 벡터를 보여준다. $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_r\rbrace`$이 $`\mathbb{R}^n`$의 기저가 되지만, 이는 $`\dim \mathbb{R}^n = n`$에 모순이다. 따라서 $`m \geq n`$.

마지막으로 $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$이 LI가 아니라고 가정하자. 그러면 $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_r\rbrace`$이 $`r < n`$으로 $`\mathbb{R}^n`$의 기저가 되는데, 이는 $`\dim \mathbb{R}^n = n`$에 모순이다. 따라서 $`\lbrace\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n\rbrace`$은 LI이다. $`\square`$

### 1.11 $`A \in \mathbb{R}^{n \times n}`$의 성질 (Properties of $`A \in \mathbb{R}^{n \times n}`$)

**i)** $`A`$의 $`n`$개 열이 LI이면, $`\mathbb{R}^n`$을 생성한다.

**ii)** $`n`$개 열이 $`\mathbb{R}^n`$을 생성하면, LI이다 (위 증명 참고).

풀 랭크 정사각 행렬 $`A \in \mathbb{R}^{n \times n}`$에 대해:
- $`A\mathbf{x} = \mathbf{b}`$는 풀 수 있다 (존재성)
- $`\mathbf{x} = A^{-1}\mathbf{b}`$는 유일하다 (유일성)

**iii)** 풀 랭크 $`A, B \in \mathbb{R}^{n \times n}`$에 대해 $`AB = I`$이면, $`BA = I`$이다.

**증명:** $`AB = I`$. 그러면 $`B(AB) = BI = B`$ (결합법칙). $`(BA)B = IB`$ (분배법칙). $`(BA - I)B = 0`$, 따라서 $`BA = I`$.

### 1.12 예제 3: 네 부분공간 (Four Subspaces)

```math
A = \begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}
```

$`R_2 - 3R_1`$:

```math
\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0, \quad \text{rank}(A) = \text{rank}(R_0) = 1 = r
```

$`n = 2`$, $`n - r = 2 - 1 = 1`$ (자유 변수 1개).

```math
C(A^T) = \text{span}\lbrace(1, 2)\rbrace, \quad C(A) = \text{span}\left\lbrace\begin{pmatrix} 1 \\ 3 \end{pmatrix}\right\rbrace
```

$`A\mathbf{x} = \mathbf{0}`$: $`x_2 = 1`$로 선택, $`x_1 + 2 = 0`$, $`x_1 = -2`$.

```math
\mathcal{N}(A) = \text{span}\left\lbrace\begin{pmatrix} -2 \\ 1 \end{pmatrix}\right\rbrace
```

$`A^T\mathbf{y} = \mathbf{0}`$: $`\begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$. 축소 후: $`\begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}`$. $`y_2 = 1`$로 취하면, $`y_1 = -3`$.

```math
\mathcal{N}(A^T) = \text{span}\left\lbrace\begin{pmatrix} -3 \\ 1 \end{pmatrix}\right\rbrace
```

$`\mathbf{b} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}`$라 하자. $`A\mathbf{x} = \mathbf{b}`$의 해 $`\mathbf{x}`$는:

```math
\mathbf{x} = \mathbf{x}_p + c\,\mathbf{x}_n
```

```math
\begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}
```

$`x_2 = 0`$으로 취하면, $`x_1 = 10`$. 따라서 $`\mathbf{x}_p = \begin{pmatrix} 10 \\ 0 \end{pmatrix}`$.

```math
\mathbf{x} = \begin{pmatrix} 10 \\ 0 \end{pmatrix} + c\begin{pmatrix} -2 \\ 1 \end{pmatrix}
```

$`c = 1`$일 때: $`\begin{pmatrix} 1 & 2 \\ 3 & 6 \end{pmatrix}\begin{pmatrix} 8 \\ 1 \end{pmatrix} = \begin{pmatrix} 10 \\ 30 \end{pmatrix}`$. 검증 완료.

**참고:** $`\mathbf{x}_p`$는 $`\mathbf{x}_n`$에 직교하지 **않는다**: $`(10 \quad 0)\begin{pmatrix} -2 \\ 1 \end{pmatrix} = -20 \neq 0`$.

$`\mathbf{x}_p`$를 $`\mathbf{x}_r`$과 $`\mathbf{x}_n`$으로 더 분해할 수 있다:

```math
\begin{pmatrix} 10 \\ 0 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix} - 4\begin{pmatrix} -2 \\ 1 \end{pmatrix}
```

여기서 $`\mathbf{x}_r = \begin{pmatrix} 2 \\ 4 \end{pmatrix}`$이고 영공간 성분은 $`-4\begin{pmatrix} -2 \\ 1 \end{pmatrix}`$이다.

그러면 $`\mathbf{x} = \mathbf{x}_r + (c - 4)\mathbf{x}_n = \mathbf{x}_r + c'\mathbf{x}_n`$.

검증: $`(2 \quad 4)\begin{pmatrix} -2 \\ 1 \end{pmatrix} = -4 + 4 = 0`$이므로 $`\mathbf{x}_r \perp \mathbf{x}_n`$.

---

<br>

## 2. 직선과 부분공간으로의 사영 (4.2)

### 2.1 핵심 사실 요약 (Key Facts Summary)

**(1)** $`\mathbf{a}`$를 지나는 직선 위로의 $`\mathbf{b}`$의 사영(projection)은 $`\mathbf{b}`$에 가장 가까운 점이다:

```math
\mathbf{p} = \mathbf{a}\frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}
```

**(2)** 오차 $`\mathbf{e} = \mathbf{b} - \mathbf{p}`$는 $`\mathbf{a}`$에 수직이다. 직각삼각형 $`\mathbf{b}`$, $`\mathbf{p}`$, $`\mathbf{e}`$에 대해:

```math
\|\mathbf{p}\|^2 + \|\mathbf{e}\|^2 = \|\mathbf{b}\|^2
```

**(3)** 부분공간 $`S`$로의 $`\mathbf{b}`$의 사영은 $`S`$에서 가장 가까운 벡터 $`\mathbf{p}`$이다; $`\mathbf{b} - \mathbf{p}`$는 $`S`$에 직교한다.

**(4)** $`A`$가 독립인 열을 가질 때 $`A^T A`$는 가역(invertible)이고 대칭이다: $`\mathcal{N}(A^T A) = \mathcal{N}(A)`$.

**(5)** 열 공간 $`C(A)`$로의 $`\mathbf{b}`$의 사영:

```math
\mathbf{p} = A(A^T A)^{-1} A^T \mathbf{b}
```

**(6)** $`C(A)`$로의 사영 행렬:

```math
P = A(A^T A)^{-1} A^T
```

$`\mathbf{p} = P\mathbf{b}`$이고 $`P^2 = P = P^T`$이다.

### 2.2 직선으로의 사영 (Projection onto a Line)

직선이 $`\mathbf{a}`$ 방향으로 원점을 지난다. $`\mathbf{b}`$를 직선 위로 사영한다. $`\mathbf{b}`$에서 $`\mathbf{p}`$까지의 직선은 벡터 $`\mathbf{a}`$에 수직이다:

```math
\mathbf{e} = \mathbf{b} - \mathbf{p} \perp \mathbf{a}
```

사영 $`\mathbf{p}`$는 $`\mathbf{a}`$의 배수이다: $`\mathbf{p} = \alpha \mathbf{a}`$.

**유도:**

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

**특수한 경우:**
- $`\mathbf{b} = \mathbf{a}`$이면, $`\alpha = 1`$, $`\mathbf{p} = \mathbf{a}`$.
- $`\mathbf{b} \perp \mathbf{a}`$ (즉, $`\mathbf{a} \cdot \mathbf{b} = 0`$)이면, $`\alpha = 0`$, $`\mathbf{p} = \mathbf{0}`$.

**예제 1:** $`\mathbf{b} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$을 $`\mathbf{a} = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}`$ 위로 사영하여 $`\mathbf{p} = \alpha\mathbf{a}`$를 구한다.

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

$`\mathbf{b}`$는 두 부분으로 분리된다: $`\mathbf{b} = \mathbf{p} + \mathbf{e}`$, 여기서 $`\mathbf{p} \perp \mathbf{e}`$.

```math
\|\mathbf{p}\| = \|\mathbf{b}\|\cos\theta, \quad \|\mathbf{e}\| = \|\mathbf{b}\|\sin\theta
```

### 2.3 사영 행렬 P (직선으로의 사영) (Projection Matrix P onto a Line)

```math
\mathbf{p} = \alpha\mathbf{a} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\,\mathbf{a} = \mathbf{a}\frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}} = \left(\frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}\right)\mathbf{b}
```

```math
P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}
```

참고: $`\mathbf{a}\mathbf{a}^T`$는 **랭크 1** 행렬이다 (열 곱하기 행). $`\mathbf{a}`$를 지나는 직선인 1차원 부분공간으로 사영하며, 이는 $`C(P)`$이다.

**예제:** $`\mathbf{a} = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}`$를 지나는 직선으로의 사영 행렬 $`P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}`$를 구하라.

```math
\mathbf{a}^T\mathbf{a} = 9
```

```math
\mathbf{a}\mathbf{a}^T = \begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix}(1\;2\;2) = \begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}
```

```math
P = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}
```

$`\mathbf{a} = \begin{pmatrix} 2 \\ 4 \\ 4 \end{pmatrix}`$이면? $`\mathbf{a}^T\mathbf{a} = 4(1\;2\;2)\begin{pmatrix} 1 \\ 2 \\ 2 \end{pmatrix} = 36`$이고, $`\mathbf{a}\mathbf{a}^T = 4\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}`$. 따라서 $`P = \frac{1}{4 \cdot 9} \cdot 4\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}`$ — **동일하다!**

**$`P^2 = P`$ 검증:**

```math
P^2 = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix}\frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = \frac{1}{81}\begin{pmatrix} 9 & 18 & 18 \\ 18 & 36 & 36 \\ 18 & 36 & 36 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 1 & 2 & 2 \\ 2 & 4 & 4 \\ 2 & 4 & 4 \end{pmatrix} = P
```

**대각합(trace):** $`\text{diag}(P) \cdot \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \frac{1}{9}(1\;4\;4)\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \frac{1}{9} \cdot 9 = 1`$.

**여사영 (complementary projection):** $`P`$가 한 부분공간으로 사영하면, $`I - P`$는 수직 부분공간(직교 여공간)으로 사영한다. $`I - P`$는 $`\mathbf{a}`$에 수직인 평면으로 사영한다.

### 2.4 $`\mathbb{R}^3`$에서의 사영 예제 (Projection in $`\mathbb{R}^3`$ Example)

$`\mathbf{b} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix} \in \mathbb{R}^3`$을 생각하자.

- $`\mathbf{p}_1 = P_1\mathbf{b} = \begin{pmatrix} 0 \\ 0 \\ 4 \end{pmatrix}`$는 $`\mathbf{b}`$의 $`z`$축으로의 사영이다.
- $`\mathbf{p}_2 = P_2\mathbf{b} = \begin{pmatrix} 2 \\ 3 \\ 0 \end{pmatrix}`$는 $`\mathbf{b}`$의 $`xy`$평면으로의 사영이다.

```math
P_1 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad P_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
```

**관찰 사항:**
- $`P_1 + P_2 = I_{3 \times 3}`$
- $`P_1 P_2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}`$ (영행렬)
- $`P_1, P_2`$는 수직이다; $`xy`$평면과 $`z`$축은 직교 부분공간이다
- 직선과 평면은 직교 여공간이다: $`\dim(\text{직선}) + \dim(\text{평면}) = 1 + 2 = 3`$
- 모든 벡터 $`\mathbf{b}`$는 두 부분공간에서의 부분의 합이다:

```math
\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} x \\ y \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \\ z \end{pmatrix} = \mathbf{p}_2 + \mathbf{p}_1
```

### 2.5 부분공간으로의 사영 (Projection onto a Subspace)

$`\mathbb{R}^m`$의 모든 부분공간은 자체의 $`m \times m`$ 사영 행렬 $`P`$를 갖는다. 사영 행렬 $`P`$는 해당 부분을 만든다: $`\mathbf{p} = P\mathbf{b}`$.

부분공간은 기저로 구성된다. 예를 들어:

```math
A_1 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}, \quad A_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}, \quad A_3 = \begin{pmatrix} 2 & 3 \\ 2 & 3 \\ 0 & 0 \end{pmatrix}
```

$`C(A_1)`$은 $`z`$축, $`C(A_2)`$는 $`xy`$평면, $`C(A_3)`$는 $`xy`$평면이다.

LI인 $`n`$개의 벡터 $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3, \ldots, \mathbf{a}_n \in \mathbb{R}^m`$으로 시작한다. 다음 결합을 찾는다:

```math
\mathbf{p} = \alpha_1\mathbf{a}_1 + \alpha_2\mathbf{a}_2 + \cdots + \alpha_n\mathbf{a}_n
```

이것이 주어진 벡터 $`\mathbf{b}`$에 가장 가깝도록 한다. 각 $`\mathbf{b} \in \mathbb{R}^m`$을 $`\mathbf{a}`$들이 생성하는 $`n`$차원 부분공간으로 사영하는 것이다:

```math
A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)
```

$`C(A)`$는 $`\mathbb{R}^m`$의 부분공간이다. $`A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n \in C(A)`$.

$`\mathbf{b}`$에 가장 가까운 특정 결합 $`\mathbf{p} = A\hat{\boldsymbol{\alpha}}`$를 찾는다. $`\hat{\boldsymbol{\alpha}}`$는 $`C(A)`$에서의 최적 벡터이다.

$`n = 1`$일 때: $`\alpha = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}`$.

$`n > 1`$일 때: $`\hat{\boldsymbol{\alpha}} = \begin{pmatrix} \alpha_1 \\ \alpha_2 \\ \vdots \\ \alpha_n \end{pmatrix}`$를 구해야 한다.

### 2.6 $`n = 2`$에서의 유도 (Derivation for $`n = 2`$)

$`n = 2`$, $`A = (\mathbf{a}_1 \; \mathbf{a}_2)`$로 하자.

부분공간 $`S`$는 $`\mathbf{a}_1, \mathbf{a}_2`$로 생성된다:

```math
\mathbf{p} = \alpha_1\mathbf{a}_1 + \alpha_2\mathbf{a}_2 = A\hat{\boldsymbol{\alpha}} \in C(A)
```

오차 벡터 $`\mathbf{e} = \mathbf{b} - \mathbf{p}`$는 부분공간 $`S`$에 수직이다:

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

$`A^T A`$가 가역이면:

```math
\hat{\boldsymbol{\alpha}} = (A^T A)^{-1} A^T \mathbf{b}
```

($`n = 1`$일 때: $`\alpha = (\mathbf{a}^T\mathbf{a})^{-1}\mathbf{a}^T\mathbf{b}`$.)

```math
\mathbf{p} = A\hat{\boldsymbol{\alpha}} = A(A^T A)^{-1} A^T \mathbf{b}
```

```math
P = A(A^T A)^{-1} A^T
```

### 2.7 $`n`$차원 부분공간으로의 확장 (Extension to $`n`$-dimensional Subspace)

$`n`$차원 부분공간으로 쉽게 확장할 수 있으며, $`n`$개의 방정식을 갖는다:

```math
\mathbf{a}_1^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \mathbf{a}_2^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0, \quad \ldots, \quad \mathbf{a}_n^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = 0
```

```math
A^T(\mathbf{b} - A\hat{\boldsymbol{\alpha}}) = \mathbf{0}
```

**비고:**
- $`\mathbf{e} = \mathbf{b} - A\hat{\boldsymbol{\alpha}}`$는 $`\mathcal{N}(A^T)`$에 있으며, 이는 $`C(A)`$에 수직이다.
- $`A`$의 좌영공간(left nullspace)이 오차 벡터를 포함한다.
- $`\mathbf{b}`$는 사영 $`\mathbf{p}`$와 오차 $`\mathbf{e}`$로 분리된다.

### 2.8 예제 2: 부분공간으로의 사영 (Projection onto a Subspace)

```math
A = \begin{pmatrix} 1 & 0 \\ 1 & 1 \\ 1 & 2 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix}
```

$`\hat{\boldsymbol{\alpha}}`$, $`\mathbf{p}`$, $`P`$를 구하라.

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

### 2.9 $`A^T A`$는 $`A`$가 LI 열을 가질 때만 가역 ($`A^T A`$ is Invertible iff $`A`$ has LI Columns)

**참고:** $`(A^T A)^{-1} \neq A^{-1}(A^T)^{-1}`$, 왜냐하면 $`A^{-1}`$은 존재하지 않는다 ($`A`$가 정사각이 아닐 때).

**정리:** $`A^T A`$가 가역일 필요충분조건은 $`A`$가 선형독립인 열을 갖는 것이다.

**증명:**

$`A \in \mathbb{R}^{m \times n}`$이라 하자.

($`\Rightarrow`$) $`A\mathbf{x} = \mathbf{0}`$이면 $`\mathbf{x} \in \mathcal{N}(A)`$. $`A^T`$를 곱하면 $`A^T A\mathbf{x} = \mathbf{0}`$이고, 이는 $`\mathbf{x} \in \mathcal{N}(A^T A)`$를 의미한다. 즉, $`\mathcal{N}(A) \ni \mathbf{x} \longrightarrow \mathbf{x} \in \mathcal{N}(A^T A)`$.

($`\Leftarrow`$) $`A^T A\mathbf{x} = \mathbf{0}`$에서, $`\mathbf{x}^T`$를 곱하면:

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

따라서 $`\mathcal{N}(A) = \lbrace\mathbf{0}\rbrace \iff \mathcal{N}(A^T A) = \lbrace\mathbf{0}\rbrace`$.

$`\iff`$ $`A`$가 가역 $`\iff`$ $`A^T A`$가 가역. $`\square`$

**크기 고려:** $`A \in \mathbb{R}^{m \times n}`$, $`A^T \in \mathbb{R}^{n \times m}`$, $`A^T A \in \mathbb{R}^{n \times n}`$.

$`A^T A`$는 **대칭** (symmetric)이다.

**반례:** $`A = \begin{pmatrix} 1 & 2 \\ 1 & 2 \\ 0 & 0 \end{pmatrix}`$는 종속인 열을 갖는다.

```math
A^T A = \begin{pmatrix} 1 & 1 & 0 \\ 2 & 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 1 & 2 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \end{pmatrix}
```

$`\det(A^T A) = 2 \cdot 8 - 4 \cdot 4 = 0`$. $`A^T A`$는 가역이 **아니며**, $`A^T A`$는 특이(singular)행렬이다.

### 2.10 풀이 예제 (4.2A) (Worked Example)

**문제:** 벡터 $`\mathbf{b} = (3, 4, 4)`$를 $`\mathbf{a} = (2, 2, 1)`$을 지나는 직선 위로, 그 다음 $`\mathbf{a}^* = (1, 0, 0)`$도 포함하는 평면 위로 사영하라. 첫 번째 오차 벡터 $`\mathbf{b} - \mathbf{p}`$가 $`\mathbf{a}`$에 수직인지, 두 번째 오차 벡터 $`\mathbf{e}^* = \mathbf{b} - \mathbf{p}^*`$가 $`\mathbf{a}`$와 $`\mathbf{a}^*`$에도 수직인지 확인하라. $`\mathbf{a}`$와 $`\mathbf{a}^*`$의 평면으로의 $`3 \times 3`$ 사영 행렬 $`P`$를 구하라. 평면으로의 사영이 $`\mathbf{p} = 0`$인 벡터를 찾아라. $`P^2 = P = P^T`$에 주목하라.

**풀이:**

**직선으로의 사영:**

```math
\mathbf{p} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\mathbf{a} = \frac{18}{9}(2, 2, 1) = (4, 4, 2) = 2\mathbf{a}
```

오차 벡터 $`\mathbf{e} = \mathbf{b} - \mathbf{p} = (-1, 0, 2)`$는 $`\mathbf{a} = (2, 2, 1)`$에 수직이다. 확인: $`\mathbf{e}^T\mathbf{a} = -2 + 0 + 2 = 0`$. 따라서 $`\mathbf{p}`$가 맞다.

**평면으로의 사영:** $`\mathbf{a} = (2, 2, 1)`$과 $`\mathbf{a}^* = (1, 0, 0)`$의 평면은 $`A = [\mathbf{a} \;\; \mathbf{a}^*]`$의 열 공간이다:

```math
A = \begin{pmatrix} 2 & 1 \\ 2 & 0 \\ 1 & 0 \end{pmatrix}, \quad A^T A = \begin{pmatrix} 9 & 2 \\ 2 & 1 \end{pmatrix}, \quad (A^T A)^{-1} = \frac{1}{5}\begin{pmatrix} 1 & -2 \\ -2 & 9 \end{pmatrix}
```

```math
P = A(A^T A)^{-1}A^T = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0.8 & 0.4 \\ 0 & 0.4 & 0.2 \end{pmatrix}
```

이제 $`\mathbf{p}^* = P\mathbf{b} = (3, 4.8, 2.4)`$. 오차 $`\mathbf{e}^* = \mathbf{b} - \mathbf{p}^* = (0, -0.8, 1.6)`$은 $`\mathbf{a}`$와 $`\mathbf{a}^*`$에 수직이다.

**검증:**

```math
(\mathbf{e}^*)^T\mathbf{a} = \begin{pmatrix} 0 \\ -0.8 \\ 1.6 \end{pmatrix}^T\begin{pmatrix} 2 \\ 2 \\ 1 \end{pmatrix} = -1.6 + 1.6 = 0
```

```math
(\mathbf{e}^*)^T\mathbf{a}^* = \begin{pmatrix} 0 \\ -0.8 \\ 1.6 \end{pmatrix}^T\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} = 0 + 0 + 0 = 0
```

$`A^T\mathbf{e}^* = \mathbf{0}`$이므로, $`A(A^T A)^{-1}A^T\mathbf{e}^* = \mathbf{0}`$, 즉 $`P\mathbf{e}^* = \mathbf{0}`$이고, $`\mathbf{e}^* \in \mathcal{N}(P)`$이다.

---

<br>
