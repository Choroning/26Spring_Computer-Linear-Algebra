# 제6장 강의 — 고유값과 고유벡터

> **최종 수정일:** 2026-03-31
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 6

> **선수 지식**: [선형대수학] 행렬식, 행렬 연산 (제1-5장).
>
> **학습 목표**:
> 1. 고유값과 고유벡터를 정의하고 계산할 수 있다
> 2. 고유분해를 사용하여 행렬을 대각화할 수 있다
> 3. 고유값을 이용하여 안정성과 동적 특성을 분석할 수 있다
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 6

> **선수 지식**: [선형대수학] 행렬식, 행렬 연산 (제1-5장).
>
> **학습 목표**:
> 1. 고유값과 고유벡터를 정의하고 계산할 수 있다
> 2. 고유분해를 사용하여 행렬을 대각화할 수 있다
> 3. 고유값을 이용하여 안정성과 동적 특성을 분석할 수 있다

> **📑 이 문서는 3개 파트로 나뉘어 있습니다.** GitHub은 파일 하나당 렌더되는 수식 개수에 한도가 있어, 모든 수식이 제대로 보이도록 절 단위로 나눴습니다.
>
> **Part 1** · [Part 2](Concepts-2.ko.md) · [Part 3](Concepts-3.ko.md)

---

<br>

## 목차

- [1. 고유값 소개 (6.1)](#1-고유값-소개-61)
  - [1.1 고유값과 고유벡터의 정의](#11-고유값과-고유벡터의-정의)
  - [1.2 A의 거듭제곱과 고유값](#12-a의-거듭제곱과-고유값)
  - [1.3 고유값의 성질](#13-고유값의-성질)
  - [1.4 고유값을 구하는 방정식](#14-고유값을-구하는-방정식)
  - [1.5 행렬식과 대각합](#15-행렬식과-대각합)
  - [1.6 풀이 예제](#16-풀이-예제)
  - [1.7 허수 고유값](#17-허수-고유값)
  - [1.8 회전 행렬의 고유값](#18-회전-행렬의-고유값)
  - [1.9 AB와 A+B의 고유값](#19-ab와-ab의-고유값)
- [2. 행렬의 대각화 (6.2)](#2-행렬의-대각화-62)
  - [2.1 대각화의 핵심 사실](#21-대각화의-핵심-사실)
  - [2.2 대각화 절차](#22-대각화-절차)
  - [2.3 풀이 예제: 대각화](#23-풀이-예제-대각화)
  - [2.4 대각화에 대한 비고](#24-대각화에-대한-비고)
  - [2.5 증명: 서로 다른 고유값의 고유벡터는 선형독립](#25-증명-서로-다른-고유값의-고유벡터는-선형독립)
  - [2.6 A의 거듭제곱 (마르코프 행렬 예제)](#26-a의-거듭제곱-마르코프-행렬-예제)
  - [2.7 닮은 행렬](#27-닮은-행렬)
  - [2.8 행렬의 거듭제곱과 피보나치 수](#28-행렬의-거듭제곱과-피보나치-수)
  - [2.9 대각화 불가능한 행렬과 중복도](#29-대각화-불가능한-행렬과-중복도)
- [3. 대칭 양정치 행렬 (6.3)](Concepts-2.ko.md#3-대칭-양정치-행렬-63)
  - [3.1 대칭 행렬의 핵심 성질](Concepts-2.ko.md#31-대칭-행렬의-핵심-성질)
  - [3.2 스펙트럼 정리](Concepts-2.ko.md#32-스펙트럼-정리)
  - [3.3 증명: 대칭 행렬은 정규직교 고유기저를 갖는다](Concepts-2.ko.md#33-증명-대칭-행렬은-정규직교-고유기저를-갖는다)
  - [3.4 양정치 행렬의 정의](Concepts-2.ko.md#34-양정치-행렬의-정의)
  - [3.5 양정치 행렬의 성질](Concepts-2.ko.md#35-양정치-행렬의-성질)
  - [3.6 양정치 행렬 판별법](Concepts-2.ko.md#36-양정치-행렬-판별법)
  - [3.7 풀이 예제: 양정치와 양반정치](Concepts-2.ko.md#37-풀이-예제-양정치와-양반정치)
  - [3.8 타원과 이차 형식](Concepts-2.ko.md#38-타원과-이차-형식)
  - [3.9 양정치 행렬과 최솟값 문제](Concepts-2.ko.md#39-양정치-행렬과-최솟값-문제)
  - [3.10 양반정치 행렬](Concepts-2.ko.md#310-양반정치-행렬)
  - [3.11 합동 행렬](Concepts-2.ko.md#311-합동-행렬)
  - [3.12 최적화와 머신 러닝](Concepts-2.ko.md#312-최적화와-머신-러닝)
- [4. 선형 미분방정식 풀기 (6.5)](Concepts-3.ko.md#4-선형-미분방정식-풀기-65)
  - [4.1 핵심 사실](Concepts-3.ko.md#41-핵심-사실)
  - [4.2 스칼라 상미분방정식 복습](Concepts-3.ko.md#42-스칼라-상미분방정식-복습)
  - [4.3 du/dt = Au의 풀이](Concepts-3.ko.md#43-dudt--au의-풀이)
  - [4.4 일반적인 n x n 풀이 절차](Concepts-3.ko.md#44-일반적인-n-x-n-풀이-절차)
  - [4.5 행렬의 지수 함수](Concepts-3.ko.md#45-행렬의-지수-함수)
  - [4.6 이계 방정식](Concepts-3.ko.md#46-이계-방정식)
  - [4.7 2 x 2 행렬의 안정성](Concepts-3.ko.md#47-2-x-2-행렬의-안정성)
  - [4.8 풀이 예제](Concepts-3.ko.md#48-풀이-예제)
- [요약](Concepts-3.ko.md#요약)

---

<br>

## 1. 고유값 소개 (6.1)

### 1.1 고유값과 고유벡터의 정의

$`A`$가 $`\mathbf{x}`$에 작용할 때, 방향을 바꾸지 않고 벡터 $`\mathbf{x}`$를 $`\lambda`$만큼 늘이거나 줄이기만 한다.

```math
A\mathbf{x} = \lambda \mathbf{x}
```

- $`\mathbf{x}`$는 **영이 아닌 벡터** 이며, **고유벡터**(eigenvector)라 한다.
- $`\lambda`$는 $`\mathbf{x}`$에 대응하는 **고유값**(eigenvalue)이다.

**공식 표현:**

1. $`A\mathbf{x} = \lambda\mathbf{x}`$이면, $`\mathbf{x} \neq \mathbf{0}`$는 $`A`$의 고유벡터이고, 수 $`\lambda`$는 고유값이다.

2. $`A^n \mathbf{x} = \lambda^n \mathbf{x}`$ (모든 $`n`$에 대해); $`(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}`$; 그리고 $`\lambda \neq 0`$이면 $`A^{-1}\mathbf{x} = \frac{1}{\lambda}\mathbf{x} = \lambda^{-1}\mathbf{x}`$.

**증명** ($`A^{-1}`$에 대해):

```math
A\mathbf{x} = \lambda\mathbf{x} \implies \mathbf{x} = A^{-1}(A\mathbf{x}) = A^{-1}(\lambda\mathbf{x}) = \lambda A^{-1}\mathbf{x} \quad \square
```

3. $`(A - \lambda I)\mathbf{x} = \mathbf{0} \implies \det(A - \lambda I) = 0`$. 이 방정식은 $`n`$개의 $`\lambda`$를 생성한다.

4. $`A \in \mathbb{R}^{n \times n}`$일 때:

```math
\det(A) = \lambda_1 \lambda_2 \cdots \lambda_n; \quad \text{trace}(A) = \lambda_1 + \lambda_2 + \cdots + \lambda_n
```

5. **사영 행렬**(projection matrix) $`P`$는 $`\lambda = 1`$ 또는 $`\lambda = 0`$을 갖는다. 정사각 행렬 $`P`$가 $`P^2 = P`$를 만족하면 "사영 행렬"이라 한다.

### 1.2 A의 거듭제곱과 고유값

관계식에 $`A`$를 곱하면 어떻게 되는가?

```math
A(A\mathbf{x}) = A(\lambda\mathbf{x}) \implies A^2\mathbf{x} = \lambda A\mathbf{x} = \lambda^2 \mathbf{x}
```

$`\mathbf{x}`$에 $`A`$를 계속 곱하면:

```math
A^3\mathbf{x} = \lambda^3\mathbf{x}, \quad A^k\mathbf{x} = \lambda^k\mathbf{x}, \quad \ldots, \quad A^{100}\mathbf{x} = \lambda^{100}\mathbf{x}
```

```math
\boxed{A^k \mathbf{x} = \lambda^k \mathbf{x}}
```

**$`A^k`$의 거동:**

- $`i = 1, 2, \ldots, n`$에 대해 $`|\lambda_i| < 1`$이면, $`A^k`$는 결국 **영** 에 수렴한다.
- 어떤 $`|\lambda_i| > 1`$이면, $`A^k`$는 결국 **증가** 한다.
- $`\lambda = 1`$이면, 시스템 상태는 시간에 따라 증가하거나 감쇠하지 않고 **일정** 하게 유지된다.

$`A\mathbf{x} = \mathbf{x}`$이면, $`\mathbf{x}`$는 시스템의 **고정점**(fixed point) 또는 **평형점**(equilibrium point)이다. 즉, 시스템이 $`\mathbf{x}`$ 방향에서 정상 상태를 유지한다. 시스템은 **정상 상태**(steady state)에 도달한다: $`A^k \mathbf{x} = \mathbf{x}`$.

### 1.3 고유값의 성질

$`A \in \mathbb{R}^{n \times n}`$ 행렬은 $`n`$개의 고유값을 갖는다. **특성 다항식**(characteristic polynomial)을 풀어 고유값을 구할 수 있다:

```math
\det(A - \lambda I) = 0
```

**행렬의 성질은 고유값에 큰 영향을 미친다:**

**(1)** 행렬 $`A`$의 **대각합**(trace)은 고유값의 **합** 과 같다.

예: $`A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}`$, $`\text{trace}(A) = 1 + 4 = 5`$

$`\det(A - \lambda I) = 0 \implies \lambda_1 = 0, \lambda_2 = 5`$이고, $`\lambda_1 + \lambda_2 = 5`$.

**(2)** $`A`$의 **행렬식**(determinant)은 고유값의 **곱** 과 같다.

예: $`\det(A) = 1 \cdot 4 - 2 \cdot 2 = 0`$이고, $`\lambda_1 \cdot \lambda_2 = 0 \cdot 5 = 0`$.

**(3)** **대칭 행렬** (즉, $`A = A^T`$)은 **실수 고유값** 만 갖는다.

예: $`A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0`$이므로, $`\lambda = \pm 1`$.

**(4)** **대칭 양정치**(SPD) 행렬은 **실수이고 양수인** 고유값을 갖는다. SPD 행렬은 이차 형식(quadratic form) $`f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A\mathbf{x}`$가 **유일한 최솟값** 을 가짐을 보장하므로, 최적화(optimization)에서 매우 중요하다. (6.3절 참고)

**(5)** **대각 행렬**(diagonal matrix)의 경우, 고유값은 **대각 원소** 이다.

예: $`A = \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} a - \lambda & 0 \\ 0 & b - \lambda \end{vmatrix} = (a - \lambda)(b - \lambda) = 0`$이므로, $`\lambda_1 = a, \lambda_2 = b`$.

**삼각 행렬**(triangular matrix)의 경우에도 고유값은 **대각 원소** 이다.

예: $`A = \begin{pmatrix} a & b \\ 0 & c \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} a - \lambda & b \\ 0 & c - \lambda \end{vmatrix} = (a - \lambda)(c - \lambda) = 0`$이므로, $`\lambda_1 = a, \lambda_2 = c`$.

**(6)** **반대칭 행렬**(skew-symmetric matrix, 즉 $`A = -A^T`$)은 **순허수**(purely imaginary) 고유값 또는 **영** 고유값을 갖는다.

예: $`A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ -1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0`$이므로, $`\lambda = \pm i`$.

### 1.4 고유값을 구하는 방정식

$`A\mathbf{x} = \lambda\mathbf{x}`$에서:

```math
A\mathbf{x} - \lambda\mathbf{x} = (A - \lambda I)\mathbf{x} = \mathbf{0}
```

고유벡터는 $`A - \lambda I`$의 **영공간**(nullspace)을 구성한다.

고유값 $`\lambda`$를 알면, $`(A - \lambda I)\mathbf{x} = \mathbf{0}`$을 풀어 고유벡터를 구한다.

$`\mathbf{x}`$가 영이 아닌 벡터이면, $`A - \lambda I`$는 **특이**(singular)하다.

```math
\boxed{\det(A - \lambda I) = 0}
```

이 특성 다항식($`\det(A - \lambda I) = 0`$)은 $`\lambda`$만을 포함하며, $`\lambda`$는 $`A - \lambda I`$의 주대각선을 따라 나타난다. 행렬식은 $`(-\lambda)^n`$을 포함한다.

- 이 방정식은 $`n`$개의 해 $`\lambda_1`$부터 $`\lambda_n`$을 갖는다.
- $`A_{n \times n}`$은 $`n`$개의 고유값을 갖는다.

**증명 (행렬식 = 고유값의 곱):**

```math
A\mathbf{x} = \lambda\mathbf{x} \implies \det(A - \lambda I) = 0
```

$`\lambda`$에 대한 특성 다항식:

```math
\det(\lambda I - A) = 0 \iff \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0 = 0
```

```math
\iff (\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n) = 0
```

$`\lambda = 0`$을 대입하면:

```math
\det(0 \cdot I - A) = \det(-A) = (-1)^n \det(A)
```

```math
(-\lambda_1)(-\lambda_2)\cdots(-\lambda_n) = (-1)^n \lambda_1\lambda_2\cdots\lambda_n
```

```math
\therefore \det(A) = \lambda_1\lambda_2\cdots\lambda_n \quad \square
```

또한, $`(\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n)`$을 전개하면:

```math
= \lambda^n - (\lambda_1 + \lambda_2 + \cdots + \lambda_n)\lambda^{n-1} + \cdots + \det(A) = 0
```

$`\lambda^{n-1}`$의 계수는 $`\text{trace}(A) = a_{11} + a_{22} + \cdots + a_{nn}`$과 같다 (증명 생략).

**행렬의 행렬식의 도함수:**

i) $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$, $`\det(A) = ad - bc`$.

$`a = a(x), b = b(x), c = c(x), d = d(x)`$로 놓으면:

```math
\frac{d}{dx}\det(A) = \frac{d}{dx}(a(x)d(x) - b(x)c(x)) = a'd + ad' - b'c - bc'
```

```math
= \underbrace{a'd - b'c} + \underbrace{ad' - bc'}
```

```math
\begin{vmatrix} a & b \\ c & d \end{vmatrix}' = \begin{vmatrix} a' & b' \\ c & d \end{vmatrix} + \begin{vmatrix} a & b \\ c' & d' \end{vmatrix}
```

행렬식을 미분할 때, 나머지 행은 그대로 두고 **한 행(또는 열)만 미분** 한다.

$`n \times n`$ 행렬의 경우: 행을 $`\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_n`$으로 쓰면:

```math
\begin{vmatrix} a_{11} & \cdots & a_{1n} \\ \vdots & & \vdots \\ a_{n1} & \cdots & a_{nn} \end{vmatrix}' = \sum_{i=1}^{n} \begin{vmatrix} \mathbf{r}_1 \\ \vdots \\ \mathbf{r}_i' \\ \vdots \\ \mathbf{r}_n \end{vmatrix}
```

### 1.5 행렬식과 대각합

**관찰 1:** 소거법(elimination)은 고유값을 **보존하지 않는다**.

```math
A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}
```

$`\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 2 & 6 - \lambda \end{vmatrix} = \lambda^2 - 7\lambda = 0`$이므로, $`\lambda_1 = 7, \lambda_2 = 0`$.

$`\det(A_0 - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 - \lambda = 0`$이므로, $`\lambda_1 = 1, \lambda_2 = 0`$. (다르다!)

**관찰 2:** 곱 $`\lambda_1 \cdot \lambda_2`$과 합 $`\lambda_1 + \lambda_2`$는 행렬로부터 구할 수 있다.

$`A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}`$에 대해:

- $`\lambda_1 \lambda_2 = 7 \cdot 0 = 0 = \det(A) = 6 - 6 = 0`$
- $`\lambda_1 + \lambda_2 = 7 + 0 = 7 = \text{trace}(A) = 1 + 6 = 7`$

$`n`$개 고유값의 곱 $`\lambda_1 \lambda_2 \cdots \lambda_n`$은 $`A`$의 **행렬식** 과 같다.

합 $`\lambda_1 + \lambda_2 + \cdots + \lambda_n`$은 $`n`$개 대각 원소의 합 = $`A`$의 **대각합**(trace)과 같다.

### 1.6 풀이 예제

**예제 1:** 마르코프 행렬(Markov matrix)

```math
A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix}
```

```math
\det(A - \lambda I) = \begin{vmatrix} 0.8 - \lambda & 0.3 \\ 0.2 & 0.7 - \lambda \end{vmatrix} = 0.56 - 1.5\lambda + \lambda^2 - 0.06 = \lambda^2 - \frac{3}{2}\lambda + \frac{1}{2} = (\lambda - \frac{1}{2})(\lambda - 1) = 0
```

```math
\therefore \lambda_1 = 1, \quad \lambda_2 = \frac{1}{2}
```

이는 $`A - \lambda_1 I = A - I`$와 $`A - \lambda_2 I = A - \frac{1}{2}I`$가 **역행렬이 존재하지 않음** 을 의미한다.

고유벡터 $`\mathbf{x}_1`$, $`\mathbf{x}_2`$는 $`(A - I)\mathbf{x}_1 = \mathbf{0}`$과 $`(A - \frac{1}{2}I)\mathbf{x}_2 = \mathbf{0}`$을 만족한다.

즉, $`\mathbf{x}_1 \in \mathcal{N}(A - I)`$이고 $`\mathbf{x}_2 \in \mathcal{N}(A - \frac{1}{2}I)`$이다.

i) $`(A - I)\mathbf{x}_1 = \begin{pmatrix} -0.2 & 0.3 \\ 0.2 & -0.3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`-2x_1 + 3x_2 = 0`$. $`x_1 = 3`$으로 선택하면, $`x_2 = 2`$: $`\mathbf{x}_1 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}`$

ii) $`(A - \frac{1}{2}I)\mathbf{x}_2 = \begin{pmatrix} 0.3 & 0.3 \\ 0.2 & 0.2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_1 + x_2 = 0`$이므로, $`x_1 = -x_2`$. $`x_1 = 1`$로 선택하면, $`x_2 = -1`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

$`A`$를 $`\mathbf{x}_1`$에 곱하면: $`A\mathbf{x}_1 = \mathbf{x}_1`$, $`A^2\mathbf{x}_1 = \mathbf{x}_1`$, ..., $`A^{100}\mathbf{x}_1 = \mathbf{x}_1`$.

$`A`$를 $`\mathbf{x}_2`$에 곱하면: $`A\mathbf{x}_2 = \frac{1}{2}\mathbf{x}_2`$, $`A^2\mathbf{x}_2 = (\frac{1}{2})^2\mathbf{x}_2`$, ..., $`A^{100}\mathbf{x}_2 = (\frac{1}{2})^{100}\mathbf{x}_2`$.

$`\mathbf{x}_1`$과 $`\mathbf{x}_2`$는 모두 자신의 방향을 유지한다.

$`A`$의 고유벡터 $`\mathbf{x}`$는 $`A^n\mathbf{x} = \lambda^n\mathbf{x}`$이므로 모든 $`A^n`$의 고유벡터이기도 하다.

고유벡터 $`\mathbf{x}_1`$과 $`\mathbf{x}_2`$는 $`\mathbb{R}^2`$를 생성(span)한다.

임의의 벡터 $`\mathbf{x}`$는 $`\mathbf{x}_1`$과 $`\mathbf{x}_2`$의 선형결합이다: $`\mathbf{x} = c\mathbf{x}_1 + d\mathbf{x}_2 = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix}`$.

```math
A\mathbf{x} = cA\mathbf{x}_1 + dA\mathbf{x}_2 = c\mathbf{x}_1 + d\tfrac{1}{2}\mathbf{x}_2
```

```math
A^2\mathbf{x} = c\mathbf{x}_1 + d(\tfrac{1}{2})^2\mathbf{x}_2
```

```math
A^{100}\mathbf{x} = c\,\underbrace{\mathbf{x}_1}_{\text{정상 상태}} + d(\tfrac{1}{2})^{100}\underbrace{\mathbf{x}_2}_{\text{감쇠 모드}} \approx c\,\mathbf{x}_1
```

$`\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$로 놓으면:

```math
\begin{pmatrix} 1 \\ 0 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix}
```

```math
\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 & 1 \\ 2 & -3 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 0.2 \\ 0.4 \end{pmatrix}
```

```math
= 0.2\begin{pmatrix} 3 \\ 2 \end{pmatrix} + 0.4\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix} + \begin{pmatrix} 0.4 \\ -0.4 \end{pmatrix}
```

```math
A^{100}\begin{pmatrix} 1 \\ 0 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}
```

마찬가지로 $`\mathbf{x} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$에 대해: $`A^{100}\begin{pmatrix} 0 \\ 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}`$.

```math
A^{100}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_1)
```

$`A`$의 거듭제곱이 높을수록, 열들이 **정상 상태** 에 더 가까워진다.

**예제 2:** 사영 행렬(projection matrix)

```math
P = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}
```

고유값: $`\lambda = 1`$과 $`\lambda = 0`$.

```math
\det(P - \lambda I) = \begin{vmatrix} \frac{1}{2} - \lambda & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} - \lambda \end{vmatrix} = (\frac{1}{2})^2 \begin{vmatrix} 1 - 2\lambda & 1 \\ 1 & 1 - 2\lambda \end{vmatrix} = \frac{1}{4}[(1 - 2\lambda)^2 - 1] = \frac{1}{4}(4\lambda^2 - 4\lambda) = \lambda(\lambda - 1) = 0
```

```math
\therefore \lambda_1 = 1, \; \lambda_2 = 0
```

i) $`\lambda_1 = 1`$: $`(P - I)\mathbf{x}_1 = \begin{pmatrix} -\frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & -\frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = x_2`$, $`x_1 = 1`$로 선택: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`\lambda_2 = 0`$: $`P\mathbf{x}_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 + x_2 = 0`$, $`x_1 = -x_2`$, $`x_1 = 1`$로 선택: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) $`P\mathbf{x}_1 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \mathbf{x}_1`$

$`P\mathbf{x}_2 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}`$

$`\therefore \mathbf{x}_1 \in C(P)`$: 열공간(column space)은 자기 자신으로 사영된다.

$`\mathbf{x}_2 \in \mathcal{N}(P)`$

iv) $`\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}`$로 놓으면

$`P\mathbf{w} = P\begin{pmatrix} 1 \\ -1 \end{pmatrix} + P\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \mathbf{0} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}`$

사영 행렬의 고유값은 오직 **0과 1** 뿐이다.

**예제 3:** 교환 행렬(exchange matrix)

```math
E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
```

고유값: 1과 $`-1`$.

```math
\det(E - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = (\lambda - 1)(\lambda + 1) = 0
```

```math
\therefore \lambda_1 = 1, \; \lambda_2 = -1
```

i) $`\lambda_1 = 1`$: $`(E - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`-x_1 + x_2 = 0`$, $`x_1 = x_2`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`\lambda_2 = -1`$: $`(E + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = -x_2`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) 고유벡터 $`\mathbf{x}_1`$과 $`\mathbf{x}_2`$는 $`P`$의 것과 동일하다.

```math
E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} - \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 2P - I
```

행렬이 $`I`$만큼 이동(shift)하면, 각 $`\lambda`$도 1만큼 이동한다.

$`\det(2P - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & 1 - \lambda \end{vmatrix} = (1 - \lambda)^2 - 1 = \lambda^2 - 2\lambda = \lambda(\lambda - 2) = 0`$

$`\therefore \lambda_1 = 2, \lambda_2 = 0`$ ($`2P`$에 대해), 즉 $`E = 2P - I`$에 대해: $`\lambda_1 = 1, \lambda_2 = -1`$.

**예제 4:** 특이 행렬(singular matrix)

```math
A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}
```

$`\lambda`$와 대응하는 $`\mathbf{x}`$를 구하라.

```math
\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 2 \\ 2 & 4 - \lambda \end{vmatrix} = (1 - \lambda)(4 - \lambda) - 4 = \lambda^2 - 5\lambda = \lambda(\lambda - 5) = 0
```

```math
\therefore \lambda_1 = 5 \text{ 이고 } \lambda_2 = 0
```

i) $`\lambda_1 = 5`$: $`(A - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`2x_1 - x_2 = 0 \implies x_2 = 2x_1`$. $`x_1 = 1`$로 놓으면 $`x_2 = 2`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$

ii) $`\lambda_2 = 0`$: $`(A - 0I)\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_1 + 2x_2 = 0 \implies x_1 = -2x_2`$. $`x_2 = 1`$로 놓으면 $`x_1 = -2`$: $`\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}`$

**비고 1:** $`\mathbf{x}_1`$과 $`\mathbf{x}_2`$는 $`(A - \lambda I)`$의 영공간에 있다.

이 예제에서 $`A`$가 특이(singular)이므로 $`\lambda_2 = 0`$이 고유값이다. $`A`$가 가역이면 "0"은 고유값이 아니다: $`A\mathbf{x} = \mathbf{0} \implies \mathbf{x} = \mathbf{0}`$.

**비고 2:** $`A \in \mathbb{R}^{2 \times 2}`$에서, $`A - \lambda I`$가 특이하면 두 행 모두 벡터 $`(a, b)`$의 배수이다:

```math
\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

고유벡터는 $`(b, -a)`$의 임의 배수이다:

```math
\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} b \\ -a \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

**비고 3:** 이 예제는 두 개의 서로 다른 고유값을 가지며, $`\mathbb{R}^2`$를 생성한다. $`2 \times 2`$ 행렬이 하나의 고유값만 가지면 $`\mathbb{R}^2`$를 생성할 수 없다.

예: $`A = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}`$, $`\det(A - \lambda I) = (3 - \lambda)^2 = 0`$이므로, $`\lambda = 3`$. 기하적 중복도(geometric multiplicity) = 1, 대수적 중복도(algebraic multiplicity) = 2. 행렬 $`A`$는 **결함**(defective)이 있으며 **대각화 불가능** 하다.

### 1.7 허수 고유값

**예제 5:** 90도 회전

```math
R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}
```

실수 고유벡터가 없다.

```math
\det(R - \lambda I) = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0 \implies \lambda = \pm i
```

$`\lambda_1 + \lambda_2 = 0 = \text{trace}(R)`$, $`\lambda_1\lambda_2 = 1 = \det(R)`$.

$`R\mathbf{x} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -y \\ x \end{pmatrix}`$

회전 후, 어떤 실수 벡터 $`R\mathbf{x}`$도 $`\mathbf{x}`$와 같은 방향에 머물지 않으므로, $`R\mathbf{x} = \lambda\mathbf{x}`$를 만족하는 실수 $`\lambda`$는 존재하지 않는다.

i) $`\lambda_1 = i`$: $`(R - iI)\mathbf{x} = \begin{pmatrix} -i & -1 \\ 1 & -i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`-ix_1 = x_2`$, $`x_1 = 1`$로 선택: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ -i \end{pmatrix}`$

ii) $`\lambda_2 = -i`$: $`(R + iI)\mathbf{x} = \begin{pmatrix} i & -1 \\ 1 & i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`ix_1 = x_2`$, $`x_1 = 1`$로 선택: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ i \end{pmatrix}`$

복소 고유벡터 $`\mathbf{x}_1`$과 $`\mathbf{x}_2`$는 복소 공간에서 회전되면서도 방향을 유지한다.

### 1.8 회전 행렬의 고유값

회전 행렬(rotation matrix)은 $`\lambda = e^{i\theta}`$와 $`e^{-i\theta}`$를 갖는다.

```math
R\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}, \quad R\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} -\sin\theta \\ \cos\theta \end{pmatrix}
```

```math
\therefore R = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}
```

```math
\det(R - \lambda I) = \begin{vmatrix} \cos\theta - \lambda & -\sin\theta \\ \sin\theta & \cos\theta - \lambda \end{vmatrix} = (\cos\theta - \lambda)^2 + \sin^2\theta
```

```math
= \lambda^2 - 2\cos\theta\,\lambda + \cos^2\theta + \sin^2\theta = \lambda^2 - 2\cos\theta\,\lambda + 1
```

```math
\therefore \lambda = \cos\theta \pm \sqrt{\cos^2\theta - 1} = \cos\theta \pm \sqrt{-\sin^2\theta} = \cos\theta \pm i\sin\theta
```

```math
\lambda_1 = \cos\theta + i\sin\theta = e^{i\theta}, \quad \lambda_2 = \cos\theta - i\sin\theta = e^{-i\theta}
```

**회전 행렬 $`R`$의 두 가지 성질:**

1. $`R`$은 **직교 행렬**(orthogonal matrix): $`R^T R = I`$, $`|\lambda| = 1`$.
2. $`|\lambda| = 1`$: $`R`$이 직교 행렬이므로 ($`R^T = R^{-1}`$), 모든 고유값은 $`|\lambda| = 1`$을 만족한다. 고유값은 $`\lambda = e^{\pm i\theta}`$이다.

### 1.9 AB와 A+B의 고유값

$`A\mathbf{x} = \lambda\mathbf{x}`$, $`B\mathbf{x} = \beta\mathbf{x}`$를 고려하자.

즉, $`\lambda`$와 $`\beta`$는 각각 $`A`$와 $`B`$의 고유값이다.

```math
AB\mathbf{x} = A(\beta\mathbf{x}) = \beta A\mathbf{x} = \beta\lambda\mathbf{x}
```

**이것은 일반적으로 참이 아니다.** 왜? 일반적으로, $`\mathbf{x}`$는 $`A`$와 $`B`$ **모두** 의 고유벡터가 아니기 때문이다.

마찬가지로, $`(A + B)\mathbf{x} \neq (\lambda + \beta)\mathbf{x}`$.

**비고:** $`A`$와 $`B`$가 동일한 $`n`$개의 선형독립 고유벡터를 공유하는 것은 $`AB = BA`$일 때 그리고 그때에만 성립한다.

---

<br>

## 2. 행렬의 대각화 (6.2)

### 2.1 대각화의 핵심 사실

**(1)** $`AX = X\Lambda`$의 열은 $`A\mathbf{x}_k = \lambda_k\mathbf{x}_k`$이다. 고유값 행렬 $`\Lambda`$는 대각 행렬이다.

**(2)** $`X`$에 있는 $`n`$개의 선형독립 고유벡터가 $`A`$를 대각화한다:

```math
A = X\Lambda X^{-1}
```

```math
AA = X\Lambda X^{-1} X\Lambda X^{-1} = X\Lambda^2 X^{-1}
```

```math
\boxed{A^k = X\Lambda^k X^{-1}}
```

**(3)** $`\mathbf{u}_{k+1} = A\mathbf{u}_k`$를 $`\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0`$로 풀 수 있다.

**(4)** 동일한 고유값이 없으면 $`\implies`$ 고유벡터 $`X`$가 가역 $`\implies`$ $`A`$를 대각화할 수 있다.

중복 고유값 $`\implies`$ $`A`$가 선형독립 고유벡터가 너무 적을 수 있다 $`\implies`$ $`X^{-1}`$이 실패한다.

**(5)** 모든 행렬 $`C = B^{-1}AB`$는 $`A`$와 동일한 고유값을 갖는다. 이러한 $`C`$들은 $`A`$에 **닮은**(similar) 행렬이다.

### 2.2 대각화 절차

$`\mathbf{x}`$가 고유벡터이면, $`A\mathbf{x} = \lambda\mathbf{x}`$이다. $`A`$를 $`\mathbf{x}`$에 적용하는 것은 단지 $`\lambda`$를 곱하는 것이다 — **매우 효율적** 이다.

$`A`$가 대각화 가능하면, $`A^{100}\mathbf{x} = X\Lambda^{100}X^{-1}\mathbf{x}`$—역시 **매우 효율적** 이다.

**대각화:** $`A_{n \times n}`$이 선형독립 고유벡터 $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$을 갖는다고 하자.

$`X = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)`$으로 놓으면:

```math
AX = (A\mathbf{x}_1 \; A\mathbf{x}_2 \; \cdots \; A\mathbf{x}_n) = (\lambda_1\mathbf{x}_1 \; \lambda_2\mathbf{x}_2 \; \cdots \; \lambda_n\mathbf{x}_n) = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix}
```

```math
AX = X\Lambda
```

```math
AXX^{-1} = X\Lambda X^{-1} \implies A = X\Lambda X^{-1}
```

```math
X^{-1}AX = X^{-1}X\Lambda = \Lambda
```

```math
\therefore \Lambda = X^{-1}AX
```

### 2.3 풀이 예제: 대각화

```math
A = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}
```

$`\det(A - \lambda I) = \begin{vmatrix} 2 - \lambda & 4 \\ 0 & 6 - \lambda \end{vmatrix} = (6 - \lambda)(2 - \lambda) = 0`$이므로, $`\lambda_1 = 6, \lambda_2 = 2`$.

i) $`(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`-x_1 + x_2 = 0 \implies \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`(A - \lambda_2 I)\mathbf{x}_2 = \begin{pmatrix} 0 & 4 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_2 = 0`$, $`x_1 = 1`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

iii) $`A\mathbf{x}_1 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = 6\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`A\mathbf{x}_2 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = 2\begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

```math
A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix}
```

```math
\begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 6 & \\ & 2 \end{pmatrix}
```

곱하면: $`X^{-1} = -\begin{pmatrix} 0 & -1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & -1 \end{pmatrix}`$

$`X^{-1}AX = X^{-1}X\Lambda = \Lambda = \begin{pmatrix} 6 & \\ & 2 \end{pmatrix}`$

$`AXX^{-1} = X\Lambda X^{-1}`$

**관찰:** $`A^2 = (X\Lambda X^{-1})(X\Lambda X^{-1}) = X\Lambda^2 X^{-1}`$

$`\implies X^{-1}A^2 X = \Lambda^2`$

$`A^2`$는 $`X`$에 있는 동일한 고유벡터를 가지며, $`\Lambda^2`$에서 제곱된 고유값 $`36, 4`$를 갖는다.

### 2.4 대각화에 대한 비고

**비고:** 행렬 $`X`$는 고유벡터가 **선형독립**(LI)이므로 역행렬이 존재한다.

**비고:** 고유값이 $`n`$개의 서로 다른 수라고 하자. 이는 $`n`$개의 고유벡터가 선형독립임을 의미한다. $`X^{-1}`$이 존재한다. 중복 고유값이 없는 행렬은 대각화할 수 있다.

**비고:** 고유벡터에 영이 아닌 임의의 상수를 곱할 수 있다: $`\alpha A\mathbf{x} = \alpha\lambda\mathbf{x} = A(\alpha\mathbf{x})`$.

**비고:** $`\Lambda`$에서의 고유값 순서는 $`X`$에서의 고유벡터 순서와 동일하다:

```math
A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix} \implies A(\mathbf{x}_2 \; \mathbf{x}_1) = (\mathbf{x}_2 \; \mathbf{x}_1)\begin{pmatrix} \lambda_2 & \\ & \lambda_1 \end{pmatrix}
```

**비고:** 일부 행렬은 고유벡터가 너무 적다. 이러한 행렬은 **대각화할 수 없다**.

예: $`A = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}`$

$`\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & -1 \\ 1 & -1 - \lambda \end{vmatrix} = -(1 - \lambda)(1 + \lambda) + 1 = -(1 - \lambda^2) + 1 = \lambda^2 = 0`$

$`\therefore \lambda_1 = \lambda_2 = 0`$ ($`\lambda`$의 중복).

$`(A - 0I)\mathbf{x} = A\mathbf{x} = \mathbf{0}`$: $`\begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = x_2`$, $`x_1 = 1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

고유벡터가 **하나** 뿐이다. $`A`$는 **대각화할 수 없다**.

**비고:** 가역성(invertibility)은 **영이 아닌** 고유값과 관련된다. $`\lambda_i = 0`$이면, $`\det(A) = \lambda_1\lambda_2\cdots\lambda_n = 0`$이므로 $`A`$는 **특이**(singular)하다.

### 2.5 증명: 서로 다른 고유값의 고유벡터는 선형독립

**명제:** $`A`$가 $`n`$개의 서로 다른 고유값을 갖는 $`n \times n`$ 행렬이면, 대응하는 고유벡터는 **선형독립**(LI)이다.

**증명:** 고유벡터 $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$이 **선형종속**(LD, linearly dependent)이라 가정하자.

$`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p`$가 선형독립이고, $`\mathbf{x}_{p+1}, \mathbf{x}_{p+2}, \ldots, \mathbf{x}_n`$이 선형종속이라 하자.

즉, **모두 영이 아닌** 상수가 존재하여:

```math
\mathbf{x}_{p+1} = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_p\mathbf{x}_p
```

선형결합에 $`A`$를 곱하면:

```math
A\mathbf{x}_{p+1} = \lambda_{p+1}\mathbf{x}_{p+1} = c_1 A\mathbf{x}_1 + c_2 A\mathbf{x}_2 + \cdots + c_p A\mathbf{x}_p = c_1\lambda_1\mathbf{x}_1 + c_2\lambda_2\mathbf{x}_2 + \cdots + c_p\lambda_p\mathbf{x}_p \quad \text{--- (1)}
```

선형결합에 $`\lambda_{p+1}`$을 곱하면:

```math
\lambda_{p+1}\mathbf{x}_{p+1} = c_1\lambda_{p+1}\mathbf{x}_1 + c_2\lambda_{p+1}\mathbf{x}_2 + \cdots + c_p\lambda_{p+1}\mathbf{x}_p \quad \text{--- (2)}
```

(1)에서 (2)를 빼면:

```math
\mathbf{0} = c_1(\lambda_1 - \lambda_{p+1})\mathbf{x}_1 + \cdots + c_p(\lambda_p - \lambda_{p+1})\mathbf{x}_p
```

$`\lambda`$가 서로 다르고 $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p`$가 선형독립이므로:

```math
c_1 = c_2 = \cdots = c_p = 0
```

이는 $`\mathbf{x}_{p+1} = \mathbf{0}`$임을 의미한다. 이는 가정에 모순이다.

따라서 $`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$은 **선형독립** 이다. $`\square`$

### 2.6 A의 거듭제곱 (마르코프 행렬 예제)

```math
A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix} \quad \text{(마르코프 행렬)}
```

$`\lambda_1 = 1, \lambda_2 = 0.5 \implies \Lambda = \begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}`$

$`\mathbf{x}_1 = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}, \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \implies X = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}`$

```math
A = X\Lambda X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}
```

```math
A^2 = X\Lambda^2 X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.25 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}
```

```math
A^k = X\Lambda^k X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1^k & \\ & (0.5)^k \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}
```

$`k \to \infty`$이면, $`(0.5)^k \to 0`$:

```math
A^{\infty} = X\Lambda^{\infty}X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix} = \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix}
```

**Q.** $`A^k`$는 언제 영행렬로 수렴하는가?

**A.** 모든 $`|\lambda_i| < 1`$일 때.

### 2.7 닮은 행렬

$`\Lambda`$가 고정되어 있다고 하자. 고유벡터 행렬 $`X`$를 바꾸면, **같은 고유값** $`\Lambda`$를 갖는 서로 다른 행렬을 얻는다.

```math
A_1 = X_1 \Lambda X_1^{-1}, \quad A_2 = X_2 \Lambda X_2^{-1}, \quad \ldots
```

이들은 **닮은 행렬**(similar matrices)이다: $`\Lambda, A_1, A_2, \ldots`$

모든 $`A_1, A_2, \ldots`$는 $`C`$에 닮았다. 이들은 모두 $`C`$의 고유값을 공유한다.

이 개념을 대각화 불가능한 행렬로 확장한다. 상수 행렬 $`C`$와 가역 행렬 $`B`$를 선택한다. $`A = BCB^{-1}`$을 구성한다. $`A`$와 $`C`$는 **닮았다**.

$`A`$와 $`C`$는 동일한 $`n`$개의 고유값을 갖는다.

**명제:** $`C\mathbf{x} = \lambda\mathbf{x}`$이면, $`BCB^{-1}`$은 새로운 고유벡터 $`B\mathbf{x}`$에 대해 동일한 고유값 $`\lambda`$를 갖는다.

**증명:**

```math
(BCB^{-1})(B\mathbf{x}) = BCI\mathbf{x} = BC\mathbf{x} = B\lambda\mathbf{x} = \lambda(B\mathbf{x}) \quad \square
```

예: $`A_1 = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}`$와 $`A_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}`$는 닮았다.

$`\det(A_1 - \lambda I) = \lambda(1 - \lambda) = 0`$이므로, $`\lambda = 0, 1`$.

$`\det(A_2 - \lambda I) = \lambda^2 - \lambda + \frac{1}{4} - \frac{1}{4} = \lambda(\lambda - 1) = 0`$이므로, $`\lambda = 0, 1`$.

### 2.8 행렬의 거듭제곱과 피보나치 수

피보나치 수는 이전 두 수의 합이다:

```math
0, 1, 1, 2, 3, 5, 8, 13, \ldots
```

```math
a_k + a_{k+1} = a_{k+2}, \quad k \geq 0
```

$`a_{k+1} = a_{k+1}`$을 도입하면:

```math
\begin{cases} a_k + a_{k+1} = a_{k+2} \\ a_{k+1} = a_{k+1} \end{cases}
```

$`\mathbf{u}_k = \begin{pmatrix} a_{k+1} \\ a_k \end{pmatrix}`$로 놓으면:

```math
\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\mathbf{u}_k = \mathbf{u}_{k+1}
```

```math
\therefore \boxed{\mathbf{u}_{k+1} = A\mathbf{u}_k}
```

$`\mathbf{u}_0 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$에서:

$`\mathbf{u}_1 = A\mathbf{u}_0 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`\mathbf{u}_2 = A^2\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{u}_3 = A^3\mathbf{u}_0 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}`$, ...

$`\mathbf{u}_k = A^k\mathbf{u}_0`$

$`A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}`$, $`\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & -\lambda \end{vmatrix} = -\lambda(1 - \lambda) - 1 = \lambda^2 - \lambda - 1 = 0`$

```math
\therefore \lambda = \frac{1 \pm \sqrt{1 + 4}}{2} = \frac{1}{2} \pm \frac{\sqrt{5}}{2}
```

두 개의 서로 다른 고유값 $`\implies`$ 두 개의 선형독립 고유벡터 $`\implies`$ $`X^{-1}`$ 존재 $`\implies`$ $`A = X\Lambda X^{-1}`$.

```math
\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0
```

i) $`\mathbf{u}_0`$을 선형결합 $`X\mathbf{c}`$로 쓰면:

```math
\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = X\mathbf{c} \implies \mathbf{c} = X^{-1}\mathbf{u}_0
```

ii) $`\Lambda^k`$를 $`\mathbf{c}`$에 곱하면:

```math
\begin{pmatrix} \lambda_1^k & \\ & \lambda_2^k \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix}
```

iii) $`X`$를 $`\Lambda^k\mathbf{c}`$에 곱하면:

```math
\mathbf{u}_k = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix} = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2
```

**$`A \in \mathbb{R}^{n \times n}`$으로 일반화:**

```math
\mathbf{u}_k = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2 + \cdots + c_n\lambda_n^k\mathbf{x}_n
```

$`\mathbf{u}_{k+1} = A\mathbf{u}_k`$의 해.

### 2.9 대각화 불가능한 행렬과 중복도

$`\lambda`$가 $`A`$의 고유값이라 하자.

1. **고유벡터 (기하적):** $`A\mathbf{x} = \lambda\mathbf{x}`$, 영이 아닌 $`\mathbf{x}`$.
2. **고유값 (대수적):** $`\det(A - \lambda I) = 0`$.

$`\lambda`$는 단순 고유값이거나, **중복** 고유값일 수 있다 (예: $`\lambda^2 = 0 \to \lambda = 0`$).

**중복도(Multiplicity):**

1. **기하적 중복도 (GM):** $`\lambda`$에 대한 독립 벡터의 수 = $`\dim \mathcal{N}(A - \lambda I)`$.
2. **대수적 중복도 (AM):** 고유값 중 $`\lambda`$의 반복 횟수, 즉 $`\det(A - \lambda I) = 0`$의 $`n`$개 근 중 $`\lambda`$의 수.

예: $`A`$가 $`\lambda = 4, 4, 4`$를 가지면: AM = 3, GM = 1, 2, 또는 3.

예: $`A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}`$, $`|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 = 0`$이므로, $`\lambda = 0, 0`$.

AM = 2, GM = 1 (고유벡터 1개).

$`A\mathbf{x} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$이므로, $`x_2 = 0`$: $`\mathbf{x} = \begin{pmatrix} c \\ 0 \end{pmatrix}`$.

**GM < AM** 이면, $`A`$는 **대각화 불가능** 하다.

---

<br>
