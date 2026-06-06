# 제1장 강의 — 벡터와 행렬

> **최종 수정일:** 2026-06-06
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 1

> **선수 지식**: [선형대수학] 고등학교 수준의 기본 대수학 (방정식, 그래프).
>
> **학습 목표**:
> 1. 벡터를 정의하고 선형 결합을 수행할 수 있다
> 2. 내적을 계산하고 직교성의 기하학적 의미를 이해할 수 있다
> 3. 기본적인 행렬-벡터 곱셈을 수행할 수 있다
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 1

> **선수 지식**: [선형대수학] 고등학교 수준의 기본 대수학 (방정식, 그래프).
>
> **학습 목표**:
> 1. 벡터를 정의하고 선형 결합을 수행할 수 있다
> 2. 내적을 계산하고 직교성의 기하학적 의미를 이해할 수 있다
> 3. 기본적인 행렬-벡터 곱셈을 수행할 수 있다

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.**
>
> **Part 1** · [Part 2](Concepts-Part2.ko.md)

---

<br>

## 목차

- [0. 강의 개요](#0-강의-개요)
- [1. 벡터와 행렬 (워밍업)](#1-벡터와-행렬-워밍업)
- [1.1 벡터와 선형결합 (Linear Combinations)](#11-벡터와-선형결합-linear-combinations)
  - [1.1.1 R^2에서의 선형결합](#111-r2에서의-선형결합)
  - [1.1.2 두 가지 핵심 질문](#112-두-가지-핵심-질문)
  - [1.1.3 고차원으로의 확장](#113-고차원으로의-확장)
  - [1.1.4 두 방정식 풀기](#114-두-방정식-풀기)
  - [1.1.5 소거법이 실패할 수 있는가?](#115-소거법이-실패할-수-있는가)
  - [1.1.6 3차원에서의 벡터](#116-3차원에서의-벡터)
- [1.2 내적으로부터의 길이와 각도](#12-내적으로부터의-길이와-각도)
  - [1.2.1 내적의 정의](#121-내적의-정의)
  - [1.2.2 벡터의 길이](#122-벡터의-길이)
  - [1.2.3 단위벡터 (Unit Vectors)](#123-단위벡터-unit-vectors)
  - [1.2.4 직교벡터 (Perpendicular Vectors)](#124-직교벡터-perpendicular-vectors)
  - [1.2.5 두 벡터 사이의 각도](#125-두-벡터-사이의-각도)
  - [1.2.6 슈바르츠 부등식 (Schwarz Inequality)](#126-슈바르츠-부등식-schwarz-inequality)
  - [1.2.7 삼각 부등식 (Triangle Inequality)](#127-삼각-부등식-triangle-inequality)
  - [1.2.8 3차원에서의 평면](#128-3차원에서의-평면)
- [1.3 행렬과 열공간 (Column Spaces)](Concepts-Part2.ko.md#13-행렬과-열공간-column-spaces)
  - [1.3.1 행렬-벡터 곱셈](Concepts-Part2.ko.md#131-행렬-벡터-곱셈)
  - [1.3.2 열공간 (Column Space)](Concepts-Part2.ko.md#132-열공간-column-space)
  - [1.3.3 독립, 종속, 그리고 열공간](Concepts-Part2.ko.md#133-독립-종속-그리고-열공간)
  - [1.3.4 생성 (Span)](Concepts-Part2.ko.md#134-생성-span)
  - [1.3.5 랭크와 기저 (Rank and Basis)](Concepts-Part2.ko.md#135-랭크와-기저-rank-and-basis)
  - [1.3.6 랭크 1인 행렬](Concepts-Part2.ko.md#136-랭크-1인-행렬)
- [1.4 행렬 곱셈 AB와 CR](Concepts-Part2.ko.md#14-행렬-곱셈-ab와-cr)
  - [1.4.1 행렬 곱셈의 규칙](Concepts-Part2.ko.md#141-행렬-곱셈의-규칙)
  - [1.4.2 AB의 열 해석](Concepts-Part2.ko.md#142-ab의-열-해석)
  - [1.4.3 계산 비용](Concepts-Part2.ko.md#143-계산-비용)
  - [1.4.4 행렬 곱셈의 성질](Concepts-Part2.ko.md#144-행렬-곱셈의-성질)
  - [1.4.5 랭크 1 행렬과 A = CR](Concepts-Part2.ko.md#145-랭크-1-행렬과-a--cr)
  - [1.4.6 C와 R 구하기](Concepts-Part2.ko.md#146-c와-r-구하기)
  - [1.4.7 A의 열 곱하기 B의 행 (외적, Outer Product)](Concepts-Part2.ko.md#147-a의-열-곱하기-b의-행-외적-outer-product)
- [요약](Concepts-Part2.ko.md#요약)

---

<br>

## 0. 강의 개요

**과목:** DCSS321 — 전산선형대수학 (Computer Linear Algebra) (강신후)

**321에서 배울 내용:**
- 벡터 공간 (Vector spaces)
- 선형 변환 (Linear Transformations)
- 연립일차방정식 풀기 (Solve System of Linear Equations)
- 고유값과 고유벡터 (Eigenvalues and Eigenvectors)

**강의 계획 (Part 1):**

| 주차 | 강의 내용 | 핵심 개념 | 자료 |
|:----:|:------------|:-------------|:----------|
| 1 | 벡터 소개 | 벡터 | Chap. 1 |
| 2 | 연립일차방정식 풀기 | $`Ax = b`$ | Chap. 2 |
| 3 | 연립일차방정식 풀기 | $`Ax = b`$ | Chap. 3 |
| 4 | 4대 기본 부분공간 | 열공간 (Column space) | Chap. 3.1--3.2 |
| 5 | 4대 기본 부분공간 | 영공간 (Null space) | Chap. 3.3--3.5 |
| 6 | 직교성 (Orthogonality) | 직교성 | Chap. 4.1--4.3 |
| 7 | 직교성 | 직교성 | Chap. 4.4--4.5 |
| 8 | **중간시험** | | |

**강의 계획 (Part 2):**

| 주차 | 강의 내용 | 핵심 개념 | 자료 |
|:----:|:------------|:-------------|:----------|
| 9 | 행렬식 (Determinants) | 행렬식 | Chap. 5 |
| 10 | 고유값과 고유벡터 | 고유분해 (Eigen decomposition) | Chap. 6.1--6.3 |
| 11 | 고유값과 고유벡터 | 고유분해 | Chap. 6.4--6.5 |
| 12 | 특이값 분해 (SVD) | SVD | Chap. 7 |
| 13 | 선형 변환 | 선형 사상 (Linear map) | Chap. 8 |
| 14 | 선형 변환 | 선형 사상 | Chap. 8 |
| 15 | 최적화에서의 선형대수 | SGD | Chap. 9 |
| 16 | **기말시험** | | |

---

<br>

## 1. 벡터와 행렬 (워밍업)

**제1장의 내용:**
- 1.1 벡터와 선형결합 (Linear Combination)
- 1.2 내적으로부터의 길이와 각도 (Lengths and Angles from Dot Products)
- 1.3 행렬과 열공간 (Matrices and their Column Spaces)
- 1.4 행렬 곱셈 $`AB`$와 $`CR`$ (Matrix Multiplication)

### 핵심 아이디어

- 선형대수는 **벡터** $`\mathbf{v}`$와 **행렬** $`A`$에 관한 것이다.
- 연산을 정의한다: $`+, -, \cdot`$ (덧셈, 뺄셈, 곱셈).
- 벡터와 스칼라를 살펴보자:

```math
\mathbf{v} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

```math
\mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 7 \end{pmatrix} \in \mathbb{R}^2
```

- 스칼라 $`c, d \in \mathbb{R}`$에 대해:

```math
c\mathbf{v} + d\mathbf{w} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix} \in \mathbb{R}^2
```

> 선형결합은 $`xy`$ 평면을 채운다.

### 벡터의 길이

벡터 $`\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} \in \mathbb{R}^2`$의 길이:

```math
\|\mathbf{v}\| = \sqrt{v_1^2 + v_2^2}
```

**예시:** $`\mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}`$, $`\|\mathbf{w}\| = \sqrt{1^2 + 3^2} = \sqrt{10}`$

### 내적 (Dot Product)

$`\mathbf{v}`$와 $`\mathbf{w}`$의 내적:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}^T \begin{pmatrix} w_1 \\ w_2 \end{pmatrix} = v_1 w_1 + v_2 w_2
```

**예시:**

```math
\begin{pmatrix} 2 \\ 4 \end{pmatrix} \cdot \begin{pmatrix} 1 \\ 3 \end{pmatrix} = (2)(1) + (4)(3) = 2 + 12 = 14
```

### 열벡터로 구성된 행렬

행렬 $`A`$는 두 개의 열을 포함한다:

```math
A = \begin{pmatrix} \mathbf{v} & \mathbf{w} \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}
```

### 행렬-벡터 곱은 선형결합

행렬 $`A`$에 벡터 $`\begin{pmatrix} c \\ d \end{pmatrix}`$를 곱하면:

```math
A \begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix} \begin{pmatrix} c \\ d \end{pmatrix} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

이것은 **$`\mathbf{v}`$와 $`\mathbf{w}`$의 선형결합** 이다.

$`\mathbf{x} = \begin{pmatrix} c \\ d \end{pmatrix}`$로 놓으면, 모든 결합 $`A\mathbf{x}`$는 행렬 $`A`$의 **열공간 (Column Space)** 을 생성한다. (열공간은 평면이다!)

### 종속벡터 (Dependent Vectors)

$`\mathbf{z} = \mathbf{v} + \mathbf{w}`$로 놓자.

```math
B = \begin{pmatrix} \mathbf{v} & \mathbf{w} & \mathbf{z} \end{pmatrix} = \begin{pmatrix} 2 & 1 & 3 \\ 4 & 3 & 7 \end{pmatrix}
```

- $`B`$의 열공간은 여전히 $`xy`$ 평면이다.
- $`\mathbf{v}`$와 $`\mathbf{w}`$는 **독립 (Independent)** 벡터이다.
- $`\mathbf{z}`$는 **종속 (Dependent)** 벡터이다.

### 행렬 곱셈 미리보기

행렬 곱셈 $`AB`$는 **$`A`$ 곱하기 $`B`$의 각 열** 로 해석할 수 있다.

---

<br>

## 1.1 벡터와 선형결합 (Linear Combinations)

### 1.1.1 R^2에서의 선형결합

**(1)** $`2\mathbf{v} - 3\mathbf{w}`$는 벡터 $`\mathbf{v}`$와 $`\mathbf{w}`$의 선형결합 $`c\mathbf{v} + d\mathbf{w}`$이다.

**(2)** $`\mathbf{v} = \begin{pmatrix} 4 \\ 1 \end{pmatrix}`$이고 $`\mathbf{w} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$이면:

```math
2\mathbf{v} - 3\mathbf{w} = 2\begin{pmatrix} 4 \\ 1 \end{pmatrix} - 3\begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ -5 \end{pmatrix}
```

**(3)** 모든 결합 $`c\begin{pmatrix} 4 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$는 $`xy`$ 평면을 채운다.

**(4)** 벡터 $`c\begin{pmatrix} 4 \\ 1 \\ 0 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}`$는 $`xyz`$ 공간에서 **평면** 을 채운다.

$`\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$은 그 평면 위에 있지 **않다**.

### 선형결합의 구성

벡터 $`\mathbf{v}`$와 $`\mathbf{w}`$의 선형결합은 두 가지 기본 연산으로 구성된다:
1. **스칼라 곱셈:** $`c\mathbf{v}`$, $`d\mathbf{w}`$
2. **벡터 덧셈:** $`\mathbf{v} + \mathbf{w}`$

이로부터 $`\mathbf{v}`$와 $`\mathbf{w}`$의 **선형결합** $`c\mathbf{v} + d\mathbf{w}`$가 만들어진다.

### 1.1.2 두 가지 핵심 질문

두 가지 질문이 생긴다:

**(1) 기술하라:** 모든 결합 $`c\mathbf{v} + d\mathbf{w}`$를 기술하라.
- 결과는 **평면** 또는 **직선** 이다.

**(2) 구하라:** $`c\mathbf{v} + d\mathbf{w} = \mathbf{x}`$를 만족하는 $`c`$와 $`d`$를 구하라.

**예시:** 다음을 만족하는 $`c`$와 $`d`$를 구하라:

```math
c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 4 \\ 3 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}
```

### 1.1.3 고차원으로의 확장

지금까지 $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$이었다.

벡터의 차원을 높이고 벡터의 개수를 늘려보자:

```math
\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \in \mathbb{R}^m
```

$`m`$차원 공간에서의 $`n`$개의 벡터.

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}
```

$`m`$개의 행과 $`n`$개의 열: $`m \times n`$ 행렬.

고차원에서도 동일한 질문이 제기된다:

**(1) 기술하라:** 모든 결합을 기술하라:

```math
A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{v}_1 + x_2\mathbf{v}_2 + \cdots + x_n\mathbf{v}_n
```

$`A`$의 열들의 결합.

**(2) 구하라:** 다음을 만족하는 $`x_1`$부터 $`x_n`$까지를 구하라:

```math
A\mathbf{x} = \mathbf{b}
```

### 선형결합 $`c\mathbf{v} + d\mathbf{w}`$

$`\mathbf{v} \in \mathbb{R}^2`$를 생각하자. 두 성분을 가진다: $`\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}`$.

기하학적으로 $`\mathbf{v} + \mathbf{w} = \mathbf{w} + \mathbf{v}`$ (교환법칙은 평행사변형 법칙으로 보여진다).

벡터 $`c\mathbf{v}`$는 $`xy`$ 평면에서 무한히 긴 직선을 채운다.

$`\mathbf{w}`$가 그 직선 위에 있지 않으면, 벡터 $`d\mathbf{w}`$는 두 번째 직선을 채운다.

**선형결합 $`c\mathbf{v} + d\mathbf{w}`$는 평면을 채운다.**

**특수한 경우:**
- $`1\mathbf{v} + 1\mathbf{w}`$ = 벡터의 합
- $`1\mathbf{v} - 1\mathbf{w}`$ = 벡터의 차
- $`0\mathbf{v} + 0\mathbf{w}`$ = 영벡터
- $`c\mathbf{v} + 0\mathbf{w}`$ = $`\mathbf{v}`$ 방향의 벡터 $`c\mathbf{v}`$

### 1.1.4 두 방정식 풀기

**풀어보자:**

```math
c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \quad \cdots (*)
```

이것은 다음 연립방정식과 동치이다:

```math
\begin{cases} 2c + 2d = 8 \\ c - d = 2 \end{cases}
```

**소거법으로 풀기:**

첫 번째 방정식을 2로 나누면 $`c + d = 4`$를 얻는다. 그런 다음 두 방정식을 더하면:

```math
c + d = 4
```

```math
c - d = 2
```

더하면: $`2c = 6`$, 따라서 $`c = 3`$.

$`3 + d = 4`$이므로 $`d = 1`$.

**검증:** $`(*)`$는 다음이 된다:

```math
3\begin{pmatrix} 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \checkmark
```

**행렬 형태로:**

```math
\begin{pmatrix} 2 & 2 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix}
```

$`\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}`$, $`\mathbf{b} = \begin{pmatrix} 8 \\ 2 \end{pmatrix}`$로 놓으면:

```math
c\mathbf{v} + d\mathbf{w} = \mathbf{b}
```

```math
\Updownarrow
```

```math
cv_1 + dw_1 = b_1
```

```math
cv_2 + dw_2 = b_2
```

```math
\Updownarrow
```

```math
\begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
```

**기하학적으로 이것은 무엇을 의미하는가?**

연립방정식 $`xv_1 + yw_1 = b_1`$과 $`xv_2 + yw_2 = b_2`$는 두 직선을 나타낸다. 해 $`(c, d)`$는 교점이다:

- 두 일차방정식은 점 $`(c, d)`$에서 만난다.
- $`\mathbf{v}`$와 $`\mathbf{w}`$는 **선형독립 (Linearly Independent)** 이다.
- $`A = \begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}`$는 **가역 (Invertible)** 이다.

### 1.1.5 소거법이 실패할 수 있는가?

$`\mathbf{v} \| \mathbf{w}`$ (즉, $`\frac{v_1}{v_2} = \frac{w_1}{w_2}`$)일 때 **실패한다**.

연립방정식:

```math
xv_1 + yw_1 = b_1
```

```math
xv_2 + yw_2 = b_2
```

두 가지 경우가 있다:
- **해 없음**: $`\mathbf{v}`$와 $`\mathbf{w}`$의 모든 결합이 같은 직선 위에 놓인다. $`\mathbf{b}`$가 그 직선 위에 있지 않으면 해가 없다. $`\mathbf{v}`$와 $`\mathbf{w}`$의 어떤 결합도 $`\mathbf{b}`$와 같지 않다.
- **무한히 많은 해**: $`\mathbf{b}`$가 그 직선 위에 있으면 무한히 많은 해가 존재한다.

### 1.1.6 3차원에서의 벡터

$`\mathbf{v}, \mathbf{w} \in \mathbb{R}^3`$ (3차원 공간)을 가정하자.

```math
\mathbf{v} = \begin{pmatrix} 2 \\ 3 \\ 1 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 4 \\ 1 \end{pmatrix}
```

```math
c\mathbf{v} + d\mathbf{w} = \begin{pmatrix} 2c + d \\ 3c + d \\ c \end{pmatrix}
```

> $`c\mathbf{v} + d\mathbf{w}`$는 전체 3차원 공간을 채우지 **못한다**. 최대 **2차원 평면** 만 채울 수 있다!

**3차원 공간을 채우려면 세 개의 독립 벡터가 필요하다.**

**예시:** $`\hat{i} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}`$, $`\hat{k} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$

```math
c\hat{i} + d\hat{j} + e\hat{k} = \begin{pmatrix} c \\ d \\ e \end{pmatrix}
```

$`\hat{i}, \hat{j}, \hat{k}`$는 3차원 공간에서 $`x, y, z`$ 축 방향의 단위벡터에 해당한다.

```math
\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix} = v_1\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + v_2\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + v_3\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
```

**벡터 형태:** $`= v_1\hat{i} + v_2\hat{j} + v_3\hat{k}`$

**행렬 형태:**

```math
= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}
```

> $`\mathbf{v} = I\mathbf{v}`$ (항등행렬, Identity matrix)

### 평면인지 어떻게 아는가?

영이 아닌 벡터 $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^3`$가 **독립** ($`\mathbf{w}`$가 $`c\mathbf{v}`$의 배수가 아닌)이라고 가정하자.

그러면 이들의 선형결합은 3차원 공간 내의 **평면** (평평한 면)을 채운다.

---

<br>

## 1.2 내적으로부터의 길이와 각도

### 1.2.1 내적의 정의

**(1)** $`\mathbf{v} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 4 \\ 6 \end{pmatrix}`$로 놓자.

$`\mathbf{v}`$와 $`\mathbf{w}`$의 내적:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}^T \begin{pmatrix} 4 \\ 6 \end{pmatrix} = 1 \cdot 4 + 2 \cdot 6 = 16
```

**$`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$에 대한 일반적 정의:**

```math
\mathbf{v} \cdot \mathbf{w} = v_1 w_1 + v_2 w_2
```

**$`\mathbf{v}, \mathbf{w} \in \mathbb{R}^n`$으로의 확장:**

```math
\mathbf{v} \cdot \mathbf{w} = v_1 w_1 + v_2 w_2 + \cdots + v_n w_n = \sum_{i=1}^{n} v_i w_i
```

내적 $`\mathbf{v} \cdot \mathbf{v}`$는 **길이의 제곱** 을 알려준다:

```math
\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + \cdots + v_n^2
```

벡터 $`\mathbf{v}`$의 길이의 제곱.

**예시:** $`\mathbf{v} \in \mathbb{R}^2`$, $`\|\mathbf{v}\|^2 = v_1^2 + v_2^2`$ (피타고라스 공식).

**예시:** $`\mathbf{v} \in \mathbb{R}^3`$, $`\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + v_3^2 = (v_1^2 + v_2^2) + v_3^2`$ 여기서 $`(v_1^2 + v_2^2)`$는 $`xy`$ 평면에서의 값.

### 1.2.2 벡터의 길이

**(2)** $`\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}`$의 길이의 제곱:

```math
\|\mathbf{v}\|^2 = \mathbf{v} \cdot \mathbf{v} = 1^2 + 3^2 + 2^2 = 1 + 9 + 4 = 14
```

```math
\therefore \|\mathbf{v}\| = \sqrt{14}
```

### 1.2.3 단위벡터 (Unit Vectors)

- $`\mathbf{v}`$가 $`\|\mathbf{v}\| = 1`$일 때 **단위벡터 (Unit Vector)** 이다.
- $`\mathbf{v} \neq \mathbf{0}`$이면, $`\frac{\mathbf{v}}{\|\mathbf{v}\|}`$는 단위벡터이다.

**예시 (ex1):** $`\mathbf{u} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}`$는 단위벡터이다.

```math
\|\mathbf{u}\|^2 = \cos^2\theta + \sin^2\theta = 1
```

### 1.2.4 직교벡터 (Perpendicular Vectors)

**(3)** $`\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}`$는 $`\mathbf{w} = \begin{pmatrix} 4 \\ -4 \\ 4 \end{pmatrix}`$에 수직이다:

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}^T \begin{pmatrix} 4 \\ -4 \\ 4 \end{pmatrix} = (1)(4) + (3)(-4) + (2)(4) = 4 - 12 + 8 = 0
```

$`\mathbf{v}`$와 $`\mathbf{w}`$ 사이의 각도가 $`90°`$라고 가정하면:

```math
\cos\theta \Rightarrow \cos 90° = 0
```

```math
\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\| \|\mathbf{w}\| \cos\theta = 0
```

**수직 벡터에 대한 피타고라스 정리:**

```math
\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w})
```

```math
= \mathbf{v} \cdot \mathbf{v} + \mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{v} + \mathbf{w} \cdot \mathbf{w}
```

```math
= \|\mathbf{v}\|^2 + 2\mathbf{v} \cdot \mathbf{w} + \|\mathbf{w}\|^2
```

$`\mathbf{v} \cdot \mathbf{w} = 0`$일 때:

```math
\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 \quad \text{(피타고라스 정리)}
```

$`\mathbf{v} \perp \mathbf{w}`$일 때 $`\|\mathbf{v} - \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2`$도 성립한다.

**예시 (ex2):** $`\mathbf{v} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$는 $`x`$축과 $`45°`$ 각도이다.

$`\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$는 $`x`$축과 $`-45°`$ 각도이다.

```math
\mathbf{v} + \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} + \begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 2 \\ 0 \end{pmatrix}
```

```math
\mathbf{v} - \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} - \begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0 \\ 2 \end{pmatrix}
```

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}^T \begin{pmatrix} 1 \\ -1 \end{pmatrix} = 1 \cdot 1 + 1 \cdot (-1) = 0
```

```math
\|\mathbf{v}\| = \sqrt{2}, \quad \|\mathbf{w}\| = \sqrt{2}
```

```math
\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = 4
```

```math
\|\mathbf{v} - \mathbf{w}\|^2 = 4
```

**예시 (ex3):** $`\mathbf{v} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} -1 \\ 2 \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}^T \begin{pmatrix} -1 \\ 2 \end{pmatrix} = -4 + 4 = 0
```

```math
\Rightarrow \mathbf{v} \perp \mathbf{w}
```

가중치와 거리의 곱 $`v_1 w_1`$과 $`v_2 w_2`$가 균형을 이룬다.

**예시 (ex4):** $`\mathbf{v} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$ (단위벡터), $`\mathbf{w} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \cos\theta
```

$`\|\mathbf{v}\| = \|\mathbf{w}\| = 1`$이면 $`\mathbf{v}`$와 $`\mathbf{w}`$ 사이의 각도는 $`\cos\theta = \mathbf{v} \cdot \mathbf{w}`$이다.

### 1.2.5 두 벡터 사이의 각도

**(4)** $`\mathbf{v} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$과 $`\mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$ 사이의 각도 $`\theta = 45°`$:

```math
\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|} = \frac{1}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}
```

**(5)** 모든 각도는 $`|\cos\theta| \leq 1`$을 만족한다.

모든 벡터는 $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$이고 $`\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|`$이다.

**예시 (ex5) — 응용: 총비용**

가격 벡터 $`\mathbf{p} = \begin{pmatrix} p_1 \\ p_2 \\ p_3 \end{pmatrix}`$, 수량 벡터 $`\mathbf{q} = \begin{pmatrix} q_1 \\ q_2 \\ q_3 \end{pmatrix}`$

$`p_i q_i \Rightarrow`$ 가격 $`p_i`$로 $`q_i`$ 단위를 구매.

```math
\mathbf{p} \cdot \mathbf{q} = p_1 q_1 + p_2 q_2 + p_3 q_3 \Rightarrow \text{총비용}
```

**내적 $`\mathbf{v} \cdot \mathbf{w}`$는 영이 아닌 두 벡터 $`\mathbf{v}`$와 $`\mathbf{w}`$ 사이의 각도를 구해준다.**

**예시 (ex6):** $`\mathbf{v} = \begin{pmatrix} \cos\alpha \\ \sin\alpha \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} \cos\beta \\ \sin\beta \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = \cos\alpha\cos\beta + \sin\alpha\sin\beta
```

삼각함수에 의해: $`= \cos(\beta - \alpha)`$

벡터 사이의 각도는 $`\theta = \beta - \alpha`$이다.

**$`\mathbf{v} \cdot \mathbf{w}`$의 부호는 직각 아래인지 위인지를 알려준다:**

| 경우 | 각도 | 내적 |
|:-----|:------|:------------|
| i) | $`\theta < 90°`$ | $`\mathbf{v} \cdot \mathbf{w} > 0`$ |
| ii) | $`\theta = 90°`$ | $`\mathbf{v} \cdot \mathbf{w} = 0`$ |
| iii) | $`\theta > 90°`$ | $`\mathbf{v} \cdot \mathbf{w} < 0`$ |

$`|\cos\theta| \leq 1`$을 알고 있다. 따라서 영이 아닌 $`\mathbf{v}`$와 $`\mathbf{w}`$에 대해 각도를 다음으로 측정할 수 있다:

```math
\frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \cos\theta
```

### 1.2.6 슈바르츠 부등식 (Schwarz Inequality)

```math
\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\|\|\mathbf{w}\|\cos\theta
```

```math
|\mathbf{v} \cdot \mathbf{w}| = \|\mathbf{v}\|\|\mathbf{w}\||\cos\theta| \leq \|\mathbf{v}\|\|\mathbf{w}\|
```

> **슈바르츠 부등식 (Cauchy--Schwarz--Bunyakovsky):**
> ```math
> |\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|
> ```

**예시 (ex7):** $`\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$에 대해 $`\cos\theta`$를 구하라

```math
\mathbf{v} \cdot \mathbf{w} = 2 \cdot 1 + 1 \cdot 2 = 4
```

```math
\|\mathbf{v}\| = \sqrt{5}, \quad \|\mathbf{w}\| = \sqrt{5}
```

```math
\cos\theta = \frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \frac{4}{5}
```

슈바르츠 부등식에 의해: $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$, 즉 $`4 < 5`$. $`\checkmark`$

**예시 (ex8):** $`\mathbf{v} = \begin{pmatrix} a \\ b \end{pmatrix}`$, $`\mathbf{w} = \begin{pmatrix} b \\ a \end{pmatrix}`$

```math
\mathbf{v} \cdot \mathbf{w} = ab + ba = 2ab
```

```math
\|\mathbf{v}\| = \|\mathbf{w}\| = \sqrt{a^2 + b^2}
```

슈바르츠 부등식 $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$:

```math
|2ab| \leq a^2 + b^2
```

```math
\Leftrightarrow 0 \leq a^2 + b^2 - 2|ab|
```

- $`ab \geq 0`$이면: $`a^2 + b^2 - 2ab = (a - b)^2 \geq 0`$ $`\checkmark`$
- $`ab < 0`$이면: $`a^2 + b^2 + 2ab = (a + b)^2 \geq 0`$ $`\checkmark`$

### 1.2.7 삼각 부등식 (Triangle Inequality)

**삼각 부등식은 슈바르츠 부등식으로부터 직접 도출된다.**

```math
\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w}) = \mathbf{v} \cdot \mathbf{v} + 2\mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{w}
```

```math
\leq \|\mathbf{v}\|^2 + 2|\mathbf{v} \cdot \mathbf{w}| + \|\mathbf{w}\|^2
```

```math
\leq \|\mathbf{v}\|^2 + 2\|\mathbf{v}\|\|\mathbf{w}\| + \|\mathbf{w}\|^2 \quad \text{(슈바르츠에 의해)}
```

```math
= (\|\mathbf{v}\| + \|\mathbf{w}\|)^2
```

> **삼각 부등식:**
> ```math
> \|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
> ```

**ex7로부터의 검증:** 삼각 부등식에 의해:

```math
\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
```

```math
\sqrt{3^2 + 3^2} = \sqrt{18} = 3\sqrt{2} < 2\sqrt{5}
```

```math
\sqrt{18} < \sqrt{20} \quad \checkmark
```

### 1.2.8 3차원에서의 평면

**3차원에서의 평면.** $`\mathbf{n}`$이 단위벡터 $`\mathbf{n} = \begin{pmatrix} n_1 \\ n_2 \\ n_3 \end{pmatrix}`$라 하자.

$`\mathbf{w} \perp \mathbf{n}`$, 즉 $`\mathbf{w} \cdot \mathbf{n} = 0`$을 만족하는 모든 벡터 $`\mathbf{w} \in \mathbb{R}^3`$를 살펴보자.

$`\mathbf{w} \cdot \mathbf{n} = 0`$을 만족하는 벡터 $`\mathbf{w}`$는 $`\mathbb{R}^3`$에서 **2차원 평면** 을 채운다.

```math
\mathbf{w} \cdot \mathbf{n} = w_1 n_1 + w_2 n_2 + w_3 n_3 = 0
```

**예시:** $`\mathbf{n} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$

```math
\mathbf{w} \cdot \mathbf{n} = w_3 = 0
```

또는 $`z = 0`$이며, 이는 $`xy`$ 평면을 나타낸다.

---

<br>
