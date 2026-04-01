# 제1장 강의 — 벡터와 행렬

> **최종 수정일:** 2026-03-31
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
- [1.3 행렬과 열공간 (Column Spaces)](#13-행렬과-열공간-column-spaces)
  - [1.3.1 행렬-벡터 곱셈](#131-행렬-벡터-곱셈)
  - [1.3.2 열공간 (Column Space)](#132-열공간-column-space)
  - [1.3.3 독립, 종속, 그리고 열공간](#133-독립-종속-그리고-열공간)
  - [1.3.4 생성 (Span)](#134-생성-span)
  - [1.3.5 랭크와 기저 (Rank and Basis)](#135-랭크와-기저-rank-and-basis)
  - [1.3.6 랭크 1인 행렬](#136-랭크-1인-행렬)
- [1.4 행렬 곱셈 AB와 CR](#14-행렬-곱셈-ab와-cr)
  - [1.4.1 행렬 곱셈의 규칙](#141-행렬-곱셈의-규칙)
  - [1.4.2 AB의 열 해석](#142-ab의-열-해석)
  - [1.4.3 계산 비용](#143-계산-비용)
  - [1.4.4 행렬 곱셈의 성질](#144-행렬-곱셈의-성질)
  - [1.4.5 랭크 1 행렬과 A = CR](#145-랭크-1-행렬과-a--cr)
  - [1.4.6 C와 R 구하기](#146-c와-r-구하기)
  - [1.4.7 A의 열 곱하기 B의 행 (외적, Outer Product)](#147-a의-열-곱하기-b의-행-외적-outer-product)
- [요약](#요약)

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
| 2 | 연립일차방정식 풀기 | $Ax = b$ | Chap. 2 |
| 3 | 연립일차방정식 풀기 | $Ax = b$ | Chap. 3 |
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
- 1.4 행렬 곱셈 $AB$와 $CR$ (Matrix Multiplication)

### 핵심 아이디어

- 선형대수는 **벡터** $\mathbf{v}$와 **행렬** $A$에 관한 것이다.
- 연산을 정의한다: $+, -, \cdot$ (덧셈, 뺄셈, 곱셈).
- 벡터와 스칼라를 살펴보자:

$$\mathbf{v} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$$

$$\mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 7 \end{pmatrix} \in \mathbb{R}^2$$

- 스칼라 $c, d \in \mathbb{R}$에 대해:

$$c\mathbf{v} + d\mathbf{w} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix} \in \mathbb{R}^2$$

> 선형결합은 $xy$ 평면을 채운다.

### 벡터의 길이

벡터 $\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} \in \mathbb{R}^2$의 길이:

$$\|\mathbf{v}\| = \sqrt{v_1^2 + v_2^2}$$

**예시:** $\mathbf{w} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$, $\|\mathbf{w}\| = \sqrt{1^2 + 3^2} = \sqrt{10}$

### 내적 (Dot Product)

$\mathbf{v}$와 $\mathbf{w}$의 내적:

$$\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}^T \begin{pmatrix} w_1 \\ w_2 \end{pmatrix} = v_1 w_1 + v_2 w_2$$

**예시:**

$$\begin{pmatrix} 2 \\ 4 \end{pmatrix} \cdot \begin{pmatrix} 1 \\ 3 \end{pmatrix} = (2)(1) + (4)(3) = 2 + 12 = 14$$

### 열벡터로 구성된 행렬

행렬 $A$는 두 개의 열을 포함한다:

$$A = \begin{pmatrix} \mathbf{v} & \mathbf{w} \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}$$

### 행렬-벡터 곱은 선형결합

행렬 $A$에 벡터 $\begin{pmatrix} c \\ d \end{pmatrix}$를 곱하면:

$$A \begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix} \begin{pmatrix} c \\ d \end{pmatrix} = c\begin{pmatrix} 2 \\ 4 \end{pmatrix} + d\begin{pmatrix} 1 \\ 3 \end{pmatrix}$$

이것은 **$\mathbf{v}$와 $\mathbf{w}$의 선형결합**이다.

$\mathbf{x} = \begin{pmatrix} c \\ d \end{pmatrix}$로 놓으면, 모든 결합 $A\mathbf{x}$는 행렬 $A$의 **열공간 (Column Space)**을 생성한다. (열공간은 평면이다!)

### 종속벡터 (Dependent Vectors)

$\mathbf{z} = \mathbf{v} + \mathbf{w}$로 놓자.

$$B = \begin{pmatrix} \mathbf{v} & \mathbf{w} & \mathbf{z} \end{pmatrix} = \begin{pmatrix} 2 & 1 & 3 \\ 4 & 3 & 7 \end{pmatrix}$$

- $B$의 열공간은 여전히 $xy$ 평면이다.
- $\mathbf{v}$와 $\mathbf{w}$는 **독립 (Independent)** 벡터이다.
- $\mathbf{z}$는 **종속 (Dependent)** 벡터이다.

### 행렬 곱셈 미리보기

행렬 곱셈 $AB$는 **$A$ 곱하기 $B$의 각 열**로 해석할 수 있다.

---

<br>

## 1.1 벡터와 선형결합 (Linear Combinations)

### 1.1.1 R^2에서의 선형결합

**(1)** $2\mathbf{v} - 3\mathbf{w}$는 벡터 $\mathbf{v}$와 $\mathbf{w}$의 선형결합 $c\mathbf{v} + d\mathbf{w}$이다.

**(2)** $\mathbf{v} = \begin{pmatrix} 4 \\ 1 \end{pmatrix}$이고 $\mathbf{w} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$이면:

$$2\mathbf{v} - 3\mathbf{w} = 2\begin{pmatrix} 4 \\ 1 \end{pmatrix} - 3\begin{pmatrix} 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ -5 \end{pmatrix}$$

**(3)** 모든 결합 $c\begin{pmatrix} 4 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \end{pmatrix}$는 $xy$ 평면을 채운다.

**(4)** 벡터 $c\begin{pmatrix} 4 \\ 1 \\ 0 \end{pmatrix} + d\begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}$는 $xyz$ 공간에서 **평면**을 채운다.

$\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$은 그 평면 위에 있지 **않다**.

### 선형결합의 구성

벡터 $\mathbf{v}$와 $\mathbf{w}$의 선형결합은 두 가지 기본 연산으로 구성된다:
1. **스칼라 곱셈:** $c\mathbf{v}$, $d\mathbf{w}$
2. **벡터 덧셈:** $\mathbf{v} + \mathbf{w}$

이로부터 $\mathbf{v}$와 $\mathbf{w}$의 **선형결합** $c\mathbf{v} + d\mathbf{w}$가 만들어진다.

### 1.1.2 두 가지 핵심 질문

두 가지 질문이 생긴다:

**(1) 기술하라:** 모든 결합 $c\mathbf{v} + d\mathbf{w}$를 기술하라.
- 결과는 **평면** 또는 **직선**이다.

**(2) 구하라:** $c\mathbf{v} + d\mathbf{w} = \mathbf{x}$를 만족하는 $c$와 $d$를 구하라.

**예시:** 다음을 만족하는 $c$와 $d$를 구하라:

$$c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 4 \\ 3 \end{pmatrix} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}$$

### 1.1.3 고차원으로의 확장

지금까지 $\mathbf{v}, \mathbf{w} \in \mathbb{R}^2$이었다.

벡터의 차원을 높이고 벡터의 개수를 늘려보자:

$$\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \in \mathbb{R}^m$$

$m$차원 공간에서의 $n$개의 벡터.

$$A = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}$$

$m$개의 행과 $n$개의 열: $m \times n$ 행렬.

고차원에서도 동일한 질문이 제기된다:

**(1) 기술하라:** 모든 결합을 기술하라:

$$A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \\ | & | & & | \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{v}_1 + x_2\mathbf{v}_2 + \cdots + x_n\mathbf{v}_n$$

$A$의 열들의 결합.

**(2) 구하라:** 다음을 만족하는 $x_1$부터 $x_n$까지를 구하라:

$$A\mathbf{x} = \mathbf{b}$$

### 선형결합 $c\mathbf{v} + d\mathbf{w}$

$\mathbf{v} \in \mathbb{R}^2$를 생각하자. 두 성분을 가진다: $\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix}$.

기하학적으로 $\mathbf{v} + \mathbf{w} = \mathbf{w} + \mathbf{v}$ (교환법칙은 평행사변형 법칙으로 보여진다).

벡터 $c\mathbf{v}$는 $xy$ 평면에서 무한히 긴 직선을 채운다.

$\mathbf{w}$가 그 직선 위에 있지 않으면, 벡터 $d\mathbf{w}$는 두 번째 직선을 채운다.

**선형결합 $c\mathbf{v} + d\mathbf{w}$는 평면을 채운다.**

**특수한 경우:**
- $1\mathbf{v} + 1\mathbf{w}$ = 벡터의 합
- $1\mathbf{v} - 1\mathbf{w}$ = 벡터의 차
- $0\mathbf{v} + 0\mathbf{w}$ = 영벡터
- $c\mathbf{v} + 0\mathbf{w}$ = $\mathbf{v}$ 방향의 벡터 $c\mathbf{v}$

### 1.1.4 두 방정식 풀기

**풀어보자:**

$$c\begin{pmatrix} 2 \\ 1 \end{pmatrix} + d\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \quad \cdots (*)$$

이것은 다음 연립방정식과 동치이다:

$$\begin{cases} 2c + 2d = 8 \\ c - d = 2 \end{cases}$$

**소거법으로 풀기:**

첫 번째 방정식을 2로 나누면 $c + d = 4$를 얻는다. 그런 다음 두 방정식을 더하면:

$$c + d = 4$$

$$c - d = 2$$

더하면: $2c = 6$, 따라서 $c = 3$.

$3 + d = 4$이므로 $d = 1$.

**검증:** $(*)$는 다음이 된다:

$$3\begin{pmatrix} 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix} \checkmark$$

**행렬 형태로:**

$$\begin{pmatrix} 2 & 2 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 8 \\ 2 \end{pmatrix}$$

$\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$, $\mathbf{w} = \begin{pmatrix} 2 \\ -1 \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 8 \\ 2 \end{pmatrix}$로 놓으면:

$$c\mathbf{v} + d\mathbf{w} = \mathbf{b}$$

$$\Updownarrow$$

$$cv_1 + dw_1 = b_1$$

$$cv_2 + dw_2 = b_2$$

$$\Updownarrow$$

$$\begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$$

**기하학적으로 이것은 무엇을 의미하는가?**

연립방정식 $xv_1 + yw_1 = b_1$과 $xv_2 + yw_2 = b_2$는 두 직선을 나타낸다. 해 $(c, d)$는 교점이다:

- 두 일차방정식은 점 $(c, d)$에서 만난다.
- $\mathbf{v}$와 $\mathbf{w}$는 **일차독립 (Linearly Independent)**이다.
- $A = \begin{pmatrix} v_1 & w_1 \\ v_2 & w_2 \end{pmatrix}$는 **가역 (Invertible)**이다.

### 1.1.5 소거법이 실패할 수 있는가?

$\mathbf{v} \| \mathbf{w}$ (즉, $\frac{v_1}{v_2} = \frac{w_1}{w_2}$)일 때 **실패한다**.

연립방정식:

$$xv_1 + yw_1 = b_1$$

$$xv_2 + yw_2 = b_2$$

두 가지 경우가 있다:
- **해 없음**: $\mathbf{v}$와 $\mathbf{w}$의 모든 결합이 같은 직선 위에 놓인다. $\mathbf{b}$가 그 직선 위에 있지 않으면 해가 없다. $\mathbf{v}$와 $\mathbf{w}$의 어떤 결합도 $\mathbf{b}$와 같지 않다.
- **무한히 많은 해**: $\mathbf{b}$가 그 직선 위에 있으면 무한히 많은 해가 존재한다.

### 1.1.6 3차원에서의 벡터

$\mathbf{v}, \mathbf{w} \in \mathbb{R}^3$ (3차원 공간)을 가정하자.

$$\mathbf{v} = \begin{pmatrix} 2 \\ 3 \\ 1 \end{pmatrix}, \quad \mathbf{w} = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad \mathbf{v} + \mathbf{w} = \begin{pmatrix} 3 \\ 4 \\ 1 \end{pmatrix}$$

$$c\mathbf{v} + d\mathbf{w} = \begin{pmatrix} 2c + d \\ 3c + d \\ c \end{pmatrix}$$

> $c\mathbf{v} + d\mathbf{w}$는 전체 3차원 공간을 채우지 **못한다**. 최대 **2차원 평면**만 채울 수 있다!

**3차원 공간을 채우려면 세 개의 독립 벡터가 필요하다.**

**예시:** $\hat{i} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$, $\hat{k} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$

$$c\hat{i} + d\hat{j} + e\hat{k} = \begin{pmatrix} c \\ d \\ e \end{pmatrix}$$

$\hat{i}, \hat{j}, \hat{k}$는 3차원 공간에서 $x, y, z$ 축 방향의 단위벡터에 해당한다.

$$\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix} = v_1\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + v_2\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + v_3\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$$

**벡터 형태:** $= v_1\hat{i} + v_2\hat{j} + v_3\hat{k}$

**행렬 형태:**

$$= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}$$

> $\mathbf{v} = I\mathbf{v}$ (항등행렬, Identity matrix)

### 평면인지 어떻게 아는가?

영이 아닌 벡터 $\mathbf{v}, \mathbf{w} \in \mathbb{R}^3$가 **독립** ($\mathbf{w}$가 $c\mathbf{v}$의 배수가 아닌)이라고 가정하자.

그러면 이들의 선형결합은 3차원 공간 내의 **평면** (평평한 면)을 채운다.

---

<br>

## 1.2 내적으로부터의 길이와 각도

### 1.2.1 내적의 정의

**(1)** $\mathbf{v} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$, $\mathbf{w} = \begin{pmatrix} 4 \\ 6 \end{pmatrix}$로 놓자.

$\mathbf{v}$와 $\mathbf{w}$의 내적:

$$\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}^T \begin{pmatrix} 4 \\ 6 \end{pmatrix} = 1 \cdot 4 + 2 \cdot 6 = 16$$

**$\mathbf{v}, \mathbf{w} \in \mathbb{R}^2$에 대한 일반적 정의:**

$$\mathbf{v} \cdot \mathbf{w} = v_1 w_1 + v_2 w_2$$

**$\mathbf{v}, \mathbf{w} \in \mathbb{R}^n$으로의 확장:**

$$\mathbf{v} \cdot \mathbf{w} = v_1 w_1 + v_2 w_2 + \cdots + v_n w_n = \sum_{i=1}^{n} v_i w_i$$

내적 $\mathbf{v} \cdot \mathbf{v}$는 **길이의 제곱**을 알려준다:

$$\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + \cdots + v_n^2$$

벡터 $\mathbf{v}$의 길이의 제곱.

**예시:** $\mathbf{v} \in \mathbb{R}^2$, $\|\mathbf{v}\|^2 = v_1^2 + v_2^2$ (피타고라스 공식).

**예시:** $\mathbf{v} \in \mathbb{R}^3$, $\|\mathbf{v}\|^2 = v_1^2 + v_2^2 + v_3^2 = (v_1^2 + v_2^2) + v_3^2$ 여기서 $(v_1^2 + v_2^2)$는 $xy$ 평면에서의 값.

### 1.2.2 벡터의 길이

**(2)** $\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}$의 길이의 제곱:

$$\|\mathbf{v}\|^2 = \mathbf{v} \cdot \mathbf{v} = 1^2 + 3^2 + 2^2 = 1 + 9 + 4 = 14$$

$$\therefore \|\mathbf{v}\| = \sqrt{14}$$

### 1.2.3 단위벡터 (Unit Vectors)

- $\mathbf{v}$가 $\|\mathbf{v}\| = 1$일 때 **단위벡터 (Unit Vector)**이다.
- $\mathbf{v} \neq \mathbf{0}$이면, $\frac{\mathbf{v}}{\|\mathbf{v}\|}$는 단위벡터이다.

**예시 (ex1):** $\mathbf{u} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}$는 단위벡터이다.

$$\|\mathbf{u}\|^2 = \cos^2\theta + \sin^2\theta = 1$$

### 1.2.4 직교벡터 (Perpendicular Vectors)

**(3)** $\mathbf{v} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}$는 $\mathbf{w} = \begin{pmatrix} 4 \\ -4 \\ 4 \end{pmatrix}$에 수직이다:

$$\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}^T \begin{pmatrix} 4 \\ -4 \\ 4 \end{pmatrix} = (1)(4) + (3)(-4) + (2)(4) = 4 - 12 + 8 = 0$$

$\mathbf{v}$와 $\mathbf{w}$ 사이의 각도가 $90°$라고 가정하면:

$$\cos\theta \Rightarrow \cos 90° = 0$$

$$\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\| \|\mathbf{w}\| \cos\theta = 0$$

**수직 벡터에 대한 피타고라스 정리:**

$$\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w})$$

$$= \mathbf{v} \cdot \mathbf{v} + \mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{v} + \mathbf{w} \cdot \mathbf{w}$$

$$= \|\mathbf{v}\|^2 + 2\mathbf{v} \cdot \mathbf{w} + \|\mathbf{w}\|^2$$

$\mathbf{v} \cdot \mathbf{w} = 0$일 때:

$$\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 \quad \text{(피타고라스 정리)}$$

$\mathbf{v} \perp \mathbf{w}$일 때 $\|\mathbf{v} - \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2$도 성립한다.

**예시 (ex2):** $\mathbf{v} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$는 $x$축과 $45°$ 각도이다.

$\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$는 $x$축과 $-45°$ 각도이다.

$$\mathbf{v} + \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} + \begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 2 \\ 0 \end{pmatrix}$$

$$\mathbf{v} - \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} - \begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0 \\ 2 \end{pmatrix}$$

$$\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}^T \begin{pmatrix} 1 \\ -1 \end{pmatrix} = 1 \cdot 1 + 1 \cdot (-1) = 0$$

$$\|\mathbf{v}\| = \sqrt{2}, \quad \|\mathbf{w}\| = \sqrt{2}$$

$$\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = 4$$

$$\|\mathbf{v} - \mathbf{w}\|^2 = 4$$

**예시 (ex3):** $\mathbf{v} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$, $\mathbf{w} = \begin{pmatrix} -1 \\ 2 \end{pmatrix}$

$$\mathbf{v} \cdot \mathbf{w} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}^T \begin{pmatrix} -1 \\ 2 \end{pmatrix} = -4 + 4 = 0$$

$$\Rightarrow \mathbf{v} \perp \mathbf{w}$$

가중치와 거리의 곱 $v_1 w_1$과 $v_2 w_2$가 균형을 이룬다.

**예시 (ex4):** $\mathbf{v} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$ (단위벡터), $\mathbf{w} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}$

$$\mathbf{v} \cdot \mathbf{w} = \cos\theta$$

$\|\mathbf{v}\| = \|\mathbf{w}\| = 1$이면 $\mathbf{v}$와 $\mathbf{w}$ 사이의 각도는 $\cos\theta = \mathbf{v} \cdot \mathbf{w}$이다.

### 1.2.5 두 벡터 사이의 각도

**(4)** $\mathbf{v} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$과 $\mathbf{w} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$ 사이의 각도 $\theta = 45°$:

$$\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|} = \frac{1}{1 \cdot \sqrt{2}} = \frac{1}{\sqrt{2}}$$

**(5)** 모든 각도는 $|\cos\theta| \leq 1$을 만족한다.

모든 벡터는 $|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|$이고 $\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|$이다.

**예시 (ex5) — 응용: 총비용**

가격 벡터 $\mathbf{p} = \begin{pmatrix} p_1 \\ p_2 \\ p_3 \end{pmatrix}$, 수량 벡터 $\mathbf{q} = \begin{pmatrix} q_1 \\ q_2 \\ q_3 \end{pmatrix}$

$p_i q_i \Rightarrow$ 가격 $p_i$로 $q_i$ 단위를 구매.

$$\mathbf{p} \cdot \mathbf{q} = p_1 q_1 + p_2 q_2 + p_3 q_3 \Rightarrow \text{총비용}$$

**내적 $\mathbf{v} \cdot \mathbf{w}$는 영이 아닌 두 벡터 $\mathbf{v}$와 $\mathbf{w}$ 사이의 각도를 구해준다.**

**예시 (ex6):** $\mathbf{v} = \begin{pmatrix} \cos\alpha \\ \sin\alpha \end{pmatrix}$, $\mathbf{w} = \begin{pmatrix} \cos\beta \\ \sin\beta \end{pmatrix}$

$$\mathbf{v} \cdot \mathbf{w} = \cos\alpha\cos\beta + \sin\alpha\sin\beta$$

삼각함수에 의해: $= \cos(\beta - \alpha)$

벡터 사이의 각도는 $\theta = \beta - \alpha$이다.

**$\mathbf{v} \cdot \mathbf{w}$의 부호는 직각 아래인지 위인지를 알려준다:**

| 경우 | 각도 | 내적 |
|:-----|:------|:------------|
| i) | $\theta < 90°$ | $\mathbf{v} \cdot \mathbf{w} > 0$ |
| ii) | $\theta = 90°$ | $\mathbf{v} \cdot \mathbf{w} = 0$ |
| iii) | $\theta > 90°$ | $\mathbf{v} \cdot \mathbf{w} < 0$ |

$|\cos\theta| \leq 1$을 알고 있다. 따라서 영이 아닌 $\mathbf{v}$와 $\mathbf{w}$에 대해 각도를 다음으로 측정할 수 있다:

$$\frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \cos\theta$$

### 1.2.6 슈바르츠 부등식 (Schwarz Inequality)

$$\mathbf{v} \cdot \mathbf{w} = \|\mathbf{v}\|\|\mathbf{w}\|\cos\theta$$

$$|\mathbf{v} \cdot \mathbf{w}| = \|\mathbf{v}\|\|\mathbf{w}\||\cos\theta| \leq \|\mathbf{v}\|\|\mathbf{w}\|$$

> **슈바르츠 부등식 (Cauchy--Schwarz--Bunyakovsky):**
> $$|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|$$

**예시 (ex7):** $\mathbf{v} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$, $\mathbf{w} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$에 대해 $\cos\theta$를 구하라

$$\mathbf{v} \cdot \mathbf{w} = 2 \cdot 1 + 1 \cdot 2 = 4$$

$$\|\mathbf{v}\| = \sqrt{5}, \quad \|\mathbf{w}\| = \sqrt{5}$$

$$\cos\theta = \frac{\mathbf{v}}{\|\mathbf{v}\|} \cdot \frac{\mathbf{w}}{\|\mathbf{w}\|} = \frac{4}{5}$$

슈바르츠 부등식에 의해: $|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|$, 즉 $4 < 5$. $\checkmark$

**예시 (ex8):** $\mathbf{v} = \begin{pmatrix} a \\ b \end{pmatrix}$, $\mathbf{w} = \begin{pmatrix} b \\ a \end{pmatrix}$

$$\mathbf{v} \cdot \mathbf{w} = ab + ba = 2ab$$

$$\|\mathbf{v}\| = \|\mathbf{w}\| = \sqrt{a^2 + b^2}$$

슈바르츠 부등식 $|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|$:

$$|2ab| \leq a^2 + b^2$$

$$\Leftrightarrow 0 \leq a^2 + b^2 - 2|ab|$$

- $ab \geq 0$이면: $a^2 + b^2 - 2ab = (a - b)^2 \geq 0$ $\checkmark$
- $ab < 0$이면: $a^2 + b^2 + 2ab = (a + b)^2 \geq 0$ $\checkmark$

### 1.2.7 삼각 부등식 (Triangle Inequality)

**삼각 부등식은 슈바르츠 부등식으로부터 직접 도출된다.**

$$\|\mathbf{v} + \mathbf{w}\|^2 = (\mathbf{v} + \mathbf{w}) \cdot (\mathbf{v} + \mathbf{w}) = \mathbf{v} \cdot \mathbf{v} + 2\mathbf{v} \cdot \mathbf{w} + \mathbf{w} \cdot \mathbf{w}$$

$$\leq \|\mathbf{v}\|^2 + 2|\mathbf{v} \cdot \mathbf{w}| + \|\mathbf{w}\|^2$$

$$\leq \|\mathbf{v}\|^2 + 2\|\mathbf{v}\|\|\mathbf{w}\| + \|\mathbf{w}\|^2 \quad \text{(슈바르츠에 의해)}$$

$$= (\|\mathbf{v}\| + \|\mathbf{w}\|)^2$$

> **삼각 부등식:**
> $$\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|$$

**ex7로부터의 검증:** 삼각 부등식에 의해:

$$\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|$$

$$\sqrt{3^2 + 3^2} = \sqrt{18} = 3\sqrt{2} < 2\sqrt{5}$$

$$\sqrt{18} < \sqrt{20} \quad \checkmark$$

### 1.2.8 3차원에서의 평면

**3차원에서의 평면.** $\mathbf{n}$이 단위벡터 $\mathbf{n} = \begin{pmatrix} n_1 \\ n_2 \\ n_3 \end{pmatrix}$라 하자.

$\mathbf{w} \perp \mathbf{n}$, 즉 $\mathbf{w} \cdot \mathbf{n} = 0$을 만족하는 모든 벡터 $\mathbf{w} \in \mathbb{R}^3$를 살펴보자.

$\mathbf{w} \cdot \mathbf{n} = 0$을 만족하는 벡터 $\mathbf{w}$는 $\mathbb{R}^3$에서 **2차원 평면**을 채운다.

$$\mathbf{w} \cdot \mathbf{n} = w_1 n_1 + w_2 n_2 + w_3 n_3 = 0$$

**예시:** $\mathbf{n} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$

$$\mathbf{w} \cdot \mathbf{n} = w_3 = 0$$

또는 $z = 0$이며, 이는 $xy$ 평면을 나타낸다.

---

<br>

## 1.3 행렬과 열공간 (Column Spaces)

### 1.3.1 행렬-벡터 곱셈

**(1)** $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}$는 **3 x 2 행렬**이다: 3개의 행과 2개의 열, 랭크 2.

**(2)** $A\mathbf{x}$의 3개의 성분은 $A$의 3개의 행과 벡터 $\mathbf{x}$의 내적이다.

**예시:** $\mathbf{x} = \begin{pmatrix} 7 \\ 8 \end{pmatrix}$

$$A\mathbf{x} = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 7 \\ 8 \end{pmatrix} = \begin{pmatrix} 7 + 16 \\ 21 + 32 \\ 35 + 48 \end{pmatrix} = \begin{pmatrix} 23 \\ 53 \\ 83 \end{pmatrix}$$

**(3)** $A\mathbf{x}$는 $A$의 **열들의 결합**이다:

$$\begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 7 \\ 8 \end{pmatrix} = 7\begin{pmatrix} 1 \\ 3 \\ 5 \end{pmatrix} + 8\begin{pmatrix} 2 \\ 4 \\ 6 \end{pmatrix}$$

**(4)** $A$의 열공간은 열들의 모든 결합 $A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2$를 포함한다.

**(5)** 랭크 1 행렬: $A$의 모든 열이 하나의 직선 위에 있다.

### $m \times n$ 행렬 살펴보기

$$A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

$m = n$이면: **정방행렬 (Square matrix)**.

**정방행렬의 예시 ($A \in \mathbb{R}^{3 \times 3}$):**

$$\text{항등행렬 (Identity): } \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \qquad \text{대각행렬 (Diagonal): } \begin{pmatrix} 2 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 0 & 5 \end{pmatrix}$$

$$\text{삼각행렬 (Triangular): } \begin{pmatrix} 2 & 1 & -3 \\ 0 & 4 & 7 \\ 0 & 0 & 5 \end{pmatrix} \qquad \text{대칭행렬 (Symmetric): } \begin{pmatrix} 2 & 1 & -3 \\ 1 & 4 & 7 \\ -3 & 7 & 5 \end{pmatrix}$$

### 열을 벡터로 해석하기

$A$의 열들을 벡터로 해석할 수 있다:

$$A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}$$

여기서 $\mathbf{a}_i \in \mathbb{R}^m$, $i = 1, 2, \ldots, n$.

**질문:** $m \times n$ 행렬 $A$는 $n \times 1$ 벡터 $\mathbf{x}$를 어떻게 곱하는가?

1. $\mathbf{x}$와 $A$의 행들의 **내적**
2. $A$의 열들의 **선형결합**

### $A\mathbf{x}$의 두 가지 관점

**(1) 행 관점 (내적):**

$$A\mathbf{x} = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}$$

3개의 행 $\Rightarrow$ 3개의 내적.

$A$의 각 행은 벡터 $\mathbf{x}$와 같은 수의 성분을 가진다.

**(2) 열 관점 (선형결합):**

$$A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n$$

$A$의 **열들의 결합**.

**예시:** $A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}$

$$A\mathbf{x} = x_1\begin{pmatrix} -1 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix} + x_4\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}$$

### 1.3.2 열공간 (Column Space)

$A = (\mathbf{a}_1 \ \mathbf{a}_2)$이고 $\mathbf{a}_1$과 $\mathbf{a}_2$가 일차독립일 때:

$$A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2$$

는 **모든 결합의 평면** (생성, Span)을 보여준다.

### 1.3.3 독립, 종속, 그리고 열공간

**예시 (ex1):**

$$A_1 = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 4 & 0 \\ 3 & 5 & 6 \end{pmatrix}$$

각 열이 새로운 방향을 제시한다. 이들의 결합은 **3차원 공간**을 채운다.

$$A_1\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = x_1\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix} + x_2\begin{pmatrix} 0 \\ 4 \\ 5 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 0 \\ 6 \end{pmatrix}$$

**독립 열:** $A_1\mathbf{x} = \mathbf{0}$이면 $\mathbf{x} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$일 때만 성립

**예시 (ex2):**

$$A_2 = \begin{pmatrix} 1 & 2 & 3 \\ 1 & 4 & 5 \\ 6 & 0 & 6 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)$$

$\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2$. $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$의 결합은 3차원 공간을 채우지 **못한다**. $\mathbf{a}_3$는 3차원에서 $\mathbf{a}_1$과 $\mathbf{a}_2$의 평면 위에 놓인다.

**예시 (ex3):**

$$A_3 = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 6 & 8 \\ 5 & 15 & 20 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)$$

$\mathbf{a}_2 = 3\mathbf{a}_1$이고 $\mathbf{a}_3 = 4\mathbf{a}_1$.

$A_3$의 세 열 모두 3차원 공간에서 **같은 직선** 위에 놓인다.

### 열공간 $C(A)$

열공간 $C(A)$는 모든 벡터 $A\mathbf{x}$를 포함한다: 열들의 모든 결합.

**$A$의 열공간에 대한 고찰:**

$$A_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix} \quad \text{(독립 열)}$$

$$A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix} \quad \text{4번째 열은 종속인가?}$$

$A_5$에 대해: $\mathbf{a}_4 = \mathbf{a}_1 - \mathbf{a}_2 + \mathbf{a}_3$.

$\mathbf{v} = A_4\mathbf{x}$를 생각하자:

$$\mathbf{v} = x_1\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + x_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$$

$A_4\mathbf{x} = \mathbf{v}$ 풀기: 방정식으로부터,

$$v_4 = x_4, \quad v_3 = x_4 + x_3, \quad v_2 = x_4 + x_3 + x_2, \quad v_1 = x_4 + x_3 + x_2 + x_1$$

따라서: $x_4 = v_4$, $x_3 = v_3 - v_4$, $x_2 = v_2 - v_3$, $x_1 = v_1 - v_2$

$$\therefore \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \\ v_4 \end{pmatrix} = (v_1 - v_2)\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + (v_2 - v_3)\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + (v_3 - v_4)\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + v_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}$$

이는 **모든 $\mathbf{v}$가 열공간에 속한다**는 것을 의미한다. 네 개의 방정식 $A_4\mathbf{x} = \mathbf{v}$가 풀렸다.

### 1.3.4 생성 (Span)

**SPAN (생성)**은 벡터 집합의 모든 선형결합을 기술한다.

**$A$의 열들의 생성은 열공간이다.**

**예시:** $\mathbf{a}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}$, $\mathbf{a}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} \in \mathbb{R}^4$

- $c_1\mathbf{a}_1$은 4차원에서 하나의 직선에 해당한다.
- $c_2\mathbf{a}_2$는 4차원에서 하나의 직선에 해당한다.
- $c_1\mathbf{a}_1 + c_2\mathbf{a}_2$는 4차원 공간에서 2차원 평면을 채운다.
- 이 평면은 열 $\mathbf{a}_1$과 $\mathbf{a}_2$의 생성이다.

**예시:** $A_4$는 4개의 독립 열을 가진다. $A_4$의 열공간은 $\mathbb{R}^4$ 전체이다.

**예시:** $A_5$는 하나의 종속 열이 있다 (독립 열은 세 개뿐이다).

$$A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix}$$

$A_5$의 열공간은 $\mathbb{R}^4$ 내의 **3차원 부분공간**이다. 4번째 열은 그 부분공간에 속한다.

$A_5\mathbf{x} = \mathbf{v}$는 $\mathbf{v} \in C(A_5)$일 때만 풀 수 있다.

### 1.3.5 랭크와 기저 (Rank and Basis)

$A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}$로 놓자

$C(A)$는 $A$의 열공간, $\mathbf{a}_i \in \mathbb{R}^m$.

열공간 $C(A)$는 $\mathbb{R}^m$ 전체를 채울 수도 있고 그렇지 않을 수도 있다.

**예시:** $m = 3$인 경우. $A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}$

$C(A)$는:
- 전체 공간 $\mathbb{R}^3$ — $A$가 3개의 독립 열 (I.C.)을 가질 때
- $\mathbb{R}^3$에서의 평면 — $A$가 2개의 I.C.를 가질 때
- $\mathbb{R}^3$에서의 직선 — $A$가 1개의 I.C.를 가질 때
- 단일 점 $\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$ — $A$가 영행렬일 때

**더 많은 예시 ($3 \times 3$ 행렬):**

$$A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}: C(A) = \mathbb{R}^3$$

$$A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}: C(A) = xy \text{ 평면}$$

$$A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}: C(A) = x \text{ 축}$$

**핵심 질문:**

- **Q1:** $A$의 몇 개의 열이 독립인가? 그 수 $r$이 $A$의 **"랭크 (Rank)"**이다.
- **Q2:** 처음 $r$개의 독립 열은 무엇인가? 이것이 열공간의 **"기저 (Basis)"**이다.
- **Q3:** 그 $r$개의 기저 벡터의 어떤 결합이 나머지 $n - r$개의 열을 만드는가?
- **Q4:** 임의의 $A$를 $m \times r$ 열 행렬 $C$ 곱하기 $r \times n$ 행렬 $R$로 쓸 수 있다: $A = CR$.
- **Q5:** $R$의 $r$개의 행은 $A$의 **행공간 (Row space)**의 기저이다. $R$의 행은 $A$에서 직접 오지 않는다.

### 1.3.6 랭크 1인 행렬

랭크 1 행렬에서는 **모든 열벡터가 같은 직선을 따라 놓인다**.

**예시:**

$$A_6 = \begin{pmatrix} 1 & 3 & -2 \\ 4 & 12 & -8 \\ 2 & 6 & -4 \end{pmatrix}$$

$\mathbf{a}_2 = 3\mathbf{a}_1$, $\mathbf{a}_3 = -2\mathbf{a}_1$.

- $A_6$의 랭크 $r = 1$.
- $C(A_6) = c_1\mathbf{a}_1$, 원점을 지나는 직선.
- $A_6$의 모든 행은 하나의 행의 배수이다.

> 열공간이 $\mathbb{R}^m$에서의 단일 직선이면, 행공간은 $\mathbb{R}^n$에서의 단일 직선이다.

**질문: 왜 이렇게 되는가?**

모든 열이 같은 방향이면, 모든 행도 같은 방향이다.

**예시:** $A = \begin{pmatrix} a & ma \\ b & mb \end{pmatrix}$

행 2는 행 1의 배수이다. 열 랭크가 1이면 행 랭크도 1이다.

**예시:** $A = \begin{pmatrix} a & ma & pa \\ b & mb & pb \\ c & mc & pc \end{pmatrix}$

행 2 = $\frac{b}{a}$ 행 1, 행 3 = $\frac{c}{a}$ 행 1.

열 랭크가 1이면 행 랭크도 1이다.

> **질문. 행 랭크와 열 랭크는 모든 행렬에서 같은가.** 답: **그렇다.**

---

<br>

## 1.4 행렬 곱셈 AB와 CR

### 1.4.1 행렬 곱셈의 규칙

**(1)** $AB$를 곱하려면: $A \in \mathbb{R}^{m \times n}$, $B \in \mathbb{R}^{n \times p}$.

**$A$의 행 길이가 $B$의 열 길이와 같아야 한다.**

**(2) 내적 관점:**

$$(AB)_{ij} = (\text{$A$의 행 } i) \cdot (\text{$B$의 열 } j)$$

**(3) 열 관점:**

$$AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p)$$

$AB$의 열 $j$ = $A\mathbf{b}_j$

**(4) 비교환성:** 일반적으로 $AB \neq BA$.

**(5)** $A$가 $C$에서 $r$개의 독립 열을 가지면: $A = CR$ 여기서 $C \in \mathbb{R}^{m \times r}$, $R \in \mathbb{R}^{r \times n}$.

### 1.4.2 AB의 열 해석

$B = \begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix}$로 놓자

$$AB = A\begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix} = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}$$

이것은 **$A$의 결합**이다.

**예시:** $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$

$$AB = (A\mathbf{b}_1 \ A\mathbf{b}_2) = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}$$

$$A\mathbf{b}_1 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}$$

$$A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$$

**예시:** $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$, $B = \begin{pmatrix} 5 & 6 \\ 7 & 8 \end{pmatrix}$

$$A\mathbf{b}_1 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 5 \\ 7 \end{pmatrix}$$

**(1) 내적 접근법:**

$$= \begin{pmatrix} 1 \cdot 5 + 2 \cdot 7 \\ 3 \cdot 5 + 4 \cdot 7 \end{pmatrix} = \begin{pmatrix} 19 \\ 43 \end{pmatrix}$$

**(2) 열 결합 접근법:**

$$= 5\begin{pmatrix} 1 \\ 3 \end{pmatrix} + 7\begin{pmatrix} 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 19 \\ 43 \end{pmatrix}$$

두 접근법 모두 "4"번의 곱셈을 사용한다.

$$A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 6 \\ 8 \end{pmatrix} = \begin{pmatrix} 6 + 16 \\ 18 + 32 \end{pmatrix} = \begin{pmatrix} 22 \\ 50 \end{pmatrix}$$

$$\therefore AB = \begin{pmatrix} 19 & 22 \\ 43 & 50 \end{pmatrix}$$

### 1.4.3 계산 비용

$A \in \mathbb{R}^{m \times n}$, $B \in \mathbb{R}^{n \times p}$로 놓자.

$$AB = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}$$

**(1) 내적 횟수:** $A\mathbf{b}_1 \Rightarrow$ $m$개의 내적. $AB \Rightarrow$ $mp$개의 내적. 각 내적은 $n$번의 곱셈 $\Rightarrow$ 총 **$mnp$번의 곱셈**.

**(2) 열 결합 횟수:** $A\mathbf{b}_1 = b_{11}\mathbf{a}_1 + b_{12}\mathbf{a}_2 + \cdots + b_{1n}\mathbf{a}_n \Rightarrow$ 각각 $m$번의 곱셈, 하나의 열에 $mn$번의 곱셈. $AB \Rightarrow$ **$mnp$번의 곱셈**.

> 행렬 곱셈의 비용은 $mnp$이다.
>
> $A, B \in \mathbb{R}^{n \times n}$이면, $AB$의 비용은 $n^3$이다.

### 1.4.4 행렬 곱셈의 성질

**비교환성:**

$$AB \neq BA \quad \text{일반적으로.}$$

행렬 곱셈은 교환법칙이 성립하지 **않는다**.

**결합법칙:**

$$(AB)C = A(BC)$$

행렬 곱셈에서 결합법칙은 **성립한다**.

**분배법칙:**

$$A(B + C) = AB + AC$$

### 1.4.5 랭크 1 행렬과 A = CR

**$AB$ 복습:**

**(1) 내적:** $(AB)_{ij} = (\text{$A$의 행 } i) \cdot (\text{$B$의 열 } j)$. $(m \times n)(n \times p) \Rightarrow mp$개의 내적.

**(2) 열 결합:** $AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p) = (A\mathbf{b}_1 \ A\mathbf{b}_2 \ \cdots \ A\mathbf{b}_p)$. $AB$의 열 $j$ = $A\mathbf{b}_j$.

**랭크 1 행렬과 $A = CR$:**

랭크 1 행렬의 모든 열은 같은 직선 위에 놓인다.

$A$의 모든 열이 같은 열 방향에 있으면, $A$의 모든 행도 같은 행 방향에 있다.

**예시:** 랭크 $r = 1$.

$$A = \begin{pmatrix} 1 & 2 & 10 & 100 \\ 3 & 6 & 30 & 300 \\ 2 & 4 & 20 & 200 \end{pmatrix}$$

하나의 독립 열, 하나의 독립 행.

$A$의 열공간이 직선이면, $A$의 행공간도 직선이다.

$A$를 다음과 같이 분해할 수 있다:

$$A = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}\begin{pmatrix} 1 & 2 & 10 & 100 \end{pmatrix} = CR$$

여기서 $C \in \mathbb{R}^{3 \times 1}$이고 $R \in \mathbb{R}^{1 \times 4}$.

### 1.4.6 C와 R 구하기

**$C$는 $A$의 처음 $r$개의 독립 열을 포함한다.**

주어진 $A$에서 **왼쪽에서 오른쪽으로** 독립 열을 찾는다:

$$A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}$$

1. $\mathbf{a}_1 \neq \mathbf{0}$이면, $\mathbf{a}_1$을 $C$에 넣는다.
2. $\mathbf{a}_2 \neq c_1\mathbf{a}_1$ ($\mathbf{a}_1$의 배수가 아니면), $\mathbf{a}_2$를 $C$에 넣는다.
3. $\mathbf{a}_3 \neq c_1\mathbf{a}_1 + c_2\mathbf{a}_2$ ($\mathbf{a}_1$과 $\mathbf{a}_2$의 결합이 아니면), $\mathbf{a}_3$를 $C$에 넣는다.
4. $C$가 $r$개의 열을 가질 때까지 계속한다.

**수 $r$이 $A$와 $C$의 랭크이다.**

$\Rightarrow C\mathbf{x} = \mathbf{0}$이면 $\mathbf{x} = \mathbf{0}$일 때만 성립. 열들의 어떤 결합도 영벡터를 만들지 않는다.

**예시:**

$$A = \begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)$$

1. $C = \begin{pmatrix} 2 \\ 4 \\ 1 \end{pmatrix}$
2. $\mathbf{a}_2 = 3\mathbf{a}_1$
3. $\mathbf{a}_3 \neq c_1\mathbf{a}_1$이므로, $C = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}$

$\therefore$ 랭크 $r = 2$.

**질문. $A = CR$에서 $R$은 무엇인가?**

$$\begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 3 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$3 \times 3 = (3 \times 2)(2 \times 3)$$

열을 재배열하면:

$$\begin{pmatrix} 2 & 4 & 6 \\ 4 & 8 & 12 \\ 1 & 5 & 3 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 0 \end{pmatrix}$$

**관찰:**

1. $C$는 $A$의 $r$개의 독립 열의 완전한 집합을 포함한다.
2. $R = (I \ F)$는 항등행렬 $I \in \mathbb{R}^{r \times r}$을 포함한다.
3. $A$의 종속 열은 $C$에 있는 독립 열들의 결합이다.
4. $A = CR = C(I \ F) = (C \ CF)$
5. $C$는 $A$와 같은 열공간을 가진다. $R$은 $A$와 같은 행공간을 가진다.

**예시 (ex9):**

$$\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 4 & 5 \\ 7 & 8 \end{pmatrix}\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}$$

랭크 $r = 2$ 행렬.

$\mathbf{a}_j = C\mathbf{r}_j$

$(1 \ 2 \ 3) = (1 \ 2)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}$

$(4 \ 5 \ 6) = (4 \ 5)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}$

$(7 \ 8 \ 9) = (7 \ 8)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}$

$A$의 행 $i$ = $C$의 행 $i$ 곱하기 $R$: $R$의 행들의 결합.

**예시:** $A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 1 & 2 & 4 & 5 \end{pmatrix}_{2 \times 4}$

$$= \begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}\begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}$$

$C \in \mathbb{R}^{2 \times 2}$, $R \in \mathbb{R}^{2 \times 4}$. 랭크 $r = 2$이며, $R$을 사용하여 $C$로부터 $A$의 모든 열을 복원한다.

### R 행렬을 구하는 방법

제3장에서 다룰 **"소거법 (Elimination)"**을 사용할 수 있다.

**예시:**

$$A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}$$

**단계 1:** $R2 - 2R1$, $R3 - 3R1$:

$$\begin{pmatrix} 1 & 3 & 4 \\ 0 & -2 & -6 \\ 0 & -2 & -6 \end{pmatrix}$$

**단계 2:** $R3 - R2$:

$$\begin{pmatrix} 1 & 3 & 4 \\ 0 & -2 & -6 \\ 0 & 0 & 0 \end{pmatrix} \quad \text{행 랭크 } r = 2$$

**단계 3:** $R2 / (-2)$:

$$\begin{pmatrix} 1 & 3 & 4 \\ 0 & 1 & 3 \\ 0 & 0 & 0 \end{pmatrix}$$

**단계 4:** $R1 - 3R2$:

$$\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \\ 0 & 0 & 0 \end{pmatrix} \quad \leftarrow R \text{ 행렬}$$

$\mathbf{a}_1$과 $\mathbf{a}_2$의 열은 독립이다.

$$A = \begin{pmatrix} 1 & 3 \\ 2 & 4 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \end{pmatrix} = CR$$

### A = CR의 핵심 성질

- $A = CR$
- $C$의 $r$개의 열은 $A$의 열공간의 **기저**이다.
- $R$의 $r$개의 행은 $A$의 행공간의 **기저**이다.
- $\Rightarrow$ $r$ 차원.

> $A = CR$에서, "$R$"은 **기약행사다리꼴 (Reduced Row Echelon Form)**이다.

**예시 (ex5):**

$$A = \begin{pmatrix} | & | & | \\ \mathbf{a}_1 & \mathbf{a}_2 & 3\mathbf{a}_1 + 4\mathbf{a}_2 \\ | & | & | \end{pmatrix}_{m \times 3}$$

$$= \begin{pmatrix} | & | \\ \mathbf{a}_1 & \mathbf{a}_2 \\ | & | \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 4 \end{pmatrix}$$

$$= \mathbf{a}_1(1 \ 0 \ 3) + \mathbf{a}_2(0 \ 1 \ 4) = (\mathbf{a}_1 \ \mathbf{a}_2 \ 3\mathbf{a}_1 + 4\mathbf{a}_2) = A$$

### 1.4.7 A의 열 곱하기 B의 행 (외적, Outer Product)

$$AB = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}_{m \times n} \begin{pmatrix} - & \mathbf{b}_1^* & - \\ - & \mathbf{b}_2^* & - \\ & \vdots & \\ - & \mathbf{b}_n^* & - \end{pmatrix}_{n \times p}$$

$$= \mathbf{a}_1\mathbf{b}_1^* + \mathbf{a}_2\mathbf{b}_2^* + \cdots + \mathbf{a}_n\mathbf{b}_n^*$$

여기서 $\mathbf{a}_k\mathbf{b}_k^* = \begin{pmatrix} | \\ \mathbf{a}_k \\ | \end{pmatrix}_{m \times 1}(- \ \mathbf{b}_k^* \ -)_{1 \times p}$

이것은 **열 곱하기 행 = 랭크 1 행렬**이며, $mp$개의 원소를 가진다.

$$AB = \sum_{k=1}^{n} \mathbf{a}_k\mathbf{b}_k^*$$

$n$개의 랭크 1 행렬의 합. 각각 $mp$번의 곱셈 $\Rightarrow$ 총 $mnp$번의 곱셈.

**참고:** 행렬 $A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}$에 대해, $\text{rank}(A) = 2$.

$\Rightarrow$ $A$는 역행렬이 없다.
$\Rightarrow$ $A$의 행렬식은 0이다.

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:--------|:---------|
| **선형결합 (Linear Combination)** | $c\mathbf{v} + d\mathbf{w}$: 스칼라 곱셈 + 벡터 덧셈 |
| **열공간 (Column Space) $C(A)$** | 모든 벡터 $A\mathbf{x}$ = $A$의 열들의 모든 결합 |
| **생성 (Span)** | 벡터 집합의 모든 선형결합 |
| **독립벡터 (Independent Vectors)** | $A\mathbf{x} = \mathbf{0}$이 $\mathbf{x} = \mathbf{0}$일 때만 성립 |
| **종속벡터 (Dependent Vector)** | 다른 벡터들의 결합으로 쓸 수 있음 |
| **랭크 (Rank)** | 독립 열의 수 (= 독립 행의 수) |
| **기저 (Basis)** | 처음 $r$개의 독립 열이 $C(A)$의 기저를 형성 |
| **내적 (Dot Product)** | $\mathbf{v} \cdot \mathbf{w} = \sum v_i w_i$; 길이와 각도 정보 제공 |
| **길이 (노름, Norm)** | $\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}}$ |
| **단위벡터 (Unit Vector)** | $\|\mathbf{v}\| = 1$; 임의의 영이 아닌 벡터의 정규화: $\mathbf{v}/\|\mathbf{v}\|$ |
| **직교성 (Perpendicularity)** | $\mathbf{v} \perp \mathbf{w} \Leftrightarrow \mathbf{v} \cdot \mathbf{w} = 0$ |
| **각도 공식 (Angle Formula)** | $\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|}$ |
| **슈바르츠 부등식 (Schwarz Inequality)** | $|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|$ |
| **삼각 부등식 (Triangle Inequality)** | $\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|$ |
| **피타고라스 정리 (Pythagorean Theorem)** | $\mathbf{v} \perp \mathbf{w}$이면: $\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2$ |
| **행렬 곱셈 $AB$** | $(AB)_{ij} = \text{row}_i(A) \cdot \text{col}_j(B)$; $AB$의 열 $j = A\mathbf{b}_j$ |
| **외적 (Outer Product)** | $AB = \sum \mathbf{a}_k \mathbf{b}_k^*$ (랭크 1 행렬의 합) |
| **$AB$의 비용** | $A \in \mathbb{R}^{m \times n}$, $B \in \mathbb{R}^{n \times p}$에 대해 $mnp$번의 곱셈 |
| **비교환성 (Non-commutativity)** | 일반적으로 $AB \neq BA$ |
| **결합법칙 (Associativity)** | $(AB)C = A(BC)$ |
| **분배법칙 (Distributivity)** | $A(B+C) = AB + AC$ |
| **$A = CR$ 분해** | $C$: $A$의 독립 열; $R$: 기약행사다리꼴 |
| **랭크 1 행렬** | 모든 열이 한 직선 위; $A = \mathbf{c}\mathbf{r}^T$ |
| **행 랭크 = 열 랭크** | 모든 행렬에 대해 항상 성립 |
| **항등행렬 (Identity Matrix)** | $\mathbf{v} = I\mathbf{v}$; 표준 기저 벡터가 열 |
| **법선벡터에 의한 평면** | $\mathbf{w} \cdot \mathbf{n} = 0$은 $\mathbf{n}$에 수직인 평면을 정의 |
| **가역성 (Invertibility)** | $A$ 가역 $\Leftrightarrow$ 열이 독립 $\Leftrightarrow$ $A\mathbf{x} = \mathbf{b}$의 유일한 해 |
| **소거법 실패 (Elimination Failure)** | $\mathbf{v} \| \mathbf{w}$일 때: 해 없음 또는 무한히 많은 해 |

---
