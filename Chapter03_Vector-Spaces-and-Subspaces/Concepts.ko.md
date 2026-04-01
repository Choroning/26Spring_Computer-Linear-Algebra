# 제3장 강의 — 벡터 공간과 부분공간

> **최종 수정일:** 2026-03-31
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 3

> **선수 지식**: [선형대수학] 소거법과 행렬 연산 (제1-2장).
>
> **학습 목표**:
> 1. 벡터 공간과 부분공간을 정의하고 닫힘 성질을 검증할 수 있다
> 2. 열 공간, 영 공간, 행 공간, 좌영 공간을 구할 수 있다
> 3. 기저, 차원, 행렬의 랭크를 결정할 수 있다

---

<br>

## 목차

- [1. 개요: 네 가지 기본 부분공간](#1-개요-네-가지-기본-부분공간)
- [2. 벡터 공간과 부분공간 (3.1)](#2-벡터-공간과-부분공간-31)
  - [2.1 벡터 공간의 요건](#21-벡터-공간의-요건)
  - [2.2 공간 R^n](#22-공간-rn)
  - [2.3 벡터 공간의 정의](#23-벡터-공간의-정의)
  - [2.4 여덟 가지 공리](#24-여덟-가지-공리)
  - [2.5 공리의 결과](#25-공리의-결과)
  - [2.6 예제: 벡터 공간인가?](#26-예제-벡터-공간인가)
  - [2.7 벡터 공간의 예시](#27-벡터-공간의-예시)
  - [2.8 일반화된 벡터 공간](#28-일반화된-벡터-공간)
  - [2.9 벡터 공간의 부분공간](#29-벡터-공간의-부분공간)
  - [2.10 부분공간의 정의](#210-부분공간의-정의)
  - [2.11 부분공간과 비부분공간의 예시](#211-부분공간과-비부분공간의-예시)
  - [2.12 열공간과 행공간](#212-열공간과-행공간)
  - [2.13 생성 (Spanning)](#213-생성-spanning)
- [3. 소거법을 이용한 영공간 계산 (3.2)](#3-소거법을-이용한-영공간-계산-32)
  - [3.1 핵심 사항: A = CR](#31-핵심-사항-a--cr)
  - [3.2 Ax = 0의 모든 해 구하기](#32-ax--0의-모든-해-구하기)
  - [3.3 기약행사다리꼴](#33-기약행사다리꼴)
  - [3.4 특수해와 영공간의 기저](#34-특수해와-영공간의-기저)
  - [3.5 영공간 행렬: (-F; I)의 열들](#35-영공간-행렬--fi의-열들)
  - [3.6 행렬 분해 A = CR과 N(A)](#36-행렬-분해-a--cr과-na)
- [4. Ax = b의 완전해 (3.3)](#4-ax--b의-완전해-33)
  - [4.1 완전해의 구조](#41-완전해의-구조)
  - [4.2 풀이 예제: 특수해 구하기](#42-풀이-예제-특수해-구하기)
  - [4.3 완전해 분해](#43-완전해-분해)
  - [4.4 열 완전 계수: r = n](#44-열-완전-계수-r--n)
  - [4.5 풀림 조건](#45-풀림-조건)
  - [4.6 행 완전 계수와 완전해](#46-행-완전-계수와-완전해)
  - [4.7 선형 방정식의 네 가지 경우](#47-선형-방정식의-네-가지-경우)
- [5. 독립, 기저, 차원 (3.4)](#5-독립-기저-차원-34)
  - [5.1 독립 벡터](#51-독립-벡터)
  - [5.2 영공간을 통한 선형 독립](#52-영공간을-통한-선형-독립)
  - [5.3 부분공간을 생성하는 벡터](#53-부분공간을-생성하는-벡터)
  - [5.4 벡터 공간의 기저](#54-벡터-공간의-기저)
  - [5.5 벡터 공간의 차원](#55-벡터-공간의-차원)
  - [5.6 행렬 공간과 함수 공간의 기저](#56-행렬-공간과-함수-공간의-기저)
- [6. 네 부분공간의 차원 (3.5)](#6-네-부분공간의-차원-35)
  - [6.1 차원 요약](#61-차원-요약)
  - [6.2 부분공간의 직교성](#62-부분공간의-직교성)
  - [6.3 R_0의 네 부분공간](#63-r_0의-네-부분공간)
  - [6.4 A와 R_0의 관계](#64-a와-r_0의-관계)
  - [6.5 선형대수학의 기본 정리](#65-선형대수학의-기본-정리)
- [요약](#요약)

---

<br>

## 1. 개요: 네 가지 기본 부분공간

제3장은 다섯 가지 주요 주제를 다룬다:

- **3.1** 벡터 공간과 부분공간 (Vector Spaces and Subspaces) — 벡터 공간을 어떻게 정의하는가? 핵심 연산은 $\mathbf{u} + \mathbf{v}$와 $c\mathbf{u}$이다. 벡터 $\mathbf{u}$와 스칼라 $c$가 만족해야 하는 8가지 규칙이 있다.
- **3.2** $A$의 영공간 (Null Space): $A\mathbf{x} = \mathbf{0}$ 풀기
- **3.3** $A\mathbf{x} = \mathbf{b}$의 완전해 (Complete Solution)
- **3.4** 독립 (Independence), **기저 (Basis)**, **차원 (Dimension)** — 공간을 설명하는 벡터 집합. $A \in \mathbb{R}^{n \times n}$이라 하자. $A$가 $r$개의 독립 열을 가지면 $C(A)$의 차원은 $r$이다. $A\mathbf{x} = \mathbf{0}$의 $(n - r)$개 특수해는 $A$의 영공간 $N(A)$의 기저가 된다.
- **3.5** 네 부분공간의 차원 (Dimensions of the Four Subspaces):

| 부분공간 | 차원 |
|:---------|:----------|
| $A$의 열공간 (Column Space) | $r$ |
| $A$의 행공간 (Row Space) | $r$ |
| $A$의 영공간 (Null Space) | $n - r$ |
| $A^T$의 영공간 (좌영공간, Left Nullspace) | $m - r$ |

---

<br>

## 2. 벡터 공간과 부분공간 (3.1)

### 2.1 벡터 공간의 요건

1. **요건:** 모든 선형 결합 $c\mathbf{u} + d\mathbf{w}$가 벡터 공간 안에 머물러야 한다.
2. $A$의 **행공간 (Row Space)** 은 $A$의 행들에 의해 "생성"된다. $A$의 열들은 열공간 $C(A)$를 생성한다.
3. **행렬 (Matrix)** $M_1$부터 $M_n$까지와 **함수 (Function)** $f_1$부터 $f_n$까지가 **행렬 공간 (Matrix Space)** 과 **함수 공간 (Function Space)** 을 생성한다.

### 2.2 공간 R^n

공간 $\mathbb{R}^n$은 길이 $n$인 모든 열벡터 $\mathbf{v}$를 포함한다.

$$
\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} \in \mathbb{R}^n \quad \text{여기서 } x_1, x_2, \ldots, x_n \in \mathbb{R} \text{ (실수)}
$$

예시:
- $x \in \mathbb{R}^1$
- $\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} \in \mathbb{R}^2$
- $\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} \in \mathbb{R}^3$

**참고:** $x_1, x_2, \ldots, x_n \in \mathbb{C}$ (복소수)이면, 공간은 $\mathbb{C}^n$이 된다.

**참고:** $c, d \in \mathbb{R}$이고 $\mathbf{u}, \mathbf{w} \in \mathbb{R}^n$이면, $c\mathbf{u} + d\mathbf{w} \in \mathbb{R}^n$이다. 선형 결합은 벡터 공간 $\mathbb{R}^n$ 안에 머문다.

### 2.3 벡터 공간의 정의

**벡터 공간 (Vector Space)** (= 선형 공간, Linear Space) $V$는 원소(벡터)들이 서로 더해지고 수와 곱해질 수 있는 집합이다.

$\mathbf{u}, \mathbf{w} \in V$ (벡터 공간)일 때:
1. $\mathbf{u} + \mathbf{w} \in V$ — 벡터 덧셈 (Vector Addition)
2. $\alpha \mathbf{u} \in V \quad \forall \alpha \in \mathbb{F}$ — 스칼라 곱 (Scalar Multiplication)

**체 (Field)** 는 덧셈, 뺄셈, 곱셈, 나눗셈이 정의된 집합이다.

체의 예: $\mathbb{R}$ (실수), $\mathbb{C}$ (복소수).

### 2.4 여덟 가지 공리

$\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$이고 $c, d \in \mathbb{F}$일 때. 벡터 공간은 다음 여덟 가지 공리를 만족해야 한다:

**(1) 벡터 덧셈의 결합법칙 (Associativity):**

$$\mathbf{u} + (\mathbf{v} + \mathbf{w}) = (\mathbf{u} + \mathbf{v}) + \mathbf{w}$$

**(2) 벡터 덧셈의 교환법칙 (Commutativity):**

$$\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}$$

**(3) 벡터 덧셈의 항등원 (영벡터, Zero Vector):**

$$\exists!\; \mathbf{0} \in V \text{ s.t. } \mathbf{v} + \mathbf{0} = \mathbf{v} \quad \forall \mathbf{v} \in V$$

**(4) 벡터 덧셈의 역원 (Inverse Element):**

$$\forall \mathbf{v} \in V, \; \exists!\; -\mathbf{v} \in V \text{ s.t. } \mathbf{v} + (-\mathbf{v}) = \mathbf{0}$$

($-\mathbf{v}$는 $\mathbf{v}$의 덧셈 역원 (Additive Inverse))

**(5) 스칼라 곱과 체 곱의 양립성 (Compatibility):**

$$c(d\mathbf{v}) = (cd)\mathbf{v}$$

**(6) 스칼라 곱의 항등원 (Identity Element):**

$$1 \cdot \mathbf{v} = \mathbf{v}$$

**(7) 벡터 덧셈에 대한 스칼라 곱의 분배법칙 (Distributivity):**

$$c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}$$

**(8) 체 덧셈에 대한 스칼라 곱의 분배법칙 (Distributivity):**

$$(c + d)\mathbf{u} = c\mathbf{u} + d\mathbf{u}$$

### 2.5 공리의 결과

여덟 가지 공리로부터 다음 성질이 성립한다:

$$0\mathbf{u} = \mathbf{0}$$

$$c\mathbf{0} = \mathbf{0}$$

$$(-1)\mathbf{u} = -\mathbf{u}$$

$$c\mathbf{v} = \mathbf{0} \implies c = 0 \text{ 또는 } \mathbf{v} = \mathbf{0}$$

### 2.6 예제: 벡터 공간인가?

**Q.** 모든 양수 벡터의 집합 $\mathbb{X}$, 즉 $\mathbb{X} \ni \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}$이고 모든 $v_i > 0$인 경우, 벡터 공간인가?

**A.** 아니다. $-\mathbf{v} \notin \mathbb{X}$.

---

**Q.** $\mathbb{X}$를 $A\mathbf{x} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$의 해 집합이라 하자. $\mathbb{X}$는 벡터 공간인가?

**A.** 아니다. $\mathbf{u}, \mathbf{w} \in \mathbb{X}$이면, $A\mathbf{u} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$이고 $A\mathbf{w} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$이다. 하지만 $\mathbf{u} + \mathbf{w} \notin \mathbb{X}$인데, $A(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w} = \begin{pmatrix} 2 \\ \vdots \\ 2 \end{pmatrix} \neq \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}$이기 때문이다.

---

**Q.** $\mathbb{R}^n$에서 직선은 벡터 공간인가?

직선은 점들의 모임이다. $\mathbf{q}$는 직선 위의 위치를 나타낸다:

$$\mathbf{q} = \mathbf{p} + t\mathbf{d}, \quad t \in \mathbb{R}$$

$\mathbf{p}$와 $\mathbf{d}$는 고정되어 있고, $t$는 $-\infty$에서 $\infty$까지 변한다.

$\mathbb{X}$를 모든 벡터 $\mathbf{q}$의 집합이라 하자: $\mathbb{X} = \{\mathbf{q} : \mathbf{p} + t\mathbf{d} \;\forall\; t \in \mathbb{R}\}$.

$\mathbf{a}, \mathbf{b} \in \mathbb{X}$를 취하자. $\mathbf{a} + \mathbf{b}$가 $\mathbb{X}$에 속하는가?

$$\mathbf{a} = \mathbf{p} + t_1 \mathbf{d}$$
$$\mathbf{b} = \mathbf{p} + t_2 \mathbf{d}$$
$$\mathbf{a} + \mathbf{b} = 2\mathbf{p} + (t_1 + t_2)\mathbf{d} \notin \mathbb{X} \quad \text{만약 } \mathbf{p} \neq \mathbf{0}$$

그러나:
- $\mathbf{a} + \mathbf{b} = (t_1 + t_2)\mathbf{d} \in \mathbb{X}$ (만약 $\mathbf{p} = \mathbf{0}$)
- $c\mathbf{a} = ct_1\mathbf{d} \in \mathbb{X}$ (만약 $\mathbf{p} = \mathbf{0}$)

**원점을 지나는 직선은 벡터 공간이다.** $\mathbb{R}^n$에서 $\mathbf{0}$을 지나는 직선은 $\mathbb{R}^n$의 **부분공간 (Subspace)** 이다 — 다른 벡터 공간 안의 벡터 공간.

### 2.7 벡터 공간의 예시

- $\mathbb{R}^n$은 벡터 공간이다.
- $\mathbb{Z} = \{\mathbf{0}\}$은 벡터 공간이다 (가장 작은 벡터 공간):
  - i) $\mathbf{0} + \mathbf{0} = \mathbf{0} \in \mathbb{Z}$
  - ii) $c\mathbf{0} = \mathbf{0} \in \mathbb{Z}$

$\mathbb{Z}$는 **가역 행렬의 영공간 (Null Space)** 이다. $A\mathbf{x} = \mathbf{0}$의 유일한 해가 영벡터 $\mathbf{x} = \mathbf{0}$이면, $A$의 열들은 선형 독립 (LI)이고 $A$의 영공간은 $\mathbb{Z}$이다.

### 2.8 일반화된 벡터 공간

**참고:** 벡터 공간 (= 선형 공간)은 원소들이 서로 더해지고 수와 곱해질 수 있는 집합이다.

즉, $\mathbf{u}, \mathbf{w} \in V \implies c\mathbf{u} + d\mathbf{w} \in V$.

**행렬과 함수도 벡터로 간주할 수 있다.**

**행렬의 벡터 공간 (Vector Space of Matrices):**

$A, B \in \mathbb{R}^{m \times n} \implies cA + dB \in \mathbb{R}^{m \times n}$

고정된 크기의 모든 행렬의 집합은 벡터 공간을 형성한다. (여덟 가지 규칙을 만족하는지 확인하라.)

**함수의 벡터 공간 (Vector Space of Functions):**

$\mathbb{F}$를 $\mathbb{R}$에서 원소를 취하여 실수로 사상하는 함수들의 집합이라 하자:

$$\mathbb{F} = \{f \mid f: \mathbb{R} \to \mathbb{R}\}$$

$f, g \in \mathbb{F}$이고 $c, d \in \mathbb{R}$이면, $cf + dg \in \mathbb{F}$이다.

정의:
- $(f + g)(x) = f(x) + g(x)$
- $(cf)(x) = c(f(x))$

**예:** $\mathbb{F}$ = 함수 $y = ce^x$의 직선.

$\mathbb{F} = \{f \mid f: \mathbb{R} \ni x \to ce^x \in \mathbb{R}, \;\forall c \in \mathbb{R}\}$

$f(x) = e^x$, $g(x) = 2e^x$.

$(f + g)(x) = f(x) + g(x) = e^x + 2e^x = 3e^x \in \mathbb{F}$

$(cf)(x) = c(f(x)) = ce^x \in \mathbb{F}$

**비고:** **집합 (Set)** 은 추가적인 구조 없이 원소들의 모임에 불과하다. **공간 (Space)** 은 추가적인 구조가 정의된 집합이다. 예: 벡터 공간은 벡터 덧셈과 스칼라 곱이 정의된 벡터들의 집합이다. 벡터 공간은 8가지 규칙을 만족해야 한다.

### 2.9 벡터 공간의 부분공간

벡터 공간 $\mathbb{R}^n$을 고려하자. 여기서 $\mathbf{v} \in \mathbb{R}^n$은 $n$개의 성분을 가진 열벡터이다.

$\mathbb{R}^n$ 안에 중요한 벡터 공간들이 있다. 이들이 $\mathbb{R}^n$의 **부분공간 (Subspace)** 이다.

평면은 벡터 공간이다. $\mathbf{v}, \mathbf{w} \in \mathbb{R}^2$이면, $c\mathbf{v} + d\mathbf{w} \in \mathbb{R}^2$이다.

$\mathbb{R}^3$에서 원점 $\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$을 지나는 평면은 벡터 공간이다. 이 평면은 $\mathbb{R}^2$가 아닌데, $\mathbf{u}, \mathbf{v} \in \mathbb{R}^3$이기 때문이다. 이 평면은 $\mathbb{R}^3$ 안의 벡터 공간이다. 이것은 전체 벡터 공간 $\mathbb{R}^3$의 **부분공간** 이다.

### 2.10 부분공간의 정의

**정의.** 벡터 공간의 **부분공간 (Subspace)** 은 다음을 만족하는 벡터들의 집합($\mathbf{0}$ 포함)이다:

> i) $\mathbf{u} + \mathbf{w}$가 부분공간에 속한다
>
> ii) $c\mathbf{u}$가 부분공간에 속한다

$\mathbf{u}, \mathbf{w} \in$ 부분공간이고 $c \in \mathbb{R}$ (또는 $\mathbb{F}$)일 때.

벡터 집합이 덧셈 $\mathbf{u} + \mathbf{w}$과 곱셈 $c\mathbf{u}$에 대해 "**닫혀 (Closed)**" 있다.

조건 i)과 ii)는 다음을 의미한다: **모든 선형 결합이 부분공간 안에 머문다.**

**참고:** 모든 부분공간은 영벡터를 포함한다. ii)에서 $c = 0$으로 놓으면 $0\mathbf{u} = \mathbf{0}$이 부분공간에 속한다.

### 2.11 부분공간과 비부분공간의 예시

**Q:** $\mathbb{R}^3$에서 평면 $z = 5$는 부분공간인가? **아니다.** (원점을 포함하지 않는다.)

원점을 지나는 직선은 부분공간이다.

$\mathbb{R}^3$은 자기 자신의 부분공간이다.

단일 벡터 $\mathbb{Z} = \{\mathbf{0}\}$은 $\mathbb{R}^3$의 부분공간이다.

---

**예 1.** $\mathbb{R}^2$는 벡터 공간이다. 제1사분면은 부분공간인가?

$$\mathcal{U} = \left\{\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0, y \geq 0 \right\}$$

$\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}$, $c = -1$을 취하자.

그러면 $c\mathbf{u} = \begin{pmatrix} -2 \\ -3 \end{pmatrix} \notin \mathcal{U}$.

이는 규칙 ii)를 위반한다. $\mathcal{U}$는 $\mathbb{R}^2$의 부분공간이 **아니다**.

---

**예 2.** $\mathcal{U} = \left\{\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0 \text{ 이고 } y \geq 0, \text{ 또는 } x \leq 0 \text{ 이고 } y \leq 0 \right\}$

$\mathcal{U}$는 부분공간인가?

$\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}$, $\mathbf{w} = \begin{pmatrix} -3 \\ -2 \end{pmatrix} \in \mathcal{U}$를 취하자.

하지만 합 $\mathbf{u} + \mathbf{w} = \begin{pmatrix} -1 \\ 1 \end{pmatrix} \notin \mathcal{U}$.

$\mathcal{U}$는 부분공간이 **아니다**.

---

**예 3.** $\mathbb{M}$은 $2 \times 2$ 행렬 $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$의 벡터 공간이다.

$\mathcal{U}$는 모든 상삼각 행렬 (Upper Triangular Matrix) $\begin{pmatrix} a & b \\ 0 & d \end{pmatrix}$의 집합이다.

$\mathbb{D}$는 모든 대각 행렬 (Diagonal Matrix) $\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix}$의 집합이다.

**$\mathcal{U}$와 $\mathbb{D}$ 모두 $\mathbb{M}$의 부분공간이다:**

- $A, B \in \mathcal{U} \implies cA + dB \in \mathcal{U}$
- $A, B \in \mathbb{D} \implies cA + dB \in \mathbb{D}$

참고:

$$\begin{pmatrix} a & b \\ c & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + c\begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$$

### 2.12 열공간과 행공간

**$A$의 열공간 (Column Space):**

$$A\mathbf{x} = \mathbf{b}$$

$$\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}$$

$$= x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n$$

**$A$의 열공간** 은 열벡터들로 구성된 벡터 공간이다.

$A\mathbf{x} = \mathbf{b}$를 푸는 것은 $\mathbf{b}$를 열들의 선형 결합으로 표현하는 것이다. 우변 벡터 $\mathbf{b}$는 $A$의 열공간에 속해야 한다.

> **방정식 $A\mathbf{x} = \mathbf{b}$가 풀리려면 $\mathbf{b}$가 $A$의 열공간에 속해야 한다.**

**$A$의 행공간 (Row Space):**

$A$의 행들은 $A^T$의 열들이다.

$$A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} \implies A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}$$

$A^T$는 $m$개의 열벡터를 가진다.

**$A$의 행공간은 $A^T$의 열공간이다.**

방정식 $A^T \mathbf{y} = \mathbf{c}$가 풀리려면 $\mathbf{c}$가 $A^T$의 열공간에 속해야 한다 ($= \mathbf{c}$가 $A$의 행공간에 속해야 한다).

**예:** 계수 1 행렬 (Rank 1 Matrix) $A = \mathbf{u}\mathbf{v}^T$를 고려하자:

$$A = \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix} \begin{pmatrix} v_1 & v_2 & \cdots & v_n \end{pmatrix} = \begin{pmatrix} v_1\mathbf{u} & v_2\mathbf{u} & \cdots & v_n\mathbf{u} \end{pmatrix}$$

$C(A)$는 $A$의 모든 열벡터의 직선: $c\mathbf{u}$.

$A^T = \mathbf{v}\mathbf{u}^T = \begin{pmatrix} u_1\mathbf{v} & u_2\mathbf{v} & \cdots & u_m\mathbf{v} \end{pmatrix}$

$C(A^T)$는 $A^T$의 모든 열벡터의 직선: $c\mathbf{v}$.

### 2.13 생성 (Spanning)

**$A$의 열들이 벡터 공간 $C(A)$를 생성 (Span)한다.**

$\mathbb{S}$를 $\mathbb{R}^m$의 벡터 집합이라 하자. $\mathbb{S} = \left\{ \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix}, \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_m \end{pmatrix} \right\}$이면, $\mathbb{S}$는 $\mathbb{R}^m$의 부분공간이 **아니다**. 왜냐하면 $\mathbf{u} + \mathbf{v} \notin \mathbb{S}$이기 때문이다.

$\mathbb{S}$에 속한 벡터들의 모든 결합을 포함하면, 벡터 공간 $V$를 얻는다.

**집합 $\mathbb{S}$가 $V$를 생성한다.** $V$는 $\mathbb{S}$를 포함하는 가장 작은 벡터 공간이다.

행렬 $A \in \mathbb{R}^{m \times n}$을 고려하자:
- $n$개의 열이 열공간 $C(A)$를 생성한다
- $A^T$의 $m$개의 열이 행공간 $C(A^T)$를 생성한다

---

<br>

## 3. 소거법을 이용한 영공간 계산 (3.2)

### 3.1 핵심 사항: A = CR

$$A = CR$$

1. $\mathbb{R}^n$에서 영공간 $N(A)$는 $A\mathbf{x} = \mathbf{0}$의 모든 해 $\mathbf{x}$를 포함한다. 여기에는 $\mathbf{x} = \mathbf{0}$도 포함된다.
2. $A$에서 $R_0$, $R$로의 소거는 영공간을 **바꾸지 않는다**.
3. 기약행사다리꼴 (Reduced Row Echelon Form) $R_0 = \text{rref}(A)$는 $r$개의 열에 $I$를 갖고, $n - r$개의 열에 $F$를 갖는다.
4. 열 $j$가 이전 열들에 종속이면, $A\mathbf{x} = \mathbf{0}$은 $x_j = 1$인 "특수해 (Special Solution)"를 가진다.
5. $A\mathbf{x} = \mathbf{0}$의 $n - r$개의 특수해는 $-F$와 $I$를 포함한다.
6. $m < n$인 모든 짧고 넓은 행렬은 영공간에 $A\mathbf{x} = \mathbf{0}$의 영이 아닌 해를 가진다.

### 3.2 Ax = 0의 모든 해 구하기

$A\mathbf{x} = \mathbf{0}$의 모든 해를 구하고 싶다.

$A \in \mathbb{R}^{n \times n}$이 가역이면 ($\text{rank}(A) = n$), 유일한 해는 $\mathbf{x} = \mathbf{0}$이다. $A$의 영공간은 영벡터만 포함한다: $N(A) = \{\mathbf{0}\}$.

일반적으로는 $\text{rank}(A) = r$이다. 즉, $A$는 $r$개의 독립 열을 가진다. 나머지 $n - r$개의 종속 열은 독립 열들의 결합이다. $N(A)$에서 $n - r$개의 벡터를 찾을 것인데, 이들이 $A\mathbf{x} = \mathbf{0}$의 특수해이다.

제2장에서 가역 행렬 $A$는 상삼각 행렬 (Upper Triangular Matrix) $U$로 축소되었다. $A \in \mathbb{R}^{m \times n}$에 대해, $A\mathbf{x} = \mathbf{0}$은 $R\mathbf{x} = \mathbf{0}$ (사다리꼴, Echelon Form)으로 단순화된다.

### 3.3 기약행사다리꼴

이 절에서는 $A \in \mathbb{R}^{m \times n}$을 고려하고, $A$를 기약행사다리꼴 (Reduced Row Echelon Form)로 소거한다: $R_0 = \text{rref}(A)$.

$R_0$는 영행을 가질 수 있다. $R_0$의 모든 영행을 제거하면 $R$이 된다.

**예 1:** $R = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = (I \quad F)$

$\text{rank}(R) = 2$, $n = 4$, $n - r = 2$개의 종속 열.

$$R\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_1 + 3x_3 + 5x_4 \\ x_2 + 4x_3 + 6x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}$$

특수해는 어떻게 구하는가? 두 변수를 고정하고 $R\mathbf{x} = \mathbf{0}$을 이용한다:

$x_3 = 1, x_4 = 0$으로 취하면: $\Rightarrow x_1 + 3 = 0, \; x_2 + 4 = 0$. $\therefore \mathbf{s}_1 = \begin{pmatrix} -3 \\ -4 \\ 1 \\ 0 \end{pmatrix}$

$x_3 = 0, x_4 = 1$로 취하면: $\Rightarrow x_1 + 5 = 0, \; x_2 + 6 = 0$. $\therefore \mathbf{s}_2 = \begin{pmatrix} -5 \\ -6 \\ 0 \\ 1 \end{pmatrix}$

### 3.4 특수해와 영공간의 기저

두 특수해 $\mathbf{s}_1, \mathbf{s}_2$는 $R$의 영공간에 속한다:

$$\mathbf{s}_1, \mathbf{s}_2 \in N(R)$$

$$R\mathbf{s}_1 = \mathbf{0}, \quad R\mathbf{s}_2 = \mathbf{0}$$

$$\Rightarrow R(c_1 \mathbf{s}_1 + c_2 \mathbf{s}_2) = \mathbf{0}$$

$\mathbf{s}_1$과 $\mathbf{s}_2$의 선형 결합은 $R$의 영공간에 속한다. **특수해 $\mathbf{s}_1$과 $\mathbf{s}_2$는 영공간의 기저이다.**

---

**예 2:** $R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix}$ (영행이 있는 기약행사다리꼴)

$$R_0 \mathbf{x} = \begin{pmatrix} x_1 + 7x_2 + 8x_4 \\ x_3 + 9x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$\text{rank}(R_0) = 2$, $n = 4$, $n - r = 2$. 두 개의 특수해를 구하자.

$x_1 = 0, x_2 = 1$로 취하면:

$7 + 8x_4 = 0 \Rightarrow x_4 = -7/8$, $x_3 = -9(-7/8) = 63/8$

$\Rightarrow \mathbf{x} = \begin{pmatrix} 0 \\ 1 \\ 63/8 \\ -7/8 \end{pmatrix} \Rightarrow \mathbf{s}_1 = \begin{pmatrix} 0 \\ 8 \\ 63 \\ -7 \end{pmatrix}$

$x_1 = 1, x_2 = 0$으로 취하면:

$1 + 8x_4 = 0 \Rightarrow x_4 = -1/8$, $x_3 = +9/8$

$\Rightarrow \mathbf{x} = \begin{pmatrix} 1 \\ 0 \\ 9/8 \\ -1/8 \end{pmatrix} \Rightarrow \mathbf{s}_2 = \begin{pmatrix} 8 \\ 0 \\ 9 \\ -1 \end{pmatrix}$

$x_2 = 1, x_4 = 0$과 $x_2 = 0, x_4 = 1$로 취하여 $x_1$과 $x_3$을 구할 수도 있다.

### 3.5 영공간 행렬: (-F; I)의 열들

이 절에서는 $A \in \mathbb{R}^{m \times n}$을 고려하고, $A$를 기약행사다리꼴 $R_0 = \text{rref}(A)$로 소거한다.

$R_0$는 영행을 가질 수 있다. $R_0$의 모든 영행을 제거하면 $R$이 된다.

$$R = (I \quad F)Q$$

여기서 $Q$는 순열 행렬 (Permutation Matrix)이다 (피벗 열이 앞에 오도록 열을 재배치).

**상기** $R\mathbf{x} = \mathbf{0}$, 즉 $(I \quad F)\mathbf{x} = \mathbf{0}$ (예 1):

$$\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

자유 변수를 $\begin{pmatrix} 1 \\ 0 \end{pmatrix}$과 $\begin{pmatrix} 0 \\ 1 \end{pmatrix}$로 설정하면:

$$I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 3 \\ 4 \end{pmatrix}$$

$$I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 5 \\ 6 \end{pmatrix}$$

동치로:

$$\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} -3 & -5 \\ -4 & -6 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$$

즉:

$$\boxed{(I \quad F) \begin{pmatrix} -F \\ I \end{pmatrix} = O}$$

**$(I \quad F)\mathbf{x} = \mathbf{0}$의 두 특수해는 $\begin{pmatrix} -F \\ I \end{pmatrix}$의 열들이다.**

### 3.6 행렬 분해 A = CR과 N(A)

**$A$에서 $\text{rref}(A)$로의 소거: 기약행사다리꼴**

소거를 적용하여 $A$를 $R_0$로 축소한다. 그러면 $R_0$에서 "$I$"가 $A$의 독립 열들의 행렬 $C$를 찾아준다. $R_0$에서 영행을 제거하면 $A = CR$에서의 행 행렬 $R$이 된다.

제2장에서 $A$는 정사각이고 가역이었다:

- $A\mathbf{x} = \mathbf{b}$ $\xrightarrow{a) \text{ 소거}}$ $U\mathbf{x} = \mathbf{c}$ $\xrightarrow{b) \text{ 후진 대입}}$ $\mathbf{x} = U^{-1}\mathbf{c}$

계수 $r$인 임의의 행렬 $A$에 대해:

$$A = CR = C(I \quad F)$$

소거는 $I_{r \times r}$ 단위 행렬에 도달할 때까지 계속된다.

---

**예 1:**

$$A = \begin{pmatrix} 1 & 2 & 11 & 17 \\ 3 & 7 & 37 & 57 \end{pmatrix}$$

$\xrightarrow{R_2 - 3R_1}$ $\begin{pmatrix} 1 & 2 & 11 & 17 \\ 0 & 1 & 4 & 6 \end{pmatrix}$ $\xrightarrow{R_1 - 2R_2}$ $\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = R$

$\text{rank}(A) = 2$, $n = 4$, $n - r = 2$.

$A = (W \quad H)$이고, $W = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}$ (독립 열), $H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix}$ (종속 열).

$R = (I \quad F)$. 소거는 $W$를 역행렬로 만들었다. 이는 $A$에 $W^{-1}$을 곱하는 것과 같다:

$$W^{-1}A = W^{-1}(W \quad H) = (I \quad W^{-1}H) = (I \quad F) = R$$

$W^{-1}H = F \implies H = WF$.

$W$는 독립 열로 구성되어 있다. $F$는 $A$의 독립 열들을 어떻게 결합하는지 알려준다.

$$H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix} = WF = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 3 & 5 \\ 4 & 6 \end{pmatrix}$$

---

**예 2:**

$$A = \begin{pmatrix} 1 & 7 & 3 & 35 \\ 2 & 14 & 6 & 70 \\ 2 & 14 & 9 & 97 \end{pmatrix}$$

$\xrightarrow{R_2 - 2R_1, \; R_3 - 2R_1}$ $\begin{pmatrix} 1 & 7 & 3 & 35 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 3 & 27 \end{pmatrix}$ $\xrightarrow{R_1 - R_3, \; R_3/3}$ $\begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 9 \end{pmatrix}$ $\xrightarrow{\text{swap } R_2 \text{ and } R_3}$

$$R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q$$

여기서 $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$, $F = \begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}$이고, $Q$는 열을 재배치하는 순열 행렬이다.

영행을 제거하면: $R = (I \quad F)Q$.

열 순열 $PAQ$ 후:

$PAQ = \begin{pmatrix} 1 & 3 & 7 & 35 \\ 2 & 9 & 14 & 97 \\ 2 & 6 & 14 & 70 \end{pmatrix}$이고, $W = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}$, $H = \begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix}$

$W^{-1}H = F \implies H = WF$:

$$\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}$$

**예 2 상기:** 단위 행렬이 $R_0$의 1번째와 3번째 열에 나타나므로, $A$의 1번째와 3번째 열이 독립이다. $A$의 2번째와 4번째 열은 독립 열들의 선형 결합이다:

$$\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}$$

**행렬 분해 $A = CR$과 $N(A)$:**

소거를 적용하여 $A$를 $R_0$로 축소한다. 그러면 $R_0$에서 "$I$"가 $A$의 독립 열들의 행렬 $C$를 찾아준다. $R_0$에서 영행을 제거하면 $A = CR$에서의 행 행렬 $R$이 된다.

예 2의 두 특수해 $\mathbf{s}_1, \mathbf{s}_2$를 구하자:

$$R = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \end{pmatrix}$$

$$R\mathbf{s} = \begin{pmatrix} s_1 + 7s_2 + 8s_4 \\ s_3 + 9s_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$s_2 = 1, s_4 = 0$으로 취하면: $\Rightarrow s_3 = 0, s_1 = -7$. $\therefore \mathbf{s}_1 = \begin{pmatrix} -7 \\ 1 \\ 0 \\ 0 \end{pmatrix}$

$s_2 = 0, s_4 = 1$로 취하면: $\Rightarrow s_1 = -8, s_3 = -9$. $\therefore \mathbf{s}_2 = \begin{pmatrix} -8 \\ 0 \\ -9 \\ 1 \end{pmatrix}$

---

<br>

## 4. Ax = b의 완전해 (3.3)

### 4.1 완전해의 구조

1. **$A\mathbf{x} = \mathbf{b}$의 완전해 (Complete Solution):**

$$\mathbf{x} = \text{하나의 특수해 } \mathbf{x}_p + \text{영공간의 임의의 } \mathbf{x}_n$$

2. $A\mathbf{x} = \mathbf{b}$에 **소거** 를 적용하면 $R_0 \mathbf{x} = \mathbf{d}$가 된다. $R_0$의 영행에 대응하는 $\mathbf{d}$의 성분이 0일 때 풀 수 있다.

$$\begin{pmatrix} R \\ 0 \; 0 \; \cdots \; 0 \\ 0 \; 0 \; \cdots \; 0 \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ 0 \\ 0 \end{pmatrix}$$

3. $R_0 \mathbf{x} = \mathbf{d}$가 풀릴 때, 하나의 특수해 $\mathbf{x}_p$는 **모든 자유 변수가 0** 인 해이다.

4. $A$가 **열 완전 계수 (Full Column Rank)** $r = n$일 때 영공간 $N(A) = \{\text{영벡터}\}$. 자유 변수가 없다.

5. $A$가 **행 완전 계수 (Full Row Rank)** $r = m$일 때 열공간 $C(A)$가 $\mathbb{R}^m$이다: $A\mathbf{x} = \mathbf{b}$는 항상 풀린다.

### 4.2 풀이 예제: 특수해 구하기

이전 절에서 $A\mathbf{x} = \mathbf{0}$의 해를 구했다. 이 절에서는 $A\mathbf{x} = \mathbf{b}$의 해를 구한다.

좌변의 행 연산은 우변에도 적용해야 한다. **첨가 행렬 (Augmented Matrix)** $(A | \mathbf{b})$를 사용한다.

$$\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 1 & 3 & 1 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 7 \end{pmatrix} \quad \Rightarrow \quad (A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & 1 \\ 0 & 0 & 1 & 4 & | & 6 \\ 1 & 3 & 1 & 6 & | & 7 \end{pmatrix}$$

소거 후 ($R_3 - R_1 - R_2$):

$$\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} \quad \Rightarrow \quad (R_0 | \mathbf{d})$$

마지막 방정식은 $0 = 0$이다. 일반적인 $\mathbf{b}$를 고려하자:

$$(A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 1 & 3 & 1 & 6 & | & b_3 \end{pmatrix} \xrightarrow{R_3 - R_1 - R_2} \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 0 & 0 & 0 & 0 & | & b_3 - b_1 - b_2 \end{pmatrix} = (R_0 | \mathbf{d})$$

$b_3 - b_1 - b_2 = 0$이면 세 번째 방정식에서 $0 = 0$을 얻을 수 있다.

**하나의 특수해 $A\mathbf{x}_p = \mathbf{b}$:**

$$R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} = \mathbf{d}$$

$\text{rank}(R_0) = 2$, $n = 4$, $n - r = 2$개의 자유 변수.

**$x_2 = 1, x_4 = 0$으로 취하면:**

$x_1 + 3 = 1 \Rightarrow x_1 = -2$, $x_3 = 6$

$\therefore \mathbf{x} = \begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix}$

**$x_2 = 0, x_4 = 1$로 취하면:**

$x_1 + 2 = 1 \Rightarrow x_1 = -1$, $x_3 + 4 = 6 \Rightarrow x_3 = 2$

$\therefore \mathbf{x} = \begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix}$

**$x_2 = 0, x_4 = 2$로 취하면:**

$x_1 + 2 \cdot 2 = 1 \Rightarrow x_1 = -3$, $x_3 + 4 \cdot 2 = 6 \Rightarrow x_3 = -2$

$\therefore \mathbf{x} = \begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix}$

자유 변수로 인해 무한히 많은 해를 구할 수 있다.

### 4.3 완전해 분해

해는 **특수해 (Particular Solution)** 와 영공간의 **특수해들 (Special Solutions)** 로 분해할 수 있다.

$$R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$R_0 = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q$이고, $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$, $F = \begin{pmatrix} 3 & 2 \\ 0 & 4 \end{pmatrix}$

영공간 벡터는 $Q^T \begin{pmatrix} -F \\ I \end{pmatrix}$의 열들이다:

$$\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} -3 & -2 \\ 0 & -4 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} -3 & -2 \\ 1 & 0 \\ 0 & -4 \\ 0 & 1 \end{pmatrix}$$

**$R_0 \mathbf{x} = \mathbf{d}$의 해 재검증:**

$$\begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix}$$

$$\begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}$$

$$\begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -4 \\ 0 \\ -8 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + 2\begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}$$

**Q. 특수해는 어떻게 구하는가?**

$x_2 = 0, x_4 = 0$ (모든 자유 변수를 0으로 설정):

$$R_0 \mathbf{x} = \begin{pmatrix} x_1 + 3x_2 + 2x_4 \\ x_3 + 4x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix}$$

$\Rightarrow x_1 = 1, \; x_3 = 6$

$$\therefore \mathbf{x}_p = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix}$$

**$A\mathbf{x} = \mathbf{b}$의 완전해** $\mathbf{x}_p + \mathbf{x}_n$:

$$\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + x_2 \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_4 \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}$$

여기서 $\mathbf{x}_p$가 특수해이고 나머지 항들은 영공간 벡터이다.

### 4.4 열 완전 계수: r = n

**Q. $m = n = r$일 때 $\mathbf{x}_p, \mathbf{x}_n$은?**

**A.** $\mathbf{x}_n = \mathbf{0}$.

$A\mathbf{x} = \mathbf{b} \iff A(\mathbf{x}_p + \mathbf{x}_n) = \mathbf{b} \implies A(\mathbf{x}_p + \mathbf{0}) = \mathbf{b} \implies \mathbf{x}_p = A^{-1}\mathbf{b}$

### 4.5 풀림 조건

**예 1.** $A\mathbf{x} = \mathbf{b}$가 풀리기 위한 $(b_1, b_2, b_3)$의 조건을 구하라:

$$A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ -2 & -3 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \\ b_3 \end{pmatrix}$$

$$(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & | & b_1 \\ 1 & 2 & | & b_2 \\ -2 & -3 & | & b_3 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 + 2R_1} \begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & -1 & | & b_3 + 2b_1 \end{pmatrix} \xrightarrow{R_3 + R_2}$$

$$\begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & | & 2b_1 - b_2 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} = (R_0, \mathbf{d})$$

$b_3 + b_2 + b_1 = 0$이면, $A\mathbf{x} = \mathbf{b}$가 풀린다. 즉, $b_3 + b_2 + b_1 = 0$이 $\mathbf{b}$를 $A$의 열공간에 넣기 위한 조건이다.

$\text{rank}(A) = 2$, $n = 2$, $n - r = 0$. 자유 변수 없음. $N(A) = \{\mathbf{0}\}$, $\mathbf{x}_n = \mathbf{0}$.

$$\mathbf{x}_p = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix}$$

$$\therefore \mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$b_3 + b_2 + b_1 \neq 0$이면, $A\mathbf{x} = \mathbf{b}$의 해가 없다.

**열 완전 계수 $r = n$일 때:**

$$R_0 = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix}$$

1. 행렬 $A$는 $n$개의 독립 열을 가진다.
2. $A$의 영공간은 $\mathbb{Z} = \{\mathbf{0}\}$이다.
3. $A\mathbf{x} = \mathbf{b}$에 해가 있으면, **유일한** 해를 가진다.

$$R_0 \mathbf{x} = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ O_{(m-n) \times 1} \end{pmatrix}$$

하단 $m - n$개의 행이 $\mathbf{b}$가 $A$의 열공간에 속하기 위한 $m - n$개의 조건을 제공한다.

### 4.6 행 완전 계수와 완전해

$A \in \mathbb{R}^{m \times n}$, $m \leq n$ (짧고 넓은 행렬).

행렬이 **행 완전 계수 (Full Row Rank)** 를 가지려면 $r = m$.

**예 2.** $A\mathbf{x} = \mathbf{b}$가 $n = 3$개의 미지수를 갖지만 $m = 2$개의 방정식만 있을 때:

$$x + y + z = 3, \quad x + 2y - z = 4$$

$$(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 1 & 2 & -1 & | & 4 \end{pmatrix} \xrightarrow{R_2 - R_1} \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & 3 & | & 2 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} = (R | \mathbf{d})$$

$\text{rank}(A) = 2$, $n = 3$, $n - r = 1$. 1개의 자유 변수, 1개의 특수해.

**i) $\mathbf{x}_p$:**

$$\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$$

$$x + 3z = 2, \quad y - 2z = 1$$

$z = 0$으로 취하면: $\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$, $\mathbf{x}_p = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}$ (특수해).

**ii) $\mathbf{x}_n$:**

$$\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$z = 1$로 취하면: $x + 3 = 0 \Rightarrow x = -3$, $y - 2 = 0 \Rightarrow y = 2$.

$\therefore \mathbf{x}_n = \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}$

**iii) 완전해:**

$$\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix} + \alpha \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}$$

이것은 $\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$을 지나는 방향의 **직선** 이며 ($A\mathbf{x}_n = \mathbf{0}$), $\mathbf{x}_p$만큼 이동한 것이다 ($A\mathbf{x} = \mathbf{b}$).

**행 완전 계수 ($r = m$)인 모든 행렬 $A$는 다음 성질을 가진다:**
1. 모든 행에 피벗이 있고, $R_0$에 영행이 없다. $R_0 = R$.
2. $A\mathbf{x} = \mathbf{b}$는 모든 $\mathbf{b}$에 대해 해를 가진다.
3. $A$의 열공간은 전체 공간 $\mathbb{R}^m$이다.

$A\mathbf{x} = \mathbf{b} \to R\mathbf{x} = \mathbf{d}$이고, $R = (I_{m \times m} \quad F_{m \times (n-m)})$.

4. $m < n$이면, $A\mathbf{x} = \mathbf{b}$는 **부정 (Underdetermined)** 이다 (많은 해).

행 완전 계수 ($r = m$)일 때, $m$개의 행이 선형 독립이다. $A^T$의 열들이 LI이다. $A^T$의 영공간은 $\mathbb{Z} = \{\mathbf{0}\}$이다.

### 4.7 선형 방정식의 네 가지 경우

선형 방정식의 네 가지 경우는 계수 $r$에 따라 달라진다:

| 경우 | 형태 | $R_0$ 형식 | $A\mathbf{x} = \mathbf{b}$의 해 |
|:-----|:------|:-----------|:----------------------------------------|
| $r = m$ 이고 $r = n$ | 정사각, 가역 | $(I)$ | 1개의 해 |
| $r = m$ 이고 $r < n$ | 짧고 넓음 | $(I \quad F)$ | $\infty$개의 해 |
| $r < m$ 이고 $r = n$ | 키 크고 좁음 | $\begin{pmatrix} I \\ 0 \end{pmatrix}$ | 0 또는 1개의 해 |
| $r < m$ 이고 $r < n$ | 완전 계수 아님 | $\begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix}$ | 0 또는 $\infty$개의 해 |

---

<br>

## 5. 독립, 기저, 차원 (3.4)

### 5.1 독립 벡터

**(1) 독립 벡터 (Independent Vectors):** 유일한 영 결합

$$c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \cdots + c_n \mathbf{v}_n = \mathbf{0}$$

은 모든 $c_1 = c_2 = \cdots = c_n = 0$을 가진다.

이는 다음을 의미한다: 적어도 하나의 스칼라가 0이 아닌 경우, 예를 들어 $c_1 \neq 0$이면:

$$\mathbf{v}_1 = -\frac{c_2}{c_1}\mathbf{v}_2 - \frac{c_3}{c_1}\mathbf{v}_3 - \cdots - \frac{c_n}{c_1}\mathbf{v}_n$$

$\mathbf{v}_1$은 다른 벡터들의 선형 결합이다.

**(2) 벡터 $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$이 공간 $\mathbb{S}$를 생성한다 (Span): $\mathbb{S}$ = $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$의 모든 결합일 때.**

예: $\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$ $\implies$ $\begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}$. $\hat{i}, \hat{j}$가 공간 $\mathbb{R}^2$를 생성한다.

**(3)** 벡터 $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$이 $\mathbb{S}$의 **기저 (Basis)** 가 되려면: **선형 독립 (Linearly Independent)** 이면서 $\mathbb{S}$를 **생성 (Span)** 해야 한다.

공간의 모든 벡터는 기저 벡터들의 **유일한** 결합이다.

예: $\mathbb{R}^2 \ni \begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}$

**(4)** 벡터 공간 $\mathbb{S}$의 **차원 (Dimension)** 은 $\mathbb{S}$의 모든 기저에 있는 벡터의 수 $n$이다.

$A \in \mathbb{R}^{m \times n}$을 고려하자. $n$개의 열 중 $r$개가 독립이고, 나머지 $n - r$개의 열이 종속이다. $C(A)$의 차원은 $r$이며, 이는 $A$의 계수 (Rank)이다.

이 절의 네 가지 핵심 개념:
1. 독립 벡터 (Independent Vectors)
2. 공간의 생성 (Spanning a Space)
3. 공간의 기저 (Basis for a Space)
4. 공간의 차원 (Dimension of a Space)

### 5.2 영공간을 통한 선형 독립

**정의.** $A$의 열들이 **선형 독립 (Linearly Independent)** 이려면 $A\mathbf{x} = \mathbf{0}$의 유일한 해가 $\mathbf{x} = \mathbf{0}$이어야 한다.

$$\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n = \mathbf{0}$$

$x_1 = x_2 = \cdots = x_n = 0$일 때만 성립.

**기하학적 해석:**

예: $\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3 \in \mathbb{R}^3$이 같은 평면에 있지 않으면 $\implies$ 이 벡터들은 독립이다. $\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$의 어떤 결합도 $\mathbf{0}$이 되지 않는다 (자명한 경우 제외).

예: $\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3$이 $\mathbb{R}^3$에서 같은 평면에 있다. 이들은 종속이다. 예를 들어, $\mathbf{w}_3 = \mathbf{w}_1 + \mathbf{w}_2$ $\iff$ $1 \cdot \mathbf{w}_1 + 1 \cdot \mathbf{w}_2 - 1 \cdot \mathbf{w}_3 = \mathbf{0}$. 결합이 $\mathbf{0}$을 주지만, 0이 아닌 계수가 있으므로 종속이다.

**종속 여부의 빠른 판별:**

- Q. $\begin{pmatrix} 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 10^{-5} \end{pmatrix}$는 종속인가? **아니다.**
- Q. $\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} -1 \\ -1 \end{pmatrix}$는 종속인가? **그렇다.**
- Q. $\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 0 \end{pmatrix}$는 종속인가? **그렇다.**
- Q. $\mathbb{R}^2$에서 임의의 세 벡터는 종속이다. **참.**

**예 1.** $A = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix}$은 종속 열을 가지는가?

$A$의 영공간 $N(A)$를 확인하자.

$$A\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$\xrightarrow{R_2 - 2R_1, \; R_3 - R_1}$ $\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -1 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$, $\text{rank}(A) = 2 < 3$.

$x_3 = 1$로 취하면: $x_1 + 3 = 0 \Rightarrow x_1 = -3$, $x_2 - 1 = 0 \Rightarrow x_2 = 1$.

0이 아닌 계수가 영벡터를 만든다:

$$-3\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + 1\begin{pmatrix} 3 \\ 5 \\ 3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$A$의 열들이 독립이면, $\text{rank}(A) = r = n$이고 $A$의 영공간은 영벡터만 포함한다: $N(A) = \{\mathbf{0}\}$.

**$\mathbb{R}^m$에서 $n > m$이면 $n$개의 벡터는 반드시 선형 종속이다.**

예: 5개의 성분을 가진 7개의 열이 있다면 ($5 \times 7$ 행렬). 7개의 열벡터는 $\mathbb{R}^5$에서 온다. 5개의 행에 5개 초과의 피벗이 있을 수 없다. $A\mathbf{x} = \mathbf{0}$은 최소 $2 \;(= 7 - 5)$개의 자유 변수를 가진다. 즉, 0이 아닌 해를 가진다.

**참고:** $n \leq m$이면 열들이 종속일 수도 독립일 수도 있다. 소거가 $r$개의 피벗을 드러낸다.

### 5.3 부분공간을 생성하는 벡터

$A$를 행렬이라 하자. $C(A)$는 $A\mathbf{x}$의 모든 결합으로 이루어진 열공간이다:

$$\begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{v}_1 + x_2 \mathbf{v}_2 + \cdots + x_n \mathbf{v}_n$$

**예:** $A = \begin{pmatrix} 1 & 4 \\ 2 & 7 \\ 3 & 5 \end{pmatrix}$, $A^T = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 7 & 5 \end{pmatrix}$

$C(A)$는 $\mathbb{R}^3$에서의 평면이다. $C(A^T)$, 즉 $A$의 행공간은 $\mathbb{R}^2$이다.

$A \in \mathbb{R}^{m \times n}$에 대해: 행벡터는 $\mathbb{R}^n$에, 열벡터는 $\mathbb{R}^m$에 속한다.

### 5.4 벡터 공간의 기저

**정의.** 벡터 공간의 **기저 (Basis)** 는 두 가지 성질을 가진 벡터 수열이다:
> i) 기저 벡터들이 **선형 독립** 이다.
> ii) 그들이 공간을 **생성** 한다.

**예:** $\begin{pmatrix} x \\ y \end{pmatrix} = x\begin{pmatrix} 1 \\ 0 \end{pmatrix} + y\begin{pmatrix} 0 \\ 1 \end{pmatrix}$이고, $\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$.

i) $\hat{i}$와 $\hat{j}$는 LI이다.
ii) $\hat{i}, \hat{j}$가 $\mathbb{R}^2$를 생성한다.

결합 $\mathbf{x} = x\hat{i} + y\hat{j}$는 $\hat{i}, \hat{j}$가 LI이므로 유일하다.

**$\mathbf{v}$를 기저 벡터들의 결합으로 쓰는 방법은 하나뿐이다.**

**증명.** $\mathbf{v} = a_1 \mathbf{v}_1 + \cdots + a_n \mathbf{v}_n$이고 $\mathbf{v} = b_1 \mathbf{v}_1 + \cdots + b_n \mathbf{v}_n$라 하자.

빼면 영벡터를 얻는다:

$$\mathbf{0} = (a_1 - b_1)\mathbf{v}_1 + \cdots + (a_n - b_n)\mathbf{v}_n$$

$\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$의 LI에 의해: $a_1 - b_1 = 0$, $a_2 - b_2 = 0$, $\ldots$, $a_n - b_n = 0$.

따라서 $a_i = b_i$ ($i = 1, 2, \ldots, n$). $\square$

---

**예 3.** $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$의 열들이 $\mathbb{R}^2$의 "표준 기저 (Standard Basis)"를 만든다.

$\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$라 하자.

i) $\hat{i}, \hat{j}$는 LI이다.
ii) $\hat{i}, \hat{j}$가 $\mathbb{R}^2$를 생성한다.

$\therefore \hat{i}, \hat{j}$는 $\mathbb{R}^2$의 기저 벡터이다.

---

**예 4.** 모든 가역 $n \times n$ 행렬의 열들이 $\mathbb{R}^n$의 기저를 제공한다.

**가역 행렬 예:**

$$A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = R$$

3개의 LI 열벡터. $C(A) = \mathbb{R}^3$. 3개의 0이 아닌 피벗. $\text{rank}(A) = 3$. $N(A) = \{\mathbf{0}\}$.

**기저는 유일하지 않다.**

벡터 $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n$이 $n \times n$ 가역 행렬의 열들일 때 $\mathbb{R}^n$의 기저이다. $\mathbb{R}^n$은 무한히 많은 서로 다른 기저를 가진다.

---

**특이 행렬 (Singular Matrix) 예:**

$$B = \begin{pmatrix} 1 & 0 & 1 \\ 1 & 1 & 2 \\ 1 & 1 & 2 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 1 & 2 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 0 & 0 \end{pmatrix} = R_0$$

$C(A) \neq \mathbb{R}^3$. $\text{rank}(A) = 2 < 3$.

열들이 종속이면, **피벗 열 (Pivot Column)** 만 취한다.

예: $B$의 1번째, 2번째 열.

- 모든 독립 벡터 집합은 기저로 **확장** 할 수 있다. (예: $B$의 1, 2번째 열이 $\mathbb{R}^3$에서 평면을 생성.)
- 모든 생성 벡터 집합은 기저로 **축소** 할 수 있다.

---

**예 5.** $A = \begin{pmatrix} 2 & 4 \\ 3 & 6 \end{pmatrix}$

$\xrightarrow{R_2 - \frac{3}{2}R_1}$ $\begin{pmatrix} 2 & 4 \\ 0 & 0 \end{pmatrix}$ $\xrightarrow{R_1 / 2}$ $\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0$

$\text{rank}(A) = 1 < 2$. 하나의 피벗 열, 하나의 피벗 행.

---

**예 6.** $R_0 = \begin{pmatrix} 1 & 2 & 0 & 3 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix}$

1번째, 3번째 열이 피벗 열이다. $\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$이 $C(R_0)$의 기저 벡터이다.

$C(R_0)$는 $\mathbb{R}^3$에서의 $xy$ 평면이다.

또한, 2번째와 3번째 열벡터도 $C(R_0)$의 기저이다.

---

**모든 벡터 공간의 기저는 같은 수의 벡터를 포함한다.** 모든 기저에 있는 벡터의 수가 그 공간의 "**차원 (Dimension)**"이다.

### 5.5 벡터 공간의 차원

**벡터 공간의 차원:**

$\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m$과 $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n$이 같은 벡터 공간의 두 기저이면, $m = n$이다.

**증명.** $n > m$이라 하자. $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m$이 기저이므로, 각 $\mathbf{w}_i$ ($i = 1, 2, \ldots, n$)는 $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m$의 결합이어야 한다:

$$\mathbf{w}_1 = a_{11}\mathbf{v}_1 + a_{21}\mathbf{v}_2 + \cdots + a_{m1}\mathbf{v}_m$$
$$\mathbf{w}_2 = a_{12}\mathbf{v}_1 + a_{22}\mathbf{v}_2 + \cdots + a_{m2}\mathbf{v}_m$$
$$\vdots$$
$$\mathbf{w}_n = a_{1n}\mathbf{v}_1 + a_{2n}\mathbf{v}_2 + \cdots + a_{mn}\mathbf{v}_m$$

이로부터:

$$W = (\mathbf{w}_1 \quad \mathbf{w}_2 \quad \cdots \quad \mathbf{w}_n) = (\mathbf{v}_1 \quad \mathbf{v}_2 \quad \cdots \quad \mathbf{v}_m) \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} = VA$$

$A$는 $m \times n$ 행렬 (짧고 넓음)이다. $n > m$이므로, $A\mathbf{x} = \mathbf{0}$은 0이 아닌 해를 가진다.

$A\mathbf{x} = \mathbf{0}$에서: $VA\mathbf{x} = V\mathbf{0} = \mathbf{0}$.

즉 $W\mathbf{x} = \mathbf{0}$, 다시 말해 $x_1\mathbf{w}_1 + x_2\mathbf{w}_2 + \cdots + x_n\mathbf{w}_n = \mathbf{0}$.

$\mathbf{x}$가 0이 아닌 벡터이므로, $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n$은 LI가 아니다. $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n$은 기저가 될 수 **없다**. 이는 $\mathbf{w}_i$ ($i = 1, 2, \ldots, n$)가 기저라는 가정에 모순이다. $\square$

---

**정의.** 공간의 **차원 (Dimension)** 은 모든 기저에 있는 벡터의 수이다.

**예:** $\mathbf{u} = \begin{pmatrix} 1 \\ 5 \\ 2 \end{pmatrix}$를 지나는 직선은 1차원이다.

그 직선에 수직인 평면: $\mathbf{u} \cdot \mathbf{x} = 0 \iff x + 5y + 2z = 0$.

행렬 $A = \begin{pmatrix} 1 & 5 & 2 \end{pmatrix}$의 평면 영공간:

$(1 \quad 5 \quad 2)\begin{pmatrix} x \\ y \\ z \end{pmatrix} = 0$

$n = 3$, $\text{rank}(A) = r = 1$, $n - r = 2$개의 자유 변수.

$n - r$개의 특수해가 영공간의 기저를 제공한다: 차원 $n - r$.

특수해를 구하면:
- $y = 1, z = 0$으로 취하면 $\implies x = -5$
- $y = 0, z = 1$로 취하면 $\implies x = -2$

$\begin{pmatrix} -5 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} -2 \\ 0 \\ 1 \end{pmatrix}$이 기저이다 (2차원).

**차원 요약:**
- 행공간의 차원은 $r$
- 열공간의 차원은 $r$
- 영공간의 차원은 $n - r$
- $N(A^T)$의 차원은 $m - r$

### 5.6 행렬 공간과 함수 공간의 기저

**행렬 공간과 함수 공간의 기저**

"독립," "기저," "차원"이라는 용어는 열벡터에만 한정되지 않는다.

**행렬 공간 (Matrix Spaces):**

벡터 공간 $\mathbb{M}$은 모든 $2 \times 2$ 행렬을 포함한다. 차원은 4이다.

예: $A_1, A_2, A_3, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

i) LI이다: $x_1 A_1 + x_2 A_2 + x_3 A_3 + x_4 A_4 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$ $\iff$ $x_1 = x_2 = x_3 = x_4 = 0$.

ii) $\mathbb{M}$을 생성한다: $\begin{pmatrix} a & b \\ c & d \end{pmatrix} = aA_1 + bA_2 + cA_3 + dA_4$. $A_1, A_2, A_3, A_4$의 선형 결합이 $\mathbb{M}$의 임의의 행렬을 만들 수 있다.

**행렬 공간의 부분공간:**

- $A_1, A_2, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$는 상삼각 행렬 (Upper Triangular Matrix) 부분공간 $\mathcal{U}$의 기저이다: $\mathcal{U} \ni \begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

- $A_1, A_4$는 대각 행렬 (Diagonal Matrix) 부분공간 $\mathbb{D}$의 기저이다: $\mathbb{D} \ni \begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

- $A_1, A_4, A_2 + A_3$는 대칭 행렬 (Symmetric Matrix)의 기저이다: $\mathbb{S} \ni \begin{pmatrix} a & b \\ b & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$

**행렬 부분공간의 차원 ($n \times n$ 행렬에 대해):**

| 공간 | 차원 |
|:------|:----------|
| 전체 $n \times n$ 행렬 공간 | $n^2$ |
| 대각 행렬 | $n$ |
| 상삼각 행렬 | $\frac{1}{2}n^2 + \frac{1}{2}n$ |
| 대칭 행렬 | $\frac{1}{2}n^2 + \frac{1}{2}n$ |

상삼각 행렬의 경우: 0의 총 개수 $= \sum_{i=1}^{n}(i-1) = \sum_{i=1}^{n} i - n = \frac{n(n+1)}{2} - n$. 0이 아닌 원소의 총 개수 $= n^2 - \frac{n(n+1)}{2} + n = \frac{n^2}{2} + \frac{n}{2}$.

상삼각 행렬과 대칭 행렬의 원소 수 (차원)은 같다: $\frac{1}{2}n^2 + \frac{1}{2}n$.

---

**함수 공간 (Function Spaces):**

$$\frac{d^2y}{dx^2} = 0, \quad \frac{d^2y}{dx^2} = -y, \quad \frac{d^2y}{dx^2} = y$$

이들은 2계 도함수를 포함한다. 미적분학에서 해 $y(x)$를 구한다:

| 상미분방정식 (ODE) | 해 | 기저 | 차원 |
|:----|:---------|:------|:----------|
| $y'' = 0$ | $y = cx + d$ | $\{1, x\}$ | 2 |
| $y'' = -y$ | $y = c\sin x + d\cos x$ | $\{\sin x, \cos x\}$ | 2 |
| $y'' = y$ | $y = ce^x + de^{-x}$ | $\{e^x, e^{-x}\}$ | 2 |

기저 벡터들은 2계 도함수의 **영공간 (Nullspace)** 에 속한다.

$y'' = 2$는 특수해 $y_p = x^2$을 가진다: $\frac{dy_p}{dx} = 2x$, $\frac{d^2 y_p}{dx^2} = 2$.

따라서 일반해 (완전해)는:

$$y(x) = y_p(x) + y_n(x) = x^2 + cx + d$$

---

<br>

## 6. 네 부분공간의 차원 (3.5)

### 6.1 차원 요약

1. 열공간 $C(A)$와 행공간 $C(A^T)$는 모두 차원 $r$ (= $A$의 계수)을 가진다.
2. $A$의 영공간 $N(A)$는 차원 $n - r$을 가진다.
3. $A$의 좌영공간 (Left Nullspace) $N(A^T)$는 차원 $m - r$을 가진다.
4. $A$에서 $R_0$로의 소거는 $C(A)$와 $N(A^T)$를 **바꾸지만**, 차원은 변하지 않는다.

### 6.2 부분공간의 직교성

"**계수 (Rank)**"와 "**차원 (Dimension)**"을 연결하려 한다:
- **계수** 는 독립 열의 수를 센다.
- **차원** 은 기저에 있는 벡터의 수이다.

$A \in \mathbb{R}^{m \times n}$의 계수는 네 가지 기본 부분공간 모두의 차원을 알려준다:

1. **행공간 (Row Space)**, $C(A^T)$는 $\mathbb{R}^n$의 부분공간, 차원 $r$. (각 행벡터 $\in \mathbb{R}^n$.)
2. **열공간 (Column Space)**, $C(A)$는 $\mathbb{R}^m$의 부분공간, 차원 $r$. (각 열벡터 $\in \mathbb{R}^m$.)
3. **영공간 (Nullspace)**, $N(A)$는 $\mathbb{R}^n$의 부분공간, 차원 $n - r$. ($A\mathbf{x} = \mathbf{0}$, $\mathbf{x} \in \mathbb{R}^n$.)
4. **좌영공간 (Left Nullspace)**, $N(A^T)$는 $\mathbb{R}^m$의 부분공간, 차원 $m - r$. ($A^T\mathbf{y} = \mathbf{0}$, $\mathbf{y} \in \mathbb{R}^m$.)

**직교 쌍 (Orthogonal Pairs):**

$$C(A) \subset \mathbb{R}^m, \quad C(A^T) \subset \mathbb{R}^n$$
$$N(A^T) \subset \mathbb{R}^m, \quad N(A) \subset \mathbb{R}^n$$

- $C(A)$와 $C(A^T)$는 $r$차원이다.
- $A$의 열공간과 행공간은 같은 차원 $r$을 가진다.
- $N(A)$는 차원 $n - r$.
- $N(A^T)$는 차원 $m - r$.

**$N(A)$는 $C(A^T)$ (행공간)에 수직이다:**

예: $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$

행벡터 $\{(1, 2), (3, 4)\}$가 $A$의 행공간을 생성한다. $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$의 해는 $A$의 영공간 $N(A)$에 속한다.

$N(A) \ni \mathbf{x}$는 내적 (Inner Product)의 의미에서 $A$의 행공간의 모든 벡터에 수직이다.

즉, $\mathbf{y} \in C(A^T)$이고 $\mathbf{x} \in N(A)$이면: $\mathbf{y} \cdot \mathbf{x} = 0$.

**$N(A^T)$는 $C(A)$ (열공간)에 수직이다:**

예: $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$

열벡터 $\{\begin{pmatrix} 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ 4 \end{pmatrix}\}$가 $A$의 열공간을 생성한다.

$N(A^T) \ni \mathbf{y}$이면 $A^T\mathbf{y} = \mathbf{0}$:

$\begin{pmatrix} 1 & 3 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$\Rightarrow C(A) \ni \mathbf{x} \perp \mathbf{y} \in N(A^T)$

즉, $\alpha\begin{pmatrix} 1 \\ 3 \end{pmatrix} + \beta\begin{pmatrix} 2 \\ 4 \end{pmatrix} \perp \mathbf{y}$

### 6.3 R_0의 네 부분공간

**$R_0 = \text{rref}(A)$라 하자.** 네 차원은 $R_0$와 $A$에 대해 같다.

**예:** $3 \times 5$ 행렬 $R_0$를 고려하자:

$$R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}$$

피벗 행: 1과 2. 피벗 열: 1과 4. $\text{rank}(R_0) = r = 2$.

**행공간:** 기저 벡터 $\{(1, 3, 5, 0, 7), \;(0, 0, 0, 1, 2)\}$로 생성된다. $\dim C(R_0^T) = 2 = r$.

**열공간:** 1번째와 4번째 열벡터 $\left\{\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\}$가 $C(R_0)$의 기저이다. $\dim C(R_0) = 2 = r$.

**영공간:** $R_0 \mathbf{x} = \mathbf{0}$:

$$\begin{pmatrix} x_1 + 3x_2 + 5x_3 + 7x_5 \\ x_4 + 2x_5 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$n - r = 5 - 2 = 3$, 3개의 자유 변수.

i) $x_2 = 1, x_3 = 0, x_5 = 0$으로 취하면: $x_1 + 3 = 0 \Rightarrow x_1 = -3$, $x_4 = 0$. $\mathbf{x} = \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}$

ii) $x_2 = 0, x_3 = 1, x_5 = 0$으로 취하면: $x_1 + 5 = 0 \Rightarrow x_1 = -5$, $x_4 = 0$. $\mathbf{x} = \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}$

iii) $x_2 = 0, x_3 = 0, x_5 = 1$로 취하면: $x_1 + 7 = 0 \Rightarrow x_1 = -7$, $x_4 + 2 = 0 \Rightarrow x_4 = -2$. $\mathbf{x} = \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}$

3개의 특수해가 $N(R_0)$의 기저 $\left\{\begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}\right\}$를 형성한다.

$\dim N(R_0) = 3 = 5 - 2 = n - r$.

**좌영공간:** $R_0^T \mathbf{y} = \mathbf{0}$:

$$\begin{pmatrix} 1 & 0 & 0 \\ 3 & 0 & 0 \\ 5 & 0 & 0 \\ 0 & 1 & 0 \\ 7 & 2 & 0 \end{pmatrix} \begin{pmatrix} y_1 \\ y_2 \\ y_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix} \iff \begin{pmatrix} y_1 \\ 3y_1 \\ 5y_1 \\ y_2 \\ 7y_1 + 2y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}$$

$m - r = 3 - 2 = 1$, 1개의 자유 변수.

$y_3 = 1$로 취하면: $\Rightarrow y_1 = y_2 = 0$. $\therefore \mathbf{y} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$

$\left\{\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\right\}$이 $N(R_0^T)$의 기저이다.

$\dim N(R_0^T) = 1 = m - r = 3 - 2$.

$R_0^T \mathbf{y} = \mathbf{0} \iff \mathbf{y}^T R_0 = \mathbf{0}$. $\mathbf{y}^T$는 $R_0$의 '왼쪽'에 있는 행벡터이다.

**직교성 요약:**

$$C(A) \subset \mathbb{R}^m \perp N(A^T) \subset \mathbb{R}^m$$
$$C(A^T) \subset \mathbb{R}^n \perp N(A) \subset \mathbb{R}^n$$

$\mathbb{R}^n$에서: 행공간과 영공간은 $r$과 $n - r$ 차원을 가진다.

$\mathbb{R}^m$에서: 열공간과 좌영공간은 $r$과 $m - r$ 차원을 가진다.

### 6.4 A와 R_0의 관계

**$A$의 네 부분공간 차원은 $R_0$와 같다.**

$$A = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 1 & 3 & 5 & 1 & 9 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}$$

$A$는 $R_0$와 같은 행공간을 가지지만, 열공간은 $R_0$와 다르다.

$C(A)$의 기저: $\left\{\begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}\right\}$ ($R_0$가 아닌 $A$의 피벗 열)

$C(R_0)$의 기저: $\left\{\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\}$

($\mathbb{R}^3$에서 서로 다른 평면이지만, 같은 차원.)

**핵심 관계:**

1. $A$는 $R_0$, $R$과 **같은 행공간** 을 가진다: $C(A^T) = C(R_0^T) = C(R^T)$, 차원 $r$.

2. $A$의 열공간은 차원 $r$을 가진다: $\dim C(A) = \dim C(A^T)$. $C(A) \neq C(R_0)$이지만, $\dim C(A) = \dim C(R_0) = r$.

3. $A$는 $R_0$와 **같은 영공간** 을 가진다: $A\mathbf{x} = \mathbf{0} \iff R_0 \mathbf{x} = \mathbf{0}$. **소거는 해를 바꾸지 않는다.** 특수해가 기저를 형성한다. $n - r$개의 자유 변수 $\implies$ $\dim N(A) = n - r$.

$$\dim C(A) + \dim N(A) = r + (n - r) = n$$

4. $A$의 좌영공간 $N(A^T)$: $\dim C(A^T) + \dim N(A^T) = r + (m - r) = m$.

### 6.5 선형대수학의 기본 정리

$$\boxed{\textbf{선형대수학의 기본 정리 (Fundamental Theorem of Linear Algebra)}}$$

> i) $\dim C(A) = \dim C(A^T) = r$
>
> ii) $\dim N(A) = n - r$, $\quad \dim N(A^T) = m - r$

**$A$의 네 부분공간:**

$$C(A^T) \subset \mathbb{R}^n \quad \perp \quad N(A) \subset \mathbb{R}^n$$
$$C(A) \subset \mathbb{R}^m \quad \perp \quad N(A^T) \subset \mathbb{R}^m$$

$\mathbb{R}^n$에서: $C(A^T)$ (차원 $r$)과 $N(A)$ (차원 $n - r$)는 직교 여공간 (Orthogonal Complements)이다.

$\mathbb{R}^m$에서: $C(A)$ (차원 $r$)과 $N(A^T)$ (차원 $m - r$)는 직교 여공간이다.

---

<br>

## 요약

| 개념 | 핵심 내용 |
|:--------|:---------|
| 벡터 공간 (Vector Space) | 덧셈과 스칼라 곱에 닫혀 있으며 8가지 공리를 만족하는 집합 |
| 체 (Field) | $+, -, \times, \div$가 정의된 집합 ($\mathbb{R}$, $\mathbb{C}$) |
| 부분공간 (Subspace) | 그 자체가 벡터 공간인 벡터 공간의 부분집합 ($\mathbf{0}$ 포함 필수) |
| 열공간 $C(A)$ | $A$의 열들의 모든 선형 결합; $A\mathbf{x} = \mathbf{b}$가 풀리려면 $\mathbf{b} \in C(A)$ |
| 행공간 $C(A^T)$ | $A^T$의 열공간; $A$의 행들이 생성 |
| 영공간 $N(A)$ | $A\mathbf{x} = \mathbf{0}$의 모든 해; $\mathbb{R}^n$의 부분공간 |
| 좌영공간 $N(A^T)$ | $A^T\mathbf{y} = \mathbf{0}$인 모든 $\mathbf{y}$; $\mathbb{R}^m$의 부분공간 |
| 기약행사다리꼴 $R_0 = \text{rref}(A)$ | 피벗 열에 $I$, 자유 열에 $F$를 포함 |
| $A = CR$ | $C$ = $A$의 독립 열; $R = (I \; F)$ = 기약행사다리꼴 (영행 제거) |
| 특수해 (Special Solutions) | $\begin{pmatrix} -F \\ I \end{pmatrix}$의 열들; $N(A)$의 기저를 형성 |
| 완전해 (Complete Solution) | $\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n$ (특수해 + 영공간) |
| 특수해 $\mathbf{x}_p$ (Particular Solution) | 모든 자유 변수를 0으로 놓고 $R_0\mathbf{x} = \mathbf{d}$ 풀기 |
| 열 완전 계수 ($r = n$) | $N(A) = \{\mathbf{0}\}$; $A\mathbf{x} = \mathbf{b}$의 해가 최대 1개; $R_0 = \begin{pmatrix} I \\ 0 \end{pmatrix}$ |
| 행 완전 계수 ($r = m$) | $A\mathbf{x} = \mathbf{b}$가 항상 풀림; $C(A) = \mathbb{R}^m$; $R = (I \; F)$ |
| 선형 독립 (Linear Independence) | $c_1\mathbf{v}_1 + \cdots + c_n\mathbf{v}_n = \mathbf{0}$이 모든 $c_i = 0$일 때만 성립 |
| 생성 (Spanning) | 벡터들이 $\mathbb{S}$를 생성: $\mathbb{S}$ = 그 벡터들의 모든 결합 |
| 기저 (Basis) | 선형 독립이면서 공간을 생성하는 벡터들; 표현이 유일 |
| 차원 (Dimension) | 모든 기저에 있는 벡터의 수; 기저에 관계없이 불변 |
| $\dim C(A) = \dim C(A^T) = r$ | 열공간과 행공간은 같은 차원 (계수) |
| $\dim N(A) = n - r$ | 영공간 차원은 자유 변수의 수 |
| $\dim N(A^T) = m - r$ | 좌영공간 차원 |
| $r + (n - r) = n$ | 행공간 + 영공간이 $\mathbb{R}^n$을 채움 |
| $r + (m - r) = m$ | 열공간 + 좌영공간이 $\mathbb{R}^m$을 채움 |
| 직교성 (Orthogonality) | $N(A) \perp C(A^T)$ ($\mathbb{R}^n$); $N(A^T) \perp C(A)$ ($\mathbb{R}^m$) |
| 행렬 공간 차원 | $n \times n$: $n^2$; 대각: $n$; 상삼각: $\frac{n^2+n}{2}$; 대칭: $\frac{n^2+n}{2}$ |
| 함수 공간 (Function Space) | $y'' = 0, -y, y$의 해가 기저 $\{1,x\}$, $\{\sin x, \cos x\}$, $\{e^x, e^{-x}\}$인 2차원 공간 형성 |
| 기본 정리 (Fundamental Theorem) | $\dim C(A) = \dim C(A^T) = r$; $\dim N(A) = n-r$; $\dim N(A^T) = m-r$ |

---
