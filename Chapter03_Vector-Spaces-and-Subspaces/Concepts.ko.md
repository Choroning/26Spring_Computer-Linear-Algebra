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

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.** GitHub은 파일 하나당 렌더되는 수식 개수에 한도가 있어, 모든 수식이 제대로 보이도록 절 단위로 나눴습니다.
>
> **Part 1** · [Part 2](Concepts-2.ko.md)

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
- [5. 독립, 기저, 차원 (3.4)](Concepts-2.ko.md#5-독립-기저-차원-34)
  - [5.1 독립 벡터](Concepts-2.ko.md#51-독립-벡터)
  - [5.2 영공간을 통한 선형 독립](Concepts-2.ko.md#52-영공간을-통한-선형-독립)
  - [5.3 부분공간을 생성하는 벡터](Concepts-2.ko.md#53-부분공간을-생성하는-벡터)
  - [5.4 벡터 공간의 기저](Concepts-2.ko.md#54-벡터-공간의-기저)
  - [5.5 벡터 공간의 차원](Concepts-2.ko.md#55-벡터-공간의-차원)
  - [5.6 행렬 공간과 함수 공간의 기저](Concepts-2.ko.md#56-행렬-공간과-함수-공간의-기저)
- [6. 네 부분공간의 차원 (3.5)](Concepts-2.ko.md#6-네-부분공간의-차원-35)
  - [6.1 차원 요약](Concepts-2.ko.md#61-차원-요약)
  - [6.2 부분공간의 직교성](Concepts-2.ko.md#62-부분공간의-직교성)
  - [6.3 R_0의 네 부분공간](Concepts-2.ko.md#63-r_0의-네-부분공간)
  - [6.4 A와 R_0의 관계](Concepts-2.ko.md#64-a와-r_0의-관계)
  - [6.5 선형대수학의 기본 정리](Concepts-2.ko.md#65-선형대수학의-기본-정리)
- [요약](Concepts-2.ko.md#요약)

---

<br>

## 1. 개요: 네 가지 기본 부분공간

제3장은 다섯 가지 주요 주제를 다룬다:

- **3.1** 벡터 공간과 부분공간 (Vector Spaces and Subspaces) — 벡터 공간을 어떻게 정의하는가? 핵심 연산은 $`\mathbf{u} + \mathbf{v}`$와 $`c\mathbf{u}`$이다. 벡터 $`\mathbf{u}`$와 스칼라 $`c`$가 만족해야 하는 8가지 규칙이 있다.
- **3.2** $`A`$의 영공간 (Null Space): $`A\mathbf{x} = \mathbf{0}`$ 풀기
- **3.3** $`A\mathbf{x} = \mathbf{b}`$의 완전해 (Complete Solution)
- **3.4** 독립 (Independence), **기저 (Basis)**, **차원 (Dimension)** — 공간을 설명하는 벡터 집합. $`A \in \mathbb{R}^{n \times n}`$이라 하자. $`A`$가 $`r`$개의 독립 열을 가지면 $`C(A)`$의 차원은 $`r`$이다. $`A\mathbf{x} = \mathbf{0}`$의 $`(n - r)`$개 특수해는 $`A`$의 영공간 $`N(A)`$의 기저가 된다.
- **3.5** 네 부분공간의 차원 (Dimensions of the Four Subspaces):

| 부분공간 | 차원 |
|:---------|:----------|
| $`A`$의 열공간 (Column Space) | $`r`$ |
| $`A`$의 행공간 (Row Space) | $`r`$ |
| $`A`$의 영공간 (Null Space) | $`n - r`$ |
| $`A^T`$의 영공간 (좌영공간, Left Nullspace) | $`m - r`$ |

---

<br>

## 2. 벡터 공간과 부분공간 (3.1)

### 2.1 벡터 공간의 요건

1. **요건:** 모든 선형 결합 $`c\mathbf{u} + d\mathbf{w}`$가 벡터 공간 안에 머물러야 한다.
2. $`A`$의 **행공간 (Row Space)** 은 $`A`$의 행들에 의해 "생성"된다. $`A`$의 열들은 열공간 $`C(A)`$를 생성한다.
3. **행렬 (Matrix)** $`M_1`$부터 $`M_n`$까지와 **함수 (Function)** $`f_1`$부터 $`f_n`$까지가 **행렬 공간 (Matrix Space)** 과 **함수 공간 (Function Space)** 을 생성한다.

### 2.2 공간 R^n

공간 $`\mathbb{R}^n`$은 길이 $`n`$인 모든 열벡터 $`\mathbf{v}`$를 포함한다.

```math
\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} \in \mathbb{R}^n \quad \text{여기서 } x_1, x_2, \ldots, x_n \in \mathbb{R} \text{ (실수)}
```

예시:
- $`x \in \mathbb{R}^1`$
- $`\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} \in \mathbb{R}^2`$
- $`\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} \in \mathbb{R}^3`$

**참고:** $`x_1, x_2, \ldots, x_n \in \mathbb{C}`$ (복소수)이면, 공간은 $`\mathbb{C}^n`$이 된다.

**참고:** $`c, d \in \mathbb{R}`$이고 $`\mathbf{u}, \mathbf{w} \in \mathbb{R}^n`$이면, $`c\mathbf{u} + d\mathbf{w} \in \mathbb{R}^n`$이다. 선형 결합은 벡터 공간 $`\mathbb{R}^n`$ 안에 머문다.

### 2.3 벡터 공간의 정의

**벡터 공간 (Vector Space)** (= 선형 공간, Linear Space) $`V`$는 원소(벡터)들이 서로 더해지고 수와 곱해질 수 있는 집합이다.

$`\mathbf{u}, \mathbf{w} \in V`$ (벡터 공간)일 때:
1. $`\mathbf{u} + \mathbf{w} \in V`$ — 벡터 덧셈 (Vector Addition)
2. $`\alpha \mathbf{u} \in V \quad \forall \alpha \in \mathbb{F}`$ — 스칼라 곱 (Scalar Multiplication)

**체 (Field)** 는 덧셈, 뺄셈, 곱셈, 나눗셈이 정의된 집합이다.

체의 예: $`\mathbb{R}`$ (실수), $`\mathbb{C}`$ (복소수).

### 2.4 여덟 가지 공리

$`\mathbf{u}, \mathbf{v}, \mathbf{w} \in V`$이고 $`c, d \in \mathbb{F}`$일 때. 벡터 공간은 다음 여덟 가지 공리를 만족해야 한다:

**(1) 벡터 덧셈의 결합법칙 (Associativity):**

```math
\mathbf{u} + (\mathbf{v} + \mathbf{w}) = (\mathbf{u} + \mathbf{v}) + \mathbf{w}
```

**(2) 벡터 덧셈의 교환법칙 (Commutativity):**

```math
\mathbf{u} + \mathbf{v} = \mathbf{v} + \mathbf{u}
```

**(3) 벡터 덧셈의 항등원 (영벡터, Zero Vector):**

```math
\exists!\; \mathbf{0} \in V \text{ s.t. } \mathbf{v} + \mathbf{0} = \mathbf{v} \quad \forall \mathbf{v} \in V
```

**(4) 벡터 덧셈의 역원 (Inverse Element):**

```math
\forall \mathbf{v} \in V, \; \exists!\; -\mathbf{v} \in V \text{ s.t. } \mathbf{v} + (-\mathbf{v}) = \mathbf{0}
```

($`-\mathbf{v}`$는 $`\mathbf{v}`$의 덧셈 역원 (Additive Inverse))

**(5) 스칼라 곱과 체 곱의 양립성 (Compatibility):**

```math
c(d\mathbf{v}) = (cd)\mathbf{v}
```

**(6) 스칼라 곱의 항등원 (Identity Element):**

```math
1 \cdot \mathbf{v} = \mathbf{v}
```

**(7) 벡터 덧셈에 대한 스칼라 곱의 분배법칙 (Distributivity):**

```math
c(\mathbf{u} + \mathbf{v}) = c\mathbf{u} + c\mathbf{v}
```

**(8) 체 덧셈에 대한 스칼라 곱의 분배법칙 (Distributivity):**

```math
(c + d)\mathbf{u} = c\mathbf{u} + d\mathbf{u}
```

### 2.5 공리의 결과

여덟 가지 공리로부터 다음 성질이 성립한다:

```math
0\mathbf{u} = \mathbf{0}
```

```math
c\mathbf{0} = \mathbf{0}
```

```math
(-1)\mathbf{u} = -\mathbf{u}
```

```math
c\mathbf{v} = \mathbf{0} \implies c = 0 \text{ 또는 } \mathbf{v} = \mathbf{0}
```

### 2.6 예제: 벡터 공간인가?

**Q.** 모든 양수 벡터의 집합 $`\mathbb{X}`$, 즉 $`\mathbb{X} \ni \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}`$이고 모든 $`v_i > 0`$인 경우, 벡터 공간인가?

**A.** 아니다. $`-\mathbf{v} \notin \mathbb{X}`$.

---

**Q.** $`\mathbb{X}`$를 $`A\mathbf{x} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$의 해 집합이라 하자. $`\mathbb{X}`$는 벡터 공간인가?

**A.** 아니다. $`\mathbf{u}, \mathbf{w} \in \mathbb{X}`$이면, $`A\mathbf{u} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$이고 $`A\mathbf{w} = \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$이다. 하지만 $`\mathbf{u} + \mathbf{w} \notin \mathbb{X}`$인데, $`A(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w} = \begin{pmatrix} 2 \\ \vdots \\ 2 \end{pmatrix} \neq \begin{pmatrix} 1 \\ \vdots \\ 1 \end{pmatrix}`$이기 때문이다.

---

**Q.** $`\mathbb{R}^n`$에서 직선은 벡터 공간인가?

직선은 점들의 모임이다. $`\mathbf{q}`$는 직선 위의 위치를 나타낸다:

```math
\mathbf{q} = \mathbf{p} + t\mathbf{d}, \quad t \in \mathbb{R}
```

$`\mathbf{p}`$와 $`\mathbf{d}`$는 고정되어 있고, $`t`$는 $`-\infty`$에서 $`\infty`$까지 변한다.

$`\mathbb{X}`$를 모든 벡터 $`\mathbf{q}`$의 집합이라 하자: $`\mathbb{X} = \lbrace\mathbf{q} : \mathbf{p} + t\mathbf{d} \;\forall\; t \in \mathbb{R}\rbrace`$.

$`\mathbf{a}, \mathbf{b} \in \mathbb{X}`$를 취하자. $`\mathbf{a} + \mathbf{b}`$가 $`\mathbb{X}`$에 속하는가?

```math
\mathbf{a} = \mathbf{p} + t_1 \mathbf{d}
```
```math
\mathbf{b} = \mathbf{p} + t_2 \mathbf{d}
```
```math
\mathbf{a} + \mathbf{b} = 2\mathbf{p} + (t_1 + t_2)\mathbf{d} \notin \mathbb{X} \quad \text{만약 } \mathbf{p} \neq \mathbf{0}
```

그러나:
- $`\mathbf{a} + \mathbf{b} = (t_1 + t_2)\mathbf{d} \in \mathbb{X}`$ (만약 $`\mathbf{p} = \mathbf{0}`$)
- $`c\mathbf{a} = ct_1\mathbf{d} \in \mathbb{X}`$ (만약 $`\mathbf{p} = \mathbf{0}`$)

**원점을 지나는 직선은 벡터 공간이다.** $`\mathbb{R}^n`$에서 $`\mathbf{0}`$을 지나는 직선은 $`\mathbb{R}^n`$의 **부분공간 (Subspace)** 이다 — 다른 벡터 공간 안의 벡터 공간.

### 2.7 벡터 공간의 예시

- $`\mathbb{R}^n`$은 벡터 공간이다.
- $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$은 벡터 공간이다 (가장 작은 벡터 공간):
  - i) $`\mathbf{0} + \mathbf{0} = \mathbf{0} \in \mathbb{Z}`$
  - ii) $`c\mathbf{0} = \mathbf{0} \in \mathbb{Z}`$

$`\mathbb{Z}`$는 **가역 행렬의 영공간 (Null Space)** 이다. $`A\mathbf{x} = \mathbf{0}`$의 유일한 해가 영벡터 $`\mathbf{x} = \mathbf{0}`$이면, $`A`$의 열들은 선형 독립 (LI)이고 $`A`$의 영공간은 $`\mathbb{Z}`$이다.

### 2.8 일반화된 벡터 공간

**참고:** 벡터 공간 (= 선형 공간)은 원소들이 서로 더해지고 수와 곱해질 수 있는 집합이다.

즉, $`\mathbf{u}, \mathbf{w} \in V \implies c\mathbf{u} + d\mathbf{w} \in V`$.

**행렬과 함수도 벡터로 간주할 수 있다.**

**행렬의 벡터 공간 (Vector Space of Matrices):**

$`A, B \in \mathbb{R}^{m \times n} \implies cA + dB \in \mathbb{R}^{m \times n}`$

고정된 크기의 모든 행렬의 집합은 벡터 공간을 형성한다. (여덟 가지 규칙을 만족하는지 확인하라.)

**함수의 벡터 공간 (Vector Space of Functions):**

$`\mathbb{F}`$를 $`\mathbb{R}`$에서 원소를 취하여 실수로 사상하는 함수들의 집합이라 하자:

```math
\mathbb{F} = \lbracef \mid f: \mathbb{R} \to \mathbb{R}\rbrace
```

$`f, g \in \mathbb{F}`$이고 $`c, d \in \mathbb{R}`$이면, $`cf + dg \in \mathbb{F}`$이다.

정의:
- $`(f + g)(x) = f(x) + g(x)`$
- $`(cf)(x) = c(f(x))`$

**예:** $`\mathbb{F}`$ = 함수 $`y = ce^x`$의 직선.

$`\mathbb{F} = \lbracef \mid f: \mathbb{R} \ni x \to ce^x \in \mathbb{R}, \;\forall c \in \mathbb{R}\rbrace`$

$`f(x) = e^x`$, $`g(x) = 2e^x`$.

$`(f + g)(x) = f(x) + g(x) = e^x + 2e^x = 3e^x \in \mathbb{F}`$

$`(cf)(x) = c(f(x)) = ce^x \in \mathbb{F}`$

**비고:** **집합 (Set)** 은 추가적인 구조 없이 원소들의 모임에 불과하다. **공간 (Space)** 은 추가적인 구조가 정의된 집합이다. 예: 벡터 공간은 벡터 덧셈과 스칼라 곱이 정의된 벡터들의 집합이다. 벡터 공간은 8가지 규칙을 만족해야 한다.

### 2.9 벡터 공간의 부분공간

벡터 공간 $`\mathbb{R}^n`$을 고려하자. 여기서 $`\mathbf{v} \in \mathbb{R}^n`$은 $`n`$개의 성분을 가진 열벡터이다.

$`\mathbb{R}^n`$ 안에 중요한 벡터 공간들이 있다. 이들이 $`\mathbb{R}^n`$의 **부분공간 (Subspace)** 이다.

평면은 벡터 공간이다. $`\mathbf{v}, \mathbf{w} \in \mathbb{R}^2`$이면, $`c\mathbf{v} + d\mathbf{w} \in \mathbb{R}^2`$이다.

$`\mathbb{R}^3`$에서 원점 $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$을 지나는 평면은 벡터 공간이다. 이 평면은 $`\mathbb{R}^2`$가 아닌데, $`\mathbf{u}, \mathbf{v} \in \mathbb{R}^3`$이기 때문이다. 이 평면은 $`\mathbb{R}^3`$ 안의 벡터 공간이다. 이것은 전체 벡터 공간 $`\mathbb{R}^3`$의 **부분공간** 이다.

### 2.10 부분공간의 정의

**정의.** 벡터 공간의 **부분공간 (Subspace)** 은 다음을 만족하는 벡터들의 집합($`\mathbf{0}`$ 포함)이다:

> i) $`\mathbf{u} + \mathbf{w}`$가 부분공간에 속한다
>
> ii) $`c\mathbf{u}`$가 부분공간에 속한다

$`\mathbf{u}, \mathbf{w} \in`$ 부분공간이고 $`c \in \mathbb{R}`$ (또는 $`\mathbb{F}`$)일 때.

벡터 집합이 덧셈 $`\mathbf{u} + \mathbf{w}`$과 곱셈 $`c\mathbf{u}`$에 대해 "**닫혀 (Closed)**" 있다.

조건 i)과 ii)는 다음을 의미한다: **모든 선형 결합이 부분공간 안에 머문다.**

**참고:** 모든 부분공간은 영벡터를 포함한다. ii)에서 $`c = 0`$으로 놓으면 $`0\mathbf{u} = \mathbf{0}`$이 부분공간에 속한다.

### 2.11 부분공간과 비부분공간의 예시

**Q:** $`\mathbb{R}^3`$에서 평면 $`z = 5`$는 부분공간인가? **아니다.** (원점을 포함하지 않는다.)

원점을 지나는 직선은 부분공간이다.

$`\mathbb{R}^3`$은 자기 자신의 부분공간이다.

단일 벡터 $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$은 $`\mathbb{R}^3`$의 부분공간이다.

---

**예 1.** $`\mathbb{R}^2`$는 벡터 공간이다. 제1사분면은 부분공간인가?

```math
\mathcal{U} = \left\lbrace\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0, y \geq 0 \right\rbrace
```

$`\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}`$, $`c = -1`$을 취하자.

그러면 $`c\mathbf{u} = \begin{pmatrix} -2 \\ -3 \end{pmatrix} \notin \mathcal{U}`$.

이는 규칙 ii)를 위반한다. $`\mathcal{U}`$는 $`\mathbb{R}^2`$의 부분공간이 **아니다**.

---

**예 2.** $`\mathcal{U} = \left\lbrace\begin{pmatrix} x \\ y \end{pmatrix} \;\middle|\; x \geq 0 \text{ 이고 } y \geq 0, \text{ 또는 } x \leq 0 \text{ 이고 } y \leq 0 \right\rbrace`$

$`\mathcal{U}`$는 부분공간인가?

$`\mathbf{u} = \begin{pmatrix} 2 \\ 3 \end{pmatrix} \in \mathcal{U}`$, $`\mathbf{w} = \begin{pmatrix} -3 \\ -2 \end{pmatrix} \in \mathcal{U}`$를 취하자.

하지만 합 $`\mathbf{u} + \mathbf{w} = \begin{pmatrix} -1 \\ 1 \end{pmatrix} \notin \mathcal{U}`$.

$`\mathcal{U}`$는 부분공간이 **아니다**.

---

**예 3.** $`\mathbb{M}`$은 $`2 \times 2`$ 행렬 $`\begin{pmatrix} a & b \\ c & d \end{pmatrix}`$의 벡터 공간이다.

$`\mathcal{U}`$는 모든 상삼각 행렬 (Upper Triangular Matrix) $`\begin{pmatrix} a & b \\ 0 & d \end{pmatrix}`$의 집합이다.

$`\mathbb{D}`$는 모든 대각 행렬 (Diagonal Matrix) $`\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix}`$의 집합이다.

**$`\mathcal{U}`$와 $`\mathbb{D}`$ 모두 $`\mathbb{M}`$의 부분공간이다:**

- $`A, B \in \mathcal{U} \implies cA + dB \in \mathcal{U}`$
- $`A, B \in \mathbb{D} \implies cA + dB \in \mathbb{D}`$

참고:

```math
\begin{pmatrix} a & b \\ c & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + c\begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
```

### 2.12 열공간과 행공간

**$`A`$의 열공간 (Column Space):**

```math
A\mathbf{x} = \mathbf{b}
```

```math
\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix}
```

```math
= x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n
```

**$`A`$의 열공간** 은 열벡터들로 구성된 벡터 공간이다.

$`A\mathbf{x} = \mathbf{b}`$를 푸는 것은 $`\mathbf{b}`$를 열들의 선형 결합으로 표현하는 것이다. 우변 벡터 $`\mathbf{b}`$는 $`A`$의 열공간에 속해야 한다.

> **방정식 $`A\mathbf{x} = \mathbf{b}`$가 풀리려면 $`\mathbf{b}`$가 $`A`$의 열공간에 속해야 한다.**

**$`A`$의 행공간 (Row Space):**

$`A`$의 행들은 $`A^T`$의 열들이다.

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} \implies A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}
```

$`A^T`$는 $`m`$개의 열벡터를 가진다.

**$`A`$의 행공간은 $`A^T`$의 열공간이다.**

방정식 $`A^T \mathbf{y} = \mathbf{c}`$가 풀리려면 $`\mathbf{c}`$가 $`A^T`$의 열공간에 속해야 한다 ($`= \mathbf{c}`$가 $`A`$의 행공간에 속해야 한다).

**예:** 계수 1 행렬 (Rank 1 Matrix) $`A = \mathbf{u}\mathbf{v}^T`$를 고려하자:

```math
A = \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix} \begin{pmatrix} v_1 & v_2 & \cdots & v_n \end{pmatrix} = \begin{pmatrix} v_1\mathbf{u} & v_2\mathbf{u} & \cdots & v_n\mathbf{u} \end{pmatrix}
```

$`C(A)`$는 $`A`$의 모든 열벡터의 직선: $`c\mathbf{u}`$.

$`A^T = \mathbf{v}\mathbf{u}^T = \begin{pmatrix} u_1\mathbf{v} & u_2\mathbf{v} & \cdots & u_m\mathbf{v} \end{pmatrix}`$

$`C(A^T)`$는 $`A^T`$의 모든 열벡터의 직선: $`c\mathbf{v}`$.

### 2.13 생성 (Spanning)

**$`A`$의 열들이 벡터 공간 $`C(A)`$를 생성 (Span)한다.**

$`\mathbb{S}`$를 $`\mathbb{R}^m`$의 벡터 집합이라 하자. $`\mathbb{S} = \left\lbrace \begin{pmatrix} u_1 \\ u_2 \\ \vdots \\ u_m \end{pmatrix}, \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_m \end{pmatrix} \right\rbrace`$이면, $`\mathbb{S}`$는 $`\mathbb{R}^m`$의 부분공간이 **아니다**. 왜냐하면 $`\mathbf{u} + \mathbf{v} \notin \mathbb{S}`$이기 때문이다.

$`\mathbb{S}`$에 속한 벡터들의 모든 결합을 포함하면, 벡터 공간 $`V`$를 얻는다.

**집합 $`\mathbb{S}`$가 $`V`$를 생성한다.** $`V`$는 $`\mathbb{S}`$를 포함하는 가장 작은 벡터 공간이다.

행렬 $`A \in \mathbb{R}^{m \times n}`$을 고려하자:
- $`n`$개의 열이 열공간 $`C(A)`$를 생성한다
- $`A^T`$의 $`m`$개의 열이 행공간 $`C(A^T)`$를 생성한다

---

<br>

## 3. 소거법을 이용한 영공간 계산 (3.2)

### 3.1 핵심 사항: A = CR

```math
A = CR
```

1. $`\mathbb{R}^n`$에서 영공간 $`N(A)`$는 $`A\mathbf{x} = \mathbf{0}`$의 모든 해 $`\mathbf{x}`$를 포함한다. 여기에는 $`\mathbf{x} = \mathbf{0}`$도 포함된다.
2. $`A`$에서 $`R_0`$, $`R`$로의 소거는 영공간을 **바꾸지 않는다**.
3. 기약행사다리꼴 (Reduced Row Echelon Form) $`R_0 = \text{rref}(A)`$는 $`r`$개의 열에 $`I`$를 갖고, $`n - r`$개의 열에 $`F`$를 갖는다.
4. 열 $`j`$가 이전 열들에 종속이면, $`A\mathbf{x} = \mathbf{0}`$은 $`x_j = 1`$인 "특수해 (Special Solution)"를 가진다.
5. $`A\mathbf{x} = \mathbf{0}`$의 $`n - r`$개의 특수해는 $`-F`$와 $`I`$를 포함한다.
6. $`m < n`$인 모든 짧고 넓은 행렬은 영공간에 $`A\mathbf{x} = \mathbf{0}`$의 영이 아닌 해를 가진다.

### 3.2 Ax = 0의 모든 해 구하기

$`A\mathbf{x} = \mathbf{0}`$의 모든 해를 구하고 싶다.

$`A \in \mathbb{R}^{n \times n}`$이 가역이면 ($`\text{rank}(A) = n`$), 유일한 해는 $`\mathbf{x} = \mathbf{0}`$이다. $`A`$의 영공간은 영벡터만 포함한다: $`N(A) = \lbrace\mathbf{0}\rbrace`$.

일반적으로는 $`\text{rank}(A) = r`$이다. 즉, $`A`$는 $`r`$개의 독립 열을 가진다. 나머지 $`n - r`$개의 종속 열은 독립 열들의 결합이다. $`N(A)`$에서 $`n - r`$개의 벡터를 찾을 것인데, 이들이 $`A\mathbf{x} = \mathbf{0}`$의 특수해이다.

제2장에서 가역 행렬 $`A`$는 상삼각 행렬 (Upper Triangular Matrix) $`U`$로 축소되었다. $`A \in \mathbb{R}^{m \times n}`$에 대해, $`A\mathbf{x} = \mathbf{0}`$은 $`R\mathbf{x} = \mathbf{0}`$ (사다리꼴, Echelon Form)으로 단순화된다.

### 3.3 기약행사다리꼴

이 절에서는 $`A \in \mathbb{R}^{m \times n}`$을 고려하고, $`A`$를 기약행사다리꼴 (Reduced Row Echelon Form)로 소거한다: $`R_0 = \text{rref}(A)`$.

$`R_0`$는 영행을 가질 수 있다. $`R_0`$의 모든 영행을 제거하면 $`R`$이 된다.

**예 1:** $`R = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = (I \quad F)`$

$`\text{rank}(R) = 2`$, $`n = 4`$, $`n - r = 2`$개의 종속 열.

```math
R\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_1 + 3x_3 + 5x_4 \\ x_2 + 4x_3 + 6x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}
```

특수해는 어떻게 구하는가? 두 변수를 고정하고 $`R\mathbf{x} = \mathbf{0}`$을 이용한다:

$`x_3 = 1, x_4 = 0`$으로 취하면: $`\Rightarrow x_1 + 3 = 0, \; x_2 + 4 = 0`$. $`\therefore \mathbf{s}_1 = \begin{pmatrix} -3 \\ -4 \\ 1 \\ 0 \end{pmatrix}`$

$`x_3 = 0, x_4 = 1`$로 취하면: $`\Rightarrow x_1 + 5 = 0, \; x_2 + 6 = 0`$. $`\therefore \mathbf{s}_2 = \begin{pmatrix} -5 \\ -6 \\ 0 \\ 1 \end{pmatrix}`$

### 3.4 특수해와 영공간의 기저

두 특수해 $`\mathbf{s}_1, \mathbf{s}_2`$는 $`R`$의 영공간에 속한다:

```math
\mathbf{s}_1, \mathbf{s}_2 \in N(R)
```

```math
R\mathbf{s}_1 = \mathbf{0}, \quad R\mathbf{s}_2 = \mathbf{0}
```

```math
\Rightarrow R(c_1 \mathbf{s}_1 + c_2 \mathbf{s}_2) = \mathbf{0}
```

$`\mathbf{s}_1`$과 $`\mathbf{s}_2`$의 선형 결합은 $`R`$의 영공간에 속한다. **특수해 $`\mathbf{s}_1`$과 $`\mathbf{s}_2`$는 영공간의 기저이다.**

---

**예 2:** $`R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix}`$ (영행이 있는 기약행사다리꼴)

```math
R_0 \mathbf{x} = \begin{pmatrix} x_1 + 7x_2 + 8x_4 \\ x_3 + 9x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`\text{rank}(R_0) = 2`$, $`n = 4`$, $`n - r = 2`$. 두 개의 특수해를 구하자.

$`x_1 = 0, x_2 = 1`$로 취하면:

$`7 + 8x_4 = 0 \Rightarrow x_4 = -7/8`$, $`x_3 = -9(-7/8) = 63/8`$

$`\Rightarrow \mathbf{x} = \begin{pmatrix} 0 \\ 1 \\ 63/8 \\ -7/8 \end{pmatrix} \Rightarrow \mathbf{s}_1 = \begin{pmatrix} 0 \\ 8 \\ 63 \\ -7 \end{pmatrix}`$

$`x_1 = 1, x_2 = 0`$으로 취하면:

$`1 + 8x_4 = 0 \Rightarrow x_4 = -1/8`$, $`x_3 = +9/8`$

$`\Rightarrow \mathbf{x} = \begin{pmatrix} 1 \\ 0 \\ 9/8 \\ -1/8 \end{pmatrix} \Rightarrow \mathbf{s}_2 = \begin{pmatrix} 8 \\ 0 \\ 9 \\ -1 \end{pmatrix}`$

$`x_2 = 1, x_4 = 0`$과 $`x_2 = 0, x_4 = 1`$로 취하여 $`x_1`$과 $`x_3`$을 구할 수도 있다.

### 3.5 영공간 행렬: (-F; I)의 열들

이 절에서는 $`A \in \mathbb{R}^{m \times n}`$을 고려하고, $`A`$를 기약행사다리꼴 $`R_0 = \text{rref}(A)`$로 소거한다.

$`R_0`$는 영행을 가질 수 있다. $`R_0`$의 모든 영행을 제거하면 $`R`$이 된다.

```math
R = (I \quad F)Q
```

여기서 $`Q`$는 순열 행렬 (Permutation Matrix)이다 (피벗 열이 앞에 오도록 열을 재배치).

**상기** $`R\mathbf{x} = \mathbf{0}`$, 즉 $`(I \quad F)\mathbf{x} = \mathbf{0}`$ (예 1):

```math
\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

자유 변수를 $`\begin{pmatrix} 1 \\ 0 \end{pmatrix}`$과 $`\begin{pmatrix} 0 \\ 1 \end{pmatrix}`$로 설정하면:

```math
I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 3 \\ 4 \end{pmatrix}
```

```math
I\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} + F\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \implies \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = -\begin{pmatrix} 5 \\ 6 \end{pmatrix}
```

동치로:

```math
\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} \begin{pmatrix} -3 & -5 \\ -4 & -6 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}
```

즉:

```math
\boxed{(I \quad F) \begin{pmatrix} -F \\ I \end{pmatrix} = O}
```

**$`(I \quad F)\mathbf{x} = \mathbf{0}`$의 두 특수해는 $`\begin{pmatrix} -F \\ I \end{pmatrix}`$의 열들이다.**

### 3.6 행렬 분해 A = CR과 N(A)

**$`A`$에서 $`\text{rref}(A)`$로의 소거: 기약행사다리꼴**

소거를 적용하여 $`A`$를 $`R_0`$로 축소한다. 그러면 $`R_0`$에서 "$`I`$"가 $`A`$의 독립 열들의 행렬 $`C`$를 찾아준다. $`R_0`$에서 영행을 제거하면 $`A = CR`$에서의 행 행렬 $`R`$이 된다.

제2장에서 $`A`$는 정사각이고 가역이었다:

- $`A\mathbf{x} = \mathbf{b}`$ $`\xrightarrow{a) \text{ 소거}}`$ $`U\mathbf{x} = \mathbf{c}`$ $`\xrightarrow{b) \text{ 후진 대입}}`$ $`\mathbf{x} = U^{-1}\mathbf{c}`$

계수 $`r`$인 임의의 행렬 $`A`$에 대해:

```math
A = CR = C(I \quad F)
```

소거는 $`I_{r \times r}`$ 단위 행렬에 도달할 때까지 계속된다.

---

**예 1:**

```math
A = \begin{pmatrix} 1 & 2 & 11 & 17 \\ 3 & 7 & 37 & 57 \end{pmatrix}
```

$`\xrightarrow{R_2 - 3R_1}`$ $`\begin{pmatrix} 1 & 2 & 11 & 17 \\ 0 & 1 & 4 & 6 \end{pmatrix}`$ $`\xrightarrow{R_1 - 2R_2}`$ $`\begin{pmatrix} 1 & 0 & 3 & 5 \\ 0 & 1 & 4 & 6 \end{pmatrix} = R`$

$`\text{rank}(A) = 2`$, $`n = 4`$, $`n - r = 2`$.

$`A = (W \quad H)`$이고, $`W = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}`$ (독립 열), $`H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix}`$ (종속 열).

$`R = (I \quad F)`$. 소거는 $`W`$를 역행렬로 만들었다. 이는 $`A`$에 $`W^{-1}`$을 곱하는 것과 같다:

```math
W^{-1}A = W^{-1}(W \quad H) = (I \quad W^{-1}H) = (I \quad F) = R
```

$`W^{-1}H = F \implies H = WF`$.

$`W`$는 독립 열로 구성되어 있다. $`F`$는 $`A`$의 독립 열들을 어떻게 결합하는지 알려준다.

```math
H = \begin{pmatrix} 11 & 17 \\ 37 & 57 \end{pmatrix} = WF = \begin{pmatrix} 1 & 2 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 3 & 5 \\ 4 & 6 \end{pmatrix}
```

---

**예 2:**

```math
A = \begin{pmatrix} 1 & 7 & 3 & 35 \\ 2 & 14 & 6 & 70 \\ 2 & 14 & 9 & 97 \end{pmatrix}
```

$`\xrightarrow{R_2 - 2R_1, \; R_3 - 2R_1}`$ $`\begin{pmatrix} 1 & 7 & 3 & 35 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 3 & 27 \end{pmatrix}`$ $`\xrightarrow{R_1 - R_3, \; R_3/3}`$ $`\begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 9 \end{pmatrix}`$ $`\xrightarrow{\text{swap } R_2 \text{ and } R_3}`$

```math
R_0 = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \\ 0 & 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q
```

여기서 $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}`$, $`F = \begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}`$이고, $`Q`$는 열을 재배치하는 순열 행렬이다.

영행을 제거하면: $`R = (I \quad F)Q`$.

열 순열 $`PAQ`$ 후:

$`PAQ = \begin{pmatrix} 1 & 3 & 7 & 35 \\ 2 & 9 & 14 & 97 \\ 2 & 6 & 14 & 70 \end{pmatrix}`$이고, $`W = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}`$, $`H = \begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix}`$

$`W^{-1}H = F \implies H = WF`$:

```math
\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}
```

**예 2 상기:** 단위 행렬이 $`R_0`$의 1번째와 3번째 열에 나타나므로, $`A`$의 1번째와 3번째 열이 독립이다. $`A`$의 2번째와 4번째 열은 독립 열들의 선형 결합이다:

```math
\begin{pmatrix} 7 & 35 \\ 14 & 97 \\ 14 & 70 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ 2 & 9 \\ 2 & 6 \end{pmatrix}\begin{pmatrix} 7 & 8 \\ 0 & 9 \end{pmatrix}
```

**행렬 분해 $`A = CR`$과 $`N(A)`$:**

소거를 적용하여 $`A`$를 $`R_0`$로 축소한다. 그러면 $`R_0`$에서 "$`I`$"가 $`A`$의 독립 열들의 행렬 $`C`$를 찾아준다. $`R_0`$에서 영행을 제거하면 $`A = CR`$에서의 행 행렬 $`R`$이 된다.

예 2의 두 특수해 $`\mathbf{s}_1, \mathbf{s}_2`$를 구하자:

```math
R = \begin{pmatrix} 1 & 7 & 0 & 8 \\ 0 & 0 & 1 & 9 \end{pmatrix}
```

```math
R\mathbf{s} = \begin{pmatrix} s_1 + 7s_2 + 8s_4 \\ s_3 + 9s_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

$`s_2 = 1, s_4 = 0`$으로 취하면: $`\Rightarrow s_3 = 0, s_1 = -7`$. $`\therefore \mathbf{s}_1 = \begin{pmatrix} -7 \\ 1 \\ 0 \\ 0 \end{pmatrix}`$

$`s_2 = 0, s_4 = 1`$로 취하면: $`\Rightarrow s_1 = -8, s_3 = -9`$. $`\therefore \mathbf{s}_2 = \begin{pmatrix} -8 \\ 0 \\ -9 \\ 1 \end{pmatrix}`$

---

<br>

## 4. Ax = b의 완전해 (3.3)

### 4.1 완전해의 구조

1. **$`A\mathbf{x} = \mathbf{b}`$의 완전해 (Complete Solution):**

```math
\mathbf{x} = \text{하나의 특수해 } \mathbf{x}_p + \text{영공간의 임의의 } \mathbf{x}_n
```

2. $`A\mathbf{x} = \mathbf{b}`$에 **소거** 를 적용하면 $`R_0 \mathbf{x} = \mathbf{d}`$가 된다. $`R_0`$의 영행에 대응하는 $`\mathbf{d}`$의 성분이 0일 때 풀 수 있다.

```math
\begin{pmatrix} R \\ 0 \; 0 \; \cdots \; 0 \\ 0 \; 0 \; \cdots \; 0 \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ 0 \\ 0 \end{pmatrix}
```

3. $`R_0 \mathbf{x} = \mathbf{d}`$가 풀릴 때, 하나의 특수해 $`\mathbf{x}_p`$는 **모든 자유 변수가 0** 인 해이다.

4. $`A`$가 **열 완전 계수 (Full Column Rank)** $`r = n`$일 때 영공간 $`N(A) = \lbrace\text{영벡터}\rbrace`$. 자유 변수가 없다.

5. $`A`$가 **행 완전 계수 (Full Row Rank)** $`r = m`$일 때 열공간 $`C(A)`$가 $`\mathbb{R}^m`$이다: $`A\mathbf{x} = \mathbf{b}`$는 항상 풀린다.

### 4.2 풀이 예제: 특수해 구하기

이전 절에서 $`A\mathbf{x} = \mathbf{0}`$의 해를 구했다. 이 절에서는 $`A\mathbf{x} = \mathbf{b}`$의 해를 구한다.

좌변의 행 연산은 우변에도 적용해야 한다. **첨가 행렬 (Augmented Matrix)** $`(A | \mathbf{b})`$를 사용한다.

```math
\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 1 & 3 & 1 & 6 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 7 \end{pmatrix} \quad \Rightarrow \quad (A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & 1 \\ 0 & 0 & 1 & 4 & | & 6 \\ 1 & 3 & 1 & 6 & | & 7 \end{pmatrix}
```

소거 후 ($`R_3 - R_1 - R_2`$):

```math
\begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} \quad \Rightarrow \quad (R_0 | \mathbf{d})
```

마지막 방정식은 $`0 = 0`$이다. 일반적인 $`\mathbf{b}`$를 고려하자:

```math
(A | \mathbf{b}) = \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 1 & 3 & 1 & 6 & | & b_3 \end{pmatrix} \xrightarrow{R_3 - R_1 - R_2} \begin{pmatrix} 1 & 3 & 0 & 2 & | & b_1 \\ 0 & 0 & 1 & 4 & | & b_2 \\ 0 & 0 & 0 & 0 & | & b_3 - b_1 - b_2 \end{pmatrix} = (R_0 | \mathbf{d})
```

$`b_3 - b_1 - b_2 = 0`$이면 세 번째 방정식에서 $`0 = 0`$을 얻을 수 있다.

**하나의 특수해 $`A\mathbf{x}_p = \mathbf{b}`$:**

```math
R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix} = \mathbf{d}
```

$`\text{rank}(R_0) = 2`$, $`n = 4`$, $`n - r = 2`$개의 자유 변수.

**$`x_2 = 1, x_4 = 0`$으로 취하면:**

$`x_1 + 3 = 1 \Rightarrow x_1 = -2`$, $`x_3 = 6`$

$`\therefore \mathbf{x} = \begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix}`$

**$`x_2 = 0, x_4 = 1`$로 취하면:**

$`x_1 + 2 = 1 \Rightarrow x_1 = -1`$, $`x_3 + 4 = 6 \Rightarrow x_3 = 2`$

$`\therefore \mathbf{x} = \begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix}`$

**$`x_2 = 0, x_4 = 2`$로 취하면:**

$`x_1 + 2 \cdot 2 = 1 \Rightarrow x_1 = -3`$, $`x_3 + 4 \cdot 2 = 6 \Rightarrow x_3 = -2`$

$`\therefore \mathbf{x} = \begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix}`$

자유 변수로 인해 무한히 많은 해를 구할 수 있다.

### 4.3 완전해 분해

해는 **특수해 (Particular Solution)** 와 영공간의 **특수해들 (Special Solutions)** 로 분해할 수 있다.

```math
R_0 \mathbf{x} = \begin{pmatrix} 1 & 3 & 0 & 2 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`R_0 = \begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix} Q`$이고, $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}`$, $`F = \begin{pmatrix} 3 & 2 \\ 0 & 4 \end{pmatrix}`$

영공간 벡터는 $`Q^T \begin{pmatrix} -F \\ I \end{pmatrix}`$의 열들이다:

```math
\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} -3 & -2 \\ 0 & -4 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} -3 & -2 \\ 1 & 0 \\ 0 & -4 \\ 0 & 1 \end{pmatrix}
```

**$`R_0 \mathbf{x} = \mathbf{d}`$의 해 재검증:**

```math
\begin{pmatrix} -2 \\ 1 \\ 6 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix}
```

```math
\begin{pmatrix} -1 \\ 0 \\ 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}
```

```math
\begin{pmatrix} -3 \\ 0 \\ -2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + \begin{pmatrix} -4 \\ 0 \\ -8 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + 2\begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}
```

**Q. 특수해는 어떻게 구하는가?**

$`x_2 = 0, x_4 = 0`$ (모든 자유 변수를 0으로 설정):

```math
R_0 \mathbf{x} = \begin{pmatrix} x_1 + 3x_2 + 2x_4 \\ x_3 + 4x_4 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 6 \\ 0 \end{pmatrix}
```

$`\Rightarrow x_1 = 1, \; x_3 = 6`$

```math
\therefore \mathbf{x}_p = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix}
```

**$`A\mathbf{x} = \mathbf{b}`$의 완전해** $`\mathbf{x}_p + \mathbf{x}_n`$:

```math
\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 1 \\ 0 \\ 6 \\ 0 \end{pmatrix} + x_2 \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_4 \begin{pmatrix} -2 \\ 0 \\ -4 \\ 1 \end{pmatrix}
```

여기서 $`\mathbf{x}_p`$가 특수해이고 나머지 항들은 영공간 벡터이다.

### 4.4 열 완전 계수: r = n

**Q. $`m = n = r`$일 때 $`\mathbf{x}_p, \mathbf{x}_n`$은?**

**A.** $`\mathbf{x}_n = \mathbf{0}`$.

$`A\mathbf{x} = \mathbf{b} \iff A(\mathbf{x}_p + \mathbf{x}_n) = \mathbf{b} \implies A(\mathbf{x}_p + \mathbf{0}) = \mathbf{b} \implies \mathbf{x}_p = A^{-1}\mathbf{b}`$

### 4.5 풀림 조건

**예 1.** $`A\mathbf{x} = \mathbf{b}`$가 풀리기 위한 $`(b_1, b_2, b_3)`$의 조건을 구하라:

```math
A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ -2 & -3 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \\ b_3 \end{pmatrix}
```

```math
(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & | & b_1 \\ 1 & 2 & | & b_2 \\ -2 & -3 & | & b_3 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 + 2R_1} \begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & -1 & | & b_3 + 2b_1 \end{pmatrix} \xrightarrow{R_3 + R_2}
```

```math
\begin{pmatrix} 1 & 1 & | & b_1 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & | & 2b_1 - b_2 \\ 0 & 1 & | & b_2 - b_1 \\ 0 & 0 & | & b_3 + b_2 + b_1 \end{pmatrix} = (R_0, \mathbf{d})
```

$`b_3 + b_2 + b_1 = 0`$이면, $`A\mathbf{x} = \mathbf{b}`$가 풀린다. 즉, $`b_3 + b_2 + b_1 = 0`$이 $`\mathbf{b}`$를 $`A`$의 열공간에 넣기 위한 조건이다.

$`\text{rank}(A) = 2`$, $`n = 2`$, $`n - r = 0`$. 자유 변수 없음. $`N(A) = \lbrace\mathbf{0}\rbrace`$, $`\mathbf{x}_n = \mathbf{0}`$.

```math
\mathbf{x}_p = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix}
```

```math
\therefore \mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2b_1 - b_2 \\ b_2 - b_1 \end{pmatrix} + \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

$`b_3 + b_2 + b_1 \neq 0`$이면, $`A\mathbf{x} = \mathbf{b}`$의 해가 없다.

**열 완전 계수 $`r = n`$일 때:**

```math
R_0 = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix}
```

1. 행렬 $`A`$는 $`n`$개의 독립 열을 가진다.
2. $`A`$의 영공간은 $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$이다.
3. $`A\mathbf{x} = \mathbf{b}`$에 해가 있으면, **유일한** 해를 가진다.

```math
R_0 \mathbf{x} = \begin{pmatrix} I_{n \times n} \\ O_{(m-n) \times n} \end{pmatrix} \mathbf{x} = \begin{pmatrix} \vdots \\ O_{(m-n) \times 1} \end{pmatrix}
```

하단 $`m - n`$개의 행이 $`\mathbf{b}`$가 $`A`$의 열공간에 속하기 위한 $`m - n`$개의 조건을 제공한다.

### 4.6 행 완전 계수와 완전해

$`A \in \mathbb{R}^{m \times n}`$, $`m \leq n`$ (짧고 넓은 행렬).

행렬이 **행 완전 계수 (Full Row Rank)** 를 가지려면 $`r = m`$.

**예 2.** $`A\mathbf{x} = \mathbf{b}`$가 $`n = 3`$개의 미지수를 갖지만 $`m = 2`$개의 방정식만 있을 때:

```math
x + y + z = 3, \quad x + 2y - z = 4
```

```math
(A | \mathbf{b}) = \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 1 & 2 & -1 & | & 4 \end{pmatrix} \xrightarrow{R_2 - R_1} \begin{pmatrix} 1 & 1 & 1 & | & 3 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} \xrightarrow{R_1 - R_2} \begin{pmatrix} 1 & 0 & 3 & | & 2 \\ 0 & 1 & -2 & | & 1 \end{pmatrix} = (R | \mathbf{d})
```

$`\text{rank}(A) = 2`$, $`n = 3`$, $`n - r = 1`$. 1개의 자유 변수, 1개의 특수해.

**i) $`\mathbf{x}_p`$:**

```math
\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}
```

```math
x + 3z = 2, \quad y - 2z = 1
```

$`z = 0`$으로 취하면: $`\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$, $`\mathbf{x}_p = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix}`$ (특수해).

**ii) $`\mathbf{x}_n`$:**

```math
\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -2 \end{pmatrix} \begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

$`z = 1`$로 취하면: $`x + 3 = 0 \Rightarrow x = -3`$, $`y - 2 = 0 \Rightarrow y = 2`$.

$`\therefore \mathbf{x}_n = \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}`$

**iii) 완전해:**

```math
\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n = \begin{pmatrix} 2 \\ 1 \\ 0 \end{pmatrix} + \alpha \begin{pmatrix} -3 \\ 2 \\ 1 \end{pmatrix}
```

이것은 $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$을 지나는 방향의 **직선** 이며 ($`A\mathbf{x}_n = \mathbf{0}`$), $`\mathbf{x}_p`$만큼 이동한 것이다 ($`A\mathbf{x} = \mathbf{b}`$).

**행 완전 계수 ($`r = m`$)인 모든 행렬 $`A`$는 다음 성질을 가진다:**
1. 모든 행에 피벗이 있고, $`R_0`$에 영행이 없다. $`R_0 = R`$.
2. $`A\mathbf{x} = \mathbf{b}`$는 모든 $`\mathbf{b}`$에 대해 해를 가진다.
3. $`A`$의 열공간은 전체 공간 $`\mathbb{R}^m`$이다.

$`A\mathbf{x} = \mathbf{b} \to R\mathbf{x} = \mathbf{d}`$이고, $`R = (I_{m \times m} \quad F_{m \times (n-m)})`$.

4. $`m < n`$이면, $`A\mathbf{x} = \mathbf{b}`$는 **부정 (Underdetermined)** 이다 (많은 해).

행 완전 계수 ($`r = m`$)일 때, $`m`$개의 행이 선형 독립이다. $`A^T`$의 열들이 LI이다. $`A^T`$의 영공간은 $`\mathbb{Z} = \lbrace\mathbf{0}\rbrace`$이다.

### 4.7 선형 방정식의 네 가지 경우

선형 방정식의 네 가지 경우는 계수 $`r`$에 따라 달라진다:

| 경우 | 형태 | $`R_0`$ 형식 | $`A\mathbf{x} = \mathbf{b}`$의 해 |
|:-----|:------|:-----------|:----------------------------------------|
| $`r = m`$ 이고 $`r = n`$ | 정사각, 가역 | $`(I)`$ | 1개의 해 |
| $`r = m`$ 이고 $`r < n`$ | 짧고 넓음 | $`(I \quad F)`$ | $`\infty`$개의 해 |
| $`r < m`$ 이고 $`r = n`$ | 키 크고 좁음 | $`\begin{pmatrix} I \\ 0 \end{pmatrix}`$ | 0 또는 1개의 해 |
| $`r < m`$ 이고 $`r < n`$ | 완전 계수 아님 | $`\begin{pmatrix} I & F \\ 0 & 0 \end{pmatrix}`$ | 0 또는 $`\infty`$개의 해 |

---

<br>
