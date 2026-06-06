# 제1장 강의 — 벡터와 행렬 — Part 2/2

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.** GitHub은 파일 하나당 렌더되는 수식 개수에 한도가 있어, 모든 수식이 제대로 보이도록 절 단위로 나눴습니다.
>
> [Part 1](Concepts.ko.md) · **Part 2**

---

<br>

## 목차

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

## 1.3 행렬과 열공간 (Column Spaces)

### 1.3.1 행렬-벡터 곱셈

**(1)** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}`$는 **3 x 2 행렬** 이다: 3개의 행과 2개의 열, 랭크 2.

**(2)** $`A\mathbf{x}`$의 3개의 성분은 $`A`$의 3개의 행과 벡터 $`\mathbf{x}`$의 내적이다.

**예시:** $`\mathbf{x} = \begin{pmatrix} 7 \\ 8 \end{pmatrix}`$

```math
A\mathbf{x} = \begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 7 \\ 8 \end{pmatrix} = \begin{pmatrix} 7 + 16 \\ 21 + 32 \\ 35 + 48 \end{pmatrix} = \begin{pmatrix} 23 \\ 53 \\ 83 \end{pmatrix}
```

**(3)** $`A\mathbf{x}`$는 $`A`$의 **열들의 결합** 이다:

```math
\begin{pmatrix} 1 & 2 \\ 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 7 \\ 8 \end{pmatrix} = 7\begin{pmatrix} 1 \\ 3 \\ 5 \end{pmatrix} + 8\begin{pmatrix} 2 \\ 4 \\ 6 \end{pmatrix}
```

**(4)** $`A`$의 열공간은 열들의 모든 결합 $`A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2`$를 포함한다.

**(5)** 랭크 1 행렬: $`A`$의 모든 열이 하나의 직선 위에 있다.

### $`m \times n`$ 행렬 살펴보기

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}
```

$`m = n`$이면: **정방행렬 (Square matrix)**.

**정방행렬의 예시 ($`A \in \mathbb{R}^{3 \times 3}`$):**

```math
\text{항등행렬 (Identity): } \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \qquad \text{대각행렬 (Diagonal): } \begin{pmatrix} 2 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 0 & 5 \end{pmatrix}
```

```math
\text{삼각행렬 (Triangular): } \begin{pmatrix} 2 & 1 & -3 \\ 0 & 4 & 7 \\ 0 & 0 & 5 \end{pmatrix} \qquad \text{대칭행렬 (Symmetric): } \begin{pmatrix} 2 & 1 & -3 \\ 1 & 4 & 7 \\ -3 & 7 & 5 \end{pmatrix}
```

### 열을 벡터로 해석하기

$`A`$의 열들을 벡터로 해석할 수 있다:

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}
```

여기서 $`\mathbf{a}_i \in \mathbb{R}^m`$, $`i = 1, 2, \ldots, n`$.

**질문:** $`m \times n`$ 행렬 $`A`$는 $`n \times 1`$ 벡터 $`\mathbf{x}`$를 어떻게 곱하는가?

1. $`\mathbf{x}`$와 $`A`$의 행들의 **내적**
2. $`A`$의 열들의 **선형결합**

### $`A\mathbf{x}`$의 두 가지 관점

**(1) 행 관점 (내적):**

```math
A\mathbf{x} = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}
```

3개의 행 $`\Rightarrow`$ 3개의 내적.

$`A`$의 각 행은 벡터 $`\mathbf{x}`$와 같은 수의 성분을 가진다.

**(2) 열 관점 (선형결합):**

```math
A\mathbf{x} = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2 + \cdots + x_n\mathbf{a}_n
```

$`A`$의 **열들의 결합**.

**예시:** $`A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}`$

```math
A\mathbf{x} = x_1\begin{pmatrix} -1 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix} + x_4\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \\ x_4 - x_3 \end{pmatrix}
```

### 1.3.2 열공간 (Column Space)

$`A = (\mathbf{a}_1 \ \mathbf{a}_2)`$이고 $`\mathbf{a}_1`$과 $`\mathbf{a}_2`$가 선형독립일 때:

```math
A\mathbf{x} = x_1\mathbf{a}_1 + x_2\mathbf{a}_2
```

는 **모든 결합의 평면** (생성, Span)을 보여준다.

### 1.3.3 독립, 종속, 그리고 열공간

**예시 (ex1):**

```math
A_1 = \begin{pmatrix} 1 & 0 & 0 \\ 2 & 4 & 0 \\ 3 & 5 & 6 \end{pmatrix}
```

각 열이 새로운 방향을 제시한다. 이들의 결합은 **3차원 공간** 을 채운다.

```math
A_1\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = x_1\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix} + x_2\begin{pmatrix} 0 \\ 4 \\ 5 \end{pmatrix} + x_3\begin{pmatrix} 0 \\ 0 \\ 6 \end{pmatrix}
```

**독립 열:** $`A_1\mathbf{x} = \mathbf{0}`$이면 $`\mathbf{x} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$일 때만 성립

**예시 (ex2):**

```math
A_2 = \begin{pmatrix} 1 & 2 & 3 \\ 1 & 4 & 5 \\ 6 & 0 & 6 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

$`\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2`$. $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3`$의 결합은 3차원 공간을 채우지 **못한다**. $`\mathbf{a}_3`$는 3차원에서 $`\mathbf{a}_1`$과 $`\mathbf{a}_2`$의 평면 위에 놓인다.

**예시 (ex3):**

```math
A_3 = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 6 & 8 \\ 5 & 15 & 20 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

$`\mathbf{a}_2 = 3\mathbf{a}_1`$이고 $`\mathbf{a}_3 = 4\mathbf{a}_1`$.

$`A_3`$의 세 열 모두 3차원 공간에서 **같은 직선** 위에 놓인다.

### 열공간 $`C(A)`$

열공간 $`C(A)`$는 모든 벡터 $`A\mathbf{x}`$를 포함한다: 열들의 모든 결합.

**$`A`$의 열공간에 대한 고찰:**

```math
A_4 = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix} \quad \text{(독립 열)}
```

```math
A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix} \quad \text{4번째 열은 종속인가?}
```

$`A_5`$에 대해: $`\mathbf{a}_4 = \mathbf{a}_1 - \mathbf{a}_2 + \mathbf{a}_3`$.

$`\mathbf{v} = A_4\mathbf{x}`$를 생각하자:

```math
\mathbf{v} = x_1\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + x_2\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + x_3\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + x_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
```

$`A_4\mathbf{x} = \mathbf{v}`$ 풀기: 방정식으로부터,

```math
v_4 = x_4, \quad v_3 = x_4 + x_3, \quad v_2 = x_4 + x_3 + x_2, \quad v_1 = x_4 + x_3 + x_2 + x_1
```

따라서: $`x_4 = v_4`$, $`x_3 = v_3 - v_4`$, $`x_2 = v_2 - v_3`$, $`x_1 = v_1 - v_2`$

```math
\therefore \mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \\ v_4 \end{pmatrix} = (v_1 - v_2)\begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix} + (v_2 - v_3)\begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} + (v_3 - v_4)\begin{pmatrix} 1 \\ 1 \\ 1 \\ 0 \end{pmatrix} + v_4\begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
```

이는 **모든 $`\mathbf{v}`$가 열공간에 속한다** 는 것을 의미한다. 네 개의 방정식 $`A_4\mathbf{x} = \mathbf{v}`$가 풀렸다.

### 1.3.4 생성 (Span)

**SPAN (생성)** 은 벡터 집합의 모든 선형결합을 기술한다.

**$`A`$의 열들의 생성은 열공간이다.**

**예시:** $`\mathbf{a}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}`$, $`\mathbf{a}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix} \in \mathbb{R}^4`$

- $`c_1\mathbf{a}_1`$은 4차원에서 하나의 직선에 해당한다.
- $`c_2\mathbf{a}_2`$는 4차원에서 하나의 직선에 해당한다.
- $`c_1\mathbf{a}_1 + c_2\mathbf{a}_2`$는 4차원 공간에서 2차원 평면을 채운다.
- 이 평면은 열 $`\mathbf{a}_1`$과 $`\mathbf{a}_2`$의 생성이다.

**예시:** $`A_4`$는 4개의 독립 열을 가진다. $`A_4`$의 열공간은 $`\mathbb{R}^4`$ 전체이다.

**예시:** $`A_5`$는 하나의 종속 열이 있다 (독립 열은 세 개뿐이다).

```math
A_5 = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 1 & 1 \\ 1 & 0 & 0 & 1 \end{pmatrix}
```

$`A_5`$의 열공간은 $`\mathbb{R}^4`$ 내의 **3차원 부분공간** 이다. 4번째 열은 그 부분공간에 속한다.

$`A_5\mathbf{x} = \mathbf{v}`$는 $`\mathbf{v} \in C(A_5)`$일 때만 풀 수 있다.

### 1.3.5 랭크와 기저 (Rank and Basis)

$`A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix} \in \mathbb{R}^{m \times n}`$로 놓자

$`C(A)`$는 $`A`$의 열공간, $`\mathbf{a}_i \in \mathbb{R}^m`$.

열공간 $`C(A)`$는 $`\mathbb{R}^m`$ 전체를 채울 수도 있고 그렇지 않을 수도 있다.

**예시:** $`m = 3`$인 경우. $`A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}`$

$`C(A)`$는:
- 전체 공간 $`\mathbb{R}^3`$ — $`A`$가 3개의 독립 열 (I.C.)을 가질 때
- $`\mathbb{R}^3`$에서의 평면 — $`A`$가 2개의 I.C.를 가질 때
- $`\mathbb{R}^3`$에서의 직선 — $`A`$가 1개의 I.C.를 가질 때
- 단일 점 $`\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$ — $`A`$가 영행렬일 때

**더 많은 예시 ($`3 \times 3`$ 행렬):**

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}: C(A) = \mathbb{R}^3
```

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}: C(A) = xy \text{ 평면}
```

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}: C(A) = x \text{ 축}
```

**핵심 질문:**

- **Q1:** $`A`$의 몇 개의 열이 독립인가? 그 수 $`r`$이 $`A`$의 **"랭크 (Rank)"** 이다.
- **Q2:** 처음 $`r`$개의 독립 열은 무엇인가? 이것이 열공간의 **"기저 (Basis)"** 이다.
- **Q3:** 그 $`r`$개의 기저 벡터의 어떤 결합이 나머지 $`n - r`$개의 열을 만드는가?
- **Q4:** 임의의 $`A`$를 $`m \times r`$ 열 행렬 $`C`$ 곱하기 $`r \times n`$ 행렬 $`R`$로 쓸 수 있다: $`A = CR`$.
- **Q5:** $`R`$의 $`r`$개의 행은 $`A`$의 **행공간 (Row space)** 의 기저이다. $`R`$의 행은 $`A`$에서 직접 오지 않는다.

### 1.3.6 랭크 1인 행렬

랭크 1 행렬에서는 **모든 열벡터가 같은 직선을 따라 놓인다**.

**예시:**

```math
A_6 = \begin{pmatrix} 1 & 3 & -2 \\ 4 & 12 & -8 \\ 2 & 6 & -4 \end{pmatrix}
```

$`\mathbf{a}_2 = 3\mathbf{a}_1`$, $`\mathbf{a}_3 = -2\mathbf{a}_1`$.

- $`A_6`$의 랭크 $`r = 1`$.
- $`C(A_6) = c_1\mathbf{a}_1`$, 원점을 지나는 직선.
- $`A_6`$의 모든 행은 하나의 행의 배수이다.

> 열공간이 $`\mathbb{R}^m`$에서의 단일 직선이면, 행공간은 $`\mathbb{R}^n`$에서의 단일 직선이다.

**질문: 왜 이렇게 되는가?**

모든 열이 같은 방향이면, 모든 행도 같은 방향이다.

**예시:** $`A = \begin{pmatrix} a & ma \\ b & mb \end{pmatrix}`$

행 2는 행 1의 배수이다. 열 랭크가 1이면 행 랭크도 1이다.

**예시:** $`A = \begin{pmatrix} a & ma & pa \\ b & mb & pb \\ c & mc & pc \end{pmatrix}`$

행 2 = $`\frac{b}{a}`$ 행 1, 행 3 = $`\frac{c}{a}`$ 행 1.

열 랭크가 1이면 행 랭크도 1이다.

> **질문. 행 랭크와 열 랭크는 모든 행렬에서 같은가.** 답: **그렇다.**

---

<br>

## 1.4 행렬 곱셈 AB와 CR

### 1.4.1 행렬 곱셈의 규칙

**(1)** $`AB`$를 곱하려면: $`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$.

**$`A`$의 행 길이가 $`B`$의 열 길이와 같아야 한다.**

**(2) 내적 관점:**

```math
(AB)_{ij} = (A\text{의 행 } i) \cdot (B\text{의 열 } j)
```

**(3) 열 관점:**

```math
AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p)
```

$`AB`$의 열 $`j`$ = $`A\mathbf{b}_j`$

**(4) 비교환성:** 일반적으로 $`AB \neq BA`$.

**(5)** $`A`$가 $`C`$에서 $`r`$개의 독립 열을 가지면: $`A = CR`$ 여기서 $`C \in \mathbb{R}^{m \times r}`$, $`R \in \mathbb{R}^{r \times n}`$.

### 1.4.2 AB의 열 해석

$`B = \begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix}`$로 놓자

```math
AB = A\begin{pmatrix} | & | & & | \\ \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_p \\ | & | & & | \end{pmatrix} = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}
```

이것은 **$`A`$의 결합** 이다.

**예시:** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$, $`B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}`$

```math
AB = (A\mathbf{b}_1 \ A\mathbf{b}_2) = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}
```

```math
A\mathbf{b}_1 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}
```

```math
A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}
```

**예시:** $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$, $`B = \begin{pmatrix} 5 & 6 \\ 7 & 8 \end{pmatrix}`$

```math
A\mathbf{b}_1 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 5 \\ 7 \end{pmatrix}
```

**(1) 내적 접근법:**

```math
= \begin{pmatrix} 1 \cdot 5 + 2 \cdot 7 \\ 3 \cdot 5 + 4 \cdot 7 \end{pmatrix} = \begin{pmatrix} 19 \\ 43 \end{pmatrix}
```

**(2) 열 결합 접근법:**

```math
= 5\begin{pmatrix} 1 \\ 3 \end{pmatrix} + 7\begin{pmatrix} 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 19 \\ 43 \end{pmatrix}
```

두 접근법 모두 "4"번의 곱셈을 사용한다.

```math
A\mathbf{b}_2 = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} 6 \\ 8 \end{pmatrix} = \begin{pmatrix} 6 + 16 \\ 18 + 32 \end{pmatrix} = \begin{pmatrix} 22 \\ 50 \end{pmatrix}
```

```math
\therefore AB = \begin{pmatrix} 19 & 22 \\ 43 & 50 \end{pmatrix}
```

### 1.4.3 계산 비용

$`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$로 놓자.

```math
AB = \begin{pmatrix} | & | & & | \\ A\mathbf{b}_1 & A\mathbf{b}_2 & \cdots & A\mathbf{b}_p \\ | & | & & | \end{pmatrix}
```

**(1) 내적 횟수:** $`A\mathbf{b}_1 \Rightarrow`$ $`m`$개의 내적. $`AB \Rightarrow`$ $`mp`$개의 내적. 각 내적은 $`n`$번의 곱셈 $`\Rightarrow`$ 총 **$`mnp`$번의 곱셈**.

**(2) 열 결합 횟수:** $`A\mathbf{b}_1 = b_{11}\mathbf{a}_1 + b_{12}\mathbf{a}_2 + \cdots + b_{1n}\mathbf{a}_n \Rightarrow`$ 각각 $`m`$번의 곱셈, 하나의 열에 $`mn`$번의 곱셈. $`AB \Rightarrow`$ **$`mnp`$번의 곱셈**.

> 행렬 곱셈의 비용은 $`mnp`$이다.
>
> $`A, B \in \mathbb{R}^{n \times n}`$이면, $`AB`$의 비용은 $`n^3`$이다.

### 1.4.4 행렬 곱셈의 성질

**비교환성:**

```math
AB \neq BA \quad \text{일반적으로.}
```

행렬 곱셈은 교환법칙이 성립하지 **않는다**.

**결합법칙:**

```math
(AB)C = A(BC)
```

행렬 곱셈에서 결합법칙은 **성립한다**.

**분배법칙:**

```math
A(B + C) = AB + AC
```

### 1.4.5 랭크 1 행렬과 A = CR

**$`AB`$ 복습:**

**(1) 내적:** $`(AB)_{ij} = (A\text{의 행 } i) \cdot (B\text{의 열 } j)`$. $`(m \times n)(n \times p) \Rightarrow mp`$개의 내적.

**(2) 열 결합:** $`AB = A(\mathbf{b}_1 \ \mathbf{b}_2 \ \cdots \ \mathbf{b}_p) = (A\mathbf{b}_1 \ A\mathbf{b}_2 \ \cdots \ A\mathbf{b}_p)`$. $`AB`$의 열 $`j`$ = $`A\mathbf{b}_j`$.

**랭크 1 행렬과 $`A = CR`$:**

랭크 1 행렬의 모든 열은 같은 직선 위에 놓인다.

$`A`$의 모든 열이 같은 열 방향에 있으면, $`A`$의 모든 행도 같은 행 방향에 있다.

**예시:** 랭크 $`r = 1`$.

```math
A = \begin{pmatrix} 1 & 2 & 10 & 100 \\ 3 & 6 & 30 & 300 \\ 2 & 4 & 20 & 200 \end{pmatrix}
```

하나의 독립 열, 하나의 독립 행.

$`A`$의 열공간이 직선이면, $`A`$의 행공간도 직선이다.

$`A`$를 다음과 같이 분해할 수 있다:

```math
A = \begin{pmatrix} 1 \\ 3 \\ 2 \end{pmatrix}\begin{pmatrix} 1 & 2 & 10 & 100 \end{pmatrix} = CR
```

여기서 $`C \in \mathbb{R}^{3 \times 1}`$이고 $`R \in \mathbb{R}^{1 \times 4}`$.

### 1.4.6 C와 R 구하기

**$`C`$는 $`A`$의 처음 $`r`$개의 독립 열을 포함한다.**

주어진 $`A`$에서 **왼쪽에서 오른쪽으로** 독립 열을 찾는다:

```math
A = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}
```

1. $`\mathbf{a}_1 \neq \mathbf{0}`$이면, $`\mathbf{a}_1`$을 $`C`$에 넣는다.
2. $`\mathbf{a}_2 \neq c_1\mathbf{a}_1`$ ($`\mathbf{a}_1`$의 배수가 아니면), $`\mathbf{a}_2`$를 $`C`$에 넣는다.
3. $`\mathbf{a}_3 \neq c_1\mathbf{a}_1 + c_2\mathbf{a}_2`$ ($`\mathbf{a}_1`$과 $`\mathbf{a}_2`$의 결합이 아니면), $`\mathbf{a}_3`$를 $`C`$에 넣는다.
4. $`C`$가 $`r`$개의 열을 가질 때까지 계속한다.

**수 $`r`$이 $`A`$와 $`C`$의 랭크이다.**

$`\Rightarrow C\mathbf{x} = \mathbf{0}`$이면 $`\mathbf{x} = \mathbf{0}`$일 때만 성립. 열들의 어떤 결합도 영벡터를 만들지 않는다.

**예시:**

```math
A = \begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = (\mathbf{a}_1 \ \mathbf{a}_2 \ \mathbf{a}_3)
```

1. $`C = \begin{pmatrix} 2 \\ 4 \\ 1 \end{pmatrix}`$
2. $`\mathbf{a}_2 = 3\mathbf{a}_1`$
3. $`\mathbf{a}_3 \neq c_1\mathbf{a}_1`$이므로, $`C = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}`$

$`\therefore`$ 랭크 $`r = 2`$.

**질문. $`A = CR`$에서 $`R`$은 무엇인가?**

```math
\begin{pmatrix} 2 & 6 & 4 \\ 4 & 12 & 8 \\ 1 & 3 & 5 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 3 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
3 \times 3 = (3 \times 2)(2 \times 3)
```

열을 재배열하면:

```math
\begin{pmatrix} 2 & 4 & 6 \\ 4 & 8 & 12 \\ 1 & 5 & 3 \end{pmatrix} = \begin{pmatrix} 2 & 4 \\ 4 & 8 \\ 1 & 5 \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 0 \end{pmatrix}
```

**관찰:**

1. $`C`$는 $`A`$의 $`r`$개의 독립 열의 완전한 집합을 포함한다.
2. $`R = (I \ F)`$는 항등행렬 $`I \in \mathbb{R}^{r \times r}`$을 포함한다.
3. $`A`$의 종속 열은 $`C`$에 있는 독립 열들의 결합이다.
4. $`A = CR = C(I \ F) = (C \ CF)`$
5. $`C`$는 $`A`$와 같은 열공간을 가진다. $`R`$은 $`A`$와 같은 행공간을 가진다.

**예시 (ex9):**

```math
\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 4 & 5 \\ 7 & 8 \end{pmatrix}\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}
```

랭크 $`r = 2`$ 행렬.

$`\mathbf{a}_j = C\mathbf{r}_j`$

$`(1 \ 2 \ 3) = (1 \ 2)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`(4 \ 5 \ 6) = (4 \ 5)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`(7 \ 8 \ 9) = (7 \ 8)\begin{pmatrix} 1 & 0 & -1 \\ 0 & 1 & 2 \end{pmatrix}`$

$`A`$의 행 $`i`$ = $`C`$의 행 $`i`$ 곱하기 $`R`$: $`R`$의 행들의 결합.

**예시:** $`A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 1 & 2 & 4 & 5 \end{pmatrix}_{2 \times 4}`$

```math
= \begin{pmatrix} 1 & 3 \\ 1 & 4 \end{pmatrix}\begin{pmatrix} 1 & 2 & 0 & 1 \\ 0 & 0 & 1 & 1 \end{pmatrix}
```

$`C \in \mathbb{R}^{2 \times 2}`$, $`R \in \mathbb{R}^{2 \times 4}`$. 랭크 $`r = 2`$이며, $`R`$을 사용하여 $`C`$로부터 $`A`$의 모든 열을 복원한다.

### R 행렬을 구하는 방법

제3장에서 다룰 **"소거법 (Elimination)"** 을 사용할 수 있다.

**예시:**

```math
A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}
```

**단계 1:** $`R2 - 2R1`$, $`R3 - 3R1`$:

```math
\begin{pmatrix} 1 & 3 & 4 \\ 0 & -2 & -6 \\ 0 & -2 & -6 \end{pmatrix}
```

**단계 2:** $`R3 - R2`$:

```math
\begin{pmatrix} 1 & 3 & 4 \\ 0 & -2 & -6 \\ 0 & 0 & 0 \end{pmatrix} \quad \text{행 랭크 } r = 2
```

**단계 3:** $`R2 / (-2)`$:

```math
\begin{pmatrix} 1 & 3 & 4 \\ 0 & 1 & 3 \\ 0 & 0 & 0 \end{pmatrix}
```

**단계 4:** $`R1 - 3R2`$:

```math
\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \\ 0 & 0 & 0 \end{pmatrix} \quad \leftarrow R \text{ 행렬}
```

$`\mathbf{a}_1`$과 $`\mathbf{a}_2`$의 열은 독립이다.

```math
A = \begin{pmatrix} 1 & 3 \\ 2 & 4 \\ 3 & 7 \end{pmatrix}\begin{pmatrix} 1 & 0 & -5 \\ 0 & 1 & 3 \end{pmatrix} = CR
```

### A = CR의 핵심 성질

- $`A = CR`$
- $`C`$의 $`r`$개의 열은 $`A`$의 열공간의 **기저** 이다.
- $`R`$의 $`r`$개의 행은 $`A`$의 행공간의 **기저** 이다.
- $`\Rightarrow`$ $`r`$ 차원.

> $`A = CR`$에서, "$`R`$"은 **기약행사다리꼴 (Reduced Row Echelon Form)** 이다.

**예시 (ex5):**

```math
A = \begin{pmatrix} | & | & | \\ \mathbf{a}_1 & \mathbf{a}_2 & 3\mathbf{a}_1 + 4\mathbf{a}_2 \\ | & | & | \end{pmatrix}_{m \times 3}
```

```math
= \begin{pmatrix} | & | \\ \mathbf{a}_1 & \mathbf{a}_2 \\ | & | \end{pmatrix}\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & 4 \end{pmatrix}
```

```math
= \mathbf{a}_1(1 \ 0 \ 3) + \mathbf{a}_2(0 \ 1 \ 4) = (\mathbf{a}_1 \ \mathbf{a}_2 \ 3\mathbf{a}_1 + 4\mathbf{a}_2) = A
```

### 1.4.7 A의 열 곱하기 B의 행 (외적, Outer Product)

```math
AB = \begin{pmatrix} | & | & & | \\ \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \\ | & | & & | \end{pmatrix}_{m \times n} \begin{pmatrix} - & \mathbf{b}_1^* & - \\ - & \mathbf{b}_2^* & - \\ & \vdots & \\ - & \mathbf{b}_n^* & - \end{pmatrix}_{n \times p}
```

```math
= \mathbf{a}_1\mathbf{b}_1^* + \mathbf{a}_2\mathbf{b}_2^* + \cdots + \mathbf{a}_n\mathbf{b}_n^*
```

여기서 $`\mathbf{a}_k\mathbf{b}_k^* = \begin{pmatrix} | \\ \mathbf{a}_k \\ | \end{pmatrix}_{m \times 1}(- \ \mathbf{b}_k^* \ -)_{1 \times p}`$

이것은 **열 곱하기 행 = 랭크 1 행렬** 이며, $`mp`$개의 원소를 가진다.

```math
AB = \sum_{k=1}^{n} \mathbf{a}_k\mathbf{b}_k^*
```

$`n`$개의 랭크 1 행렬의 합. 각각 $`mp`$번의 곱셈 $`\Rightarrow`$ 총 $`mnp`$번의 곱셈.

**참고:** 행렬 $`A = \begin{pmatrix} 1 & 3 & 4 \\ 2 & 4 & 2 \\ 3 & 7 & 6 \end{pmatrix}`$에 대해, $`\text{rank}(A) = 2`$.

$`\Rightarrow`$ $`A`$는 역행렬이 없다.
$`\Rightarrow`$ $`A`$의 행렬식은 0이다.

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:--------|:---------|
| **선형결합 (Linear Combination)** | $`c\mathbf{v} + d\mathbf{w}`$: 스칼라 곱셈 + 벡터 덧셈 |
| **열공간 (Column Space) $`C(A)`$** | 모든 벡터 $`A\mathbf{x}`$ = $`A`$의 열들의 모든 결합 |
| **생성 (Span)** | 벡터 집합의 모든 선형결합 |
| **독립벡터 (Independent Vectors)** | $`A\mathbf{x} = \mathbf{0}`$이 $`\mathbf{x} = \mathbf{0}`$일 때만 성립 |
| **종속벡터 (Dependent Vector)** | 다른 벡터들의 결합으로 쓸 수 있음 |
| **랭크 (Rank)** | 독립 열의 수 (= 독립 행의 수) |
| **기저 (Basis)** | 처음 $`r`$개의 독립 열이 $`C(A)`$의 기저를 형성 |
| **내적 (Dot Product)** | $`\mathbf{v} \cdot \mathbf{w} = \sum v_i w_i`$; 길이와 각도 정보 제공 |
| **길이 (노름, Norm)** | $`\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}}`$ |
| **단위벡터 (Unit Vector)** | $`\|\mathbf{v}\| = 1`$; 임의의 영이 아닌 벡터의 정규화: $`\mathbf{v}/\|\mathbf{v}\|`$ |
| **직교성 (Perpendicularity)** | $`\mathbf{v} \perp \mathbf{w} \Leftrightarrow \mathbf{v} \cdot \mathbf{w} = 0`$ |
| **각도 공식 (Angle Formula)** | $`\cos\theta = \frac{\mathbf{v} \cdot \mathbf{w}}{\|\mathbf{v}\|\|\mathbf{w}\|}`$ |
| **슈바르츠 부등식 (Schwarz Inequality)** | $`|\mathbf{v} \cdot \mathbf{w}| \leq \|\mathbf{v}\|\|\mathbf{w}\|`$ |
| **삼각 부등식 (Triangle Inequality)** | $`\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|`$ |
| **피타고라스 정리 (Pythagorean Theorem)** | $`\mathbf{v} \perp \mathbf{w}`$이면: $`\|\mathbf{v} + \mathbf{w}\|^2 = \|\mathbf{v}\|^2 + \|\mathbf{w}\|^2`$ |
| **행렬 곱셈 $`AB`$** | $`(AB)_{ij} = \text{row}_i(A) \cdot \text{col}_j(B)`$; $`AB`$의 열 $`j = A\mathbf{b}_j`$ |
| **외적 (Outer Product)** | $`AB = \sum \mathbf{a}_k \mathbf{b}_k^*`$ (랭크 1 행렬의 합) |
| **$`AB`$의 비용** | $`A \in \mathbb{R}^{m \times n}`$, $`B \in \mathbb{R}^{n \times p}`$에 대해 $`mnp`$번의 곱셈 |
| **비교환성 (Non-commutativity)** | 일반적으로 $`AB \neq BA`$ |
| **결합법칙 (Associativity)** | $`(AB)C = A(BC)`$ |
| **분배법칙 (Distributivity)** | $`A(B+C) = AB + AC`$ |
| **$`A = CR`$ 분해** | $`C`$: $`A`$의 독립 열; $`R`$: 기약행사다리꼴 |
| **랭크 1 행렬** | 모든 열이 한 직선 위; $`A = \mathbf{c}\mathbf{r}^T`$ |
| **행 랭크 = 열 랭크** | 모든 행렬에 대해 항상 성립 |
| **항등행렬 (Identity Matrix)** | $`\mathbf{v} = I\mathbf{v}`$; 표준 기저 벡터가 열 |
| **법선벡터에 의한 평면** | $`\mathbf{w} \cdot \mathbf{n} = 0`$은 $`\mathbf{n}`$에 수직인 평면을 정의 |
| **가역성 (Invertibility)** | $`A`$ 가역 $`\Leftrightarrow`$ 열이 독립 $`\Leftrightarrow`$ $`A\mathbf{x} = \mathbf{b}`$의 유일한 해 |
| **소거법 실패 (Elimination Failure)** | $`\mathbf{v} \| \mathbf{w}`$일 때: 해 없음 또는 무한히 많은 해 |

---
