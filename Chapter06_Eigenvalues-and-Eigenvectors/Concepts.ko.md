# 제6장 강의 — 고유값과 고유벡터

> **Last Updated:** 2026-03-31

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
  - [2.5 증명: 서로 다른 고유값의 고유벡터는 일차독립](#25-증명-서로-다른-고유값의-고유벡터는-일차독립)
  - [2.6 A의 거듭제곱 (마르코프 행렬 예제)](#26-a의-거듭제곱-마르코프-행렬-예제)
  - [2.7 닮은 행렬](#27-닮은-행렬)
  - [2.8 행렬의 거듭제곱과 피보나치 수](#28-행렬의-거듭제곱과-피보나치-수)
  - [2.9 대각화 불가능한 행렬과 중복도](#29-대각화-불가능한-행렬과-중복도)
- [3. 대칭 양정치 행렬 (6.3)](#3-대칭-양정치-행렬-63)
  - [3.1 대칭 행렬의 핵심 성질](#31-대칭-행렬의-핵심-성질)
  - [3.2 스펙트럼 정리](#32-스펙트럼-정리)
  - [3.3 증명: 대칭 행렬은 정규직교 고유기저를 갖는다](#33-증명-대칭-행렬은-정규직교-고유기저를-갖는다)
  - [3.4 양정치 행렬의 정의](#34-양정치-행렬의-정의)
  - [3.5 양정치 행렬의 성질](#35-양정치-행렬의-성질)
  - [3.6 양정치 행렬 판별법](#36-양정치-행렬-판별법)
  - [3.7 풀이 예제: 양정치와 양반정치](#37-풀이-예제-양정치와-양반정치)
  - [3.8 타원과 이차 형식](#38-타원과-이차-형식)
  - [3.9 양정치 행렬과 최솟값 문제](#39-양정치-행렬과-최솟값-문제)
  - [3.10 양반정치 행렬](#310-양반정치-행렬)
  - [3.11 합동 행렬](#311-합동-행렬)
  - [3.12 최적화와 머신 러닝](#312-최적화와-머신-러닝)
- [4. 선형 미분방정식 풀기 (6.5)](#4-선형-미분방정식-풀기-65)
  - [4.1 핵심 사실](#41-핵심-사실)
  - [4.2 스칼라 상미분방정식 복습](#42-스칼라-상미분방정식-복습)
  - [4.3 du/dt = Au의 풀이](#43-dudt--au의-풀이)
  - [4.4 일반적인 n x n 풀이 절차](#44-일반적인-n-x-n-풀이-절차)
  - [4.5 행렬의 지수 함수](#45-행렬의-지수-함수)
  - [4.6 이계 방정식](#46-이계-방정식)
  - [4.7 2 x 2 행렬의 안정성](#47-2-x-2-행렬의-안정성)
  - [4.8 풀이 예제](#48-풀이-예제)
- [요약](#요약)

---

<br>

## 1. 고유값 소개 (6.1)

### 1.1 고유값과 고유벡터의 정의

$A$가 $\mathbf{x}$에 작용할 때, 방향을 바꾸지 않고 벡터 $\mathbf{x}$를 $\lambda$만큼 늘이거나 줄이기만 한다.

$$A\mathbf{x} = \lambda \mathbf{x}$$

- $\mathbf{x}$는 **영이 아닌 벡터**이며, **고유벡터**(eigenvector)라 한다.
- $\lambda$는 $\mathbf{x}$에 대응하는 **고유값**(eigenvalue)이다.

**공식 표현:**

1. $A\mathbf{x} = \lambda\mathbf{x}$이면, $\mathbf{x} \neq \mathbf{0}$는 $A$의 고유벡터이고, 수 $\lambda$는 고유값이다.

2. $A^n \mathbf{x} = \lambda^n \mathbf{x}$ (모든 $n$에 대해); $(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}$; 그리고 $\lambda \neq 0$이면 $A^{-1}\mathbf{x} = \frac{1}{\lambda}\mathbf{x} = \lambda^{-1}\mathbf{x}$.

**증명** ($A^{-1}$에 대해):

$$A\mathbf{x} = \lambda\mathbf{x} \implies \mathbf{x} = A^{-1}(A\mathbf{x}) = A^{-1}(\lambda\mathbf{x}) = \lambda A^{-1}\mathbf{x} \quad \square$$

3. $(A - \lambda I)\mathbf{x} = \mathbf{0} \implies \det(A - \lambda I) = 0$. 이 방정식은 $n$개의 $\lambda$를 생성한다.

4. $A \in \mathbb{R}^{n \times n}$일 때:

$$\det(A) = \lambda_1 \lambda_2 \cdots \lambda_n; \quad \text{trace}(A) = \lambda_1 + \lambda_2 + \cdots + \lambda_n$$

5. **사영 행렬**(projection matrix) $P$는 $\lambda = 1$ 또는 $\lambda = 0$을 갖는다. 정사각 행렬 $P$가 $P^2 = P$를 만족하면 "사영 행렬"이라 한다.

### 1.2 A의 거듭제곱과 고유값

관계식에 $A$를 곱하면 어떻게 되는가?

$$A(A\mathbf{x}) = A(\lambda\mathbf{x}) \implies A^2\mathbf{x} = \lambda A\mathbf{x} = \lambda^2 \mathbf{x}$$

$\mathbf{x}$에 $A$를 계속 곱하면:

$$A^3\mathbf{x} = \lambda^3\mathbf{x}, \quad A^k\mathbf{x} = \lambda^k\mathbf{x}, \quad \ldots, \quad A^{100}\mathbf{x} = \lambda^{100}\mathbf{x}$$

$$\boxed{A^k \mathbf{x} = \lambda^k \mathbf{x}}$$

**$A^k$의 거동:**

- $i = 1, 2, \ldots, n$에 대해 $|\lambda_i| < 1$이면, $A^k$는 결국 **영**에 수렴한다.
- 어떤 $|\lambda_i| > 1$이면, $A^k$는 결국 **증가**한다.
- $\lambda = 1$이면, 시스템 상태는 시간에 따라 증가하거나 감쇠하지 않고 **일정**하게 유지된다.

$A\mathbf{x} = \mathbf{x}$이면, $\mathbf{x}$는 시스템의 **고정점**(fixed point) 또는 **평형점**(equilibrium point)이다. 즉, 시스템이 $\mathbf{x}$ 방향에서 정상 상태를 유지한다. 시스템은 **정상 상태**(steady state)에 도달한다: $A^k \mathbf{x} = \mathbf{x}$.

### 1.3 고유값의 성질

$A \in \mathbb{R}^{n \times n}$ 행렬은 $n$개의 고유값을 갖는다. **특성 다항식**(characteristic polynomial)을 풀어 고유값을 구할 수 있다:

$$\det(A - \lambda I) = 0$$

**행렬의 성질은 고유값에 큰 영향을 미친다:**

**(1)** 행렬 $A$의 **대각합**(trace)은 고유값의 **합**과 같다.

예: $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$, $\text{trace}(A) = 1 + 4 = 5$

$\det(A - \lambda I) = 0 \implies \lambda_1 = 0, \lambda_2 = 5$이고, $\lambda_1 + \lambda_2 = 5$.

**(2)** $A$의 **행렬식**(determinant)은 고유값의 **곱**과 같다.

예: $\det(A) = 1 \cdot 4 - 2 \cdot 2 = 0$이고, $\lambda_1 \cdot \lambda_2 = 0 \cdot 5 = 0$.

**(3)** **대칭 행렬** (즉, $A = A^T$)은 **실수 고유값**만 갖는다.

예: $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0$이므로, $\lambda = \pm 1$.

**(4)** **대칭 양정치**(SPD) 행렬은 **실수이고 양수인** 고유값을 갖는다. SPD 행렬은 이차 형식(quadratic form) $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A\mathbf{x}$가 **유일한 최솟값**을 가짐을 보장하므로, 최적화(optimization)에서 매우 중요하다. (6.3절 참고)

**(5)** **대각 행렬**(diagonal matrix)의 경우, 고유값은 **대각 원소**이다.

예: $A = \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} a - \lambda & 0 \\ 0 & b - \lambda \end{vmatrix} = (a - \lambda)(b - \lambda) = 0$이므로, $\lambda_1 = a, \lambda_2 = b$.

**삼각 행렬**(triangular matrix)의 경우에도 고유값은 **대각 원소**이다.

예: $A = \begin{pmatrix} a & b \\ 0 & c \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} a - \lambda & b \\ 0 & c - \lambda \end{vmatrix} = (a - \lambda)(c - \lambda) = 0$이므로, $\lambda_1 = a, \lambda_2 = c$.

**(6)** **반대칭 행렬**(skew-symmetric matrix, 즉 $A = -A^T$)은 **순허수**(purely imaginary) 고유값 또는 **영** 고유값을 갖는다.

예: $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ -1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0$이므로, $\lambda = \pm i$.

### 1.4 고유값을 구하는 방정식

$A\mathbf{x} = \lambda\mathbf{x}$에서:

$$A\mathbf{x} - \lambda\mathbf{x} = (A - \lambda I)\mathbf{x} = \mathbf{0}$$

고유벡터는 $A - \lambda I$의 **영공간**(nullspace)을 구성한다.

고유값 $\lambda$를 알면, $(A - \lambda I)\mathbf{x} = \mathbf{0}$을 풀어 고유벡터를 구한다.

$\mathbf{x}$가 영이 아닌 벡터이면, $A - \lambda I$는 **특이**(singular)하다.

$$\boxed{\det(A - \lambda I) = 0}$$

이 특성 다항식($\det(A - \lambda I) = 0$)은 $\lambda$만을 포함하며, $\lambda$는 $A - \lambda I$의 주대각선을 따라 나타난다. 행렬식은 $(-\lambda)^n$을 포함한다.

- 이 방정식은 $n$개의 해 $\lambda_1$부터 $\lambda_n$을 갖는다.
- $A_{n \times n}$은 $n$개의 고유값을 갖는다.

**증명 (행렬식 = 고유값의 곱):**

$$A\mathbf{x} = \lambda\mathbf{x} \implies \det(A - \lambda I) = 0$$

$\lambda$에 대한 특성 다항식:

$$\det(\lambda I - A) = 0 \iff \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0 = 0$$

$$\iff (\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n) = 0$$

$\lambda = 0$을 대입하면:

$$\det(0 \cdot I - A) = \det(-A) = (-1)^n \det(A)$$

$$(-\lambda_1)(-\lambda_2)\cdots(-\lambda_n) = (-1)^n \lambda_1\lambda_2\cdots\lambda_n$$

$$\therefore \det(A) = \lambda_1\lambda_2\cdots\lambda_n \quad \square$$

또한, $(\lambda - \lambda_1)(\lambda - \lambda_2)\cdots(\lambda - \lambda_n)$을 전개하면:

$$= \lambda^n - (\lambda_1 + \lambda_2 + \cdots + \lambda_n)\lambda^{n-1} + \cdots + \det(A) = 0$$

$\lambda^{n-1}$의 계수는 $\text{trace}(A) = a_{11} + a_{22} + \cdots + a_{nn}$과 같다 (증명 생략).

**행렬의 행렬식의 도함수:**

i) $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, $\det(A) = ad - bc$.

$a = a(x), b = b(x), c = c(x), d = d(x)$로 놓으면:

$$\frac{d}{dx}\det(A) = \frac{d}{dx}(a(x)d(x) - b(x)c(x)) = a'd + ad' - b'c - bc'$$

$$= \underbrace{a'd - b'c} + \underbrace{ad' - bc'}$$

$$\begin{vmatrix} a & b \\ c & d \end{vmatrix}' = \begin{vmatrix} a' & b' \\ c & d \end{vmatrix} + \begin{vmatrix} a & b \\ c' & d' \end{vmatrix}$$

행렬식을 미분할 때, 나머지 행은 그대로 두고 **한 행(또는 열)만 미분**한다.

$n \times n$ 행렬의 경우: 행을 $\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_n$으로 쓰면:

$$\begin{vmatrix} a_{11} & \cdots & a_{1n} \\ \vdots & & \vdots \\ a_{n1} & \cdots & a_{nn} \end{vmatrix}' = \sum_{i=1}^{n} \begin{vmatrix} \mathbf{r}_1 \\ \vdots \\ \mathbf{r}_i' \\ \vdots \\ \mathbf{r}_n \end{vmatrix}$$

### 1.5 행렬식과 대각합

**관찰 1:** 소거법(elimination)은 고유값을 **보존하지 않는다**.

$$A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 \\ 0 & 0 \end{pmatrix}$$

$\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 2 & 6 - \lambda \end{vmatrix} = \lambda^2 - 7\lambda = 0$이므로, $\lambda_1 = 7, \lambda_2 = 0$.

$\det(A_0 - \lambda I) = \begin{vmatrix} 1 - \lambda & 3 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 - \lambda = 0$이므로, $\lambda_1 = 1, \lambda_2 = 0$. (다르다!)

**관찰 2:** 곱 $\lambda_1 \cdot \lambda_2$과 합 $\lambda_1 + \lambda_2$는 행렬로부터 구할 수 있다.

$A = \begin{pmatrix} 1 & 3 \\ 2 & 6 \end{pmatrix}$에 대해:

- $\lambda_1 \lambda_2 = 7 \cdot 0 = 0 = \det(A) = 6 - 6 = 0$
- $\lambda_1 + \lambda_2 = 7 + 0 = 7 = \text{trace}(A) = 1 + 6 = 7$

$n$개 고유값의 곱 $\lambda_1 \lambda_2 \cdots \lambda_n$은 $A$의 **행렬식**과 같다.

합 $\lambda_1 + \lambda_2 + \cdots + \lambda_n$은 $n$개 대각 원소의 합 = $A$의 **대각합**(trace)과 같다.

### 1.6 풀이 예제

**예제 1:** 마르코프 행렬(Markov matrix)

$$A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix}$$

$$\det(A - \lambda I) = \begin{vmatrix} 0.8 - \lambda & 0.3 \\ 0.2 & 0.7 - \lambda \end{vmatrix} = 0.56 - 1.5\lambda + \lambda^2 - 0.06 = \lambda^2 - \frac{3}{2}\lambda + \frac{1}{2} = (\lambda - \frac{1}{2})(\lambda - 1) = 0$$

$$\therefore \lambda_1 = 1, \quad \lambda_2 = \frac{1}{2}$$

이는 $A - \lambda_1 I = A - I$와 $A - \lambda_2 I = A - \frac{1}{2}I$가 **역행렬이 존재하지 않음**을 의미한다.

고유벡터 $\mathbf{x}_1$, $\mathbf{x}_2$는 $(A - I)\mathbf{x}_1 = \mathbf{0}$과 $(A - \frac{1}{2}I)\mathbf{x}_2 = \mathbf{0}$을 만족한다.

즉, $\mathbf{x}_1 \in \mathcal{N}(A - I)$이고 $\mathbf{x}_2 \in \mathcal{N}(A - \frac{1}{2}I)$이다.

i) $(A - I)\mathbf{x}_1 = \begin{pmatrix} -0.2 & 0.3 \\ 0.2 & -0.3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$-2x_1 + 3x_2 = 0$. $x_1 = 3$으로 선택하면, $x_2 = 2$: $\mathbf{x}_1 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}$

ii) $(A - \frac{1}{2}I)\mathbf{x}_2 = \begin{pmatrix} 0.3 & 0.3 \\ 0.2 & 0.2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_1 + x_2 = 0$이므로, $x_1 = -x_2$. $x_1 = 1$로 선택하면, $x_2 = -1$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

$A$를 $\mathbf{x}_1$에 곱하면: $A\mathbf{x}_1 = \mathbf{x}_1$, $A^2\mathbf{x}_1 = \mathbf{x}_1$, ..., $A^{100}\mathbf{x}_1 = \mathbf{x}_1$.

$A$를 $\mathbf{x}_2$에 곱하면: $A\mathbf{x}_2 = \frac{1}{2}\mathbf{x}_2$, $A^2\mathbf{x}_2 = (\frac{1}{2})^2\mathbf{x}_2$, ..., $A^{100}\mathbf{x}_2 = (\frac{1}{2})^{100}\mathbf{x}_2$.

$\mathbf{x}_1$과 $\mathbf{x}_2$는 모두 자신의 방향을 유지한다.

$A$의 고유벡터 $\mathbf{x}$는 $A^n\mathbf{x} = \lambda^n\mathbf{x}$이므로 모든 $A^n$의 고유벡터이기도 하다.

고유벡터 $\mathbf{x}_1$과 $\mathbf{x}_2$는 $\mathbb{R}^2$를 생성(span)한다.

임의의 벡터 $\mathbf{x}$는 $\mathbf{x}_1$과 $\mathbf{x}_2$의 일차결합이다: $\mathbf{x} = c\mathbf{x}_1 + d\mathbf{x}_2 = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix}$.

$$A\mathbf{x} = cA\mathbf{x}_1 + dA\mathbf{x}_2 = c\mathbf{x}_1 + d\tfrac{1}{2}\mathbf{x}_2$$

$$A^2\mathbf{x} = c\mathbf{x}_1 + d(\tfrac{1}{2})^2\mathbf{x}_2$$

$$A^{100}\mathbf{x} = c\,\underbrace{\mathbf{x}_1}_{\text{정상 상태}} + d(\tfrac{1}{2})^{100}\underbrace{\mathbf{x}_2}_{\text{감쇠 모드}} \approx c\,\mathbf{x}_1$$

$\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$로 놓으면:

$$\begin{pmatrix} 1 \\ 0 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} c \\ d \end{pmatrix}$$

$$\begin{pmatrix} c \\ d \end{pmatrix} = \begin{pmatrix} 3 & 1 \\ 2 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 & 1 \\ 2 & -3 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 0.2 \\ 0.4 \end{pmatrix}$$

$$= 0.2\begin{pmatrix} 3 \\ 2 \end{pmatrix} + 0.4\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix} + \begin{pmatrix} 0.4 \\ -0.4 \end{pmatrix}$$

$$A^{100}\begin{pmatrix} 1 \\ 0 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}$$

마찬가지로 $\mathbf{x} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$에 대해: $A^{100}\begin{pmatrix} 0 \\ 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}$.

$$A^{100}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \approx \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_1)$$

$A$의 거듭제곱이 높을수록, 열들이 **정상 상태**에 더 가까워진다.

**예제 2:** 사영 행렬(projection matrix)

$$P = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$$

고유값: $\lambda = 1$과 $\lambda = 0$.

$$\det(P - \lambda I) = \begin{vmatrix} \frac{1}{2} - \lambda & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} - \lambda \end{vmatrix} = (\frac{1}{2})^2 \begin{vmatrix} 1 - 2\lambda & 1 \\ 1 & 1 - 2\lambda \end{vmatrix} = \frac{1}{4}[(1 - 2\lambda)^2 - 1] = \frac{1}{4}(4\lambda^2 - 4\lambda) = \lambda(\lambda - 1) = 0$$

$$\therefore \lambda_1 = 1, \; \lambda_2 = 0$$

i) $\lambda_1 = 1$: $(P - I)\mathbf{x}_1 = \begin{pmatrix} -\frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & -\frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = x_2$, $x_1 = 1$로 선택: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $\lambda_2 = 0$: $P\mathbf{x}_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 + x_2 = 0$, $x_1 = -x_2$, $x_1 = 1$로 선택: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) $P\mathbf{x}_1 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \mathbf{x}_1$

$P\mathbf{x}_2 = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 0 \\ 0 \end{pmatrix} = \mathbf{0}$

$\therefore \mathbf{x}_1 \in C(P)$: 열공간(column space)은 자기 자신으로 사영된다.

$\mathbf{x}_2 \in \mathcal{N}(P)$

iv) $\mathbf{w} = \begin{pmatrix} 1 \\ -1 \end{pmatrix} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$로 놓으면

$P\mathbf{w} = P\begin{pmatrix} 1 \\ -1 \end{pmatrix} + P\begin{pmatrix} 2 \\ 2 \end{pmatrix} = \mathbf{0} + \begin{pmatrix} 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}$

사영 행렬의 고유값은 오직 **0과 1**뿐이다.

**예제 3:** 교환 행렬(exchange matrix)

$$E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$$

고유값: 1과 $-1$.

$$\det(E - \lambda I) = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = (\lambda - 1)(\lambda + 1) = 0$$

$$\therefore \lambda_1 = 1, \; \lambda_2 = -1$$

i) $\lambda_1 = 1$: $(E - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $-x_1 + x_2 = 0$, $x_1 = x_2$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $\lambda_2 = -1$: $(E + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = -x_2$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) 고유벡터 $\mathbf{x}_1$과 $\mathbf{x}_2$는 $P$의 것과 동일하다.

$$E = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} - \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 2P - I$$

행렬이 $I$만큼 이동(shift)하면, 각 $\lambda$도 1만큼 이동한다.

$\det(2P - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & 1 - \lambda \end{vmatrix} = (1 - \lambda)^2 - 1 = \lambda^2 - 2\lambda = \lambda(\lambda - 2) = 0$

$\therefore \lambda_1 = 2, \lambda_2 = 0$ ($2P$에 대해), 즉 $E = 2P - I$에 대해: $\lambda_1 = 1, \lambda_2 = -1$.

**예제 4:** 특이 행렬(singular matrix)

$$A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$$

$\lambda$와 대응하는 $\mathbf{x}$를 구하라.

$$\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 2 \\ 2 & 4 - \lambda \end{vmatrix} = (1 - \lambda)(4 - \lambda) - 4 = \lambda^2 - 5\lambda = \lambda(\lambda - 5) = 0$$

$$\therefore \lambda_1 = 5 \text{ 이고 } \lambda_2 = 0$$

i) $\lambda_1 = 5$: $(A - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$2x_1 - x_2 = 0 \implies x_2 = 2x_1$. $x_1 = 1$로 놓으면 $x_2 = 2$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$

ii) $\lambda_2 = 0$: $(A - 0I)\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_1 + 2x_2 = 0 \implies x_1 = -2x_2$. $x_2 = 1$로 놓으면 $x_1 = -2$: $\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}$

**비고 1:** $\mathbf{x}_1$과 $\mathbf{x}_2$는 $(A - \lambda I)$의 영공간에 있다.

이 예제에서 $A$가 특이(singular)이므로 $\lambda_2 = 0$이 고유값이다. $A$가 가역이면 "0"은 고유값이 아니다: $A\mathbf{x} = \mathbf{0} \implies \mathbf{x} = \mathbf{0}$.

**비고 2:** $A \in \mathbb{R}^{2 \times 2}$에서, $A - \lambda I$가 특이하면 두 행 모두 벡터 $(a, b)$의 배수이다:

$$\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

고유벡터는 $(b, -a)$의 임의 배수이다:

$$\begin{pmatrix} a & b \\ ka & kb \end{pmatrix}\begin{pmatrix} b \\ -a \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

**비고 3:** 이 예제는 두 개의 서로 다른 고유값을 가지며, $\mathbb{R}^2$를 생성한다. $2 \times 2$ 행렬이 하나의 고유값만 가지면 $\mathbb{R}^2$를 생성할 수 없다.

예: $A = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$, $\det(A - \lambda I) = (3 - \lambda)^2 = 0$이므로, $\lambda = 3$. 기하적 중복도(geometric multiplicity) = 1, 대수적 중복도(algebraic multiplicity) = 2. 행렬 $A$는 **결함**(defective)이 있으며 **대각화 불가능**하다.

### 1.7 허수 고유값

**예제 5:** 90도 회전

$$R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$$

실수 고유벡터가 없다.

$$\det(R - \lambda I) = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0 \implies \lambda = \pm i$$

$\lambda_1 + \lambda_2 = 0 = \text{trace}(R)$, $\lambda_1\lambda_2 = 1 = \det(R)$.

$R\mathbf{x} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -y \\ x \end{pmatrix}$

회전 후, 어떤 실수 벡터 $R\mathbf{x}$도 $\mathbf{x}$와 같은 방향에 머물지 않으므로, $R\mathbf{x} = \lambda\mathbf{x}$를 만족하는 실수 $\lambda$는 존재하지 않는다.

i) $\lambda_1 = i$: $(R - iI)\mathbf{x} = \begin{pmatrix} -i & -1 \\ 1 & -i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$-ix_1 = x_2$, $x_1 = 1$로 선택: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ -i \end{pmatrix}$

ii) $\lambda_2 = -i$: $(R + iI)\mathbf{x} = \begin{pmatrix} i & -1 \\ 1 & i \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$ix_1 = x_2$, $x_1 = 1$로 선택: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ i \end{pmatrix}$

복소 고유벡터 $\mathbf{x}_1$과 $\mathbf{x}_2$는 복소 공간에서 회전되면서도 방향을 유지한다.

### 1.8 회전 행렬의 고유값

회전 행렬(rotation matrix)은 $\lambda = e^{i\theta}$와 $e^{-i\theta}$를 갖는다.

$$R\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}, \quad R\begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} -\sin\theta \\ \cos\theta \end{pmatrix}$$

$$\therefore R = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

$$\det(R - \lambda I) = \begin{vmatrix} \cos\theta - \lambda & -\sin\theta \\ \sin\theta & \cos\theta - \lambda \end{vmatrix} = (\cos\theta - \lambda)^2 + \sin^2\theta$$

$$= \lambda^2 - 2\cos\theta\,\lambda + \cos^2\theta + \sin^2\theta = \lambda^2 - 2\cos\theta\,\lambda + 1$$

$$\therefore \lambda = \cos\theta \pm \sqrt{\cos^2\theta - 1} = \cos\theta \pm \sqrt{-\sin^2\theta} = \cos\theta \pm i\sin\theta$$

$$\lambda_1 = \cos\theta + i\sin\theta = e^{i\theta}, \quad \lambda_2 = \cos\theta - i\sin\theta = e^{-i\theta}$$

**회전 행렬 $R$의 두 가지 성질:**

1. $R$은 **직교 행렬**(orthogonal matrix): $R^T R = I$, $|\lambda| = 1$.
2. $R$은 **반대칭 행렬**(skew-symmetric matrix): $R = -R^T$, $\lambda$는 순허수.

### 1.9 AB와 A+B의 고유값

$A\mathbf{x} = \lambda\mathbf{x}$, $B\mathbf{x} = \beta\mathbf{x}$를 고려하자.

즉, $\lambda$와 $\beta$는 각각 $A$와 $B$의 고유값이다.

$$AB\mathbf{x} = A(\beta\mathbf{x}) = \beta A\mathbf{x} = \beta\lambda\mathbf{x}$$

**이것은 일반적으로 참이 아니다.** 왜? 일반적으로, $\mathbf{x}$는 $A$와 $B$ **모두**의 고유벡터가 아니기 때문이다.

마찬가지로, $(A + B)\mathbf{x} \neq (\lambda + \beta)\mathbf{x}$.

**비고:** $A$와 $B$가 동일한 $n$개의 일차독립 고유벡터를 공유하는 것은 $AB = BA$일 때 그리고 그때에만 성립한다.

---

<br>

## 2. 행렬의 대각화 (6.2)

### 2.1 대각화의 핵심 사실

**(1)** $AX = X\Lambda$의 열은 $A\mathbf{x}_k = \lambda_k\mathbf{x}_k$이다. 고유값 행렬 $\Lambda$는 대각 행렬이다.

**(2)** $X$에 있는 $n$개의 일차독립 고유벡터가 $A$를 대각화한다:

$$A = X\Lambda X^{-1}$$

$$AA = X\Lambda X^{-1} X\Lambda X^{-1} = X\Lambda^2 X^{-1}$$

$$\boxed{A^k = X\Lambda^k X^{-1}}$$

**(3)** $\mathbf{u}_{k+1} = A\mathbf{u}_k$를 $\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0$로 풀 수 있다.

**(4)** 동일한 고유값이 없으면 $\implies$ 고유벡터 $X$가 가역 $\implies$ $A$를 대각화할 수 있다.

중복 고유값 $\implies$ $A$가 일차독립 고유벡터가 너무 적을 수 있다 $\implies$ $X^{-1}$이 실패한다.

**(5)** 모든 행렬 $C = B^{-1}AB$는 $A$와 동일한 고유값을 갖는다. 이러한 $C$들은 $A$에 **닮은**(similar) 행렬이다.

### 2.2 대각화 절차

$\mathbf{x}$가 고유벡터이면, $A\mathbf{x} = \lambda\mathbf{x}$이다. $A$를 $\mathbf{x}$에 적용하는 것은 단지 $\lambda$를 곱하는 것이다 — **매우 효율적**이다.

$A$가 대각화 가능하면, $A^{100}\mathbf{x} = X\Lambda^{100}X^{-1}\mathbf{x}$—역시 **매우 효율적**이다.

**대각화:** $A_{n \times n}$이 일차독립 고유벡터 $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$을 갖는다고 하자.

$X = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)$으로 놓으면:

$$AX = (A\mathbf{x}_1 \; A\mathbf{x}_2 \; \cdots \; A\mathbf{x}_n) = (\lambda_1\mathbf{x}_1 \; \lambda_2\mathbf{x}_2 \; \cdots \; \lambda_n\mathbf{x}_n) = (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix}$$

$$AX = X\Lambda$$

$$AXX^{-1} = X\Lambda X^{-1} \implies A = X\Lambda X^{-1}$$

$$X^{-1}AX = X^{-1}X\Lambda = \Lambda$$

$$\therefore \Lambda = X^{-1}AX$$

### 2.3 풀이 예제: 대각화

$$A = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}$$

$\det(A - \lambda I) = \begin{vmatrix} 2 - \lambda & 4 \\ 0 & 6 - \lambda \end{vmatrix} = (6 - \lambda)(2 - \lambda) = 0$이므로, $\lambda_1 = 6, \lambda_2 = 2$.

i) $(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$-x_1 + x_2 = 0 \implies \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $(A - \lambda_2 I)\mathbf{x}_2 = \begin{pmatrix} 0 & 4 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_2 = 0$, $x_1 = 1$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

iii) $A\mathbf{x}_1 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = 6\begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $A\mathbf{x}_2 = \begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = 2\begin{pmatrix} 1 \\ 0 \end{pmatrix}$

$$A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix}$$

$$\begin{pmatrix} 2 & 4 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 6 & \\ & 2 \end{pmatrix}$$

곱하면: $X^{-1} = -\begin{pmatrix} 0 & -1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & -1 \end{pmatrix}$

$X^{-1}AX = X^{-1}X\Lambda = \Lambda = \begin{pmatrix} 6 & \\ & 2 \end{pmatrix}$

$AXX^{-1} = X\Lambda X^{-1}$

**관찰:** $A^2 = (X\Lambda X^{-1})(X\Lambda X^{-1}) = X\Lambda^2 X^{-1}$

$\implies X^{-1}A^2 X = \Lambda^2$

$A^2$는 $X$에 있는 동일한 고유벡터를 가지며, $\Lambda^2$에서 제곱된 고유값 $36, 4$를 갖는다.

### 2.4 대각화에 대한 비고

**비고:** 행렬 $X$는 고유벡터가 **일차독립**(LI)이므로 역행렬이 존재한다.

**비고:** 고유값이 $n$개의 서로 다른 수라고 하자. 이는 $n$개의 고유벡터가 일차독립임을 의미한다. $X^{-1}$이 존재한다. 중복 고유값이 없는 행렬은 대각화할 수 있다.

**비고:** 고유벡터에 영이 아닌 임의의 상수를 곱할 수 있다: $\alpha A\mathbf{x} = \alpha\lambda\mathbf{x} = A(\alpha\mathbf{x})$.

**비고:** $\Lambda$에서의 고유값 순서는 $X$에서의 고유벡터 순서와 동일하다:

$$A(\mathbf{x}_1 \; \mathbf{x}_2) = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} \lambda_1 & \\ & \lambda_2 \end{pmatrix} \implies A(\mathbf{x}_2 \; \mathbf{x}_1) = (\mathbf{x}_2 \; \mathbf{x}_1)\begin{pmatrix} \lambda_2 & \\ & \lambda_1 \end{pmatrix}$$

**비고:** 일부 행렬은 고유벡터가 너무 적다. 이러한 행렬은 **대각화할 수 없다**.

예: $A = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}$

$\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & -1 \\ 1 & -1 - \lambda \end{vmatrix} = -(1 - \lambda)(1 + \lambda) + 1 = -(1 - \lambda^2) + 1 = \lambda^2 = 0$

$\therefore \lambda_1 = \lambda_2 = 0$ ($\lambda$의 중복).

$(A - 0I)\mathbf{x} = A\mathbf{x} = \mathbf{0}$: $\begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = x_2$, $x_1 = 1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

고유벡터가 **하나**뿐이다. $A$는 **대각화할 수 없다**.

**비고:** 가역성(invertibility)은 **영이 아닌** 고유값과 관련된다. $\lambda_i = 0$이면, $\det(A) = \lambda_1\lambda_2\cdots\lambda_n = 0$이므로 $A$는 **특이**(singular)하다.

### 2.5 증명: 서로 다른 고유값의 고유벡터는 일차독립

**명제:** $A$가 $n$개의 서로 다른 고유값을 갖는 $n \times n$ 행렬이면, 대응하는 고유벡터는 **일차독립**(LI)이다.

**증명:** 고유벡터 $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$이 **일차종속**(LD, linearly dependent)이라 가정하자.

$\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p$가 일차독립이고, $\mathbf{x}_{p+1}, \mathbf{x}_{p+2}, \ldots, \mathbf{x}_n$이 일차종속이라 하자.

즉, **모두 영이 아닌** 상수가 존재하여:

$$\mathbf{x}_{p+1} = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_p\mathbf{x}_p$$

일차결합에 $A$를 곱하면:

$$A\mathbf{x}_{p+1} = \lambda_{p+1}\mathbf{x}_{p+1} = c_1 A\mathbf{x}_1 + c_2 A\mathbf{x}_2 + \cdots + c_p A\mathbf{x}_p = c_1\lambda_1\mathbf{x}_1 + c_2\lambda_2\mathbf{x}_2 + \cdots + c_p\lambda_p\mathbf{x}_p \quad \text{--- (1)}$$

일차결합에 $\lambda_{p+1}$을 곱하면:

$$\lambda_{p+1}\mathbf{x}_{p+1} = c_1\lambda_{p+1}\mathbf{x}_1 + c_2\lambda_{p+1}\mathbf{x}_2 + \cdots + c_p\lambda_{p+1}\mathbf{x}_p \quad \text{--- (2)}$$

(1)에서 (2)를 빼면:

$$\mathbf{0} = c_1(\lambda_1 - \lambda_{p+1})\mathbf{x}_1 + \cdots + c_p(\lambda_p - \lambda_{p+1})\mathbf{x}_p$$

$\lambda$가 서로 다르고 $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p$가 일차독립이므로:

$$c_1 = c_2 = \cdots = c_p = 0$$

이는 $\mathbf{x}_{p+1} = \mathbf{0}$임을 의미한다. 이는 가정에 모순이다.

따라서 $\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$은 **일차독립**이다. $\square$

### 2.6 A의 거듭제곱 (마르코프 행렬 예제)

$$A = \begin{pmatrix} 0.8 & 0.3 \\ 0.2 & 0.7 \end{pmatrix} \quad \text{(마르코프 행렬)}$$

$\lambda_1 = 1, \lambda_2 = 0.5 \implies \Lambda = \begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}$

$\mathbf{x}_1 = \begin{pmatrix} 0.6 \\ 0.4 \end{pmatrix}, \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \implies X = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}$

$$A = X\Lambda X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.5 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}$$

$$A^2 = X\Lambda^2 X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0.25 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}$$

$$A^k = X\Lambda^k X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1^k & \\ & (0.5)^k \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix}$$

$k \to \infty$이면, $(0.5)^k \to 0$:

$$A^{\infty} = X\Lambda^{\infty}X^{-1} = \begin{pmatrix} 0.6 & 1 \\ 0.4 & -1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0.4 & -0.6 \end{pmatrix} = \begin{pmatrix} 0.6 & 0.6 \\ 0.4 & 0.4 \end{pmatrix}$$

**Q.** $A^k$는 언제 영행렬로 수렴하는가?

**A.** 모든 $|\lambda_i| < 1$일 때.

### 2.7 닮은 행렬

$\Lambda$가 고정되어 있다고 하자. 고유벡터 행렬 $X$를 바꾸면, **같은 고유값** $\Lambda$를 갖는 서로 다른 행렬을 얻는다.

$$A_1 = X_1 \Lambda X_1^{-1}, \quad A_2 = X_2 \Lambda X_2^{-1}, \quad \ldots$$

이들은 **닮은 행렬**(similar matrices)이다: $\Lambda, A_1, A_2, \ldots$

모든 $A_1, A_2, \ldots$는 $C$에 닮았다. 이들은 모두 $C$의 고유값을 공유한다.

이 개념을 대각화 불가능한 행렬로 확장한다. 상수 행렬 $C$와 가역 행렬 $B$를 선택한다. $A = BCB^{-1}$을 구성한다. $A$와 $C$는 **닮았다**.

$A$와 $C$는 동일한 $n$개의 고유값을 갖는다.

**명제:** $C\mathbf{x} = \lambda\mathbf{x}$이면, $BCB^{-1}$은 새로운 고유벡터 $B\mathbf{x}$에 대해 동일한 고유값 $\lambda$를 갖는다.

**증명:**

$$(BCB^{-1})(B\mathbf{x}) = BCI\mathbf{x} = BC\mathbf{x} = B\lambda\mathbf{x} = \lambda(B\mathbf{x}) \quad \square$$

예: $A_1 = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$와 $A_2 = \begin{pmatrix} \frac{1}{2} & \frac{1}{2} \\ \frac{1}{2} & \frac{1}{2} \end{pmatrix}$는 닮았다.

$\det(A_1 - \lambda I) = \lambda(1 - \lambda) = 0$이므로, $\lambda = 0, 1$.

$\det(A_2 - \lambda I) = \lambda^2 - \lambda + \frac{1}{4} - \frac{1}{4} = \lambda(\lambda - 1) = 0$이므로, $\lambda = 0, 1$.

### 2.8 행렬의 거듭제곱과 피보나치 수

피보나치 수는 이전 두 수의 합이다:

$$0, 1, 1, 2, 3, 5, 8, 13, \ldots$$

$$a_k + a_{k+1} = a_{k+2}, \quad k \geq 0$$

$a_{k+1} = a_{k+1}$을 도입하면:

$$\begin{cases} a_k + a_{k+1} = a_{k+2} \\ a_{k+1} = a_{k+1} \end{cases}$$

$\mathbf{u}_k = \begin{pmatrix} a_{k+1} \\ a_k \end{pmatrix}$로 놓으면:

$$\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}\mathbf{u}_k = \mathbf{u}_{k+1}$$

$$\therefore \boxed{\mathbf{u}_{k+1} = A\mathbf{u}_k}$$

$\mathbf{u}_0 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$에서:

$\mathbf{u}_1 = A\mathbf{u}_0 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $\mathbf{u}_2 = A^2\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$, $\mathbf{u}_3 = A^3\mathbf{u}_0 = \begin{pmatrix} 3 \\ 2 \end{pmatrix}$, ...

$\mathbf{u}_k = A^k\mathbf{u}_0$

$A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}$, $\det(A - \lambda I) = \begin{vmatrix} 1 - \lambda & 1 \\ 1 & -\lambda \end{vmatrix} = -\lambda(1 - \lambda) - 1 = \lambda^2 - \lambda - 1 = 0$

$$\therefore \lambda = \frac{1 \pm \sqrt{1 + 4}}{2} = \frac{1}{2} \pm \frac{\sqrt{5}}{2}$$

두 개의 서로 다른 고유값 $\implies$ 두 개의 일차독립 고유벡터 $\implies$ $X^{-1}$ 존재 $\implies$ $A = X\Lambda X^{-1}$.

$$\mathbf{u}_k = A^k \mathbf{u}_0 = X\Lambda^k X^{-1}\mathbf{u}_0$$

i) $\mathbf{u}_0$을 일차결합 $X\mathbf{c}$로 쓰면:

$$\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = X\mathbf{c} \implies \mathbf{c} = X^{-1}\mathbf{u}_0$$

ii) $\Lambda^k$를 $\mathbf{c}$에 곱하면:

$$\begin{pmatrix} \lambda_1^k & \\ & \lambda_2^k \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix}$$

iii) $X$를 $\Lambda^k\mathbf{c}$에 곱하면:

$$\mathbf{u}_k = (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1\lambda_1^k \\ c_2\lambda_2^k \end{pmatrix} = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2$$

**$A \in \mathbb{R}^{n \times n}$으로 일반화:**

$$\mathbf{u}_k = c_1\lambda_1^k\mathbf{x}_1 + c_2\lambda_2^k\mathbf{x}_2 + \cdots + c_n\lambda_n^k\mathbf{x}_n$$

$\mathbf{u}_{k+1} = A\mathbf{u}_k$의 해.

### 2.9 대각화 불가능한 행렬과 중복도

$\lambda$가 $A$의 고유값이라 하자.

1. **고유벡터 (기하적):** $A\mathbf{x} = \lambda\mathbf{x}$, 영이 아닌 $\mathbf{x}$.
2. **고유값 (대수적):** $\det(A - \lambda I) = 0$.

$\lambda$는 단순 고유값이거나, **중복** 고유값일 수 있다 (예: $\lambda^2 = 0 \to \lambda = 0$).

**중복도(Multiplicity):**

1. **기하적 중복도 (GM):** $\lambda$에 대한 독립 벡터의 수 = $\dim \mathcal{N}(A - \lambda I)$.
2. **대수적 중복도 (AM):** 고유값 중 $\lambda$의 반복 횟수, 즉 $\det(A - \lambda I) = 0$의 $n$개 근 중 $\lambda$의 수.

예: $A$가 $\lambda = 4, 4, 4$를 가지면: AM = 3, GM = 1, 2, 또는 3.

예: $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 0 & -\lambda \end{vmatrix} = \lambda^2 = 0$이므로, $\lambda = 0, 0$.

AM = 2, GM = 1 (고유벡터 1개).

$A\mathbf{x} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$이므로, $x_2 = 0$: $\mathbf{x} = \begin{pmatrix} c \\ 0 \end{pmatrix}$.

**GM < AM**이면, $A$는 **대각화 불가능**하다.

---

<br>

## 3. 대칭 양정치 행렬 (6.3)

### 3.1 대칭 행렬의 핵심 성질

**(1)** 대칭 행렬 $S$는 $n$개의 **실수** 고유값 $\lambda_i$와, $n$개의 **정규직교**(orthonormal) 고유벡터 $\mathbf{q}_i$를 갖는다.

**(2)** $S$는 직교 고유벡터 행렬 $Q$로 대각화된다:

$$S = Q\Lambda Q^{-1} = Q\Lambda Q^T$$

**(3a)** **양정치**(positive definite): 모든 $\lambda_i > 0$, 모든 피벗(pivot) $> 0$, 모든 왼쪽 상단 행렬식(leading determinant) $> 0$.

**(3b)** 에너지 판정법은 $\mathbf{x}^T S\mathbf{x} > 0$ (모든 $\mathbf{x} \neq \mathbf{0}$에 대해). 그러면 $S = A^T A$이고 $A$의 열이 독립이다.

**(4)** **양반정치**(positive semidefinite)는 $\lambda = 0$, 피벗 = 0, 행렬식 = 0, 에너지 $\mathbf{x}^T S\mathbf{x} = 0$, 임의의 $A$를 허용한다.

**대칭 행렬 ($S = S^T$)은 특별하다:**

1. 모든 $n$개의 고유값 $\lambda$는 **실수**이다.
2. $n$개의 고유벡터 $\mathbf{q}$는 **직교**하도록 선택할 수 있다.

예: $I = I^T = \begin{pmatrix} 1 & & \\ & 1 & \\ & & \ddots \end{pmatrix}$, $\lambda = 1, 1, \ldots, 1$. $I\mathbf{x} = 1 \cdot \mathbf{x}$—모든 영이 아닌 벡터 $\mathbf{x}$가 고유벡터이다.

직교하도록 선택할 수 있고, **단위 벡터**로 재조정(rescale)할 수 있다 $\implies$ **정규직교**(orthonormal).

### 3.2 스펙트럼 정리

$Q = (\mathbf{q}_1 \; \mathbf{q}_2 \; \cdots \; \mathbf{q}_n)$이고, $\|\mathbf{q}_i\| = 1$, $\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 1 & i = j \\ 0 & i \neq j \end{cases}$로 놓으면:

$$Q^T Q = \begin{pmatrix} - \mathbf{q}_1^T - \\ - \mathbf{q}_2^T - \\ \vdots \\ - \mathbf{q}_n^T - \end{pmatrix}\begin{pmatrix} | & | & & | \\ \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_n \\ | & | & & | \end{pmatrix} = I$$

정사각 행렬 $Q$에 대해: $Q^T Q = I \implies Q^T = Q^{-1}$.

$AX = X\Lambda$을 상기하자. 대칭 행렬 $S$에 대해, $X$ 대신 정규직교 $Q$를 사용한다:

$$SQ = Q\Lambda$$

$$\boxed{SQQ^T = Q\Lambda Q^{-1} = Q\Lambda Q^T}$$

**스펙트럼 정리(Spectral Theorem):** 모든 실수 대칭 행렬 $S$는 다음 형태를 갖는다

$$S = Q\Lambda Q^T$$

대칭성 증명: $S^T = (Q\Lambda Q^T)^T = Q\Lambda^T Q^T = Q\Lambda Q^T = S$.

$S$가 정규직교 고유벡터를 가지면, $S$는 대칭이다.

### 3.3 증명: 대칭 행렬은 정규직교 고유기저를 갖는다

**명제:** $S = S^T$이면, $S$는 정규직교 고유기저(orthonormal eigenbasis)를 갖는다.

**증명:** 두 고유벡터 $\mathbf{u}, \mathbf{v}$를 고려하자.

i) $\mathbf{v} \cdot (S\mathbf{u}) = \mathbf{v}^T S\mathbf{u} = \mathbf{v}^T S^T\mathbf{u}$ (대칭성에 의해) $= (S\mathbf{v})^T\mathbf{u} = (S\mathbf{v}) \cdot \mathbf{u}$

ii) $\alpha, \beta$를 대응하는 고유값이라 하자: $S\mathbf{u} = \alpha\mathbf{u}$, $S\mathbf{v} = \beta\mathbf{v}$.

$\mathbf{v} \cdot (S\mathbf{u}) = (S\mathbf{v}) \cdot \mathbf{u}$로부터:

$$\mathbf{v} \cdot (\alpha\mathbf{u}) = (\beta\mathbf{v}) \cdot \mathbf{u}$$

$$(\alpha - \beta)\mathbf{v} \cdot \mathbf{u} = 0$$

$\alpha \neq \beta$이면: $\mathbf{v} \cdot \mathbf{u} = 0$.

$\therefore \mathbf{v} \perp \mathbf{u}$ — **직교**(orthogonal).

$\mathbf{u}, \mathbf{v}$를 $\|\mathbf{u}\| = \|\mathbf{v}\| = 1$로 재조정하면, $\mathbf{u}, \mathbf{v}$는 **정규직교** 고유벡터이다. $\square$

### 3.4 양정치 행렬의 정의

**정의.** $n \times n$ 대칭 실수 행렬 $S$가 **양정치**(positive definite)라 함은

$$\mathbf{x}^T S\mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}, \; \mathbf{x} \in \mathbb{R}^n$$

(동치: $\forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$)

**정의.** $n \times n$ 대칭 실수 행렬 $S$가 **양반정치**(positive semidefinite, = 비음정치)라 함은

$$\mathbf{x}^T S\mathbf{x} \geq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n$$

**정의.** $n \times n$ 대칭 실수 행렬 $S$가 **음정치**(negative definite)라 함은

$$\mathbf{x}^T S\mathbf{x} < 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$$

**정의.** $n \times n$ 대칭 실수 행렬 $S$가 **음반정치**(negative semidefinite)라 함은

$$\mathbf{x}^T S\mathbf{x} \leq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n$$

### 3.5 양정치 행렬의 성질

**성질 1.** 양정치 행렬 $S$는 **모든 양의 고유값**을 갖는다.

**증명:** $S$는 대칭 $\implies S = Q\Lambda Q^T$, $S\mathbf{x} = Q\Lambda Q^T\mathbf{x}$.

$$\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\Lambda Q^T\mathbf{x} = (Q^T\mathbf{x})^T \Lambda (Q^T\mathbf{x})$$

$\mathbf{y} = Q^T\mathbf{x}$로 놓으면:

$$= \mathbf{y}^T \Lambda \mathbf{y} = \mathbf{y}^T\begin{pmatrix} \lambda_1 y_1 \\ \lambda_2 y_2 \\ \vdots \\ \lambda_n y_n \end{pmatrix} = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2$$

$S$가 양정치이므로, $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$.

$\mathbf{x} = Q\mathbf{e}_i$로 선택하면 ($\mathbf{e}_i$는 $I$의 $i$번째 열), $\mathbf{y} = Q^T\mathbf{x} = Q^T Q\mathbf{e}_i = \mathbf{e}_i$이므로, $\mathbf{x}^T S\mathbf{x} = \lambda_i > 0$.

따라서 모든 $\lambda_1, \lambda_2, \ldots, \lambda_n > 0$. $\square$

**비고:** 모든 고유값이 양수라 해서 행렬이 양정치인 것은 아니다.

예: $A = \begin{pmatrix} 1 & -3 \\ 0 & 1 \end{pmatrix} \to \lambda_1 = \lambda_2 = 1 > 0$이지만, $\mathbf{x}^T A\mathbf{x} = (x_1, x_2)\begin{pmatrix} x_1 - 3x_2 \\ x_2 \end{pmatrix} = x_1^2 - 3x_1 x_2 + x_2^2 < 0$ ($x_1 = x_2 = 1$일 때).

(동치 관계가 성립하려면 행렬이 **대칭**이어야 한다.)

**양의 에너지는 양의 고유값 ($\lambda > 0$)과 밀접한 관련이 있다:**

$S\mathbf{x} = \lambda\mathbf{x}$이면, $\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T\lambda\mathbf{x} = \lambda\mathbf{x}^T\mathbf{x} = \lambda\|\mathbf{x}\|^2$.

따라서 $\lambda > 0 \implies \mathbf{x}^T S\mathbf{x} > 0$ ($\mathbf{x} \neq \mathbf{0}$).

**명제:** $S$의 고유벡터에 대해 $\mathbf{x}^T S\mathbf{x} > 0$이면, 모든 $\mathbf{x} \neq \mathbf{0}$에 대해 $\mathbf{x}^T S\mathbf{x} > 0$이다.

**증명:** $\mathbf{x} = Q\mathbf{c} = c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n$으로 놓자. 여기서 $\mathbf{q}_i$는 $S$의 $i$번째 고유벡터이자 정규직교 벡터이다.

$$\mathbf{x}^T S\mathbf{x} = (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T S(c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)$$

$$= (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T(c_1\lambda_1\mathbf{q}_1 + c_2\lambda_2\mathbf{q}_2 + \cdots + c_n\lambda_n\mathbf{q}_n)$$

$$= c_1^2\lambda_1\mathbf{q}_1^T\mathbf{q}_1 + c_2^2\lambda_2\mathbf{q}_2^T\mathbf{q}_2 + \cdots + c_n^2\lambda_n\mathbf{q}_n^T\mathbf{q}_n > 0$$

$\lambda_1, \lambda_2, \ldots, \lambda_n > 0$이면 성립. $\square$

**성질 2.** $S$가 양정치이면, **가역**이고, $\det(S) > 0$이며, $S^{-1}$도 양정치이다.

**증명:**

i) $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$

$\implies S\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$

$\implies \mathcal{N}(S) = \{\mathbf{0}\}$

$\implies S$는 full rank $\iff \dim C(S) = n$

$\therefore S$는 **가역**이다.

ii) $S$가 양정치이므로, $S$는 모든 양의 고유값을 갖는다.

따라서 $0 < \lambda_1\lambda_2\cdots\lambda_n = \det(S)$.

iii) $(S^{-1})^T = (S^T)^{-1} = S^{-1}$이므로, $S^{-1}$은 **대칭**이다.

$\mathbf{x}^T S^{-1}\mathbf{x} = \mathbf{x}^T S^{-1}SS^{-1}\mathbf{x} = (S^{-T}\mathbf{x})^T S(S^{-1}\mathbf{x}) = (S^{-1}\mathbf{x})^T S(S^{-1}\mathbf{x})$

$\mathbf{z} = S^{-1}\mathbf{x}$로 놓으면: $= \mathbf{z}^T S\mathbf{z} > 0$. $\square$

**성질 3.** $S$가 양정치이면, $S = A^T A$이고 $A$의 열은 독립이다.

**증명:** 대칭 행렬 $S = Q\Lambda Q^T$. $S$가 양정치이므로 모든 고유값이 양수이다. 대각 행렬 $\Lambda$는 제곱근 $\sqrt{\Lambda}$를 갖는다:

$$\Lambda = \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix} = \begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix}\begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix} = \sqrt{\Lambda}\sqrt{\Lambda}$$

$A = Q\sqrt{\Lambda}Q^T$로 놓으면:

$$A^T A = (Q\sqrt{\Lambda}Q^T)(Q\sqrt{\Lambda}Q^T) = Q\sqrt{\Lambda}\underbrace{Q^T Q}_{I}\sqrt{\Lambda}Q^T = Q\sqrt{\Lambda}\sqrt{\Lambda}Q^T = Q\Lambda Q^T = S$$

**비고:** 에너지 $= \mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2$.

$\|A\mathbf{x}\| > 0$ ($A\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$인 경우) $\iff$ $A$가 full rank.

### 3.6 양정치 행렬 판별법

**(1) 양의 고유값:** 양정치 행렬 $S$는 모든 양의 고유값을 갖는다.

예: $S = \mathbf{u}\mathbf{v}^T$ (랭크 1 행렬)은 양정치가 아니다.

$A = \begin{pmatrix} a \\ b \end{pmatrix}(1 \; k) = \begin{pmatrix} a & ka \\ b & kb \end{pmatrix}$

랭크 1 행렬 $\to$ 1개 일차독립 벡터 $\to$ $(n-1)$개 종속 벡터.

$A\mathbf{u} = \mathbf{u}\mathbf{v}^T\mathbf{u} = (\mathbf{v}^T\mathbf{u})\mathbf{u}$이므로, $\lambda = \mathbf{v}^T\mathbf{u} = (1, k)\begin{pmatrix} a \\ b \end{pmatrix} = a + kb$.

$\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_{n-1}$를 $\mathcal{N}(A)$의 기저 벡터라 하자. $A\mathbf{w}_i = \mathbf{0} = 0 \cdot \mathbf{w}_i$이므로, $i = 1, 2, \ldots, n-1$에 대해 $\lambda = 0$.

$\therefore (n-1)$개의 영 고유값.

**(2) 양의 에너지 판정법:** $\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$.

**(3) $S = A^T A$에 대한 양의 에너지 판정법:**

$$\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2$$

예: $S = A^T A$ ($A = \begin{pmatrix} 1 & 1 & 2 \\ 1 & 2 & 3 \end{pmatrix}$, 열 $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$).

참고: $\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2$이므로, $\text{rank}(A) = 2 < 3$, $\dim C(A^T) = 2$, $\dim \mathcal{N}(A) = 1$.

$\dim \mathcal{N}(A) \neq 0$이므로, $\mathbf{x}^T S\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0$. $S$는 양정치가 아니다. $S$는 **양반정치**이다.

**(4) 행렬식 판정법:** $\det(S) > 0$인지 확인 — 더 정확하게는, 모든 **왼쪽 상단 행렬식**(leading determinant)을 확인한다.

$$S = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & -1 & 2 & -1 \\ & & -1 & 2 \end{pmatrix}$$

$D_1 = |2| = 2 > 0$

$D_2 = \begin{vmatrix} 2 & -1 \\ -1 & 2 \end{vmatrix} = 3 > 0$

$D_3 = \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{vmatrix} = 8 - 4 = 4 > 0$

$D_4 = |S| = 2D_3 + \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & -1 & -1 \end{vmatrix} = 8 - 4 + 1 = 5 > 0$

모든 왼쪽 상단 행렬식이 양수이므로, $S$는 양정치이다.

**(5) 피벗 판정법:** **피벗**이 양수인지 확인한다.

$S = \begin{pmatrix} 2 & -1 \\ -1 & 2 & -1 \\ & -1 & 2 \end{pmatrix}$의 경우:

소거 후: 피벗은 $2, \frac{3}{2}, \frac{4}{3}$ — 모두 양수.

$S = LU$, $S = LDL^T$, $S = A^T A$.

선행 행렬식(leading determinant)은 피벗과 밀접한 관련이 있다: $D_2/D_1$, $D_3/D_1$ 등.

**$2 \times 2$ SPD 행렬** $S = \begin{pmatrix} a & b \\ b & d \end{pmatrix}$에 대해:

- 행렬식 판정: $D_1 = a > 0$, $D_2 = ad - b^2 > 0$.
- 피벗 판정: $d_1 = a > 0$, $d_2 = \frac{ad - b^2}{a} > 0$.
- 고유값: $\lambda_1 > 0$, $\lambda_2 > 0$.
- 에너지: $(x \; y)\begin{pmatrix} a & b \\ b & d \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = ax^2 + bxy + byx + dy^2 = ax^2 + 2bxy + dy^2 > 0$.

### 3.7 풀이 예제: 양정치와 양반정치

**예제:** $S = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix} \implies \lambda_1 = 2 > 0, \lambda_2 = 6 > 0$

$S\mathbf{x} = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix}$

$\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix} = 2x_1^2 + 6x_2^2 > 0$. $\therefore S$는 양정치이다.

**예제:** $S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}$

$|S - \lambda I| = \begin{vmatrix} 5 - \lambda & 4 \\ 4 & 5 - \lambda \end{vmatrix} = \lambda^2 - 10\lambda + 25 - 16 = (\lambda - 9)(\lambda - 1) = 0$

$\therefore \lambda_1 = 9 > 0$, $\lambda_2 = 1 > 0$.

i) $(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 = x_2$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 + x_2 = 0$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) $Q = (\mathbf{x}_1 \; \mathbf{x}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$, $Q^{-1} = Q = Q^T$

iv) $\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T(Q\Lambda Q^T)\mathbf{x}$, $\mathbf{y} = Q^T\mathbf{x}$로 놓으면:

$= \lambda_1 y_1^2 + \lambda_2 y_2^2 = 9y_1^2 + y_2^2 > 0$

$\therefore S$는 양정치이다.

**대안적 접근법 (에너지 판정법):**

$\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = (x_1, x_2)\begin{pmatrix} 5x_1 + 4x_2 \\ 4x_1 + 5x_2 \end{pmatrix}$

$= 5x_1^2 + 4x_1x_2 + 4x_1x_2 + 5x_2^2 = x_1^2 + 4x_1^2 + 8x_1x_2 + 4x_2^2 + x_1^2$

$= x_1^2 + \underbrace{4(x_1^2 + 2x_1x_2 + x_2^2)}_{4(x_1 + x_2)^2} > 0$

**예제:** $S = \begin{pmatrix} 4 & 5 \\ 5 & 4 \end{pmatrix}$

$|S - \lambda I| = \lambda^2 - 8\lambda + 16 - 25 = \lambda^2 - 8\lambda - 9 = (\lambda - 9)(\lambda + 1) = 0$

$\therefore \lambda_1 = 9$이고 $\lambda_2 = -1$.

$\lambda_2 < 0$이므로, $S$는 양정치가 **아니다**.

### 3.8 타원과 이차 형식

**타원**(The Ellipse) $ax^2 + 2bxy + cy^2 = 1$.

**예제 1:** 타원 $5x^2 + 8xy + 5y^2 = 1$을 고려하자.

$$(x \; y)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 1 \implies \mathbf{x}^T S\mathbf{x} = 1$$

$S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}$, $\lambda^2 - 10\lambda + 9 = (\lambda - 9)(\lambda - 1) = 0$이므로, $\lambda = 9, 1$.

i) $(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) $\mathbf{q}_1 = \frac{1}{\sqrt{2}}\mathbf{x}_1$, $\mathbf{q}_2 = \frac{1}{\sqrt{2}}\mathbf{x}_2$

$Q = (\mathbf{q}_1 \; \mathbf{q}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$, $Q^{-1} = Q^T = Q$

$S = Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T$

iv) $\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T\mathbf{x}$

$\mathbf{z} = Q^T\mathbf{x}$로 놓으면: $= \mathbf{z}^T\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}\mathbf{z} = (z_1, z_2)\begin{pmatrix} 9z_1 \\ z_2 \end{pmatrix} = 9z_1^2 + z_2^2$

$\mathbf{z} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix} x_1 + x_2 \\ x_1 - x_2 \end{pmatrix}$

$\therefore z_1 = \frac{x_1 + x_2}{\sqrt{2}}, \; z_2 = \frac{x_1 - x_2}{\sqrt{2}}$

$\mathbf{x}^T S\mathbf{x} = 9z_1^2 + z_2^2 = 1$

$z_2 = 0$일 때: $z_1^2 = \frac{1}{9}$이므로, $z_1 = \pm \frac{1}{3}$ ($\mathbf{q}_1$ 방향의 반축).

$z_1 = 0$일 때: $z_2^2 = 1$이므로, $z_2 = \pm 1$ ($\mathbf{q}_2$ 방향의 반축).

**예제 2:** 양반정치

$$T = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}\begin{pmatrix} 3 & 1 \end{pmatrix} = A^T A$$

행렬식: $D_1 = 9 > 0$, $D_2 = 9 - 9 = 0 \to$ 역행렬 없음.

$\lambda^2 - \text{trace}(T)\lambda + \det(T) = \lambda^2 - 10\lambda = 0$이므로, $\lambda_1 = 10, \lambda_2 = 0$.

$(T - 10I)\mathbf{x}_1 = \begin{pmatrix} -1 & 3 \\ 3 & -9 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$

$(T - 0I)\mathbf{x}_2 = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -3 \end{pmatrix}$

에너지: $ax^2 + 2bxy + dy^2 = 9x^2 + 6xy + y^2 = (3x + y)^2 \geq 0$.

$E(x,y) = 1$은 띠(band)이다: $3x + y = \pm 1$. 에너지 $E(x,y)$의 그래프는 계곡(valley)이다. 띠의 축은 $T$의 고유벡터 방향을 따른다.

### 3.9 양정치 행렬과 최솟값 문제

에너지 $E = \mathbf{x}^T S\mathbf{x} = \begin{pmatrix} x \\ y \end{pmatrix}^T\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 5x^2 + 8xy + 5y^2 > 0$을 고려하자.

(그릇 모양, **볼록**(convex))

$E(x,y) = 0$일 때 $x = y = 0$이다. 이는 **최솟값 문제**와 연결된다.

이계 도함수의 행렬은 모든 점에서 양정치이다.

1차 도함수: $\frac{\partial E}{\partial x} = 10x + 8y$, $\frac{\partial E}{\partial y} = 8x + 10y$.

2차 도함수: $\frac{\partial^2 E}{\partial x^2} = 10 > 0$, $\frac{\partial^2 E}{\partial x \partial y} = 8 > 0$, $\frac{\partial^2 E}{\partial y^2} = 10 > 0$.

$$\nabla E = \left(\frac{\partial E}{\partial x}, \frac{\partial E}{\partial y}\right) \quad \xrightarrow{(x,y) = (0,0)} \quad \nabla E = (0, 0)$$

$$J(\nabla E^T) = \begin{pmatrix} \frac{\partial^2 E}{\partial x^2} & \frac{\partial^2 E}{\partial y \partial x} \\ \frac{\partial^2 E}{\partial x \partial y} & \frac{\partial^2 E}{\partial y^2} \end{pmatrix} = H \quad \text{(헤시안 행렬, Hessian matrix)}$$

$H_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}$

에너지를 $\frac{1}{2}\mathbf{x}^T S\mathbf{x}$로 정의하면, $H = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix} = S$.

예: $E = x^2 + y^2$는 $x = y = 0$에서 최솟값을 갖는다.

$f = e^{x^2 + y^2}$는 어떠한가?

$\nabla f = \left(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}\right)$이고, $\frac{\partial f}{\partial x} = e^{x^2+y^2} \cdot 2x$, $\frac{\partial f}{\partial y} = e^{x^2+y^2} \cdot 2y$.

$J(\nabla f)^T = \begin{pmatrix} (2 + 4x^2)e^{x^2+y^2} & 4xy \, e^{x^2+y^2} \\ 4xy \, e^{x^2+y^2} & (2 + 4y^2)e^{x^2+y^2} \end{pmatrix} = S$

$S$가 양정치이므로, $f$는 **순볼록**(strictly convex)이다.

**예제 1 (양정치):**

$S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix}$

행렬식: $D_1 = 9 > 0$, $D_2 = 27 - 9 = 18 > 0$.

피벗: $\begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} \xrightarrow{R_2 - \frac{1}{3}R_1} \begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix}$. 피벗: 9, 2 (둘 다 양수).

에너지: $E(x,y) = 9x^2 + 6xy + 3y^2 = (3x + y)^2 + 2y^2 > 0$. (순볼록 함수.)

대각합과 행렬식: $\text{trace}(S) = 12$, $\det(S) = 18$.

$\lambda^2 - 12\lambda + 18 = 0 \implies \lambda = 6 \pm \sqrt{36 - 18} = 6 \pm 3\sqrt{2}$

i) $\lambda_1 = 6 + 3\sqrt{2}$: $(S - \lambda_1 I)\mathbf{x} = \begin{pmatrix} 3 - 3\sqrt{2} & 3 \\ 3 & 3 - 3\sqrt{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$

$x_2 = (-1 + \sqrt{2})x_1$, $x_1 = 1$로 선택: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ -1 + \sqrt{2} \end{pmatrix}$

ii) $\lambda_2 = 6 - 3\sqrt{2}$: $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 - \sqrt{2} \end{pmatrix}$

**분해:**

$S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} = LU = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix} = LDL^T = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & \frac{1}{3} \\ 0 & 1 \end{pmatrix}$

$= (L\sqrt{D})(\sqrt{D}L^T) = A^T A$ ($A = \begin{pmatrix} 3 & 1 \\ 0 & \sqrt{2} \end{pmatrix}$)

### 3.10 양반정치 행렬

$\mathbf{x}^T S\mathbf{x} \geq 0$

$$S = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$$

$\det(S) = 4 - 4 = 0$ — 특이 행렬.

$\text{trace}(S) = 1 + 4 = 5$. $\lambda^2 - 5\lambda + 0 = \lambda(\lambda - 5) = 0$이므로, $\lambda_1 = 5, \lambda_2 = 0$.

i) $(S - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_2 = 2x_1$: $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$

ii) $S\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 + 2x_2 = 0$: $\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}$

$E_{21}S = U$: $S = E_{21}^{-1}U = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = LDL^T$

$= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$

$= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}$

$S$는 $A$의 **종속** 열을 갖는 $A^T A$로 분해된다.

$\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0$

### 3.11 합동 행렬

$S$가 양반정치이면, 모든 행렬 $A^T SA$도 양반정치이다:

$$\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} \geq 0 \quad \forall \mathbf{x}$$

여기서 $\mathbf{y} = A\mathbf{x}$. $\square$

$\mathbf{x}^T S\mathbf{x} > 0$이고, 모든 $\mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\}$에 대해 $A\mathbf{x} \neq \mathbf{0}$이면, $A^T SA$는 **양정치**이다.

$$\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} > 0 \quad \forall \mathbf{x} \in \mathbb{R}^n \setminus \{\mathbf{0}\} \quad \square$$

행렬 $A^T SA$를 $S$에 **"합동"**(congruent)이라 한다.

$S^T = S$가 $P$개의 양의 고유값, $N$개의 음의 고유값, $Z$개의 영 고유값을 가지면, $A$가 가역인 경우 $A^T SA$에 대해서도 동일하다.

**반대칭 행렬의 고유값 성질 증명:**

$A = -\bar{A}^T$. $A\mathbf{x} = \lambda\mathbf{x}$, $\overline{A\mathbf{x}} = \bar{\lambda}\bar{\mathbf{x}}$.

$\bar{\mathbf{x}}^T A\mathbf{x} = \bar{\mathbf{x}}^T(-\bar{A}^T)\mathbf{x} = -(\bar{A}\bar{\mathbf{x}})^T\mathbf{x}$

$= \bar{\mathbf{x}}^T\lambda\mathbf{x}$ 이고 $= -(\bar{\lambda}\bar{\mathbf{x}})^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}$

$\lambda\bar{\mathbf{x}}^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}$

$(\lambda + \bar{\lambda})\bar{\mathbf{x}}^T\mathbf{x} = 0$

$\therefore \lambda + \bar{\lambda} = (r + is) + (r - is) = 2r = 0$

$\lambda$의 실수 부분은 **영**이다. 반대칭 행렬은 영 또는 **순허수** 고유값을 갖는다.

### 3.12 최적화와 머신 러닝

$f(x)$를 최소화하기 위한 **경사 하강법**(gradient descent).

i) $f(x) = x^2 + 4x + 4 \implies f'(x) = 2x + 4$

ii) $x_0 = 10$, $\eta = 0.1$ (학습률, learning rate), 정지 기준 $|f'(x)| < 0.01$

iii) $x_{k+1} = x_k - \eta f'(x_k)$ — **최급강하 방향**(steepest direction).

$x_k \to x^*$ (최솟값 점)까지 근사를 반복한다.

$f(\mathbf{x})$에 대해: $\mathbf{x}_{k+1} = \mathbf{x}_k - \eta \nabla f(\mathbf{x}_k)$

$\mathbf{x}_k \to \mathbf{x}^*$이면, $\mathbf{x}^*$에서 $\nabla f$가 영 $\implies \mathbf{x}_{k+1} \approx \mathbf{x}_k$.

$J(\nabla f)^T = S$가 양정치이면, $f$는 **볼록**(convex) — $\mathbf{x}^*$를 찾기 쉽다.

---

<br>

## 4. 선형 미분방정식 풀기 (6.5)

### 4.1 핵심 사실

**(1)** $A\mathbf{x} = \lambda\mathbf{x}$이면, $\mathbf{u}(t) = e^{\lambda t}\mathbf{x}$는 $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$를 풀어준다. 각 $\lambda$와 $\mathbf{x}$는 해 $e^{\lambda t}\mathbf{x}$를 준다.

**(2)** $A = X\Lambda X^{-1}$이면,

$$\mathbf{u}(t) = e^{At}\mathbf{u}(0) = Xe^{\Lambda t}X^{-1}\mathbf{u}(0) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n$$

**(3)** 행렬 지수 함수: $e^{At} = I + At + \cdots + \frac{(At)^n}{n!} + \cdots = Xe^{\Lambda t}X^{-1}$ ($A = X\Lambda X^{-1}$인 경우).

**(4)** $A$의 모든 고유값의 **실수 부분이 $< 0$**이면, $A$는 **안정**(stable)이고 $\mathbf{u}(t) \to \mathbf{0}$이며 $e^{At} \to 0$.

**(5)** $\frac{d^2\mathbf{u}}{dt^2} + B\frac{d\mathbf{u}}{dt} + C\mathbf{u} = 0$은 $\frac{d}{dt}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix} = \begin{pmatrix} 0 & I \\ -C & -B \end{pmatrix}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix}$를 의미한다 ($\mathbf{v} = \frac{d\mathbf{u}}{dt}$).

### 4.2 스칼라 상미분방정식 복습

상미분방정식(ordinary differential equation)을 고려하자:

$$\frac{du}{dt} = u \implies u(t) = Ce^t$$

$$\boxed{\frac{du}{dt} = \lambda u} \implies \boxed{u(t) = Ce^{\lambda t}}$$

확인: $\frac{du}{dt} = \lambda Ce^{\lambda t} = \lambda u$.

$t = 0$일 때: $u(0) = Ce^{\lambda \cdot 0} = C \cdot 1 = C$. 초기값은 $C$이다.

$\therefore u(t) = u(0)e^{\lambda t} = u_0 e^{\lambda t}$

**거동:**

- $\lambda > 0$: **불안정**(unstable) (증가)
- $\lambda = 0$: **정상 상태**(steady state) (일정)
- $\lambda < 0$: **안정**(stable) (감쇠)

**$\lambda \in \mathbb{C}$이면?** 예: $\lambda = -1 + 2i$

$e^{\lambda t} = e^{-t + 2it} = e^{-t}e^{2it}$

$|e^{\lambda t}| = |e^{-t}||e^{2it}| = e^{-t} \cdot \underbrace{1}_{|e^{i\theta}| = 1}$

**관찰:** 안정성은 $\lambda$의 **실수 부분**에 의존한다.

- $\text{Re}(\lambda) > 0$: 불안정
- $\text{Re}(\lambda) = 0$: 정상 (진동)
- $\text{Re}(\lambda) < 0$: 안정

### 4.3 du/dt = Au의 풀이

**예제 1:** 연립 상미분방정식(coupled ODE)을 고려하자:

$$\frac{dy}{dt} = z, \quad \frac{dz}{dt} = y$$

$$\implies \frac{d}{dt}\begin{pmatrix} y \\ z \end{pmatrix} = \begin{pmatrix} z \\ y \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y \\ z \end{pmatrix}$$

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u}$$

초기 조건 $\mathbf{u}_0 = \begin{pmatrix} y(0) \\ z(0) \end{pmatrix}$가 주어질 때, $\mathbf{u}(t) = ?$

$\lambda$를 $A$의 고유값, $\mathbf{x}$를 대응하는 고유벡터라 하자.

$\mathbf{u} = e^{\lambda t}\mathbf{x}$로 선택하고 연립 ODE에 대입하면:

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u} = e^{\lambda t}A\mathbf{x} = e^{\lambda t}\lambda\mathbf{x} = \lambda e^{\lambda t}\mathbf{x}$$

$\mathbf{u} = e^{\lambda t}\mathbf{x}$가 연립 ODE를 만족하므로, $e^{\lambda t}\mathbf{x}$는 $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$의 해이다.

$$\implies \mathbf{u} = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2$$

$A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$의 고유값은 무엇인가?

$|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0$이므로, $\lambda = \pm 1$.

i) $\lambda_1 = 1$: $(A - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$이므로, $\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

ii) $\lambda_2 = -1$: $(A + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$이므로, $\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}$

iii) 완전한 해: $\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2$

iv) 상수 $c_1$, $c_2$는 초기 조건 $\mathbf{u}_0 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$로 결정할 수 있다:

$$\mathbf{u}(t=0) = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$$

$\mathbf{u}_0$을 $A$의 고유벡터 결합으로 쓰면:

$$= (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}$$

$$\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = -\frac{1}{2}\begin{pmatrix} -1 & -1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$$

### 4.4 일반적인 n x n 풀이 절차

$n \times n$ 행렬 $A$로 일반화:

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u}$$

$A$는 $n$개의 고유값과 $n$개의 고유벡터 ($\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n$)를 갖는다.

1. $\mathbf{u}_0$을 $\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_n\mathbf{x}_n$으로 쓴다.
2. 고유벡터 $\mathbf{x}_i$에 $e^{\lambda_i t}$를 곱한다.
3. $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$의 해는:

$$\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}$$

**예제 2:**

$$\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 1 \\ 0 & 0 & 3 \end{pmatrix}\mathbf{u}, \quad \mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}$$

$\lambda_1 = 1, \lambda_2 = 2, \lambda_3 = 3$ (상삼각 행렬 — 대각 원소가 고유값).

i) $(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$

ii) $(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}$

iii) $(A - 3I)\mathbf{x}_3 = \begin{pmatrix} -2 & 1 & 1 \\ 0 & -1 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_3 = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$

iv) $\mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2 \; \mathbf{x}_3)\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix}$

$$\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}^{-1}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}$$

$|X| = 1$, $X^{-1} = \frac{1}{|X|}C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}$

$$\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix}$$

v) $\therefore \mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + c_3 e^{\lambda_3 t}\mathbf{x}_3$

$$= 2e^t\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + 3e^{2t}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} + 4e^{3t}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$$

### 4.5 행렬의 지수 함수

상기: $e^x = 1 + x + \frac{1}{2}x^2 + \frac{1}{6}x^3 + \cdots + \frac{x^n}{n!} + \cdots = \sum_{i=0}^{\infty}\frac{x^i}{i!}$

$\frac{1}{1-x} = 1 + x + x^2 + x^3 + \cdots + x^n + \cdots = \sum_{i=0}^{\infty}x^i$ (등비급수, $|x| < 1$)

**행렬 지수 함수(matrix exponential):**

$$e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots + \frac{1}{n!}(At)^n + \cdots$$

$(I - At)^{-1} = I + At + (At)^2 + (At)^3 + \cdots + (At)^n + \cdots$

$\mathbf{u}(t) = e^{At}\mathbf{u}(0) = e^{At}\mathbf{u}_0$으로 놓자.

$\mathbf{u}$가 $\frac{d\mathbf{u}}{dt} = A\mathbf{u}$의 해인지 확인:

$$\frac{d\mathbf{u}}{dt} = \frac{d}{dt}(e^{At}\mathbf{u}_0) = \frac{de^{At}}{dt}\mathbf{u}_0 = \left(A + \frac{1}{2} \cdot 2 \cdot (At)A + \frac{3}{6}(At)^2 A + \cdots + \frac{n}{n!}(At)^{n-1}A + \cdots\right)\mathbf{u}_0$$

$$= A\left(I + At + \frac{1}{2}(At)^2 + \cdots + \frac{1}{(n-1)!}(At)^{n-1} + \cdots\right)\mathbf{u}_0 = Ae^{At}\mathbf{u}_0 = A\mathbf{u}$$

**$A$가 대각화 가능하면,** $A = X\Lambda X^{-1}$이고:

$$e^{At} = Xe^{\Lambda t}X^{-1}$$

여기서

$$e^{\Lambda t} = I + \Lambda t + \frac{1}{2}(\Lambda t)^2 + \frac{1}{6}(\Lambda t)^3 + \cdots = \begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}$$

**증명:**

$$e^{At} = I + At + \frac{1}{2}(At)^2 + \cdots$$

$$= I + X\Lambda X^{-1}t + \frac{1}{2}t^2(X\Lambda X^{-1})(X\Lambda X^{-1}) + \frac{1}{6}t^3(X\Lambda X^{-1})(X\Lambda X^{-1})(X\Lambda X^{-1}) + \cdots$$

$$= I + t(X\Lambda X^{-1}) + \frac{1}{2}t^2(X\Lambda^2 X^{-1}) + \frac{1}{6}t^3(X\Lambda^3 X^{-1}) + \cdots$$

$$= (XIX^{-1}) + (Xt\Lambda X^{-1}) + (X\frac{t^2}{2}\Lambda^2 X^{-1}) + (X\frac{t^3}{6}\Lambda^3 X^{-1}) + \cdots$$

$$= X\left(I + t\Lambda + \frac{t^2}{2}\Lambda^2 + \frac{t^3}{6}\Lambda^3 + \cdots\right)X^{-1}$$

$$\therefore e^{At} = Xe^{\Lambda t}X^{-1}$$

**관찰:** $A = X\Lambda X^{-1}$이고 $e^{At} = Xe^{\Lambda t}X^{-1}$ $\implies$ $e^{At}$와 $A$는 **동일한 고유벡터**를 공유한다.

**해:**

$$\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\mathbf{u}_0$$

$\mathbf{u}_0 = X\mathbf{c}$ (고유벡터의 일차결합)이므로:

$$= Xe^{\Lambda t}\underbrace{X^{-1}X}_{I}\mathbf{c} = Xe^{\Lambda t}\mathbf{c} = X\begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{pmatrix} = X\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}$$

$$= (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}$$

$$\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}$$

**성질:** $e^{A(s+t)} = e^{As} \cdot e^{At}$

**증명:** (이중 급수와 이항 정리를 사용한 상세 계산)

$$e^{As} \cdot e^{At} = \left(\sum_{j=0}^{\infty}\frac{(As)^j}{j!}\right)\left(\sum_{k=0}^{\infty}\frac{(At)^k}{k!}\right) = \sum_{j=0}^{\infty}\sum_{k=0}^{\infty}\frac{A^{j+k}s^j t^k}{j!\,k!}$$

$n = j + k$, $k = n - j$로 놓으면:

$$= \sum_{n=0}^{\infty}\frac{A^n}{n!}\sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^j t^{n-j} = \sum_{n=0}^{\infty}\frac{A^n}{n!}(s+t)^n = \sum_{n=0}^{\infty}\frac{(A(s+t))^n}{n!} = e^{A(s+t)} \quad \square$$

이항 정리(binomial theorem)를 사용: $(s + t)^n = \sum_{j=0}^{n}\binom{n}{j}s^{n-j}t^j = \sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^{n-j}t^j$.

$e^{A(s+t)} = e^{As} \cdot e^{At}$에서, $s = -t$로 놓으면:

$$e^{A(-t+t)} = e^0 = I = e^{-At} \cdot e^{At}$$

이는 $e^{-At}$가 $e^{At}$의 **역행렬**임을 의미한다.

**예제 4:** $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$, $e^{At} = ?$

$|A - \lambda I| = \lambda^2 + 1 = 0$이므로, $\lambda = \pm i$ (반대칭 행렬).

$A^2 = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}$

$A^3 = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$

$A^4 = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$

$$e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots$$

$$= \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}t + \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\frac{t^2}{2} + \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\frac{t^3}{6} + \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\frac{t^4}{24} + \cdots$$

$$= \begin{pmatrix} 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots & t - \frac{t^3}{6} + \cdots \\ -t + \frac{t^3}{6} - \cdots & 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots \end{pmatrix}$$

$$e^{At} = \begin{pmatrix} \cos(t) & \sin(t) \\ -\sin(t) & \cos(t) \end{pmatrix} \quad \text{(반대칭 행렬은 회전을 준다)}$$

**예제 5:** $\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}\mathbf{u}$를 $\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$ ($t = 0$)에서 풀어라.

상삼각 행렬: $\lambda_1 = 1, \lambda_2 = 2$.

i) $(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

ii) $(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

iii) $X = (\mathbf{x}_1 \; \mathbf{x}_2) = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \to X^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$

$A = X\Lambda X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$

$e^{At} = Xe^{\Lambda t}X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} e^t & \\ & e^{2t} \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$

iv) $\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\begin{pmatrix} 2 \\ 1 \end{pmatrix}$

$= Xe^{\Lambda t}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 2 \\ 1 \end{pmatrix} = Xe^{\Lambda t}\begin{pmatrix} 1 \\ 1 \end{pmatrix}$

$= e^t\begin{pmatrix} 1 \\ 0 \end{pmatrix} + e^{2t}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} e^t + e^{2t} \\ e^{2t} \end{pmatrix}$

### 4.6 이계 방정식

$my'' + by' + ky = 0$ (질량-스프링-댐퍼 시스템, mass-spring-damper system)을 고려하자.

$y = e^{\lambda t}$로 놓으면: $y' = \lambda e^{\lambda t} = \lambda y$, $y'' = (y')' = \lambda^2 y$.

스프링 방정식은 다음이 된다: $m\lambda^2 y + b\lambda y + ky = (m\lambda^2 + b\lambda + k)e^{\lambda t} = 0$

$$\implies m\lambda^2 + b\lambda + k = 0$$

$$\iff \lambda^2 + \frac{b}{m}\lambda + \frac{k}{m} = 0 \quad \text{---(*)}$$

$\lambda_1 + \lambda_2 = -\frac{b}{m}$, $\lambda_1\lambda_2 = \frac{k}{m}$

$\lambda = \frac{1}{2}\left(-\frac{b}{m} \pm \sqrt{\frac{b^2}{m^2} - \frac{4k}{m}}\right)$

$y_1 = e^{\lambda_1 t}$, $y_2 = e^{\lambda_2 t}$로 놓으면, 완전한 해는:

$$y(t) = c_1 y_1 + c_2 y_2 \quad \text{(} \lambda_1 \neq \lambda_2 \text{인 경우)} \quad \text{---(**)}$$

**이계 ODE를 일계 ODE로 변환:**

$my'' + by' + ky = 0 \iff y'' + \tilde{b}y' + \tilde{k}y = 0$ ($\tilde{b} = \frac{b}{m}$, $\tilde{k} = \frac{k}{m}$).

$\mathbf{y} = \begin{pmatrix} y' \\ y \end{pmatrix}$로 놓으면:

$y'' = -\tilde{b}y' - \tilde{k}y = (-\tilde{b} \; -\tilde{k})\begin{pmatrix} y' \\ y \end{pmatrix}$

$y' = y' = (1 \; 0)\begin{pmatrix} y' \\ y \end{pmatrix}$

$$\implies \frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -\tilde{b} & -\tilde{k} \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix} \implies \frac{d\mathbf{y}}{dt} = A\mathbf{y}$$

i) $A$의 고유값:

$|A - \lambda I| = \begin{vmatrix} -\tilde{b} - \lambda & -\tilde{k} \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + \tilde{b}\lambda + \tilde{k} = 0$

$\lambda$에 대한 방정식은 (*)와 동일하다.

ii) $\lambda_1, \lambda_2$를 $A$의 고유값이라 하면:

$(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -\tilde{b} - \lambda_1 & -\tilde{k} \\ 1 & -\lambda_1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$, $x_1 - \lambda_1 x_2 = 0$: $\mathbf{x}_1 = \begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix}$

$(A - \lambda_2 I)\mathbf{x}_2$: $x_1 - \lambda_2 x_2 = 0$: $\mathbf{x}_2 = \begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}$

iii) $\mathbf{y}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2$

$$\begin{pmatrix} y'(t) \\ y(t) \end{pmatrix} = c_1 e^{\lambda_1 t}\begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix} + c_2 e^{\lambda_2 t}\begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}$$

$$\therefore y(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t} \quad \text{(**)와 동일}$$

$m = 1, b = 0, k = 1$일 때: $y'' + y = 0 \iff y'' = -y$

$$\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -y \\ y' \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}$$

$|A - \lambda I| = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0$이므로, $\lambda = \pm i$ — **진동**(oscillation).

### 4.7 2 x 2 행렬의 안정성

$\frac{d\mathbf{u}}{dt} = A\mathbf{u}$의 해에 대해 근본적인 질문이 있다:

$t \to \infty$일 때 해가 $\mathbf{u} = \mathbf{0}$에 접근하는가?

완전한 해 $\mathbf{u}(t)$는 순수 해 $e^{\lambda t}\mathbf{x}$로 구성되므로, 안정성은 $A$의 고유값에 의존한다.

$\lambda = r + is$

$e^{\lambda t} = e^{rt}e^{ist}$

$|e^{\lambda t}| = |e^{rt}|$

$\lambda$의 **실수 부분**이 증가 ($r > 0$) 또는 감쇠 ($r < 0$)를 제어한다.

임의의 $2 \times 2$ 행렬 $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$에 대해:

$\lambda_1$과 $\lambda_2$가 음수이려면:

i) $\lambda_1 + \lambda_2 = a + d < 0$ (대각합이 음수)

ii) $\lambda_1 \lambda_2 > 0$ (행렬식이 양수)

### 4.8 풀이 예제

**예제 6:** $y'' + 4y' + 3y = 0$을 풀어라.

$y$를 $e^{\lambda t}$로 대체: $(\lambda^2 + 4\lambda + 3)e^{\lambda t} = 0$

$\implies \lambda^2 + 4\lambda + 3 = 0 \implies (\lambda + 3)(\lambda + 1) = 0$

$\therefore \lambda = -1, -3$

$\implies y(t) = c_1 e^{-t} + c_2 e^{-3t}$ — 감쇠하는 해 $\to$ **안정** 해.

$\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}$를 도입하면: $y'' = -4y' - 3y$, $y' = y'$

$$\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -4 & -3 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad \frac{d\mathbf{u}}{dt} = A\mathbf{u}$$

$|A - \lambda I| = \begin{vmatrix} -4 - \lambda & -3 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 4\lambda + 3 = 0$이므로, $\lambda = -1, -3$.

대응하는 고유벡터:

$\begin{pmatrix} -4+1 & -3 \\ 1 & +1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}$

$\begin{pmatrix} -4+3 & -3 \\ 1 & +3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} -3 \\ 1 \end{pmatrix}$

완전한 해는:

$$\begin{pmatrix} y' \\ y \end{pmatrix} = \mathbf{u}(t) = c_1 e^{-t}\begin{pmatrix} -1 \\ 1 \end{pmatrix} + c_2 e^{-3t}\begin{pmatrix} -3 \\ 1 \end{pmatrix}$$

이로부터 $y(t) = c_1 e^{-t} + c_2 e^{-3t}$.

**예제 7:** $y'' - 2y' + y = 0$

$y$를 $e^{\lambda t}$로 대체: $(\lambda^2 - 2\lambda + 1) = (\lambda - 1)^2 = 0$이므로, $\lambda = 1, 1$ (중근).

$\implies y(t) = c_1 e^t + ?$

근이 중복이므로, 두 번째 해는 $c_2 \cdot t e^t$이다.

$$y(t) = c_1 e^t + c_2 t e^t$$

$\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}$를 도입하면: $y'' = 2y' - y$, $y' = y'$

$$\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad A = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}$$

$|A - \lambda I| = \begin{vmatrix} 2 - \lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 2\lambda + 1 = 0$이므로, $\lambda = 1, 1$.

$(A - I)\mathbf{x}_1 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$

고유벡터의 수가 2보다 적으므로, 대각화는 **불가능**하다.

따라서, 급수의 정의로부터 $e^{At}$를 계산한다:

$$e^{At} = I + (At) + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots$$

$A = I + (A - I)$로 쓰면:

$$e^{At} = e^{It + (A-I)t} = e^{It} \cdot e^{(A-I)t}$$

$(A - I)^2 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$이므로, $k \geq 2$에 대해 $(A - I)^k = 0$.

$$e^{(A-I)t} = I + (A-I)t = I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t$$

$$e^{At} = e^t\left(I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t\right) = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}$$

$$\mathbf{u}(t) = e^{At}\mathbf{u}_0$$

$$\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}\begin{pmatrix} y_0' \\ y_0 \end{pmatrix}$$

$$\therefore y(t) = (e^t - te^t)y_0 + te^t y_0'$$

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:--------|:---------|
| 고유값 방정식 | $A\mathbf{x} = \lambda\mathbf{x}$, $\mathbf{x} \neq \mathbf{0}$ |
| 특성 다항식 | $\det(A - \lambda I) = 0$은 $n$개의 고유값을 준다 |
| 대각합과 행렬식 | $\text{trace}(A) = \sum \lambda_i$, $\det(A) = \prod \lambda_i$ |
| 고유값의 거듭제곱 | $A^k\mathbf{x} = \lambda^k\mathbf{x}$ |
| 역행렬의 고유값 | $A^{-1}\mathbf{x} = \lambda^{-1}\mathbf{x}$ ($\lambda \neq 0$인 경우) |
| 이동 규칙 | $(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}$ |
| 대칭 행렬 | 실수 고유값, 직교 고유벡터 |
| 반대칭 행렬 | 순허수 또는 영 고유값 |
| 회전 행렬 | $\lambda = e^{\pm i\theta}$ |
| 사영 행렬 | $\lambda = 0$ 또는 $1$만 가능 |
| 대각화 | $n$개 일차독립 고유벡터가 존재하면 $A = X\Lambda X^{-1}$ |
| 대각화를 통한 행렬 거듭제곱 | $A^k = X\Lambda^k X^{-1}$ |
| 서로 다른 고유값 | 고유벡터가 일차독립 $\implies$ 대각화 가능 |
| 닮은 행렬 | $C = B^{-1}AB$는 $A$와 같은 고유값을 가짐 |
| GM vs AM | GM $\leq$ AM; GM $<$ AM이면 대각화 불가능 |
| 스펙트럼 정리 | 대칭 $S$에 대해 $S = Q\Lambda Q^T$ |
| 양정치 | 모든 $\mathbf{x} \neq \mathbf{0}$에 대해 $\mathbf{x}^T S\mathbf{x} > 0$; 모든 $\lambda_i > 0$; 모든 피벗 $> 0$; 모든 선행 행렬식 $> 0$ |
| 양반정치 | $\mathbf{x}^T S\mathbf{x} \geq 0$; $\lambda = 0$ 허용 |
| $A^T A$의 에너지 판정 | $\mathbf{x}^T A^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0$ |
| 양정치 $\implies$ 가역 | $\det(S) > 0$, $S^{-1}$도 양정치 |
| 양정치 $\implies S = A^T A$ | $A$의 열이 독립 |
| 합동 행렬 | $A^T SA$는 정치성(definiteness) 유형을 보존 |
| 타원 방정식 | $\mathbf{x}^T S\mathbf{x} = 1$; 축은 고유벡터 방향, 길이 $1/\sqrt{\lambda}$ |
| 헤시안과 볼록성 | $H$ 양정치 $\implies$ $f$는 순볼록 |
| 경사 하강법 | $\mathbf{x}_{k+1} = \mathbf{x}_k - \eta\nabla f(\mathbf{x}_k)$ |
| ODE 해 | $\frac{d\mathbf{u}}{dt} = A\mathbf{u} \implies \mathbf{u}(t) = \sum c_i e^{\lambda_i t}\mathbf{x}_i$ |
| 행렬 지수 함수 | $e^{At} = I + At + \frac{(At)^2}{2!} + \cdots = Xe^{\Lambda t}X^{-1}$ |
| 안정성 (연속) | 모든 $\text{Re}(\lambda_i) < 0$이면 안정 |
| 안정성 (이산) | 모든 $|\lambda_i| < 1$이면 $A^k \to 0$ |
| $e^{A(s+t)} = e^{As}e^{At}$ | 행렬 지수 함수의 가법 성질 |
| $(e^{At})^{-1} = e^{-At}$ | 행렬 지수 함수의 역행렬 |
| 이계 ODE | $my'' + by' + ky = 0 \iff$ 일계 시스템 $\frac{d\mathbf{y}}{dt} = A\mathbf{y}$ |
| 중복 고유값 ODE | 대각화 불가 경우: $e^{At} = e^{It} \cdot e^{(A-I)t}$ 사용 |

---
