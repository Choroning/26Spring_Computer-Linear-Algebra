# 제6장 강의 — 고유값과 고유벡터 — Part 2/3

> **📑 이 문서는 3개 파트로 나뉘어 있습니다.** GitHub은 파일 하나당 렌더되는 수식 개수에 한도가 있어, 모든 수식이 제대로 보이도록 절 단위로 나눴습니다.
>
> [Part 1](Concepts.ko.md) · **Part 2** · [Part 3](Concepts-3.ko.md)

---

<br>

## 목차

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

---

<br>

## 3. 대칭 양정치 행렬 (6.3)

### 3.1 대칭 행렬의 핵심 성질

**(1)** 대칭 행렬 $`S`$는 $`n`$개의 **실수** 고유값 $`\lambda_i`$와, $`n`$개의 **정규직교**(orthonormal) 고유벡터 $`\mathbf{q}_i`$를 갖는다.

**(2)** $`S`$는 직교 고유벡터 행렬 $`Q`$로 대각화된다:

```math
S = Q\Lambda Q^{-1} = Q\Lambda Q^T
```

**(3a)** **양정치**(positive definite): 모든 $`\lambda_i > 0`$, 모든 피벗(pivot) $`> 0`$, 모든 왼쪽 상단 행렬식(leading determinant) $`> 0`$.

**(3b)** 에너지 판정법은 $`\mathbf{x}^T S\mathbf{x} > 0`$ (모든 $`\mathbf{x} \neq \mathbf{0}`$에 대해). 그러면 $`S = A^T A`$이고 $`A`$의 열이 독립이다.

**(4)** **양반정치**(positive semidefinite)는 $`\lambda = 0`$, 피벗 = 0, 행렬식 = 0, 에너지 $`\mathbf{x}^T S\mathbf{x} = 0`$, 임의의 $`A`$를 허용한다.

**대칭 행렬 ($`S = S^T`$)은 특별하다:**

1. 모든 $`n`$개의 고유값 $`\lambda`$는 **실수** 이다.
2. $`n`$개의 고유벡터 $`\mathbf{q}`$는 **직교** 하도록 선택할 수 있다.

예: $`I = I^T = \begin{pmatrix} 1 & & \\ & 1 & \\ & & \ddots \end{pmatrix}`$, $`\lambda = 1, 1, \ldots, 1`$. $`I\mathbf{x} = 1 \cdot \mathbf{x}`$—모든 영이 아닌 벡터 $`\mathbf{x}`$가 고유벡터이다.

직교하도록 선택할 수 있고, **단위 벡터** 로 재조정(rescale)할 수 있다 $`\implies`$ **정규직교**(orthonormal).

### 3.2 스펙트럼 정리

$`Q = (\mathbf{q}_1 \; \mathbf{q}_2 \; \cdots \; \mathbf{q}_n)`$이고, $`\|\mathbf{q}_i\| = 1`$, $`\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 1 & i = j \\ 0 & i \neq j \end{cases}`$로 놓으면:

```math
Q^T Q = \begin{pmatrix} - \mathbf{q}_1^T - \\ - \mathbf{q}_2^T - \\ \vdots \\ - \mathbf{q}_n^T - \end{pmatrix}\begin{pmatrix} | & | & & | \\ \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_n \\ | & | & & | \end{pmatrix} = I
```

정사각 행렬 $`Q`$에 대해: $`Q^T Q = I \implies Q^T = Q^{-1}`$.

$`AX = X\Lambda`$을 상기하자. 대칭 행렬 $`S`$에 대해, $`X`$ 대신 정규직교 $`Q`$를 사용한다:

```math
SQ = Q\Lambda
```

```math
\boxed{SQQ^T = Q\Lambda Q^{-1} = Q\Lambda Q^T}
```

**스펙트럼 정리(Spectral Theorem):** 모든 실수 대칭 행렬 $`S`$는 다음 형태를 갖는다

```math
S = Q\Lambda Q^T
```

대칭성 증명: $`S^T = (Q\Lambda Q^T)^T = Q\Lambda^T Q^T = Q\Lambda Q^T = S`$.

$`S`$가 정규직교 고유벡터를 가지면, $`S`$는 대칭이다.

### 3.3 증명: 대칭 행렬은 정규직교 고유기저를 갖는다

**명제:** $`S = S^T`$이면, $`S`$는 정규직교 고유기저(orthonormal eigenbasis)를 갖는다.

**증명:** 두 고유벡터 $`\mathbf{u}, \mathbf{v}`$를 고려하자.

i) $`\mathbf{v} \cdot (S\mathbf{u}) = \mathbf{v}^T S\mathbf{u} = \mathbf{v}^T S^T\mathbf{u}`$ (대칭성에 의해) $`= (S\mathbf{v})^T\mathbf{u} = (S\mathbf{v}) \cdot \mathbf{u}`$

ii) $`\alpha, \beta`$를 대응하는 고유값이라 하자: $`S\mathbf{u} = \alpha\mathbf{u}`$, $`S\mathbf{v} = \beta\mathbf{v}`$.

$`\mathbf{v} \cdot (S\mathbf{u}) = (S\mathbf{v}) \cdot \mathbf{u}`$로부터:

```math
\mathbf{v} \cdot (\alpha\mathbf{u}) = (\beta\mathbf{v}) \cdot \mathbf{u}
```

```math
(\alpha - \beta)\mathbf{v} \cdot \mathbf{u} = 0
```

$`\alpha \neq \beta`$이면: $`\mathbf{v} \cdot \mathbf{u} = 0`$.

$`\therefore \mathbf{v} \perp \mathbf{u}`$ — **직교**(orthogonal).

$`\mathbf{u}, \mathbf{v}`$를 $`\|\mathbf{u}\| = \|\mathbf{v}\| = 1`$로 재조정하면, $`\mathbf{u}, \mathbf{v}`$는 **정규직교** 고유벡터이다. $`\square`$

### 3.4 양정치 행렬의 정의

**정의.** $`n \times n`$ 대칭 실수 행렬 $`S`$가 **양정치**(positive definite)라 함은

```math
\mathbf{x}^T S\mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}, \; \mathbf{x} \in \mathbb{R}^n
```

(동치: $`\forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$)

**정의.** $`n \times n`$ 대칭 실수 행렬 $`S`$가 **양반정치**(positive semidefinite, = 비음정치)라 함은

```math
\mathbf{x}^T S\mathbf{x} \geq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n
```

**정의.** $`n \times n`$ 대칭 실수 행렬 $`S`$가 **음정치**(negative definite)라 함은

```math
\mathbf{x}^T S\mathbf{x} < 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace
```

**정의.** $`n \times n`$ 대칭 실수 행렬 $`S`$가 **음반정치**(negative semidefinite)라 함은

```math
\mathbf{x}^T S\mathbf{x} \leq 0 \quad \forall \; \mathbf{x} \in \mathbb{R}^n
```

### 3.5 양정치 행렬의 성질

**성질 1.** 양정치 행렬 $`S`$는 **모든 양의 고유값** 을 갖는다.

**증명:** $`S`$는 대칭 $`\implies S = Q\Lambda Q^T`$, $`S\mathbf{x} = Q\Lambda Q^T\mathbf{x}`$.

```math
\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\Lambda Q^T\mathbf{x} = (Q^T\mathbf{x})^T \Lambda (Q^T\mathbf{x})
```

$`\mathbf{y} = Q^T\mathbf{x}`$로 놓으면:

```math
= \mathbf{y}^T \Lambda \mathbf{y} = \mathbf{y}^T\begin{pmatrix} \lambda_1 y_1 \\ \lambda_2 y_2 \\ \vdots \\ \lambda_n y_n \end{pmatrix} = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2
```

$`S`$가 양정치이므로, $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$.

$`\mathbf{x} = Q\mathbf{e}_i`$로 선택하면 ($`\mathbf{e}_i`$는 $`I`$의 $`i`$번째 열), $`\mathbf{y} = Q^T\mathbf{x} = Q^T Q\mathbf{e}_i = \mathbf{e}_i`$이므로, $`\mathbf{x}^T S\mathbf{x} = \lambda_i > 0`$.

따라서 모든 $`\lambda_1, \lambda_2, \ldots, \lambda_n > 0`$. $`\square`$

**비고:** 모든 고유값이 양수라 해서 행렬이 양정치인 것은 아니다.

예: $`A = \begin{pmatrix} 1 & -3 \\ 0 & 1 \end{pmatrix} \to \lambda_1 = \lambda_2 = 1 > 0`$이지만, $`\mathbf{x}^T A\mathbf{x} = (x_1, x_2)\begin{pmatrix} x_1 - 3x_2 \\ x_2 \end{pmatrix} = x_1^2 - 3x_1 x_2 + x_2^2 < 0`$ ($`x_1 = x_2 = 1`$일 때).

(동치 관계가 성립하려면 행렬이 **대칭** 이어야 한다.)

**양의 에너지는 양의 고유값 ($`\lambda > 0`$)과 밀접한 관련이 있다:**

$`S\mathbf{x} = \lambda\mathbf{x}`$이면, $`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T\lambda\mathbf{x} = \lambda\mathbf{x}^T\mathbf{x} = \lambda\|\mathbf{x}\|^2`$.

따라서 $`\lambda > 0 \implies \mathbf{x}^T S\mathbf{x} > 0`$ ($`\mathbf{x} \neq \mathbf{0}`$).

**명제:** $`S`$의 고유벡터에 대해 $`\mathbf{x}^T S\mathbf{x} > 0`$이면, 모든 $`\mathbf{x} \neq \mathbf{0}`$에 대해 $`\mathbf{x}^T S\mathbf{x} > 0`$이다.

**증명:** $`\mathbf{x} = Q\mathbf{c} = c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n`$으로 놓자. 여기서 $`\mathbf{q}_i`$는 $`S`$의 $`i`$번째 고유벡터이자 정규직교 벡터이다.

```math
\mathbf{x}^T S\mathbf{x} = (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T S(c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)
```

```math
= (c_1\mathbf{q}_1 + c_2\mathbf{q}_2 + \cdots + c_n\mathbf{q}_n)^T(c_1\lambda_1\mathbf{q}_1 + c_2\lambda_2\mathbf{q}_2 + \cdots + c_n\lambda_n\mathbf{q}_n)
```

```math
= c_1^2\lambda_1\mathbf{q}_1^T\mathbf{q}_1 + c_2^2\lambda_2\mathbf{q}_2^T\mathbf{q}_2 + \cdots + c_n^2\lambda_n\mathbf{q}_n^T\mathbf{q}_n > 0
```

$`\lambda_1, \lambda_2, \ldots, \lambda_n > 0`$이면 성립. $`\square`$

**성질 2.** $`S`$가 양정치이면, **가역** 이고, $`\det(S) > 0`$이며, $`S^{-1}`$도 양정치이다.

**증명:**

i) $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$

$`\implies S\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$

$`\implies \mathcal{N}(S) = \lbrace\mathbf{0}\rbrace`$

$`\implies S`$는 full rank $`\iff \dim C(S) = n`$

$`\therefore S`$는 **가역** 이다.

ii) $`S`$가 양정치이므로, $`S`$는 모든 양의 고유값을 갖는다.

따라서 $`0 < \lambda_1\lambda_2\cdots\lambda_n = \det(S)`$.

iii) $`(S^{-1})^T = (S^T)^{-1} = S^{-1}`$이므로, $`S^{-1}`$은 **대칭** 이다.

$`\mathbf{x}^T S^{-1}\mathbf{x} = \mathbf{x}^T S^{-1}SS^{-1}\mathbf{x} = (S^{-T}\mathbf{x})^T S(S^{-1}\mathbf{x}) = (S^{-1}\mathbf{x})^T S(S^{-1}\mathbf{x})`$

$`\mathbf{z} = S^{-1}\mathbf{x}`$로 놓으면: $`= \mathbf{z}^T S\mathbf{z} > 0`$. $`\square`$

**성질 3.** $`S`$가 양정치이면, $`S = A^T A`$이고 $`A`$의 열은 독립이다.

**증명:** 대칭 행렬 $`S = Q\Lambda Q^T`$. $`S`$가 양정치이므로 모든 고유값이 양수이다. 대각 행렬 $`\Lambda`$는 제곱근 $`\sqrt{\Lambda}`$를 갖는다:

```math
\Lambda = \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots \\ & & & \lambda_n \end{pmatrix} = \begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix}\begin{pmatrix} \sqrt{\lambda_1} & & \\ & \sqrt{\lambda_2} & \\ & & \ddots \\ & & & \sqrt{\lambda_n} \end{pmatrix} = \sqrt{\Lambda}\sqrt{\Lambda}
```

$`A = Q\sqrt{\Lambda}Q^T`$로 놓으면:

```math
A^T A = (Q\sqrt{\Lambda}Q^T)(Q\sqrt{\Lambda}Q^T) = Q\sqrt{\Lambda}\underbrace{Q^T Q}_{I}\sqrt{\Lambda}Q^T = Q\sqrt{\Lambda}\sqrt{\Lambda}Q^T = Q\Lambda Q^T = S
```

**비고:** 에너지 $`= \mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2`$.

$`\|A\mathbf{x}\| > 0`$ ($`A\mathbf{x} \neq \mathbf{0} \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$인 경우) $`\iff`$ $`A`$가 full rank.

### 3.6 양정치 행렬 판별법

**(1) 양의 고유값:** 양정치 행렬 $`S`$는 모든 양의 고유값을 갖는다.

예: $`S = \mathbf{u}\mathbf{v}^T`$ (랭크 1 행렬)은 양정치가 아니다.

$`A = \begin{pmatrix} a \\ b \end{pmatrix}(1 \; k) = \begin{pmatrix} a & ka \\ b & kb \end{pmatrix}`$

랭크 1 행렬 $`\to`$ 1개 일차독립 벡터 $`\to`$ $`(n-1)`$개 종속 벡터.

$`A\mathbf{u} = \mathbf{u}\mathbf{v}^T\mathbf{u} = (\mathbf{v}^T\mathbf{u})\mathbf{u}`$이므로, $`\lambda = \mathbf{v}^T\mathbf{u} = (1, k)\begin{pmatrix} a \\ b \end{pmatrix} = a + kb`$.

$`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_{n-1}`$를 $`\mathcal{N}(A)`$의 기저 벡터라 하자. $`A\mathbf{w}_i = \mathbf{0} = 0 \cdot \mathbf{w}_i`$이므로, $`i = 1, 2, \ldots, n-1`$에 대해 $`\lambda = 0`$.

$`\therefore (n-1)`$개의 영 고유값.

**(2) 양의 에너지 판정법:** $`\mathbf{x}^T S\mathbf{x} > 0 \; \forall \; \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$.

**(3) $`S = A^T A`$에 대한 양의 에너지 판정법:**

```math
\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2
```

예: $`S = A^T A`$ ($`A = \begin{pmatrix} 1 & 1 & 2 \\ 1 & 2 & 3 \end{pmatrix}`$, 열 $`\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3`$).

참고: $`\mathbf{a}_3 = \mathbf{a}_1 + \mathbf{a}_2`$이므로, $`\text{rank}(A) = 2 < 3`$, $`\dim C(A^T) = 2`$, $`\dim \mathcal{N}(A) = 1`$.

$`\dim \mathcal{N}(A) \neq 0`$이므로, $`\mathbf{x}^T S\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0`$. $`S`$는 양정치가 아니다. $`S`$는 **양반정치** 이다.

**(4) 행렬식 판정법:** $`\det(S) > 0`$인지 확인 — 더 정확하게는, 모든 **왼쪽 상단 행렬식**(leading determinant)을 확인한다.

```math
S = \begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & -1 & 2 & -1 \\ & & -1 & 2 \end{pmatrix}
```

$`D_1 = |2| = 2 > 0`$

$`D_2 = \begin{vmatrix} 2 & -1 \\ -1 & 2 \end{vmatrix} = 3 > 0`$

$`D_3 = \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{vmatrix} = 8 - 4 = 4 > 0`$

$`D_4 = |S| = 2D_3 + \begin{vmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & -1 & -1 \end{vmatrix} = 8 - 4 + 1 = 5 > 0`$

모든 왼쪽 상단 행렬식이 양수이므로, $`S`$는 양정치이다.

**(5) 피벗 판정법:** **피벗** 이 양수인지 확인한다.

$`S = \begin{pmatrix} 2 & -1 \\ -1 & 2 & -1 \\ & -1 & 2 \end{pmatrix}`$의 경우:

소거 후: 피벗은 $`2, \frac{3}{2}, \frac{4}{3}`$ — 모두 양수.

$`S = LU`$, $`S = LDL^T`$, $`S = A^T A`$.

선행 행렬식(leading determinant)은 피벗과 밀접한 관련이 있다: $`D_2/D_1`$, $`D_3/D_1`$ 등.

**$`2 \times 2`$ SPD 행렬** $`S = \begin{pmatrix} a & b \\ b & d \end{pmatrix}`$에 대해:

- 행렬식 판정: $`D_1 = a > 0`$, $`D_2 = ad - b^2 > 0`$.
- 피벗 판정: $`d_1 = a > 0`$, $`d_2 = \frac{ad - b^2}{a} > 0`$.
- 고유값: $`\lambda_1 > 0`$, $`\lambda_2 > 0`$.
- 에너지: $`(x \; y)\begin{pmatrix} a & b \\ b & d \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = ax^2 + bxy + byx + dy^2 = ax^2 + 2bxy + dy^2 > 0`$.

### 3.7 풀이 예제: 양정치와 양반정치

**예제:** $`S = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix} \implies \lambda_1 = 2 > 0, \lambda_2 = 6 > 0`$

$`S\mathbf{x} = \begin{pmatrix} 2 & 0 \\ 0 & 6 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix}`$

$`\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 2x_1 \\ 6x_2 \end{pmatrix} = 2x_1^2 + 6x_2^2 > 0`$. $`\therefore S`$는 양정치이다.

**예제:** $`S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}`$

$`|S - \lambda I| = \begin{vmatrix} 5 - \lambda & 4 \\ 4 & 5 - \lambda \end{vmatrix} = \lambda^2 - 10\lambda + 25 - 16 = (\lambda - 9)(\lambda - 1) = 0`$

$`\therefore \lambda_1 = 9 > 0`$, $`\lambda_2 = 1 > 0`$.

i) $`(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 = x_2`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 + x_2 = 0`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \to \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) $`Q = (\mathbf{x}_1 \; \mathbf{x}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}`$, $`Q^{-1} = Q = Q^T`$

iv) $`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T(Q\Lambda Q^T)\mathbf{x}`$, $`\mathbf{y} = Q^T\mathbf{x}`$로 놓으면:

$`= \lambda_1 y_1^2 + \lambda_2 y_2^2 = 9y_1^2 + y_2^2 > 0`$

$`\therefore S`$는 양정치이다.

**대안적 접근법 (에너지 판정법):**

$`\mathbf{x}^T S\mathbf{x} = (x_1, x_2)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = (x_1, x_2)\begin{pmatrix} 5x_1 + 4x_2 \\ 4x_1 + 5x_2 \end{pmatrix}`$

$`= 5x_1^2 + 4x_1x_2 + 4x_1x_2 + 5x_2^2 = x_1^2 + 4x_1^2 + 8x_1x_2 + 4x_2^2 + x_1^2`$

$`= x_1^2 + \underbrace{4(x_1^2 + 2x_1x_2 + x_2^2)}_{4(x_1 + x_2)^2} > 0`$

**예제:** $`S = \begin{pmatrix} 4 & 5 \\ 5 & 4 \end{pmatrix}`$

$`|S - \lambda I| = \lambda^2 - 8\lambda + 16 - 25 = \lambda^2 - 8\lambda - 9 = (\lambda - 9)(\lambda + 1) = 0`$

$`\therefore \lambda_1 = 9`$이고 $`\lambda_2 = -1`$.

$`\lambda_2 < 0`$이므로, $`S`$는 양정치가 **아니다**.

### 3.8 타원과 이차 형식

**타원**(The Ellipse) $`ax^2 + 2bxy + cy^2 = 1`$.

**예제 1:** 타원 $`5x^2 + 8xy + 5y^2 = 1`$을 고려하자.

```math
(x \; y)\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 1 \implies \mathbf{x}^T S\mathbf{x} = 1
```

$`S = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}`$, $`\lambda^2 - 10\lambda + 9 = (\lambda - 9)(\lambda - 1) = 0`$이므로, $`\lambda = 9, 1`$.

i) $`(S - 9I)\mathbf{x}_1 = \begin{pmatrix} -4 & 4 \\ 4 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`(S - I)\mathbf{x}_2 = \begin{pmatrix} 4 & 4 \\ 4 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) $`\mathbf{q}_1 = \frac{1}{\sqrt{2}}\mathbf{x}_1`$, $`\mathbf{q}_2 = \frac{1}{\sqrt{2}}\mathbf{x}_2`$

$`Q = (\mathbf{q}_1 \; \mathbf{q}_2) = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}`$, $`Q^{-1} = Q^T = Q`$

$`S = Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T`$

iv) $`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T Q\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}Q^T\mathbf{x}`$

$`\mathbf{z} = Q^T\mathbf{x}`$로 놓으면: $`= \mathbf{z}^T\begin{pmatrix} 9 & \\ & 1 \end{pmatrix}\mathbf{z} = (z_1, z_2)\begin{pmatrix} 9z_1 \\ z_2 \end{pmatrix} = 9z_1^2 + z_2^2`$

$`\mathbf{z} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix} x_1 + x_2 \\ x_1 - x_2 \end{pmatrix}`$

$`\therefore z_1 = \frac{x_1 + x_2}{\sqrt{2}}, \; z_2 = \frac{x_1 - x_2}{\sqrt{2}}`$

$`\mathbf{x}^T S\mathbf{x} = 9z_1^2 + z_2^2 = 1`$

$`z_2 = 0`$일 때: $`z_1^2 = \frac{1}{9}`$이므로, $`z_1 = \pm \frac{1}{3}`$ ($`\mathbf{q}_1`$ 방향의 반축).

$`z_1 = 0`$일 때: $`z_2^2 = 1`$이므로, $`z_2 = \pm 1`$ ($`\mathbf{q}_2`$ 방향의 반축).

**예제 2:** 양반정치

```math
T = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}\begin{pmatrix} 3 & 1 \end{pmatrix} = A^T A
```

행렬식: $`D_1 = 9 > 0`$, $`D_2 = 9 - 9 = 0 \to`$ 역행렬 없음.

$`\lambda^2 - \text{trace}(T)\lambda + \det(T) = \lambda^2 - 10\lambda = 0`$이므로, $`\lambda_1 = 10, \lambda_2 = 0`$.

$`(T - 10I)\mathbf{x}_1 = \begin{pmatrix} -1 & 3 \\ 3 & -9 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 3 \\ 1 \end{pmatrix}`$

$`(T - 0I)\mathbf{x}_2 = \begin{pmatrix} 9 & 3 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ -3 \end{pmatrix}`$

에너지: $`ax^2 + 2bxy + dy^2 = 9x^2 + 6xy + y^2 = (3x + y)^2 \geq 0`$.

$`E(x,y) = 1`$은 띠(band)이다: $`3x + y = \pm 1`$. 에너지 $`E(x,y)`$의 그래프는 계곡(valley)이다. 띠의 축은 $`T`$의 고유벡터 방향을 따른다.

### 3.9 양정치 행렬과 최솟값 문제

에너지 $`E = \mathbf{x}^T S\mathbf{x} = \begin{pmatrix} x \\ y \end{pmatrix}^T\begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = 5x^2 + 8xy + 5y^2 > 0`$을 고려하자.

(그릇 모양, **볼록**(convex))

$`E(x,y) = 0`$일 때 $`x = y = 0`$이다. 이는 **최솟값 문제** 와 연결된다.

이계 도함수의 행렬은 모든 점에서 양정치이다.

1차 도함수: $`\frac{\partial E}{\partial x} = 10x + 8y`$, $`\frac{\partial E}{\partial y} = 8x + 10y`$.

2차 도함수: $`\frac{\partial^2 E}{\partial x^2} = 10 > 0`$, $`\frac{\partial^2 E}{\partial x \partial y} = 8 > 0`$, $`\frac{\partial^2 E}{\partial y^2} = 10 > 0`$.

```math
\nabla E = \left(\frac{\partial E}{\partial x}, \frac{\partial E}{\partial y}\right) \quad \xrightarrow{(x,y) = (0,0)} \quad \nabla E = (0, 0)
```

```math
J(\nabla E^T) = \begin{pmatrix} \frac{\partial^2 E}{\partial x^2} & \frac{\partial^2 E}{\partial y \partial x} \\ \frac{\partial^2 E}{\partial x \partial y} & \frac{\partial^2 E}{\partial y^2} \end{pmatrix} = H \quad \text{(헤시안 행렬, Hessian matrix)}
```

$`H_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}`$

에너지를 $`\frac{1}{2}\mathbf{x}^T S\mathbf{x}`$로 정의하면, $`H = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix} = S`$.

예: $`E = x^2 + y^2`$는 $`x = y = 0`$에서 최솟값을 갖는다.

$`f = e^{x^2 + y^2}`$는 어떠한가?

$`\nabla f = \left(\frac{\partial f}{\partial x}, \frac{\partial f}{\partial y}\right)`$이고, $`\frac{\partial f}{\partial x} = e^{x^2+y^2} \cdot 2x`$, $`\frac{\partial f}{\partial y} = e^{x^2+y^2} \cdot 2y`$.

$`J(\nabla f)^T = \begin{pmatrix} (2 + 4x^2)e^{x^2+y^2} & 4xy \, e^{x^2+y^2} \\ 4xy \, e^{x^2+y^2} & (2 + 4y^2)e^{x^2+y^2} \end{pmatrix} = S`$

$`S`$가 양정치이므로, $`f`$는 **순볼록**(strictly convex)이다.

**예제 1 (양정치):**

$`S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix}`$

행렬식: $`D_1 = 9 > 0`$, $`D_2 = 27 - 9 = 18 > 0`$.

피벗: $`\begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} \xrightarrow{R_2 - \frac{1}{3}R_1} \begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix}`$. 피벗: 9, 2 (둘 다 양수).

에너지: $`E(x,y) = 9x^2 + 6xy + 3y^2 = (3x + y)^2 + 2y^2 > 0`$. (순볼록 함수.)

대각합과 행렬식: $`\text{trace}(S) = 12`$, $`\det(S) = 18`$.

$`\lambda^2 - 12\lambda + 18 = 0 \implies \lambda = 6 \pm \sqrt{36 - 18} = 6 \pm 3\sqrt{2}`$

i) $`\lambda_1 = 6 + 3\sqrt{2}`$: $`(S - \lambda_1 I)\mathbf{x} = \begin{pmatrix} 3 - 3\sqrt{2} & 3 \\ 3 & 3 - 3\sqrt{2} \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`x_2 = (-1 + \sqrt{2})x_1`$, $`x_1 = 1`$로 선택: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ -1 + \sqrt{2} \end{pmatrix}`$

ii) $`\lambda_2 = 6 - 3\sqrt{2}`$: $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 - \sqrt{2} \end{pmatrix}`$

**분해:**

$`S = \begin{pmatrix} 9 & 3 \\ 3 & 3 \end{pmatrix} = LU = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & 3 \\ 0 & 2 \end{pmatrix} = LDL^T = \begin{pmatrix} 1 & 0 \\ \frac{1}{3} & 1 \end{pmatrix}\begin{pmatrix} 9 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & \frac{1}{3} \\ 0 & 1 \end{pmatrix}`$

$`= (L\sqrt{D})(\sqrt{D}L^T) = A^T A`$ ($`A = \begin{pmatrix} 3 & 1 \\ 0 & \sqrt{2} \end{pmatrix}`$)

### 3.10 양반정치 행렬

$`\mathbf{x}^T S\mathbf{x} \geq 0`$

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}
```

$`\det(S) = 4 - 4 = 0`$ — 특이 행렬.

$`\text{trace}(S) = 1 + 4 = 5`$. $`\lambda^2 - 5\lambda + 0 = \lambda(\lambda - 5) = 0`$이므로, $`\lambda_1 = 5, \lambda_2 = 0`$.

i) $`(S - 5I)\mathbf{x}_1 = \begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_2 = 2x_1`$: $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 2 \end{pmatrix}`$

ii) $`S\mathbf{x}_2 = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 + 2x_2 = 0`$: $`\mathbf{x}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}`$

$`E_{21}S = U`$: $`S = E_{21}^{-1}U = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = LDL^T`$

$`= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}`$

$`= \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}`$

$`S`$는 $`A`$의 **종속** 열을 갖는 $`A^T A`$로 분해된다.

$`\mathbf{x}^T S\mathbf{x} = \mathbf{x}^T A^T A\mathbf{x} = (A\mathbf{x})^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0`$

### 3.11 합동 행렬

$`S`$가 양반정치이면, 모든 행렬 $`A^T SA`$도 양반정치이다:

```math
\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} \geq 0 \quad \forall \mathbf{x}
```

여기서 $`\mathbf{y} = A\mathbf{x}`$. $`\square`$

$`\mathbf{x}^T S\mathbf{x} > 0`$이고, 모든 $`\mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace`$에 대해 $`A\mathbf{x} \neq \mathbf{0}`$이면, $`A^T SA`$는 **양정치** 이다.

```math
\mathbf{x}^T(A^T SA)\mathbf{x} = (A\mathbf{x})^T S(A\mathbf{x}) = \mathbf{y}^T S\mathbf{y} > 0 \quad \forall \mathbf{x} \in \mathbb{R}^n \setminus \lbrace\mathbf{0}\rbrace \quad \square
```

행렬 $`A^T SA`$를 $`S`$에 **"합동"**(congruent)이라 한다.

$`S^T = S`$가 $`P`$개의 양의 고유값, $`N`$개의 음의 고유값, $`Z`$개의 영 고유값을 가지면, $`A`$가 가역인 경우 $`A^T SA`$에 대해서도 동일하다.

**반대칭 행렬의 고유값 성질 증명:**

$`A = -\bar{A}^T`$. $`A\mathbf{x} = \lambda\mathbf{x}`$, $`\overline{A\mathbf{x}} = \bar{\lambda}\bar{\mathbf{x}}`$.

$`\bar{\mathbf{x}}^T A\mathbf{x} = \bar{\mathbf{x}}^T(-\bar{A}^T)\mathbf{x} = -(\bar{A}\bar{\mathbf{x}})^T\mathbf{x}`$

$`= \bar{\mathbf{x}}^T\lambda\mathbf{x}`$ 이고 $`= -(\bar{\lambda}\bar{\mathbf{x}})^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}`$

$`\lambda\bar{\mathbf{x}}^T\mathbf{x} = -\bar{\lambda}\bar{\mathbf{x}}^T\mathbf{x}`$

$`(\lambda + \bar{\lambda})\bar{\mathbf{x}}^T\mathbf{x} = 0`$

$`\therefore \lambda + \bar{\lambda} = (r + is) + (r - is) = 2r = 0`$

$`\lambda`$의 실수 부분은 **영** 이다. 반대칭 행렬은 영 또는 **순허수** 고유값을 갖는다.

### 3.12 최적화와 머신 러닝

$`f(x)`$를 최소화하기 위한 **경사 하강법**(gradient descent).

i) $`f(x) = x^2 + 4x + 4 \implies f'(x) = 2x + 4`$

ii) $`x_0 = 10`$, $`\eta = 0.1`$ (학습률, learning rate), 정지 기준 $`|f'(x)| < 0.01`$

iii) $`x_{k+1} = x_k - \eta f'(x_k)`$ — **최급강하 방향**(steepest direction).

$`x_k \to x^*`$ (최솟값 점)까지 근사를 반복한다.

$`f(\mathbf{x})`$에 대해: $`\mathbf{x}_{k+1} = \mathbf{x}_k - \eta \nabla f(\mathbf{x}_k)`$

$`\mathbf{x}_k \to \mathbf{x}^*`$이면, $`\mathbf{x}^*`$에서 $`\nabla f`$가 영 $`\implies \mathbf{x}_{k+1} \approx \mathbf{x}_k`$.

$`J(\nabla f)^T = S`$가 양정치이면, $`f`$는 **볼록**(convex) — $`\mathbf{x}^*`$를 찾기 쉽다.

---

<br>
