# 제4장 강의 — 직교성 — Part 2/2

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.**
>
> [Part 1](Concepts.ko.md) · **Part 2**

---

<br>

## 목차

- [3. 최소제곱 근사 (4.3)](#3-최소제곱-근사-43)
- [4. 직교 기저와 그람-슈미트 (4.4)](#4-직교-기저와-그람-슈미트-44)
- [5. 행렬의 유사역행렬 (4.5)](#5-행렬의-유사역행렬-45)
- [요약](#요약)

---

<br>

## 3. 최소제곱 근사 (4.3)

### 3.1 핵심 사실 요약 (Key Facts Summary)

**(1)** $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$를 풀면 $`\mathbf{b}`$의 $`A`$의 열 공간으로의 사영 $`\mathbf{p} = A\hat{\boldsymbol{\alpha}}`$를 얻는다.

**(2)** $`A\mathbf{x} = \mathbf{b}`$에 해가 없을 때, $`\hat{\boldsymbol{\alpha}}`$는 "최소제곱 해"이다: $`\|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2 = \text{최솟값}`$.

**(3)** $`E = \|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2`$의 편미분을 0으로 놓으면 ($`\frac{\partial E}{\partial \alpha_i} = 0`$) 역시 $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$를 만든다.

**(4)** 점 $`(t_1, b_1), (t_2, b_2), \ldots, (t_m, b_m)`$을 직선으로 적합할 때, $`A`$는 열 $`(1, 1, \ldots, 1)`$과 $`(t_1, t_2, \ldots, t_m)`$을 갖는다.

**(5)** 이 경우 $`A^T A`$는 $`2 \times 2`$ 행렬 $`\begin{pmatrix} m & \sum t_i \\ \sum t_i & \sum t_i^2 \end{pmatrix}`$이고 $`A^T\mathbf{b}`$는 $`\begin{pmatrix} \sum b_i \\ \sum t_i b_i \end{pmatrix}`$이다.

### 3.2 과잉결정 시스템 (Overdetermined Systems)

$`A\mathbf{x} = \mathbf{b}`$에 해가 없는 경우가 자주 발생한다. 행렬 $`A`$는 열보다 행이 더 많다 ($`m > n`$). 방정식 수가 미지수 수보다 많다. $`n`$개의 열은 $`m`$차원 공간의 작은 부분만 생성한다.

$`\mathbf{b}`$는 $`C(A)`$ 밖에 있을 수 있다. 이 경우에도 $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$를 풀어 $`A\mathbf{x} = \mathbf{b}`$의 근사를 찾을 수 있다.

$`\hat{\boldsymbol{\alpha}}`$는 다음을 최소화하므로 **최소제곱 해** 이다:

```math
E = \|A\hat{\boldsymbol{\alpha}} - \mathbf{b}\|^2
```

### 3.3 예제 1: 직선 적합 (Fitting a Line)

점 $`\begin{pmatrix} 0 \\ 6 \end{pmatrix}`$, $`\begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\begin{pmatrix} 2 \\ 0 \end{pmatrix}`$에 가장 가까운 직선을 찾아라.

```math
y = ax + b
```

```math
y_1 = ax_1 + b \Rightarrow 6 = a \cdot 0 + b
```

```math
y_2 = ax_2 + b \Rightarrow 0 = a \cdot 1 + b
```

```math
y_3 = ax_3 + b \Rightarrow 0 = a \cdot 2 + b
```

```math
\begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 1 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}
```

불행히도, $`\mathbf{b} \notin C(A)`$. $`A\mathbf{x} = \mathbf{b}`$는 풀 수 없다.

대신, $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$를 풀어 $`\hat{\mathbf{x}}`$에 대한 근사 $`\hat{\boldsymbol{\alpha}}`$를 찾는다. 실제로, $`\hat{\boldsymbol{\alpha}}`$는 $`E`$의 최솟값을 준다.

### 3.4 오차 최소화: 미적분 접근 (Minimizing the Error: Calculus Approach)

```math
E = \|\mathbf{y} - \hat{\mathbf{y}}\|^2 = (y_1 - \hat{y}_1)^2 + (y_2 - \hat{y}_2)^2 + (y_3 - \hat{y}_3)^2 = \sum_{i=1}^{3}(y_i - \hat{y}_i)^2
```

여기서 $`\hat{y} = ax + b`$.

```math
\frac{\partial E}{\partial a} = 2\sum_{i=1}^{3}(y_i - \hat{y}_i)x_i = 0
```

```math
\frac{\partial E}{\partial b} = 2\sum_{i=1}^{3}(y_i - \hat{y}_i) \cdot 1 = 0
```

$`\sum(y_i - ax_i - b)x_i = 0`$으로부터: $`\sum y_i x_i = a\sum x_i^2 + b\sum x_i`$

$`\sum(y_i - \hat{y}_i) = 0`$으로부터: $`\sum y_i = a\sum x_i + b\sum 1`$

```math
\begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix} = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}
```

계산:
- $`\sum x_i = 0 + 1 + 2 = 3`$
- $`\sum x_i^2 = 0 + 1 + 4 = 5`$
- $`\sum y_i = 6 + 0 + 0 = 6`$
- $`\sum x_i y_i = 6 \cdot 0 + 0 \cdot 1 + 0 \cdot 2 = 0`$

```math
\begin{pmatrix} 0 \\ 6 \end{pmatrix} = \begin{pmatrix} 5 & 3 \\ 3 & 3 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix}
```

```math
\begin{pmatrix} a \\ b \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 3 & -3 \\ -3 & 5 \end{pmatrix}\begin{pmatrix} 0 \\ 6 \end{pmatrix} = \begin{pmatrix} -3 \\ 5 \end{pmatrix}
```

따라서 최적 직선은 $`\hat{y} = -3x + 5`$이다.

### 3.5 직선 적합의 일반 정규방정식 (General Normal Equations for Line Fitting)

```math
A^T A = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}, \quad A^T\mathbf{b} = \begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix}
```

정규방정식 $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$:

```math
\begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & \sum 1 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} \sum x_i y_i \\ \sum y_i \end{pmatrix}
```

### 3.6 오차 최소화: 세 가지 접근 (Minimizing the Error: Three Approaches)

**오차 $`\mathbf{e} = \mathbf{b} - A\hat{\boldsymbol{\alpha}}`$를 어떻게 가능한 한 작게 만들 수 있는가?**

**(1) 기하학적:** $`\mathbf{b}`$에 가장 가까운 점을 찾는다. $`\|\mathbf{e}\|`$는 $`\mathbf{e} \perp \mathbf{a}`$, 즉 $`\mathbf{e} \perp C(A)`$일 때 최솟값이 된다.

**(2) 대수적:** $`A\mathbf{x} = \mathbf{b} = \mathbf{p} + \mathbf{e}`$는 풀 수 없다. $`A\hat{\boldsymbol{\alpha}} = \mathbf{p}`$는 풀 수 있다. $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{p}`$이므로, $`\hat{\boldsymbol{\alpha}} = (A^T A)^{-1}A^T\mathbf{p}`$.

**(3) 임의의 $`\mathbf{x}`$에 대한 제곱 오차:**

```math
\|A\mathbf{x} - \mathbf{b}\|^2 = \|A\mathbf{x} - \mathbf{p} - \mathbf{e}\|^2
```

```math
= (A\mathbf{x} - \mathbf{p} - \mathbf{e})^T(A\mathbf{x} - \mathbf{p} - \mathbf{e})
```

$`A\mathbf{x} - \mathbf{p} \in C(A)`$이고 $`\mathbf{e} \in \mathcal{N}(A^T)`$이므로, 교차항이 사라진다:

```math
= \|A\mathbf{x} - \mathbf{p}\|^2 + \|\mathbf{e}\|^2
```

$`\mathbf{x} = \hat{\boldsymbol{\alpha}}`$로 선택하여 $`A\hat{\boldsymbol{\alpha}} - \mathbf{p} = \mathbf{0}`$이 되면. $`A\mathbf{x} - \mathbf{b}`$의 제곱 길이가 최소화된다:

```math
\|A\mathbf{x} - \mathbf{b}\|^2 = \|\mathbf{e}\|^2
```

최소제곱 해 $`\hat{\boldsymbol{\alpha}}`$는 $`E = \|A\mathbf{x} - \mathbf{b}\|^2`$를 가능한 한 작게 만든다.

**(4) 미적분:**

```math
E = \|A\mathbf{x} - \mathbf{b}\|^2 = (A\mathbf{x} - \mathbf{b})^T(A\mathbf{x} - \mathbf{b})
```

```math
= (A\mathbf{x})^T A\mathbf{x} - (A\mathbf{x})^T\mathbf{b} - \mathbf{b}^T(A\mathbf{x}) + \mathbf{b}^T\mathbf{b}
```

```math
= \mathbf{x}^T A^T A\mathbf{x} - 2\mathbf{x}^T A^T\mathbf{b} + \mathbf{b}^T\mathbf{b}
```

첨자 표기법 사용: $`E = x_i(A^T A)_{ij}x_j - 2x_i(A^T\mathbf{b})_i + b_ib_i`$.

```math
\frac{\partial E}{\partial x_i} = 2(A^T A)_{ij}x_j - 2(A^T\mathbf{b})_i = 0
```

이것은 $`A^T A\mathbf{x} = A^T\mathbf{b}`$를 준다.

$`E = \|A\mathbf{x} - \mathbf{b}\|^2`$의 편미분은 $`A^T A\mathbf{x} = A^T\mathbf{b}`$일 때 0이다.

### 3.7 최소제곱의 큰 그림 (The Big Picture for Least Squares)

$`\mathbf{b}`$를 $`\mathbf{p}`$와 $`\mathbf{e}`$로 분리:

- $`\hat{\boldsymbol{\alpha}} \in C(A^T)`$ (행 공간), $`A\hat{\boldsymbol{\alpha}} = \mathbf{p}`$는 풀 수 있다
- $`\mathbf{p} = P\mathbf{b}`$는 $`C(A)`$에서 $`\mathbf{b}`$에 가장 가까운 점
- $`\mathbf{e}`$는 $`\mathcal{N}(A^T)`$에서의 최소 오차
- $`\mathbf{b} \notin C(A)`$이면, $`A\mathbf{x} = \mathbf{b}`$는 풀 수 없다

### 3.8 예제 2: 또 다른 직선 적합 (Another Line Fit)

점 $`(-2, 1)`$, $`(0, 2)`$, $`(2, 4)`$가 주어졌다. 최소제곱 오차를 최소화하는 직선을 찾아라.

```math
\hat{y} = ax + b
```

```math
1 = a(-2) + b, \quad 2 = a(0) + b, \quad 4 = a(2) + b
```

```math
\begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 1 & -2 \\ 1 & 0 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} b \\ a \end{pmatrix}
```

```math
A^T A = \begin{pmatrix} 1 & 1 & 1 \\ -2 & 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & -2 \\ 1 & 0 \\ 1 & 2 \end{pmatrix} = \begin{pmatrix} 3 & 0 \\ 0 & 8 \end{pmatrix}
```

```math
A^T\mathbf{b} = \begin{pmatrix} 1 & 1 & 1 \\ -2 & 0 & 2 \end{pmatrix}\begin{pmatrix} 1 \\ 2 \\ 4 \end{pmatrix} = \begin{pmatrix} 7 \\ 6 \end{pmatrix}
```

```math
\begin{pmatrix} 3 & 0 \\ 0 & 8 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} 7 \\ 6 \end{pmatrix}
```

```math
\alpha_1 = 7/3, \quad \alpha_2 = 3/4
```

참고: $`\mathbf{a}_2 \perp \mathbf{a}_1`$ (직교 열 벡터)이므로 $`A^T A`$가 대각 행렬이 된다.

### 3.9 $`A`$에 종속 열이 있을 때: $`\hat{\boldsymbol{\alpha}}`$는? (Dependent Columns in $`A`$: What is $`\hat{\boldsymbol{\alpha}}`$?)

$`A`$에 종속 열이 있으면 어떤 $`\hat{\boldsymbol{\alpha}}`$가 최적인가?

```math
A\mathbf{x} = \mathbf{b}: \quad \begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 3 \\ 3 \end{pmatrix}, \quad \mathbf{b} \notin C(A)
```

```math
\mathbf{p} = \begin{pmatrix} 2 \\ 2 \\ 2 \end{pmatrix} \in C(A)
```

$`A\hat{\boldsymbol{\alpha}} = \mathbf{p}`$는 **많은 해** 를 갖는다: $`\begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} = \begin{pmatrix} 2 \\ 2 \\ 2 \end{pmatrix}`$.

문제는 $`A`$가 종속 열을 가지고 $`\mathcal{N}(A) \ni \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$이라는 것이다.

최적 해를 어떻게 찾을 수 있는가? 4.5절에서 $`A`$의 "유사역행렬" (pseudoinverse)을 배운다.

### 3.10 포물선 적합 (Fitting by a Parabola)

높이 $`b_1, b_2, \ldots, b_m`$을 시간 $`t_1, t_2, \ldots, t_m`$에서 포물선으로 적합하려 한다:

```math
\hat{b} = \alpha_1 + \alpha_2 t + \alpha_3 t^2
```

$`m > 3`$일 때:

```math
b_1 = \alpha_1 + \alpha_2 t_1 + \alpha_3 t_1^2, \quad b_2 = \alpha_1 + \alpha_2 t_2 + \alpha_3 t_2^2, \quad \ldots, \quad b_m = \alpha_1 + \alpha_2 t_m + \alpha_3 t_m^2
```

```math
\begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{pmatrix} = \begin{pmatrix} 1 & t_1 & t_1^2 \\ 1 & t_2 & t_2^2 \\ \vdots & \vdots & \vdots \\ 1 & t_m & t_m^2 \end{pmatrix}\begin{pmatrix} \alpha_1 \\ \alpha_2 \\ \alpha_3 \end{pmatrix}
```

최소제곱을 사용하여 $`\hat{\boldsymbol{\alpha}}`$를 찾을 수 있다: $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$.

### 3.11 예제 3: 포물선 적합 (Parabola Fit)

$`\hat{b}(t) = \alpha_1 + \alpha_2 t + \alpha_3 t^2`$. $`t = 0, 1, 2`$에서 세 높이 $`b_1, b_2, b_3`$이 측정되었으며, $`b_1 = 6`$, $`b_2 = 0`$, $`b_3 = 0`$이다.

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 1 \\ 1 & 2 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 6 \\ 0 \\ 0 \end{pmatrix}
```

$`A\hat{\boldsymbol{\alpha}} = \mathbf{b}`$. $`\text{rank}(A) = 3`$이므로, $`\hat{\boldsymbol{\alpha}} = A^{-1}\mathbf{b}`$.

$`(A|\mathbf{b})`$의 행 축소:

```math
\begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 1 & 1 & 1 & | & 0 \\ 1 & 2 & 4 & | & 0 \end{pmatrix} \xrightarrow{R_2 - R_1, \, R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 2 & 4 & | & -6 \end{pmatrix}
```

```math
\xrightarrow{R_3/2} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 1 & 2 & | & -3 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 1 & | & -6 \\ 0 & 0 & 1 & | & 3 \end{pmatrix}
```

```math
\xrightarrow{R_2 - R_3} \begin{pmatrix} 1 & 0 & 0 & | & 6 \\ 0 & 1 & 0 & | & -9 \\ 0 & 0 & 1 & | & 3 \end{pmatrix} = (I | A^{-1}\mathbf{b})
```

```math
\therefore \hat{b}(t) = 6 - 9t + 3t^2
```

---

<br>

## 4. 직교 기저와 그람-슈미트 (4.4)

### 4.1 핵심 사실 요약 (Key Facts Summary)

**(1)** 열 $`\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n`$이 **정규직교** (orthonormal)일 조건:

```math
\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 0 & i \neq j\text{일 때 (직교 벡터)} \\ 1 & i = j\text{일 때 (단위 벡터: } \|\mathbf{q}_i\| = 1) \end{cases}
```

이때 $`Q^T Q = I_{n \times n}`$.

**(2)** $`Q`$가 정사각이면, $`QQ^T = I`$이고 $`Q^T = Q^{-1}`$이다. 이때 $`Q`$는 **"직교 행렬"** (orthogonal matrix)이다.

**(3)** $`Q\mathbf{x} = \mathbf{b}`$의 최소제곱 해는 $`\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}`$이다. $`\mathbf{b}`$의 사영:

```math
\mathbf{p} = QQ^T\mathbf{b} = P\mathbf{b} = \mathbf{q}_1 d_1 + \mathbf{q}_2 d_2 + \cdots + \mathbf{q}_n d_n
```

**(4)** 그람-슈미트(Gram-Schmidt) 과정은 독립인 $`\mathbf{a}_i`$를 직교 $`\mathbf{q}_i`$로 변환한다. $`\mathbf{q}_1 = \mathbf{r}_1/\|\mathbf{r}_1\|`$로 시작한다.

**(5)** $`\mathbf{q}_i`$는 $`(\mathbf{a}_i - \text{사영 } \mathbf{p}_i) / \|\mathbf{a}_i - \mathbf{p}_i\|`$이며, 사영 $`\mathbf{p}_i = (\mathbf{a}_i^T\mathbf{q}_1)\mathbf{q}_1 + \cdots + (\mathbf{a}_i^T\mathbf{q}_{i-1})\mathbf{q}_{i-1}`$이다.

**(6)** 각 $`\mathbf{a}_i`$는 $`\mathbf{q}_1`$에서 $`\mathbf{q}_n`$까지의 결합이 된다. 그러면 $`A = QR`$: 직교 $`Q`$와 삼각 $`R`$.

### 4.2 목표: 직교 열 (Goal: Orthogonal Columns)

**i)** $`A`$에서 직교 열이 좋다.

```math
A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \cdots \; \mathbf{a}_n)
```

```math
\mathbf{a}_1^T\mathbf{a}_2 = 0, \quad \ldots, \quad \mathbf{a}_1^T\mathbf{a}_n = 0
```

```math
A^T A = \begin{pmatrix} \mathbf{a}_1^T\mathbf{a}_1 & 0 & \cdots & 0 \\ 0 & \mathbf{a}_2^T\mathbf{a}_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \mathbf{a}_n^T\mathbf{a}_n \end{pmatrix}
```

이는 **대각 행렬** 이다.

**ii)** 그람-슈미트 과정을 통해 직교 벡터 $`\mathbf{q}_i`$를 구성한다.

### 4.3 정규직교 벡터와 행렬 (Orthonormal Vectors and Matrices)

**정의:** $`n`$개의 벡터가 **정규직교** 일 조건:

```math
\mathbf{q}_i^T\mathbf{q}_j = \begin{cases} 0 & i \neq j\text{일 때} \\ 1 & i = j\text{일 때} \end{cases}
```

정규직교 열을 가진 행렬 $`Q`$는 $`Q^T Q = I`$를 만족한다. 일반적으로 $`m > n`$.

참고: $`Q`$는 정사각일 필요가 없다. $`QQ^T \neq I`$인 것이 일반적이다.

**예제:** $`Q \in \mathbb{R}^{3 \times 2}`$ 정규직교 행렬:

```math
Q = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix}
```

```math
Q^T Q = \frac{1}{2}\begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = I_{2 \times 2}
```

```math
QQ^T = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & -1 & 0 \\ 1 & 1 & 0 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{pmatrix} \neq I_{3 \times 3}
```

$`Q`$가 정사각일 때, $`Q^T Q = I \Longrightarrow Q^T = Q^{-1}`$ (역행렬의 정의에 의해).

### 4.4 직교 행렬의 예 (Examples of Orthogonal Matrices)

**예제 1: 회전 (Rotation).** $`Q`$는 평면의 모든 벡터를 각도 $`\theta`$만큼 회전시킨다.

```math
Q = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}, \quad Q^T = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}
```

```math
Q^T Q = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
```

$`\mathbf{a}_1`$과 $`\mathbf{a}_2`$는 평면 $`\mathbb{R}^2`$의 정규직교 기저이다. 표준 기저 벡터 $`\mathbf{i}, \mathbf{j}`$는 $`\theta`$만큼 회전된다. $`Q^{-1}`$는 벡터를 $`-\theta`$만큼 역회전시킨다.

```math
Q^{-1} = \begin{pmatrix} \cos(-\theta) & -\sin(-\theta) \\ \sin(-\theta) & \cos(-\theta) \end{pmatrix}
```

**예제 2: 순열 (Permutation).**

```math
\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} y \\ z \\ x \end{pmatrix}, \quad \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} y \\ x \end{pmatrix}
```

모든 열은 단위 벡터이다. 또한 직교이다. 순열 행렬의 역은 전치이다: $`Q^{-1} = Q^T`$. **모든 순열 행렬은 직교 행렬이다.**

**예제 3: 반사 (Reflection).**

단위 법선 벡터 $`\mathbf{n}`$ ($`\|\mathbf{n}\| = 1`$)이 주어지면, $`\mathbf{u}`$를 분해한다:

```math
\mathbf{u} = \mathbf{u}_{\text{normal}} + \mathbf{u}_{\text{tangential}} = \mathbf{u}_\perp + \mathbf{u}_\parallel
```

```math
\mathbf{u}_\perp = (\mathbf{u} \cdot \mathbf{n})\mathbf{n} = \mathbf{n}(\mathbf{n}^T\mathbf{u}) = (\mathbf{n}\mathbf{n}^T)\mathbf{u}
```

```math
\mathbf{u}_\parallel = \mathbf{u} - \mathbf{u}_\perp = I\mathbf{u} - (\mathbf{n}\mathbf{n}^T)\mathbf{u} = (I - \mathbf{n}\mathbf{n}^T)\mathbf{u}
```

반사 $`\mathbf{v}`$:

```math
\mathbf{v} = \mathbf{u}_\parallel - \mathbf{u}_\perp = (I - \mathbf{n}\mathbf{n}^T)\mathbf{u} - (\mathbf{n}\mathbf{n}^T)\mathbf{u} = (I - 2\mathbf{n}\mathbf{n}^T)\mathbf{u}
```

```math
Q = I - 2\mathbf{n}\mathbf{n}^T
```

**$`Q`$가 직교임을 검증:**

```math
Q^T Q = (I - 2\mathbf{n}\mathbf{n}^T)^T(I - 2\mathbf{n}\mathbf{n}^T) = (I - 2\mathbf{n}\mathbf{n}^T)(I - 2\mathbf{n}\mathbf{n}^T)
```

```math
= I - 4\mathbf{n}\mathbf{n}^T + 4(\mathbf{n}\mathbf{n}^T)(\mathbf{n}\mathbf{n}^T) = I - 4\mathbf{n}\mathbf{n}^T + 4\mathbf{n}\mathbf{n}^T = I
```

($`\mathbf{n}^T\mathbf{n} = \|\mathbf{n}\|^2 = 1`$이므로.)

**핵심 성질:** 회전, 순열, 반사 행렬은 모든 벡터의 **길이와 각도를 보존** 한다:

```math
\|Q\mathbf{x}\|^2 = (Q\mathbf{x})^T(Q\mathbf{x}) = \mathbf{x}^T Q^T Q\mathbf{x} = \mathbf{x}^T I\mathbf{x} = \|\mathbf{x}\|^2
```

$`Q`$가 정규직교 열을 가지면 ($`Q^T Q = I`$), 길이를 변하지 않게 둔다.

### 4.5 정규직교 기저를 이용한 사영 (Projections Using Orthonormal Bases)

$`Q`$가 $`A`$를 대체한다.

$`C(A)`$로의 $`\mathbf{b}`$의 사영은 $`\mathbf{p} = A(A^T A)^{-1}A^T\mathbf{b}`$임을 상기하자.

$`A = Q`$이고 $`A^T A = Q^T Q = I`$이면:

```math
\mathbf{p} = QQ^T\mathbf{b}
```

```math
\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}
```

```math
\mathbf{p} = Q\hat{\boldsymbol{\alpha}} = \mathbf{q}_1 d_1 + \mathbf{q}_2 d_2 + \cdots + \mathbf{q}_n d_n
```

```math
= \mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \cdots + \mathbf{q}_n(\mathbf{q}_n^T\mathbf{b})
```

**비고:** $`Q \in \mathbb{R}^{n \times n}`$ (정사각)이면, $`Q^T = Q^{-1}`$이므로, $`\hat{\boldsymbol{\alpha}} = Q^{-1}\mathbf{b}`$이고 $`\mathbf{p} = QQ^{-1}\mathbf{b} = \mathbf{b}`$, $`P = I`$.

```math
\mathbf{b} = \mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \cdots + \mathbf{q}_n(\mathbf{q}_n^T\mathbf{b}) = QQ^T\mathbf{b}
```

$`QQ^T = I`$는 **푸리에 급수** (Fourier Series)의 기초이며, $`\mathbf{b}`$를 수직 조각으로 분해한다. 그 다음 기저 벡터의 선형결합이 $`\mathbf{b}`$를 다시 합친다.

### 4.6 예제 4: 정사각 직교 행렬 (Square Orthogonal Matrix)

```math
Q = \frac{1}{3}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix}, \quad Q^T = Q
```

```math
Q^T Q = \frac{1}{9}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix}\begin{pmatrix} -1 & 2 & 2 \\ 2 & -1 & 2 \\ 2 & 2 & -1 \end{pmatrix} = \frac{1}{9}\begin{pmatrix} 9 & 0 & 0 \\ 0 & 9 & 0 \\ 0 & 0 & 9 \end{pmatrix} = I = QQ^T
```

$`\mathbf{b} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$로 하면:

```math
QQ^T\mathbf{b} = I\mathbf{b} = \mathbf{b}
```

```math
\mathbf{q}_1(\mathbf{q}_1^T\mathbf{b}) + \mathbf{q}_2(\mathbf{q}_2^T\mathbf{b}) + \mathbf{q}_3(\mathbf{q}_3^T\mathbf{b})
```

```math
= \frac{1}{3}\begin{pmatrix} -1 \\ 2 \\ 2 \end{pmatrix} \cdot 2 + \frac{1}{3}\begin{pmatrix} 2 \\ -1 \\ 2 \end{pmatrix} \cdot 2 + \frac{1}{3}\begin{pmatrix} 2 \\ 2 \\ -1 \end{pmatrix} \cdot (-1) = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
```

### 4.7 그람-슈미트 과정 (The Gram-Schmidt Process)

벡터 $`\mathbf{a}, \mathbf{b}, \mathbf{c}`$가 주어지면, 세 직교 벡터 $`\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3`$를 만들고, 이로부터 세 정규직교 벡터를 생성한다:

```math
\mathbf{q}_1 = \frac{\mathbf{r}_1}{\|\mathbf{r}_1\|}, \quad \mathbf{q}_2 = \frac{\mathbf{r}_2}{\|\mathbf{r}_2\|}, \quad \mathbf{q}_3 = \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|}
```

**i)** $`\mathbf{r}_1 = \mathbf{a}`$이므로, $`\mathbf{q}_1 = \mathbf{r}_1/\|\mathbf{r}_1\|`$.

모든 벡터 $`\mathbf{a}`$, $`\mathbf{r}_1`$, $`\mathbf{q}_1`$은 한 직선 위에 있다.

**ii)** $`\mathbf{r}_2`$는 $`\mathbf{r}_1`$에 수직이어야 한다.

$`\mathbf{b}`$를 $`\mathbf{q}_1`$에 대해 $`\mathbf{b}_\parallel`$과 $`\mathbf{b}_\perp`$로 분리:

```math
\mathbf{b}_\parallel = (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1 = \left(\mathbf{b} \cdot \frac{\mathbf{r}_1}{\|\mathbf{r}_1\|}\right)\frac{\mathbf{r}_1}{\|\mathbf{r}_1\|} = \frac{1}{\|\mathbf{r}_1\|^2}(\mathbf{r}_1^T\mathbf{b})\,\mathbf{r}_1 = \left(\frac{\mathbf{r}_1^T\mathbf{b}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1
```

```math
\mathbf{b}_\perp = \mathbf{b} - \mathbf{b}_\parallel = \mathbf{b} - \left(\frac{\mathbf{r}_1^T\mathbf{b}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1 = \mathbf{r}_2
```

```math
\mathbf{q}_2 = \frac{\mathbf{r}_2}{\|\mathbf{r}_2\|}
```

**iii)** $`\mathbf{r}_3`$는 $`\mathbf{r}_1, \mathbf{r}_2`$에 수직이어야 한다.

모든 새 벡터에서 이미 설정된 방향으로의 사영을 뺀다:

```math
\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2
```

```math
= \mathbf{c} - \left(\frac{\mathbf{r}_1^T\mathbf{c}}{\mathbf{r}_1^T\mathbf{r}_1}\right)\mathbf{r}_1 - \left(\frac{\mathbf{r}_2^T\mathbf{c}}{\mathbf{r}_2^T\mathbf{r}_2}\right)\mathbf{r}_2
```

```math
\mathbf{q}_3 = \frac{\mathbf{r}_3}{\|\mathbf{r}_3\|}
```

모든 벡터 $`\mathbf{a}, \mathbf{b}, \mathbf{c}, \mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3, \mathbf{q}_1, \mathbf{q}_2, \mathbf{q}_3`$는 하나의 부분공간 ($`\mathbb{R}^3`$)에 있다.

### 4.8 그람-슈미트 예제 (Gram-Schmidt Example)

```math
\mathbf{a} = \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix}, \quad \mathbf{c} = \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix}
```

**i)** $`\mathbf{r}_1 = \mathbf{a} = \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}`$, $`\|\mathbf{r}_1\| = \sqrt{2}`$:

```math
\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}
```

**ii)** $`\mathbf{r}_2 = \mathbf{b} - (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1`$:

```math
= \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix} - \frac{1}{\sqrt{2}} \cdot 2 \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} = \begin{pmatrix} 2 \\ 0 \\ -2 \end{pmatrix} - \begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}
```

```math
\mathbf{q}_2 = \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}
```

**iii)** $`\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2`$:

```math
= \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix} - \frac{6}{\sqrt{2}} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} - \frac{-6}{\sqrt{6}} \cdot \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix}
```

```math
= \begin{pmatrix} 3 \\ -3 \\ 3 \end{pmatrix} - 3\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix} + \begin{pmatrix} 1 \\ 1 \\ -2 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}
```

```math
\mathbf{q}_3 = \frac{1}{\sqrt{3}}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}
```

### 4.9 분해 $`A = QR`$ (Factorization $`A = QR`$)

```math
(\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = (\mathbf{q}_1 \;\; \mathbf{q}_2 \;\; \mathbf{q}_3)\,R
```

벡터 $`\mathbf{a}, \mathbf{b}, \mathbf{c}`$는 $`\mathbf{q}_1, \mathbf{q}_2, \mathbf{q}_3`$의 결합이다:

```math
\mathbf{a} = (\mathbf{q}_1 \cdot \mathbf{a})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{a})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{a})\mathbf{q}_3
```

```math
\mathbf{b} = (\mathbf{q}_1 \cdot \mathbf{b})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{b})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{b})\mathbf{q}_3
```

```math
\mathbf{c} = (\mathbf{q}_1 \cdot \mathbf{c})\mathbf{q}_1 + (\mathbf{q}_2 \cdot \mathbf{c})\mathbf{q}_2 + (\mathbf{q}_3 \cdot \mathbf{c})\mathbf{q}_3
```

참고: 나중의 $`\mathbf{q}`$들은 이전의 $`\mathbf{a}`$들에 직교하므로, $`R`$의 대각선 아래 많은 항이 0이다:

```math
(\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = (\mathbf{q}_1 \;\; \mathbf{q}_2 \;\; \mathbf{q}_3)\begin{pmatrix} \mathbf{q}_1^T\mathbf{a} & \mathbf{q}_1^T\mathbf{b} & \mathbf{q}_1^T\mathbf{c} \\ 0 & \mathbf{q}_2^T\mathbf{b} & \mathbf{q}_2^T\mathbf{c} \\ 0 & 0 & \mathbf{q}_3^T\mathbf{c} \end{pmatrix}
```

```math
A = QR
```

$`Q^T A = Q^T Q R = IR = R`$.

LI인 $`\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_n`$으로부터 그람-슈미트가 정규직교 벡터 $`\mathbf{q}_1, \mathbf{q}_2, \ldots, \mathbf{q}_n`$을 만든다. 이들로 구성된 행렬은 $`A = QR`$을 만족한다. $`R = Q^T A`$는 나중의 $`\mathbf{q}`$들이 이전의 $`\mathbf{a}`$들에 직교하므로 **상삼각 행렬** 이다.

**예제로부터:**

```math
A = (\mathbf{a} \;\; \mathbf{b} \;\; \mathbf{c}) = \begin{pmatrix} 1 & 2 & 3 \\ -1 & 0 & -3 \\ 0 & -2 & 3 \end{pmatrix}
```

```math
= \frac{1}{\sqrt{6}}\begin{pmatrix} \sqrt{3} & 1 & \sqrt{2} \\ -\sqrt{3} & 1 & \sqrt{2} \\ 0 & -2 & \sqrt{2} \end{pmatrix} \begin{pmatrix} \sqrt{2} & \sqrt{2} & \sqrt{18} \\ 0 & \sqrt{6} & \sqrt{6} \\ 0 & 0 & \sqrt{3} \end{pmatrix}
```

```math
Q \qquad \qquad R
```

$`\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3`$의 **길이** 는 $`R`$의 대각선이며, 양수이다.

독립 열을 가진 임의의 $`m \times n`$ 행렬 $`A`$는 $`A = QR`$로 분해할 수 있다. $`Q \in \mathbb{R}^{m \times n}`$은 직교 열을 갖고, $`R \in \mathbb{R}^{n \times n}`$은 양의 대각선을 가진 상삼각 행렬이다.

### 4.10 응용: QR을 통한 최소제곱 (Application: Least Squares via QR)

```math
A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}
```

$`A = QR`$를 대입: $`A^T A = (QR)^T QR = R^T Q^T QR = R^T R`$.

```math
R^T R\hat{\boldsymbol{\alpha}} = R^T Q^T\mathbf{b}
```

$`R`$이 가역이므로:

```math
R\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b} \quad \text{(매우 빠르다)}
```

```math
\hat{\boldsymbol{\alpha}} = R^{-1}Q^T\mathbf{b}
```

---

<br>

## 5. 행렬의 유사역행렬 (4.5)

### 5.1 양측 역행렬과 편측 역행렬 (Two-Sided and One-Sided Inverses)

**(1)** **양측 역행렬 (two-sided inverse):**

```math
A^{-1}A = AA^{-1} = I
```

**(2)** **편측 역행렬 (one-sided inverse):**

```math
A^+ A = I \quad (A^+ \text{는 } A\text{의 좌역행렬})
```

```math
AA^+ = I \quad (A^+ \text{는 } A\text{의 우역행렬})
```

모든 행렬 $`A`$는 유사역행렬(pseudoinverse) $`A^+`$를 갖는다.

### 5.2 항등-유사 행렬의 세 경우 (Three Cases with Identity-like Matrices)

**(1)** $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}_{2 \times 2}`$, $`m = n = 2`$, $`\text{rank}(I) = 2 = r`$.

$`\dim C(I) = 2`$, $`\dim \mathcal{N}(I^T) = 0`$, $`\dim C(I^T) = 2`$, $`\dim \mathcal{N}(I) = 0`$.

네 부분공간 모두 전체이거나 자명하다.

**(2)** $`I_L = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}_{2 \times 3}`$, $`m = 2`$, $`n = 3`$, $`\text{rank}(I_L) = 2 = r = m`$.

모든 행이 LI이지만, $`r = m < n`$. 영공간은 비자명 원소를 갖는다.

$`\dim C(I_L) = 2`$, $`\dim \mathcal{N}(I_L^T) = 0`$, $`\dim C(I_L^T) = 2`$, $`\dim \mathcal{N}(I_L) = 1`$.

$`I_L I_R = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I_{2 \times 2}`$

$`I_L`$ = $`I_R`$의 좌역행렬, $`I_R`$ = $`I_L`$의 우역행렬.

**(3)** $`I_R = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}_{3 \times 2}`$, $`m = 3`$, $`n = 2`$, $`\text{rank}(I_R) = 2 = r = n`$.

모든 열이 LI이지만, $`r < m`$.

$`\dim C(I_R) = 2`$, $`\dim \mathcal{N}(I_R^T) = 1`$, $`\dim C(I_R^T) = 2`$, $`\dim \mathcal{N}(I_R) = 0`$.

### 5.3 좌역행렬과 우역행렬 (Left Inverse and Right Inverse)

| | $`A^+A = I_{n \times n}`$ (좌역행렬) | $`AA^+ = I_{m \times m}`$ (우역행렬) |
|:---|:---|:---|
| **조건** | 풀 열 랭크: $`r = n < m`$ | 풀 행 랭크: $`r = m < n`$ |
| **해석** | 미지수 수 < 방정식 수 | 방정식 수 < 미지수 수 |
| **$`A\mathbf{x} = \mathbf{b}`$의 해** | 0개 또는 1개 | 무한히 많은 해 |
| **영공간** | $`\mathcal{N}(A) = \lbrace\mathbf{0}\rbrace`$ | $`\mathcal{N}(A^T) = \lbrace\mathbf{0}\rbrace`$ |
| **$`A^T A`$ 또는 $`AA^T`$** | $`A^T A`$는 $`n \times n`$이고 가역 | $`AA^T`$는 $`m \times m`$이고 가역 |
| **공식** | $`A^+ = (A^T A)^{-1}A^T`$ | $`A^+ = A^T(AA^T)^{-1}`$ |

**좌역행렬의 경우:** $`A^+A = I`$는 최소제곱에서의 행렬을 설명한다. $`\hat{\boldsymbol{\alpha}} = A^+\mathbf{b}`$는 $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$의 해이며, $`\mathbf{b}`$는 $`C(A)`$에 없을 수 있고 $`\hat{\boldsymbol{\alpha}}`$는 $`C(A^T)`$에 있다.

**우역행렬의 경우:** $`AA^+ = I`$. $`\mathbf{x}^+ = A^+\mathbf{b}`$는 $`A\mathbf{x} = \mathbf{b}`$의 최소 길이 해이며, $`\mathbf{b}`$는 $`C(A)`$에 있고 $`\mathbf{x}^+`$는 $`C(A^T)`$에 있다.

### 5.4 일반 행렬 $`A_{m \times n}`$의 유사역행렬 $`A^+`$ (The Pseudoinverse $`A^+`$ of a General Matrix)

**단계 1:** 모든 벡터 $`\mathbf{b} \in \mathbb{R}^m`$은 두 수직 부분 $`\mathbf{p}`$와 $`\mathbf{z}`$를 가지며, $`\mathbf{b} = \mathbf{p} + \mathbf{z}`$이다.

**단계 2:** $`\mathbf{p} \in C(A)`$이고 $`\mathbf{z} \in \mathcal{N}(A^T)`$ ($`\Leftrightarrow A^T\mathbf{z} = \mathbf{0}`$, $`A^+\mathbf{z} = \mathbf{0}`$).

**단계 3:** $`C(A^T) \ni \mathbf{x}^+ \longrightarrow A\mathbf{x}^+ = \mathbf{p}`$. 이 부분을 역전: $`\mathbf{x}^+ = A^+\mathbf{p}`$.

**단계 4:** $`A^+\mathbf{b} = A^+(\mathbf{p} + \mathbf{z}) = A^+\mathbf{p} + A^+\mathbf{z} = \mathbf{x}^+`$.

### 5.5 $`A^+`$의 큰 그림 (Big Picture for $`A^+`$)

- $`A^+`$는 $`A^T`$와 동일한 네 부분공간을 공유한다.
- $`A^+`$는 가능할 때 $`A`$를 역전한다: 열 공간에서 행 공간으로.

```math
C(A^T) \xrightarrow{A} C(A) \xrightarrow{A^+} C(A^T)
```

- $`A^+A`$는 $`\mathbf{x}^+ \in C(A^T)`$를 같은 $`\mathbf{x}^+`$로 되돌린다.
- $`A^+A`$는 $`A`$의 행 공간 $`C(A^T)`$로의 $`n \times n`$ 사영 행렬 $`P_{\text{row}}`$이다.

```math
P_{\text{row}} = A^+A = (A^+A)^2 = (A^+A)^T
```

- $`AA^+`$는 $`A`$의 열 공간 $`C(A)`$로의 $`m \times m`$ 사영 행렬 $`P_{\text{col}}`$이다.

```math
P_{\text{col}} = AA^+ = (AA^+)^2 = (AA^+)^T
```

### 5.6 예제 1: $`A^+`$ 구하기 (Finding $`A^+`$)

```math
A = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ b_2 \end{pmatrix} = \mathbf{p} + \mathbf{z}
```

$`r = 1`$.

**i)** $`\mathbf{p} \in C(A)`$: $`A\mathbf{x} = \mathbf{p}`$: $`\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ 0 \end{pmatrix}`$

$`2x_1 = b_1`$, $`x_1 = b_1/2`$.

```math
\mathbf{x} = \begin{pmatrix} b_1/2 \\ x_2 \end{pmatrix} = \underbrace{\begin{pmatrix} b_1/2 \\ 0 \end{pmatrix}}_{\mathbf{x}^+} + \underbrace{\begin{pmatrix} 0 \\ x_2 \end{pmatrix}}_{\mathbf{x}_n}
```

**ii)** $`\mathbf{z} \in \mathcal{N}(A^T)`$: $`A^T\mathbf{z} = \mathbf{0}`$: $`\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} z_1 \\ z_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$이므로 $`z_1 = 0`$, $`\mathbf{z} = \begin{pmatrix} 0 \\ z_2 \end{pmatrix}`$.

**iii)** $`A^+\mathbf{b} = A^+(\mathbf{p} + \mathbf{z}) = A^+\mathbf{p} + A^+\mathbf{z} = \mathbf{x}^+ = \begin{pmatrix} b_1/2 \\ 0 \end{pmatrix}`$

```math
= \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
```

```math
A^+ = \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}
```

**iv)** 검증:

```math
AA^+ = \begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = P_{\text{col}}
```

```math
A^+A = \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = P_{\text{row}}
```

### 5.7 대각 행렬의 유사역행렬 (Pseudoinverse of a Diagonal Matrix)

```math
D = \begin{pmatrix} 2 & & \\ & 3 & \\ & & 0 \end{pmatrix}, \quad D^+ = \begin{pmatrix} 1/2 & & \\ & 1/3 & \\ & & 0 \end{pmatrix}
```

0이 아닌 대각 성분을 역수로 바꾸고; 0은 0으로 둔다.

### 5.8 분해를 통한 유사역행렬 (Pseudoinverse via Factorization)

$`U, V`$를 가역 행렬이라 하자. $`A = UDV^T`$이면:

```math
A^+ = (UDV^T)^+ = V^{-T}D^+U^{-1}
```

**예제 (풀 열 랭크):** $`A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`\text{rank}(A) = 1 = r = n < m`$. 풀 열 랭크, 좌역행렬이 존재한다.

$`A^+A = 1`$. $`A^+ = (A^T A)^{-1}A^T = ((1\;1)\begin{pmatrix} 1 \\ 1 \end{pmatrix})^{-1}(1\;1) = \frac{1}{2}(1\;1)`$.

**예제 (풀 행 랭크):** $`A = (1\;1)`$, $`\text{rank}(A) = 1 = m < n`$. 풀 행 랭크, 우역행렬이 존재한다.

$`AA^+ = 1`$. $`A^+ = A^T(AA^T)^{-1} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}((1\;1)\begin{pmatrix} 1 \\ 1 \end{pmatrix})^{-1} = \frac{1}{2}\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$.

### 5.9 $`A`$의 중요한 작용: 행 공간에서 열 공간으로 (Important Action of $`A`$: Row Space to Column Space)

```math
C(A^T) \xrightarrow{A} C(A)
```

```math
C(A) \xrightarrow{A^+} C(A^T)
```

**i)** $`\mathbf{x}_1, \mathbf{x}_2 \in C(A^T)`$, $`\mathbf{x}_1 \neq \mathbf{x}_2 \Rightarrow A\mathbf{x}_1 \neq A\mathbf{x}_2`$.

**증명:** $`A\mathbf{x}_1 = \mathbf{b}`$이고 $`A\mathbf{x}_2 = \mathbf{b}`$라 가정하자. 빼면: $`A(\mathbf{x}_1 - \mathbf{x}_2) = \mathbf{0}`$. 영공간의 정의에 의해, $`\mathbf{x}_1 - \mathbf{x}_2 \in \mathcal{N}(A)`$. 그런데 $`\mathbf{x}_1 - \mathbf{x}_2 \in C(A^T)`$. $`\mathcal{N}(A)`$가 $`C(A^T)`$에 직교하므로, $`\mathbf{x}_1 - \mathbf{x}_2 = \mathbf{0}`$. 이는 $`\mathbf{x}_1 \neq \mathbf{x}_2`$에 모순이다. 따라서 $`A\mathbf{x}_1 \neq A\mathbf{x}_2`$. $`\square`$

**ii)** $`\forall\,\mathbf{b} \in C(A)`$, $`\exists!\,\mathbf{x}^+ \in C(A^T)`$ s.t. $`A^+\mathbf{b} = \mathbf{x}^+`$.

### 5.10 $`A = CR`$의 유사역행렬 $`A^+ = R^+C^+`$ (Pseudoinverse $`A^+ = R^+C^+`$ of $`A = CR`$)

$`A = CR`$이 주어지고, $`C`$는 $`r`$개의 독립 열(풀 열 랭크), $`R`$는 $`r`$개의 독립 행(풀 행 랭크)을 가지면:

- $`C`$ 풀 열 랭크: $`\exists (C^T C)^{-1}`$, $`C`$의 좌역행렬 존재: $`C^+C = I`$, $`C^+ = (C^T C)^{-1}C^T`$.
- $`R`$ 풀 행 랭크: $`\exists (RR^T)^{-1}`$, $`R`$의 우역행렬 존재: $`RR^+ = I`$, $`R^+ = R^T(RR^T)^{-1}`$.

```math
A^+ = R^+C^+ = R^T(RR^T)^{-1}(C^T C)^{-1}C^T
```

```math
= R^T(C^T A R^T)^{-1}C^T
```

($`C^T C R R^T = C^T(CR)R^T = C^T A R^T`$이므로.)

### 5.11 예제: $`5 \times 4`$ 행렬의 유사역행렬 (Pseudoinverse of a $`5 \times 4`$ Matrix)

```math
A = \begin{pmatrix} -1 & 1 & 0 & 0 \\ -1 & 0 & 1 & 0 \\ 0 & -1 & 1 & 0 \\ -1 & 0 & 0 & 1 \\ 0 & -1 & 0 & 1 \end{pmatrix}_{5 \times 4}
```

행 축소 단계: $`R_2 - R_1`$, $`R_4 - R_1`$, 그 다음 $`R_3 - R_2`$, $`R_4 - R_2`$, $`R_5 - R_2`$, $`R_3`$와 $`R_5`$ 교환, $`R_4 - R_3`$, 그 다음 $`R_1 + R_2`$, $`R_2 + R_3`$, $`\times(-1)`$:

```math
R_0 = \begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}
```

열 1, 2, 3이 LI이다. $`\text{rank}(A) = 3`$.

```math
R = \begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{pmatrix}
```

```math
A = CR = \begin{pmatrix} -1 & 1 & 0 \\ -1 & 0 & 1 \\ 0 & -1 & 1 \\ -1 & 0 & 0 \\ 0 & -1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & -1 \end{pmatrix}
```

$`A^+ = (CR)^+ = R^+C^+ = R^T(C^T AR^T)^{-1}C^T`$.

$`C^T AR^T = \begin{pmatrix} 4 & 0 & 0 \\ 0 & 4 & 0 \\ -1 & -1 & 2 \end{pmatrix}`$임이 밝혀진다.

```math
A^+ = \frac{1}{8}\begin{pmatrix} -2 & -2 & 0 & -2 & 0 \\ 2 & 0 & -2 & 0 & -2 \\ 0 & 3 & 3 & -1 & -1 \\ 0 & -1 & -1 & 3 & 3 \end{pmatrix}
```

$`A = CR`$를 통한 사상:

```math
C(A^T) \ni \mathbf{x} \xrightarrow{R} R\mathbf{x} \xrightarrow{C} CR\mathbf{x} = A\mathbf{x} \in C(A)
```

$`R`$은 풀 행 랭크: $`\mathcal{N}(R^T) = \lbrace\mathbf{0}\rbrace`$. $`C`$는 풀 열 랭크: $`\mathcal{N}(C) = \lbrace\mathbf{0}\rbrace`$.

---

<br>

## 요약

| 개념 | 핵심 내용 |
|:--------|:---------|
| 직교 벡터 (orthogonal vectors) | $`\mathbf{v}^T\mathbf{w} = 0`$이면 $`\|\mathbf{v}\|^2 + \|\mathbf{w}\|^2 = \|\mathbf{v} + \mathbf{w}\|^2`$ |
| 기본 부분공간의 직교성 | $`C(A^T) \perp \mathcal{N}(A)`$ ($`\mathbb{R}^n`$); $`C(A) \perp \mathcal{N}(A^T)`$ ($`\mathbb{R}^m`$) |
| 직교 여공간 (orthogonal complement) | $`\dim V + \dim V^\perp = n`$; 모든 $`\mathbf{x} = \mathbf{x}_r + \mathbf{x}_n`$ |
| 직선으로의 사영 | $`\mathbf{p} = \frac{\mathbf{a}^T\mathbf{b}}{\mathbf{a}^T\mathbf{a}}\mathbf{a}`$; 사영 행렬 $`P = \frac{\mathbf{a}\mathbf{a}^T}{\mathbf{a}^T\mathbf{a}}`$ |
| 부분공간으로의 사영 | $`\mathbf{p} = A(A^T A)^{-1}A^T\mathbf{b}`$; $`P = A(A^T A)^{-1}A^T`$; $`P^2 = P = P^T`$ |
| $`A^T A`$ 가역성 | $`A^T A`$ 가역 $`\iff`$ $`A`$가 LI 열 $`\iff`$ $`\mathcal{N}(A) = \lbrace\mathbf{0}\rbrace`$ |
| 최소제곱 해 | $`\hat{\boldsymbol{\alpha}} = (A^T A)^{-1}A^T\mathbf{b}`$는 $`E = \|A\mathbf{x} - \mathbf{b}\|^2`$를 최소화 |
| 정규방정식 (normal equations) | $`A^T A\hat{\boldsymbol{\alpha}} = A^T\mathbf{b}`$ (기하, 대수, 미적분 유도) |
| 직선 적합 (line fitting) | $`A = \begin{pmatrix} x_1 & 1 \\ \vdots & \vdots \\ x_m & 1 \end{pmatrix}`$; $`A^T A = \begin{pmatrix} \sum x_i^2 & \sum x_i \\ \sum x_i & m \end{pmatrix}`$ |
| 포물선 적합 (parabola fitting) | $`A`$의 열은 $`(1, 1, \ldots)`$, $`(t_1, t_2, \ldots)`$, $`(t_1^2, t_2^2, \ldots)`$ |
| 정규직교 벡터 (orthonormal vectors) | $`\mathbf{q}_i^T\mathbf{q}_j = \delta_{ij}`$; $`Q^T Q = I`$ |
| 직교 행렬 (정사각 $`Q`$) | $`Q^T Q = QQ^T = I`$; $`Q^T = Q^{-1}`$; 길이와 각도 보존 |
| 직교 행렬의 예 | 회전, 순열, 반사 ($`Q = I - 2\mathbf{n}\mathbf{n}^T`$) |
| 정규직교 기저로의 사영 | $`\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}`$; $`\mathbf{p} = QQ^T\mathbf{b}`$; $`Q`$ 정사각이면 $`P = I`$ |
| 그람-슈미트 과정 | $`\mathbf{r}_1 = \mathbf{a}`$; $`\mathbf{r}_2 = \mathbf{b} - (\mathbf{b} \cdot \mathbf{q}_1)\mathbf{q}_1`$; $`\mathbf{r}_3 = \mathbf{c} - (\mathbf{c} \cdot \mathbf{q}_1)\mathbf{q}_1 - (\mathbf{c} \cdot \mathbf{q}_2)\mathbf{q}_2`$; 각각 정규화 |
| QR 분해 (QR factorization) | $`A = QR`$; $`Q`$ 직교 열, $`R`$ 양의 대각선 상삼각 행렬 |
| QR을 통한 최소제곱 | $`R\hat{\boldsymbol{\alpha}} = Q^T\mathbf{b}`$ (정규방정식보다 훨씬 빠름) |
| 유사역행렬 $`A^+`$ (pseudoinverse) | $`C(A) \to C(A^T)`$로 $`A`$를 역전; $`\mathbf{z} \in \mathcal{N}(A^T)`$이면 $`A^+\mathbf{z} = \mathbf{0}`$ |
| 좌역행렬 ($`r = n < m`$) | $`A^+ = (A^T A)^{-1}A^T`$; $`A^+A = I`$; 최소제곱 해 |
| 우역행렬 ($`r = m < n`$) | $`A^+ = A^T(AA^T)^{-1}`$; $`AA^+ = I`$; 최소 길이 해 |
| $`A^+A`$와 $`AA^+`$ | $`A^+A = P_{\text{row}}`$ (행 공간으로의 사영); $`AA^+ = P_{\text{col}}`$ (열 공간으로의 사영) |
| $`A = CR`$을 통한 유사역행렬 | $`A^+ = R^+C^+ = R^T(C^T AR^T)^{-1}C^T`$ |

---
