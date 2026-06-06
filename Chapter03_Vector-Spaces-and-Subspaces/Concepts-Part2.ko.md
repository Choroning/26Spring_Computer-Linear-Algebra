# 제3장 강의 — 벡터 공간과 부분공간 — Part 2/2

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.**
>
> [Part 1](Concepts-Part1.ko.md) · **Part 2**

---

<br>

## 목차

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

## 5. 독립, 기저, 차원 (3.4)

### 5.1 독립 벡터

**(1) 독립 벡터 (Independent Vectors):** 유일한 영 결합

```math
c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \cdots + c_n \mathbf{v}_n = \mathbf{0}
```

은 모든 $`c_1 = c_2 = \cdots = c_n = 0`$을 가진다.

이는 다음을 의미한다: 적어도 하나의 스칼라가 0이 아닌 경우, 예를 들어 $`c_1 \neq 0`$이면:

```math
\mathbf{v}_1 = -\frac{c_2}{c_1}\mathbf{v}_2 - \frac{c_3}{c_1}\mathbf{v}_3 - \cdots - \frac{c_n}{c_1}\mathbf{v}_n
```

$`\mathbf{v}_1`$은 다른 벡터들의 선형 결합이다.

**(2) 벡터 $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$이 공간 $`\mathbb{S}`$를 생성한다 (Span): $`\mathbb{S}`$ = $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$의 모든 결합일 때.**

예: $`\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$ $`\implies`$ $`\begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}`$. $`\hat{i}, \hat{j}`$가 공간 $`\mathbb{R}^2`$를 생성한다.

**(3)** 벡터 $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$이 $`\mathbb{S}`$의 **기저 (Basis)** 가 되려면: **선형 독립 (Linearly Independent)** 이면서 $`\mathbb{S}`$를 **생성 (Span)** 해야 한다.

공간의 모든 벡터는 기저 벡터들의 **유일한** 결합이다.

예: $`\mathbb{R}^2 \ni \begin{pmatrix} x \\ y \end{pmatrix} = x\hat{i} + y\hat{j}`$

**(4)** 벡터 공간 $`\mathbb{S}`$의 **차원 (Dimension)** 은 $`\mathbb{S}`$의 모든 기저에 있는 벡터의 수 $`n`$이다.

$`A \in \mathbb{R}^{m \times n}`$을 고려하자. $`n`$개의 열 중 $`r`$개가 독립이고, 나머지 $`n - r`$개의 열이 종속이다. $`C(A)`$의 차원은 $`r`$이며, 이는 $`A`$의 계수 (Rank)이다.

이 절의 네 가지 핵심 개념:
1. 독립 벡터 (Independent Vectors)
2. 공간의 생성 (Spanning a Space)
3. 공간의 기저 (Basis for a Space)
4. 공간의 차원 (Dimension of a Space)

### 5.2 영공간을 통한 선형 독립

**정의.** $`A`$의 열들이 **선형 독립 (Linearly Independent)** 이려면 $`A\mathbf{x} = \mathbf{0}`$의 유일한 해가 $`\mathbf{x} = \mathbf{0}`$이어야 한다.

```math
\begin{pmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \cdots & \mathbf{a}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 + \cdots + x_n \mathbf{a}_n = \mathbf{0}
```

$`x_1 = x_2 = \cdots = x_n = 0`$일 때만 성립.

**기하학적 해석:**

예: $`\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3 \in \mathbb{R}^3`$이 같은 평면에 있지 않으면 $`\implies`$ 이 벡터들은 독립이다. $`\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3`$의 어떤 결합도 $`\mathbf{0}`$이 되지 않는다 (자명한 경우 제외).

예: $`\mathbf{w}_1, \mathbf{w}_2, \mathbf{w}_3`$이 $`\mathbb{R}^3`$에서 같은 평면에 있다. 이들은 종속이다. 예를 들어, $`\mathbf{w}_3 = \mathbf{w}_1 + \mathbf{w}_2`$ $`\iff`$ $`1 \cdot \mathbf{w}_1 + 1 \cdot \mathbf{w}_2 - 1 \cdot \mathbf{w}_3 = \mathbf{0}`$. 결합이 $`\mathbf{0}`$을 주지만, 0이 아닌 계수가 있으므로 종속이다.

**종속 여부의 빠른 판별:**

- Q. $`\begin{pmatrix} 1 \\ 0 \end{pmatrix}, \begin{pmatrix} 1 \\ 10^{-5} \end{pmatrix}`$는 종속인가? **아니다.**
- Q. $`\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} -1 \\ -1 \end{pmatrix}`$는 종속인가? **그렇다.**
- Q. $`\begin{pmatrix} 1 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$는 종속인가? **그렇다.**
- Q. $`\mathbb{R}^2`$에서 임의의 세 벡터는 종속이다. **참.**

**예 1.** $`A = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix}`$은 종속 열을 가지는가?

$`A`$의 영공간 $`N(A)`$를 확인하자.

```math
A\mathbf{x} = \begin{pmatrix} 1 & 0 & 3 \\ 2 & 1 & 5 \\ 1 & 0 & 3 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`\xrightarrow{R_2 - 2R_1, \; R_3 - R_1}`$ $`\begin{pmatrix} 1 & 0 & 3 \\ 0 & 1 & -1 \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}`$, $`\text{rank}(A) = 2 < 3`$.

$`x_3 = 1`$로 취하면: $`x_1 + 3 = 0 \Rightarrow x_1 = -3`$, $`x_2 - 1 = 0 \Rightarrow x_2 = 1`$.

0이 아닌 계수가 영벡터를 만든다:

```math
-3\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} + 1\begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix} + 1\begin{pmatrix} 3 \\ 5 \\ 3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`A`$의 열들이 독립이면, $`\text{rank}(A) = r = n`$이고 $`A`$의 영공간은 영벡터만 포함한다: $`N(A) = \lbrace\mathbf{0}\rbrace`$.

**$`\mathbb{R}^m`$에서 $`n > m`$이면 $`n`$개의 벡터는 반드시 선형 종속이다.**

예: 5개의 성분을 가진 7개의 열이 있다면 ($`5 \times 7`$ 행렬). 7개의 열벡터는 $`\mathbb{R}^5`$에서 온다. 5개의 행에 5개 초과의 피벗이 있을 수 없다. $`A\mathbf{x} = \mathbf{0}`$은 최소 $`2 \;(= 7 - 5)`$개의 자유 변수를 가진다. 즉, 0이 아닌 해를 가진다.

**참고:** $`n \leq m`$이면 열들이 종속일 수도 독립일 수도 있다. 소거가 $`r`$개의 피벗을 드러낸다.

### 5.3 부분공간을 생성하는 벡터

$`A`$를 행렬이라 하자. $`C(A)`$는 $`A\mathbf{x}`$의 모든 결합으로 이루어진 열공간이다:

```math
\begin{pmatrix} \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = x_1 \mathbf{v}_1 + x_2 \mathbf{v}_2 + \cdots + x_n \mathbf{v}_n
```

**예:** $`A = \begin{pmatrix} 1 & 4 \\ 2 & 7 \\ 3 & 5 \end{pmatrix}`$, $`A^T = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 7 & 5 \end{pmatrix}`$

$`C(A)`$는 $`\mathbb{R}^3`$에서의 평면이다. $`C(A^T)`$, 즉 $`A`$의 행공간은 $`\mathbb{R}^2`$이다.

$`A \in \mathbb{R}^{m \times n}`$에 대해: 행벡터는 $`\mathbb{R}^n`$에, 열벡터는 $`\mathbb{R}^m`$에 속한다.

### 5.4 벡터 공간의 기저

**정의.** 벡터 공간의 **기저 (Basis)** 는 두 가지 성질을 가진 벡터 수열이다:
> i) 기저 벡터들이 **선형 독립** 이다.
> ii) 그들이 공간을 **생성** 한다.

**예:** $`\begin{pmatrix} x \\ y \end{pmatrix} = x\begin{pmatrix} 1 \\ 0 \end{pmatrix} + y\begin{pmatrix} 0 \\ 1 \end{pmatrix}`$이고, $`\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$.

i) $`\hat{i}`$와 $`\hat{j}`$는 LI이다.
ii) $`\hat{i}, \hat{j}`$가 $`\mathbb{R}^2`$를 생성한다.

결합 $`\mathbf{x} = x\hat{i} + y\hat{j}`$는 $`\hat{i}, \hat{j}`$가 LI이므로 유일하다.

**$`\mathbf{v}`$를 기저 벡터들의 결합으로 쓰는 방법은 하나뿐이다.**

**증명.** $`\mathbf{v} = a_1 \mathbf{v}_1 + \cdots + a_n \mathbf{v}_n`$이고 $`\mathbf{v} = b_1 \mathbf{v}_1 + \cdots + b_n \mathbf{v}_n`$라 하자.

빼면 영벡터를 얻는다:

```math
\mathbf{0} = (a_1 - b_1)\mathbf{v}_1 + \cdots + (a_n - b_n)\mathbf{v}_n
```

$`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$의 LI에 의해: $`a_1 - b_1 = 0`$, $`a_2 - b_2 = 0`$, $`\ldots`$, $`a_n - b_n = 0`$.

따라서 $`a_i = b_i`$ ($`i = 1, 2, \ldots, n`$). $`\square`$

---

**예 3.** $`I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}`$의 열들이 $`\mathbb{R}^2`$의 "표준 기저 (Standard Basis)"를 만든다.

$`\hat{i} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$, $`\hat{j} = \begin{pmatrix} 0 \\ 1 \end{pmatrix}`$라 하자.

i) $`\hat{i}, \hat{j}`$는 LI이다.
ii) $`\hat{i}, \hat{j}`$가 $`\mathbb{R}^2`$를 생성한다.

$`\therefore \hat{i}, \hat{j}`$는 $`\mathbb{R}^2`$의 기저 벡터이다.

---

**예 4.** 모든 가역 $`n \times n`$ 행렬의 열들이 $`\mathbb{R}^n`$의 기저를 제공한다.

**가역 행렬 예:**

```math
A = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = R
```

3개의 LI 열벡터. $`C(A) = \mathbb{R}^3`$. 3개의 0이 아닌 피벗. $`\text{rank}(A) = 3`$. $`N(A) = \lbrace\mathbf{0}\rbrace`$.

**기저는 유일하지 않다.**

벡터 $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`$이 $`n \times n`$ 가역 행렬의 열들일 때 $`\mathbb{R}^n`$의 기저이다. $`\mathbb{R}^n`$은 무한히 많은 서로 다른 기저를 가진다.

---

**특이 행렬 (Singular Matrix) 예:**

```math
B = \begin{pmatrix} 1 & 0 & 1 \\ 1 & 1 & 2 \\ 1 & 1 & 2 \end{pmatrix} \xrightarrow{R_2 - R_1, \; R_3 - R_1} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 1 & 2 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 0 & 0 \end{pmatrix} = R_0
```

$`C(A) \neq \mathbb{R}^3`$. $`\text{rank}(A) = 2 < 3`$.

열들이 종속이면, **피벗 열 (Pivot Column)** 만 취한다.

예: $`B`$의 1번째, 2번째 열.

- 모든 독립 벡터 집합은 기저로 **확장** 할 수 있다. (예: $`B`$의 1, 2번째 열이 $`\mathbb{R}^3`$에서 평면을 생성.)
- 모든 생성 벡터 집합은 기저로 **축소** 할 수 있다.

---

**예 5.** $`A = \begin{pmatrix} 2 & 4 \\ 3 & 6 \end{pmatrix}`$

$`\xrightarrow{R_2 - \frac{3}{2}R_1}`$ $`\begin{pmatrix} 2 & 4 \\ 0 & 0 \end{pmatrix}`$ $`\xrightarrow{R_1 / 2}`$ $`\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix} = R_0`$

$`\text{rank}(A) = 1 < 2`$. 하나의 피벗 열, 하나의 피벗 행.

---

**예 6.** $`R_0 = \begin{pmatrix} 1 & 2 & 0 & 3 \\ 0 & 0 & 1 & 4 \\ 0 & 0 & 0 & 0 \end{pmatrix}`$

1번째, 3번째 열이 피벗 열이다. $`\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}`$이 $`C(R_0)`$의 기저 벡터이다.

$`C(R_0)`$는 $`\mathbb{R}^3`$에서의 $`xy`$ 평면이다.

또한, 2번째와 3번째 열벡터도 $`C(R_0)`$의 기저이다.

---

**모든 벡터 공간의 기저는 같은 수의 벡터를 포함한다.** 모든 기저에 있는 벡터의 수가 그 공간의 "**차원 (Dimension)**"이다.

### 5.5 벡터 공간의 차원

**벡터 공간의 차원:**

$`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m`$과 $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n`$이 같은 벡터 공간의 두 기저이면, $`m = n`$이다.

**증명.** $`n > m`$이라 하자. $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m`$이 기저이므로, 각 $`\mathbf{w}_i`$ ($`i = 1, 2, \ldots, n`$)는 $`\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_m`$의 결합이어야 한다:

```math
\mathbf{w}_1 = a_{11}\mathbf{v}_1 + a_{21}\mathbf{v}_2 + \cdots + a_{m1}\mathbf{v}_m
```
```math
\mathbf{w}_2 = a_{12}\mathbf{v}_1 + a_{22}\mathbf{v}_2 + \cdots + a_{m2}\mathbf{v}_m
```
```math
\vdots
```
```math
\mathbf{w}_n = a_{1n}\mathbf{v}_1 + a_{2n}\mathbf{v}_2 + \cdots + a_{mn}\mathbf{v}_m
```

이로부터:

```math
W = (\mathbf{w}_1 \quad \mathbf{w}_2 \quad \cdots \quad \mathbf{w}_n) = (\mathbf{v}_1 \quad \mathbf{v}_2 \quad \cdots \quad \mathbf{v}_m) \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & & & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix} = VA
```

$`A`$는 $`m \times n`$ 행렬 (짧고 넓음)이다. $`n > m`$이므로, $`A\mathbf{x} = \mathbf{0}`$은 0이 아닌 해를 가진다.

$`A\mathbf{x} = \mathbf{0}`$에서: $`VA\mathbf{x} = V\mathbf{0} = \mathbf{0}`$.

즉 $`W\mathbf{x} = \mathbf{0}`$, 다시 말해 $`x_1\mathbf{w}_1 + x_2\mathbf{w}_2 + \cdots + x_n\mathbf{w}_n = \mathbf{0}`$.

$`\mathbf{x}`$가 0이 아닌 벡터이므로, $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n`$은 LI가 아니다. $`\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_n`$은 기저가 될 수 **없다**. 이는 $`\mathbf{w}_i`$ ($`i = 1, 2, \ldots, n`$)가 기저라는 가정에 모순이다. $`\square`$

---

**정의.** 공간의 **차원 (Dimension)** 은 모든 기저에 있는 벡터의 수이다.

**예:** $`\mathbf{u} = \begin{pmatrix} 1 \\ 5 \\ 2 \end{pmatrix}`$를 지나는 직선은 1차원이다.

그 직선에 수직인 평면: $`\mathbf{u} \cdot \mathbf{x} = 0 \iff x + 5y + 2z = 0`$.

행렬 $`A = \begin{pmatrix} 1 & 5 & 2 \end{pmatrix}`$의 평면 영공간:

$`(1 \quad 5 \quad 2)\begin{pmatrix} x \\ y \\ z \end{pmatrix} = 0`$

$`n = 3`$, $`\text{rank}(A) = r = 1`$, $`n - r = 2`$개의 자유 변수.

$`n - r`$개의 특수해가 영공간의 기저를 제공한다: 차원 $`n - r`$.

특수해를 구하면:
- $`y = 1, z = 0`$으로 취하면 $`\implies x = -5`$
- $`y = 0, z = 1`$로 취하면 $`\implies x = -2`$

$`\begin{pmatrix} -5 \\ 1 \\ 0 \end{pmatrix}, \begin{pmatrix} -2 \\ 0 \\ 1 \end{pmatrix}`$이 기저이다 (2차원).

**차원 요약:**
- 행공간의 차원은 $`r`$
- 열공간의 차원은 $`r`$
- 영공간의 차원은 $`n - r`$
- $`N(A^T)`$의 차원은 $`m - r`$

### 5.6 행렬 공간과 함수 공간의 기저

**행렬 공간과 함수 공간의 기저**

"독립," "기저," "차원"이라는 용어는 열벡터에만 한정되지 않는다.

**행렬 공간 (Matrix Spaces):**

벡터 공간 $`\mathbb{M}`$은 모든 $`2 \times 2`$ 행렬을 포함한다. 차원은 4이다.

예: $`A_1, A_2, A_3, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

i) LI이다: $`x_1 A_1 + x_2 A_2 + x_3 A_3 + x_4 A_4 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}`$ $`\iff`$ $`x_1 = x_2 = x_3 = x_4 = 0`$.

ii) $`\mathbb{M}`$을 생성한다: $`\begin{pmatrix} a & b \\ c & d \end{pmatrix} = aA_1 + bA_2 + cA_3 + dA_4`$. $`A_1, A_2, A_3, A_4`$의 선형 결합이 $`\mathbb{M}`$의 임의의 행렬을 만들 수 있다.

**행렬 공간의 부분공간:**

- $`A_1, A_2, A_4 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$는 상삼각 행렬 (Upper Triangular Matrix) 부분공간 $`\mathcal{U}`$의 기저이다: $`\mathcal{U} \ni \begin{pmatrix} a & b \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

- $`A_1, A_4`$는 대각 행렬 (Diagonal Matrix) 부분공간 $`\mathbb{D}`$의 기저이다: $`\mathbb{D} \ni \begin{pmatrix} a & 0 \\ 0 & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

- $`A_1, A_4, A_2 + A_3`$는 대칭 행렬 (Symmetric Matrix)의 기저이다: $`\mathbb{S} \ni \begin{pmatrix} a & b \\ b & d \end{pmatrix} = a\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} + b\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} + d\begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}`$

**행렬 부분공간의 차원 ($`n \times n`$ 행렬에 대해):**

| 공간 | 차원 |
|:------|:----------|
| 전체 $`n \times n`$ 행렬 공간 | $`n^2`$ |
| 대각 행렬 | $`n`$ |
| 상삼각 행렬 | $`\frac{1}{2}n^2 + \frac{1}{2}n`$ |
| 대칭 행렬 | $`\frac{1}{2}n^2 + \frac{1}{2}n`$ |

상삼각 행렬의 경우: 0의 총 개수 $`= \sum_{i=1}^{n}(i-1) = \sum_{i=1}^{n} i - n = \frac{n(n+1)}{2} - n`$. 0이 아닌 원소의 총 개수 $`= n^2 - \frac{n(n+1)}{2} + n = \frac{n^2}{2} + \frac{n}{2}`$.

상삼각 행렬과 대칭 행렬의 원소 수 (차원)은 같다: $`\frac{1}{2}n^2 + \frac{1}{2}n`$.

---

**함수 공간 (Function Spaces):**

```math
\frac{d^2y}{dx^2} = 0, \quad \frac{d^2y}{dx^2} = -y, \quad \frac{d^2y}{dx^2} = y
```

이들은 2계 도함수를 포함한다. 미적분학에서 해 $`y(x)`$를 구한다:

| 상미분방정식 (ODE) | 해 | 기저 | 차원 |
|:----|:---------|:------|:----------|
| $`y'' = 0`$ | $`y = cx + d`$ | $`\lbrace1, x\rbrace`$ | 2 |
| $`y'' = -y`$ | $`y = c\sin x + d\cos x`$ | $`\lbrace\sin x, \cos x\rbrace`$ | 2 |
| $`y'' = y`$ | $`y = ce^x + de^{-x}`$ | $`\lbracee^x, e^{-x}\rbrace`$ | 2 |

기저 벡터들은 2계 도함수의 **영공간 (Nullspace)** 에 속한다.

$`y'' = 2`$는 특수해 $`y_p = x^2`$을 가진다: $`\frac{dy_p}{dx} = 2x`$, $`\frac{d^2 y_p}{dx^2} = 2`$.

따라서 일반해 (완전해)는:

```math
y(x) = y_p(x) + y_n(x) = x^2 + cx + d
```

---

<br>

## 6. 네 부분공간의 차원 (3.5)

### 6.1 차원 요약

1. 열공간 $`C(A)`$와 행공간 $`C(A^T)`$는 모두 차원 $`r`$ (= $`A`$의 계수)을 가진다.
2. $`A`$의 영공간 $`N(A)`$는 차원 $`n - r`$을 가진다.
3. $`A`$의 좌영공간 (Left Nullspace) $`N(A^T)`$는 차원 $`m - r`$을 가진다.
4. $`A`$에서 $`R_0`$로의 소거는 $`C(A)`$와 $`N(A^T)`$를 **바꾸지만**, 차원은 변하지 않는다.

### 6.2 부분공간의 직교성

"**계수 (Rank)**"와 "**차원 (Dimension)**"을 연결하려 한다:
- **계수** 는 독립 열의 수를 센다.
- **차원** 은 기저에 있는 벡터의 수이다.

$`A \in \mathbb{R}^{m \times n}`$의 계수는 네 가지 기본 부분공간 모두의 차원을 알려준다:

1. **행공간 (Row Space)**, $`C(A^T)`$는 $`\mathbb{R}^n`$의 부분공간, 차원 $`r`$. (각 행벡터 $`\in \mathbb{R}^n`$.)
2. **열공간 (Column Space)**, $`C(A)`$는 $`\mathbb{R}^m`$의 부분공간, 차원 $`r`$. (각 열벡터 $`\in \mathbb{R}^m`$.)
3. **영공간 (Nullspace)**, $`N(A)`$는 $`\mathbb{R}^n`$의 부분공간, 차원 $`n - r`$. ($`A\mathbf{x} = \mathbf{0}`$, $`\mathbf{x} \in \mathbb{R}^n`$.)
4. **좌영공간 (Left Nullspace)**, $`N(A^T)`$는 $`\mathbb{R}^m`$의 부분공간, 차원 $`m - r`$. ($`A^T\mathbf{y} = \mathbf{0}`$, $`\mathbf{y} \in \mathbb{R}^m`$.)

**직교 쌍 (Orthogonal Pairs):**

```math
C(A) \subset \mathbb{R}^m, \quad C(A^T) \subset \mathbb{R}^n
```
```math
N(A^T) \subset \mathbb{R}^m, \quad N(A) \subset \mathbb{R}^n
```

- $`C(A)`$와 $`C(A^T)`$는 $`r`$차원이다.
- $`A`$의 열공간과 행공간은 같은 차원 $`r`$을 가진다.
- $`N(A)`$는 차원 $`n - r`$.
- $`N(A^T)`$는 차원 $`m - r`$.

**$`N(A)`$는 $`C(A^T)`$ (행공간)에 수직이다:**

예: $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$

행벡터 $`\lbrace(1, 2), (3, 4)\rbrace`$가 $`A`$의 행공간을 생성한다. $`\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$의 해는 $`A`$의 영공간 $`N(A)`$에 속한다.

$`N(A) \ni \mathbf{x}`$는 내적 (Inner Product)의 의미에서 $`A`$의 행공간의 모든 벡터에 수직이다.

즉, $`\mathbf{y} \in C(A^T)`$이고 $`\mathbf{x} \in N(A)`$이면: $`\mathbf{y} \cdot \mathbf{x} = 0`$.

**$`N(A^T)`$는 $`C(A)`$ (열공간)에 수직이다:**

예: $`A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}`$

열벡터 $`\lbrace\begin{pmatrix} 1 \\ 3 \end{pmatrix}, \begin{pmatrix} 2 \\ 4 \end{pmatrix}\rbrace`$가 $`A`$의 열공간을 생성한다.

$`N(A^T) \ni \mathbf{y}`$이면 $`A^T\mathbf{y} = \mathbf{0}`$:

$`\begin{pmatrix} 1 & 3 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$

$`\Rightarrow C(A) \ni \mathbf{x} \perp \mathbf{y} \in N(A^T)`$

즉, $`\alpha\begin{pmatrix} 1 \\ 3 \end{pmatrix} + \beta\begin{pmatrix} 2 \\ 4 \end{pmatrix} \perp \mathbf{y}`$

### 6.3 R_0의 네 부분공간

**$`R_0 = \text{rref}(A)`$라 하자.** 네 차원은 $`R_0`$와 $`A`$에 대해 같다.

**예:** $`3 \times 5`$ 행렬 $`R_0`$를 고려하자:

```math
R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}
```

피벗 행: 1과 2. 피벗 열: 1과 4. $`\text{rank}(R_0) = r = 2`$.

**행공간:** 기저 벡터 $`\lbrace(1, 3, 5, 0, 7), \;(0, 0, 0, 1, 2)\rbrace`$로 생성된다. $`\dim C(R_0^T) = 2 = r`$.

**열공간:** 1번째와 4번째 열벡터 $`\left\lbrace\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\rbrace`$가 $`C(R_0)`$의 기저이다. $`\dim C(R_0) = 2 = r`$.

**영공간:** $`R_0 \mathbf{x} = \mathbf{0}`$:

```math
\begin{pmatrix} x_1 + 3x_2 + 5x_3 + 7x_5 \\ x_4 + 2x_5 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`n - r = 5 - 2 = 3`$, 3개의 자유 변수.

i) $`x_2 = 1, x_3 = 0, x_5 = 0`$으로 취하면: $`x_1 + 3 = 0 \Rightarrow x_1 = -3`$, $`x_4 = 0`$. $`\mathbf{x} = \begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}`$

ii) $`x_2 = 0, x_3 = 1, x_5 = 0`$으로 취하면: $`x_1 + 5 = 0 \Rightarrow x_1 = -5`$, $`x_4 = 0`$. $`\mathbf{x} = \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}`$

iii) $`x_2 = 0, x_3 = 0, x_5 = 1`$로 취하면: $`x_1 + 7 = 0 \Rightarrow x_1 = -7`$, $`x_4 + 2 = 0 \Rightarrow x_4 = -2`$. $`\mathbf{x} = \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}`$

3개의 특수해가 $`N(R_0)`$의 기저 $`\left\lbrace\begin{pmatrix} -3 \\ 1 \\ 0 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -5 \\ 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} -7 \\ 0 \\ 0 \\ -2 \\ 1 \end{pmatrix}\right\rbrace`$를 형성한다.

$`\dim N(R_0) = 3 = 5 - 2 = n - r`$.

**좌영공간:** $`R_0^T \mathbf{y} = \mathbf{0}`$:

```math
\begin{pmatrix} 1 & 0 & 0 \\ 3 & 0 & 0 \\ 5 & 0 & 0 \\ 0 & 1 & 0 \\ 7 & 2 & 0 \end{pmatrix} \begin{pmatrix} y_1 \\ y_2 \\ y_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix} \iff \begin{pmatrix} y_1 \\ 3y_1 \\ 5y_1 \\ y_2 \\ 7y_1 + 2y_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}
```

$`m - r = 3 - 2 = 1`$, 1개의 자유 변수.

$`y_3 = 1`$로 취하면: $`\Rightarrow y_1 = y_2 = 0`$. $`\therefore \mathbf{y} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}`$

$`\left\lbrace\begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\right\rbrace`$이 $`N(R_0^T)`$의 기저이다.

$`\dim N(R_0^T) = 1 = m - r = 3 - 2`$.

$`R_0^T \mathbf{y} = \mathbf{0} \iff \mathbf{y}^T R_0 = \mathbf{0}`$. $`\mathbf{y}^T`$는 $`R_0`$의 '왼쪽'에 있는 행벡터이다.

**직교성 요약:**

```math
C(A) \subset \mathbb{R}^m \perp N(A^T) \subset \mathbb{R}^m
```
```math
C(A^T) \subset \mathbb{R}^n \perp N(A) \subset \mathbb{R}^n
```

$`\mathbb{R}^n`$에서: 행공간과 영공간은 $`r`$과 $`n - r`$ 차원을 가진다.

$`\mathbb{R}^m`$에서: 열공간과 좌영공간은 $`r`$과 $`m - r`$ 차원을 가진다.

### 6.4 A와 R_0의 관계

**$`A`$의 네 부분공간 차원은 $`R_0`$와 같다.**

```math
A = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 1 & 3 & 5 & 1 & 9 \end{pmatrix} \longrightarrow R_0 = \begin{pmatrix} 1 & 3 & 5 & 0 & 7 \\ 0 & 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 & 0 \end{pmatrix}
```

$`A`$는 $`R_0`$와 같은 행공간을 가지지만, 열공간은 $`R_0`$와 다르다.

$`C(A)`$의 기저: $`\left\lbrace\begin{pmatrix} 1 \\ 0 \\ 1 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}\right\rbrace`$ ($`R_0`$가 아닌 $`A`$의 피벗 열)

$`C(R_0)`$의 기저: $`\left\lbrace\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}, \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}\right\rbrace`$

($`\mathbb{R}^3`$에서 서로 다른 평면이지만, 같은 차원.)

**핵심 관계:**

1. $`A`$는 $`R_0`$, $`R`$과 **같은 행공간** 을 가진다: $`C(A^T) = C(R_0^T) = C(R^T)`$, 차원 $`r`$.

2. $`A`$의 열공간은 차원 $`r`$을 가진다: $`\dim C(A) = \dim C(A^T)`$. $`C(A) \neq C(R_0)`$이지만, $`\dim C(A) = \dim C(R_0) = r`$.

3. $`A`$는 $`R_0`$와 **같은 영공간** 을 가진다: $`A\mathbf{x} = \mathbf{0} \iff R_0 \mathbf{x} = \mathbf{0}`$. **소거는 해를 바꾸지 않는다.** 특수해가 기저를 형성한다. $`n - r`$개의 자유 변수 $`\implies`$ $`\dim N(A) = n - r`$.

```math
\dim C(A) + \dim N(A) = r + (n - r) = n
```

4. $`A`$의 좌영공간 $`N(A^T)`$: $`\dim C(A^T) + \dim N(A^T) = r + (m - r) = m`$.

### 6.5 선형대수학의 기본 정리

```math
\boxed{\textbf{선형대수학의 기본 정리 (Fundamental Theorem of Linear Algebra)}}
```

> i) $`\dim C(A) = \dim C(A^T) = r`$
>
> ii) $`\dim N(A) = n - r`$, $`\quad \dim N(A^T) = m - r`$

**$`A`$의 네 부분공간:**

```math
C(A^T) \subset \mathbb{R}^n \quad \perp \quad N(A) \subset \mathbb{R}^n
```
```math
C(A) \subset \mathbb{R}^m \quad \perp \quad N(A^T) \subset \mathbb{R}^m
```

$`\mathbb{R}^n`$에서: $`C(A^T)`$ (차원 $`r`$)과 $`N(A)`$ (차원 $`n - r`$)는 직교 여공간 (Orthogonal Complements)이다.

$`\mathbb{R}^m`$에서: $`C(A)`$ (차원 $`r`$)과 $`N(A^T)`$ (차원 $`m - r`$)는 직교 여공간이다.

---

<br>

## 요약

| 개념 | 핵심 내용 |
|:--------|:---------|
| 벡터 공간 (Vector Space) | 덧셈과 스칼라 곱에 닫혀 있으며 8가지 공리를 만족하는 집합 |
| 체 (Field) | $`+, -, \times, \div`$가 정의된 집합 ($`\mathbb{R}`$, $`\mathbb{C}`$) |
| 부분공간 (Subspace) | 그 자체가 벡터 공간인 벡터 공간의 부분집합 ($`\mathbf{0}`$ 포함 필수) |
| 열공간 $`C(A)`$ | $`A`$의 열들의 모든 선형 결합; $`A\mathbf{x} = \mathbf{b}`$가 풀리려면 $`\mathbf{b} \in C(A)`$ |
| 행공간 $`C(A^T)`$ | $`A^T`$의 열공간; $`A`$의 행들이 생성 |
| 영공간 $`N(A)`$ | $`A\mathbf{x} = \mathbf{0}`$의 모든 해; $`\mathbb{R}^n`$의 부분공간 |
| 좌영공간 $`N(A^T)`$ | $`A^T\mathbf{y} = \mathbf{0}`$인 모든 $`\mathbf{y}`$; $`\mathbb{R}^m`$의 부분공간 |
| 기약행사다리꼴 $`R_0 = \text{rref}(A)`$ | 피벗 열에 $`I`$, 자유 열에 $`F`$를 포함 |
| $`A = CR`$ | $`C`$ = $`A`$의 독립 열; $`R = (I \; F)`$ = 기약행사다리꼴 (영행 제거) |
| 특수해 (Special Solutions) | $`\begin{pmatrix} -F \\ I \end{pmatrix}`$의 열들; $`N(A)`$의 기저를 형성 |
| 완전해 (Complete Solution) | $`\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n`$ (특수해 + 영공간) |
| 특수해 $`\mathbf{x}_p`$ (Particular Solution) | 모든 자유 변수를 0으로 놓고 $`R_0\mathbf{x} = \mathbf{d}`$ 풀기 |
| 열 완전 계수 ($`r = n`$) | $`N(A) = \lbrace\mathbf{0}\rbrace`$; $`A\mathbf{x} = \mathbf{b}`$의 해가 최대 1개; $`R_0 = \begin{pmatrix} I \\ 0 \end{pmatrix}`$ |
| 행 완전 계수 ($`r = m`$) | $`A\mathbf{x} = \mathbf{b}`$가 항상 풀림; $`C(A) = \mathbb{R}^m`$; $`R = (I \; F)`$ |
| 선형 독립 (Linear Independence) | $`c_1\mathbf{v}_1 + \cdots + c_n\mathbf{v}_n = \mathbf{0}`$이 모든 $`c_i = 0`$일 때만 성립 |
| 생성 (Spanning) | 벡터들이 $`\mathbb{S}`$를 생성: $`\mathbb{S}`$ = 그 벡터들의 모든 결합 |
| 기저 (Basis) | 선형 독립이면서 공간을 생성하는 벡터들; 표현이 유일 |
| 차원 (Dimension) | 모든 기저에 있는 벡터의 수; 기저에 관계없이 불변 |
| $`\dim C(A) = \dim C(A^T) = r`$ | 열공간과 행공간은 같은 차원 (계수) |
| $`\dim N(A) = n - r`$ | 영공간 차원은 자유 변수의 수 |
| $`\dim N(A^T) = m - r`$ | 좌영공간 차원 |
| $`r + (n - r) = n`$ | 행공간 + 영공간이 $`\mathbb{R}^n`$을 채움 |
| $`r + (m - r) = m`$ | 열공간 + 좌영공간이 $`\mathbb{R}^m`$을 채움 |
| 직교성 (Orthogonality) | $`N(A) \perp C(A^T)`$ ($`\mathbb{R}^n`$); $`N(A^T) \perp C(A)`$ ($`\mathbb{R}^m`$) |
| 행렬 공간 차원 | $`n \times n`$: $`n^2`$; 대각: $`n`$; 상삼각: $`\frac{n^2+n}{2}`$; 대칭: $`\frac{n^2+n}{2}`$ |
| 함수 공간 (Function Space) | $`y'' = 0, -y, y`$의 해가 기저 $`\lbrace1,x\rbrace`$, $`\lbrace\sin x, \cos x\rbrace`$, $`\lbracee^x, e^{-x}\rbrace`$인 2차원 공간 형성 |
| 기본 정리 (Fundamental Theorem) | $`\dim C(A) = \dim C(A^T) = r`$; $`\dim N(A) = n-r`$; $`\dim N(A^T) = m-r`$ |

---
