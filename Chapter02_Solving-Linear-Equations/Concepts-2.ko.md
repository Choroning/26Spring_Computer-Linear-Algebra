# 제2장 강의 — 연립일차방정식 풀기 — Part 2/2

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.**
>
> [Part 1](Concepts.ko.md) · **Part 2**

---

<br>

## 목차

- [5. 치환과 전치 (2.4)](#5-치환과-전치-24)
  - [5.1 치환 행렬](#51-치환-행렬)
  - [5.2 치환 행렬의 성질](#52-치환-행렬의-성질)
  - [5.3 PA = LU 분해](#53-pa--lu-분해)
  - [5.4 부분 피벗팅](#54-부분-피벗팅)
  - [5.5 PAQ: 행과 열 치환](#55-paq-행과-열-치환)
  - [5.6 A의 전치](#56-a의-전치)
  - [5.7 내적과 전치](#57-내적과-전치)
  - [5.8 대칭 행렬](#58-대칭-행렬)
  - [5.9 대칭 곱과 LDL^T](#59-대칭-곱과-ldlt)
- [6. 도함수와 유한 차분 행렬 (2.5)](#6-도함수와-유한-차분-행렬-25)
  - [6.1 테일러 급수와 근사](#61-테일러-급수와-근사)
  - [6.2 차분으로부터의 도함수](#62-차분으로부터의-도함수)
  - [6.3 이차 차분 행렬 K, T, B](#63-이차-차분-행렬-k-t-b)
  - [6.4 K의 성질](#64-k의-성질)
  - [6.5 자유-고정 행렬 T](#65-자유-고정-행렬-t)
  - [6.6 자유-자유 행렬 B](#66-자유-자유-행렬-b)
- [요약](#요약)

---

<br>

## 5. 치환과 전치 (2.4)

### 5.1 치환 행렬

**치환 행렬** (permutation matrix) $`P`$는 $`I \in \mathbb{R}^{n \times n}`$과 같은 행을 가진다.

$`n!`$가지 서로 다른 순서가 있다.

**예제:** $`P \in \mathbb{R}^{3 \times 3}`$: 3개의 행. 1행에: 3가지; 2행에: 2가지; 3행에: 1가지 $`\implies 3! = 6`$가지 순서.

$`P`$ 곱하기 $`\mathbf{x}`$는 성분 $`x_1`$부터 $`x_n`$까지를 새로운 순서로 배치한다.

그리고 $`P^T`$는 $`P^{-1}`$과 같다.

**예제:**

```math
P\mathbf{x} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \end{pmatrix}
```

```math
P^T = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}
```

$`P^{-1} = ?`$ 가우스-조르단 방법 $`(A|I) \Rightarrow (I|A^{-1})`$ 사용:

```math
\begin{pmatrix} 0 & 0 & 1 & | & 1 & 0 & 0 \\ 1 & 0 & 0 & | & 0 & 1 & 0 \\ 0 & 1 & 0 & | & 0 & 0 & 1 \end{pmatrix} \xrightarrow{\text{행 교환}} \begin{pmatrix} 1 & 0 & 0 & | & 0 & 1 & 0 \\ 0 & 1 & 0 & | & 0 & 0 & 1 \\ 0 & 0 & 1 & | & 1 & 0 & 0 \end{pmatrix}
```

```math
P^{-1} = P^T
```

**전치의 성질:**

- $`A`$의 열은 $`A^T`$의 행이다.
- $`A\mathbf{x}`$와 $`AB`$의 전치는 $`\mathbf{x}^T A^T`$와 $`B^T A^T`$이다.

**내적의 성질:**

```math
A\mathbf{x} \cdot \mathbf{y} = \mathbf{x} \cdot A^T\mathbf{y}
```

왜냐하면 $`(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T A^T \mathbf{y} = \mathbf{x}^T(A^T\mathbf{y})`$

**대칭 행렬:** $`S^T = S`$. 곱 $`A^T A`$와 $`AA^T`$는 항상 대칭이다.

---

### 5.2 치환 행렬의 성질

치환 행렬은 모든 행에 정확히 하나의 1을, 모든 열에 정확히 하나의 1을 가진다.

$`P`$를 벡터 $`\mathbf{x}`$에 곱하면, 성분의 순서가 바뀐다:

```math
P\mathbf{x} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \end{pmatrix}
```

$`P`$는 $`x_1`$을 두 번째 위치로 이동시킨다.

```math
PP\mathbf{x} = P^2\mathbf{x} = \begin{pmatrix} x_2 \\ x_3 \\ x_1 \end{pmatrix}
```

$`P^2`$는 $`x_1`$을 세 번째 위치로 이동시킨다.

```math
PPP\mathbf{x} = P^3\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = I\mathbf{x}
```

$`P^3`$은 $`x_1`$을 원래 위치로 되돌린다. $`P^3 = I`$.

**4x4 치환 행렬을 고려하자:**

**(a)** $`P`$는 $`\mathbf{x}`$의 순서를 뒤집는다:

```math
\begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_4 \\ x_3 \\ x_2 \\ x_1 \end{pmatrix}
```

$`P(P\mathbf{x}) = \mathbf{x} \implies P^2 = I`$

**(b)** $`P`$는 $`x_4`$ 위치를 바꾸지 않는다:

```math
\begin{pmatrix} 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \\ x_4 \end{pmatrix}
```

$`P(P(P\mathbf{x})) = \mathbf{x} \implies P^3 = I`$

**(c)** $`P`$는 원소를 순환적으로 이동시킨다:

```math
\begin{pmatrix} 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_4 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix}
```

$`PPPP\mathbf{x} = P^4\mathbf{x} = I\mathbf{x} \implies P^4 = I`$

**(d)** 짝수-홀수 분리:

```math
\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_0 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_0 \\ x_2 \\ x_1 \\ x_3 \end{pmatrix} \quad \text{(짝수, 홀수 분리)}
```

$`P^2 = I`$

이것은 짝수 인덱스와 홀수 인덱스 항목을 분리하는 8x8 치환 행렬을 사용하여 $`\mathbf{x} \in \mathbb{R}^8`$ 벡터로 확장된다.

**$`P^T P = I`$의 증명:**

임의의 $`P`$의 행은 $`P^{-1} = P^T`$의 열이다.

```math
P^T P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}\begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}
```

```math
= \mathbf{h}_1\mathbf{h}_1^T + \mathbf{h}_2\mathbf{h}_2^T + \mathbf{h}_3\mathbf{h}_3^T
```

$`\mathbf{h}_i`$가 표준 기저 벡터 (canonical unit vector)이므로:

```math
= \begin{pmatrix} 1 & & \\ & 1 & \\ & & 1 \end{pmatrix} = I
```

**치환 행렬의 성질:**

1. 치환 행렬 $`P`$는 각 행에 정확히 하나의 1을, 각 열에 정확히 하나의 1을 가진다.

2. $`P`$의 열들은 **직교** (orthogonal)한다. (열 사이의 내적이 모두 0이다.)

3. 치환의 곱 $`P_1 P_2`$는 치환이다. $`P`$의 역행렬도 마찬가지이다.

4. $`A`$가 가역이면, 행을 미리 정렬하는 치환 $`P`$가 존재하여, $`PA`$에 대한 소거가 피벗 위치에서 0을 만나지 않는다:

```math
PA = LU
```

---

### 5.3 PA = LU 분해

**$`P`$에 의한 행 교환:**

행렬 $`A`$를 고려하자:

```math
A = \begin{pmatrix} 1 & 2 & a \\ 2 & 4 & b \\ 3 & 7 & c \end{pmatrix}
```

1행에서 1을 피벗으로 취한다:

```math
R_2 - 2R_1, \quad R_3 - 3R_1
```

```math
EA = \begin{pmatrix} 1 & 2 & a \\ 0 & 0 & b - 2a \\ 0 & 1 & c - 3a \end{pmatrix}
```

0 피벗으로 인해, 2행과 3행을 교환한다:

```math
PEA = \begin{pmatrix} 1 & 2 & a \\ 0 & 1 & c - 3a \\ 0 & 0 & b - 2a \end{pmatrix} = U
```

$`A`$는 $`b - 2a \neq 0`$일 때에만 가역이다. $`b = 2a`$이면, $`\text{rank}(A) = 2 < 3`$이고 $`A`$는 가역이 아니다.

**먼저 2행과 3행을 교환할 수 있다:**

```math
PA = \begin{pmatrix} 1 & 2 & a \\ 3 & 7 & c \\ 2 & 4 & b \end{pmatrix}
```

```math
\xrightarrow{R_2 - 3R_1, R_3 - 2R_1}
```

```math
EPA = \begin{pmatrix} 1 & 2 & a \\ 0 & 1 & c - 3a \\ 0 & 0 & b - 2a \end{pmatrix} = U
```

```math
\therefore PA = E^{-1}U = LU
```

**Daniel Drucker의 $`P`$ 추적 방법:**

행 인덱스를 추적하는 열로 행렬을 확대한다:

```math
\begin{pmatrix} 1 & 2 & a & | & 1 \\ 2 & 4 & b & | & 2 \\ 3 & 7 & c & | & 3 \end{pmatrix} \xrightarrow{\text{소거}} \begin{pmatrix} 1 & 2 & a & | & 1 \\ 0 & 0 & b-2a & | & 2 \\ 0 & 1 & c-3a & | & 3 \end{pmatrix} \xrightarrow{\text{교환}} \begin{pmatrix} 1 & 2 & a & | & 1 \\ 0 & 1 & c-3a & | & 3 \\ 0 & 0 & b-2a & | & 2 \end{pmatrix}
```

최종 열은 $`P_{132}`$를 준다:

```math
P_{132} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}
```

---

### 5.4 부분 피벗팅

반올림 오차를 줄이기 위한 **"부분 피벗팅"** (Partial Pivoting).

피벗에 **가능한 가장 큰 수** 를 생성하도록 행을 교환하면 계산이 더 안정적이다.

**예제:**

```math
\begin{pmatrix} 1 & 2 & a & | & 1 \\ 2 & 4 & b & | & 2 \\ 3 & 7 & c & | & 3 \end{pmatrix}
```

```math
\xrightarrow{R_3 \leftrightarrow R_1} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 2 & 4 & b & | & 2 \\ 1 & 2 & a & | & 1 \end{pmatrix}
```

```math
\xrightarrow{R_2 - \frac{2}{3}R_1, R_3 - \frac{1}{3}R_1} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 0 & -\frac{2}{3} & b - \frac{2}{3}c & | & 2 \\ 0 & -\frac{1}{3} & a - \frac{1}{3}c & | & 1 \end{pmatrix}
```

```math
\xrightarrow{R_3' - R_2'(\frac{1}{2})} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 0 & -\frac{2}{3} & b - \frac{2}{3}c & | & 2 \\ 0 & 0 & a - \frac{1}{2}b & | & 1 \end{pmatrix} = U
```

각 피벗을 그 아래의 모든 수보다 크게 만들면 $`L`$의 모든 원소가 $`\leq 1`$이 된다:

```math
L = \begin{pmatrix} 1 & 0 & 0 \\ 2/3 & 1 & 0 \\ 1/3 & 1/2 & 1 \end{pmatrix}
```

```math
P_{321} = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}
```

---

### 5.5 PAQ: 행과 열 치환

$`PAQ`$는 행 치환 $`P`$와 열 치환 $`Q`$를 가진다.

$`A \in \mathbb{R}^{3 \times 3}`$에서 시작한다. $`P \in \mathbb{R}^{3 \times 3}`$로 행을 재정렬한다.

```math
P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}, \quad PA = \begin{pmatrix} a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \\ a_{11} & a_{12} & a_{13} \end{pmatrix}
```

오른쪽에서 $`Q = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}`$을 곱한다:

```math
(PA)Q = \begin{pmatrix} a_{23} & a_{22} & a_{21} \\ a_{33} & a_{32} & a_{31} \\ a_{13} & a_{12} & a_{11} \end{pmatrix}
```

열 치환 $`Q`$는 열을 재정렬한다.

$`A`$의 열공간은 $`PA`$의 열공간과 같은가?

```math
C(A) \stackrel{?}{=} C(PA)
```

**그렇다**, $`P`$는 일차 관계를 바꾸지 않는다. 따라서 $`C(A) = C(PA)`$.

**Q:** 행렬 $`A`$는 9개의 수를 가진다. $`A`$의 9개의 수를 몇 가지 다른 방법으로 배열할 수 있는가? **A:** $`9!`$

$`P, Q`$는 $`PAQ`$에 대해 $`6 \times 6 = 36`$가지로 수를 배열할 수 있다. $`PAQ`$는 $`C(A) = C(PAQ)`$를 만족하는 매우 특별한 것이다.

---

### 5.6 A의 전치

$`A`$의 **전치** (transpose)는 $`A^T`$로 표기한다. $`A^T`$의 열은 $`A`$의 행이다.

```math
A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & \vdots & \ddots & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}_{n \times m}
```

```math
A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}_{m \times n}
```

**예제:**

```math
A = \begin{pmatrix} 1 & 2 & 3 \\ 0 & 0 & 4 \end{pmatrix}, \quad A^T = \begin{pmatrix} 1 & 0 \\ 2 & 0 \\ 3 & 4 \end{pmatrix}
```

행렬은 주 대각선을 기준으로 "뒤집힌다".

```math
(A^T)_{ij} = A_{ji}
```

**전치의 규칙:**

- **합:** $`(A + B)^T = A^T + B^T`$
- **곱:** $`(AB)^T = B^T A^T`$ (역순)
- **역:** $`(A^{-1})^T = (A^T)^{-1}`$

**$`(A\mathbf{x})^T = \mathbf{x}^T A^T`$의 증명:**

```math
A\mathbf{x} = \begin{pmatrix} \sum_{j=1}^n a_{1j}x_j \\ \sum_{j=1}^n a_{2j}x_j \\ \vdots \\ \sum_{j=1}^n a_{mj}x_j \end{pmatrix}
```

```math
(A\mathbf{x})^T = \left(\sum_{j=1}^n a_{1j}x_j \quad \sum_{j=1}^n a_{2j}x_j \quad \cdots \quad \sum_{j=1}^n a_{mj}x_j\right)
```

```math
= (x_1, x_2, \ldots, x_n)\begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}_{n \times m}
```

```math
= \mathbf{x}^T A^T
```

**$`(AB)^T = B^T A^T`$의 증명:**

$`B = \begin{pmatrix} | & | & & | \\ \mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_p \\ | & | & & | \end{pmatrix}_{n \times p}`$으로 해석한다.

```math
(AB)^T = (A\mathbf{x}_1 \quad A\mathbf{x}_2 \quad \cdots \quad A\mathbf{x}_p)^T = \begin{pmatrix} \mathbf{x}_1^T A^T \\ \mathbf{x}_2^T A^T \\ \vdots \\ \mathbf{x}_p^T A^T \end{pmatrix} = \begin{pmatrix} \mathbf{x}_1^T \\ \mathbf{x}_2^T \\ \vdots \\ \mathbf{x}_p^T \end{pmatrix} A^T = B^T A^T
```

**예제:**

```math
AB = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 5 & 0 \\ 4 & 1 \end{pmatrix} = \begin{pmatrix} 5 & 0 \\ 9 & 1 \end{pmatrix}
```

```math
B^T A^T = \begin{pmatrix} 5 & 4 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 5 & 9 \\ 0 & 1 \end{pmatrix} = (AB)^T
```

또한: $`(ABC)^T = (AB \cdot C)^T = C^T(AB)^T = C^T B^T A^T`$. 역순 규칙이 여전히 성립한다.

**$`(A^{-1})^T = (A^T)^{-1}`$의 증명:**

$`A^{-1}A = I`$를 고려하면:

```math
\Rightarrow (A^{-1}A)^T = I^T = I
```
```math
\Leftrightarrow A^T(A^{-1})^T = I
```

$`\therefore (A^{-1})^T`$는 $`A^T`$의 역행렬이다:

```math
(A^{-1})^T = (A^T)^{-1}
```

이는 $`A`$가 가역일 때 정확히 $`A^T`$도 가역임을 의미한다.

**예제:**

```math
A = \begin{pmatrix} 1 & 0 \\ 6 & 1 \end{pmatrix}, \quad A^{-1} = \begin{pmatrix} 1 & 0 \\ -6 & 1 \end{pmatrix} \implies (A^{-1})^T = \begin{pmatrix} 1 & -6 \\ 0 & 1 \end{pmatrix}
```

```math
A^T = \begin{pmatrix} 1 & 6 \\ 0 & 1 \end{pmatrix} \implies (A^T)^{-1} = \begin{pmatrix} 1 & -6 \\ 0 & 1 \end{pmatrix}
```

일치한다: $`(A^{-1})^T = (A^T)^{-1}`$.

---

### 5.7 내적과 전치

**내적의 의미.**

$`\mathbf{x}`$와 $`\mathbf{y}`$의 점곱(dot product) (내적, inner product)은 수 $`x_i y_i`$의 합이다:

```math
\mathbf{x} \cdot \mathbf{y} = \sum x_i y_i
```

$`\mathbf{x}, \mathbf{y} \in \mathbb{R}^n`$라 하자:

```math
\mathbf{x} \cdot \mathbf{y} = \sum_{i=1}^n x_i y_i = \mathbf{x}^T \mathbf{y} \quad (1 \times n)(n \times 1) = (1 \times 1)
```

```math
\mathbf{x}\mathbf{y}^T = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix}(y_1, y_2, \ldots, y_n) = \mathbf{x} \otimes \mathbf{y} \quad (n \times 1)(1 \times n) = (n \times n)
```

이것은 **외적** (outer product) (랭크 1 곱, 랭크 1 행렬)이다.

**내적의 예:**

```math
\text{일} = \text{힘} \cdot \text{거리} = \mathbf{f}^T \mathbf{d} \quad [J] = [N] \cdot [m]
```

```math
\text{수입} = \text{수량} \cdot \text{가격} = \mathbf{q}^T \mathbf{p}
```

**$`A^T`$의 더 나은 정의:**

$`A^T`$를 행렬 $`A`$를 주 대각선을 기준으로 뒤집어 정의하지만, $`A^T`$의 더 나은 정의는 다음 내적을 같게 만드는 행렬이라는 것이다:

```math
(A\mathbf{x}) \cdot \mathbf{y} = \mathbf{x} \cdot (A^T\mathbf{y})
```

즉,

```math
(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T A^T \mathbf{y}
```

**예제 1:**

```math
A = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}, \quad \mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix} y_1 \\ y_2 \end{pmatrix}
```

```math
A\mathbf{x} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \end{pmatrix}
```

```math
(A\mathbf{x})^T\mathbf{y} = (x_2 - x_1)y_1 + (x_3 - x_2)y_2
```

```math
A^T\mathbf{y} = \begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} -y_1 \\ y_1 - y_2 \\ y_2 \end{pmatrix}
```

```math
\mathbf{x}^T(A^T\mathbf{y}) = (x_1, x_2, x_3)\begin{pmatrix} -y_1 \\ y_1 - y_2 \\ y_2 \end{pmatrix} = -x_1y_1 + x_2(y_1 - y_2) + x_3y_2
```

양쪽이 같다. $`\checkmark`$

**예제 2: 함수의 내적.**

```math
\mathbf{x}^T\mathbf{y} = x_1 y_1 + x_2 y_2 + \cdots + x_n y_n
```

연속 세계에서:

```math
(x, y) := \int_{-\infty}^{\infty} x(t) \, y(t) \, dt
```

마찬가지로, $`(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T(A^T\mathbf{y})`$:

```math
(Ax, y) = (x, A^T y)
```

$`A = \frac{d}{dt}`$, $`A^T = -\frac{d}{dt} = -A`$로 놓자.

```math
\left(\frac{dx}{dt}, y\right) = \left(x, -\frac{dy}{dt}\right)
```

```math
\int_{-\infty}^{\infty} \frac{dx}{dt} y \, dt = \int_{-\infty}^{\infty} x \left(-\frac{dy}{dt}\right) dt
```

**부분적분 (IBP)을 통한 증명:**

```math
(f(t)g(t))' = f'(t)g(t) + f(t)g'(t)
```

```math
\int f'g \, dt = \int (fg)' \, dt - \int fg' \, dt = (fg)\Big|_{-\infty}^{\infty} - \int fg' \, dt
```

$`f(\infty) = f(-\infty) = 0`$이라 가정하면:

```math
\boxed{\int_{-\infty}^{\infty} f'g \, dt = -\int_{-\infty}^{\infty} fg' \, dt}
```

미분은 **반대칭** (anti-symmetric)이다: $`A^T = -A`$.

대칭 행렬은 $`A^T = A`$를 가진다.

---

### 5.8 대칭 행렬

**대칭 행렬** (symmetric matrix)은 $`S^T = S`$를 만족한다:

```math
\implies (S^T)_{ij} = S_{ij} = S_{ji}
```

**예제:**

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix} = S^T, \quad D = \begin{pmatrix} 1 & 0 \\ 0 & 10 \end{pmatrix} = D^T
```

**대칭 행렬의 역행렬은 대칭 행렬이다:**

```math
(S^{-1})^T = (S^T)^{-1} = S^{-1}
```

$`\Rightarrow`$ $`S`$가 가역일 때, $`S^{-1}`$은 대칭이다.

**예제:**

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}, \quad S^{-1} = \begin{pmatrix} 5 & -2 \\ -2 & 1 \end{pmatrix} = (S^{-1})^T
```

---

### 5.9 대칭 곱과 LDL^T

**대칭 곱 $`A^T A`$, $`AA^T`$, $`LDL^T`$:**

$`A \in \mathbb{R}^{m \times n}`$에 대해:

- $`A^T A \in \mathbb{R}^{n \times n}`$: $`(A^T A)^T = A^T A \implies`$ **대칭**
- $`AA^T \in \mathbb{R}^{m \times m}`$: $`(AA^T)^T = AA^T \implies`$ **대칭**

**예제 3:**

```math
A = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}
```

```math
AA^T = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}
```

```math
A^T A = \begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}
```

**소거에서의 대칭 행렬:** $`S^T = S`$는 소거를 **두 배 빠르게** 만든다.

```math
S = \begin{pmatrix} 1 & 2 \\ 2 & 7 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix} = U
```

```math
L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}
```

```math
S = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}
```

$`LU`$ 분해에서는 $`S`$의 대칭성이 보이지 않는다. 대칭 행렬 $`S`$에 대해, $`U`$를 $`D`$와 $`L^T`$로 더 분해할 수 있다:

```math
U = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix} = \begin{pmatrix} d_1 & \\ & d_2 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & \\ & 3 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}
```

```math
S = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 3 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = LDL^T
```

**예제 4: 안장점 행렬 (Saddle-point matrix).**

직사각 행렬 $`A \in \mathbb{R}^{m \times n}`$에 대해, **안장점 행렬** $`S`$는 대칭이며 중요하다:

```math
S = \begin{pmatrix} I_{m \times m} & A_{m \times n} \\ A^T_{n \times m} & 0_{n \times n} \end{pmatrix} = S^T \quad (m+n) \times (m+n)
```

블록 소거: $`R_2 - A^T R_1`$:

```math
ES = \begin{pmatrix} I & A \\ 0 & -A^T A \end{pmatrix} = U
```

**블록 분해:**

```math
S = \begin{pmatrix} I & 0 \\ A^T & I \end{pmatrix}\begin{pmatrix} I & 0 \\ 0 & -A^T A \end{pmatrix}\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} = LDL^T
```

$`S`$가 가역 $`\iff`$ $`A^T A`$가 가역 $`\iff`$ $`\text{rank}(A^T A) = n`$ $`\iff`$ $`\mathbf{x} \neq \mathbf{0}`$일 때 $`A\mathbf{x} \neq \mathbf{0}`$ $`\iff`$ $`A`$의 열들이 선형독립.

---

<br>

## 6. 도함수와 유한 차분 행렬 (2.5)

### 6.1 테일러 급수와 근사

행렬은 **도함수** (derivatives)를 모방한다. 도함수는 공간의 한 점 $`x`$에서 또는 시간의 한 순간 $`t`$에서 무슨 일이 일어나고 있는지 알려준다.

**예제:** $`y = x^2 + 2`$

```math
\frac{dy}{dx} = 2x \xrightarrow{x=1} 2 > 0
```

```math
\frac{d^2y}{dx^2} = 2 > 0
```

$`y`$의 그래프는 위로 굽는다 (기울기가 증가한다).

**$`\Delta y = y(x+h) - y(x)`$를 고려하자:**

1. 차이는 대략: $`\Delta y \approx h \cdot \frac{dy}{dx}\Big|_{x=x_1}`$

2. 더 나은 근사: $`\Delta y \approx h \cdot \frac{dy}{dx}\Big|_{x=x_1} + \frac{h^2}{2}\frac{d^2y}{dx^2}\Big|_{x=x_1}`$ (접선 + 포물선)

3. 정확한 $`\Delta y`$는 적분이다: $`\Delta y = y(x+h) - y(x) = \int_x^{x+h} \frac{dy}{dx} \, dx`$

$`\Delta y`$의 정확도는 도함수 항을 추가함으로써 증가한다.

**테일러 급수 (Taylor Series):**

```math
y(x+h) = y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + \cdots + \frac{h^n}{n!}\frac{d^{(n)}y}{dx^n} + \cdots
```

**예제:** $`e^x`$ (전해석 함수, entire analytic function):

```math
e^{x+h} = e^x + h \cdot e^x + \frac{h^2}{2} \cdot e^x + \cdots + \frac{h^n}{n!} e^x + \cdots
```

```math
= e^x\left(1 + h + \frac{h^2}{2} + \cdots + \frac{h^n}{n!} + \cdots\right)
```

```math
\therefore e^h = 1 + h + \frac{h^2}{2} + \cdots + \frac{h^n}{n!} + \cdots
```

---

### 6.2 차분으로부터의 도함수

**공식 뒤집기: 차분으로부터의 도함수**

접선 포물선에서 시작한다:

```math
y(x+h) \approx y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2}
```

**Q:** $`y(x)`$와 $`y(x+h)`$를 알면, $`\frac{dy}{dx}`$를 어떻게 추정하는가?

**Q:** $`y(x-h)`$를 알면, $`\frac{dy}{dx}`$, $`\frac{d^2y}{dx^2}`$를 추정할 수 있는가?

**유한 차분법** (Finite difference method)을 사용하여 도함수를 근사할 수 있다.

**전진 차분 (Forward Difference):**

```math
y(x+h) \approx y(x) + h\frac{dy}{dx}
```

```math
\implies \frac{dy}{dx} \approx \frac{y(x+h) - y(x)}{h} \quad \text{(전진 차분)}
```

이것은 **1차** 근사이다:

```math
y(x+h) = y(x) + h\frac{dy}{dx} + O(h^2)
```

```math
\implies \frac{dy}{dx} = \frac{y(x+h) - y(x)}{h} + O(h) \quad \text{(절단 오차)}
```

**후진 차분 (Backward Difference):**

```math
y(x-h) \approx y(x) - h\frac{dy}{dx}
```

```math
\implies \frac{dy}{dx} \approx \frac{y(x) - y(x-h)}{h}
```

**Q:** $`\frac{dy}{dx}`$에 대한 근사의 정확도를 높일 수 있는가?

**중심 차분 (Centered Difference) (2차 정확도):**

```math
y(x+h) = y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + O(h^3)
```

```math
y(x-h) = y(x) - h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + O(h^3)
```

빼면:

```math
y(x+h) - y(x-h) = 2h\frac{dy}{dx} + O(h^3)
```

```math
\implies \frac{dy}{dx} = \frac{y(x+h) - y(x-h)}{2h} + O(h^2) \approx \frac{y(x+h) - y(x-h)}{2h}
```

이 중심 차분 공식은 **2차 정확도** 를 가진다.

**이차 차분 ($`\frac{d^2y}{dx^2}`$의 근사):**

두 테일러 전개를 더하면:

```math
y(x+h) + y(x-h) = 2y(x) + h^2\frac{d^2y}{dx^2} + O(h^4)
```

```math
\implies \frac{d^2y}{dx^2} = \frac{y(x+h) - 2y(x) + y(x-h)}{h^2} + O(h^2)
```

```math
\approx \frac{y(x+h) - 2y(x) + y(x-h)}{h^2} \quad \text{(이차 차분)}
```

```math
\frac{d^2y}{dx^2} \approx \frac{1}{h^2}(1 \quad -2 \quad 1)\begin{pmatrix} y(x-h) \\ y(x) \\ y(x+h) \end{pmatrix}
```

---

### 6.3 이차 차분 행렬 K, T, B

**1차원 영역을 고려하자:**

```math
x_0, x_1, x_2, \ldots, x_{N-1}, x_N, x_{N+1}
```

전체 영역을 $`N+1`$개의 겹치지 않는 요소로 균일 간격 $`h`$로 분할한다.

여기서 $`x_0 = 0`$, $`x_{N+1} = 1`$은 경계점이다.

```math
h = \frac{1}{N+1} \implies x_i = ih = \frac{i}{N+1}
```

방정식 $`-\frac{d^2u}{dx^2} = f(x)`$를 $`u(0) = u(1) = 0`$ (경계 조건)으로 **이산화** (discretize)한다.

경계 조건에서: $`u_0 = u_{N+1} = 0`$.

**Q:** $`u_1, u_2, \ldots, u_N`$은 무엇인가?

```math
\left.\frac{d^2u}{dx^2}\right|_{x=x_1} \approx \frac{u_0 - 2u_1 + u_2}{h^2} = -f(x_1)
```

```math
\left.\frac{d^2u}{dx^2}\right|_{x=x_2} \approx \frac{u_1 - 2u_2 + u_3}{h^2} = -f(x_2)
```

```math
\vdots
```

```math
\left.\frac{d^2u}{dx^2}\right|_{x=x_N} \approx \frac{u_{N-1} - 2u_N + u_{N+1}}{h^2} = -f(x_N)
```

행렬 형태:

```math
\frac{1}{h^2}\begin{pmatrix} 2 & -1 & 0 & \cdots & 0 \\ -1 & 2 & -1 & \cdots & 0 \\ 0 & -1 & 2 & \cdots & 0 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & -1 & 2 \end{pmatrix}\begin{pmatrix} u_1 \\ u_2 \\ u_3 \\ \vdots \\ u_{N-1} \\ u_N \end{pmatrix} = \begin{pmatrix} f_1 \\ f_2 \\ f_3 \\ \vdots \\ f_{N-1} \\ f_N \end{pmatrix}
```

여기서 $`f_i := f(x_i)`$.

```math
\boxed{\frac{1}{h^2}K\mathbf{u} = \mathbf{f}}
```

행렬 $`K`$는 **고정-고정 경계 조건** (fixed-fixed BC): $`u_0 = 0`$이고 $`u_{N+1} = 0`$에서 $`-\frac{d^2u}{dx^2}`$의 자연스러운 근사를 제공한다.

---

### 6.4 K의 성질

$`N = 4`$로 놓자:

```math
K_4 = \begin{pmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}
```

**성질:**

**1. $`K`$는 대칭이다.** $`K^T = K`$.

**2. $`K`$는 띠행렬이다.** $`K`$의 모든 비영 원소는 주 대각선 주위의 "띠"에 놓인다. 좁은 띠를 가진 행렬은 **희소** (sparse)하다 (대부분 0). 이것은 **삼중대각 행렬** (tridiagonal matrix)이다.

예시: $`N = 100`$:
- 2의 개수: 100
- $`-1`$의 개수: $`99 + 99 = 198`$
- 비영 원소의 수 $`\approx 300`$. $`10000`$개 원소 중: $`\frac{300}{10000} = 3\%`$.

**3. $`K`$는 상수 대각선을 가진다.** 이는 푸리에 변환, 필터, 합성곱 행렬, 토에플리츠 행렬 (Toeplitz matrix)과 관련된다. $`K`$는 $`(-1, 2, -1)`$ 패턴이 각 행에 나타나므로 **이동 불변** (shift-invariant)이다.

**4. $`K`$는 가역이다.** 소거로 확인할 수 있다. 결과 상삼각 행렬 $`U`$에 0 피벗이 없으면, $`K`$는 가역이다.

**5. 대칭 행렬 $`K_n`$은 양의 정부호 (positive definite)이다:**

```math
\mathbf{x}^T K \mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}
```

**예제:** $`K = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}`$

```math
K\mathbf{x} = \begin{pmatrix} 2x - y \\ -x + 2y \end{pmatrix}
```

```math
\mathbf{x}^T(K\mathbf{x}) = 2x^2 - 2xy + 2y^2 = x^2 + (x-y)^2 + y^2 > 0
```

$`\mathbf{x}^T K \mathbf{x} \geq 0 \; \forall \; \mathbf{x} \neq \mathbf{0}`$일 때, $`K`$는 **양의 준정부호** (positive semi-definite)라 한다.

**피벗:**
- 가역 행렬은 $`n`$개의 비영 피벗을 가진다.
- 양의 정부호 대칭 행렬은 $`n`$개의 **양의** 피벗을 가진다.
- 양의 준정부호 대칭 행렬은 $`n`$개의 **비음의** 피벗을 가진다.

---

### 6.5 자유-고정 행렬 T

```math
-\frac{d^2u}{dx^2} = f(x), \quad \text{여기서 } \frac{du}{dx} = 0 \text{ (} x = 0 \text{에서)}, \quad u(1) = 0
```

$`x = 0`$에서의 노이만 경계 조건 (Neumann BC): $`\frac{u_1 - u_0}{h} = 0 \implies \frac{1}{h^2}(u_0 - u_1) = 0`$. 여기서 $`u_0`$은 **미지수** 이다.

이것은 행렬 $`T`$를 준다:

```math
\frac{1}{h^2}\begin{pmatrix} 1 & -1 & 0 & 0 & \cdots & 0 \\ -1 & 2 & -1 & 0 & \cdots & 0 \\ 0 & -1 & 2 & -1 & \cdots & 0 \\ \vdots & & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & -1 & 2 \end{pmatrix}\begin{pmatrix} u_0 \\ u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N \end{pmatrix} = \begin{pmatrix} 0 \\ f_1 \\ f_2 \\ \vdots \\ f_{N-1} \\ f_N \end{pmatrix}
```

$`N = 4`$일 때:

```math
T_4 = \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}
```

**$`T_4`$의 소거:**

```math
\xrightarrow{R_2 + R_1} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix} \xrightarrow{R_3' + R_2'} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix} \xrightarrow{R_4'' + R_3''} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 1 \end{pmatrix} = U = L^T
```

```math
\therefore T_4 = LL^T
```

```math
L = \begin{pmatrix} 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}
```

**$`L^{-1}`$은 무엇인가?** 가우스-조르단 방법 사용: $`(L|I) \Rightarrow (I|L^{-1})`$

```math
\begin{pmatrix} 1 & 0 & 0 & 0 & | & 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 & | & 0 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 & | & 0 & 0 & 1 & 0 \\ 0 & 0 & -1 & 1 & | & 0 & 0 & 0 & 1 \end{pmatrix} \implies \begin{pmatrix} 1 & 0 & 0 & 0 & | & 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & | & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 & | & 1 & 1 & 1 & 0 \\ 0 & 0 & 0 & 1 & | & 1 & 1 & 1 & 1 \end{pmatrix}
```

```math
L^{-1} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 1 & 1 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix}
```

```math
T_4^{-1} = (LL^T)^{-1} = L^{-T}L^{-1} = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 1 & 1 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix} = \begin{pmatrix} 4 & 3 & 2 & 1 \\ 3 & 3 & 2 & 1 \\ 2 & 2 & 2 & 1 \\ 1 & 1 & 1 & 1 \end{pmatrix}
```

---

### 6.6 자유-자유 행렬 B

**자유-자유 행렬 $`B`$는 특이 (singular)하다.**

- 가역이 아니다
- 영이 아닌 $`\mathbf{x}`$가 존재하여 $`B\mathbf{x} = \mathbf{0}`$

방정식 $`-\frac{d^2u}{dx^2} = f(x)`$에 대해:

```math
\frac{du}{dx} = 0 \text{ (} x = 0 \text{에서)}, \quad \frac{du}{dx} = 0 \text{ (} x = 1 \text{에서)}
```

```math
\frac{u_1 - u_0}{h} = 0 \implies \frac{1}{h^2}(u_0 - u_1) = 0 \quad (u_0 \text{는 미지수})
```

```math
\frac{u_{N+1} - u_N}{h} = 0 \implies \frac{1}{h^2}(-u_N + u_{N+1}) = 0 \quad (u_{N+1} \text{는 미지수})
```

```math
\frac{1}{h^2}\begin{pmatrix} 1 & -1 & 0 & 0 & \cdots & 0 & 0 \\ -1 & 2 & -1 & 0 & \cdots & 0 & 0 \\ 0 & -1 & 2 & -1 & \cdots & 0 & 0 \\ \vdots & & & \ddots & & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 & 0 \\ 0 & 0 & \cdots & 0 & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} u_0 \\ u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N \\ u_{N+1} \end{pmatrix} = \begin{pmatrix} 0 \\ f_1 \\ f_2 \\ \vdots \\ f_{N-1} \\ f_N \\ 0 \end{pmatrix}
```

**$`B_3`$을 고려하자:**

```math
B_3 = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}
```

```math
B_3 \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```

$`\mathbf{x} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$은 $`B_3`$의 **영공간** (null space)에 속한다.

$`\mathbf{x} = c\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$ 벡터들은 $`B_3`$의 영공간에 속한다.

$`\Rightarrow`$ $`B_3`$은 **가역이 될 수 없다.**

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:-----|:-----------|
| $`A\mathbf{x} = \mathbf{b}`$ | $`n`$개의 미지수를 가진 $`n`$개의 연립일차방정식 |
| 소거법 (Elimination) | 피벗 아래에 0을 만들기 위해 행 $`j`$의 $`l_{ij}`$배를 행 $`i`$에서 뺌 |
| 후진 대입 (Back substitution) | $`U\mathbf{x} = \mathbf{c}`$를 마지막 행부터 위로 풂 |
| 세 가지 경우 | 유일한 해 ($`\text{rank} = n`$), 해 없음, 무한히 많은 해 |
| 확대 행렬 (Augmented matrix) | $`(A \mid \mathbf{b})`$는 소거 중 양변을 추적 |
| 동차 연립방정식 (Homogeneous system) | $`A\mathbf{x} = \mathbf{0}`$; $`\text{rank}(A) < n`$일 때 비자명 해 존재 |
| 피벗 (Pivots) | $`U`$의 대각 원소; 가역성을 위해 영이 아니어야 함 |
| 행 교환 (Row exchange) | 피벗 위치에 0이 나타날 때 행을 교환 |
| 소거 행렬 $`E_{ij}`$ (Elimination matrix) | $`(i,j)`$ 위치에 $`-l_{ij}`$를 가진 단위 행렬 |
| $`EA = U`$ | 모든 소거 행렬의 곱이 $`A`$를 상삼각 $`U`$로 변환 |
| $`A = LU`$ | $`L = E^{-1}`$은 승수 $`l_{ij}`$가 올바른 위치에 있는 하삼각 행렬 |
| 역행렬 $`A^{-1}`$ (Inverse) | $`A`$가 $`n`$개의 독립 열을 가질 때 존재 ($`\text{rank}(A) = n`$) |
| 역행렬의 유일성 | 좌역행렬은 우역행렬과 같음; $`BA = I`$이고 $`AC = I \implies B = C`$ |
| $`(AB)^{-1} = B^{-1}A^{-1}`$ | 역행렬은 역순으로 나옴 |
| $`2 \times 2`$ 역행렬 | $`A^{-1} = \frac{1}{ad-bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}`$; $`\det(A) \neq 0`$ 필요 |
| 가우스-조르단 (Gauss-Jordan) | $`(A \mid I) \Rightarrow (I \mid A^{-1})`$로 역행렬을 명시적으로 구함 |
| 소거의 비용 | $`A \to U`$: $`\frac{1}{3}n^3`$; $`\mathbf{b} \to \mathbf{c} \to \mathbf{x}`$: $`n^2`$ |
| $`A = LU`$의 두 번째 증명 | 열 곱하기 행: $`A = \sum \mathbf{l}_k \mathbf{u}_k`$ |
| 행 교환 없는 $`A = LU`$ | 모든 좌상단 $`k \times k`$ 부분행렬이 가역이어야 함 |
| 치환 행렬 $`P`$ (Permutation matrix) | $`I`$의 행을 다른 순서로; $`P^{-1} = P^T`$; $`n!`$개 가능 |
| $`PA = LU`$ | 행 교환이 $`P`$에 기록된 일반 분해 |
| 부분 피벗팅 (Partial pivoting) | 가장 큰 원소를 피벗으로 만들기 위한 행 교환; 반올림 오차 감소 |
| $`PAQ`$ | 행 치환 $`P`$, 열 치환 $`Q`$; $`C(A) = C(PA)`$ |
| 전치 $`A^T`$ (Transpose) | $`(A^T)_{ij} = A_{ji}`$; $`A^T`$의 열은 $`A`$의 행 |
| 전치 규칙 | $`(A+B)^T = A^T + B^T`$; $`(AB)^T = B^T A^T`$; $`(A^{-1})^T = (A^T)^{-1}`$ |
| 내적 (Inner product) | $`\mathbf{x} \cdot \mathbf{y} = \mathbf{x}^T\mathbf{y}`$; $`(A\mathbf{x}) \cdot \mathbf{y} = \mathbf{x} \cdot (A^T\mathbf{y})`$ |
| 외적 (Outer product) | $`\mathbf{x}\mathbf{y}^T`$는 랭크 1인 $`n \times n`$ 행렬 |
| 대칭 행렬 $`S`$ (Symmetric matrix) | $`S^T = S`$; $`S^{-1}`$도 대칭 |
| $`A^T A`$와 $`AA^T`$ | 항상 대칭; $`A`$의 열이 독립일 때 $`A^T A`$가 가역 |
| $`S = LDL^T`$ | 대칭 분해; $`LU`$가 보여주지 못하는 대칭성을 드러냄 |
| 부분적분과 전치 (IBP and transpose) | $`A = d/dt \implies A^T = -d/dt`$ (반대칭); $`(Ax, y) = (x, A^T y)`$ |
| 전진 차분 (Forward difference) | $`\frac{dy}{dx} \approx \frac{y(x+h)-y(x)}{h}`$; 1차 정확도 $`O(h)`$ |
| 중심 차분 (Centered difference) | $`\frac{dy}{dx} \approx \frac{y(x+h)-y(x-h)}{2h}`$; 2차 정확도 $`O(h^2)`$ |
| 이차 차분 (Second difference) | $`\frac{d^2y}{dx^2} \approx \frac{y(x+h)-2y(x)+y(x-h)}{h^2}`$; 2차 정확도 |
| 행렬 $`K`$ (고정-고정) | 삼중대각 $`(-1, 2, -1)`$; 대칭, 띠, 가역, 양의 정부호 |
| 행렬 $`T`$ (자유-고정) | 첫 행이 $`(1, -1, 0, \ldots)`$; $`T = LL^T`$; 가역 |
| 행렬 $`B`$ (자유-자유) | 특이; $`B\mathbf{1} = \mathbf{0}`$; 상수 벡터가 영공간에 속함 |
| 양의 정부호 (Positive definite) | $`\mathbf{x}^T K\mathbf{x} > 0`$, 모든 $`\mathbf{x} \neq \mathbf{0}`$에 대해; 모든 피벗이 양수 |
| 양의 준정부호 (Positive semi-definite) | $`\mathbf{x}^T K\mathbf{x} \geq 0`$, 모든 $`\mathbf{x} \neq \mathbf{0}`$에 대해; 모든 피벗이 비음수 |

---
