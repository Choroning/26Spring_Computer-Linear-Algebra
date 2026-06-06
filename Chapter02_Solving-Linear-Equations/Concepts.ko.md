# 제2장 강의 — 연립일차방정식 풀기

> **최종 수정일:** 2026-06-06
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 2

> **선수 지식**: [선형대수학] 벡터와 선형 결합 (제1장).
>
> **학습 목표**:
> 1. 가우스 소거법을 적용하여 연립방정식을 풀 수 있다
> 2. 소거법을 행렬 분해(A = LU)로 표현할 수 있다
> 3. 시스템의 해가 없는 경우, 유일한 경우, 무한히 많은 경우를 판별할 수 있다
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 2

> **선수 지식**: [선형대수학] 벡터와 선형 결합 (제1장).
>
> **학습 목표**:
> 1. 가우스 소거법을 적용하여 연립방정식을 풀 수 있다
> 2. 소거법을 행렬 분해(A = LU)로 표현할 수 있다
> 3. 시스템의 해가 없는 경우, 유일한 경우, 무한히 많은 경우를 판별할 수 있다

> **📑 이 문서는 2개 파트로 나뉘어 있습니다.**
>
> **Part 1** · [Part 2](Concepts-2.ko.md)

---

<br>

## 목차

- [1. 서론](#1-서론)
- [2. 소거법과 후진 대입 (2.1)](#2-소거법과-후진-대입-21)
  - [2.1 소거 과정](#21-소거-과정)
  - [2.2 해의 세 가지 경우](#22-해의-세-가지-경우)
  - [2.3 2x2 연립방정식 예제](#23-2x2-연립방정식-예제)
  - [2.4 동차 연립방정식](#24-동차-연립방정식)
  - [2.5 후진 대입 예제](#25-후진-대입-예제)
  - [2.6 각 열에 대한 소거](#26-각-열에-대한-소거)
  - [2.7 소거법의 실패 가능성](#27-소거법의-실패-가능성)
  - [2.8 종속 열과 독립 열](#28-종속-열과-독립-열)
  - [2.9 행 관점과 열 관점](#29-행-관점과-열-관점)
- [3. 소거 행렬과 역행렬 (2.2)](#3-소거-행렬과-역행렬-22)
  - [3.1 소거와 치환의 예제](#31-소거와-치환의-예제)
  - [3.2 소거 행렬과 A = LU](#32-소거-행렬과-a--lu)
  - [3.3 역행렬에 관한 사실들](#33-역행렬에-관한-사실들)
  - [3.4 곱 AB의 역행렬](#34-곱-ab의-역행렬)
  - [3.5 소거 행렬의 역행렬](#35-소거-행렬의-역행렬)
  - [3.6 L은 E의 역행렬](#36-l은-e의-역행렬)
- [4. 행렬 연산과 A = LU (2.3)](#4-행렬-연산과-a--lu-23)
  - [4.1 핵심 사실](#41-핵심-사실)
  - [4.2 역행렬의 명시적 계산](#42-역행렬의-명시적-계산)
  - [4.3 가우스-조르단 소거법](#43-가우스-조르단-소거법)
  - [4.4 소거의 비용](#44-소거의-비용)
  - [4.5 위대한 분해 A = LU](#45-위대한-분해-a--lu)
  - [4.6 A = LU의 두 번째 증명](#46-a--lu의-두-번째-증명)
  - [4.7 행 교환 없는 소거](#47-행-교환-없는-소거)
- [5. 치환과 전치 (2.4)](Concepts-2.ko.md#5-치환과-전치-24)
  - [5.1 치환 행렬](Concepts-2.ko.md#51-치환-행렬)
  - [5.2 치환 행렬의 성질](Concepts-2.ko.md#52-치환-행렬의-성질)
  - [5.3 PA = LU 분해](Concepts-2.ko.md#53-pa--lu-분해)
  - [5.4 부분 피벗팅](Concepts-2.ko.md#54-부분-피벗팅)
  - [5.5 PAQ: 행과 열 치환](Concepts-2.ko.md#55-paq-행과-열-치환)
  - [5.6 A의 전치](Concepts-2.ko.md#56-a의-전치)
  - [5.7 내적과 전치](Concepts-2.ko.md#57-내적과-전치)
  - [5.8 대칭 행렬](Concepts-2.ko.md#58-대칭-행렬)
  - [5.9 대칭 곱과 LDL^T](Concepts-2.ko.md#59-대칭-곱과-ldlt)
- [6. 도함수와 유한 차분 행렬 (2.5)](Concepts-2.ko.md#6-도함수와-유한-차분-행렬-25)
  - [6.1 테일러 급수와 근사](Concepts-2.ko.md#61-테일러-급수와-근사)
  - [6.2 차분으로부터의 도함수](Concepts-2.ko.md#62-차분으로부터의-도함수)
  - [6.3 이차 차분 행렬 K, T, B](Concepts-2.ko.md#63-이차-차분-행렬-k-t-b)
  - [6.4 K의 성질](Concepts-2.ko.md#64-k의-성질)
  - [6.5 자유-고정 행렬 T](Concepts-2.ko.md#65-자유-고정-행렬-t)
  - [6.6 자유-자유 행렬 B](Concepts-2.ko.md#66-자유-자유-행렬-b)
- [요약](Concepts-2.ko.md#요약)

---

<br>

## 1. 서론

제2장은 연립일차방정식을 푸는 것에 집중한다:

```math
A\mathbf{x} = \mathbf{b}
```

**다루는 절:**

1. **2.1** 소거법과 후진 대입 (Elimination and Back Substitution)
2. **2.2** 소거 행렬과 역행렬 (Elimination Matrices and Inverse Matrix)
3. **2.3** 행렬 연산과 $`A = LU`$ (Matrix Computations and $`A = LU`$)
4. **2.4** 치환과 전치 (Permutations and Transposes)
5. **2.5** 도함수와 유한 차분 행렬 (Derivatives and Finite Difference Matrices)

**정사각 행렬** $`A \in \mathbb{R}^{n \times n}`$에 집중한다.

$`A\mathbf{x} = \mathbf{b}`$는 $`n`$개의 방정식을 제공한다.

**예제 (2x2 연립방정식):**

```math
a_{11}x_1 + a_{12}x_2 = b_1
```
```math
a_{21}x_1 + a_{22}x_2 = b_2
```

행렬 형태:

```math
\begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}
```

**일반 형태** ($`n`$개의 방정식, $`n`$개의 미지수):

```math
\begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{n1} & a_{n2} & \cdots & a_{nn} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_n \end{pmatrix}
```

주어진 $`\mathbf{b}`$에 대해 유일한 해 $`\mathbf{x} \in \mathbb{R}^n`$이 존재하면, 역행렬 $`A^{-1}`$이 존재하여 다음을 만족한다:

```math
A^{-1}A = AA^{-1} = I
```

이 경우 해는 다음과 같다:

```math
A\mathbf{x} = \mathbf{b} \implies A^{-1}A\mathbf{x} = A^{-1}\mathbf{b} \implies \boxed{\mathbf{x} = A^{-1}\mathbf{b}}
```

이 장의 목표는 $`A^{-1}`$을 명시적으로 계산하지 **않고** 해 $`\mathbf{x}`$를 구하는 것이다. 다음 두 가지 방법을 사용한다:

1. **소거법** (Elimination)
2. **후진 대입** (Back substitution)

**과정 개요:**

```math
\begin{pmatrix} A & | & \mathbf{b} \end{pmatrix} \xrightarrow{\text{elimination}} \begin{pmatrix} U & | & \mathbf{c} \end{pmatrix} \xrightarrow{\text{back substitution}} \mathbf{x}
```

여기서 $`U`$는 **상삼각 행렬** (upper triangular matrix)이다.

오른쪽에서:

```math
A\mathbf{x} = \mathbf{b} \implies U\mathbf{x} = \mathbf{c} \implies \mathbf{x} = U^{-1}\mathbf{c} = A^{-1}\mathbf{b}
```

---

<br>

## 2. 소거법과 후진 대입 (2.1)

### 2.1 소거 과정

**1단계:** 소거법은 행 $`j`$의 $`l_{ij}`$배를 행 $`i`$에서 빼서, 행 $`i`$에 0을 만든다.

```math
\text{row } i \leftarrow \text{row } i - l_{ij} \cdot \text{row } j
```

**예제:**

```math
\begin{pmatrix} 2 & 3 \\ 4 & 2 \end{pmatrix} \xrightarrow{R_2 - 2 \cdot R_1} \begin{pmatrix} 2 & 3 \\ 0 & -4 \end{pmatrix}
```

**2단계:** $`A\mathbf{x} = \mathbf{b}`$는 다음 중 하나가 된다:
- $`U\mathbf{x} = \mathbf{c}`$ (유일한 해)
- 해 없음
- 무한히 많은 해

**3단계:** $`U\mathbf{x} = \mathbf{c}`$는 **후진 대입** 으로 풀리며, $`U`$는 상삼각 행렬이다.

---

### 2.2 해의 세 가지 경우

$`A\mathbf{x} = \mathbf{b}`$를 고려하자. 여기서 $`A \in \mathbb{R}^{n \times n}`$, $`\mathbf{x}, \mathbf{b} \in \mathbb{R}^{n \times 1}`$이다.

**세 가지 경우** 가 있다:

**경우 1: 유일한 해** — $`\exists! \; \mathbf{x}`$ s.t. $`A\mathbf{x} = \mathbf{b}`$

- $`A`$는 **독립인 열** 을 가진다
- $`A\mathbf{x} = \mathbf{0}`$의 유일한 해는 $`\mathbf{x} = \mathbf{0}`$이다
- $`A`$는 역행렬 $`A^{-1}`$을 가진다

**경우 2: 해 없음** — $`A\mathbf{x} = \mathbf{b}`$에 대해

- $`\mathbf{b}`$가 $`A`$의 열공간(column space)에 **속하지 않는다**
- $`\mathbf{b}`$가 $`A`$의 열들의 선형결합이 아니다

**경우 3: 무한히 많은 해** — $`A\mathbf{x} = \mathbf{b}`$에 대해

- $`A`$의 열들이 **독립이 아니다**
- $`\mathbf{b}`$가 $`A`$의 열공간에 속한다
- $`n > \text{rank}(A) = \text{rank}(A|\mathbf{b})`$

---

### 2.3 2x2 연립방정식 예제

**예제 1: 유일한 해**

```math
x + 2y = 1, \quad 3x + y = -2
```

```math
\begin{pmatrix} 1 & 2 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 1 \\ -2 \end{pmatrix}
```

확대 행렬(augmented matrix):

```math
\begin{pmatrix} 1 & 2 & | & 1 \\ 3 & 1 & | & -2 \end{pmatrix} \xrightarrow{R_2 - 3R_1} \begin{pmatrix} 1 & 2 & | & 1 \\ 0 & -5 & | & -5 \end{pmatrix} \implies \begin{pmatrix} 1 & 2 & | & 1 \\ 0 & 1 & | & 1 \end{pmatrix} \xrightarrow{R_1 - 2R_2} \begin{pmatrix} 1 & 0 & | & -1 \\ 0 & 1 & | & 1 \end{pmatrix}
```

$`\text{rank}(A) = 2`$, $`\text{rank}(A|\mathbf{b}) = 2`$, 미지수의 수 $`= 2`$.

$`\Rightarrow \exists! \; \mathbf{x}`$. 해: $`x = -1, y = 1`$.

**예제 2: 해 없음**

```math
3x + 2y = 3, \quad -6x - 4y = 0
```

```math
\begin{pmatrix} 3 & 2 \\ -6 & -4 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 3 \\ 0 \end{pmatrix}
```

확대 연립방정식:

```math
\begin{pmatrix} 3 & 2 & | & 3 \\ -6 & -4 & | & 0 \end{pmatrix} \xrightarrow{R_2 + 2R_1} \begin{pmatrix} 3 & 2 & | & 3 \\ 0 & 0 & | & 6 \end{pmatrix} \implies \begin{pmatrix} 1 & 2/3 & | & 1 \\ 0 & 0 & | & 1 \end{pmatrix}
```

$`\text{rank}(A) = 1`$, $`\text{rank}(A|\mathbf{b}) = 2`$.

$`\Rightarrow`$ **해가 존재하지 않는다.**

**예제 3: 무한히 많은 해**

```math
3x + 2y = 3, \quad -6x - 4y = -6
```

```math
\begin{pmatrix} 3 & 2 \\ -6 & -4 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 3 \\ -6 \end{pmatrix}
```

```math
\begin{pmatrix} 3 & 2 & | & 3 \\ -6 & -4 & | & -6 \end{pmatrix} \xrightarrow{R_2 + 2R_1} \begin{pmatrix} 3 & 2 & | & 3 \\ 0 & 0 & | & 0 \end{pmatrix} \implies \begin{pmatrix} 1 & 2/3 & | & 1 \\ 0 & 0 & | & 0 \end{pmatrix}
```

$`\text{rank}(A) = 1`$, $`\text{rank}(A|\mathbf{b}) = 1`$, 미지수의 수는 2.

$`\Rightarrow`$ **무한히 많은 해.** 하나의 방정식에 두 개의 미지수, 자유 매개변수 1개.

---

### 2.4 동차 연립방정식

$`\mathbf{b} = \mathbf{0}`$일 때, **동차 연립방정식** (homogeneous system)이 된다:

```math
A\mathbf{x} = \mathbf{0}
```

$`\mathbf{x} = \mathbf{0}`$은 **자명한 해** (trivial solution)이다. 다른 해가 있는가?

$`\text{rank}(A) < n`$일 때 **있다**. $`A\mathbf{x} = \mathbf{0}`$을 만족하는 영이 아닌 벡터 $`\mathbf{x}`$를 $`X`$ (영공간 벡터)로 표기한다.

**핵심 성질:** $`A\mathbf{x} = \mathbf{b}`$의 해 $`\mathbf{x}`$가 하나 존재하면, $`AX = \mathbf{0}`$의 임의의 해를 더할 수 있다:

```math
\mathbf{x} + \alpha X \text{ 는 같은 방정식을 만족한다.}
```

**증명:** $`\alpha \in \mathbb{R}`$에 대해,

```math
A(\mathbf{x} + \alpha X) = A\mathbf{x} + \alpha AX = \mathbf{b} + \mathbf{0} = \mathbf{b}
```

**예제:**

```math
\begin{pmatrix} 2 & 3 \\ 4 & 6 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 6 \\ 12 \end{pmatrix}
```

$`\text{rank}(A) = 1`$, $`\text{rank}(A|\mathbf{b}) = 1 < n = 2`$. $`\Rightarrow`$ 무한히 많은 해.

특수해 선택: $`\mathbf{x} = \begin{pmatrix} 3 \\ 0 \end{pmatrix}`$

$`\begin{pmatrix} 2 & 3 \\ 4 & 6 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$의 비자명 동차해 선택:

```math
X = \begin{pmatrix} 3 \\ -2 \end{pmatrix}
```

모든 벡터 $`\alpha X`$를 특수해 $`\mathbf{x}`$에 더할 수 있다:

```math
\mathbf{x} + \alpha X = \begin{pmatrix} 3 \\ 0 \end{pmatrix} + \alpha \begin{pmatrix} 3 \\ -2 \end{pmatrix}
```

이것은 $`A\mathbf{x} = \mathbf{b}`$에 대한 해의 **직선** (line)을 형성한다.

---

### 2.5 후진 대입 예제

$`A\mathbf{x} = \mathbf{b}`$에 소거법을 적용하면 $`U\mathbf{x} = \mathbf{c}`$가 된다.

**예제:**

```math
\begin{pmatrix} 2 & 3 & 4 \\ 0 & 5 & 6 \\ 0 & 0 & 7 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 19 \\ 17 \\ 14 \end{pmatrix}
```

후진 대입으로 $`\mathbf{x}`$를 구한다:

1. 3번째 행에서: $`7x_3 = 14 \implies x_3 = 2`$
2. 2번째 행에서: $`5x_2 + 6x_3 = 17 \implies 5x_2 + 12 = 17 \implies x_2 = 1`$
3. 1번째 행에서: $`2x_1 + 3x_2 + 4x_3 = 19 \implies 2x_1 + 3 + 8 = 19 \implies x_1 = 4`$

```math
\therefore \mathbf{x} = \begin{pmatrix} 4 \\ 1 \\ 2 \end{pmatrix}
```

**참고 사항:**
- **피벗** (pivots) $`2, 5, 7`$로 나누어야 했다.
- 피벗은 소거 후에 발견된다.
- 0이 피벗이 되는 것은 허용하지 **않는다**. 필요하면 **행 교환** 을 수행한다.
- 독립인 열을 가진 모든 정사각 행렬 $`A`$는 **영이 아닌 피벗** 을 가진 삼각 행렬로 축소될 수 있다.

---

### 2.6 각 열에 대한 소거

**예제:**

```math
A = \begin{pmatrix} 2 & 3 & 4 \\ 4 & 11 & 14 \\ 2 & 8 & 17 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 19 \\ 55 \\ 50 \end{pmatrix}
```

**1단계:** 확대 행렬 $`(A|\mathbf{b})`$를 만들고 $`R_2 - 2R_1`$을 적용:

```math
\begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 4 & 11 & 14 & | & 55 \\ 2 & 8 & 17 & | & 50 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 2 & 8 & 17 & | & 50 \end{pmatrix}
```

이 연산은 다음에 대응된다:

```math
-2R_1 + R_2 + 0 \cdot R_3 = (-2 \quad 1 \quad 0)\begin{pmatrix} R_1 \\ R_2 \\ R_3 \end{pmatrix}
```

소거 행렬(elimination matrix)을 도입한다:

```math
E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

결과는 $`(E_{21}A \mid E_{21}\mathbf{b})`$이다.

**2단계:** $`R_3 - R_1`$ 적용:

```math
\begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 2 & 8 & 17 & | & 50 \end{pmatrix} \xrightarrow{R_3 - R_1} \begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 0 & 5 & 13 & | & 31 \end{pmatrix}
```

```math
E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -1 & 0 & 1 \end{pmatrix}
```

결과: $`(E_{31}E_{21}A \mid E_{31}E_{21}\mathbf{b})`$

**3단계:** $`R_3' - R_2'`$ 적용:

```math
\begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 0 & 5 & 13 & | & 31 \end{pmatrix} \xrightarrow{R_3' - R_2'} \begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 0 & 0 & 7 & | & 14 \end{pmatrix}
```

```math
E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}
```

최종 결과: $`(E_{32}E_{31}E_{21}A \mid E_{32}E_{31}E_{21}\mathbf{b}) = (U \mid \mathbf{c})`$

---

### 2.7 소거법의 실패 가능성

**피벗 위치에 0이 나타날 때** 발생한다.

**예제 (행 교환으로 해결 가능):**

```math
\begin{pmatrix} 2 & 3 & 4 \\ 4 & 6 & 14 \\ 2 & 8 & 17 \end{pmatrix} \rightarrow \begin{pmatrix} 2 & 3 & 4 \\ 0 & 0 & 6 \\ 0 & 5 & 13 \end{pmatrix}
```

치환 행렬 $`P`$를 사용하여 **2행과 3행을 교환**:

```math
\begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 2 & 3 & 4 \\ 0 & 0 & 6 \\ 0 & 5 & 13 \end{pmatrix} = \begin{pmatrix} 2 & 3 & 4 \\ 0 & 5 & 13 \\ 0 & 0 & 6 \end{pmatrix}
```

**예제 (해결 불가 — 특이):**

```math
A^* = \begin{pmatrix} 2 & 3 & 4 \\ 4 & 6 & 14 \\ 2 & 3 & 17 \end{pmatrix} \rightarrow \begin{pmatrix} 2 & 3 & 4 \\ 0 & 0 & 6 \\ 0 & 0 & 13 \end{pmatrix} = U^*
```

두 번째 열에서 피벗을 찾을 수 없다.

- $`A^*`$는 완전 계수(full rank)를 가지지 **않는다**.
- $`A^*`$와 $`U^*`$는 가역이 **아니다**.
- 1열과 2열이 같은 방향이다.
- $`A^* X = \mathbf{0}`$은 영이 아닌 해 $`X`$를 가진다.

---

### 2.8 종속 열과 독립 열

삼각 행렬 $`U`$가 **완전 계수** 를 가지려면 주 대각선에 **0이 없어야** 한다.

완전 계수인 $`A`$에 대해:
- $`U`$의 열들은 독립이다.
- $`U`$의 행들은 독립이다.

$`U`$의 **대각선에 0** 이 있을 때:
- $`U`$는 **특이** (singular) 행렬이다.
- $`U^{-1}`$이 존재하지 않는다.
- $`A^{-1}`$이 존재하지 않는다.
- $`A`$는 **특이** 행렬이다.

---

### 2.9 행 관점과 열 관점

**행 관점 (The Row Picture):**

각 방정식은 직선(2D), 평면(3D), 또는 초평면을 나타낸다. 해는 이들이 교차하는 곳이다.

**예제 1 — 해 없음 (평행선):**

```math
x - 2y = -1, \quad x - 2y = 1 \implies \begin{pmatrix} 1 & -2 \\ 1 & -2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -1 \\ 1 \end{pmatrix}
```

**예제 2 — 무한히 많은 해 (같은 직선):**

```math
x - 2y = 1, \quad x - 2y = 1 \implies \begin{pmatrix} 1 & -2 \\ 1 & -2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}
```

해의 직선.

**예제 3 — 유일한 해 (교차하는 직선):**

```math
x - 2y = 7, \quad x + y = 2 \implies \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 7 \\ 2 \end{pmatrix}
```

교차점에서의 유일한 해.

**열 관점 (The Column Picture):**

```math
A = \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 7 \\ 2 \end{pmatrix}
```

열벡터: $`\mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$, $`\mathbf{a}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}`$

```math
x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 = \mathbf{b}
```

$`\mathbf{b}`$는 $`\mathbf{a}_1`$과 $`\mathbf{a}_2`$의 **선형결합** (linear combination)이다.

$`\Leftrightarrow`$ $`\mathbf{b}`$는 $`A`$의 **열공간** (column space)에 속한다.

---

<br>

## 3. 소거 행렬과 역행렬 (2.2)

### 3.1 소거와 치환의 예제

**예제 1 (치환 불필요):**

```math
A = \begin{pmatrix} 2 & 4 & -2 \\ 4 & 9 & -3 \\ -2 & -3 & 7 \end{pmatrix}
```

```math
E_{21}: R_2 - 2R_1 \implies E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 2 & 4 & -2 \\ 0 & 1 & 1 \\ -2 & -3 & 7 \end{pmatrix}
```

```math
E_{31}: R_3 + R_1 \implies E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 1 & 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 2 & 4 & -2 \\ 0 & 1 & 1 \\ 0 & 1 & 5 \end{pmatrix}
```

```math
E_{32}: R_3' - R_2' \implies E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 2 & 4 & -2 \\ 0 & 1 & 1 \\ 0 & 0 & 4 \end{pmatrix} = U
```

같은 연산을 $`\mathbf{b} = \begin{pmatrix} 2 \\ 8 \\ 10 \end{pmatrix}`$에 적용:

```math
\mathbf{b} = \begin{pmatrix} 2 \\ 8 \\ 10 \end{pmatrix} \xrightarrow{E_{21}} \begin{pmatrix} 2 \\ 4 \\ 10 \end{pmatrix} \xrightarrow{E_{31}} \begin{pmatrix} 2 \\ 4 \\ 12 \end{pmatrix} \xrightarrow{E_{32}} \begin{pmatrix} 2 \\ 4 \\ 8 \end{pmatrix} = \mathbf{c}
```

```math
(A|\mathbf{b}) \xrightarrow{E} (U|\mathbf{c})
```

**예제 2 (치환 필요):**

```math
A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 2 & 3 \\ 0 & 4 & 5 \end{pmatrix}
```

```math
E_{21}: R_2 - 2R_1 \implies \begin{pmatrix} 1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 4 & 5 \end{pmatrix}
```

0 피벗 발생! $`P = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}`$를 사용하여 행 교환:

```math
\begin{pmatrix} 1 & 1 & 1 \\ 0 & 4 & 5 \\ 0 & 0 & 1 \end{pmatrix} = U
```

**$`P`$를 $`A`$에 먼저 적용할 수 있는가?**

```math
PA = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 4 & 5 \\ 2 & 2 & 3 \end{pmatrix}
```

```math
E_{31}: R_3 - 2R_1 \implies E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 1 & 1 & 1 \\ 0 & 4 & 5 \\ 0 & 0 & 1 \end{pmatrix} = U \quad \text{가능하다!}
```

**전체 방정식:**

```math
EPA = U \iff \boxed{PA = E^{-1}U = LU}
```

---

### 3.2 소거 행렬과 A = LU

1. 소거는 $`A`$에 $`E_{21}, E_{31}, \ldots, E_{n1}`$을 곱하고, 그 다음 $`E_{32}, E_{42}, \ldots, E_{n2}`$를 곱하여 $`A`$가 $`EA = U`$가 된다.

2. 역순으로, $`E`$들의 역행렬이 $`U`$에 곱해져 $`A = E^{-1}U`$를 복원한다. 이것이 $`A = LU`$이다.

3. $`A^{-1}A = I`$이고 $`(LU)^{-1} = U^{-1}L^{-1}`$이다. 그러면 $`A\mathbf{x} = \mathbf{b}`$는 다음이 된다:

```math
\mathbf{x} = A^{-1}\mathbf{b} = U^{-1}L^{-1}\mathbf{b}
```

소거의 모든 단계는 행렬로 수행될 수 있다:

```math
E_{n2} \cdots E_{42} E_{32} A = C
```

이 단계들은 행렬로 되돌릴 수 있다:

```math
A = E_{32}^{-1} E_{42}^{-1} \cdots E_{n2}^{-1} C
```

**예제:**

$`A = \begin{pmatrix} 3 & 1 & 0 \\ -3 & 1 & 1 \\ 6 & 8 & 4 \end{pmatrix}`$라 하자.

```math
E_{21}: R_2 + R_1 \implies E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
E_{31}: R_3 - 2R_1 \implies E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 6 & 4 \end{pmatrix}
```

```math
E_{32}: R_3' - 3R_2' \implies E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -3 & 1 \end{pmatrix}
```

```math
\begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 1 \end{pmatrix} = U
```

$`E = E_{32}E_{31}E_{21}`$로 놓으면:

```math
EA = U \implies A = E^{-1}U = LU
```

---

### 3.3 역행렬에 관한 사실들

**정의:** 행렬 $`A`$가 **가역** (invertible)이란 $`A`$를 역변환하는 행렬 $`A^{-1}`$이 존재하는 것이다:

```math
A^{-1}A = AA^{-1} = I
```

$`A \in \mathbb{R}^{n \times n}`$라 하자. $`A`$가 $`n`$개의 독립인 열을 가지면, $`A`$는 가역이다. 이는 $`\text{rank}(A) = n`$을 의미한다.

**참고 1:** 역행렬은 소거가 (행 교환과 함께) $`n`$개의 피벗을 생성할 때에만 존재한다.

소거는 명시적으로 $`A^{-1}`$을 구하지 않고 $`A\mathbf{x} = \mathbf{b}`$를 푼다.

**참고 2:** 행렬 $`A`$는 **두 개의 서로 다른 역행렬을 가질 수 없다**.

$`BA = I`$이고 $`AC = I`$라 가정하자. 그러면 $`B = C`$이다.

*증명:* 결합법칙에 의해,

```math
BAC = B(AC) = (BA)C
```

에서 $`B = C`$가 된다. $`\square`$

좌역행렬 (left inverse) $`B`$와 우역행렬 (right inverse) $`C`$는 같아야 한다.

**참고 3:** $`A`$가 가역이면, $`A\mathbf{x} = \mathbf{b}`$에 대해 $`\exists! \; \mathbf{x}`$ s.t. $`\mathbf{x} = A^{-1}\mathbf{b}`$.

*증명:*

```math
A\mathbf{x} = \mathbf{b} \implies A^{-1}A\mathbf{x} = A^{-1}\mathbf{b} \implies \mathbf{x} = A^{-1}\mathbf{b}. \quad \square
```

**참고 4:** 영이 아닌 벡터 $`\mathbf{x}`$가 존재하여 $`A\mathbf{x} = \mathbf{0}`$이라 하자. 그러면:

- $`A`$는 **종속인 열** 을 가진다
- $`A`$는 역행렬을 가질 수 없다
- 어떤 행렬도 $`\mathbf{0}`$을 $`\mathbf{x}`$로 되돌릴 수 없다

**예제:**

```math
\begin{pmatrix} 1 & 2 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} -2 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```

$`A`$가 가역이면, $`A\mathbf{x} = \mathbf{0}`$은 $`\mathbf{x} = A^{-1}\mathbf{0} = \mathbf{0}`$을 의미한다.

**참고 5:** 정사각 행렬은 열들이 독립일 때에만 가역이다.

**참고 6:** $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$는 $`ad - bc \neq 0`$ ($`A`$의 **행렬식** (determinant))일 때에만 가역이다.

```math
A^{-1} = \frac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}
```

**참고 7:** 삼각 행렬은 **영이 아닌 대각 원소** 를 가지면 역행렬을 가진다.

$`A`$가 대각 원소 $`d_1, d_2, \ldots, d_n`$ (모두 영이 아닌)을 가진 상삼각 행렬이면, $`A^{-1}`$도 대각 원소 $`\frac{1}{d_1}, \frac{1}{d_2}, \ldots, \frac{1}{d_n}`$을 가진 상삼각 행렬이다.

**예제 2:**

```math
A = \begin{pmatrix} 1 & 2 \\ 1 & 2 \end{pmatrix} \xrightarrow{R_2 - R_1} \begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}
```

$`\text{rank}(A) = 1 < 2`$. $`A`$는 1개의 피벗을 가진다. $`\det(A) = 1 \cdot 2 - 2 \cdot 1 = 0`$. $`A`$는 종속인 열을 가진다.

**예제 3:**

```math
A = \begin{pmatrix} 4 & 3 \\ 8 & 6 \end{pmatrix} \quad \text{rank}(A) = 1 < 2 \quad \text{(가역이 아님)}
```

```math
B = \begin{pmatrix} 4 & 3 \\ 8 & 7 \end{pmatrix} \quad \det(B) = 4 \cdot 7 - 3 \cdot 8 = 4 \neq 0 \quad B^{-1} = \frac{1}{4}\begin{pmatrix} 7 & -3 \\ -8 & 4 \end{pmatrix}
```

```math
C = \begin{pmatrix} 6 & 6 \\ 6 & 0 \end{pmatrix} \quad \det(C) = 6 \cdot 0 - 6 \cdot 6 = -36 \neq 0 \quad C^{-1} = -\frac{1}{36}\begin{pmatrix} 0 & -6 \\ -6 & 6 \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 0 & 1 \\ 1 & -1 \end{pmatrix}
```

```math
D = \begin{pmatrix} 6 & 6 \\ 6 & 6 \end{pmatrix} \quad \text{rank}(D) = 1 < 2 \quad \text{(가역이 아님)}
```

```math
S = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_2 - R_1, R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \xrightarrow{R_3' - R_2'} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

3개의 피벗. 가역.

$`E = E_{32}E_{31}E_{21} = S^{-1}`$을 계산:

```math
E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -1 & 0 & 1 \end{pmatrix}, \quad E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}
```

```math
E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ -1 & 0 & 1 \end{pmatrix}
```

```math
E_{32}(E_{31}E_{21}) = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix} = E = S^{-1}
```

```math
T = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}
```

$`\text{Rank}(T) = 2 < 3`$. **가역이 아니다.**

---

### 3.4 곱 AB의 역행렬

영이 아닌 두 값 $`a`$와 $`b`$에 대해, 합 $`(a + b)`$는 가역이 아닐 수 있다.

**예제:** $`a = 3 \implies a^{-1} = 1/3`$, $`b = -3 \implies b^{-1} = -1/3`$, $`a + b = 0 \implies (a+b)^{-1}`$은 존재하지 않는다.

그러나: $`ab = -9 \implies (ab)^{-1} = -1/9 = a^{-1}b^{-1}`$.

**정리:** $`A, B \in \mathbb{R}^{n \times n}`$이 가역이면, $`AB`$의 역행렬은:

```math
\boxed{(AB)^{-1} = B^{-1}A^{-1}}
```

*증명:*

```math
(AB)^{-1}AB = I
```
```math
(AB)^{-1}ABB^{-1} = IB^{-1} = B^{-1}
```
```math
(AB)^{-1}AA^{-1} = B^{-1}A^{-1}
```

```math
\therefore (AB)^{-1} = B^{-1}A^{-1} \quad \square
```

**역행렬은 역순으로 나온다!**

세 행렬에 대해:

```math
(ABC)^{-1} = C^{-1}B^{-1}A^{-1}
```

*증명:*

```math
(ABC)^{-1}ABC = I \implies (ABC)^{-1}ABCC^{-1} = C^{-1} \implies (ABC)^{-1}ABB^{-1} = C^{-1}B^{-1}
```

```math
\implies (ABC)^{-1}AA^{-1} = C^{-1}B^{-1}A^{-1} \implies (ABC)^{-1} = C^{-1}B^{-1}A^{-1}
```

---

### 3.5 소거 행렬의 역행렬

**예제 4:**

```math
E = \begin{pmatrix} 1 & 0 & 0 \\ -5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \quad \text{(1행의 5배를 2행에서 뺌)}
```

$`(-5 \quad 1 \quad 0)\begin{pmatrix} R_1 \\ R_2 \\ R_3 \end{pmatrix} = R_2 - 5R_1 = R_2'`$

되돌리기: $`R_2 = R_2' + 5R_1`$, 즉 $`(5 \quad 1 \quad 0)\begin{pmatrix} R_1 \\ R_2' \\ R_3 \end{pmatrix}`$

```math
E^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \quad \text{(1행의 5배를 2행에 더함)}
```

```math
EE^{-1} = E^{-1}E = I
```

**중요:** 정사각 행렬 $`A, C`$에 대해 $`AC = I`$이면, $`CA = I`$이다.

**예제 5:**

```math
F = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -4 & 1 \end{pmatrix} \quad (R_3' = R_3 - 4R_2)
```

```math
F^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 4 & 1 \end{pmatrix} \quad (R_3 = R_3' + 4R_2)
```

곱 $`FE`$ 계산:

```math
FE = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -4 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ -5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ -5 & 1 & 0 \\ 20 & -4 & 1 \end{pmatrix}
```

"20"은 1행에서 연산의 연쇄로 나온다:

```math
R_3'' = R_3' - 4R_2' = R_3 - 4(R_2 - 5R_1) = R_3 - 4R_2 + 20R_1
```

```math
(FE)^{-1} = E^{-1}F^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 4 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 5 & 1 & 0 \\ 0 & 4 & 1 \end{pmatrix}
```

승수(multiplier) 5와 4가 $`L = (FE)^{-1}`$의 대각선 아래에 제자리에 놓인다.

---

### 3.6 L은 E의 역행렬

**L은 E의 역행렬이다.**

예제 1을 상기하자: $`A = \begin{pmatrix} 3 & 1 & 0 \\ -3 & 1 & 1 \\ 6 & 8 & 4 \end{pmatrix}`$

```math
E = E_{32}E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -3 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
EA = \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 1 \end{pmatrix} = U
```

**E와 L에 대한 일반 공식 (3x3 경우):**

$`l_{32} = 3`$, $`l_{31} = 2`$, $`l_{21} = -1`$로 놓자. 그러면:

```math
E = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -l_{32} & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -l_{31} & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ -l_{21} & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

```math
= \begin{pmatrix} 1 & 0 & 0 \\ -l_{21} & 1 & 0 \\ l_{32}l_{21} - l_{31} & -l_{32} & 1 \end{pmatrix}
```

(3,1) 위치에 교차곱 항 $`l_{32}l_{21} - l_{31}`$이 있음에 주목하라.

```math
E^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ l_{31} & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & l_{32} & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{pmatrix} = L
```

**역행렬 $`E^{-1}`$이 아름답게 된다!** 인수 $`l_{21}, l_{31}, l_{32}`$를 볼 수 있다.

**모든 승수 $`l_{ij}`$가 $`L`$에서 올바른 위치에 나타난다.**

---

<br>

## 4. 행렬 연산과 A = LU (2.3)

### 4.1 핵심 사실

1. $`A`$에서 $`U`$로의 소거 단계는 $`\frac{1}{3}n^3`$번의 곱셈과 뺄셈이 소요된다.
2. 각 우변 $`\mathbf{b}`$는 $`n^2`$만 소요된다:
   - $`U\mathbf{x} = \mathbf{c}`$로 전진
   - 그 다음 $`\mathbf{x}`$에 대한 후진 대입
3. 행 교환 없는 소거는 $`A`$를 $`LU`$로 분해한다.

$`A\mathbf{x} = \mathbf{b}`$의 해는 $`\mathbf{x} = A^{-1}\mathbf{b}`$로 주어진다.

**Q:** $`A^{-1}`$을 명시적으로 알아야 하는가? $`A^{-1}`$을 계산하고 $`A^{-1}\mathbf{b}`$를 곱하는 것은 $`\mathbf{x}`$를 구하는 매우 **느린** 방법이다.

---

### 4.2 역행렬의 명시적 계산

**Q:** $`A^{-1}`$을 명시적으로 어떻게 구하는가?

$`AA^{-1} = I \in \mathbb{R}^{n \times n}`$에서 시작한다.

```math
I = \begin{pmatrix} | & | & & | \\ \hat{e}_1 & \hat{e}_2 & \cdots & \hat{e}_n \\ | & | & & | \end{pmatrix}
```

여기서 $`\hat{e}_1, \hat{e}_2, \ldots, \hat{e}_n`$은 **표준 기저 벡터** (단위 벡터)이다.

$`AA^{-1} = I`$를 다음과 같이 본다:

```math
A\begin{pmatrix} | & | & & | \\ \mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_n \\ | & | & & | \end{pmatrix} = \begin{pmatrix} | & | & & | \\ \hat{e}_1 & \hat{e}_2 & \cdots & \hat{e}_n \\ | & | & & | \end{pmatrix}
```

즉 $`n`$개의 방정식이다:

```math
A\mathbf{x}_1 = \hat{e}_1, \quad A\mathbf{x}_2 = \hat{e}_2, \quad \ldots, \quad A\mathbf{x}_n = \hat{e}_n
```

같은 계수 행렬 $`A`$, 서로 다른 우변 벡터.

---

### 4.3 가우스-조르단 소거법

$`n`$개의 방정식을 **가우스-조르단 소거법** (Gauss-Jordan elimination)을 사용하여 함께 푼다:

```math
(A \mid I) \implies (I \mid A^{-1})
```

**예제:**

```math
\begin{pmatrix} 1 & 0 & 0 & | & 1 & 0 & 0 \\ -1 & 1 & 0 & | & 0 & 1 & 0 \\ 0 & -1 & 1 & | & 0 & 0 & 1 \end{pmatrix} = (A|I)
```

```math
\xrightarrow{R_2 + R_1} \begin{pmatrix} 1 & 0 & 0 & | & 1 & 0 & 0 \\ 0 & 1 & 0 & | & 1 & 1 & 0 \\ 0 & -1 & 1 & | & 0 & 0 & 1 \end{pmatrix}
```

```math
\xrightarrow{R_3 + R_2} \begin{pmatrix} 1 & 0 & 0 & | & 1 & 0 & 0 \\ 0 & 1 & 0 & | & 1 & 1 & 0 \\ 0 & 0 & 1 & | & 1 & 1 & 1 \end{pmatrix} = (I|A^{-1})
```

$`A`$에 대한 소거 단계는 **한 번만** 수행하면 된다!

---

### 4.4 소거의 비용

**$`A`$를 $`U`$로 축소:**

1단계 (1열 소거): $`(n-1)`$개 행, $`n`$개 열 $`\implies (n-1)n`$번의 곱셈과 $`(n-1)n`$번의 뺄셈.

2단계 (2열 소거): $`(n-2)`$개 행, $`(n-1)`$개 열 $`\implies (n-2)(n-1)`$번의 곱셈과 뺄셈.

$`\vdots`$

$`(n-1)`$단계 ($`n-1`$열 소거): 1개 행, 2개 열 $`\implies 1 \cdot 2`$번의 곱셈과 뺄셈.

**총 곱셈 횟수:**

```math
(n-1)n + (n-2)(n-1) + \cdots + 1 \cdot 2 = \sum_{i=1}^{n-1} i(i+1) = \sum_{i=1}^{n-1} i^2 + \sum_{i=1}^{n-1} i
```

```math
= \frac{n(n+1)(2n+1)}{6} - n^2 + \frac{n^2}{2} - \frac{n}{2} + \frac{(n-1)n}{2}
```

```math
= \frac{1}{3}n^3 - \frac{n}{3}
```

$`n \to \infty`$일 때: $`\approx \dfrac{1}{3}n^3`$.

**$`A`$를 $`U`$로 축소하는 데 약 $`\frac{1}{3}n^3`$번의 곱셈과 $`\frac{1}{3}n^3`$번의 뺄셈이 필요하다.**

**$`\mathbf{b}`$를 $`\mathbf{c}`$로 축소:**

$`A`$와 유사하지만 1개 열만:

- 1단계: $`(n-1)`$번의 곱셈과 뺄셈
- 2단계: $`(n-2)`$번의 곱셈과 뺄셈
- ...
- $`(n-1)`$단계: 1번의 곱셈과 뺄셈

```math
\sum_{i=1}^{n-1} i = \frac{(n-1)n}{2} = \frac{n^2}{2} - \frac{n}{2}
```

$`n \to \infty`$일 때: $`\approx \dfrac{n^2}{2}`$.

**$`\mathbf{b}`$를 $`\mathbf{c}`$로 축소하는 데 $`\frac{n^2}{2}`$번의 곱셈과 $`\frac{n^2}{2}`$번의 뺄셈이 필요하다.**

**후진 대입** ($`U\mathbf{x} = \mathbf{c}`$):

```math
x_n = c_n / u_{nn}
```
```math
x_{n-1} = (c_{n-1} - u_{(n-1)n}x_n) / u_{(n-1)(n-1)}
```
```math
\vdots
```
```math
x_1 = (c_1 - u_{12}x_2 - u_{13}x_3 - \cdots - u_{1n}x_n) / u_{11}
```

$`n`$번의 나눗셈, $`1 + 2 + \cdots + (n-1)`$번의 곱셈, $`1 + 2 + \cdots + (n-1)`$번의 뺄셈.

```math
n + \frac{(n+1)n}{2} + \frac{(n-1)n}{2} = n^2
```

**$`\mathbf{b}`$에서 $`\mathbf{c}`$를 거쳐 $`\mathbf{x}`$까지 우변의 총 연산 횟수는 $`n^2`$이다:**
- $`n^2`$번의 곱셈
- $`n^2`$번의 뺄셈

---

### 4.5 위대한 분해 A = LU

하나의 소거 단계 $`E_{ij}`$를 역변환하려면 (행 $`i`$에서 행 $`j`$의 $`l_{ij}`$배를 뺀 것을), 빼는 대신 **더한다**.

```math
E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -l_{31} & 0 & 1 \end{pmatrix} \implies L_{31} = E_{31}^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ l_{31} & 0 & 1 \end{pmatrix}
```

**예제 (재방문):**

$`A = \begin{pmatrix} 3 & 1 & 0 \\ -3 & 1 & 1 \\ 6 & 8 & 4 \end{pmatrix}`$로 놓자.

```math
E_{21}: R_2 + R_1 \implies \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 6 & 8 & 4 \end{pmatrix}
```

```math
E_{31}: R_3 - 2R_1 \implies \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 6 & 4 \end{pmatrix}
```

```math
E_{32}: R_3' - 3R_2' \implies \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 1 \end{pmatrix} = U
```

$`U`$의 1행 = $`A`$의 1행임에 주목하라. $`U`$의 2행은 $`E_{31}`$과 $`E_{32}`$ 적용 후에도 변하지 않는다.

**$`A`$와 $`U`$의 행 사이의 관계:**

$`U`$의 3행 = $`A`$의 3행 $`-`$ 2 ($`U`$의 1행) $`-`$ 3 ($`U`$의 2행)

동치로:

$`A`$의 3행 = $`U`$의 3행 $`+`$ $`l_{31}`$ ($`U`$의 1행) $`+`$ $`l_{32}`$ ($`U`$의 2행)

여기서 $`l_{31} = 2`$이고 $`l_{32} = 3`$이다.

```math
A\text{의 3행} = (l_{31} \quad l_{32} \quad 1) \begin{pmatrix} U_1 \\ U_2 \\ U_3 \end{pmatrix}
```

```math
\implies A = LU
```

---

### 4.6 A = LU의 두 번째 증명

**열 곱하기 행 (columns times rows).**

$`A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}`$을 고려하자.

**1단계:** 1행을 피벗 행으로 취한다. 1행에 $`l_{21}, l_{31}`$을 곱하여 2행, 3행에서 뺀다.

$`l_{21} = a_{21}/a_{11}`$, $`l_{31} = a_{31}/a_{11}`$으로 선택한다.

```math
R_2 - l_{21}R_1, \quad R_3 - l_{31}R_1
```

```math
A' = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ 0 & a'_{22} & a'_{23} \\ 0 & a'_{32} & a'_{33} \end{pmatrix}
```

이 단계는 다음과 같이 볼 수 있다:

```math
A = A' + \begin{pmatrix} - & 0 & - \\ - & l_{21}R_1 & - \\ - & l_{31}R_1 & - \end{pmatrix}
```

뺀 부분은 **랭크 1 행렬** (rank 1 matrix)이다:

```math
\begin{pmatrix} 1 \\ l_{21} \\ l_{31} \end{pmatrix}\begin{pmatrix} a_{11} & a_{12} & a_{13} \end{pmatrix} = \mathbf{l}_1 \mathbf{u}_1 = \mathbf{l}_1 \otimes \mathbf{u}_1
```

```math
= \begin{pmatrix} 1 \cdot a_{11} & 1 \cdot a_{12} & 1 \cdot a_{13} \\ l_{21}a_{11} & l_{21}a_{12} & l_{21}a_{13} \\ l_{31}a_{11} & l_{31}a_{12} & l_{31}a_{13} \end{pmatrix}
```

**2단계:** $`A_2`$ (나머지 부분)의 2행을 피벗 행으로 취한다.

```math
A_2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & a'_{22} & a'_{23} \\ 0 & a'_{32} & a'_{33} \end{pmatrix}
```

$`R_3' - l_{32}R_2'`$ 적용:

```math
A'' = \begin{pmatrix} 0 & 0 & 0 \\ 0 & a'_{22} & a'_{23} \\ 0 & 0 & a''_{33} \end{pmatrix} = A_2 - \begin{pmatrix} 0 \\ 1 \\ l_{32} \end{pmatrix}(R_2')
```

이것은 또 다른 랭크 1 행렬을 생성한다:

```math
\begin{pmatrix} 0 \\ 1 \\ l_{32} \end{pmatrix}\begin{pmatrix} 0 & a'_{22} & a'_{23} \end{pmatrix} = \mathbf{l}_2 \mathbf{u}_2
```

$`\mathbf{u}_2`$는 $`U`$의 2행임에 주목하라.

나머지 부분:

```math
A_3 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\begin{pmatrix} 0 & 0 & a''_{33} \end{pmatrix} = \mathbf{l}_3 \mathbf{u}_3
```

이제 $`A`$를 다음과 같이 표현한다:

```math
A = \mathbf{l}_1 \mathbf{u}_1 + \mathbf{l}_2 \mathbf{u}_2 + \mathbf{l}_3 \mathbf{u}_3
```

```math
= \begin{pmatrix} | & | & | \\ \mathbf{l}_1 & \mathbf{l}_2 & \mathbf{l}_3 \\ | & | & | \end{pmatrix}\begin{pmatrix} - & \mathbf{u}_1 & - \\ - & \mathbf{u}_2 & - \\ - & \mathbf{u}_3 & - \end{pmatrix}
```

```math
= \begin{pmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{pmatrix}\begin{pmatrix} a_{11} & a_{12} & a_{13} \\ 0 & a'_{22} & a'_{23} \\ 0 & 0 & a''_{33} \end{pmatrix} = LU
```

**일반 확장** — $`A \in \mathbb{R}^{n \times n}`$에 대해:

```math
A = \mathbf{l}_1\mathbf{u}_1 + \mathbf{l}_2\mathbf{u}_2 + \cdots + \mathbf{l}_n\mathbf{u}_n
```

```math
= \begin{pmatrix} 1 & 0 & 0 & \cdots & 0 \\ l_{21} & 1 & 0 & \cdots & 0 \\ l_{31} & l_{32} & 1 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ l_{n1} & l_{n2} & l_{n3} & \cdots & 1 \end{pmatrix}\begin{pmatrix} a_{11} & a_{12} & a_{13} & \cdots & a_{1n} \\ 0 & a'_{22} & a'_{23} & \cdots & a'_{2n} \\ 0 & 0 & a''_{33} & \cdots & a''_{3n} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & a_{nn}^{(n-1)} \end{pmatrix} = LU
```

$`\mathbf{l}_k`$는 $`(k-1)`$개의 0으로 시작한다. $`\mathbf{u}_k`$는 $`(k-1)`$개의 0으로 시작한다.

---

### 4.7 행 교환 없는 소거

**Q:** $`A = LU`$가 **행 교환 없이** **피벗에 0 없이** 가능한 것은 언제인가?

**A:** $`A`$의 모든 좌상단 $`k \times k`$ 부분행렬이 가역이어야 한다.

3x3 행렬 $`A`$에 대해:
- $`A_1 = (a_{11})`$: $`A_1 = L_1 U_1`$ (1x1 부분행렬이 가역이어야 함)
- $`A_2 = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}`$: $`A_2 = L_2 U_2`$ (2x2 부분행렬이 가역이어야 함)
- $`A_3 = A`$ (전체 행렬): $`A_3 = L_3 U_3`$

---

<br>
