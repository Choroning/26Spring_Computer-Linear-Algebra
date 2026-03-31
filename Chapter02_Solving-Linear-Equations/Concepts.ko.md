# 제2장 강의 — 연립일차방정식 풀기

> **Last Updated:** 2026-03-31

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

## 1. 서론

제2장은 연립일차방정식을 푸는 것에 집중한다:

$$A\mathbf{x} = \mathbf{b}$$

**다루는 절:**

1. **2.1** 소거법과 후진 대입 (Elimination and Back Substitution)
2. **2.2** 소거 행렬과 역행렬 (Elimination Matrices and Inverse Matrix)
3. **2.3** 행렬 연산과 $A = LU$ (Matrix Computations and $A = LU$)
4. **2.4** 치환과 전치 (Permutations and Transposes)
5. **2.5** 도함수와 유한 차분 행렬 (Derivatives and Finite Difference Matrices)

**정사각 행렬** $A \in \mathbb{R}^{n \times n}$에 집중한다.

$A\mathbf{x} = \mathbf{b}$는 $n$개의 방정식을 제공한다.

**예제 (2x2 연립방정식):**

$$a_{11}x_1 + a_{12}x_2 = b_1$$
$$a_{21}x_1 + a_{22}x_2 = b_2$$

행렬 형태:

$$\begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$$

**일반 형태** ($n$개의 방정식, $n$개의 미지수):

$$\begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{n1} & a_{n2} & \cdots & a_{nn} \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \\ \vdots \\ b_n \end{pmatrix}$$

주어진 $\mathbf{b}$에 대해 유일한 해 $\mathbf{x} \in \mathbb{R}^n$이 존재하면, 역행렬 $A^{-1}$이 존재하여 다음을 만족한다:

$$A^{-1}A = AA^{-1} = I$$

이 경우 해는 다음과 같다:

$$A\mathbf{x} = \mathbf{b} \implies A^{-1}A\mathbf{x} = A^{-1}\mathbf{b} \implies \boxed{\mathbf{x} = A^{-1}\mathbf{b}}$$

이 장의 목표는 $A^{-1}$을 명시적으로 계산하지 **않고** 해 $\mathbf{x}$를 구하는 것이다. 다음 두 가지 방법을 사용한다:

1. **소거법** (Elimination)
2. **후진 대입** (Back substitution)

**과정 개요:**

$$\begin{pmatrix} A & | & \mathbf{b} \end{pmatrix} \xrightarrow{\text{elimination}} \begin{pmatrix} U & | & \mathbf{c} \end{pmatrix} \xrightarrow{\text{back substitution}} \mathbf{x}$$

여기서 $U$는 **상삼각 행렬** (upper triangular matrix)이다.

오른쪽에서:

$$A\mathbf{x} = \mathbf{b} \implies U\mathbf{x} = \mathbf{c} \implies \mathbf{x} = U^{-1}\mathbf{c} = A^{-1}\mathbf{b}$$

---

<br>

## 2. 소거법과 후진 대입 (2.1)

### 2.1 소거 과정

**1단계:** 소거법은 행 $j$의 $l_{ij}$배를 행 $i$에서 빼서, 행 $i$에 0을 만든다.

$$\text{row } i \leftarrow \text{row } i - l_{ij} \cdot \text{row } j$$

**예제:**

$$\begin{pmatrix} 2 & 3 \\ 4 & 2 \end{pmatrix} \xrightarrow{R_2 - 2 \cdot R_1} \begin{pmatrix} 2 & 3 \\ 0 & -4 \end{pmatrix}$$

**2단계:** $A\mathbf{x} = \mathbf{b}$는 다음 중 하나가 된다:
- $U\mathbf{x} = \mathbf{c}$ (유일한 해)
- 해 없음
- 무한히 많은 해

**3단계:** $U\mathbf{x} = \mathbf{c}$는 **후진 대입**으로 풀리며, $U$는 상삼각 행렬이다.

---

### 2.2 해의 세 가지 경우

$A\mathbf{x} = \mathbf{b}$를 고려하자. 여기서 $A \in \mathbb{R}^{n \times n}$, $\mathbf{x}, \mathbf{b} \in \mathbb{R}^{n \times 1}$이다.

**세 가지 경우**가 있다:

**경우 1: 유일한 해** — $\exists! \; \mathbf{x}$ s.t. $A\mathbf{x} = \mathbf{b}$

- $A$는 **독립인 열**을 가진다
- $A\mathbf{x} = \mathbf{0}$의 유일한 해는 $\mathbf{x} = \mathbf{0}$이다
- $A$는 역행렬 $A^{-1}$을 가진다

**경우 2: 해 없음** — $A\mathbf{x} = \mathbf{b}$에 대해

- $\mathbf{b}$가 $A$의 열공간(column space)에 **속하지 않는다**
- $\mathbf{b}$가 $A$의 열들의 일차결합이 아니다

**경우 3: 무한히 많은 해** — $A\mathbf{x} = \mathbf{b}$에 대해

- $A$의 열들이 **독립이 아니다**
- $\mathbf{b}$가 $A$의 열공간에 속한다
- $n > \text{rank}(A) = \text{rank}(A|\mathbf{b})$

---

### 2.3 2x2 연립방정식 예제

**예제 1: 유일한 해**

$$x + 2y = 1, \quad 3x + y = -2$$

$$\begin{pmatrix} 1 & 2 \\ 3 & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 1 \\ -2 \end{pmatrix}$$

확대 행렬(augmented matrix):

$$\begin{pmatrix} 1 & 2 & | & 1 \\ 3 & 1 & | & -2 \end{pmatrix} \xrightarrow{R_2 - 3R_1} \begin{pmatrix} 1 & 2 & | & 1 \\ 0 & -5 & | & -5 \end{pmatrix} \implies \begin{pmatrix} 1 & 2 & | & 1 \\ 0 & 1 & | & 1 \end{pmatrix} \xrightarrow{R_1 - 2R_2} \begin{pmatrix} 1 & 0 & | & -1 \\ 0 & 1 & | & 1 \end{pmatrix}$$

$\text{rank}(A) = 2$, $\text{rank}(A|\mathbf{b}) = 2$, 미지수의 수 $= 2$.

$\Rightarrow \exists! \; \mathbf{x}$. 해: $x = -1, y = 1$.

**예제 2: 해 없음**

$$3x + 2y = 3, \quad -6x - 4y = 0$$

$$\begin{pmatrix} 3 & 2 \\ -6 & -4 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 3 \\ 0 \end{pmatrix}$$

확대 연립방정식:

$$\begin{pmatrix} 3 & 2 & | & 3 \\ -6 & -4 & | & 0 \end{pmatrix} \xrightarrow{R_2 + 2R_1} \begin{pmatrix} 3 & 2 & | & 3 \\ 0 & 0 & | & 6 \end{pmatrix} \implies \begin{pmatrix} 1 & 2/3 & | & 1 \\ 0 & 0 & | & 1 \end{pmatrix}$$

$\text{rank}(A) = 1$, $\text{rank}(A|\mathbf{b}) = 2$.

$\Rightarrow$ **해가 존재하지 않는다.**

**예제 3: 무한히 많은 해**

$$3x + 2y = 3, \quad -6x - 4y = -6$$

$$\begin{pmatrix} 3 & 2 \\ -6 & -4 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 3 \\ -6 \end{pmatrix}$$

$$\begin{pmatrix} 3 & 2 & | & 3 \\ -6 & -4 & | & -6 \end{pmatrix} \xrightarrow{R_2 + 2R_1} \begin{pmatrix} 3 & 2 & | & 3 \\ 0 & 0 & | & 0 \end{pmatrix} \implies \begin{pmatrix} 1 & 2/3 & | & 1 \\ 0 & 0 & | & 0 \end{pmatrix}$$

$\text{rank}(A) = 1$, $\text{rank}(A|\mathbf{b}) = 1$, 미지수의 수는 2.

$\Rightarrow$ **무한히 많은 해.** 하나의 방정식에 두 개의 미지수, 자유 매개변수 1개.

---

### 2.4 동차 연립방정식

$\mathbf{b} = \mathbf{0}$일 때, **동차 연립방정식** (homogeneous system)이 된다:

$$A\mathbf{x} = \mathbf{0}$$

$\mathbf{x} = \mathbf{0}$은 **자명한 해** (trivial solution)이다. 다른 해가 있는가?

$\text{rank}(A) < n$일 때 **있다**. $A\mathbf{x} = \mathbf{0}$을 만족하는 영이 아닌 벡터 $\mathbf{x}$를 $X$ (영공간 벡터)로 표기한다.

**핵심 성질:** $A\mathbf{x} = \mathbf{b}$의 해 $\mathbf{x}$가 하나 존재하면, $AX = \mathbf{0}$의 임의의 해를 더할 수 있다:

$$\mathbf{x} + \alpha X \text{ 는 같은 방정식을 만족한다.}$$

**증명:** $\alpha \in \mathbb{R}$에 대해,

$$A(\mathbf{x} + \alpha X) = A\mathbf{x} + \alpha AX = \mathbf{b} + \mathbf{0} = \mathbf{b}$$

**예제:**

$$\begin{pmatrix} 2 & 3 \\ 4 & 6 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 6 \\ 12 \end{pmatrix}$$

$\text{rank}(A) = 1$, $\text{rank}(A|\mathbf{b}) = 1 < n = 2$. $\Rightarrow$ 무한히 많은 해.

특수해 선택: $\mathbf{x} = \begin{pmatrix} 3 \\ 0 \end{pmatrix}$

$\begin{pmatrix} 2 & 3 \\ 4 & 6 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$의 비자명 동차해 선택:

$$X = \begin{pmatrix} 3 \\ -2 \end{pmatrix}$$

모든 벡터 $\alpha X$를 특수해 $\mathbf{x}$에 더할 수 있다:

$$\mathbf{x} + \alpha X = \begin{pmatrix} 3 \\ 0 \end{pmatrix} + \alpha \begin{pmatrix} 3 \\ -2 \end{pmatrix}$$

이것은 $A\mathbf{x} = \mathbf{b}$에 대한 해의 **직선** (line)을 형성한다.

---

### 2.5 후진 대입 예제

$A\mathbf{x} = \mathbf{b}$에 소거법을 적용하면 $U\mathbf{x} = \mathbf{c}$가 된다.

**예제:**

$$\begin{pmatrix} 2 & 3 & 4 \\ 0 & 5 & 6 \\ 0 & 0 & 7 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 19 \\ 17 \\ 14 \end{pmatrix}$$

후진 대입으로 $\mathbf{x}$를 구한다:

1. 3번째 행에서: $7x_3 = 14 \implies x_3 = 2$
2. 2번째 행에서: $5x_2 + 6x_3 = 17 \implies 5x_2 + 12 = 17 \implies x_2 = 1$
3. 1번째 행에서: $2x_1 + 3x_2 + 4x_3 = 19 \implies 2x_1 + 3 + 8 = 19 \implies x_1 = 4$

$$\therefore \mathbf{x} = \begin{pmatrix} 4 \\ 1 \\ 2 \end{pmatrix}$$

**참고 사항:**
- **피벗** (pivots) $2, 5, 7$로 나누어야 했다.
- 피벗은 소거 후에 발견된다.
- 0이 피벗이 되는 것은 허용하지 **않는다**. 필요하면 **행 교환**을 수행한다.
- 독립인 열을 가진 모든 정사각 행렬 $A$는 **영이 아닌 피벗**을 가진 삼각 행렬로 축소될 수 있다.

---

### 2.6 각 열에 대한 소거

**예제:**

$$A = \begin{pmatrix} 2 & 3 & 4 \\ 4 & 11 & 14 \\ 2 & 8 & 17 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 19 \\ 55 \\ 50 \end{pmatrix}$$

**1단계:** 확대 행렬 $(A|\mathbf{b})$를 만들고 $R_2 - 2R_1$을 적용:

$$\begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 4 & 11 & 14 & | & 55 \\ 2 & 8 & 17 & | & 50 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 2 & 8 & 17 & | & 50 \end{pmatrix}$$

이 연산은 다음에 대응된다:

$$-2R_1 + R_2 + 0 \cdot R_3 = (-2 \quad 1 \quad 0)\begin{pmatrix} R_1 \\ R_2 \\ R_3 \end{pmatrix}$$

소거 행렬(elimination matrix)을 도입한다:

$$E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

결과는 $(E_{21}A \mid E_{21}\mathbf{b})$이다.

**2단계:** $R_3 - R_1$ 적용:

$$\begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 2 & 8 & 17 & | & 50 \end{pmatrix} \xrightarrow{R_3 - R_1} \begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 0 & 5 & 13 & | & 31 \end{pmatrix}$$

$$E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -1 & 0 & 1 \end{pmatrix}$$

결과: $(E_{31}E_{21}A \mid E_{31}E_{21}\mathbf{b})$

**3단계:** $R_3' - R_2'$ 적용:

$$\begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 0 & 5 & 13 & | & 31 \end{pmatrix} \xrightarrow{R_3' - R_2'} \begin{pmatrix} 2 & 3 & 4 & | & 19 \\ 0 & 5 & 6 & | & 17 \\ 0 & 0 & 7 & | & 14 \end{pmatrix}$$

$$E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}$$

최종 결과: $(E_{32}E_{31}E_{21}A \mid E_{32}E_{31}E_{21}\mathbf{b}) = (U \mid \mathbf{c})$

---

### 2.7 소거법의 실패 가능성

**피벗 위치에 0이 나타날 때** 발생한다.

**예제 (행 교환으로 해결 가능):**

$$\begin{pmatrix} 2 & 3 & 4 \\ 4 & 6 & 14 \\ 2 & 8 & 17 \end{pmatrix} \rightarrow \begin{pmatrix} 2 & 3 & 4 \\ 0 & 0 & 6 \\ 0 & 5 & 13 \end{pmatrix}$$

치환 행렬 $P$를 사용하여 **2행과 3행을 교환**:

$$\begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 2 & 3 & 4 \\ 0 & 0 & 6 \\ 0 & 5 & 13 \end{pmatrix} = \begin{pmatrix} 2 & 3 & 4 \\ 0 & 5 & 13 \\ 0 & 0 & 6 \end{pmatrix}$$

**예제 (해결 불가 — 특이):**

$$A^* = \begin{pmatrix} 2 & 3 & 4 \\ 4 & 6 & 14 \\ 2 & 3 & 17 \end{pmatrix} \rightarrow \begin{pmatrix} 2 & 3 & 4 \\ 0 & 0 & 6 \\ 0 & 0 & 13 \end{pmatrix} = U^*$$

두 번째 열에서 피벗을 찾을 수 없다.

- $A^*$는 완전 계수(full rank)를 가지지 **않는다**.
- $A^*$와 $U^*$는 가역이 **아니다**.
- 1열과 2열이 같은 방향이다.
- $A^* X = \mathbf{0}$은 영이 아닌 해 $X$를 가진다.

---

### 2.8 종속 열과 독립 열

삼각 행렬 $U$가 **완전 계수**를 가지려면 주 대각선에 **0이 없어야** 한다.

완전 계수인 $A$에 대해:
- $U$의 열들은 독립이다.
- $U$의 행들은 독립이다.

$U$의 **대각선에 0**이 있을 때:
- $U$는 **특이** (singular) 행렬이다.
- $U^{-1}$이 존재하지 않는다.
- $A^{-1}$이 존재하지 않는다.
- $A$는 **특이** 행렬이다.

---

### 2.9 행 관점과 열 관점

**행 관점 (The Row Picture):**

각 방정식은 직선(2D), 평면(3D), 또는 초평면을 나타낸다. 해는 이들이 교차하는 곳이다.

**예제 1 — 해 없음 (평행선):**

$$x - 2y = -1, \quad x - 2y = 1 \implies \begin{pmatrix} 1 & -2 \\ 1 & -2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} -1 \\ 1 \end{pmatrix}$$

**예제 2 — 무한히 많은 해 (같은 직선):**

$$x - 2y = 1, \quad x - 2y = 1 \implies \begin{pmatrix} 1 & -2 \\ 1 & -2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$$

해의 직선.

**예제 3 — 유일한 해 (교차하는 직선):**

$$x - 2y = 7, \quad x + y = 2 \implies \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 7 \\ 2 \end{pmatrix}$$

교차점에서의 유일한 해.

**열 관점 (The Column Picture):**

$$A = \begin{pmatrix} 1 & -2 \\ 1 & 1 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 7 \\ 2 \end{pmatrix}$$

열벡터: $\mathbf{a}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $\mathbf{a}_2 = \begin{pmatrix} -2 \\ 1 \end{pmatrix}$

$$x_1 \mathbf{a}_1 + x_2 \mathbf{a}_2 = \mathbf{b}$$

$\mathbf{b}$는 $\mathbf{a}_1$과 $\mathbf{a}_2$의 **일차결합** (linear combination)이다.

$\Leftrightarrow$ $\mathbf{b}$는 $A$의 **열공간** (column space)에 속한다.

---

<br>

## 3. 소거 행렬과 역행렬 (2.2)

### 3.1 소거와 치환의 예제

**예제 1 (치환 불필요):**

$$A = \begin{pmatrix} 2 & 4 & -2 \\ 4 & 9 & -3 \\ -2 & -3 & 7 \end{pmatrix}$$

$$E_{21}: R_2 - 2R_1 \implies E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 2 & 4 & -2 \\ 0 & 1 & 1 \\ -2 & -3 & 7 \end{pmatrix}$$

$$E_{31}: R_3 + R_1 \implies E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 1 & 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 2 & 4 & -2 \\ 0 & 1 & 1 \\ 0 & 1 & 5 \end{pmatrix}$$

$$E_{32}: R_3' - R_2' \implies E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 2 & 4 & -2 \\ 0 & 1 & 1 \\ 0 & 0 & 4 \end{pmatrix} = U$$

같은 연산을 $\mathbf{b} = \begin{pmatrix} 2 \\ 8 \\ 10 \end{pmatrix}$에 적용:

$$\mathbf{b} = \begin{pmatrix} 2 \\ 8 \\ 10 \end{pmatrix} \xrightarrow{E_{21}} \begin{pmatrix} 2 \\ 4 \\ 10 \end{pmatrix} \xrightarrow{E_{31}} \begin{pmatrix} 2 \\ 4 \\ 12 \end{pmatrix} \xrightarrow{E_{32}} \begin{pmatrix} 2 \\ 4 \\ 8 \end{pmatrix} = \mathbf{c}$$

$$(A|\mathbf{b}) \xrightarrow{E} (U|\mathbf{c})$$

**예제 2 (치환 필요):**

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 2 & 3 \\ 0 & 4 & 5 \end{pmatrix}$$

$$E_{21}: R_2 - 2R_1 \implies \begin{pmatrix} 1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 4 & 5 \end{pmatrix}$$

0 피벗 발생! $P = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}$를 사용하여 행 교환:

$$\begin{pmatrix} 1 & 1 & 1 \\ 0 & 4 & 5 \\ 0 & 0 & 1 \end{pmatrix} = U$$

**$P$를 $A$에 먼저 적용할 수 있는가?**

$$PA = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 4 & 5 \\ 2 & 2 & 3 \end{pmatrix}$$

$$E_{31}: R_3 - 2R_1 \implies E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 1 & 1 & 1 \\ 0 & 4 & 5 \\ 0 & 0 & 1 \end{pmatrix} = U \quad \text{가능하다!}$$

**전체 방정식:**

$$EPA = U \iff \boxed{PA = E^{-1}U = LU}$$

---

### 3.2 소거 행렬과 A = LU

1. 소거는 $A$에 $E_{21}, E_{31}, \ldots, E_{n1}$을 곱하고, 그 다음 $E_{32}, E_{42}, \ldots, E_{n2}$를 곱하여 $A$가 $EA = U$가 된다.

2. 역순으로, $E$들의 역행렬이 $U$에 곱해져 $A = E^{-1}U$를 복원한다. 이것이 $A = LU$이다.

3. $A^{-1}A = I$이고 $(LU)^{-1} = U^{-1}L^{-1}$이다. 그러면 $A\mathbf{x} = \mathbf{b}$는 다음이 된다:

$$\mathbf{x} = A^{-1}\mathbf{b} = U^{-1}L^{-1}\mathbf{b}$$

소거의 모든 단계는 행렬로 수행될 수 있다:

$$E_{n2} \cdots E_{42} E_{32} A = C$$

이 단계들은 행렬로 되돌릴 수 있다:

$$A = E_{32}^{-1} E_{42}^{-1} \cdots E_{n2}^{-1} C$$

**예제:**

$A = \begin{pmatrix} 3 & 1 & 0 \\ -3 & 1 & 1 \\ 6 & 8 & 4 \end{pmatrix}$라 하자.

$$E_{21}: R_2 + R_1 \implies E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$E_{31}: R_3 - 2R_1 \implies E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 6 & 4 \end{pmatrix}$$

$$E_{32}: R_3' - 3R_2' \implies E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -3 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 1 \end{pmatrix} = U$$

$E = E_{32}E_{31}E_{21}$로 놓으면:

$$EA = U \implies A = E^{-1}U = LU$$

---

### 3.3 역행렬에 관한 사실들

**정의:** 행렬 $A$가 **가역** (invertible)이란 $A$를 역변환하는 행렬 $A^{-1}$이 존재하는 것이다:

$$A^{-1}A = AA^{-1} = I$$

$A \in \mathbb{R}^{n \times n}$라 하자. $A$가 $n$개의 독립인 열을 가지면, $A$는 가역이다. 이는 $\text{rank}(A) = n$을 의미한다.

**참고 1:** 역행렬은 소거가 (행 교환과 함께) $n$개의 피벗을 생성할 때에만 존재한다.

소거는 명시적으로 $A^{-1}$을 구하지 않고 $A\mathbf{x} = \mathbf{b}$를 푼다.

**참고 2:** 행렬 $A$는 **두 개의 서로 다른 역행렬을 가질 수 없다**.

$BA = I$이고 $AC = I$라 가정하자. 그러면 $B = C$이다.

*증명:* 결합법칙에 의해,

$$BAC = B(AC) = (BA)C$$

에서 $B = C$가 된다. $\square$

좌역행렬 (left inverse) $B$와 우역행렬 (right inverse) $C$는 같아야 한다.

**참고 3:** $A$가 가역이면, $A\mathbf{x} = \mathbf{b}$에 대해 $\exists! \; \mathbf{x}$ s.t. $\mathbf{x} = A^{-1}\mathbf{b}$.

*증명:*

$$A\mathbf{x} = \mathbf{b} \implies A^{-1}A\mathbf{x} = A^{-1}\mathbf{b} \implies \mathbf{x} = A^{-1}\mathbf{b}. \quad \square$$

**참고 4:** 영이 아닌 벡터 $\mathbf{x}$가 존재하여 $A\mathbf{x} = \mathbf{0}$이라 하자. 그러면:

- $A$는 **종속인 열**을 가진다
- $A$는 역행렬을 가질 수 없다
- 어떤 행렬도 $\mathbf{0}$을 $\mathbf{x}$로 되돌릴 수 없다

**예제:**

$$\begin{pmatrix} 1 & 2 \\ 1 & 2 \end{pmatrix}\begin{pmatrix} -2 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$A$가 가역이면, $A\mathbf{x} = \mathbf{0}$은 $\mathbf{x} = A^{-1}\mathbf{0} = \mathbf{0}$을 의미한다.

**참고 5:** 정사각 행렬은 열들이 독립일 때에만 가역이다.

**참고 6:** $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$는 $ad - bc \neq 0$ ($A$의 **행렬식** (determinant))일 때에만 가역이다.

$$A^{-1} = \frac{1}{ad - bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$$

**참고 7:** 삼각 행렬은 **영이 아닌 대각 원소**를 가지면 역행렬을 가진다.

$A$가 대각 원소 $d_1, d_2, \ldots, d_n$ (모두 영이 아닌)을 가진 상삼각 행렬이면, $A^{-1}$도 대각 원소 $\frac{1}{d_1}, \frac{1}{d_2}, \ldots, \frac{1}{d_n}$을 가진 상삼각 행렬이다.

**예제 2:**

$$A = \begin{pmatrix} 1 & 2 \\ 1 & 2 \end{pmatrix} \xrightarrow{R_2 - R_1} \begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}$$

$\text{rank}(A) = 1 < 2$. $A$는 1개의 피벗을 가진다. $\det(A) = 1 \cdot 2 - 2 \cdot 1 = 0$. $A$는 종속인 열을 가진다.

**예제 3:**

$$A = \begin{pmatrix} 4 & 3 \\ 8 & 6 \end{pmatrix} \quad \text{rank}(A) = 1 < 2 \quad \text{(가역이 아님)}$$

$$B = \begin{pmatrix} 4 & 3 \\ 8 & 7 \end{pmatrix} \quad \det(B) = 4 \cdot 7 - 3 \cdot 8 = 4 \neq 0 \quad B^{-1} = \frac{1}{4}\begin{pmatrix} 7 & -3 \\ -8 & 4 \end{pmatrix}$$

$$C = \begin{pmatrix} 6 & 6 \\ 6 & 0 \end{pmatrix} \quad \det(C) = 6 \cdot 0 - 6 \cdot 6 = -36 \neq 0 \quad C^{-1} = -\frac{1}{36}\begin{pmatrix} 0 & -6 \\ -6 & 6 \end{pmatrix} = \frac{1}{6}\begin{pmatrix} 0 & 1 \\ 1 & -1 \end{pmatrix}$$

$$D = \begin{pmatrix} 6 & 6 \\ 6 & 6 \end{pmatrix} \quad \text{rank}(D) = 1 < 2 \quad \text{(가역이 아님)}$$

$$S = \begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_2 - R_1, R_3 - R_1} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix} \xrightarrow{R_3' - R_2'} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

3개의 피벗. 가역.

$E = E_{32}E_{31}E_{21} = S^{-1}$을 계산:

$$E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -1 & 0 & 1 \end{pmatrix}, \quad E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}$$

$$E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ -1 & 0 & 1 \end{pmatrix}$$

$$E_{32}(E_{31}E_{21}) = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix} = E = S^{-1}$$

$$T = \begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 0 \\ 1 & 1 & 1 \end{pmatrix} \xrightarrow{R_3 - R_2} \begin{pmatrix} 1 & 1 & 1 \\ 1 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

$\text{Rank}(T) = 2 < 3$. **가역이 아니다.**

---

### 3.4 곱 AB의 역행렬

영이 아닌 두 값 $a$와 $b$에 대해, 합 $(a + b)$는 가역이 아닐 수 있다.

**예제:** $a = 3 \implies a^{-1} = 1/3$, $b = -3 \implies b^{-1} = -1/3$, $a + b = 0 \implies (a+b)^{-1}$은 존재하지 않는다.

그러나: $ab = -9 \implies (ab)^{-1} = -1/9 = a^{-1}b^{-1}$.

**정리:** $A, B \in \mathbb{R}^{n \times n}$이 가역이면, $AB$의 역행렬은:

$$\boxed{(AB)^{-1} = B^{-1}A^{-1}}$$

*증명:*

$$(AB)^{-1}AB = I$$
$$(AB)^{-1}ABB^{-1} = IB^{-1} = B^{-1}$$
$$(AB)^{-1}AA^{-1} = B^{-1}A^{-1}$$

$$\therefore (AB)^{-1} = B^{-1}A^{-1} \quad \square$$

**역행렬은 역순으로 나온다!**

세 행렬에 대해:

$$(ABC)^{-1} = C^{-1}B^{-1}A^{-1}$$

*증명:*

$$(ABC)^{-1}ABC = I \implies (ABC)^{-1}ABCC^{-1} = C^{-1} \implies (ABC)^{-1}ABB^{-1} = C^{-1}B^{-1}$$

$$\implies (ABC)^{-1}AA^{-1} = C^{-1}B^{-1}A^{-1} \implies (ABC)^{-1} = C^{-1}B^{-1}A^{-1}$$

---

### 3.5 소거 행렬의 역행렬

**예제 4:**

$$E = \begin{pmatrix} 1 & 0 & 0 \\ -5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \quad \text{(1행의 5배를 2행에서 뺌)}$$

$(-5 \quad 1 \quad 0)\begin{pmatrix} R_1 \\ R_2 \\ R_3 \end{pmatrix} = R_2 - 5R_1 = R_2'$

되돌리기: $R_2 = R_2' + 5R_1$, 즉 $(5 \quad 1 \quad 0)\begin{pmatrix} R_1 \\ R_2' \\ R_3 \end{pmatrix}$

$$E^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \quad \text{(1행의 5배를 2행에 더함)}$$

$$EE^{-1} = E^{-1}E = I$$

**중요:** 정사각 행렬 $A, C$에 대해 $AC = I$이면, $CA = I$이다.

**예제 5:**

$$F = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -4 & 1 \end{pmatrix} \quad (R_3' = R_3 - 4R_2)$$

$$F^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 4 & 1 \end{pmatrix} \quad (R_3 = R_3' + 4R_2)$$

곱 $FE$ 계산:

$$FE = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -4 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ -5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ -5 & 1 & 0 \\ 20 & -4 & 1 \end{pmatrix}$$

"20"은 1행에서 연산의 연쇄로 나온다:

$$R_3'' = R_3' - 4R_2' = R_3 - 4(R_2 - 5R_1) = R_3 - 4R_2 + 20R_1$$

$$(FE)^{-1} = E^{-1}F^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 5 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 4 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 5 & 1 & 0 \\ 0 & 4 & 1 \end{pmatrix}$$

승수(multiplier) 5와 4가 $L = (FE)^{-1}$의 대각선 아래에 제자리에 놓인다.

---

### 3.6 L은 E의 역행렬

**L은 E의 역행렬이다.**

예제 1을 상기하자: $A = \begin{pmatrix} 3 & 1 & 0 \\ -3 & 1 & 1 \\ 6 & 8 & 4 \end{pmatrix}$

$$E = E_{32}E_{31}E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -3 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$EA = \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 1 \end{pmatrix} = U$$

**E와 L에 대한 일반 공식 (3x3 경우):**

$l_{32} = 3$, $l_{31} = 2$, $l_{21} = -1$로 놓자. 그러면:

$$E = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -l_{32} & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -l_{31} & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ -l_{21} & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$= \begin{pmatrix} 1 & 0 & 0 \\ -l_{21} & 1 & 0 \\ l_{32}l_{21} - l_{31} & -l_{32} & 1 \end{pmatrix}$$

(3,1) 위치에 교차곱 항 $l_{32}l_{21} - l_{31}$이 있음에 주목하라.

$$E^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ l_{31} & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & l_{32} & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{pmatrix} = L$$

**역행렬 $E^{-1}$이 아름답게 된다!** 인수 $l_{21}, l_{31}, l_{32}$를 볼 수 있다.

**모든 승수 $l_{ij}$가 $L$에서 올바른 위치에 나타난다.**

---

<br>

## 4. 행렬 연산과 A = LU (2.3)

### 4.1 핵심 사실

1. $A$에서 $U$로의 소거 단계는 $\frac{1}{3}n^3$번의 곱셈과 뺄셈이 소요된다.
2. 각 우변 $\mathbf{b}$는 $n^2$만 소요된다:
   - $U\mathbf{x} = \mathbf{c}$로 전진
   - 그 다음 $\mathbf{x}$에 대한 후진 대입
3. 행 교환 없는 소거는 $A$를 $LU$로 분해한다.

$A\mathbf{x} = \mathbf{b}$의 해는 $\mathbf{x} = A^{-1}\mathbf{b}$로 주어진다.

**Q:** $A^{-1}$을 명시적으로 알아야 하는가? $A^{-1}$을 계산하고 $A^{-1}\mathbf{b}$를 곱하는 것은 $\mathbf{x}$를 구하는 매우 **느린** 방법이다.

---

### 4.2 역행렬의 명시적 계산

**Q:** $A^{-1}$을 명시적으로 어떻게 구하는가?

$AA^{-1} = I \in \mathbb{R}^{n \times n}$에서 시작한다.

$$I = \begin{pmatrix} | & | & & | \\ \hat{e}_1 & \hat{e}_2 & \cdots & \hat{e}_n \\ | & | & & | \end{pmatrix}$$

여기서 $\hat{e}_1, \hat{e}_2, \ldots, \hat{e}_n$은 **표준 기저 벡터** (단위 벡터)이다.

$AA^{-1} = I$를 다음과 같이 본다:

$$A\begin{pmatrix} | & | & & | \\ \mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_n \\ | & | & & | \end{pmatrix} = \begin{pmatrix} | & | & & | \\ \hat{e}_1 & \hat{e}_2 & \cdots & \hat{e}_n \\ | & | & & | \end{pmatrix}$$

즉 $n$개의 방정식이다:

$$A\mathbf{x}_1 = \hat{e}_1, \quad A\mathbf{x}_2 = \hat{e}_2, \quad \ldots, \quad A\mathbf{x}_n = \hat{e}_n$$

같은 계수 행렬 $A$, 서로 다른 우변 벡터.

---

### 4.3 가우스-조르단 소거법

$n$개의 방정식을 **가우스-조르단 소거법** (Gauss-Jordan elimination)을 사용하여 함께 푼다:

$$(A \mid I) \implies (I \mid A^{-1})$$

**예제:**

$$\begin{pmatrix} 1 & 0 & 0 & | & 1 & 0 & 0 \\ -1 & 1 & 0 & | & 0 & 1 & 0 \\ 0 & -1 & 1 & | & 0 & 0 & 1 \end{pmatrix} = (A|I)$$

$$\xrightarrow{R_2 + R_1} \begin{pmatrix} 1 & 0 & 0 & | & 1 & 0 & 0 \\ 0 & 1 & 0 & | & 1 & 1 & 0 \\ 0 & -1 & 1 & | & 0 & 0 & 1 \end{pmatrix}$$

$$\xrightarrow{R_3 + R_2} \begin{pmatrix} 1 & 0 & 0 & | & 1 & 0 & 0 \\ 0 & 1 & 0 & | & 1 & 1 & 0 \\ 0 & 0 & 1 & | & 1 & 1 & 1 \end{pmatrix} = (I|A^{-1})$$

$A$에 대한 소거 단계는 **한 번만** 수행하면 된다!

---

### 4.4 소거의 비용

**$A$를 $U$로 축소:**

1단계 (1열 소거): $(n-1)$개 행, $n$개 열 $\implies (n-1)n$번의 곱셈과 $(n-1)n$번의 뺄셈.

2단계 (2열 소거): $(n-2)$개 행, $(n-1)$개 열 $\implies (n-2)(n-1)$번의 곱셈과 뺄셈.

$\vdots$

$(n-1)$단계 ($n-1$열 소거): 1개 행, 2개 열 $\implies 1 \cdot 2$번의 곱셈과 뺄셈.

**총 곱셈 횟수:**

$$(n-1)n + (n-2)(n-1) + \cdots + 1 \cdot 2 = \sum_{i=1}^{n-1} i(i+1) = \sum_{i=1}^{n-1} i^2 + \sum_{i=1}^{n-1} i$$

$$= \frac{n(n+1)(2n+1)}{6} - n^2 + \frac{n^2}{2} - \frac{n}{2} + \frac{(n-1)n}{2}$$

$$= \frac{1}{3}n^3 - \frac{n}{3}$$

$n \to \infty$일 때: $\approx \dfrac{1}{3}n^3$.

**$A$를 $U$로 축소하는 데 약 $\frac{1}{3}n^3$번의 곱셈과 $\frac{1}{3}n^3$번의 뺄셈이 필요하다.**

**$\mathbf{b}$를 $\mathbf{c}$로 축소:**

$A$와 유사하지만 1개 열만:

- 1단계: $(n-1)$번의 곱셈과 뺄셈
- 2단계: $(n-2)$번의 곱셈과 뺄셈
- ...
- $(n-1)$단계: 1번의 곱셈과 뺄셈

$$\sum_{i=1}^{n-1} i = \frac{(n-1)n}{2} = \frac{n^2}{2} - \frac{n}{2}$$

$n \to \infty$일 때: $\approx \dfrac{n^2}{2}$.

**$\mathbf{b}$를 $\mathbf{c}$로 축소하는 데 $\frac{n^2}{2}$번의 곱셈과 $\frac{n^2}{2}$번의 뺄셈이 필요하다.**

**후진 대입** ($U\mathbf{x} = \mathbf{c}$):

$$x_n = c_n / u_{nn}$$
$$x_{n-1} = (c_{n-1} - u_{(n-1)n}x_n) / u_{(n-1)(n-1)}$$
$$\vdots$$
$$x_1 = (c_1 - u_{12}x_2 - u_{13}x_3 - \cdots - u_{1n}x_n) / u_{11}$$

$n$번의 나눗셈, $1 + 2 + \cdots + (n-1)$번의 곱셈, $1 + 2 + \cdots + (n-1)$번의 뺄셈.

$$n + \frac{(n+1)n}{2} + \frac{(n-1)n}{2} = n^2$$

**$\mathbf{b}$에서 $\mathbf{c}$를 거쳐 $\mathbf{x}$까지 우변의 총 연산 횟수는 $n^2$이다:**
- $n^2$번의 곱셈
- $n^2$번의 뺄셈

---

### 4.5 위대한 분해 A = LU

하나의 소거 단계 $E_{ij}$를 역변환하려면 (행 $i$에서 행 $j$의 $l_{ij}$배를 뺀 것을), 빼는 대신 **더한다**.

$$E_{31} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -l_{31} & 0 & 1 \end{pmatrix} \implies L_{31} = E_{31}^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ l_{31} & 0 & 1 \end{pmatrix}$$

**예제 (재방문):**

$A = \begin{pmatrix} 3 & 1 & 0 \\ -3 & 1 & 1 \\ 6 & 8 & 4 \end{pmatrix}$로 놓자.

$$E_{21}: R_2 + R_1 \implies \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 6 & 8 & 4 \end{pmatrix}$$

$$E_{31}: R_3 - 2R_1 \implies \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 6 & 4 \end{pmatrix}$$

$$E_{32}: R_3' - 3R_2' \implies \begin{pmatrix} 3 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 1 \end{pmatrix} = U$$

$U$의 1행 = $A$의 1행임에 주목하라. $U$의 2행은 $E_{31}$과 $E_{32}$ 적용 후에도 변하지 않는다.

**$A$와 $U$의 행 사이의 관계:**

$U$의 3행 = $A$의 3행 $-$ 2 ($U$의 1행) $-$ 3 ($U$의 2행)

동치로:

$A$의 3행 = $U$의 3행 $+$ $l_{31}$ ($U$의 1행) $+$ $l_{32}$ ($U$의 2행)

여기서 $l_{31} = 2$이고 $l_{32} = 3$이다.

$$\text{$A$의 3행} = (l_{31} \quad l_{32} \quad 1) \begin{pmatrix} U_1 \\ U_2 \\ U_3 \end{pmatrix}$$

$$\implies A = LU$$

---

### 4.6 A = LU의 두 번째 증명

**열 곱하기 행 (columns times rows).**

$A = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{pmatrix}$을 고려하자.

**1단계:** 1행을 피벗 행으로 취한다. 1행에 $l_{21}, l_{31}$을 곱하여 2행, 3행에서 뺀다.

$l_{21} = a_{21}/a_{11}$, $l_{31} = a_{31}/a_{11}$으로 선택한다.

$$R_2 - l_{21}R_1, \quad R_3 - l_{31}R_1$$

$$A' = \begin{pmatrix} a_{11} & a_{12} & a_{13} \\ 0 & a'_{22} & a'_{23} \\ 0 & a'_{32} & a'_{33} \end{pmatrix}$$

이 단계는 다음과 같이 볼 수 있다:

$$A = A' + \begin{pmatrix} - & 0 & - \\ - & l_{21}R_1 & - \\ - & l_{31}R_1 & - \end{pmatrix}$$

뺀 부분은 **랭크 1 행렬** (rank 1 matrix)이다:

$$\begin{pmatrix} 1 \\ l_{21} \\ l_{31} \end{pmatrix}\begin{pmatrix} a_{11} & a_{12} & a_{13} \end{pmatrix} = \mathbf{l}_1 \mathbf{u}_1 = \mathbf{l}_1 \otimes \mathbf{u}_1$$

$$= \begin{pmatrix} 1 \cdot a_{11} & 1 \cdot a_{12} & 1 \cdot a_{13} \\ l_{21}a_{11} & l_{21}a_{12} & l_{21}a_{13} \\ l_{31}a_{11} & l_{31}a_{12} & l_{31}a_{13} \end{pmatrix}$$

**2단계:** $A_2$ (나머지 부분)의 2행을 피벗 행으로 취한다.

$$A_2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & a'_{22} & a'_{23} \\ 0 & a'_{32} & a'_{33} \end{pmatrix}$$

$R_3' - l_{32}R_2'$ 적용:

$$A'' = \begin{pmatrix} 0 & 0 & 0 \\ 0 & a'_{22} & a'_{23} \\ 0 & 0 & a''_{33} \end{pmatrix} = A_2 - \begin{pmatrix} 0 \\ 1 \\ l_{32} \end{pmatrix}(R_2')$$

이것은 또 다른 랭크 1 행렬을 생성한다:

$$\begin{pmatrix} 0 \\ 1 \\ l_{32} \end{pmatrix}\begin{pmatrix} 0 & a'_{22} & a'_{23} \end{pmatrix} = \mathbf{l}_2 \mathbf{u}_2$$

$\mathbf{u}_2$는 $U$의 2행임에 주목하라.

나머지 부분:

$$A_3 = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}\begin{pmatrix} 0 & 0 & a''_{33} \end{pmatrix} = \mathbf{l}_3 \mathbf{u}_3$$

이제 $A$를 다음과 같이 표현한다:

$$A = \mathbf{l}_1 \mathbf{u}_1 + \mathbf{l}_2 \mathbf{u}_2 + \mathbf{l}_3 \mathbf{u}_3$$

$$= \begin{pmatrix} | & | & | \\ \mathbf{l}_1 & \mathbf{l}_2 & \mathbf{l}_3 \\ | & | & | \end{pmatrix}\begin{pmatrix} - & \mathbf{u}_1 & - \\ - & \mathbf{u}_2 & - \\ - & \mathbf{u}_3 & - \end{pmatrix}$$

$$= \begin{pmatrix} 1 & 0 & 0 \\ l_{21} & 1 & 0 \\ l_{31} & l_{32} & 1 \end{pmatrix}\begin{pmatrix} a_{11} & a_{12} & a_{13} \\ 0 & a'_{22} & a'_{23} \\ 0 & 0 & a''_{33} \end{pmatrix} = LU$$

**일반 확장** — $A \in \mathbb{R}^{n \times n}$에 대해:

$$A = \mathbf{l}_1\mathbf{u}_1 + \mathbf{l}_2\mathbf{u}_2 + \cdots + \mathbf{l}_n\mathbf{u}_n$$

$$= \begin{pmatrix} 1 & 0 & 0 & \cdots & 0 \\ l_{21} & 1 & 0 & \cdots & 0 \\ l_{31} & l_{32} & 1 & \cdots & 0 \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ l_{n1} & l_{n2} & l_{n3} & \cdots & 1 \end{pmatrix}\begin{pmatrix} a_{11} & a_{12} & a_{13} & \cdots & a_{1n} \\ 0 & a'_{22} & a'_{23} & \cdots & a'_{2n} \\ 0 & 0 & a''_{33} & \cdots & a''_{3n} \\ \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & a_{nn}^{(n-1)} \end{pmatrix} = LU$$

$\mathbf{l}_k$는 $(k-1)$개의 0으로 시작한다. $\mathbf{u}_k$는 $(k-1)$개의 0으로 시작한다.

---

### 4.7 행 교환 없는 소거

**Q:** $A = LU$가 **행 교환 없이** **피벗에 0 없이** 가능한 것은 언제인가?

**A:** $A$의 모든 좌상단 $k \times k$ 부분행렬이 가역이어야 한다.

3x3 행렬 $A$에 대해:
- $A_1 = (a_{11})$: $A_1 = L_1 U_1$ (1x1 부분행렬이 가역이어야 함)
- $A_2 = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}$: $A_2 = L_2 U_2$ (2x2 부분행렬이 가역이어야 함)
- $A_3 = A$ (전체 행렬): $A_3 = L_3 U_3$

---

<br>

## 5. 치환과 전치 (2.4)

### 5.1 치환 행렬

**치환 행렬** (permutation matrix) $P$는 $I \in \mathbb{R}^{n \times n}$과 같은 행을 가진다.

$n!$가지 서로 다른 순서가 있다.

**예제:** $P \in \mathbb{R}^{3 \times 3}$: 3개의 행. 1행에: 3가지; 2행에: 2가지; 3행에: 1가지 $\implies 3! = 6$가지 순서.

$P$ 곱하기 $\mathbf{x}$는 성분 $x_1$부터 $x_n$까지를 새로운 순서로 배치한다.

그리고 $P^T$는 $P^{-1}$과 같다.

**예제:**

$$P\mathbf{x} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \end{pmatrix}$$

$$P^T = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$$

$P^{-1} = ?$ 가우스-조르단 방법 $(A|I) \Rightarrow (I|A^{-1})$ 사용:

$$\begin{pmatrix} 0 & 0 & 1 & | & 1 & 0 & 0 \\ 1 & 0 & 0 & | & 0 & 1 & 0 \\ 0 & 1 & 0 & | & 0 & 0 & 1 \end{pmatrix} \xrightarrow{\text{행 교환}} \begin{pmatrix} 1 & 0 & 0 & | & 0 & 1 & 0 \\ 0 & 1 & 0 & | & 0 & 0 & 1 \\ 0 & 0 & 1 & | & 1 & 0 & 0 \end{pmatrix}$$

$$P^{-1} = P^T$$

**전치의 성질:**

- $A$의 열은 $A^T$의 행이다.
- $A\mathbf{x}$와 $AB$의 전치는 $\mathbf{x}^T A^T$와 $B^T A^T$이다.

**내적의 성질:**

$$A\mathbf{x} \cdot \mathbf{y} = \mathbf{x} \cdot A^T\mathbf{y}$$

왜냐하면 $(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T A^T \mathbf{y} = \mathbf{x}^T(A^T\mathbf{y})$

**대칭 행렬:** $S^T = S$. 곱 $A^T A$와 $AA^T$는 항상 대칭이다.

---

### 5.2 치환 행렬의 성질

치환 행렬은 모든 행에 정확히 하나의 1을, 모든 열에 정확히 하나의 1을 가진다.

$P$를 벡터 $\mathbf{x}$에 곱하면, 성분의 순서가 바뀐다:

$$P\mathbf{x} = \begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \end{pmatrix}$$

$P$는 $x_1$을 두 번째 위치로 이동시킨다.

$$PP\mathbf{x} = P^2\mathbf{x} = \begin{pmatrix} x_2 \\ x_3 \\ x_1 \end{pmatrix}$$

$P^2$는 $x_1$을 세 번째 위치로 이동시킨다.

$$PPP\mathbf{x} = P^3\mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = I\mathbf{x}$$

$P^3$은 $x_1$을 원래 위치로 되돌린다. $P^3 = I$.

**4x4 치환 행렬을 고려하자:**

**(a)** $P$는 $\mathbf{x}$의 순서를 뒤집는다:

$$\begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_4 \\ x_3 \\ x_2 \\ x_1 \end{pmatrix}$$

$P(P\mathbf{x}) = \mathbf{x} \implies P^2 = I$

**(b)** $P$는 $x_4$ 위치를 바꾸지 않는다:

$$\begin{pmatrix} 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_3 \\ x_1 \\ x_2 \\ x_4 \end{pmatrix}$$

$P(P(P\mathbf{x})) = \mathbf{x} \implies P^3 = I$

**(c)** $P$는 원소를 순환적으로 이동시킨다:

$$\begin{pmatrix} 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix} = \begin{pmatrix} x_4 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix}$$

$PPPP\mathbf{x} = P^4\mathbf{x} = I\mathbf{x} \implies P^4 = I$

**(d)** 짝수-홀수 분리:

$$\begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_0 \\ x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} x_0 \\ x_2 \\ x_1 \\ x_3 \end{pmatrix} \quad \text{(짝수, 홀수 분리)}$$

$P^2 = I$

이것은 짝수 인덱스와 홀수 인덱스 항목을 분리하는 8x8 치환 행렬을 사용하여 $\mathbf{x} \in \mathbb{R}^8$ 벡터로 확장된다.

**$P^T P = I$의 증명:**

임의의 $P$의 행은 $P^{-1} = P^T$의 열이다.

$$P^T P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}\begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$$

$$= \mathbf{h}_1\mathbf{h}_1^T + \mathbf{h}_2\mathbf{h}_2^T + \mathbf{h}_3\mathbf{h}_3^T$$

$\mathbf{h}_i$가 표준 기저 벡터 (canonical unit vector)이므로:

$$= \begin{pmatrix} 1 & & \\ & 1 & \\ & & 1 \end{pmatrix} = I$$

**치환 행렬의 성질:**

1. 치환 행렬 $P$는 각 행에 정확히 하나의 1을, 각 열에 정확히 하나의 1을 가진다.

2. $P$의 열들은 **직교** (orthogonal)한다. (열 사이의 내적이 모두 0이다.)

3. 치환의 곱 $P_1 P_2$는 치환이다. $P$의 역행렬도 마찬가지이다.

4. $A$가 가역이면, 행을 미리 정렬하는 치환 $P$가 존재하여, $PA$에 대한 소거가 피벗 위치에서 0을 만나지 않는다:

$$PA = LU$$

---

### 5.3 PA = LU 분해

**$P$에 의한 행 교환:**

행렬 $A$를 고려하자:

$$A = \begin{pmatrix} 1 & 2 & a \\ 2 & 4 & b \\ 3 & 7 & c \end{pmatrix}$$

1행에서 1을 피벗으로 취한다:

$$R_2 - 2R_1, \quad R_3 - 3R_1$$

$$EA = \begin{pmatrix} 1 & 2 & a \\ 0 & 0 & b - 2a \\ 0 & 1 & c - 3a \end{pmatrix}$$

0 피벗으로 인해, 2행과 3행을 교환한다:

$$PEA = \begin{pmatrix} 1 & 2 & a \\ 0 & 1 & c - 3a \\ 0 & 0 & b - 2a \end{pmatrix} = U$$

$A$는 $b - 2a \neq 0$일 때에만 가역이다. $b = 2a$이면, $\text{rank}(A) = 2 < 3$이고 $A$는 가역이 아니다.

**먼저 2행과 3행을 교환할 수 있다:**

$$PA = \begin{pmatrix} 1 & 2 & a \\ 3 & 7 & c \\ 2 & 4 & b \end{pmatrix}$$

$$\xrightarrow{R_2 - 3R_1, R_3 - 2R_1}$$

$$EPA = \begin{pmatrix} 1 & 2 & a \\ 0 & 1 & c - 3a \\ 0 & 0 & b - 2a \end{pmatrix} = U$$

$$\therefore PA = E^{-1}U = LU$$

**Daniel Drucker의 $P$ 추적 방법:**

행 인덱스를 추적하는 열로 행렬을 확대한다:

$$\begin{pmatrix} 1 & 2 & a & | & 1 \\ 2 & 4 & b & | & 2 \\ 3 & 7 & c & | & 3 \end{pmatrix} \xrightarrow{\text{소거}} \begin{pmatrix} 1 & 2 & a & | & 1 \\ 0 & 0 & b-2a & | & 2 \\ 0 & 1 & c-3a & | & 3 \end{pmatrix} \xrightarrow{\text{교환}} \begin{pmatrix} 1 & 2 & a & | & 1 \\ 0 & 1 & c-3a & | & 3 \\ 0 & 0 & b-2a & | & 2 \end{pmatrix}$$

최종 열은 $P_{132}$를 준다:

$$P_{132} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}$$

---

### 5.4 부분 피벗팅

반올림 오차를 줄이기 위한 **"부분 피벗팅"** (Partial Pivoting).

피벗에 **가능한 가장 큰 수**를 생성하도록 행을 교환하면 계산이 더 안정적이다.

**예제:**

$$\begin{pmatrix} 1 & 2 & a & | & 1 \\ 2 & 4 & b & | & 2 \\ 3 & 7 & c & | & 3 \end{pmatrix}$$

$$\xrightarrow{R_3 \leftrightarrow R_1} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 2 & 4 & b & | & 2 \\ 1 & 2 & a & | & 1 \end{pmatrix}$$

$$\xrightarrow{R_2 - \frac{2}{3}R_1, R_3 - \frac{1}{3}R_1} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 0 & -\frac{2}{3} & b - \frac{2}{3}c & | & 2 \\ 0 & -\frac{1}{3} & a - \frac{1}{3}c & | & 1 \end{pmatrix}$$

$$\xrightarrow{R_3' - R_2'(\frac{1}{2})} \begin{pmatrix} 3 & 7 & c & | & 3 \\ 0 & -\frac{2}{3} & b - \frac{2}{3}c & | & 2 \\ 0 & 0 & a - \frac{1}{2}b & | & 1 \end{pmatrix} = U$$

각 피벗을 그 아래의 모든 수보다 크게 만들면 $L$의 모든 원소가 $\leq 1$이 된다:

$$L = \begin{pmatrix} 1 & 0 & 0 \\ 2/3 & 1 & 0 \\ 1/3 & 1/2 & 1 \end{pmatrix}$$

$$P_{321} = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}$$

---

### 5.5 PAQ: 행과 열 치환

$PAQ$는 행 치환 $P$와 열 치환 $Q$를 가진다.

$A \in \mathbb{R}^{3 \times 3}$에서 시작한다. $P \in \mathbb{R}^{3 \times 3}$로 행을 재정렬한다.

$$P = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}, \quad PA = \begin{pmatrix} a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \\ a_{11} & a_{12} & a_{13} \end{pmatrix}$$

오른쪽에서 $Q = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix}$을 곱한다:

$$(PA)Q = \begin{pmatrix} a_{23} & a_{22} & a_{21} \\ a_{33} & a_{32} & a_{31} \\ a_{13} & a_{12} & a_{11} \end{pmatrix}$$

열 치환 $Q$는 열을 재정렬한다.

$A$의 열공간은 $PA$의 열공간과 같은가?

$$C(A) \stackrel{?}{=} C(PA)$$

**그렇다**, $P$는 일차 관계를 바꾸지 않는다. 따라서 $C(A) = C(PA)$.

**Q:** 행렬 $A$는 9개의 수를 가진다. $A$의 9개의 수를 몇 가지 다른 방법으로 배열할 수 있는가? **A:** $9!$

$P, Q$는 $PAQ$에 대해 $6 \times 6 = 36$가지로 수를 배열할 수 있다. $PAQ$는 $C(A) = C(PAQ)$를 만족하는 매우 특별한 것이다.

---

### 5.6 A의 전치

$A$의 **전치** (transpose)는 $A^T$로 표기한다. $A^T$의 열은 $A$의 행이다.

$$A^T = \begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & \vdots & \ddots & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}_{n \times m}$$

$$A = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}_{m \times n}$$

**예제:**

$$A = \begin{pmatrix} 1 & 2 & 3 \\ 0 & 0 & 4 \end{pmatrix}, \quad A^T = \begin{pmatrix} 1 & 0 \\ 2 & 0 \\ 3 & 4 \end{pmatrix}$$

행렬은 주 대각선을 기준으로 "뒤집힌다".

$$(A^T)_{ij} = A_{ji}$$

**전치의 규칙:**

- **합:** $(A + B)^T = A^T + B^T$
- **곱:** $(AB)^T = B^T A^T$ (역순)
- **역:** $(A^{-1})^T = (A^T)^{-1}$

**$(A\mathbf{x})^T = \mathbf{x}^T A^T$의 증명:**

$$A\mathbf{x} = \begin{pmatrix} \sum_{j=1}^n a_{1j}x_j \\ \sum_{j=1}^n a_{2j}x_j \\ \vdots \\ \sum_{j=1}^n a_{mj}x_j \end{pmatrix}$$

$$(A\mathbf{x})^T = \left(\sum_{j=1}^n a_{1j}x_j \quad \sum_{j=1}^n a_{2j}x_j \quad \cdots \quad \sum_{j=1}^n a_{mj}x_j\right)$$

$$= (x_1, x_2, \ldots, x_n)\begin{pmatrix} a_{11} & a_{21} & \cdots & a_{m1} \\ a_{12} & a_{22} & \cdots & a_{m2} \\ \vdots & & & \vdots \\ a_{1n} & a_{2n} & \cdots & a_{mn} \end{pmatrix}_{n \times m}$$

$$= \mathbf{x}^T A^T$$

**$(AB)^T = B^T A^T$의 증명:**

$B = \begin{pmatrix} | & | & & | \\ \mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_p \\ | & | & & | \end{pmatrix}_{n \times p}$으로 해석한다.

$$(AB)^T = (A\mathbf{x}_1 \quad A\mathbf{x}_2 \quad \cdots \quad A\mathbf{x}_p)^T = \begin{pmatrix} \mathbf{x}_1^T A^T \\ \mathbf{x}_2^T A^T \\ \vdots \\ \mathbf{x}_p^T A^T \end{pmatrix} = \begin{pmatrix} \mathbf{x}_1^T \\ \mathbf{x}_2^T \\ \vdots \\ \mathbf{x}_p^T \end{pmatrix} A^T = B^T A^T$$

**예제:**

$$AB = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} 5 & 0 \\ 4 & 1 \end{pmatrix} = \begin{pmatrix} 5 & 0 \\ 9 & 1 \end{pmatrix}$$

$$B^T A^T = \begin{pmatrix} 5 & 4 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 5 & 9 \\ 0 & 1 \end{pmatrix} = (AB)^T$$

또한: $(ABC)^T = (AB \cdot C)^T = C^T(AB)^T = C^T B^T A^T$. 역순 규칙이 여전히 성립한다.

**$(A^{-1})^T = (A^T)^{-1}$의 증명:**

$A^{-1}A = I$를 고려하면:

$$\Rightarrow (A^{-1}A)^T = I^T = I$$
$$\Leftrightarrow A^T(A^{-1})^T = I$$

$\therefore (A^{-1})^T$는 $A^T$의 역행렬이다:

$$(A^{-1})^T = (A^T)^{-1}$$

이는 $A$가 가역일 때 정확히 $A^T$도 가역임을 의미한다.

**예제:**

$$A = \begin{pmatrix} 1 & 0 \\ 6 & 1 \end{pmatrix}, \quad A^{-1} = \begin{pmatrix} 1 & 0 \\ -6 & 1 \end{pmatrix} \implies (A^{-1})^T = \begin{pmatrix} 1 & -6 \\ 0 & 1 \end{pmatrix}$$

$$A^T = \begin{pmatrix} 1 & 6 \\ 0 & 1 \end{pmatrix} \implies (A^T)^{-1} = \begin{pmatrix} 1 & -6 \\ 0 & 1 \end{pmatrix}$$

일치한다: $(A^{-1})^T = (A^T)^{-1}$.

---

### 5.7 내적과 전치

**내적의 의미.**

$\mathbf{x}$와 $\mathbf{y}$의 점곱(dot product) (내적, inner product)은 수 $x_i y_i$의 합이다:

$$\mathbf{x} \cdot \mathbf{y} = \sum x_i y_i$$

$\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$라 하자:

$$\mathbf{x} \cdot \mathbf{y} = \sum_{i=1}^n x_i y_i = \mathbf{x}^T \mathbf{y} \quad (1 \times n)(n \times 1) = (1 \times 1)$$

$$\mathbf{x}\mathbf{y}^T = \begin{pmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{pmatrix}(y_1, y_2, \ldots, y_n) = \mathbf{x} \otimes \mathbf{y} \quad (n \times 1)(1 \times n) = (n \times n)$$

이것은 **외적** (outer product) (랭크 1 곱, 랭크 1 행렬)이다.

**내적의 예:**

$$\text{일} = \text{힘} \cdot \text{거리} = \mathbf{f}^T \mathbf{d} \quad [J] = [N] \cdot [m]$$

$$\text{수입} = \text{수량} \cdot \text{가격} = \mathbf{q}^T \mathbf{p}$$

**$A^T$의 더 나은 정의:**

$A^T$를 행렬 $A$를 주 대각선을 기준으로 뒤집어 정의하지만, $A^T$의 더 나은 정의는 다음 내적을 같게 만드는 행렬이라는 것이다:

$$(A\mathbf{x}) \cdot \mathbf{y} = \mathbf{x} \cdot (A^T\mathbf{y})$$

즉,

$$(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T A^T \mathbf{y}$$

**예제 1:**

$$A = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}, \quad \mathbf{x} = \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix} y_1 \\ y_2 \end{pmatrix}$$

$$A\mathbf{x} = \begin{pmatrix} x_2 - x_1 \\ x_3 - x_2 \end{pmatrix}$$

$$(A\mathbf{x})^T\mathbf{y} = (x_2 - x_1)y_1 + (x_3 - x_2)y_2$$

$$A^T\mathbf{y} = \begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} y_1 \\ y_2 \end{pmatrix} = \begin{pmatrix} -y_1 \\ y_1 - y_2 \\ y_2 \end{pmatrix}$$

$$\mathbf{x}^T(A^T\mathbf{y}) = (x_1, x_2, x_3)\begin{pmatrix} -y_1 \\ y_1 - y_2 \\ y_2 \end{pmatrix} = -x_1y_1 + x_2(y_1 - y_2) + x_3y_2$$

양쪽이 같다. $\checkmark$

**예제 2: 함수의 내적.**

$$\mathbf{x}^T\mathbf{y} = x_1 y_1 + x_2 y_2 + \cdots + x_n y_n$$

연속 세계에서:

$$(x, y) := \int_{-\infty}^{\infty} x(t) \, y(t) \, dt$$

마찬가지로, $(A\mathbf{x})^T\mathbf{y} = \mathbf{x}^T(A^T\mathbf{y})$:

$$(Ax, y) = (x, A^T y)$$

$A = \frac{d}{dt}$, $A^T = -\frac{d}{dt} = -A$로 놓자.

$$\left(\frac{dx}{dt}, y\right) = \left(x, -\frac{dy}{dt}\right)$$

$$\int_{-\infty}^{\infty} \frac{dx}{dt} y \, dt = \int_{-\infty}^{\infty} x \left(-\frac{dy}{dt}\right) dt$$

**부분적분 (IBP)을 통한 증명:**

$$(f(t)g(t))' = f'(t)g(t) + f(t)g'(t)$$

$$\int f'g \, dt = \int (fg)' \, dt - \int fg' \, dt = (fg)\Big|_{-\infty}^{\infty} - \int fg' \, dt$$

$f(\infty) = f(-\infty) = 0$이라 가정하면:

$$\boxed{\int_{-\infty}^{\infty} f'g \, dt = -\int_{-\infty}^{\infty} fg' \, dt}$$

미분은 **반대칭** (anti-symmetric)이다: $A^T = -A$.

대칭 행렬은 $A^T = A$를 가진다.

---

### 5.8 대칭 행렬

**대칭 행렬** (symmetric matrix)은 $S^T = S$를 만족한다:

$$\implies (S^T)_{ij} = S_{ij} = S_{ji}$$

**예제:**

$$S = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix} = S^T, \quad D = \begin{pmatrix} 1 & 0 \\ 0 & 10 \end{pmatrix} = D^T$$

**대칭 행렬의 역행렬은 대칭 행렬이다:**

$$(S^{-1})^T = (S^T)^{-1} = S^{-1}$$

$\Rightarrow$ $S$가 가역일 때, $S^{-1}$은 대칭이다.

**예제:**

$$S = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}, \quad S^{-1} = \begin{pmatrix} 5 & -2 \\ -2 & 1 \end{pmatrix} = (S^{-1})^T$$

---

### 5.9 대칭 곱과 LDL^T

**대칭 곱 $A^T A$, $AA^T$, $LDL^T$:**

$A \in \mathbb{R}^{m \times n}$에 대해:

- $A^T A \in \mathbb{R}^{n \times n}$: $(A^T A)^T = A^T A \implies$ **대칭**
- $AA^T \in \mathbb{R}^{m \times m}$: $(AA^T)^T = AA^T \implies$ **대칭**

**예제 3:**

$$A = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}$$

$$AA^T = \begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$$

$$A^T A = \begin{pmatrix} -1 & 0 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}$$

**소거에서의 대칭 행렬:** $S^T = S$는 소거를 **두 배 빠르게** 만든다.

$$S = \begin{pmatrix} 1 & 2 \\ 2 & 7 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix} = U$$

$$L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}$$

$$S = LU = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$$

$LU$ 분해에서는 $S$의 대칭성이 보이지 않는다. 대칭 행렬 $S$에 대해, $U$를 $D$와 $L^T$로 더 분해할 수 있다:

$$U = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix} = \begin{pmatrix} d_1 & \\ & d_2 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & \\ & 3 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$$

$$S = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 3 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix} = LDL^T$$

**예제 4: 안장점 행렬 (Saddle-point matrix).**

직사각 행렬 $A \in \mathbb{R}^{m \times n}$에 대해, **안장점 행렬** $S$는 대칭이며 중요하다:

$$S = \begin{pmatrix} I_{m \times m} & A_{m \times n} \\ A^T_{n \times m} & 0_{n \times n} \end{pmatrix} = S^T \quad (m+n) \times (m+n)$$

블록 소거: $R_2 - A^T R_1$:

$$ES = \begin{pmatrix} I & A \\ 0 & -A^T A \end{pmatrix} = U$$

**블록 분해:**

$$S = \begin{pmatrix} I & 0 \\ A^T & I \end{pmatrix}\begin{pmatrix} I & 0 \\ 0 & -A^T A \end{pmatrix}\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} = LDL^T$$

$S$가 가역 $\iff$ $A^T A$가 가역 $\iff$ $\text{rank}(A^T A) = n$ $\iff$ $\mathbf{x} \neq \mathbf{0}$일 때 $A\mathbf{x} \neq \mathbf{0}$ $\iff$ $A$의 열들이 일차독립.

---

<br>

## 6. 도함수와 유한 차분 행렬 (2.5)

### 6.1 테일러 급수와 근사

행렬은 **도함수** (derivatives)를 모방한다. 도함수는 공간의 한 점 $x$에서 또는 시간의 한 순간 $t$에서 무슨 일이 일어나고 있는지 알려준다.

**예제:** $y = x^2 + 2$

$$\frac{dy}{dx} = 2x \xrightarrow{x=1} 2 > 0$$

$$\frac{d^2y}{dx^2} = 2 > 0$$

$y$의 그래프는 위로 굽는다 (기울기가 증가한다).

**$\Delta y = y(x+h) - y(x)$를 고려하자:**

1. 차이는 대략: $\Delta y \approx h \cdot \frac{dy}{dx}\Big|_{x=x_1}$

2. 더 나은 근사: $\Delta y \approx h \cdot \frac{dy}{dx}\Big|_{x=x_1} + \frac{h^2}{2}\frac{d^2y}{dx^2}\Big|_{x=x_1}$ (접선 + 포물선)

3. 정확한 $\Delta y$는 적분이다: $\Delta y = y(x+h) - y(x) = \int_x^{x+h} \frac{dy}{dx} \, dx$

$\Delta y$의 정확도는 도함수 항을 추가함으로써 증가한다.

**테일러 급수 (Taylor Series):**

$$y(x+h) = y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + \cdots + \frac{h^n}{n!}\frac{d^{(n)}y}{dx^n} + \cdots$$

**예제:** $e^x$ (전해석 함수, entire analytic function):

$$e^{x+h} = e^x + h \cdot e^x + \frac{h^2}{2} \cdot e^x + \cdots + \frac{h^n}{n!} e^x + \cdots$$

$$= e^x\left(1 + h + \frac{h^2}{2} + \cdots + \frac{h^n}{n!} + \cdots\right)$$

$$\therefore e^h = 1 + h + \frac{h^2}{2} + \cdots + \frac{h^n}{n!} + \cdots$$

---

### 6.2 차분으로부터의 도함수

**공식 뒤집기: 차분으로부터의 도함수**

접선 포물선에서 시작한다:

$$y(x+h) \approx y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2}$$

**Q:** $y(x)$와 $y(x+h)$를 알면, $\frac{dy}{dx}$를 어떻게 추정하는가?

**Q:** $y(x-h)$를 알면, $\frac{dy}{dx}$, $\frac{d^2y}{dx^2}$를 추정할 수 있는가?

**유한 차분법** (Finite difference method)을 사용하여 도함수를 근사할 수 있다.

**전진 차분 (Forward Difference):**

$$y(x+h) \approx y(x) + h\frac{dy}{dx}$$

$$\implies \frac{dy}{dx} \approx \frac{y(x+h) - y(x)}{h} \quad \text{(전진 차분)}$$

이것은 **1차** 근사이다:

$$y(x+h) = y(x) + h\frac{dy}{dx} + O(h^2)$$

$$\implies \frac{dy}{dx} = \frac{y(x+h) - y(x)}{h} + O(h) \quad \text{(절단 오차)}$$

**후진 차분 (Backward Difference):**

$$y(x-h) \approx y(x) - h\frac{dy}{dx}$$

$$\implies \frac{dy}{dx} \approx \frac{y(x) - y(x-h)}{h}$$

**Q:** $\frac{dy}{dx}$에 대한 근사의 정확도를 높일 수 있는가?

**중심 차분 (Centered Difference) (2차 정확도):**

$$y(x+h) = y(x) + h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + O(h^3)$$

$$y(x-h) = y(x) - h\frac{dy}{dx} + \frac{h^2}{2}\frac{d^2y}{dx^2} + O(h^3)$$

빼면:

$$y(x+h) - y(x-h) = 2h\frac{dy}{dx} + O(h^3)$$

$$\implies \frac{dy}{dx} = \frac{y(x+h) - y(x-h)}{2h} + O(h^2) \approx \frac{y(x+h) - y(x-h)}{2h}$$

이 중심 차분 공식은 **2차 정확도**를 가진다.

**이차 차분 ($\frac{d^2y}{dx^2}$의 근사):**

두 테일러 전개를 더하면:

$$y(x+h) + y(x-h) = 2y(x) + h^2\frac{d^2y}{dx^2} + O(h^4)$$

$$\implies \frac{d^2y}{dx^2} = \frac{y(x+h) - 2y(x) + y(x-h)}{h^2} + O(h^2)$$

$$\approx \frac{y(x+h) - 2y(x) + y(x-h)}{h^2} \quad \text{(이차 차분)}$$

$$\frac{d^2y}{dx^2} \approx \frac{1}{h^2}(1 \quad -2 \quad 1)\begin{pmatrix} y(x-h) \\ y(x) \\ y(x+h) \end{pmatrix}$$

---

### 6.3 이차 차분 행렬 K, T, B

**1차원 영역을 고려하자:**

$$x_0, x_1, x_2, \ldots, x_{N-1}, x_N, x_{N+1}$$

전체 영역을 $N+1$개의 겹치지 않는 요소로 균일 간격 $h$로 분할한다.

여기서 $x_0 = 0$, $x_{N+1} = 1$은 경계점이다.

$$h = \frac{1}{N+1} \implies x_i = ih = \frac{i}{N+1}$$

방정식 $-\frac{d^2u}{dx^2} = f(x)$를 $u(0) = u(1) = 0$ (경계 조건)으로 **이산화** (discretize)한다.

경계 조건에서: $u_0 = u_{N+1} = 0$.

**Q:** $u_1, u_2, \ldots, u_N$은 무엇인가?

$$\left.\frac{d^2u}{dx^2}\right|_{x=x_1} \approx \frac{u_0 - 2u_1 + u_2}{h^2} = -f(x_1)$$

$$\left.\frac{d^2u}{dx^2}\right|_{x=x_2} \approx \frac{u_1 - 2u_2 + u_3}{h^2} = -f(x_2)$$

$$\vdots$$

$$\left.\frac{d^2u}{dx^2}\right|_{x=x_N} \approx \frac{u_{N-1} - 2u_N + u_{N+1}}{h^2} = -f(x_N)$$

행렬 형태:

$$\frac{1}{h^2}\begin{pmatrix} 2 & -1 & 0 & \cdots & 0 \\ -1 & 2 & -1 & \cdots & 0 \\ 0 & -1 & 2 & \cdots & 0 \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & -1 & 2 \end{pmatrix}\begin{pmatrix} u_1 \\ u_2 \\ u_3 \\ \vdots \\ u_{N-1} \\ u_N \end{pmatrix} = \begin{pmatrix} f_1 \\ f_2 \\ f_3 \\ \vdots \\ f_{N-1} \\ f_N \end{pmatrix}$$

여기서 $f_i := f(x_i)$.

$$\boxed{\frac{1}{h^2}K\mathbf{u} = \mathbf{f}}$$

행렬 $K$는 **고정-고정 경계 조건** (fixed-fixed BC): $u_0 = 0$이고 $u_{N+1} = 0$에서 $-\frac{d^2u}{dx^2}$의 자연스러운 근사를 제공한다.

---

### 6.4 K의 성질

$N = 4$로 놓자:

$$K_4 = \begin{pmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}$$

**성질:**

**1. $K$는 대칭이다.** $K^T = K$.

**2. $K$는 띠행렬이다.** $K$의 모든 비영 원소는 주 대각선 주위의 "띠"에 놓인다. 좁은 띠를 가진 행렬은 **희소** (sparse)하다 (대부분 0). 이것은 **삼중대각 행렬** (tridiagonal matrix)이다.

예시: $N = 100$:
- 2의 개수: 100
- $-1$의 개수: $99 + 99 = 198$
- 비영 원소의 수 $\approx 300$. $10000$개 원소 중: $\frac{300}{10000} = 3\%$.

**3. $K$는 상수 대각선을 가진다.** 이는 푸리에 변환, 필터, 합성곱 행렬, 토에플리츠 행렬 (Toeplitz matrix)과 관련된다. $K$는 $(-1, 2, -1)$ 패턴이 각 행에 나타나므로 **이동 불변** (shift-invariant)이다.

**4. $K$는 가역이다.** 소거로 확인할 수 있다. 결과 상삼각 행렬 $U$에 0 피벗이 없으면, $K$는 가역이다.

**5. 대칭 행렬 $K_n$은 양의 정부호 (positive definite)이다:**

$$\mathbf{x}^T K \mathbf{x} > 0 \quad \forall \; \mathbf{x} \neq \mathbf{0}$$

**예제:** $K = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$

$$K\mathbf{x} = \begin{pmatrix} 2x - y \\ -x + 2y \end{pmatrix}$$

$$\mathbf{x}^T(K\mathbf{x}) = 2x^2 - 2xy + 2y^2 = x^2 + (x-y)^2 + y^2 > 0$$

$\mathbf{x}^T K \mathbf{x} \geq 0 \; \forall \; \mathbf{x} \neq \mathbf{0}$일 때, $K$는 **양의 준정부호** (positive semi-definite)라 한다.

**피벗:**
- 가역 행렬은 $n$개의 비영 피벗을 가진다.
- 양의 정부호 대칭 행렬은 $n$개의 **양의** 피벗을 가진다.
- 양의 준정부호 대칭 행렬은 $n$개의 **비음의** 피벗을 가진다.

---

### 6.5 자유-고정 행렬 T

$$-\frac{d^2u}{dx^2} = f(x), \quad \text{여기서 } \frac{du}{dx} = 0 \text{ (} x = 0 \text{에서)}, \quad u(1) = 0$$

$x = 0$에서의 노이만 경계 조건 (Neumann BC): $\frac{u_1 - u_0}{h} = 0 \implies \frac{1}{h^2}(u_0 - u_1) = 0$. 여기서 $u_0$은 **미지수**이다.

이것은 행렬 $T$를 준다:

$$\frac{1}{h^2}\begin{pmatrix} 1 & -1 & 0 & 0 & \cdots & 0 \\ -1 & 2 & -1 & 0 & \cdots & 0 \\ 0 & -1 & 2 & -1 & \cdots & 0 \\ \vdots & & & \ddots & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & -1 & 2 \end{pmatrix}\begin{pmatrix} u_0 \\ u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N \end{pmatrix} = \begin{pmatrix} 0 \\ f_1 \\ f_2 \\ \vdots \\ f_{N-1} \\ f_N \end{pmatrix}$$

$N = 4$일 때:

$$T_4 = \begin{pmatrix} 1 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}$$

**$T_4$의 소거:**

$$\xrightarrow{R_2 + R_1} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix} \xrightarrow{R_3' + R_2'} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix} \xrightarrow{R_4'' + R_3''} \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 1 \end{pmatrix} = U = L^T$$

$$\therefore T_4 = LL^T$$

$$L = \begin{pmatrix} 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 0 & -1 & 1 \end{pmatrix}$$

**$L^{-1}$은 무엇인가?** 가우스-조르단 방법 사용: $(L|I) \Rightarrow (I|L^{-1})$

$$\begin{pmatrix} 1 & 0 & 0 & 0 & | & 1 & 0 & 0 & 0 \\ -1 & 1 & 0 & 0 & | & 0 & 1 & 0 & 0 \\ 0 & -1 & 1 & 0 & | & 0 & 0 & 1 & 0 \\ 0 & 0 & -1 & 1 & | & 0 & 0 & 0 & 1 \end{pmatrix} \implies \begin{pmatrix} 1 & 0 & 0 & 0 & | & 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & | & 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 & | & 1 & 1 & 1 & 0 \\ 0 & 0 & 0 & 1 & | & 1 & 1 & 1 & 1 \end{pmatrix}$$

$$L^{-1} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 1 & 1 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix}$$

$$T_4^{-1} = (LL^T)^{-1} = L^{-T}L^{-1} = \begin{pmatrix} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 1 & 1 & 1 & 0 \\ 1 & 1 & 1 & 1 \end{pmatrix} = \begin{pmatrix} 4 & 3 & 2 & 1 \\ 3 & 3 & 2 & 1 \\ 2 & 2 & 2 & 1 \\ 1 & 1 & 1 & 1 \end{pmatrix}$$

---

### 6.6 자유-자유 행렬 B

**자유-자유 행렬 $B$는 특이 (singular)하다.**

- 가역이 아니다
- 영이 아닌 $\mathbf{x}$가 존재하여 $B\mathbf{x} = \mathbf{0}$

방정식 $-\frac{d^2u}{dx^2} = f(x)$에 대해:

$$\frac{du}{dx} = 0 \text{ (} x = 0 \text{에서)}, \quad \frac{du}{dx} = 0 \text{ (} x = 1 \text{에서)}$$

$$\frac{u_1 - u_0}{h} = 0 \implies \frac{1}{h^2}(u_0 - u_1) = 0 \quad (u_0 \text{는 미지수})$$

$$\frac{u_{N+1} - u_N}{h} = 0 \implies \frac{1}{h^2}(-u_N + u_{N+1}) = 0 \quad (u_{N+1} \text{는 미지수})$$

$$\frac{1}{h^2}\begin{pmatrix} 1 & -1 & 0 & 0 & \cdots & 0 & 0 \\ -1 & 2 & -1 & 0 & \cdots & 0 & 0 \\ 0 & -1 & 2 & -1 & \cdots & 0 & 0 \\ \vdots & & & \ddots & & & \vdots \\ 0 & 0 & \cdots & -1 & 2 & -1 & 0 \\ 0 & 0 & \cdots & 0 & -1 & 2 & -1 \\ 0 & 0 & \cdots & 0 & 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} u_0 \\ u_1 \\ u_2 \\ \vdots \\ u_{N-1} \\ u_N \\ u_{N+1} \end{pmatrix} = \begin{pmatrix} 0 \\ f_1 \\ f_2 \\ \vdots \\ f_{N-1} \\ f_N \\ 0 \end{pmatrix}$$

**$B_3$을 고려하자:**

$$B_3 = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}$$

$$B_3 \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$\mathbf{x} = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$은 $B_3$의 **영공간** (null space)에 속한다.

$\mathbf{x} = c\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$ 벡터들은 $B_3$의 영공간에 속한다.

$\Rightarrow$ $B_3$은 **가역이 될 수 없다.**

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:-----|:-----------|
| $A\mathbf{x} = \mathbf{b}$ | $n$개의 미지수를 가진 $n$개의 연립일차방정식 |
| 소거법 (Elimination) | 피벗 아래에 0을 만들기 위해 행 $j$의 $l_{ij}$배를 행 $i$에서 뺌 |
| 후진 대입 (Back substitution) | $U\mathbf{x} = \mathbf{c}$를 마지막 행부터 위로 풂 |
| 세 가지 경우 | 유일한 해 ($\text{rank} = n$), 해 없음, 무한히 많은 해 |
| 확대 행렬 (Augmented matrix) | $(A \mid \mathbf{b})$는 소거 중 양변을 추적 |
| 동차 연립방정식 (Homogeneous system) | $A\mathbf{x} = \mathbf{0}$; $\text{rank}(A) < n$일 때 비자명 해 존재 |
| 피벗 (Pivots) | $U$의 대각 원소; 가역성을 위해 영이 아니어야 함 |
| 행 교환 (Row exchange) | 피벗 위치에 0이 나타날 때 행을 교환 |
| 소거 행렬 $E_{ij}$ (Elimination matrix) | $(i,j)$ 위치에 $-l_{ij}$를 가진 단위 행렬 |
| $EA = U$ | 모든 소거 행렬의 곱이 $A$를 상삼각 $U$로 변환 |
| $A = LU$ | $L = E^{-1}$은 승수 $l_{ij}$가 올바른 위치에 있는 하삼각 행렬 |
| 역행렬 $A^{-1}$ (Inverse) | $A$가 $n$개의 독립 열을 가질 때 존재 ($\text{rank}(A) = n$) |
| 역행렬의 유일성 | 좌역행렬은 우역행렬과 같음; $BA = I$이고 $AC = I \implies B = C$ |
| $(AB)^{-1} = B^{-1}A^{-1}$ | 역행렬은 역순으로 나옴 |
| $2 \times 2$ 역행렬 | $A^{-1} = \frac{1}{ad-bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$; $\det(A) \neq 0$ 필요 |
| 가우스-조르단 (Gauss-Jordan) | $(A \mid I) \Rightarrow (I \mid A^{-1})$로 역행렬을 명시적으로 구함 |
| 소거의 비용 | $A \to U$: $\frac{1}{3}n^3$; $\mathbf{b} \to \mathbf{c} \to \mathbf{x}$: $n^2$ |
| $A = LU$의 두 번째 증명 | 열 곱하기 행: $A = \sum \mathbf{l}_k \mathbf{u}_k$ |
| 행 교환 없는 $A = LU$ | 모든 좌상단 $k \times k$ 부분행렬이 가역이어야 함 |
| 치환 행렬 $P$ (Permutation matrix) | $I$의 행을 다른 순서로; $P^{-1} = P^T$; $n!$개 가능 |
| $PA = LU$ | 행 교환이 $P$에 기록된 일반 분해 |
| 부분 피벗팅 (Partial pivoting) | 가장 큰 원소를 피벗으로 만들기 위한 행 교환; 반올림 오차 감소 |
| $PAQ$ | 행 치환 $P$, 열 치환 $Q$; $C(A) = C(PA)$ |
| 전치 $A^T$ (Transpose) | $(A^T)_{ij} = A_{ji}$; $A^T$의 열은 $A$의 행 |
| 전치 규칙 | $(A+B)^T = A^T + B^T$; $(AB)^T = B^T A^T$; $(A^{-1})^T = (A^T)^{-1}$ |
| 내적 (Inner product) | $\mathbf{x} \cdot \mathbf{y} = \mathbf{x}^T\mathbf{y}$; $(A\mathbf{x}) \cdot \mathbf{y} = \mathbf{x} \cdot (A^T\mathbf{y})$ |
| 외적 (Outer product) | $\mathbf{x}\mathbf{y}^T$는 랭크 1인 $n \times n$ 행렬 |
| 대칭 행렬 $S$ (Symmetric matrix) | $S^T = S$; $S^{-1}$도 대칭 |
| $A^T A$와 $AA^T$ | 항상 대칭; $A$의 열이 독립일 때 $A^T A$가 가역 |
| $S = LDL^T$ | 대칭 분해; $LU$가 보여주지 못하는 대칭성을 드러냄 |
| 부분적분과 전치 (IBP and transpose) | $A = d/dt \implies A^T = -d/dt$ (반대칭); $(Ax, y) = (x, A^T y)$ |
| 전진 차분 (Forward difference) | $\frac{dy}{dx} \approx \frac{y(x+h)-y(x)}{h}$; 1차 정확도 $O(h)$ |
| 중심 차분 (Centered difference) | $\frac{dy}{dx} \approx \frac{y(x+h)-y(x-h)}{2h}$; 2차 정확도 $O(h^2)$ |
| 이차 차분 (Second difference) | $\frac{d^2y}{dx^2} \approx \frac{y(x+h)-2y(x)+y(x-h)}{h^2}$; 2차 정확도 |
| 행렬 $K$ (고정-고정) | 삼중대각 $(-1, 2, -1)$; 대칭, 띠, 가역, 양의 정부호 |
| 행렬 $T$ (자유-고정) | 첫 행이 $(1, -1, 0, \ldots)$; $T = LL^T$; 가역 |
| 행렬 $B$ (자유-자유) | 특이; $B\mathbf{1} = \mathbf{0}$; 상수 벡터가 영공간에 속함 |
| 양의 정부호 (Positive definite) | $\mathbf{x}^T K\mathbf{x} > 0$, 모든 $\mathbf{x} \neq \mathbf{0}$에 대해; 모든 피벗이 양수 |
| 양의 준정부호 (Positive semi-definite) | $\mathbf{x}^T K\mathbf{x} \geq 0$, 모든 $\mathbf{x} \neq \mathbf{0}$에 대해; 모든 피벗이 비음수 |

---
