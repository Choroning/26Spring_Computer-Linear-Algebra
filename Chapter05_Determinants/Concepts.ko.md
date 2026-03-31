# 제5장 강의 — 행렬식

> **최종 수정일:** 2026-03-31

---

<br>

## 목차

- [1. 행렬식 (개요)](#1-행렬식-개요)
- [2. 3x3 행렬식과 여인수 (N5-1)](#2-3x3-행렬식과-여인수-n5-1)
  - [2.1 행렬식의 정의](#21-행렬식의-정의)
  - [2.2 2x2 행렬식](#22-2x2-행렬식)
  - [2.3 행 교환은 부호를 바꾼다](#23-행-교환은-부호를-바꾼다)
  - [2.4 각 행에 대한 행렬식의 선형성](#24-각-행에-대한-행렬식의-선형성)
  - [2.5 3x3 행렬식](#25-3x3-행렬식)
  - [2.6 3x3 행렬식의 6개 항 (순열 접근법)](#26-3x3-행렬식의-6개-항-순열-접근법)
  - [2.7 여인수와 역행렬 공식](#27-여인수와-역행렬-공식)
  - [2.8 행 i를 따른 여인수 전개 공식](#28-행-i를-따른-여인수-전개-공식)
  - [2.9 여인수를 이용한 역행렬 공식](#29-여인수를-이용한-역행렬-공식)
  - [2.10 예제: 여인수를 이용한 3x3 역행렬](#210-예제-여인수를-이용한-3x3-역행렬)
  - [2.11 4x4 행렬식](#211-4x4-행렬식)
  - [2.12 스칼라 곱과 행렬식](#212-스칼라-곱과-행렬식)
- [3. 행렬식의 계산과 활용 (N5-2)](#3-행렬식의-계산과-활용-n5-2)
  - [3.1 유용한 성질들](#31-유용한-성질들)
  - [3.2 삼각행렬과 대각행렬의 행렬식](#32-삼각행렬과-대각행렬의-행렬식)
  - [3.3 det(A^T) = det(A)](#33-detat--deta)
  - [3.4 소거 행렬과 행렬식](#34-소거-행렬과-행렬식)
  - [3.5 스케일링 행렬](#35-스케일링-행렬)
  - [3.6 기본 행렬](#36-기본-행렬)
  - [3.7 곱의 법칙: det(AB) = det(A) det(B)](#37-곱의-법칙-detab--deta-detb)
  - [3.8 직교 행렬](#38-직교-행렬)
  - [3.9 피벗을 이용한 행렬식](#39-피벗을-이용한-행렬식)
  - [3.10 역행렬의 행렬식](#310-역행렬의-행렬식)
  - [3.11 det(A+B)는 det(A) + det(B)가 아니다](#311-detab는-deta--detb가-아니다)
  - [3.12 단일 행에 대한 선형성](#312-단일-행에-대한-선형성)
  - [3.13 det(A) = det(A^T)의 증명](#313-deta--detat의-증명)
  - [3.14 크래머 법칙](#314-크래머-법칙)
- [4. 행렬식에 의한 넓이와 부피 (N5-3)](#4-행렬식에-의한-넓이와-부피-n5-3)
  - [4.1 2D에서의 평행사변형](#41-2d에서의-평행사변형)
  - [4.2 평행사변형의 넓이](#42-평행사변형의-넓이)
  - [4.3 3D에서의 기울어진 상자](#43-3d에서의-기울어진-상자)
  - [4.4 삼각형의 넓이](#44-삼각형의-넓이)
  - [4.5 외적 (Cross Product)](#45-외적-cross-product)
  - [4.6 평행사변형 넓이의 기하학적 증명](#46-평행사변형-넓이의-기하학적-증명)
  - [4.7 상자의 부피](#47-상자의-부피)
  - [4.8 예제: 단위 정육면체와 축 정렬 상자](#48-예제-단위-정육면체와-축-정렬-상자)
- [요약](#요약)

---

<br>

## 1. 행렬식 (개요)

### 장 구성

- **5.1** &mdash; $3 \times 3$ 행렬식과 여인수 (cofactor)
- **5.2** &mdash; 행렬식의 계산과 활용
- **5.3** &mdash; 행렬식에 의한 넓이와 부피 (기하학과 선형대수의 연결)

### 핵심 개념

$A$를 정사각행렬이라 하자.

- $\det(A) = 0$이면, $A$는 **역행렬이 존재하지 않는다** (not invertible)
  - $\Longleftrightarrow$ $A$는 **특이행렬** (singular)이다
  - $\Longleftrightarrow$ $A$는 **종속인 열** (dependent columns)을 가진다

### 유용한 공식 (미리보기)

1. $PA = LU$

$$\det(PA) = \det(P)\det(A) = \pm\det(A)$$

$$\det(LU) = \underbrace{\det(L)}_{=1}\det(U) = \text{U에서 피벗들의 곱}$$

2. 여인수를 이용한 역행렬 공식:

$$(A^{-1})_{ij} = \frac{1}{\det(A)}\bigl(\text{cofactor of } A_{ij}\bigr)^T$$

---

<br>

## 2. 3x3 행렬식과 여인수 (N5-1)

### 2.1 행렬식의 정의

**행렬식** (determinant)은 정사각행렬에 대응하는 스칼라 값으로, $\det(A)$ 또는 $|A|$로 표기한다.

### 2.2 2x2 행렬식

**공식:** $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$의 행렬식은:

$$\det(A) = ad - bc$$

**예제:**

$$\det\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = 1, \qquad \det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc$$

$$\det\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = -1, \qquad \det\begin{pmatrix} c & d \\ a & b \end{pmatrix} = bc - ad$$

$$\det\begin{pmatrix} b & a \\ d & c \end{pmatrix} = bc - ad$$

**관찰:** 두 행을 교환하면 행렬식의 부호가 바뀐다 (두 열을 교환하는 경우에도 마찬가지).

**특이행렬의 경우:**

$$\det\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = 0$$

$$\det\begin{pmatrix} a & b \\ 2a & 2b \end{pmatrix} = 2ab - 2ab = 0$$

$ax + by = e$의 기울기가 $-a/b$이고 $2ax + 2by = f$의 기울기도 $-a/b$ (같은 기울기)이면, 두 직선은 평행하다.

$\det A = 0$은 $A$의 열들이 **선형독립이 아님** (NOT linearly independent)을 의미하며, 즉 $N(A) \neq \{\mathbf{0}\}$이다.

**특이행렬 예제:**

$$\text{특이행렬 } \begin{pmatrix} a & 2a \\ c & 2c \end{pmatrix} \text{의 행렬식} = 0$$

### 2.3 행 교환은 부호를 바꾼다

**성질:** 행 (또는 열) 교환은 행렬식의 부호를 바꾼다.

$$PA = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \begin{pmatrix} c & d \\ a & b \end{pmatrix}$$

$$\det(PA) = bc - ad = -\det(A)$$

### 2.4 각 행에 대한 행렬식의 선형성

$\begin{pmatrix} xa + ye & xb + yf \\ c & d \end{pmatrix}$의 행렬식은:

$$d(xa + ye) - c(xb + yf) = x(ad - bc) + y(de - cf)$$

이는 첫째 행에 대한 선형성을 보여준다.

**비고:** $3 \times 3$ 행렬식은 $3! = 6$개의 항을 가지며, 이들은 3개의 항(여인수)으로 분류된다. 이 여인수들을 $\det(A)$로 나누면 $A^{-1}$을 얻는다.

### 2.5 3x3 행렬식

$3 \times 3$ 행렬의 행렬식은 $3! = 6$개의 항을 가지며, 각 행(과 각 열)에서 하나의 원소를 선택하여 $3 \times 2 = 6$개의 항을 만든다.

**순열 행렬(permutation matrix)과 그 행렬식:**

$I_{3 \times 3}$에서 시작:

$$\det\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = 1 \xrightarrow{C_1 \leftrightarrow C_2} \det\begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix} = -1$$

$$\xrightarrow{R_2 \leftrightarrow R_3} \det\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix} = 1 \xrightarrow{C_2 \leftrightarrow C_3} \det\begin{pmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{pmatrix} = -1$$

$$\xrightarrow{R_1 \leftrightarrow R_2} \det\begin{pmatrix} 0 & 0 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = 1 \xrightarrow{R_1 \leftrightarrow R_2} \det\begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix} = -1$$

**핵심 성질:** 행 (또는 열) 교환은 $\det(A)$에 $-1$을 곱한다.

### 2.6 3x3 행렬식의 6개 항 (순열 접근법)

$A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}$에서, 6개의 순열 행렬이 각각 하나의 항을 선택한다:

| 순열 행렬 | $\det(P)$ | 선택된 원소 | 항 |
|:---|:---:|:---|:---:|
| $\begin{pmatrix} 1&0&0\\0&1&0\\0&0&1 \end{pmatrix}$ | $+1$ | $a, q, z$ | $+aqz$ |
| $\begin{pmatrix} 0&1&0\\1&0&0\\0&0&1 \end{pmatrix}$ | $-1$ | $b, p, z$ | $-bpz$ |
| $\begin{pmatrix} 0&1&0\\0&0&1\\1&0&0 \end{pmatrix}$ | $+1$ | $b, r, x$ | $+brx$ |
| $\begin{pmatrix} 0&0&1\\0&1&0\\1&0&0 \end{pmatrix}$ | $-1$ | $c, q, x$ | $-cqx$ |
| $\begin{pmatrix} 0&0&1\\1&0&0\\0&1&0 \end{pmatrix}$ | $+1$ | $c, p, y$ | $+cpy$ |
| $\begin{pmatrix} 1&0&0\\0&0&1\\0&1&0 \end{pmatrix}$ | $-1$ | $a, r, y$ | $-ary$ |

**대각선 방법 (사뤼스 법칙, Sarrus' Rule):** 행렬 옆에 처음 두 열을 다시 적고, 대각선을 따라 곱을 취한다 (아래로 = 양, 위로 = 음).

**6개 항의 결합:**

$$\det A = aqz - bpz + brx - cqx + cpy - ary$$

**첫째 행 원소** $a, b, c$**로 모으기:**

$$\det A = a(qz - ry) - b(pz - rx) + c(py - qx)$$

**공식에 사용되는 부분행렬의 행렬식:**

$$A_1 = (2), \quad \det(A_1) = 2$$

$$A_2 = \begin{pmatrix} 2 & 7 \\ -1 & 2 \end{pmatrix}, \quad \det(A_2) = 4 - (-7) = 3$$

$$A_3 = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}, \quad \det(A_3) = 8 - 2 - 2 = 4$$

### 2.7 여인수와 역행렬 공식

$3 \times 3$ 행렬 $A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}$에 대해:

$$\det(A) = a(qz - ry) - b(pz - rx) + c(py - qx)$$

$\det(A)$는 $2 \times 2$ 행렬식을 이용하여 계산됨을 관찰할 수 있다:

$$\det\begin{pmatrix} q & r \\ y & z \end{pmatrix}, \quad \det\begin{pmatrix} p & r \\ x & z \end{pmatrix}, \quad \det\begin{pmatrix} p & q \\ x & y \end{pmatrix}$$

따라서 행렬식은 다음과 같이 쓸 수 있다:

$$\det(A) = a \det\begin{pmatrix} q & r \\ y & z \end{pmatrix} + b\left(-\det\begin{pmatrix} p & r \\ x & z \end{pmatrix}\right) + c \det\begin{pmatrix} p & q \\ x & y \end{pmatrix}$$

여기서 $a, b, c$는 **인수** (1행의 원소)이고, $2 \times 2$ 행렬식은 $a, b, c$에 대한 **여인수** (cofactor)이다:

| 원소 | 여인수 |
|:------|:---------|
| $A_{11} = a$ | $C_{11} = +\det\begin{pmatrix} q & r \\ y & z \end{pmatrix}$ |
| $A_{12} = b$ | $C_{12} = -\det\begin{pmatrix} p & r \\ x & z \end{pmatrix}$ |
| $A_{13} = c$ | $C_{13} = +\det\begin{pmatrix} p & q \\ x & y \end{pmatrix}$ |

### 2.8 행 i를 따른 여인수 전개 공식

$(i,j)$ 여인수 $C_{ij}$를 구하려면: **$A$에서 행 $i$와 열 $j$를 제거**한다.

$$C_{ij} = (-1)^{i+j} \det\bigl(\text{남은 } (n{-}1) \times (n{-}1) \text{ 행렬}\bigr)$$

**행 $i$를 따른 여인수 전개 공식:**

$$\det A = A_{i1}C_{i1} + A_{i2}C_{i2} + \cdots + A_{in}C_{in}$$

$C_{ij}$는 $\det(A)$에서 $A_{ij}$가 곱해지는 모든 항의 모음이다.

**예제: $2 \times 2$ 행렬의 여인수**

$A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$가 주어지면, 여인수는:

$$A_{11} = a \to C_{11} = d, \quad A_{12} = b \to C_{12} = -c$$

$$A_{21} = c \to C_{21} = -b, \quad A_{22} = d \to C_{22} = a$$

$A$의 **여인수 행렬** (cofactor matrix)은:

$$C = \begin{pmatrix} d & -c \\ -b & a \end{pmatrix}$$

### 2.9 여인수를 이용한 역행렬 공식

$2 \times 2$ 경우의 검증:

$$AC^T = \begin{pmatrix} a & b \\ c & d \end{pmatrix}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix} = \begin{pmatrix} ad - bc & 0 \\ 0 & ad - bc \end{pmatrix} = (ad - bc)\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = \det(A) \cdot I$$

따라서:

$$A \cdot \frac{C^T}{\det(A)} = I$$

$$\boxed{A^{-1} = \frac{1}{\det(A)} C^T}$$

**중요 사항:**

- $A^{-1}$이 존재하려면, $\det(A)$는 **0이 아니어야** 한다.
- $A^{-1}$의 모든 원소는 두 행렬식의 비(ratio)이다:

$$(A^{-1})_{ij} = \frac{(C^T)_{ij}}{\det(A)} = \frac{C_{ji}}{\det(A)}$$

- $C^T = \text{adj}(A)$는 $A$의 **수반 행렬** (adjugate matrix)이다.

**일반 공식:**

$$AC^T = \begin{pmatrix} \det(A) & & 0 \\ & \det(A) & \\ 0 & & \ddots & \det(A) \end{pmatrix} = \det(A) \cdot I$$

### 2.10 예제: 여인수를 이용한 3x3 역행렬

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}, \quad \det(A) = 1$$

**행 축소** $(A | I) \to (I | A^{-1})$ 사용:

$$\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 1 & 1 & | & 0 & 1 & 0 \\ 0 & 0 & 1 & | & 0 & 0 & 1 \end{pmatrix} \xrightarrow{R_1 - R_2,\; R_1 - R_3} \begin{pmatrix} 1 & 0 & 0 & | & 1 & -1 & 0 \\ 0 & 1 & 0 & | & 0 & 1 & -1 \\ 0 & 0 & 1 & | & 0 & 0 & 1 \end{pmatrix}$$

$$A^{-1} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}$$

**여인수 공식 사용:**

$$A^{-1} = \frac{1}{\det(A)} C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T$$

모든 여인수 계산:

$$C_{11} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1, \quad C_{12} = -\det\begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix} = 0, \quad C_{13} = \det\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = 0$$

$$C_{21} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = -1, \quad C_{22} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1, \quad C_{23} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = 0$$

$$C_{31} = \det\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = 0, \quad C_{32} = -\det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = -1, \quad C_{33} = \det\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = 1$$

### 2.11 4x4 행렬식

$A_4 \in \mathbb{R}^{4 \times 4}$는 $4! = 4 \cdot 3 \cdot 2 \cdot 1 = 24$개의 항을 가진다.

**예제:**

$$A_4 = \begin{pmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{pmatrix}$$

1행을 따른 여인수 전개:

$$\det(A_4) = 2 \det\begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix} + (-1)(-1)\det\begin{pmatrix} -1 & -1 & 0 \\ 0 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}$$

$$= 2 \cdot 4 + 1 \cdot (-3) = 5$$

- 값 $2$는 $(A_4)_{11}$에서 나온다
- 값 $-1$은 $(A_4)_{12}$에서 나온다

**부호 규칙:** $(-1)^{i+j}$를 곱한다: 원소 $A_{ij}$에 대해 $i + j$가 홀수이면 음부호가 붙는다.

**$2 \times 2$ 경우:**

$$\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = \det\begin{pmatrix} a & \\ & d \end{pmatrix} + \det\begin{pmatrix} & b \\ c & \end{pmatrix} = ad - bc$$

### 2.12 스칼라 곱과 행렬식

**$2 \times 2$ 경우:**

$$A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}, \quad \det(A) = ad - bc$$

$$2A = \begin{pmatrix} 2a & 2b \\ 2c & 2d \end{pmatrix}, \quad \det(2A) = (2a)(2d) - (2b)(2c) = 2^2(ad - bc) = 4\det(A)$$

$$3A = \begin{pmatrix} 3a & 3b \\ 3c & 3d \end{pmatrix}, \quad \det(3A) = 3^2(ad - bc) = 9\det(A)$$

$$\alpha A = \begin{pmatrix} \alpha a & \alpha b \\ \alpha c & \alpha d \end{pmatrix}, \quad \det(\alpha A) = \alpha^2 \det(A)$$

**$3 \times 3$ 경우:**

$$A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}, \quad \det(A) = a(qz - ry) - b(pz - rx) + c(py - qx)$$

$$\det(2A) = 2a(2q \cdot 2z - 2r \cdot 2y) - 2b(2p \cdot 2z - 2r \cdot 2x) + 2c(2p \cdot 2y - 2q \cdot 2x) = 2^3 \det(A)$$

$$\det(\alpha A) = \alpha^3 \det(A)$$

**$A \in \mathbb{R}^{n \times n}$에 대한 일반 공식:**

$$\boxed{\det(\alpha A) = \alpha^n \det(A)}$$

---

<br>

## 3. 행렬식의 계산과 활용 (N5-2)

### 3.1 유용한 성질들

1. $\det(A^T) = \det(A)$
2. $\det(AB) = \det(A)\det(B)$
3. 직교 행렬 $Q$에 대해 $|\det(Q)| = 1$
4. 소거 행렬은 $\det(E) = 1$이므로, $\det(EA) = \det(E)\det(A) = \det(A)$

   예) $E = \begin{pmatrix} 1 & & \\ a & 1 & \\ b & c & 1 \end{pmatrix}$, $|E| = 1$

5. **크래머 법칙** (Cramer's Rule)은 행렬식의 비를 이용하여 $\mathbf{x} = A^{-1}\mathbf{b}$를 구한다 (느린 방법)
6. $\det(A) = \pm$ $A = LU$에서 피벗들의 곱
7. **큰 공식** (big formula)은 $n!$개의 순열로부터 $n!$개의 항을 가진다 ($n > 3$이면 매우 느림)

### 행렬식이 알려주는 것

1. 가역 행렬은 $\det(A) \neq 0$이다
2. 특이 행렬은 $\det(A) = 0$이다
3. $\lambda$를 고유값 (eigenvalue), $\mathbf{x}$를 고유벡터 (eigenvector)라 하면, $A\mathbf{x} = \lambda \mathbf{x}$에서

$$\Rightarrow (A - \lambda I)\mathbf{x} = \mathbf{0}$$

$\mathbf{x}$가 영벡터가 아니므로, $A - \lambda I$는 특이행렬이고:

$$\det(A - \lambda I) = 0$$

이것이 $\lambda$에 대한 방정식이 된다.

### 3.2 삼각행렬과 대각행렬의 행렬식

**상삼각행렬 (upper triangular):**

$$\det\begin{pmatrix} a & b & c \\ 0 & q & r \\ 0 & 0 & z \end{pmatrix} = aqz$$

**대각행렬 (diagonal):**

$$\det\begin{pmatrix} a & & \\ & q & \\ & & z \end{pmatrix} = aqz$$

대각 원소들을 곱하기만 하면 행렬식을 구할 수 있다.

### 3.3 det(A^T) = det(A)

$A$를 전치하면, 행렬식 공식은 같은 결과를 준다.

**예제:** $A \in \mathbb{R}^{3 \times 3}$

$$A = \begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix} \Rightarrow \det(A) = aqz + brx + cpy - cqx - ary - bpz$$

$$= a(qz - ry) - b(pz - rx) + c(py - qx)$$

$$A^T = \begin{pmatrix} a & p & x \\ b & q & y \\ c & r & z \end{pmatrix} \Rightarrow \det(A^T) = aqz + pyc + xbr - xqc - ayr - pbz$$

$$= a(qz - yr) - b(pz - xr) + c(py - xq)$$

두 식이 동일하므로, $\det(A^T) = \det(A)$가 확인된다.

### 3.4 소거 행렬과 행렬식

**예제:**

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix}$$

$R_2 - 2R_1$ 적용:

$$\begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 4 & 0 \end{pmatrix} = \underbrace{\begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}}_{E_{21}} \underbrace{\begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix}}_{A}$$

행렬식 계산:

$$\det(A) = 8 - 20 = -12$$

$$\det(E_{21}A) = -12$$

$$\det(E_{21}) = 1$$

$\det(E_{21}A) = \det(E_{21})\det(A) = \det(A) = -12$임을 관찰할 수 있다.

**$U$까지의 완전 소거:**

$$A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 4 & 5 \\ 0 & 4 & 0 \end{pmatrix} \xrightarrow{R_2 - 2R_1} \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 4 & 0 \end{pmatrix} \xrightarrow{R_3 - 2R_2} \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 3 \\ 0 & 0 & -6 \end{pmatrix} = U$$

$$E_{21} = \begin{pmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{32} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -2 & 1 \end{pmatrix}$$

$$\det(E_{21}) = 1, \quad \det(E_{32}) = 1, \quad \det(U) = -12$$

$$\det(U) = \det(E_{32}E_{21}A) = \det(E_{32})\det(E_{21}A) = \det(E_{32})\det(E_{21})\det(A) = \det(A)$$

**$U^{-1}$를 구하기 위한 역대입:**

$(U|I)$에서 시작:

$$\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 2 & 3 & | & 0 & 1 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}$$

1단계: $R_2' = R_2/2$, $S_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1/2 & 0 \\ 0 & 0 & 1 \end{pmatrix}$ 사용

$$\begin{pmatrix} 1 & 1 & 1 & | & 1 & 0 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}$$

2단계: $R_1 - R_2$, $E_{12} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$ 사용

$$\begin{pmatrix} 1 & 0 & -1/2 & | & 1 & -1/2 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & -6 & | & 0 & 0 & 1 \end{pmatrix}$$

3단계: $R_3' = -R_3/6$, $S_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & -1/6 \end{pmatrix}$ 사용

$$\begin{pmatrix} 1 & 0 & -1/2 & | & 1 & -1/2 & 0 \\ 0 & 1 & 3/2 & | & 0 & 1/2 & 0 \\ 0 & 0 & 1 & | & 0 & 0 & -1/6 \end{pmatrix}$$

4단계: $R_1' = R_1 + \frac{1}{2}R_3$ 및 $R_2' = R_2 - \frac{3}{2}R_3$, $E_{13}$과 $E_{23}$ 사용:

$$E_{13} = \begin{pmatrix} 1 & 0 & 1/2 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}, \quad E_{23} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -3/2 \\ 0 & 0 & 1 \end{pmatrix}$$

$$\begin{pmatrix} 1 & 0 & 0 & | & 1 & -1/2 & -1/12 \\ 0 & 1 & 0 & | & 0 & 1/2 & +1/4 \\ 0 & 0 & 1 & | & 0 & 0 & -1/6 \end{pmatrix} = (I | U^{-1})$$

따라서: $E_{23}E_{13}S_3E_{12}S_2 \cdot U = I$

$$U = S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1} \cdot I$$

여기서 $S_2^{-1} = \begin{pmatrix} 1 & \\ & 2 \\ & & 1 \end{pmatrix}$, $S_3^{-1} = \begin{pmatrix} 1 & \\ & 1 \\ & & -6 \end{pmatrix}$, $E_{12}^{-1} = \begin{pmatrix} 1 & 1 \\ & 1 \\ & & 1 \end{pmatrix}$, $E_{13}^{-1} = \begin{pmatrix} 1 & & -1/2 \\ & 1 & \\ & & 1 \end{pmatrix}$, $E_{23}^{-1} = \begin{pmatrix} 1 & & \\ & 1 & 3/2 \\ & & 1 \end{pmatrix}$

**$\det(U)$ 계산:**

$$\det(U) = \det(S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1}) = 2 \cdot (-6) \det(E_{32}E_{21} I E_{12}^{-1} I E_{13}^{-1}E_{23}^{-1}) = -12$$

**핵심 결과:**

$$\det(A) = \det(E_{32}^{-1}E_{21}^{-1} U) = \det(E_{32}^{-1}E_{21}^{-1} \cdot S_2^{-1}E_{12}^{-1}S_3^{-1}E_{13}^{-1}E_{23}^{-1} I)$$

$$= 2 \cdot (-6) \det(\cdots) = -12$$

### 3.5 스케일링 행렬

$$I = \begin{pmatrix} 1 & \\ & 1 \end{pmatrix}, \quad \det(I) = 1$$

**$x$-방향 스케일링:**

$$S_x = \begin{pmatrix} 2 & \\ & 1 \end{pmatrix}, \quad \det(S_x) = 2 = 2\det(I)$$

$$\begin{pmatrix} 2 & \\ & 1 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} 2x \\ y \end{pmatrix}$$

단위 정사각형이 $x$-방향으로 폭 2로 늘어난다.

**$y$-방향 스케일링:**

$$S_y = \begin{pmatrix} 1 & \\ & 2 \end{pmatrix}, \quad \det(S_y) = 2 = 2\det(I)$$

$$\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} x \\ 2y \end{pmatrix}$$

단위 정사각형이 $y$-방향으로 높이 2로 늘어난다.

**일반 스케일링 행렬:**

$$S = \begin{pmatrix} a & & \\ & b & \\ & & c \end{pmatrix} \Rightarrow \det(S) = abc = abc \cdot \det(I)$$

### 3.6 기본 행렬

$A$가 역행렬 $A^{-1}$을 가질 때, $A$는 **기본 행렬** (elementary matrix)들의 곱으로 분해할 수 있다 (소거 행렬, 순열 행렬, 스케일링 행렬).

**정의:** **기본 행렬**이란 단위 행렬에 하나의 기본 행 연산을 수행하여 얻은 정사각행렬이다. 소거 행렬, 순열 행렬, 스케일링 행렬 등이 이에 해당한다.

### 3.7 곱의 법칙: det(AB) = det(A) det(B)

$$\det(AB) = \det(A)\det(B)$$

**증명.** $A, B \in \mathbb{R}^{n \times n}$이라 하자.

**경우 i)** $A$ 또는 $B$가 가역이 아닌 경우.

- $B$가 가역이 아니라 하자. 그러면 $B\mathbf{x} = \mathbf{0}$을 만족하는 영이 아닌 벡터 $\mathbf{x}$가 존재한다.
- 양변에 $A$를 곱하면: $AB\mathbf{x} = \mathbf{0}$.
- 이는 $AB$가 가역이 아님을 의미한다.
- 따라서: $\det(AB) = 0$이고 $\det(A)\det(B) = 0$.

**경우 ii)** $A, B$ 모두 가역인 경우.

$$\text{rref}(A) = I \Rightarrow A = E_p E_{p-1} \cdots E_2 E_1 I$$

$$AB = E_p E_{p-1} \cdots E_2 E_1 B$$

$$\det(AB) = \det(E_p)\det(E_{p-1} \cdots E_2 E_1 B)$$

$$= \det(E_p)\det(E_{p-1})\det(E_{p-2} \cdots E_2 E_1 B)$$

$$= \cdots$$

$$= \det(E_p)\det(E_{p-1}) \cdots \det(E_2)\det(E_1)\det(B)$$

$$= \underbrace{\det(E_p E_{p-1} \cdots E_2 E_1)}_{\det(A)} \det(B)$$

$$= \det(A)\det(B) \qquad \square$$

### 3.8 직교 행렬

직교 행렬 $Q$는 행렬식이 $\pm 1$이다.

$$Q^T Q = I$$

$$\det(Q^T Q) = \det(I) = 1$$

$$\det(Q^T)\det(Q) = 1$$

$$(\det Q)^2 = 1$$

$$\therefore \det(Q) = \pm 1$$

### 3.9 피벗을 이용한 행렬식

가역 행렬의 경우:

$$\det(A) = \pm (\text{피벗들의 곱})$$

**유도:**

$$\det(A) = \det(LU) = \underbrace{\det(L)}_{= 1}\det(U) = \begin{vmatrix} * & * & * \\ & * & * \\ & & * \end{vmatrix} = \text{피벗들의 곱}$$

순열이 있는 경우:

$$\det(A) = \det(PLU) = \underbrace{\det(P)}_{\pm 1} \underbrace{\det(L)}_{1} \det(U) = \pm \text{피벗들의 곱}$$

### 3.10 역행렬의 행렬식

$A$가 가역이면:

$$AA^{-1} = I$$

$$\det(AA^{-1}) = \det(I) = 1$$

$$\det(A) \cdot \det(A^{-1}) = 1$$

$$\boxed{\det(A^{-1}) = \frac{1}{\det(A)}}$$

또한 다음이 성립한다:

$$\det(AB) = \det(A)\det(B) = \det(B)\det(A) = \det(BA)$$

### 3.11 det(A+B)는 det(A) + det(B)가 아니다

$$\det(A + B) \stackrel{?}{=} \det(A) + \det(B) \quad \text{(일반적으로 거짓!)}$$

**반례:**

$$A = \begin{pmatrix} 5 & -6 \\ 0 & -12 \end{pmatrix}, \quad B = \begin{pmatrix} -3 & 0 \\ 1 & 9 \end{pmatrix}$$

$$A + B = \begin{pmatrix} 2 & -6 \\ 1 & -3 \end{pmatrix}$$

$$\det(A) = -60, \quad \det(B) = -27, \quad \det(A+B) = -6 + 6 = 0$$

$$0 \neq -60 + (-27) = -87$$

### 3.12 단일 행에 대한 선형성

행렬식은 다른 행이 고정되어 있을 때, 단일 행에 대해 **선형**이다.

**$2 \times 2$ 경우 — 1행에 대한 선형성:**

$$A = \begin{pmatrix} a_1 + a_2 & b_1 + b_2 \\ c & d \end{pmatrix}$$

$$\det(A) = (a_1 + a_2)d - c(b_1 + b_2) = (a_1 d - cb_1) + (a_2 d - cb_2)$$

$$= \begin{vmatrix} a_1 & b_1 \\ c & d \end{vmatrix} + \begin{vmatrix} a_2 & b_2 \\ c & d \end{vmatrix}$$

**$2 \times 2$ 경우 — 2열에 대한 선형성 (2행 원소 분리):**

$$A = \begin{pmatrix} a & b \\ c_1 + c_2 & d_1 + d_2 \end{pmatrix}$$

$$\det(A) = a(d_1 + d_2) - b(c_1 + c_2) = (ad_1 - bc_1) + (ad_2 - bc_2)$$

$$= \begin{vmatrix} a & b \\ c_1 & d_1 \end{vmatrix} + \begin{vmatrix} a & b \\ c_2 & d_2 \end{vmatrix}$$

**일반 $n \times n$ 경우:** $A$의 행 $i$의 원소가 $a_{ij} = \alpha_{ij} + \beta_{ij}$이면:

$$\det(A) = a_{i1}C_{i1} + \cdots + a_{in}C_{in}$$

$$= (\alpha_{i1} + \beta_{i1})C_{i1} + \cdots + (\alpha_{in} + \beta_{in})C_{in}$$

$$= (\alpha_{i1}C_{i1} + \alpha_{i2}C_{i2} + \cdots + \alpha_{in}C_{in}) + (\beta_{i1}C_{i1} + \beta_{i2}C_{i2} + \cdots + \beta_{in}C_{in})$$

$$= \det(A_\alpha) + \det(A_\beta)$$

여기서 $A_\alpha$는 행 $i$를 $(\alpha_{i1}, \ldots, \alpha_{in})$으로 대체한 행렬, $A_\beta$는 행 $i$를 $(\beta_{i1}, \ldots, \beta_{in})$으로 대체한 행렬이다.

### 3.13 det(A) = det(A^T)의 증명

$$\det(A) = \det(A^T)$$

**증명.** $R_0 = \text{rref}(A)$ ("기약 행사다리꼴", reduced row-echelon form)이라 하자.

$A = E_p E_{p-1} \cdots E_2 E_1 R$이라 표기하자 (여기서 $E_i$는 기본 행렬).

전치: $A^T = R^T E_1^T E_2^T \cdots E_{p-1}^T E_p^T$

$$\det(A^T) = \det(R^T)\det(E_1^T) \cdots \det(E_p^T)$$

**$R$에 대한 고찰:**
- $A$가 가역이면, $R = I$이므로 $\det(R) = 1 = \det(R^T)$
- 그렇지 않으면, $R$은 영행을 가지므로 $\det(R) = 0 = \det(R^T)$
- 따라서 모든 경우에: $\det(R) = \det(R^T)$

임의의 기본 행렬에 대해 $\det(E_i^T) = \det(E_i)$이므로:

$$\det(A^T) = \det(R)\det(E_1) \cdots \det(E_p) = \det(A) \qquad \square$$

### 3.14 크래머 법칙

**$A\mathbf{x} = \mathbf{b}$를 풀기 위한 크래머 법칙 (Cramer's Rule):**

$A = (\mathbf{a}_1 \; \mathbf{a}_2 \; \mathbf{a}_3)$으로, 열벡터가 $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$이라 하자.

**$x_1$ 구하기:** $I$의 1열을 $\mathbf{x}$로 대체:

$$M_1 = \begin{pmatrix} x_1 & 0 & 0 \\ x_2 & 1 & 0 \\ x_3 & 0 & 1 \end{pmatrix}$$

$M_1$에 $A$를 적용:

$$AM_1 = \bigl(A\mathbf{x} \;\; A\begin{pmatrix}0\\1\\0\end{pmatrix} \;\; A\begin{pmatrix}0\\0\\1\end{pmatrix}\bigr) = (\mathbf{b} \;\; \mathbf{a}_2 \;\; \mathbf{a}_3) = B_1$$

행렬식을 취하면:

$$\det(AM_1) = \det(B_1)$$

$$\det(A)\det(M_1) = \det(B_1)$$

$$\det(A) \cdot x_1 = \det(B_1)$$

$$\therefore x_1 = \frac{\det(B_1)}{\det(A)}$$

**$x_2$ 구하기:** $M_2$를 도입:

$$AM_2 = A\begin{pmatrix} 1 & x_1 & 0 \\ 0 & x_2 & 0 \\ 0 & x_3 & 1 \end{pmatrix} = (\mathbf{a}_1 \;\; \mathbf{b} \;\; \mathbf{a}_3) = B_2$$

$$\therefore x_2 = \frac{\det(B_2)}{\det(A)}$$

**$x_3$ 구하기:** $M_3$를 도입:

$$AM_3 = A\begin{pmatrix} 1 & 0 & x_1 \\ 0 & 1 & x_2 \\ 0 & 0 & x_3 \end{pmatrix} = (\mathbf{a}_1 \;\; \mathbf{a}_2 \;\; \mathbf{b}) = B_3$$

$$\therefore x_3 = \frac{\det(B_3)}{\det(A)}$$

**일반 크래머 법칙:** $\det(A) \neq 0$이면, $A\mathbf{x} = \mathbf{b}$의 해는:

$$\boxed{x_j = \frac{\det(B_j)}{\det(A)}}$$

여기서 $B_j$는 $A$의 $j$번째 열을 벡터 $\mathbf{b}$로 대체한 행렬이다.

**예제 1:**

$$3x_1 + 4x_2 = 2, \quad 5x_1 + 6x_2 = 4$$

$$A = \begin{pmatrix} 3 & 4 \\ 5 & 6 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 2 \\ 4 \end{pmatrix}$$

$$\det(A) = 18 - 20 = -2$$

$$B_1 = \begin{pmatrix} 2 & 4 \\ 4 & 6 \end{pmatrix}, \quad \det(B_1) = 12 - 16 = -4$$

$$B_2 = \begin{pmatrix} 3 & 2 \\ 5 & 4 \end{pmatrix}, \quad \det(B_2) = 12 - 10 = 2$$

$$x_1 = \frac{-4}{-2} = 2, \quad x_2 = \frac{2}{-2} = -1$$

**검증:**

$$\begin{pmatrix} 3 & 4 \\ 5 & 6 \end{pmatrix}\begin{pmatrix} 2 \\ -1 \end{pmatrix} = \begin{pmatrix} 6 - 4 \\ 10 - 6 \end{pmatrix} = \begin{pmatrix} 2 \\ 4 \end{pmatrix} \checkmark$$

**예제 2: 크래머 법칙을 이용한 $A^{-1}$ 유도**

$$A = \begin{pmatrix} a & b \\ c & d \end{pmatrix} \in \mathbb{R}^{2 \times 2}, \quad A^{-1} = (\mathbf{x} \;\; \mathbf{y})$$

$$AA^{-1} = A(\mathbf{x} \;\; \mathbf{y}) = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$$

**i)** $A\mathbf{x} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

$$A\begin{pmatrix} x_1 & 0 \\ x_2 & 1 \end{pmatrix} = \bigl(\mathbf{b} \;\; A\begin{pmatrix}0\\1\end{pmatrix}\bigr) = \begin{pmatrix} 1 & b \\ 0 & d \end{pmatrix} = B_1$$

$$x_1 = \frac{\det(B_1)}{\det(A)} = \frac{d}{ad - bc}$$

$$A\begin{pmatrix} 1 & x_1 \\ 0 & x_2 \end{pmatrix} = \bigl(A\begin{pmatrix}1\\0\end{pmatrix} \;\; \mathbf{b}\bigr) = \begin{pmatrix} a & 1 \\ c & 0 \end{pmatrix} = B_2$$

$$x_2 = \frac{\det(B_2)}{\det(A)} = \frac{-c}{ad - bc}$$

이것으로 표준 역행렬 공식을 복원할 수 있다.

---

<br>

## 4. 행렬식에 의한 넓이와 부피 (N5-3)

### 4.1 2D에서의 평행사변형

2D에서 평행사변형은 $(0, 0)$에서 시작하며, 변이 다음과 같다:

$$\mathbf{e}_1 = (a, b), \quad \mathbf{e}_2 = (c, d)$$

네 꼭짓점은 $(0,0)$, $(a,b)$, $(c,d)$, $(a+c, b+d)$이다.

### 4.2 평행사변형의 넓이

$$\text{평행사변형의 넓이} = |ad - bc|$$

$$= \left|\begin{array}{cc} a & b \\ c & d \end{array}\right| = \left|\det\begin{pmatrix} a & b \\ c & d \end{pmatrix}\right| = \left|\det\begin{pmatrix} a & c \\ b & d \end{pmatrix}\right| = |\det(\mathbf{e}_1 \;\; \mathbf{e}_2)|$$

### 4.3 3D에서의 기울어진 상자

3D에서 기울어진 상자는 세 변 $\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3$이 $(0, 0, 0)$에서 출발한다.

$$\text{기울어진 상자의 부피} = \left|\det\begin{pmatrix} \mathbf{e}_1 & \mathbf{e}_2 & \mathbf{e}_3 \end{pmatrix}\right|$$

### 4.4 삼각형의 넓이

$$\text{삼각형의 넓이} = \triangle = \frac{1}{2} b \cdot h$$

여기서 $b = \|\mathbf{a}\|$는 밑변, $h = \|\mathbf{b}\| \sin\theta$는 높이이다.

$$\triangle = \frac{1}{2}\|\mathbf{a}\| \|\mathbf{b}\| \sin\theta$$

**방법 i)** 외적 사용:

$$\triangle = \frac{1}{2}\|\mathbf{a} \times \mathbf{b}\|$$

### 4.5 외적 (Cross Product)

3D에서 **외적** (cross product)은:

$$\mathbf{a} \times \mathbf{b} = \|\mathbf{a}\| \|\mathbf{b}\| \sin\theta \; \hat{n}$$

여기서 $\hat{n}$은 $\mathbf{a}$와 $\mathbf{b}$ 모두에 수직인 단위 법선벡터 (unit normal vector)이다 (즉, $\hat{n} \cdot \mathbf{a} = 0$이고 $\hat{n} \cdot \mathbf{b} = 0$).

**2D 벡터 (3D에 $z = 0$으로 포함)의 경우:**

$$\mathbf{a} = \begin{pmatrix} x_1 \\ y_1 \\ 0 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} x_2 \\ y_2 \\ 0 \end{pmatrix}, \quad \hat{n} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}$$

단위벡터: $\hat{i} = \begin{pmatrix}1\\0\\0\end{pmatrix}$, $\hat{j} = \begin{pmatrix}0\\1\\0\end{pmatrix}$, $\hat{k} = \begin{pmatrix}0\\0\\1\end{pmatrix}$

$$\mathbf{a} \times \mathbf{b} = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ x_1 & y_1 & 0 \\ x_2 & y_2 & 0 \end{vmatrix} = \hat{k}\begin{vmatrix} x_1 & y_1 \\ x_2 & y_2 \end{vmatrix} = \hat{k}(x_1 y_2 - x_2 y_1)$$

$$\|\mathbf{a} \times \mathbf{b}\| = |x_1 y_2 - x_2 y_1|$$

**변 벡터의 경우:**

$\mathbf{e}_1 = \begin{pmatrix} a \\ b \\ 0 \end{pmatrix}$, $\mathbf{e}_2 = \begin{pmatrix} c \\ d \\ 0 \end{pmatrix}$이면:

$$\mathbf{e}_1 \times \mathbf{e}_2 = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ a & b & 0 \\ c & d & 0 \end{vmatrix}$$

$$\boxed{\|\mathbf{e}_1 \times \mathbf{e}_2\| = |ad - bc| = 2\triangle = \text{평행사변형의 넓이}}$$

### 4.6 평행사변형 넓이의 기하학적 증명

꼭짓점이 $(0,0)$, $(a,b)$, $(c,d)$, $(a+c, b+d)$인 평행사변형을 고려하자.

이것을 $(a+c) \times (b+d)$ 크기의 직사각형으로 둘러싼다:

$$\text{직사각형 넓이} = (a+c)(b+d) = ab + ad + bc + cd$$

모서리 영역을 빼면:
- 노란 삼각형 2개: $2 \times \frac{1}{2}ab = ab$
- 초록 직사각형 2개: $2 \times bc = 2bc$
- 분홍 삼각형 2개: $2 \times \frac{1}{2}cd = cd$

$$\text{평행사변형 넓이} = (ab + ad + bc + cd) - ab - 2bc - cd = ad - bc$$

절댓값을 취하면: $\text{넓이} = |ad - bc|$

### 4.7 상자의 부피

**부피 = 밑면 넓이 $\times$ 높이**

$$\|\mathbf{e}_1 \times \mathbf{e}_2\| = \text{밑면 넓이 (평행사변형)}$$

$$h = \|\mathbf{e}_3 \cdot \hat{n}\| = \text{높이 (법선 방향으로의 사영)}$$

$$\text{부피} = h \cdot \|\mathbf{e}_1 \times \mathbf{e}_2\| = \|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\|$$

**$\mathbf{e}_1 \times \mathbf{e}_2$ 계산:**

$\mathbf{e}_1 = \begin{pmatrix} a \\ b \\ c \end{pmatrix}$, $\mathbf{e}_2 = \begin{pmatrix} p \\ q \\ r \end{pmatrix}$, $\mathbf{e}_3 = \begin{pmatrix} x \\ y \\ z \end{pmatrix}$으로 놓으면

$$\mathbf{e}_1 \times \mathbf{e}_2 = \begin{vmatrix} \hat{i} & \hat{j} & \hat{k} \\ a & b & c \\ p & q & r \end{vmatrix} = \hat{i}\begin{vmatrix}b&c\\q&r\end{vmatrix} - \hat{j}\begin{vmatrix}a&c\\p&r\end{vmatrix} + \hat{k}\begin{vmatrix}a&b\\p&q\end{vmatrix} = \begin{pmatrix} br - cq \\ cp - ar \\ aq - bp \end{pmatrix}$$

**스칼라 삼중곱 (triple scalar product):**

$$\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2) = x(br - cq) + y(cp - ar) + z(aq - bp)$$

$$= \begin{vmatrix} a & b & c \\ p & q & r \\ x & y & z \end{vmatrix}$$

따라서:

$$\boxed{\text{부피} = \left|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)\right| = \left|\det\begin{pmatrix} a & b & c \\ p & q & r \\ x & y & z \end{pmatrix}\right|}$$

### 4.8 예제: 단위 정육면체와 축 정렬 상자

**단위 정육면체:** 변이 표준 기저벡터 방향.

$$\text{부피} = \det\begin{pmatrix} 1 & & \\ & 1 & \\ & & 1 \end{pmatrix} = 1$$

**축 정렬 상자**, 변의 길이 $a, b, c$:

$$\text{부피} = \det\begin{pmatrix} a & & \\ & b & \\ & & c \end{pmatrix} = abc$$

---

<br>

## 요약

| 개념 | 핵심 내용 |
|:--------|:---------|
| 행렬식 정의 | 정사각행렬에 대응하는 스칼라 값; $\det(A)$ 또는 $|A|$ |
| $2 \times 2$ 행렬식 | $\det\begin{pmatrix}a&b\\c&d\end{pmatrix} = ad - bc$ |
| $3 \times 3$ 행렬식 | $3! = 6$ 순열로부터 6개의 항; $a(qz-ry) - b(pz-rx) + c(py-qx)$ |
| 행/열 교환 | 부호가 바뀜: $\det \to -\det$ |
| 특이행렬 | $\det(A) = 0$ iff $A$가 가역이 아님 iff 열이 종속 |
| 여인수 $C_{ij}$ | $(-1)^{i+j}\det(\text{행 } i, \text{열 } j \text{ 제거 행렬})$ |
| 여인수 전개 | $\det A = \sum_j A_{ij}C_{ij}$ (임의의 행 $i$를 따라) |
| 여인수를 이용한 역행렬 | $A^{-1} = \frac{1}{\det(A)}C^T$, 여기서 $C^T = \text{adj}(A)$ |
| 스칼라 곱 | $\det(\alpha A) = \alpha^n \det(A)$, $A \in \mathbb{R}^{n \times n}$ |
| 전치 | $\det(A^T) = \det(A)$ |
| 곱의 법칙 | $\det(AB) = \det(A)\det(B)$ |
| 소거 행렬 | $\det(E) = 1$이므로 $\det(EA) = \det(A)$ |
| 직교 행렬 | $\det(Q) = \pm 1$ |
| 삼각/대각 행렬 | 행렬식 = 대각 원소들의 곱 |
| 피벗 공식 | $\det(A) = \pm(\text{피벗들의 곱})$, $PA = LU$ 이용 |
| 역행렬의 행렬식 | $\det(A^{-1}) = 1/\det(A)$ |
| $\det(A+B)$ | 일반적으로 $\det(A) + \det(B)$가 아님 |
| 행의 선형성 | 행렬식은 각 행에 대해 개별적으로 선형 |
| 기본 행렬 | $A$는 소거, 순열, 스케일링 행렬로 분해됨 |
| 크래머 법칙 | $x_j = \det(B_j)/\det(A)$, $B_j$는 열 $j$를 $\mathbf{b}$로 대체 |
| 평행사변형 넓이 (2D) | $|ad - bc| = |\det(\mathbf{e}_1\;\mathbf{e}_2)|$ |
| 삼각형 넓이 | $\frac{1}{2}\|\mathbf{a} \times \mathbf{b}\|$ |
| 외적 (Cross product) | $\mathbf{a} \times \mathbf{b} = \|\mathbf{a}\|\|\mathbf{b}\|\sin\theta\;\hat{n}$; $3 \times 3$ 행렬식으로 계산 |
| 상자 부피 (3D) | $|\mathbf{e}_3 \cdot (\mathbf{e}_1 \times \mathbf{e}_2)| = |\det(\mathbf{e}_1\;\mathbf{e}_2\;\mathbf{e}_3)|$ (스칼라 삼중곱) |

---
