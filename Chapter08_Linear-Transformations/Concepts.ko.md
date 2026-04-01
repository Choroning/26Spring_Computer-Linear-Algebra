# 제8장 강의 — 선형 변환

> **최종 수정일:** 2026-03-31
>
> Strang, Introduction to Linear Algebra 6th Ed. Ch 8

> **선수 지식**: [선형대수학] 벡터 공간, 행렬 연산, 고유값 (제1-7장).
>
> **학습 목표**:
> 1. 선형 변환을 정의하고 선형성을 검증할 수 있다
> 2. 선형 변환의 행렬 표현을 구성할 수 있다
> 3. 변환의 핵(kernel)과 치역(range)을 분석할 수 있다

---

<br>

## 목차

- [1. 선형 변환의 개념 (8.1)](#1-선형-변환의-개념-81)
  - [1.1 정의와 선형성 조건](#11-정의와-선형성-조건)
  - [1.2 선형 변환의 형식적 정의](#12-선형-변환의-형식적-정의)
  - [1.3 비선형 예제: 이동 변환과 아핀 사상](#13-비선형-예제-이동-변환과-아핀-사상)
  - [1.4 선형 및 비선형 변환의 예제들](#14-선형-및-비선형-변환의-예제들)
  - [1.5 기저에 의해 결정되는 선형 변환](#15-기저에-의해-결정되는-선형-변환)
  - [1.6 기하학적 해석: 직선을 직선으로](#16-기하학적-해석-직선을-직선으로)
  - [1.7 미적분에서의 선형 변환 (미분)](#17-미적분에서의-선형-변환-미분)
  - [1.8 예제 5: 적분도 선형이다](#18-예제-5-적분도-선형이다)
  - [1.9 예제 6: z=1 평면으로의 사영 (비선형)](#19-예제-6-z1-평면으로의-사영-비선형)
  - [1.10 예제 7: 가역 행렬과 치역/핵](#110-예제-7-가역-행렬과-치역핵)
- [2. 선형 변환의 행렬 (8.2)](#2-선형-변환의-행렬-82)
  - [2.1 핵심 아이디어 요약](#21-핵심-아이디어-요약)
  - [2.2 기저의 선택과 표준 행렬](#22-기저의-선택과-표준-행렬)
  - [2.3 예제 1: R^2에서 R^3으로의 표준 기저](#23-예제-1-r2에서-r3으로의-표준-기저)
  - [2.4 T에 대한 행렬의 구성](#24-t에-대한-행렬의-구성)
  - [2.5 기저 변환: 행렬 B](#25-기저-변환-행렬-b)
  - [2.6 예제 3: 다항식의 미분 행렬](#26-예제-3-다항식의-미분-행렬)
  - [2.7 예제 4: 적분 행렬 (미분의 유사역행렬)](#27-예제-4-적분-행렬-미분의-유사역행렬)
  - [2.8 행렬곱 AB는 변환 TS에 대응](#28-행렬곱-ab는-변환-ts에-대응)
  - [2.9 예제 5: 회전의 합성](#29-예제-5-회전의-합성)
  - [2.10 예제 6: 역회전](#210-예제-6-역회전)
- [3. 좋은 기저 찾기 (8.3)](#3-좋은-기저-찾기-83)
  - [3.1 핵심 아이디어 요약](#31-핵심-아이디어-요약)
  - [3.2 기저 변환 공식](#32-기저-변환-공식)
  - [3.3 최적 기저 1: 고유벡터 (대각화)](#33-최적-기저-1-고유벡터-대각화)
  - [3.4 최적 기저 2: 특이벡터 (SVD)](#34-최적-기저-2-특이벡터-svd)
  - [3.5 최적 기저 3: 일반화 고유벡터 (Jordan 형식)](#35-최적-기저-3-일반화-고유벡터-jordan-형식)
  - [3.6 Jordan 형식: 구조와 정의](#36-jordan-형식-구조와-정의)
  - [3.7 예제: 2x2 Jordan 형식](#37-예제-2x2-jordan-형식)
  - [3.8 예제: 3x3 Jordan 형식](#38-예제-3x3-jordan-형식)
  - [3.9 함수 공간의 기저](#39-함수-공간의-기저)
  - [3.10 함수 공간의 직교 기저](#310-함수-공간의-직교-기저)
  - [3.11 Gram-Schmidt를 이용한 Legendre 기저 구성](#311-gram-schmidt를-이용한-legendre-기저-구성)
- [요약](#요약)

---

<br>

## 1. 선형 변환의 개념 (8.1)

### 1.1 정의와 선형성 조건

**장 개요:**

- **8.1** 선형 변환의 개념 (Idea of a Linear Transformation)
- **8.2** 선형 변환의 행렬 (The Matrix of a Linear Transformation)
- **8.3** 좋은 기저 찾기 (The Search for a Good Basis)

**선형 변환**(linear transformation) $T$는 벡터 $\mathbf{u}$를 벡터 $T(\mathbf{u})$로 보낸다:

$$T: \mathbf{u} \longrightarrow T(\mathbf{u})$$

**선형성 조건:**

$$T(c\mathbf{u} + d\mathbf{w}) = c\,T(\mathbf{u}) + d\,T(\mathbf{w})$$

**참고:** $T(\mathbf{0}) = \mathbf{0}$.

따라서 $T(\mathbf{u}) = \mathbf{u} + \mathbf{u}_0$은 선형이 **아니다** ($T(\mathbf{0}) = \mathbf{u}_0 \neq \mathbf{0}$이므로).

입력 벡터 $\mathbf{u}$와 출력 $T(\mathbf{u})$는 $\mathbb{R}^n$ 또는 행렬 공간(matrix space) 또는 함수 공간(function space)에 속할 수 있다.

$A$가 $m \times n$이면, $T(\mathbf{u}) = A\mathbf{u}$는 선형이며, 입력 공간 $\mathbb{R}^n$에서 출력 공간 $\mathbb{R}^m$으로의 변환이다.

**미분**(derivative) $T(f) = \dfrac{df}{dx}$는 선형이다.

**적분**(integral) $T^+(f) = \displaystyle\int_0^x f(t)\,dt$은 그 유사역(pseudoinverse)이다.

두 선형 변환의 **곱**(product) $ST$도 여전히 선형이다:

$$(ST)(\mathbf{u}) = S(T(\mathbf{u}))$$

---

### 1.2 선형 변환의 형식적 정의

$A: \mathbf{u} \longrightarrow A\mathbf{u}$

$A$는 $\mathbf{u}$를 다른 벡터 $A\mathbf{u}$로 변환한다.

$f: x \longrightarrow f(x)$와 유사하게, 변환 $T$는 벡터 $\mathbf{u}$를 다른 벡터 $T(\mathbf{u})$로 보낸다:

$$T: \mathbf{u} \longrightarrow T(\mathbf{u})$$

$A$로 시작하면:

$$A: \mathbf{u} \longrightarrow A\mathbf{u}$$
$$A: \mathbf{w} \longrightarrow A\mathbf{w}$$

$\mathbf{y} = \mathbf{u} + \mathbf{w}$에 대해 무슨 일이 일어나는가?

$$A\mathbf{y} = A(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w}$$

행렬 곱셈 $T(\mathbf{u}) = A\mathbf{u}$는 선형 변환을 제공한다.

변환 $T$는 $V$에 속하는 각 입력 벡터 $\mathbf{u}$에 대해 출력 $T(\mathbf{u})$를 할당한다.

**변환이 선형이 되려면 모든 $\mathbf{u}$와 $\mathbf{w} \in V$에 대해 다음 조건을 만족해야 한다:**

**(a)** $T(\mathbf{u} + \mathbf{w}) = T(\mathbf{u}) + T(\mathbf{w})$

**(b)** $T(c\mathbf{u}) = c\,T(\mathbf{u}) \quad \forall\, c \in \mathbb{C}$

$T$가 **가산적**(additive)이고 **동차적**(homogeneous)이라고 말한다.

$$T(\mathbf{0}) = T(0\mathbf{u}) = 0\,T(\mathbf{u}) = \mathbf{0}$$

(a)와 (b)를 결합하면:

$$T(c\mathbf{u} + d\mathbf{w}) = c\,T(\mathbf{u}) + d\,T(\mathbf{w})$$

**일반화:**

$$T(c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_n\mathbf{x}_n) = c_1\,T(\mathbf{x}_1) + c_2\,T(\mathbf{x}_2) + \cdots + c_n\,T(\mathbf{x}_n)$$

---

### 1.3 비선형 예제: 이동 변환과 아핀 사상

**예제.** 벡터 $\mathbf{u} \in V$를 받아 $\mathbf{u}_0 \in V$를 더하는 $T$를 정의하자:

$$T: \mathbf{u} \longrightarrow \mathbf{u} + \mathbf{u}_0$$
$$T: \mathbf{w} \longrightarrow \mathbf{w} + \mathbf{u}_0$$

가산성 검사:

$$T(\mathbf{u} + \mathbf{w}) \stackrel{?}{=} T(\mathbf{u}) + T(\mathbf{w})$$

$$\mathbf{u} + \mathbf{w} + \mathbf{u}_0 \neq \mathbf{u} + \mathbf{u}_0 + \mathbf{w} + \mathbf{u}_0$$

따라서 $T$는 선형이 **아니다**.

**선형 + 이동 변환:**

$$T(\mathbf{u}) = A\mathbf{u} + \mathbf{u}_0$$

이를 **"아핀"**(affine)이라고 부른다.

컴퓨터 그래픽스에서 아핀 사상(affine mapping)이 사용된다:
- 도형 $\mathbf{x}$로 시작
- 회전: $\mathbf{y} = R\mathbf{x}$ (선형)
- 이동: $\mathbf{z} = \mathbf{y} + \begin{pmatrix} 0 \\ u_0 \end{pmatrix}$ (이동)

---

### 1.4 선형 및 비선형 변환의 예제들

**예제 1 (선형 — 내적):**

$$\mathbf{a} = \begin{pmatrix} 1 \\ 3 \\ 4 \end{pmatrix}, \quad \mathbf{u} = \begin{pmatrix} u_1 \\ u_2 \\ u_3 \end{pmatrix}$$

$$T(\mathbf{u}) = \mathbf{a} \cdot \mathbf{u} = \mathbf{a}^T\mathbf{u} = (\mathbf{a}^T)\mathbf{u} = u_1 + 3u_2 + 4u_3$$

**내적(dot product)은 선형이다.**

---

**예제 2 (비선형 — 길이/노름):**

$T(\mathbf{u}) = \|\mathbf{u}\|$는 선형이 **아니다**.

**검사 (i) — 가산성:**

$$T(\mathbf{u} + \mathbf{w}) \stackrel{?}{=} T(\mathbf{u}) + T(\mathbf{w})$$

$$\|\mathbf{u} + \mathbf{w}\| \neq \|\mathbf{u}\| + \|\mathbf{w}\|$$

**검사 (ii) — 동차성:**

$$T(c\mathbf{u}) \stackrel{?}{=} c\,T(\mathbf{u})$$

$$T(-\mathbf{u}) \stackrel{?}{=} -T(\mathbf{u})$$

$$\|-\mathbf{u}\| = \|\mathbf{u}\| \neq -\|\mathbf{u}\|$$

동차성에 실패한다 (노름은 항상 비음수).

---

**예제 3 (선형 — 회전):**

$T$는 모든 벡터를 $30°$만큼 회전시키는 회전 행렬이다.

$$T: \mathbf{u} \longrightarrow \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\mathbf{u} \quad \text{( } \theta = 30° \text{)}$$

여기서 회전 행렬은 $R$이다.

**검사 (i):**

$$T(\mathbf{u} + \mathbf{w}) = R(\mathbf{u} + \mathbf{w}) = R\mathbf{u} + R\mathbf{w} = T(\mathbf{u}) + T(\mathbf{w})$$

**검사 (ii):**

$$T(c\mathbf{u}) = R(c\mathbf{u}) = cR\mathbf{u} = c\,T(\mathbf{u})$$

**일반화를 상기하면:**

$$T(c_1\mathbf{u}_1 + c_2\mathbf{u}_2 + \cdots + c_n\mathbf{u}_n) = c_1\,T(\mathbf{u}_1) + c_2\,T(\mathbf{u}_2) + \cdots + c_n\,T(\mathbf{u}_n)$$

---

### 1.5 기저에 의해 결정되는 선형 변환

기저(basis)에 속하는 모든 벡터 $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n$에 대해 $T(\mathbf{u})$를 안다고 하자.

$\mathbf{y} = c_1\mathbf{u}_1 + \cdots + c_n\mathbf{u}_n$으로 놓으면,

공간의 모든 $\mathbf{y}$에 대해 $T(\mathbf{y})$를 알 수 있다.

> 선형 변환은 기저에 대한 작용에 의해 완전히 결정된다.

---

### 1.6 기하학적 해석: 직선을 직선으로

선형 변환은 다음을 보존한다:
- **직선을 직선으로** (lines to lines)
- **삼각형을 삼각형으로** (triangles to triangles)
- **등간격 점들은 등간격 점들로** (equally spaced points go to equally spaced points)

---

### 1.7 미적분에서의 선형 변환 (미분)

**예제 4.** $T(u) = \dfrac{d}{dt}(u)$

**검사 (i):**

$$T(cu + dw) = \frac{d}{dt}(cu + dw) = c\frac{du}{dt} + d\frac{dw}{dt} = c\,T(u) + d\,T(w) \quad \text{선형}$$

여기서 $c, d \in \mathbb{R}$.

**$T$의 영공간 (nullspace):**

$$T(u) = \frac{d}{dt}u = 0$$

이는 $u$가 상수 $c$일 때만 성립한다.

$$dc \in \mathcal{N}(T), \quad d \in \mathbb{R}$$

영공간은 함수 공간에서의 **직선**이다.

**$T$의 열공간 (column space):**

$u = a + bt + ct^2$이면:

$$T(u) = b + 2ct, \quad \text{일차 함수이다.}$$

**차원:**

$$\dim C(T) + \dim \mathcal{N}(T) = 2 + 1 = 3$$

**미분 $T = \dfrac{d}{dt}$에 대한 행렬:**

$u = (1 \;\; t \;\; t^2) \begin{pmatrix} a \\ b \\ c \end{pmatrix}$로 보자.

$\mathbf{u}_1 = 1$, $\mathbf{u}_2 = t$, $\mathbf{u}_3 = t^2$로 놓으면:

$$\frac{d}{dt}\mathbf{u}_1 = \frac{d}{dt}(1) = 0$$

$$\frac{d}{dt}\mathbf{u}_2 = \frac{d}{dt}(t) = 1$$

$$\frac{d}{dt}\mathbf{u}_3 = \frac{d}{dt}(t^2) = 2t$$

$$T: \begin{pmatrix} a \\ b \\ c \end{pmatrix} \longrightarrow \begin{pmatrix} b \\ 2c \end{pmatrix}$$

$$a + bx + cx^2 \longrightarrow b + 2cx$$

선형 변환 $\dfrac{dy}{dx}$는 $A\mathbf{u}$에 대응된다:

$$\begin{pmatrix} b \\ 2c \end{pmatrix} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} a \\ b \\ c \end{pmatrix}$$

행렬은 $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}$이다.

---

### 1.8 예제 5: 적분도 선형이다

적분 $T^+$도 선형이다:

$$\int_0^x (D + Ex)\,dx = Dx + \frac{1}{2}Ex^2$$

$$T^+: \begin{pmatrix} D \\ E \end{pmatrix} \longrightarrow \begin{pmatrix} 0 \\ D \\ \frac{1}{2}E \end{pmatrix}$$

입력 $\mathbf{u} = D + Ex$이고 출력은 $Dx + \frac{1}{2}Ex^2$이다.

$$\begin{pmatrix} 0 \\ D \\ \frac{1}{2}E \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix}\begin{pmatrix} D \\ E \end{pmatrix}$$

적분 행렬은 $A^+$이다.

**곱 $A^+A$와 $AA^+$:**

$$A^+A = \begin{pmatrix} 0 & 0 \\ 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix}\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$$AA^+ = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} 0 & 0 \\ 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$$

---

### 1.9 예제 6: z=1 평면으로의 사영 (비선형)

벡터 $\mathbf{u} \in \mathbb{R}^3$를 수평면 $z = 1$로 사영(projection)하자:

$$T: \begin{pmatrix} x \\ y \\ z \end{pmatrix} \longrightarrow \begin{pmatrix} x \\ y \\ 1 \end{pmatrix}$$

$T(\mathbf{u})$는 선형이 **아니다**:

$$T(\mathbf{0}) = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \longrightarrow \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix} \neq \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

---

### 1.10 예제 7: 가역 행렬과 치역/핵

$A$가 가역(invertible)이라고 하자:

$$T(\mathbf{u} + \mathbf{w}) = A\mathbf{u} + A\mathbf{w} = T(\mathbf{u}) + T(\mathbf{w})$$

역변환: $T^{-1}(\mathbf{u}) = A^{-1}(\mathbf{u})$

$$T^{-1}(T(\mathbf{u})) = A^{-1}(A\mathbf{u}) = \mathbf{u}$$

$T(\mathbf{u}) = A\mathbf{u}$, $S(\mathbf{u}) = B\mathbf{u}$이면, $T \circ S(\mathbf{u})$는 $AB\mathbf{u}$에 대응된다.

**질문:** $V = \mathbb{R}^n$에서 $W = \mathbb{R}^m$으로의 모든 선형 변환은 행렬에 의해 생성되는가? **그렇다.**

$A\mathbf{u}$는 **열공간**(column space)이다.

$A$의 **영공간**(null space)은 $A\mathbf{u} = \mathbf{0}$이 되는 모든 입력을 포함한다.

$$\text{$T$의 치역 (Range)} = T(\mathbf{u}) \quad \longleftrightarrow \quad \text{열공간 (column space) } A\mathbf{u}$$

$$\text{$T$의 핵 (Kernel)} = T(\mathbf{u}) = \mathbf{0}\text{인 모든 입력} \quad \longleftrightarrow \quad \text{$A$의 영공간 (nullspace)}$$

---

<br>

## 2. 선형 변환의 행렬 (8.2)

### 2.1 핵심 아이디어 요약

1. **선형성**은 입력 기저 $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n$에 대해 $T(\mathbf{u}_1), \ldots, T(\mathbf{u}_n)$을 알면 모든 $T(\mathbf{u})$를 알 수 있게 한다.

2. $T$의 행렬에서 **열 $j$**는 입력 기저 벡터 $\mathbf{u}_j$에 $T$를 적용하여 얻는다.

3. $T(\mathbf{u}_j) = a_{1j}\,\mathbf{w}_1 + \cdots + a_{mj}\,\mathbf{w}_m = \displaystyle\sum_{i=1}^{m} a_{ij}\,\mathbf{w}_i$를 출력 기저 $\mathbf{w}$로 표현한다. 그 $a_{ij}$가 열 $j$에 들어간다.

4. 입력과 출력 기저가 $I_{n \times n}$과 $I_{m \times m}$의 열이면, 행렬 $T(\mathbf{x}) = A\mathbf{x}$는 $A$이다.

5. 기저가 $\mathbf{v}$와 $\mathbf{w}$로 바뀌면, 같은 $T$에 대한 행렬이 $A$에서 $W^{-1}AV$로 변한다.

6. **최적 기저:** $V = W$ = 고유벡터(eigenvectors)와 $V, W$ = 특이벡터(singular vectors): $A \to \Lambda$ 그리고 $\Sigma$.

---

### 2.2 기저의 선택과 표준 행렬

변환 $T$는 입력 공간 $V$에서 출력 공간 $W$로 사상한다:

$$V \xrightarrow{T} W$$

$V = \mathbb{R}^n$, $W = \mathbb{R}^m$을 선택하면, 이 변환의 행렬 $A$는 $m \times n$이 된다.

$V$와 $W$에서의 기저 선택이 $A$를 결정한다.

$\mathbb{R}^n$과 $\mathbb{R}^m$의 표준 기저 벡터는 $I$의 열이다.

이 선택은 **표준 행렬**(standard matrix)로 이어진다: $T(\mathbf{u}) = A\mathbf{u}$.

$V$와 $W$에 대해 다른 기저를 선택할 수 있다.
$\Rightarrow$ 같은 변환 $T$가 다른 행렬로 표현된다.

**질문.** $T$에 대해 **최적의 행렬**을 주는 기저를 어떻게 선택해야 하는가?

---

### 2.3 예제 1: R^2에서 R^3으로의 표준 기저

$T: \mathbb{R}^2 \to \mathbb{R}^3$

$$\mathbf{u}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix} \longrightarrow T(\mathbf{u}_1) = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix}$$

$$\mathbf{u}_2 = \begin{pmatrix} 0 \\ 1 \end{pmatrix} \longrightarrow T(\mathbf{u}_2) = \begin{pmatrix} 5 \\ 5 \\ 5 \end{pmatrix}$$

이 선형 변환을 행렬 $A$로 표현하면:

$$A = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) \end{pmatrix} = \begin{pmatrix} 2 & 5 \\ 3 & 5 \\ 4 & 5 \end{pmatrix}$$

출력이 **$A$의 열**로 들어간다. 그러면 $T(\mathbf{u}) = A\mathbf{u}$.

**검증:**

$$T(\mathbf{u}_1 + \mathbf{u}_2) = A(\mathbf{u}_1 + \mathbf{u}_2) = A\mathbf{u}_1 + A\mathbf{u}_2 = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix} + \begin{pmatrix} 5 \\ 5 \\ 5 \end{pmatrix} = A\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 7 \\ 8 \\ 9 \end{pmatrix}$$

**일반 공식:**

$$T(c_1\mathbf{u}_1 + c_2\mathbf{u}_2) = c_1\,T(\mathbf{u}_1) + c_2\,T(\mathbf{u}_2) = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = A\mathbf{c}$$

---

### 2.4 T에 대한 행렬의 구성

**임의의** 선형 변환에 대한 행렬을 구성한다.

$$V \xrightarrow{T} W, \quad \dim V = n, \quad \dim W = m$$

$V$의 기저 $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_n$을 선택한다.

$W$의 기저 $\mathbf{w}_1, \mathbf{w}_2, \ldots, \mathbf{w}_m$을 선택한다.

$A$는 $m \times n$이 된다.

$T(\mathbf{u}_1)$은 $W$의 출력 기저의 결합이다:

$$T(\mathbf{u}_1) = d_1\mathbf{w}_1 + d_2\mathbf{w}_2 + \cdots + d_m\mathbf{w}_m = a_{11}\mathbf{w}_1 + a_{21}\mathbf{w}_2 + \cdots + a_{m1}\mathbf{w}_m = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{11} \\ a_{21} \\ \vdots \\ a_{m1} \end{pmatrix}$$

$$T(\mathbf{u}_2) = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{12} \\ a_{22} \\ \vdots \\ a_{m2} \end{pmatrix}$$

$$\vdots$$

$$T(\mathbf{u}_n) = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{1n} \\ a_{2n} \\ \vdots \\ a_{mn} \end{pmatrix}$$

이를 결합하면:

$$\begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) & \cdots & T(\mathbf{u}_n) \end{pmatrix} = A$$

$$= (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \cdots \;\; \mathbf{w}_m)\begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}$$

> **$A$의 $j$번째 열은 $j$번째 기저 벡터 $\mathbf{u}_j$에 $T$를 적용하여 얻으며, 이는 출력 기저 벡터의 선형 결합이다.**

입력 기저 벡터 $\mathbf{u}_1$부터 $\mathbf{u}_n$에 대해 출력 $T(\mathbf{u})$를 안다고 하면:

$$T(\mathbf{u}_1),\; T(\mathbf{u}_2),\; \ldots,\; T(\mathbf{u}_n)$$

$$A = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) & \cdots & T(\mathbf{u}_n) \end{pmatrix}$$

$$A\mathbf{c} = c_1\,T(\mathbf{u}_1) + c_2\,T(\mathbf{u}_2) + \cdots + c_n\,T(\mathbf{u}_n)$$

---

### 2.5 기저 변환: 행렬 B

**예제 2.** $V = W = \mathbb{R}^2$

$T(\mathbf{u}) = \mathbf{u}$ (항등 변환, identity transformation)이라 하자.

**입력 기저** ($V$): 표준 기저 $\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$.

**출력 기저** ($W$): 회전된 기저 $\mathbf{w}_1 = \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}$, $\mathbf{w}_2 = \begin{pmatrix} -\sin\theta \\ \cos\theta \end{pmatrix}$.

$\mathbf{u}$를 두 기저로 표현하면:

$$\mathbf{u} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 = (\mathbf{v}_1 \;\; \mathbf{v}_2)\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = V\mathbf{c}$$

$$\mathbf{u} = d_1\mathbf{w}_1 + d_2\mathbf{w}_2 = (\mathbf{w}_1 \;\; \mathbf{w}_2)\begin{pmatrix} d_1 \\ d_2 \end{pmatrix} = W\mathbf{d}$$

둘 다 같은 벡터를 나타내므로:

$$\boxed{V\mathbf{c} = W\mathbf{d}}$$

$$\mathbf{c} = V^{-1}W\mathbf{d} \qquad \mathbf{d} = W^{-1}V\mathbf{c}$$

$$B = W^{-1}V$$

**이 예제에서:**

$$\mathbf{d} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}^{-1}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\mathbf{c}$$

$$= \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\mathbf{c}$$

$$= \begin{pmatrix} \cos(-\theta) & -\sin(-\theta) \\ \sin(-\theta) & \cos(-\theta) \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\mathbf{c}$$

$V = I$일 때, $W$에서의 계수 $\mathbf{d}$는:

$$\mathbf{d} = W^{-1}\begin{pmatrix} x \\ y \end{pmatrix}$$

---

### 2.6 예제 3: 다항식의 미분 행렬

$T(u) = \dfrac{du}{dx}$

입력 공간: $V = \text{span}\{1, x, x^2, x^3\}$, 기저 $\mathbf{u}_1 = 1,\; \mathbf{u}_2 = x,\; \mathbf{u}_3 = x^2,\; \mathbf{u}_4 = x^3$.

출력 공간: $W = \text{span}\{1, x, x^2\}$, 기저 $\mathbf{w}_1 = 1,\; \mathbf{w}_2 = x,\; \mathbf{w}_3 = x^2$.

각 기저 벡터에 $T$를 적용하면:

$$T(\mathbf{u}_1) = T(1) = 0$$
$$T(\mathbf{u}_2) = T(x) = 1 = \mathbf{w}_1$$
$$T(\mathbf{u}_3) = T(x^2) = 2x = 2\mathbf{w}_2$$
$$T(\mathbf{u}_4) = T(x^3) = 3x^2 = 3\mathbf{w}_3$$

$$A = \begin{pmatrix} T(\mathbf{u}_1) & T(\mathbf{u}_2) & T(\mathbf{u}_3) & T(\mathbf{u}_4) \end{pmatrix}$$

$$= (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \mathbf{w}_3)\begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$$

**미분 행렬**(derivative matrix)은:

$$D = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$$

입력 $u = c_1 + c_2 x + c_3 x^2 + c_4 x^3$이고 계수 벡터 $\mathbf{c} = (c_1, c_2, c_3, c_4)^T$이면:

$$T(u) = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \mathbf{w}_3) \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \\ c_3 \\ c_4 \end{pmatrix} = (\mathbf{w}_1 \;\; \mathbf{w}_2 \;\; \mathbf{w}_3)\begin{pmatrix} c_2 \\ 2c_3 \\ 3c_4 \end{pmatrix}$$

$\mathbf{d} = D\mathbf{c}$는 $T(u)$ 결합의 계수를 생성한다.

---

### 2.7 예제 4: 적분 행렬 (미분의 유사역행렬)

$d_1 + d_2 x + d_3 x^2$의 적분:

$$(1 \;\; x \;\; x^2)\begin{pmatrix} d_1 \\ d_2 \\ d_3 \end{pmatrix}$$

은 $d_1 x + \frac{1}{2}d_2 x^2 + \frac{1}{3}d_3 x^3$이다:

$$(1 \;\; x \;\; x^2 \;\; x^3)\begin{pmatrix} 0 \\ d_1 \\ \frac{1}{2}d_2 \\ \frac{1}{3}d_3 \end{pmatrix}$$

$$= (1 \;\; x \;\; x^2 \;\; x^3)\begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix}\begin{pmatrix} d_1 \\ d_2 \\ d_3 \end{pmatrix}$$

**적분 행렬**(integral matrix)은:

$$D^+ = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix}$$

**곱 $DD^+$와 $D^+D$:**

$$DD^+ = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}\begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = I$$

$$D^+D = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 0 & \frac{1}{2} & 0 \\ 0 & 0 & \frac{1}{3} \end{pmatrix}\begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix} = \begin{pmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

미분 $T$는 **핵**(kernel)을 갖는다 (상수 함수): $T(1) = 0$.

그 행렬 $D$는 **영공간**(nullspace)을 갖는다.

---

### 2.8 행렬곱 AB는 변환 TS에 대응

$TS$를 생각하자:

$$(TS)(\mathbf{u}) := T \circ S(\mathbf{u}) = T(S(\mathbf{u}))$$

$$(AB)(\mathbf{x}) := A(B\mathbf{x})$$

행렬 곱셈은 $TS$를 표현하는 올바른 행렬 $AB$를 제공한다.

$$U \xrightarrow{S} V \xrightarrow{T} W$$

$$\dim U = p, \quad \dim V = n, \quad \dim W = m$$

$$B: n \times p, \quad A: m \times n$$

선형 변환 $TS$:
- $U$의 임의의 벡터 $\mathbf{u}$로 시작
- $V$의 $S(\mathbf{u})$로 이동
- 그런 다음 $W$의 $T(S(\mathbf{u}))$로 이동

행렬 $AB$:
- $\mathbb{R}^p$의 임의의 $\mathbf{x}$로 시작
- $\mathbb{R}^n$의 $B\mathbf{x}$로 이동
- 그런 다음 $\mathbb{R}^m$의 $AB\mathbf{x}$로 이동

행렬 $AB$는 $TS$를 올바르게 표현한다:

$$TS: U \to V \to W$$
$$AB: (m \times n)(n \times p) = (m \times p)$$

---

### 2.9 예제 5: 회전의 합성

$S$는 평면을 $\theta$만큼 회전시킨다.

$T$도 $\theta$만큼 회전시킨다.

그러면 $TS$는 $2\theta$만큼 회전시킨다.

이 변환 $T^2$는 $2\theta$를 통한 회전 행렬 $A^2$에 대응된다:

$$T = S, \quad A = B, \quad T^2 = \text{$2\theta$만큼 회전}$$

$$A^2 = \begin{pmatrix} \cos 2\theta & -\sin 2\theta \\ \sin 2\theta & \cos 2\theta \end{pmatrix}$$

**행렬 곱셈으로 검증:**

$$A^2 = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$$

$$= \begin{pmatrix} \cos^2\theta - \sin^2\theta & -2\sin\theta\cos\theta \\ 2\sin\theta\cos\theta & \cos^2\theta - \sin^2\theta \end{pmatrix}$$

이배각 공식(double-angle formulas)을 사용하면:
- $\cos^2\theta - \sin^2\theta = \cos 2\theta$
- $2\sin\theta\cos\theta = \sin 2\theta$

---

### 2.10 예제 6: 역회전

$S$는 각도 $\theta$만큼 회전시킨다.

$T$는 $-\theta$만큼 회전시킨다.

$$TS = T \circ S = I \longrightarrow AB = I$$

$$AB\mathbf{x} = I\mathbf{x} = \mathbf{x}$$

$$AB = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix}\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} = \begin{pmatrix} \cos^2\theta + \sin^2\theta & 0 \\ 0 & \cos^2\theta + \sin^2\theta \end{pmatrix} = I$$

---

<br>

## 3. 좋은 기저 찾기 (8.3)

### 3.1 핵심 아이디어 요약

1. 새로운 입력 기저 $B_{\text{in}}$과 출력 기저 $B_{\text{out}}$에 대해, 모든 행렬 $A$는 $B_{\text{out}}^{-1} A B_{\text{in}}$이 된다.

2. $B_{\text{in}} = B_{\text{out}}$ = $A$의 일반화 고유벡터(generalized eigenvectors)를 사용하면 **Jordan 형식** $J = B^{-1}AB$를 얻는다.

3. **Fourier 행렬** $F = B_{\text{in}} = B_{\text{out}}$은 모든 **순환 행렬**(circulant matrix)을 대각화한다.

4. 사인, 코사인, $e^{ikx}$, **Legendre** 및 **Chebyshev**: 이들은 함수 공간에 대한 훌륭한 기저이다.

---

### 3.2 기저 변환 공식

입력 기저 벡터는 $B_{\text{in}} = (\mathbf{b}_1 \;\; \mathbf{b}_2 \;\; \cdots \;\; \mathbf{b}_n)$을 형성한다.

출력 기저 벡터는 $B_{\text{out}} = (\mathbf{b}_1' \;\; \mathbf{b}_2' \;\; \cdots \;\; \mathbf{b}_m')$을 형성한다.

항상 $B_{\text{in}}$과 $B_{\text{out}}$은 **가역**(invertible)이다.

새로운 행렬 표현은:

$$B_{\text{out}}^{-1} \; A \; B_{\text{in}}$$

여기서 $B_{\text{out}}$은 $m \times m$, $A$는 $m \times n$, $B_{\text{in}}$은 $n \times n$이다.

**예시:** $B_{\text{in}} = I_{n \times n}$이고 $B_{\text{out}} = I_{m \times m}$이면, 행렬은 $A$로 유지된다.

**$B = B_{\text{in}} = B_{\text{out}}$일 때:**

$$B^{-1}AB \quad \text{는 } A \text{와 닮음(similar)이다}$$

닮음 행렬은 **같은 고유값**(eigenvalues)을 갖는다.

---

### 3.3 최적 기저 1: 고유벡터 (대각화)

$B_{\text{in}} = B_{\text{out}}$ = 고유벡터 행렬(eigenvector matrix) $X$.

$$A = X\Lambda X^{-1} \quad \Longleftrightarrow \quad X^{-1}AX = \Lambda$$

$A$는 $n$개의 독립 고유벡터를 갖는 정사각 행렬이다. $A$는 **대각화 가능**(diagonalizable)해야 한다.

$\Lambda$는 $B_{\text{in}} = B_{\text{out}} = X$ (고유벡터)일 때 얻어진다.

---

### 3.4 최적 기저 2: 특이벡터 (SVD)

$B_{\text{in}} = V$, $B_{\text{out}} = U$: $A$의 특이벡터(singular vectors).

$$A = U\Sigma V^T$$

$$U^{-1}AV = \Sigma \quad \text{(특이값, singular values)}$$

$U, V$는 $A^TA$와 $AA^T$의 **정규직교 고유벡터**(orthonormal eigenvectors)이다.

---

### 3.5 최적 기저 3: 일반화 고유벡터 (Jordan 형식)

$B_{\text{in}} = B_{\text{out}}$ = $A$의 일반화 고유벡터(generalized eigenvectors).

$$B^{-1}AB = \text{Jordan 형식 } J$$

$A$는 정사각 행렬이지만, "$s$개의 독립 고유벡터"만 가질 수 있다:

- $s = n$이면: $B = X$ (고유벡터), $J = \Lambda$
- $s < n$이면: "$n - s$"개의 종속 열이 있다 $\to$ "$n - s$"개의 추가 **일반화 고유벡터**를 구성한다

Jordan 형식의 성질:
- (i) $J$의 대각선을 따라 $s$개의 정사각 블록이 있다
- (ii) 각 블록은 하나의 고유값 $\lambda$, 하나의 고유벡터, 그리고 대각선 위에 1을 갖는다

---

### 3.6 Jordan 형식: 구조와 정의

모든 고유벡터가 독립일 때 ($s = n$):

$$J = \Lambda = \begin{pmatrix} \lambda_1 & & \\ & \lambda_2 & \\ & & \ddots & \\ & & & \lambda_n \end{pmatrix} \quad \text{($n$개의 1x1 블록)}$$

**예제 1:**

$$J = \begin{pmatrix} \boxed{2} & & & \\ & \boxed{2} & & \\ & & \boxed{3} & 1 \\ & & & \boxed{3} \end{pmatrix}$$

이것은 **2**개의 1x1 블록과 **1**개의 2x2 블록을 갖는다.

$B^{-1}AB = J$는 **거의 대각**(nearly diagonal)이다.

**Jordan 형식 — 일반적 진술:**

모든 $A$에 대해, $B^{-1}AB$가 가능한 한 거의 대각이 되도록 $B$를 선택하고자 한다.

- $A$가 대각화 가능할 때 (즉, $n$개의 독립 고유벡터), $B = X$ (고유벡터): $X^{-1}AX = \Lambda$.
- $A$가 대각화 불가능할 때 (즉, $s < n$개의 독립 고유벡터):

$$B^{-1}AB = J$$

$A$의 Jordan 형식은:

$$J = \begin{pmatrix} J_1 & & \\ & J_2 & \\ & & \ddots & \\ & & & J_s \end{pmatrix}$$

각 **Jordan 블록** $J_i$는:

$$J_i = \begin{pmatrix} \lambda_i & 1 & 0 & \cdots & 0 \\ 0 & \lambda_i & 1 & \cdots & 0 \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & \cdots & 0 & \lambda_i & 1 \\ 0 & \cdots & 0 & 0 & \lambda_i \end{pmatrix}$$

> **최적 기저 $B$는 $B^{-1}AB = J$를 준다.**

---

### 3.7 예제: 2x2 Jordan 형식

$$A = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, \quad \lambda^2 - 2\lambda + 1 = (\lambda - 1)^2 = 0$$

**(i)** $\lambda_1 = 1$

$$(A - \lambda I)\mathbf{x} = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$x_2 = 0$, 자유 변수 $x_1 \Rightarrow \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

**(ii)** 일반화 고유벡터를 구할 수 있는가?

$$(A - \lambda_1 I)\mathbf{x}_2 = \mathbf{x}_1$$

(참고: $(A - \lambda_1 I)^2\mathbf{x}_2 = (A - \lambda_1 I)\mathbf{x}_1 = \mathbf{0}$)

$$(A - \lambda I)\mathbf{x}_2 = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$$

$$2x_2 = 1 \quad \Rightarrow \quad x_2 = \frac{1}{2}$$

자유 변수 $x_1 \to x_1 = 0$

$$\therefore \mathbf{x}_2 = \begin{pmatrix} 0 \\ \frac{1}{2} \end{pmatrix} \perp \mathbf{x}_1$$

**(iii)** 기저 행렬을 형성하면:

$$B = (\mathbf{x}_1 \;\; \mathbf{x}_2) = \begin{pmatrix} 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix}, \quad B^{-1} = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$$

$$B^{-1}AB = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & \frac{1}{2} \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & \frac{1}{2} \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = J$$

이것은 **1개의 Jordan 블록** (2x2)이다.

---

### 3.8 예제: 3x3 Jordan 형식

$$A = \begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}, \quad \lambda^3 = 0 \;\Rightarrow\; \lambda = 0$$

**(i)** 고유벡터 찾기: $(A - \lambda I)\mathbf{x} = \mathbf{0}$

$$(A - 0 \cdot I)\mathbf{x} = \begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}$$

$x_3 = 0$, $x_2 + 2x_3 = 0 \Rightarrow x_2 = 0$, 자유 변수 $x_1$

$$\Rightarrow \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$$

**(ii)** 첫 번째 일반화 고유벡터 찾기: $(A - \lambda I)\mathbf{x}_2 = \mathbf{x}_1$

$$\begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$$

$x_3 = 0$, $x_2 + 2x_3 = 1 \Rightarrow x_2 = 1$, $x_1$은 자유 $\to x_1 = 0$

$$\Rightarrow \mathbf{x}_2 = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$$

**(iii)** 두 번째 일반화 고유벡터 찾기: $(A - \lambda I)\mathbf{x}_3 = \mathbf{x}_2$

$$\begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 1 \\ 0 \end{pmatrix}$$

$x_3 = 1$, $x_2 + 2x_3 = 0 \Rightarrow x_2 = -2$, $x_1 = 0$ 선택

$$\Rightarrow \mathbf{x}_3 = \begin{pmatrix} 0 \\ -2 \\ 1 \end{pmatrix}$$

**(iv)** Jordan 연쇄 (Jordan chains):

$$(A - \lambda I)\mathbf{x}_1 = \mathbf{0}$$
$$(A - \lambda I)^2\mathbf{x}_2 = (A - \lambda I)\mathbf{x}_1 = \mathbf{0}$$
$$(A - \lambda I)^3\mathbf{x}_3 = (A - \lambda I)^2\mathbf{x}_2 = \mathbf{0}$$

기저를 형성하면:

$$B = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -2 \\ 0 & 0 & 1 \end{pmatrix}, \quad B^{-1} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}$$

$$B^{-1}AB = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 0 & 1 & 2 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & -2 \\ 0 & 0 & 1 \end{pmatrix}$$

$$= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix} = J$$

이것은 고유값 $\lambda = 0$인 **1개의 Jordan 블록** (3x3)이다.

---

### 3.9 함수 공간의 기저

표준 다항식 기저를 생각하자: $1, x, x^2, x^3, \ldots$

$$x^{10} \approx \text{span}\{1, x, x^2, \ldots, x^9\}$$

이것은 **나쁜 조건수의 기저**(ill-conditioned basis)이다.

기저가 좋은지 아닌지 어떻게 확인할 수 있는가?

$$B = \begin{pmatrix} \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_n \end{pmatrix} \quad \text{(기저 벡터)}$$

$$B^TB = I \quad \text{기저가 정규직교(orthonormal)일 때. 이것이 최적이다.}$$

**두 벡터의 내적:** $\mathbf{b}_i^T\mathbf{b}_j$

**두 함수의 내적:** $\displaystyle\int_0^1 x^i \cdot x^j \, dx$

$$\int_0^1 x^i \cdot x^j \, dx = \frac{1}{i + j + 1} \longrightarrow 0 \quad \text{$i, j \to \infty$일 때}$$

**$f(x)$와 $g(x)$의 내적:**

$$(f, g) := \int f(x)\,g(x)\,dx$$

**복소 내적:**

$$(f, g) := \int \overline{f(x)}\,g(x)\,dx$$

**가중 내적:**

$$(f, g)_w := \int w(x)\,\overline{f(x)}\,g(x)\,dx$$

---

### 3.10 함수 공간의 직교 기저

**Fourier 기저:**

$$1, \;\sin x, \;\cos x, \;\sin 2x, \;\ldots$$

**Legendre 기저:**

$$1, \;x, \;x^2 - \frac{1}{3}, \;x^3 - \frac{3}{5}x, \;\ldots$$

**Chebyshev 기저:**

$$1, \;x, \;2x^2 - 1, \;4x^3 - 3x, \;\ldots$$

$x$의 거듭제곱 $1, x, x^2, \ldots$으로부터 어떻게 Legendre 기저를 구성할 수 있는가?

---

### 3.11 Gram-Schmidt를 이용한 Legendre 기저 구성

$[-1, 1]$에서의 내적을 계산하면:

$$(1, 1) = \int_{-1}^{1} 1\,dx = 2$$

$$(1, x) = \int_{-1}^{1} 1 \cdot x\,dx = 0 \quad \text{(직교)}$$

$$(1, x^2) = \int_{-1}^{1} 1 \cdot x^2\,dx = 2\int_0^1 x^2\,dx = \frac{2}{3}$$

$x^2$은 $1$의 영이 아닌 성분을 갖는다.

$$(x, x^2) = \int_{-1}^{1} x^3\,dx = 0 \quad \text{(직교)}$$

**Gram-Schmidt 과정**으로 직교 기저를 구성하면:

$$b_1 = 1$$
$$b_2 = x \quad \perp\; b_1$$
$$b_3 = x^2 - \frac{1}{3} \quad \perp\; b_1,\; b_2$$
$$\vdots$$

**Gram-Schmidt 과정**을 사용하여 **Legendre 기저**를 얻을 수 있다.

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:-----|:-----------|
| 선형 변환 (Linear Transformation) | $T(c\mathbf{u} + d\mathbf{w}) = cT(\mathbf{u}) + dT(\mathbf{w})$; 가산성과 동차성을 만족해야 함 |
| $T(\mathbf{0}) = \mathbf{0}$ | 모든 선형 변환은 영벡터를 영벡터로 보냄 |
| 비선형 예제 | 노름 $\|\mathbf{u}\|$, 이동 $\mathbf{u} + \mathbf{u}_0$은 선형이 아님 |
| 아핀 변환 (Affine) | $T(\mathbf{u}) = A\mathbf{u} + \mathbf{u}_0$ (선형 + 이동); 컴퓨터 그래픽스에서 사용 |
| 내적 (Dot product) | $T(\mathbf{u}) = \mathbf{a}^T\mathbf{u}$는 선형 |
| 미분 (Derivative) | $T(u) = du/dt$는 함수 공간에서의 선형 변환 |
| 적분 (Integration) | $T^+(f) = \int_0^x f(t)\,dt$도 선형; 행렬 $A^+$은 $AA^+ = I$를 만족 |
| $z=1$으로의 사영 | $T(\mathbf{0}) \neq \mathbf{0}$이므로 선형이 아님 |
| 치역과 핵 (Range and Kernel) | $T$의 치역 = $A$의 열공간; $T$의 핵 = $A$의 영공간 |
| 기저에 의한 결정 | 선형 변환은 기저 벡터에 대한 값에 의해 완전히 결정됨 |
| $T$의 행렬 | 열 $j$ = $T(\mathbf{u}_j)$를 출력 기저로 표현; $A = (T(\mathbf{u}_1) \;\cdots\; T(\mathbf{u}_n))$ |
| 표준 행렬 (Standard matrix) | $I$의 열을 입출력 기저로 사용: $T(\mathbf{u}) = A\mathbf{u}$ |
| 기저 변환 (Change of basis) | 새 행렬 = $B_{\text{out}}^{-1} A B_{\text{in}}$; 항등 변환은 $B = W^{-1}V$ |
| 미분 행렬 $D$ | $\{1, x, x^2, x^3\} \to \{1, x, x^2\}$: $D = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 3 \end{pmatrix}$ |
| 적분 행렬 $D^+$ | $D$의 유사역행렬; $DD^+ = I$이지만 $D^+D \neq I$ ($D$의 영공간) |
| 합성 $TS$ | 행렬곱 $AB$는 $TS$를 표현; $(m \times n)(n \times p) = (m \times p)$ |
| 회전 합성 | $\theta$만큼 회전 두 번 $= A^2 =$ $2\theta$만큼 회전 |
| 최적 기저: 고유벡터 | $B = X$ (고유벡터) $\Rightarrow$ $X^{-1}AX = \Lambda$ (대각) |
| 최적 기저: 특이벡터 | $B_{\text{in}} = V$, $B_{\text{out}} = U$ $\Rightarrow$ $U^{-1}AV = \Sigma$ |
| 최적 기저: 일반화 고유벡터 | $B^{-1}AB = J$ (Jordan 형식); $A$의 독립 고유벡터가 $n$개 미만일 때 사용 |
| Jordan 블록 $J_i$ | 대각에 $\lambda_i$, 상대각에 1; 독립 고유벡터 $s$개에 대해 $s$개 블록 |
| Jordan 연쇄 (chains) | $\mathbf{x}_1$: 고유벡터; $(A - \lambda I)\mathbf{x}_2 = \mathbf{x}_1$; $(A - \lambda I)\mathbf{x}_3 = \mathbf{x}_2$; 등 |
| Fourier 행렬 | 모든 순환 행렬(circulant matrix)을 대각화 |
| 함수 공간 기저 | Fourier ($1, \sin x, \cos x, \ldots$), Legendre ($1, x, x^2 - 1/3, \ldots$), Chebyshev ($1, x, 2x^2-1, \ldots$) |
| 함수의 Gram-Schmidt | 내적 $(f,g) = \int f\,g\,dx$; $1, x, x^2, \ldots$로부터 Legendre 다항식 생성 |
| 정규직교 기저 (Orthonormal) | $B^TB = I$; 계산에 가장 좋은 조건수의 기저 |

---
