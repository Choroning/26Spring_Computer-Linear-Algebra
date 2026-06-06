# 제6장 강의 — 고유값과 고유벡터 — Part 3/3

> **📑 이 문서는 3개 파트로 나뉘어 있습니다.**
>
> [Part 1](Concepts-Part1.ko.md) · [Part 2](Concepts-Part2.ko.md) · **Part 3**

---

<br>

## 목차

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

## 4. 선형 미분방정식 풀기 (6.5)

### 4.1 핵심 사실

**(1)** $`A\mathbf{x} = \lambda\mathbf{x}`$이면, $`\mathbf{u}(t) = e^{\lambda t}\mathbf{x}`$는 $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$를 풀어준다. 각 $`\lambda`$와 $`\mathbf{x}`$는 해 $`e^{\lambda t}\mathbf{x}`$를 준다.

**(2)** $`A = X\Lambda X^{-1}`$이면,

```math
\mathbf{u}(t) = e^{At}\mathbf{u}(0) = Xe^{\Lambda t}X^{-1}\mathbf{u}(0) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n
```

**(3)** 행렬 지수 함수: $`e^{At} = I + At + \cdots + \frac{(At)^n}{n!} + \cdots = Xe^{\Lambda t}X^{-1}`$ ($`A = X\Lambda X^{-1}`$인 경우).

**(4)** $`A`$의 모든 고유값의 **실수 부분이 $`< 0`$** 이면, $`A`$는 **안정**(stable)이고 $`\mathbf{u}(t) \to \mathbf{0}`$이며 $`e^{At} \to 0`$.

**(5)** $`\frac{d^2\mathbf{u}}{dt^2} + B\frac{d\mathbf{u}}{dt} + C\mathbf{u} = 0`$은 $`\frac{d}{dt}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix} = \begin{pmatrix} 0 & I \\ -C & -B \end{pmatrix}\begin{pmatrix} \mathbf{u} \\ \mathbf{v} \end{pmatrix}`$를 의미한다 ($`\mathbf{v} = \frac{d\mathbf{u}}{dt}`$).

### 4.2 스칼라 상미분방정식 복습

상미분방정식(ordinary differential equation)을 고려하자:

```math
\frac{du}{dt} = u \implies u(t) = Ce^t
```

```math
\boxed{\frac{du}{dt} = \lambda u} \implies \boxed{u(t) = Ce^{\lambda t}}
```

확인: $`\frac{du}{dt} = \lambda Ce^{\lambda t} = \lambda u`$.

$`t = 0`$일 때: $`u(0) = Ce^{\lambda \cdot 0} = C \cdot 1 = C`$. 초기값은 $`C`$이다.

$`\therefore u(t) = u(0)e^{\lambda t} = u_0 e^{\lambda t}`$

**거동:**

- $`\lambda > 0`$: **불안정**(unstable) (증가)
- $`\lambda = 0`$: **정상 상태**(steady state) (일정)
- $`\lambda < 0`$: **안정**(stable) (감쇠)

**$`\lambda \in \mathbb{C}`$이면?** 예: $`\lambda = -1 + 2i`$

$`e^{\lambda t} = e^{-t + 2it} = e^{-t}e^{2it}`$

$`|e^{\lambda t}| = |e^{-t}||e^{2it}| = e^{-t} \cdot \underbrace{1}_{|e^{i\theta}| = 1}`$

**관찰:** 안정성은 $`\lambda`$의 **실수 부분** 에 의존한다.

- $`\text{Re}(\lambda) > 0`$: 불안정
- $`\text{Re}(\lambda) = 0`$: 정상 (진동)
- $`\text{Re}(\lambda) < 0`$: 안정

### 4.3 du/dt = Au의 풀이

**예제 1:** 연립 상미분방정식(coupled ODE)을 고려하자:

```math
\frac{dy}{dt} = z, \quad \frac{dz}{dt} = y
```

```math
\implies \frac{d}{dt}\begin{pmatrix} y \\ z \end{pmatrix} = \begin{pmatrix} z \\ y \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y \\ z \end{pmatrix}
```

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u}
```

초기 조건 $`\mathbf{u}_0 = \begin{pmatrix} y(0) \\ z(0) \end{pmatrix}`$가 주어질 때, $`\mathbf{u}(t) = ?`$

$`\lambda`$를 $`A`$의 고유값, $`\mathbf{x}`$를 대응하는 고유벡터라 하자.

$`\mathbf{u} = e^{\lambda t}\mathbf{x}`$로 선택하고 연립 ODE에 대입하면:

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u} = e^{\lambda t}A\mathbf{x} = e^{\lambda t}\lambda\mathbf{x} = \lambda e^{\lambda t}\mathbf{x}
```

$`\mathbf{u} = e^{\lambda t}\mathbf{x}`$가 연립 ODE를 만족하므로, $`e^{\lambda t}\mathbf{x}`$는 $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$의 해이다.

```math
\implies \mathbf{u} = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2
```

$`A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}`$의 고유값은 무엇인가?

$`|A - \lambda I| = \begin{vmatrix} -\lambda & 1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0`$이므로, $`\lambda = \pm 1`$.

i) $`\lambda_1 = 1`$: $`(A - I)\mathbf{x}_1 = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$이므로, $`\mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

ii) $`\lambda_2 = -1`$: $`(A + I)\mathbf{x}_2 = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$이므로, $`\mathbf{x}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix}`$

iii) 완전한 해: $`\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2`$

iv) 상수 $`c_1`$, $`c_2`$는 초기 조건 $`\mathbf{u}_0 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}`$로 결정할 수 있다:

```math
\mathbf{u}(t=0) = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 = \begin{pmatrix} 4 \\ 2 \end{pmatrix}
```

$`\mathbf{u}_0`$을 $`A`$의 고유벡터 결합으로 쓰면:

```math
= (\mathbf{x}_1 \; \mathbf{x}_2)\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 4 \\ 2 \end{pmatrix}
```

```math
\begin{pmatrix} c_1 \\ c_2 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}^{-1}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = -\frac{1}{2}\begin{pmatrix} -1 & -1 \\ -1 & 1 \end{pmatrix}\begin{pmatrix} 4 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 1 \end{pmatrix}
```

### 4.4 일반적인 n x n 풀이 절차

$`n \times n`$ 행렬 $`A`$로 일반화:

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u}
```

$`A`$는 $`n`$개의 고유값과 $`n`$개의 고유벡터 ($`\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_n`$)를 갖는다.

1. $`\mathbf{u}_0`$을 $`\mathbf{u}_0 = c_1\mathbf{x}_1 + c_2\mathbf{x}_2 + \cdots + c_n\mathbf{x}_n`$으로 쓴다.
2. 고유벡터 $`\mathbf{x}_i`$에 $`e^{\lambda_i t}`$를 곱한다.
3. $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$의 해는:

```math
\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}
```

**예제 2:**

```math
\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 2 & 1 \\ 0 & 0 & 3 \end{pmatrix}\mathbf{u}, \quad \mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}
```

$`\lambda_1 = 1, \lambda_2 = 2, \lambda_3 = 3`$ (상삼각 행렬 — 대각 원소가 고유값).

i) $`(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 2 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}`$

ii) $`(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 & 1 \\ 0 & 0 & 1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}`$

iii) $`(A - 3I)\mathbf{x}_3 = \begin{pmatrix} -2 & 1 & 1 \\ 0 & -1 & 1 \\ 0 & 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix} \to \mathbf{x}_3 = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}`$

iv) $`\mathbf{u}_0 = \begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = (\mathbf{x}_1 \; \mathbf{x}_2 \; \mathbf{x}_3)\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix}`$

```math
\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{pmatrix}^{-1}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix}
```

$`|X| = 1`$, $`X^{-1} = \frac{1}{|X|}C^T = \begin{pmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & -1 & 1 \end{pmatrix}^T = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}`$

```math
\begin{pmatrix} c_1 \\ c_2 \\ c_3 \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{pmatrix}\begin{pmatrix} 9 \\ 7 \\ 4 \end{pmatrix} = \begin{pmatrix} 2 \\ 3 \\ 4 \end{pmatrix}
```

v) $`\therefore \mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + c_3 e^{\lambda_3 t}\mathbf{x}_3`$

```math
= 2e^t\begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + 3e^{2t}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix} + 4e^{3t}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}
```

### 4.5 행렬의 지수 함수

상기: $`e^x = 1 + x + \frac{1}{2}x^2 + \frac{1}{6}x^3 + \cdots + \frac{x^n}{n!} + \cdots = \sum_{i=0}^{\infty}\frac{x^i}{i!}`$

$`\frac{1}{1-x} = 1 + x + x^2 + x^3 + \cdots + x^n + \cdots = \sum_{i=0}^{\infty}x^i`$ (등비급수, $`|x| < 1`$)

**행렬 지수 함수(matrix exponential):**

```math
e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots + \frac{1}{n!}(At)^n + \cdots
```

$`(I - At)^{-1} = I + At + (At)^2 + (At)^3 + \cdots + (At)^n + \cdots`$

$`\mathbf{u}(t) = e^{At}\mathbf{u}(0) = e^{At}\mathbf{u}_0`$으로 놓자.

$`\mathbf{u}`$가 $`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$의 해인지 확인:

```math
\frac{d\mathbf{u}}{dt} = \frac{d}{dt}(e^{At}\mathbf{u}_0) = \frac{de^{At}}{dt}\mathbf{u}_0 = \left(A + \frac{1}{2} \cdot 2 \cdot (At)A + \frac{3}{6}(At)^2 A + \cdots + \frac{n}{n!}(At)^{n-1}A + \cdots\right)\mathbf{u}_0
```

```math
= A\left(I + At + \frac{1}{2}(At)^2 + \cdots + \frac{1}{(n-1)!}(At)^{n-1} + \cdots\right)\mathbf{u}_0 = Ae^{At}\mathbf{u}_0 = A\mathbf{u}
```

**$`A`$가 대각화 가능하면,** $`A = X\Lambda X^{-1}`$이고:

```math
e^{At} = Xe^{\Lambda t}X^{-1}
```

여기서

```math
e^{\Lambda t} = I + \Lambda t + \frac{1}{2}(\Lambda t)^2 + \frac{1}{6}(\Lambda t)^3 + \cdots = \begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}
```

**증명:**

```math
e^{At} = I + At + \frac{1}{2}(At)^2 + \cdots
```

```math
= I + X\Lambda X^{-1}t + \frac{1}{2}t^2(X\Lambda X^{-1})(X\Lambda X^{-1}) + \frac{1}{6}t^3(X\Lambda X^{-1})(X\Lambda X^{-1})(X\Lambda X^{-1}) + \cdots
```

```math
= I + t(X\Lambda X^{-1}) + \frac{1}{2}t^2(X\Lambda^2 X^{-1}) + \frac{1}{6}t^3(X\Lambda^3 X^{-1}) + \cdots
```

```math
= (XIX^{-1}) + (Xt\Lambda X^{-1}) + (X\frac{t^2}{2}\Lambda^2 X^{-1}) + (X\frac{t^3}{6}\Lambda^3 X^{-1}) + \cdots
```

```math
= X\left(I + t\Lambda + \frac{t^2}{2}\Lambda^2 + \frac{t^3}{6}\Lambda^3 + \cdots\right)X^{-1}
```

```math
\therefore e^{At} = Xe^{\Lambda t}X^{-1}
```

**관찰:** $`A = X\Lambda X^{-1}`$이고 $`e^{At} = Xe^{\Lambda t}X^{-1}`$ $`\implies`$ $`e^{At}`$와 $`A`$는 **동일한 고유벡터** 를 공유한다.

**해:**

```math
\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\mathbf{u}_0
```

$`\mathbf{u}_0 = X\mathbf{c}`$ (고유벡터의 선형결합)이므로:

```math
= Xe^{\Lambda t}\underbrace{X^{-1}X}_{I}\mathbf{c} = Xe^{\Lambda t}\mathbf{c} = X\begin{pmatrix} e^{\lambda_1 t} & & \\ & e^{\lambda_2 t} & \\ & & \ddots \\ & & & e^{\lambda_n t} \end{pmatrix}\begin{pmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{pmatrix} = X\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}
```

```math
= (\mathbf{x}_1 \; \mathbf{x}_2 \; \cdots \; \mathbf{x}_n)\begin{pmatrix} c_1 e^{\lambda_1 t} \\ c_2 e^{\lambda_2 t} \\ \vdots \\ c_n e^{\lambda_n t} \end{pmatrix}
```

```math
\boxed{\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{x}_n}
```

**성질:** $`e^{A(s+t)} = e^{As} \cdot e^{At}`$

**증명:** (이중 급수와 이항 정리를 사용한 상세 계산)

```math
e^{As} \cdot e^{At} = \left(\sum_{j=0}^{\infty}\frac{(As)^j}{j!}\right)\left(\sum_{k=0}^{\infty}\frac{(At)^k}{k!}\right) = \sum_{j=0}^{\infty}\sum_{k=0}^{\infty}\frac{A^{j+k}s^j t^k}{j!\,k!}
```

$`n = j + k`$, $`k = n - j`$로 놓으면:

```math
= \sum_{n=0}^{\infty}\frac{A^n}{n!}\sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^j t^{n-j} = \sum_{n=0}^{\infty}\frac{A^n}{n!}(s+t)^n = \sum_{n=0}^{\infty}\frac{(A(s+t))^n}{n!} = e^{A(s+t)} \quad \square
```

이항 정리(binomial theorem)를 사용: $`(s + t)^n = \sum_{j=0}^{n}\binom{n}{j}s^{n-j}t^j = \sum_{j=0}^{n}\frac{n!}{j!(n-j)!}s^{n-j}t^j`$.

$`e^{A(s+t)} = e^{As} \cdot e^{At}`$에서, $`s = -t`$로 놓으면:

```math
e^{A(-t+t)} = e^0 = I = e^{-At} \cdot e^{At}
```

이는 $`e^{-At}`$가 $`e^{At}`$의 **역행렬** 임을 의미한다.

**예제 4:** $`A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}`$, $`e^{At} = ?`$

$`|A - \lambda I| = \lambda^2 + 1 = 0`$이므로, $`\lambda = \pm i`$ (반대칭 행렬).

$`A^2 = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}`$

$`A^3 = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}`$

$`A^4 = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I`$

```math
e^{At} = I + At + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots
```

```math
= \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}t + \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix}\frac{t^2}{2} + \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\frac{t^3}{6} + \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\frac{t^4}{24} + \cdots
```

```math
= \begin{pmatrix} 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots & t - \frac{t^3}{6} + \cdots \\ -t + \frac{t^3}{6} - \cdots & 1 - \frac{t^2}{2} + \frac{t^4}{4!} - \cdots \end{pmatrix}
```

```math
e^{At} = \begin{pmatrix} \cos(t) & \sin(t) \\ -\sin(t) & \cos(t) \end{pmatrix} \quad \text{(반대칭 행렬은 회전을 준다)}
```

**예제 5:** $`\frac{d\mathbf{u}}{dt} = A\mathbf{u} = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}\mathbf{u}`$를 $`\mathbf{u}_0 = \begin{pmatrix} 2 \\ 1 \end{pmatrix}`$ ($`t = 0`$)에서 풀어라.

상삼각 행렬: $`\lambda_1 = 1, \lambda_2 = 2`$.

i) $`(A - I)\mathbf{x}_1 = \begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}`$

ii) $`(A - 2I)\mathbf{x}_2 = \begin{pmatrix} -1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

iii) $`X = (\mathbf{x}_1 \; \mathbf{x}_2) = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \to X^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}`$

$`A = X\Lambda X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & \\ & 2 \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}`$

$`e^{At} = Xe^{\Lambda t}X^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} e^t & \\ & e^{2t} \end{pmatrix}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}`$

iv) $`\mathbf{u}(t) = e^{At}\mathbf{u}_0 = Xe^{\Lambda t}X^{-1}\begin{pmatrix} 2 \\ 1 \end{pmatrix}`$

$`= Xe^{\Lambda t}\begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 2 \\ 1 \end{pmatrix} = Xe^{\Lambda t}\begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

$`= e^t\begin{pmatrix} 1 \\ 0 \end{pmatrix} + e^{2t}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \begin{pmatrix} e^t + e^{2t} \\ e^{2t} \end{pmatrix}`$

### 4.6 이계 방정식

$`my'' + by' + ky = 0`$ (질량-스프링-댐퍼 시스템, mass-spring-damper system)을 고려하자.

$`y = e^{\lambda t}`$로 놓으면: $`y' = \lambda e^{\lambda t} = \lambda y`$, $`y'' = (y')' = \lambda^2 y`$.

스프링 방정식은 다음이 된다: $`m\lambda^2 y + b\lambda y + ky = (m\lambda^2 + b\lambda + k)e^{\lambda t} = 0`$

```math
\implies m\lambda^2 + b\lambda + k = 0
```

```math
\iff \lambda^2 + \frac{b}{m}\lambda + \frac{k}{m} = 0 \quad \text{---(*)}
```

$`\lambda_1 + \lambda_2 = -\frac{b}{m}`$, $`\lambda_1\lambda_2 = \frac{k}{m}`$

$`\lambda = \frac{1}{2}\left(-\frac{b}{m} \pm \sqrt{\frac{b^2}{m^2} - \frac{4k}{m}}\right)`$

$`y_1 = e^{\lambda_1 t}`$, $`y_2 = e^{\lambda_2 t}`$로 놓으면, 완전한 해는:

```math
y(t) = c_1 y_1 + c_2 y_2 \quad \text{(} \lambda_1 \neq \lambda_2 \text{인 경우)} \quad \text{---(**)}
```

**이계 ODE를 일계 ODE로 변환:**

$`my'' + by' + ky = 0 \iff y'' + \tilde{b}y' + \tilde{k}y = 0`$ ($`\tilde{b} = \frac{b}{m}`$, $`\tilde{k} = \frac{k}{m}`$).

$`\mathbf{y} = \begin{pmatrix} y' \\ y \end{pmatrix}`$로 놓으면:

$`y'' = -\tilde{b}y' - \tilde{k}y = (-\tilde{b} \; -\tilde{k})\begin{pmatrix} y' \\ y \end{pmatrix}`$

$`y' = y' = (1 \; 0)\begin{pmatrix} y' \\ y \end{pmatrix}`$

```math
\implies \frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -\tilde{b} & -\tilde{k} \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix} \implies \frac{d\mathbf{y}}{dt} = A\mathbf{y}
```

i) $`A`$의 고유값:

$`|A - \lambda I| = \begin{vmatrix} -\tilde{b} - \lambda & -\tilde{k} \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + \tilde{b}\lambda + \tilde{k} = 0`$

$`\lambda`$에 대한 방정식은 (*)와 동일하다.

ii) $`\lambda_1, \lambda_2`$를 $`A`$의 고유값이라 하면:

$`(A - \lambda_1 I)\mathbf{x}_1 = \begin{pmatrix} -\tilde{b} - \lambda_1 & -\tilde{k} \\ 1 & -\lambda_1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}`$, $`x_1 - \lambda_1 x_2 = 0`$: $`\mathbf{x}_1 = \begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix}`$

$`(A - \lambda_2 I)\mathbf{x}_2`$: $`x_1 - \lambda_2 x_2 = 0`$: $`\mathbf{x}_2 = \begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}`$

iii) $`\mathbf{y}(t) = c_1 e^{\lambda_1 t}\mathbf{x}_1 + c_2 e^{\lambda_2 t}\mathbf{x}_2`$

```math
\begin{pmatrix} y'(t) \\ y(t) \end{pmatrix} = c_1 e^{\lambda_1 t}\begin{pmatrix} \lambda_1 \\ 1 \end{pmatrix} + c_2 e^{\lambda_2 t}\begin{pmatrix} \lambda_2 \\ 1 \end{pmatrix}
```

```math
\therefore y(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t} \quad \text{(**)와 동일}
```

$`m = 1, b = 0, k = 1`$일 때: $`y'' + y = 0 \iff y'' = -y`$

```math
\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -y \\ y' \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}
```

$`|A - \lambda I| = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1 = 0`$이므로, $`\lambda = \pm i`$ — **진동**(oscillation).

### 4.7 2 x 2 행렬의 안정성

$`\frac{d\mathbf{u}}{dt} = A\mathbf{u}`$의 해에 대해 근본적인 질문이 있다:

$`t \to \infty`$일 때 해가 $`\mathbf{u} = \mathbf{0}`$에 접근하는가?

완전한 해 $`\mathbf{u}(t)`$는 순수 해 $`e^{\lambda t}\mathbf{x}`$로 구성되므로, 안정성은 $`A`$의 고유값에 의존한다.

$`\lambda = r + is`$

$`e^{\lambda t} = e^{rt}e^{ist}`$

$`|e^{\lambda t}| = |e^{rt}|`$

$`\lambda`$의 **실수 부분** 이 증가 ($`r > 0`$) 또는 감쇠 ($`r < 0`$)를 제어한다.

임의의 $`2 \times 2`$ 행렬 $`A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}`$에 대해:

$`\lambda_1`$과 $`\lambda_2`$가 음수이려면:

i) $`\lambda_1 + \lambda_2 = a + d < 0`$ (대각합이 음수)

ii) $`\lambda_1 \lambda_2 > 0`$ (행렬식이 양수)

### 4.8 풀이 예제

**예제 6:** $`y'' + 4y' + 3y = 0`$을 풀어라.

$`y`$를 $`e^{\lambda t}`$로 대체: $`(\lambda^2 + 4\lambda + 3)e^{\lambda t} = 0`$

$`\implies \lambda^2 + 4\lambda + 3 = 0 \implies (\lambda + 3)(\lambda + 1) = 0`$

$`\therefore \lambda = -1, -3`$

$`\implies y(t) = c_1 e^{-t} + c_2 e^{-3t}`$ — 감쇠하는 해 $`\to`$ **안정** 해.

$`\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}`$를 도입하면: $`y'' = -4y' - 3y`$, $`y' = y'`$

```math
\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} -4 & -3 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad \frac{d\mathbf{u}}{dt} = A\mathbf{u}
```

$`|A - \lambda I| = \begin{vmatrix} -4 - \lambda & -3 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 4\lambda + 3 = 0`$이므로, $`\lambda = -1, -3`$.

대응하는 고유벡터:

$`\begin{pmatrix} -4+1 & -3 \\ 1 & +1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} -1 \\ 1 \end{pmatrix}`$

$`\begin{pmatrix} -4+3 & -3 \\ 1 & +3 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_2 = \begin{pmatrix} -3 \\ 1 \end{pmatrix}`$

완전한 해는:

```math
\begin{pmatrix} y' \\ y \end{pmatrix} = \mathbf{u}(t) = c_1 e^{-t}\begin{pmatrix} -1 \\ 1 \end{pmatrix} + c_2 e^{-3t}\begin{pmatrix} -3 \\ 1 \end{pmatrix}
```

이로부터 $`y(t) = c_1 e^{-t} + c_2 e^{-3t}`$.

**예제 7:** $`y'' - 2y' + y = 0`$

$`y`$를 $`e^{\lambda t}`$로 대체: $`(\lambda^2 - 2\lambda + 1) = (\lambda - 1)^2 = 0`$이므로, $`\lambda = 1, 1`$ (중근).

$`\implies y(t) = c_1 e^t + ?`$

근이 중복이므로, 두 번째 해는 $`c_2 \cdot t e^t`$이다.

```math
y(t) = c_1 e^t + c_2 t e^t
```

$`\mathbf{u} = \begin{pmatrix} y' \\ y \end{pmatrix}`$를 도입하면: $`y'' = 2y' - y`$, $`y' = y'`$

```math
\frac{d}{dt}\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} y' \\ y \end{pmatrix}, \quad A = \begin{pmatrix} 2 & -1 \\ 1 & 0 \end{pmatrix}
```

$`|A - \lambda I| = \begin{vmatrix} 2 - \lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 - 2\lambda + 1 = 0`$이므로, $`\lambda = 1, 1`$.

$`(A - I)\mathbf{x}_1 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \end{pmatrix} \to \mathbf{x}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix}`$

고유벡터의 수가 2보다 적으므로, 대각화는 **불가능** 하다.

따라서, 급수의 정의로부터 $`e^{At}`$를 계산한다:

```math
e^{At} = I + (At) + \frac{1}{2}(At)^2 + \frac{1}{6}(At)^3 + \cdots
```

$`A = I + (A - I)`$로 쓰면:

```math
e^{At} = e^{It + (A-I)t} = e^{It} \cdot e^{(A-I)t}
```

$`(A - I)^2 = \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}`$이므로, $`k \geq 2`$에 대해 $`(A - I)^k = 0`$.

```math
e^{(A-I)t} = I + (A-I)t = I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t
```

```math
e^{At} = e^t\left(I + \begin{pmatrix} 1 & -1 \\ 1 & -1 \end{pmatrix}t\right) = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}
```

```math
\mathbf{u}(t) = e^{At}\mathbf{u}_0
```

```math
\begin{pmatrix} y' \\ y \end{pmatrix} = \begin{pmatrix} e^t + te^t & -te^t \\ te^t & e^t - te^t \end{pmatrix}\begin{pmatrix} y_0' \\ y_0 \end{pmatrix}
```

```math
\therefore y(t) = (e^t - te^t)y_0 + te^t y_0'
```

---

<br>

## 요약

| 개념 | 핵심 아이디어 |
|:--------|:---------|
| 고유값 방정식 | $`A\mathbf{x} = \lambda\mathbf{x}`$, $`\mathbf{x} \neq \mathbf{0}`$ |
| 특성 다항식 | $`\det(A - \lambda I) = 0`$은 $`n`$개의 고유값을 준다 |
| 대각합과 행렬식 | $`\text{trace}(A) = \sum \lambda_i`$, $`\det(A) = \prod \lambda_i`$ |
| 고유값의 거듭제곱 | $`A^k\mathbf{x} = \lambda^k\mathbf{x}`$ |
| 역행렬의 고유값 | $`A^{-1}\mathbf{x} = \lambda^{-1}\mathbf{x}`$ ($`\lambda \neq 0`$인 경우) |
| 이동 규칙 | $`(A + cI)\mathbf{x} = (\lambda + c)\mathbf{x}`$ |
| 대칭 행렬 | 실수 고유값, 직교 고유벡터 |
| 반대칭 행렬 | 순허수 또는 영 고유값 |
| 회전 행렬 | $`\lambda = e^{\pm i\theta}`$ |
| 사영 행렬 | $`\lambda = 0`$ 또는 $`1`$만 가능 |
| 대각화 | $`n`$개 선형독립 고유벡터가 존재하면 $`A = X\Lambda X^{-1}`$ |
| 대각화를 통한 행렬 거듭제곱 | $`A^k = X\Lambda^k X^{-1}`$ |
| 서로 다른 고유값 | 고유벡터가 선형독립 $`\implies`$ 대각화 가능 |
| 닮은 행렬 | $`C = B^{-1}AB`$는 $`A`$와 같은 고유값을 가짐 |
| GM vs AM | GM $`\leq`$ AM; GM $`<`$ AM이면 대각화 불가능 |
| 스펙트럼 정리 | 대칭 $`S`$에 대해 $`S = Q\Lambda Q^T`$ |
| 양정치 | 모든 $`\mathbf{x} \neq \mathbf{0}`$에 대해 $`\mathbf{x}^T S\mathbf{x} > 0`$; 모든 $`\lambda_i > 0`$; 모든 피벗 $`> 0`$; 모든 선행 행렬식 $`> 0`$ |
| 양반정치 | $`\mathbf{x}^T S\mathbf{x} \geq 0`$; $`\lambda = 0`$ 허용 |
| $`A^T A`$의 에너지 판정 | $`\mathbf{x}^T A^T A\mathbf{x} = \|A\mathbf{x}\|^2 \geq 0`$ |
| 양정치 $`\implies`$ 가역 | $`\det(S) > 0`$, $`S^{-1}`$도 양정치 |
| 양정치 $`\implies S = A^T A`$ | $`A`$의 열이 독립 |
| 합동 행렬 | $`A^T SA`$는 정치성(definiteness) 유형을 보존 |
| 타원 방정식 | $`\mathbf{x}^T S\mathbf{x} = 1`$; 축은 고유벡터 방향, 길이 $`1/\sqrt{\lambda}`$ |
| 헤시안과 볼록성 | $`H`$ 양정치 $`\implies`$ $`f`$는 순볼록 |
| 경사 하강법 | $`\mathbf{x}_{k+1} = \mathbf{x}_k - \eta\nabla f(\mathbf{x}_k)`$ |
| ODE 해 | $`\frac{d\mathbf{u}}{dt} = A\mathbf{u} \implies \mathbf{u}(t) = \sum c_i e^{\lambda_i t}\mathbf{x}_i`$ |
| 행렬 지수 함수 | $`e^{At} = I + At + \frac{(At)^2}{2!} + \cdots = Xe^{\Lambda t}X^{-1}`$ |
| 안정성 (연속) | 모든 $`\text{Re}(\lambda_i) < 0`$이면 안정 |
| 안정성 (이산) | 모든 $`|\lambda_i| < 1`$이면 $`A^k \to 0`$ |
| $`e^{A(s+t)} = e^{As}e^{At}`$ | 행렬 지수 함수의 가법 성질 |
| $`(e^{At})^{-1} = e^{-At}`$ | 행렬 지수 함수의 역행렬 |
| 이계 ODE | $`my'' + by' + ky = 0 \iff`$ 일계 시스템 $`\frac{d\mathbf{y}}{dt} = A\mathbf{y}`$ |
| 중복 고유값 ODE | 대각화 불가 경우: $`e^{At} = e^{It} \cdot e^{(A-I)t}`$ 사용 |

---
