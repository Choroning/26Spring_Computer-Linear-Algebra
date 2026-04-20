# Chapter 3 복습 시트 — 벡터 공간과 부분공간 (Vector Spaces and Subspaces)

> Strang, *Introduction to Linear Algebra* 6th Ed., Chapter 3
> 다루는 절: 3.1 벡터 공간·부분공간 · 3.2 영공간 계산 · 3.3 완전해 · 3.4 독립·기저·차원 · 3.5 4대 부분공간

---

## 공통 기호 (Ch3 관련 발췌)

| 기호 | 뜻 | 어느 공간에 사는가 |
|:---|:---|:---|
| $C(A)$ | 열공간 (column space) | $\mathbb{R}^m$ |
| $C(A^T)$ | 행공간 (row space) | $\mathbb{R}^n$ |
| $N(A)$ | 영공간 (null space) | $\mathbb{R}^n$ |
| $N(A^T)$ | 좌영공간 (left null space) | $\mathbb{R}^m$ |
| $r$ | 랭크 (rank) | — |
| $\mathbf{x}_p$ | **특수해** (particular solution) | $\mathbb{R}^n$ |
| $\mathbf{x}_n$ | **영공간 해** (null-space part) | $N(A)$ |
| $\mathbf{s}_j$ | **특수해** (special solution) — $N(A)$의 기저 | $N(A)$ |
| $R_0 = \mathrm{rref}(A)$ | 기약행사다리꼴 (reduced row echelon form) | — |
| 자유변수 (free variable) | 피벗이 없는 열의 변수 | — |

> ⚠ 주의: "**특수해**"가 두 가지 뜻 — 한글은 같지만 영어는 다름.
> - $\mathbf{x}_p$ = **particular solution** ($A\mathbf{x}=\mathbf{b}$의 한 해)
> - $\mathbf{s}_j$ = **special solution** ($A\mathbf{x}=\mathbf{0}$의 기저 벡터)

---

## 3-A. 벡터 공간 8공리 (the 8 axioms of a vector space)

$V$가 **벡터 공간 (vector space)** 이려면, $\mathbf{u},\mathbf{v},\mathbf{w}\in V$, $c,d\in\mathbb{F}$에 대해:

**덧셈 (addition)**:
1. 결합법칙 (associativity): $\mathbf{u}+(\mathbf{v}+\mathbf{w})=(\mathbf{u}+\mathbf{v})+\mathbf{w}$
2. 교환법칙 (commutativity): $\mathbf{u}+\mathbf{v}=\mathbf{v}+\mathbf{u}$
3. 영벡터 존재 (zero vector): $\exists\,\mathbf{0}$ s.t. $\mathbf{v}+\mathbf{0}=\mathbf{v}$
4. 덧셈 역원 존재 (additive inverse): $\exists\,-\mathbf{v}$ s.t. $\mathbf{v}+(-\mathbf{v})=\mathbf{0}$

**스칼라곱 (scalar multiplication)**:
5. $c(d\mathbf{v})=(cd)\mathbf{v}$
6. $1\mathbf{v}=\mathbf{v}$ (항등원, identity)
7. $c(\mathbf{u}+\mathbf{v})=c\mathbf{u}+c\mathbf{v}$ (벡터 덧셈에 대한 분배)
8. $(c+d)\mathbf{v}=c\mathbf{v}+d\mathbf{v}$ (스칼라 덧셈에 대한 분배)

## 3-B. 부분공간 판정 (subspace test) — 두 조건만

- (i) $\mathbf{u}+\mathbf{w}\in V$ — 덧셈에 닫힘 (closed under addition)
- (ii) $c\mathbf{u}\in V$ — 스칼라곱에 닫힘 (closed under scalar multiplication)
- → **영벡터 $\mathbf{0}$ 포함은 자동** (스칼라곱 조건에서 $c=0$)

### ⭐ 실전 빠른 판별 패턴 (문제 3.1.9 같은 거 30초 컷)

| 집합 정의 방식 | 부분공간? | 이유 |
|:---|:---:|:---|
| **등식 $= 0$** ($b_1 = b_2$, $b_1+b_2+b_3=0$) | ✅ 거의 YES | 원점 포함 + 선형성 보존 |
| **등식 $=$ 상수** ($b_1 = 1$, $z = 5$) | ❌ NO | **원점 미포함** |
| **부등식** ($b_1 \leq b_2$, $x \geq 0$) | ❌ NO | 음수배로 부등호 뒤집힘 |
| **곱 $= 0$** ($b_1 b_2 b_3 = 0$) | ❌ NO | 덧셈 깨짐 (축들의 합집합) |
| **"모든 선형결합"** (span) | ✅ 항상 YES | 정의상 닫힘 |

> 🔑 **원점 체크 = 첫 관문**. $\mathbf{0}$이 안 속하면 바로 NO. 시간 절약.

**부분공간인 것 / 아닌 것 예시**:
- ✅ 원점 지나는 직선/평면
- ✅ $\{\mathbf{0}\}$, $\mathbb{R}^n$ 자체
- ❌ $\mathbb{R}^3$의 평면 $z=5$ (원점 미포함)
- ❌ 제1사분면 (스칼라곱 $c<0$에서 닫히지 않음)
- ❌ 양수 벡터의 집합 (덧셈 역원 없음)
- ❌ 정수 벡터 집합 ($0.5 \times$ 정수 = 정수 아님)
- ❌ 두 축의 합집합 ($(1,0)+(0,1)=(1,1)$ 어느 축도 아님)

**$\mathbb{R}^n$ 외의 벡터공간 예**:
- **행렬 공간 (matrix space)**: 모든 $n\times n$ 행렬 → 차원 $n^2$
- **함수 공간 (function space)**: 예) $y''=0$의 해공간, 기저 $\{1, x\}$, 차원 2

### $C(A)$, $N(A)$가 어디 사는지 — 절대 헷갈리지 말 것

$A \in \mathbb{R}^{m \times n}$일 때:

$$\boxed{C(A) \subseteq \mathbb{R}^m, \quad N(A) \subseteq \mathbb{R}^n}$$

- $C(A)$: 열벡터 길이 = $m$ → $\mathbb{R}^m$
- $N(A)$: $x \in Ax=0$에서 $x$ 성분 수 = $n$ → $\mathbb{R}^n$

## 3-C. 영공간 $N(A)$ 계산 — 기약행사다리꼴 (rref)

**$A\mathbf{x}=\mathbf{0}$의 모든 해**를 구하는 게 목표.

### rref의 정확한 4가지 조건 (기억할 것)

1. 각 피벗이 **1**
2. 피벗 **위아래 전부 0** (echelon form보다 강한 조건)
3. 피벗이 왼쪽→오른쪽, 위→아래 **계단형** 배치
4. **영행은 맨 아래**

> ⚠ Ch2의 $U$ (echelon form) 는 "피벗 아래만 0". rref는 "위아래 다 0" + "피벗 = 1".

### 핵심 불변량 ⭐

**소거해도 $N(A)$는 안 바뀐다**: $N(A) = N(R_0) = N(R)$.
- 증명 핵심: $R = EA$, $E$ 가역 → $Ax=0 \iff EAx=0$
- 그래서 편한 $R$에서 $N$을 구해도 됨

### 계산 절차

1. $A$를 소거해서 $R_0 = \mathrm{rref}(A)$로 만듦
2. $R_0$은 $r$개의 **피벗 열 (pivot columns)** 과 $n-r$개의 **자유 열 (free columns)** 을 가짐
3. 순열하면 $R = (I\ F)$ 형태

### 특수해 (special solutions) 구하는 법

- **자유변수 중 하나만 1**, 나머지 자유변수 0으로 놓고 $R\mathbf{x}=\mathbf{0}$ 풀기
- 이렇게 $n-r$개 나옴 → 이들이 **$N(A)$의 기저**
- 한꺼번에 표현: 특수해의 행렬은 $\binom{-F}{I}$의 열들

$$\boxed{(I\ F)\binom{-F}{I} = O}$$

> 💡 **왜 자유변수에 "1"을 넣나?** 아무 값이나 넣어도 되지만, 1을 넣으면 패턴이 $\binom{-F}{I}$로 깔끔해지고 기저 벡터들이 **자동으로 독립** (자유변수 자리에 단위벡터 꼴이라서). 계산 편의 + 독립 보장.

### Rank-Nullity 관계

$$\boxed{\text{rank}(A) + \text{nullity}(A) = n}$$

- $\text{rank}(A) = r$ = 피벗 수 = 독립 열 수
- $\text{nullity}(A) = n - r$ = 자유변수 수 = 특수해 수 = $\dim N(A)$
- **열 개수 = 피벗 열 + 자유 열 = $r + (n-r) = n$**. 항상 성립.

## 3-D. $A\mathbf{x}=\mathbf{b}$ 완전해 (complete solution)

**완전해 공식 — 항상 이 구조**:

$$\boxed{\mathbf{x} = \mathbf{x}_p + \mathbf{x}_n,\qquad A\mathbf{x}_p=\mathbf{b},\ A\mathbf{x}_n=\mathbf{0}}$$

- **$\mathbf{x}_p$ (특수해, particular solution)**: 모든 자유변수를 **0으로** 놓고 $R_0\mathbf{x}=\mathbf{d}$ 풀기
- **$\mathbf{x}_n$ (영공간 부분, null-space part)**: 특수해 (special solutions)들의 결합, 즉 $N(A)$의 임의 벡터
- **풀림 조건 (solvability)**: $\mathbf{b} \in C(A)$ ⟺ $R_0$의 영행에 대응하는 $\mathbf{d}$ 성분이 모두 0

### ⭐ 6단계 표준 절차 (암기 후 기계적으로)

1. **Augmented matrix** $(A \mid \mathbf{b})$ 만들기
2. **rref로 소거** → $(R_0 \mid \mathbf{d})$
3. **해 존재 확인**: $R_0$의 영행 자리에 $\mathbf{d}$도 0인지 체크 (아니면 **해 없음**, 끝)
4. **$\mathbf{x}_p$ 구하기**: 자유변수 = 0 대입 → 피벗변수 값 계산
5. **$\mathbf{x}_n$ 구하기**: 특수해들 (3.2 방식) 찾기
6. **완전해 조합**: $\mathbf{x} = \mathbf{x}_p + c_1\mathbf{s}_1 + c_2\mathbf{s}_2 + \cdots$

> 💡 **왜 $\mathbf{x}_p$는 자유변수에 0을 넣나?** 아무 값이나 가능. 0이 계산 제일 깔끔.

### 해의 존재·유일성 — 두 축 분리

| 질문 | 무엇으로 판단? |
|:---|:---|
| 해가 **존재**하는가? | $\mathbf{b} \in C(A)$ 인가? |
| 해가 **유일**한가? | $N(A) = \{\mathbf{0}\}$ 인가? |
| 해 집합의 **차원**? | $\dim N(A) = n - r$ |

- **Full column rank** ($r = n$) → 해 있을 때 유일
- **Full row rank** ($r = m$) → 모든 $\mathbf{b}$에 해 존재

### 4가지 경우 판정표 (Ch3 일반판)

| 경우 | 형태 | $R_0$ 모양 | 해의 수 |
|:---|:---|:---|:---|
| $r=m=n$ | 정방·가역 (square, invertible) | $I$ | **유일 (unique)** |
| $r=m<n$ | 짧고 넓음 (short-wide, full row rank) | $(I\ F)$ | 항상 존재, **무한히 많음 ($\infty$)** |
| $r=n<m$ | 크고 좁음 (tall-thin, full column rank) | $\binom{I}{0}$ | **0개 또는 1개** |
| $r<m,\,r<n$ | 랭크 부족 (rank deficient) | $\binom{I\ F}{0\ 0}$ | **0개 또는 무한히 많음** |

## 3-E. 선형독립·생성·기저·차원 (independence · span · basis · dimension)

- **선형독립 (linearly independent, LI)**: $c_1\mathbf{v}_1+\cdots+c_n\mathbf{v}_n=\mathbf{0}$의 유일 해가 모든 $c_i=0$
  → 동치: $A\mathbf{x}=\mathbf{0}$의 유일 해가 $\mathbf{x}=\mathbf{0}$, 즉 $N(A)=\{\mathbf{0}\}$
- **생성 (span)**: 주어진 벡터들의 모든 선형결합의 집합
- **기저 (basis)**: 독립 + 생성 → 모든 벡터가 **유일하게** 표현됨
- **차원 (dimension)**: 기저의 크기 (어떤 기저를 뽑든 불변)

**중요 사실**:
- $\mathbb{R}^m$에서 $n > m$이면 $n$개 벡터는 **반드시 종속 (dependent)**
- 모든 가역 $n\times n$ 행렬의 열들은 $\mathbb{R}^n$의 기저

## 3-F. 4대 부분공간 & 기본 정리 ★ (Fundamental Theorem) ★

$A \in \mathbb{R}^{m\times n}$, 랭크 $r$일 때:

| 부분공간 (subspace) | 정의 | 사는 곳 | 차원 (dimension) |
|:---|:---|:---|:---:|
| **열공간 (column space) $C(A)$** | $\{A\mathbf{x}\}$ | $\mathbb{R}^m$ | $r$ |
| **행공간 (row space) $C(A^T)$** | $\{A^T\mathbf{y}\}$ | $\mathbb{R}^n$ | $r$ |
| **영공간 (null space) $N(A)$** | $A\mathbf{x}=\mathbf{0}$ | $\mathbb{R}^n$ | $n-r$ |
| **좌영공간 (left null space) $N(A^T)$** | $A^T\mathbf{y}=\mathbf{0}$ | $\mathbb{R}^m$ | $m-r$ |

**직교 여공간 쌍 (orthogonal complement pairs)** — 기본 정리의 핵심:

$$\underbrace{C(A^T) \perp N(A)}_{\mathbb{R}^n\text{에서 직교 여공간}}\qquad \underbrace{C(A) \perp N(A^T)}_{\mathbb{R}^m\text{에서 직교 여공간}}$$

$$r + (n-r) = n,\qquad r + (m-r) = m$$

### ⭐ 소거의 불변량 — 왜 기저 뽑는 법이 공간마다 다른가

**핵심 원리**: 소거는 **행 연산**만 함. 그래서:

| 공간 | 소거 후 보존? | 이유 |
|:---|:---:|:---|
| 행공간 $C(A^T)$ | ✅ 보존 | 행끼리 조작해도 span은 그대로 |
| 영공간 $N(A)$ | ✅ 보존 | $Ax=0 \iff EAx=0$ (가역) |
| 열공간 $C(A)$ | ❌ **바뀜** | 열의 실제 벡터가 변함 |
| 좌영공간 $N(A^T)$ | ❌ **바뀜** | $A^T$ 기준으로는 행 연산 아님 |

**구체 예시**:
$A = \begin{pmatrix}1&2&4\\2&4&8\end{pmatrix}$ 의 열 $(1,2)$ 방향 직선 $\ne$ $R_0 = \begin{pmatrix}1&2&4\\0&0&0\end{pmatrix}$ 의 열 $(1,0)$ 방향 = $x$축.

### 기저 뽑는 법 (how to extract a basis) — 기계적 규칙

- **$C(A)$** → **원래 $A$의 피벗 열** (⚠ $R_0$에서 위치만 확인, 벡터는 $A$에서!)
- **$N(A)$** → $R_0$에서 **특수해 (special solutions)**, $\binom{-F}{I}$의 열들
- **$C(A^T)$** → $R_0$의 **영이 아닌 행** (보존되니 편한 쪽에서)
- **$N(A^T)$** → $A^T$를 따로 소거해서 특수해

> 🔑 **한 줄 요약**: "보존되는 공간은 $R$에서, 변하는 공간은 $A$ (혹은 $A^T$) 에서 뽑는다."

### $A$와 $R_0$의 관계 요약

| 공유하는 것 | 안 공유하는 것 |
|:---|:---|
| 행공간 $C(A^T) = C(R_0^T)$ | 열공간 $C(A) \ne C(R_0)$ (차원만 같음) |
| 영공간 $N(A) = N(R_0)$ | 좌영공간 $N(A^T) \ne N(R_0^T)$ |
| rank $r$ | |

## 3-G. Ch3 자주 헷갈리는 것

| ❌ 헷갈림 | ✅ 바로잡기 |
|:---|:---|
| 제1사분면, $z=5$ 평면은 부분공간 | **아니다** (원점 미포함 / 스칼라곱에 안 닫힘) |
| 소거하면 열공간도 보존됨 | **아니다**. 행공간·영공간만 보존. 열공간은 **바뀜** |
| $C(A)$ 기저를 $R_0$의 피벗 열로 뽑기 | **$A$ 자체**의 피벗 열로 뽑아야 함 |
| $\mathbf{x}_p$와 특수해 (special sol.)가 같음 | 한글은 같지만 다른 개념! $\mathbf{x}_p$는 $A\mathbf{x}=\mathbf{b}$의 한 해, special은 $N(A)$의 기저 벡터 |
| $A\mathbf{x}=\mathbf{0}$의 유일 해는 항상 $\mathbf{0}$ | $r<n$이면 **비자명 해 (nontrivial)** 있음 |
| 기저가 하나뿐 | **무한히 많음**. 차원만 불변 |
| $n > m$인 $n$개 벡터도 독립일 수 있음 | **항상 종속** |
| $b$를 $2\times 3$ 행렬 같은 걸로 생각 | $\mathbf{b}$는 **항상 벡터 하나** ($m \times 1$). $\mathbf{x}$도 벡터 ($n \times 1$) |
| "rank 2면 평면이니 해 무한" (해 개수를 $C(A)$로 판단) | 해 개수는 **$N(A)$로** 판단. $C(A)$는 해 **존재** 판정용 |
| "왜 어떨 땐 행, 어떨 땐 열로 뽑지?" | 소거는 **행 연산**이라 행 관련 공간만 보존. 열 관련은 원래 $A$ 필요 |
| 독립 = "두 벡터 사이 스칼라배 아님" | 벡터 3개 이상이면 부족. **$c_1v_1+\cdots+c_nv_n=0$의 유일 해가 모두 0** 이 정확한 정의 |
| 영벡터 포함한 집합도 독립일 수 있음 | **무조건 종속** ($1\cdot\mathbf{0} + 0\cdot v = \mathbf{0}$) |
| 완전해에서 해 존재 판단을 원래 $A$로 | 소거 후 $(R_0 \mid \mathbf{d})$에서 영행 대응 $\mathbf{d}$ 성분이 0인지 보면 됨 |

## 3-H. 용어 미묘한 차이

- **피벗 (pivot)** vs **피벗 열 (pivot column)**: 피벗 = 수 (숫자), 피벗 열 = 열 (벡터)
- **$\mathbf{x}_p$ (particular solution)** vs **$\mathbf{s}_j$ (special solution)**:
  - $\mathbf{x}_p$: $A\mathbf{x}=\mathbf{b}$ 만족, 자유변수 전부 0
  - $\mathbf{s}_j$: $A\mathbf{x}=\mathbf{0}$의 기저, 자유변수 하나만 1
- **특이 (singular)** = 비가역 / **정칙 (nonsingular)** = 가역
- **동차 (homogeneous)** = 우변 $\mathbf{0}$ / **비동차 (nonhomogeneous)** = 우변 $\mathbf{0}$ 아님
- **부정 (underdetermined)** = $n>m$ (보통 해 $\infty$) / **과결정 (overdetermined)** = $m>n$ (보통 해 없음)
- **echelon form ($U$)** (Ch2) vs **reduced row echelon form ($R_0$)** (Ch3):
  - $U$: 피벗 **아래**만 0. 피벗이 1일 필요 없음
  - $R_0$: 피벗 **위아래** 전부 0. 피벗 **= 1** 필수
- **$C(A)$ vs $C(R_0)$**: 차원은 같지만 **공간 자체는 다름** (소거로 열이 변하니까)

## 3-I. Ch3 체크리스트

- [ ] 벡터 공간 8공리를 떠올릴 수 있는가?
- [ ] 부분공간 판정 두 조건?
- [ ] 부분공간 빠른 판별 패턴 (등식/상수/부등식/곱/span)?
- [ ] $R_0 = \mathrm{rref}(A)$를 직접 구할 수 있는가? (rref 4조건)
- [ ] 왜 자유변수에 "1"을 넣는지 설명할 수 있는가?
- [ ] 특수해 (special solutions)로 $N(A)$의 기저 뽑기?
- [ ] $\mathbf{x}_p + \mathbf{x}_n$ 형태로 완전해 쓰기 (6단계 절차)?
- [ ] $A\mathbf{x}=\mathbf{b}$ 해 존재 조건 ($R_0$의 영행에 $\mathbf{d}$도 0)?
- [ ] 4가지 경우 판정표?
- [ ] 4대 부분공간의 차원 $(r, r, n-r, m-r)$?
- [ ] 4대 부분공간의 기저 각각 뽑는 법? (왜 $C(A)$만 원래 $A$에서?)
- [ ] 소거의 불변량 (보존: $C(A^T), N(A)$ / 변함: $C(A), N(A^T)$)?
- [ ] 직교 여공간 2쌍?
- [ ] 기본 정리 (Fundamental Theorem of Linear Algebra)를 한 문장으로?

---

## 🎯 한 줄 마무리 (the one-liner)

> **선형대수는 결국 두 질문으로 귀결된다:**
> **"$\mathbf{b}$가 $C(A)$에 있는가?"** 와 **"$N(A)$는 얼마나 큰가?"**
> (Linear algebra boils down to: "Is $\mathbf{b}$ in $C(A)$?" and "How big is $N(A)$?")
> 나머지 — 소거법, $LU$, 전치, 직교성 — 는 전부 이 두 질문에 답하기 위한 도구.

---

## 부록. Ch3 영문 용어 색인 (English Glossary)

| English | 한글 | 비고 / 기호 |
|:---|:---|:---|
| basis | 기저 | 독립 + 생성 |
| column space | 열공간 | $C(A)$ |
| complete solution | 완전해 | $\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n$ |
| dimension | 차원 | 기저의 크기 |
| echelon form | 사다리꼴 | |
| equivalent | 동치 | $\iff$ |
| free variable | 자유변수 | 피벗 없는 변수 |
| full column rank | 열 완전 계수 | $r=n$ |
| full row rank | 행 완전 계수 | $r=m$ |
| Fundamental Theorem of Linear Algebra | 선형대수학의 기본 정리 | |
| homogeneous | 동차 | $A\mathbf{x}=\mathbf{0}$ |
| independent (linearly) | 선형독립 | LI |
| left null space | 좌영공간 | $N(A^T)$ |
| null space | 영공간 | $N(A)$ |
| orthogonal complement | 직교 여공간 | $N(A)\perp C(A^T)$ |
| particular solution | 특수해 (해 하나) | $\mathbf{x}_p$ |
| rank | 랭크 (계수) | $r$ |
| reduced row echelon form (rref) | 기약행사다리꼴 | $R_0$ |
| row space | 행공간 | $C(A^T)$ |
| solvable | 풀림 가능 | $\mathbf{b}\in C(A)$ |
| span | 생성 | 모든 선형 결합 |
| special solution | 특수해 ($N(A)$의 기저) | $N(A)$의 기저 벡터 |
| subspace | 부분공간 | |
| trivial solution | 자명해 | $\mathbf{x}=\mathbf{0}$ |
| unique solution | 유일 해 | $\exists!$ |
| vector space | 벡터 공간 | 8공리 |
| zero vector | 영벡터 | $\mathbf{0}$ |
