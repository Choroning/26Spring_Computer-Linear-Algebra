# [2026학년도 봄학기] 전산선형대수학

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Computer-Linear-Algebra)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Computer-Linear-Algebra)

이 레포지토리는 대학 강의 및 과제를 위해 작성된 강의 노트와 학습 자료를 체계적으로 정리하고 보관하며, 벡터, 연립일차방정식, 벡터 공간, 직교성, 행렬식, 고유값, SVD, 선형 변환을 다룹니다.

*작성자: 박철원 (고려대학교(세종), 컴퓨터소프트웨어학과) - 2026년 기준 3학년*
<br><br>

## 📑 목차

- [레포지토리 소개](#about-this-repository)
- [강의 정보](#course-information)
- [사전 요구사항](#prerequisites)
- [주차별 일정](#weekly-schedule)
- [레포지토리 구조](#repository-structure)
- [라이선스](#license)

---


<br><a name="about-this-repository"></a>
## 📝 레포지토리 소개

이 레포지토리에는 대학 수준의 전산선형대수학 과목을 위해 작성된 이중 언어 학습 자료와 과제 풀이가 포함되어 있습니다:

- 매 강의 주제별 이중 언어 개념 정리 노트 (한국어 `.ko.md` + 영어 `.md`)
- 주간 과제 풀이
- Gilbert Strang 교재 기반 Chapter 디렉토리 구조

> **🤖 AI 에이전트 활용**
> 본 과목은 AI 에이전트 사용을 권장합니다.
> 수업 전반에 걸쳐 [Claude Code](https://claude.ai/download)와 [Gemini CLI](https://github.com/google-gemini/gemini-cli)를 코딩 어시스턴트로 활용하였습니다.

<br><a name="course-information"></a>
## 📚 강의 정보

- **학기:** 2026학년도 봄학기 (3월 - 6월)
- **소속:** 고려대학교(세종)

|학수번호      |강의명    |이수구분|교수자|개설학과|
|:----------:|:-------|:----:|:------:|:----------------|
|`DCSS321-01`|전산선형대수학(영강)|전공선택|강신후 교수|컴퓨터소프트웨어학과|

- **📖 참고 자료**

| 유형 | 내용 |
|:----:|:---------|
|교재|"Introduction to Linear Algebra" by Gilbert Strang (6th Edition, Wellesley-Cambridge Press, 2023)|

<br><a name="prerequisites"></a>
## ✅ 사전 요구사항

- 이산수학에 대한 이해
- Python 인터프리터 설치
- 과학 계산 라이브러리 숙지

- **💻 개발 환경**

| 도구 | 회사 |  운영체제  | 비고 |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|수치 연산|
|Matplotlib|NumFOCUS|macOS|시각화|
|JupyterLab|Project Jupyter|macOS|대화형 노트북|

<br><a name="weekly-schedule"></a>
## 📅 주차별 일정

| 주차 | 강의 자료 | 챕터 | 주제 |
|:---:|:---|:---|:---|
| 1 | Introduction, N1 | Ch. 1 | 소개, 벡터와 행렬 |
| 2 | — | Ch. 1 | 벡터와 행렬 (계속) |
| 3 | N2 | Ch. 2 | 연립일차방정식 풀기 |
| 4 | — | Ch. 2 | 연립일차방정식 풀기 (계속) |
| 5 | N3 (A) | Ch. 3 | 벡터 공간과 부분공간 |
| 6 | N3 (B) | Ch. 3 | 벡터 공간과 부분공간 (계속) |
| 7 | — | — | 복습 |
| 8 | 중간고사 | — | — |
| 9 | N4-1, N4-2 | Ch. 4 | 직교성, 사영 |
| 10 | N4-3, N4-4, N4-5 | Ch. 4 | 최소제곱, 직교 기저, 의사역행렬 |
| 11 | N5-1, N5-2, N5-3 | Ch. 5 | 행렬식 |
| 12 | N6-1, N6-2 | Ch. 6 | 고유값, 대각화 |
| 13 | N6-3, N6-5 | Ch. 6 | 대칭 양정치 행렬, 선형 미분방정식 |
| 14 | N7-1, N7-2, N7-3 | Ch. 7 | SVD, 이미지 처리, PCA |
| 15 | N8-1, N8-2, N8-3 | Ch. 8 | 선형 변환 |
| 16 | 기말고사 | — | — |

<br><a name="repository-structure"></a>
## 🗂 레포지토리 구조

```plaintext
26Spring_Computer-Linear-Algebra
├── Chapter01_Introduction-to-Vectors
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter02_Solving-Linear-Equations
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter03_Vector-Spaces-and-Subspaces
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter04_Orthogonality
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter05_Determinants
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter06_Eigenvalues-and-Eigenvectors
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter07_Singular-Value-Decomposition
│   ├── Concepts.md
│   └── Concepts.ko.md
├── Chapter08_Linear-Transformations
│   ├── Concepts.md
│   └── Concepts.ko.md
├── docs/
├── .gitignore
├── LICENSE
├── README.ko.md
└── README.md
```

<br><a name="license"></a>
## 🤝 라이선스

이 레포지토리는 [MIT License](LICENSE) 하에 배포됩니다.

---
