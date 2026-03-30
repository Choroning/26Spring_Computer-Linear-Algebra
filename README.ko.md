# [2026학년도 봄학기] 전산선형대수학

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Computer-Linear-Algebra)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Computer-Linear-Algebra)

이 레포지토리는 대학 강의 및 과제를 위해 작성된 과제 풀이와 학습 자료를 체계적으로 정리하고 보관합니다.

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

| 주차 | 챕터 | 주제 | 핵심 개념 |
|:---:|:---|:---|:---|
| 1 | Chap. 1 | 벡터 소개 | 벡터 |
| 2 | Chap. 2 | 연립일차방정식 풀기 | Ax=b |
| 3 | Chap. 3 | 연립일차방정식 풀기 | Ax=b |
| 4 | Chap. 3.1-3.2 | 네 가지 기본 부분공간 | 열공간 |
| 5 | Chap. 3.3-3.5 | 네 가지 기본 부분공간 | 영공간 |
| 6 | Chap. 4.1-4.3 | 직교성 | 직교성 |
| 7 | Chap. 4.4-4.5 | 직교성 | 직교성 |
| 8 | 중간고사 | — | — |
| 9 | Chap. 5 | 행렬식 | 행렬식 |
| 10 | Chap. 6.1-6.3 | 고유값과 고유벡터 | 고유 분해 |
| 11 | Chap. 6.4-6.5 | 고유값과 고유벡터 | 고유 분해 |
| 12 | Chap. 7 | 특이값 분해 | SVD |
| 13 | Chap. 8 | 선형 변환 | 선형 사상 |
| 14 | Chap. 8 | 선형 변환 | 선형 사상 |
| 15 | Chap. 9 | 최적화에서의 선형대수 | SGD |
| 16 | 기말고사 | — | — |

<br><a name="repository-structure"></a>
## 🗂 레포지토리 구조

```plaintext
26Spring_Computer-Linear-Algebra
├── .gitignore
├── LICENSE
├── README.ko.md
└── README.md
```

<br><a name="license"></a>
## 🤝 라이선스

이 레포지토리는 [MIT License](LICENSE) 하에 배포됩니다.

---
