# [2026학년도 봄학기] 전산선형대수학

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Computer-Linear-Algebra)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Computer-Linear-Algebra)

이 레포지토리는 전산선형대수학 과목의 이중 언어 강의 노트와 수업 결과물을 정리한 아카이브입니다.

*작성자: 박철원 (고려대학교(세종), 컴퓨터소프트웨어학과) - 2026년 기준 3학년*
<br><br>

## 📑 목차

- [레포지토리 소개](#about-this-repository)
- [강의 정보](#course-information)
- [사전 요구사항](#prerequisites)
- [레포지토리 구조](#repository-structure)
- [라이선스](#license)

---


<br><a name="about-this-repository"></a>
## 📝 레포지토리 소개

이 레포지토리에는 Gilbert Strang의 *Introduction to Linear Algebra*를 따라 벡터, 연립일차방정식, 벡터 공간, 직교성, 행렬식, 고유값, SVD, 선형 변환을 다루는 대학 수준의 전산선형대수학 과목의 이중 언어 강의 노트가 담겨 있습니다.

- `Lecture-Notes/`에는 강의 주제마다 개념을 정리한 강의 노트가 있습니다. 강의자료의 번호(N1\~N8)가 교재의 장 번호와 같아 단원도 교재의 장 번호(`Ch##`)로 나눴으며, 1\~8장을 다룹니다.
- 모든 강의 노트는 영어 파일(`.md`)과 한국어 파일(`.ko.md`)이 같은 구성으로 짝을 이룹니다.
- 매주 나온 과제의 풀이는 올리지 않았으므로 `Coursework/` 폴더는 없습니다.

> **🤖 AI 에이전트 활용**
> 이 과목은 AI 에이전트 사용을 권장합니다.
> 수업 전반에 걸쳐 [Claude Code](https://claude.ai/download)와 [Gemini CLI](https://github.com/google-gemini/gemini-cli)를 코딩 어시스턴트로 활용했고, 강의 노트의 정리와 번역, 레포지토리 구조 정리에도 Claude Code를 활용했습니다.

<br><a name="course-information"></a>
## 📚 강의 정보

- **학기:** 2026학년도 봄학기 (3월 - 6월)
- **소속:** 고려대학교(세종)

|학수번호      |강의명    |이수구분|교수자|개설학과|
|:----------:|:-------|:----:|:------:|:----------------|
|`DCSS321-01`|전산선형대수학(영강)|전공선택|강신후 교수|컴퓨터소프트웨어학과|

### 과목 개요

벡터 공간, 선형 변환, 행렬, 연립일차방정식, 고유값 같은 선형대수학의 기본 개념을 다루는 영어 강의 과목이다. Gilbert Strang의 교재를 따라 벡터와 행렬, 연립일차방정식, 벡터 공간과 부분공간, 직교성, 행렬식, 고유값과 고유벡터, 특이값 분해(SVD), 선형 변환을 차례로 다루며, 추상적인 이론을 행렬 분해와 고유값 분해 같은 계산 기법과 연결해 과학, 공학, 데이터 분석의 문제를 푸는 능력을 기르는 것을 목표로 한다. 교수자, 수업 일정과 형식, 평가, 수업 운영 규정, 주차별 계획 및 진도, 참고 자료는 [강의 정보 문서](Course-Information.ko.md)에 정리했다.

<br><a name="prerequisites"></a>
## ✅ 사전 요구사항

- 이산수학을 이해하고 있어야 합니다.
- Python 인터프리터가 설치되어 있어야 합니다.
- 과학 계산 라이브러리에 익숙해야 합니다.

- **💻 개발 환경**

| 도구 | 회사 |  운영체제  | 비고 |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|수치 연산|
|Matplotlib|NumFOCUS|macOS|시각화|
|JupyterLab|Project Jupyter|macOS|대화형 노트북|

<br><a name="repository-structure"></a>
## 🗂 레포지토리 구조

```plaintext
26Spring_Computer-Linear-Algebra
├── Lecture-Notes               교재 장별 강의 노트
├── Course-Information.ko.md    수업 운영 정보
├── Course-Information.md
├── LICENSE                     문서 라이선스 (CC BY-NC-SA 4.0)
├── LICENSE-CODE                코드 라이선스 (PolyForm Noncommercial 1.0.0)
├── README.ko.md
└── README.md
```

단원 목록은 [`Lecture-Notes` 색인](Lecture-Notes/README.ko.md)에 있습니다.

<br><a name="license"></a>
## 🤝 라이선스

이 레포지토리의 문서(강의 노트, 해설, 요약지)는 [CC BY-NC-SA 4.0](LICENSE), 코드는 [PolyForm Noncommercial 1.0.0](LICENSE-CODE)에 따라 배포됩니다. 강의자료와 교재에서 가져온 그림과 예제 코드, 강의가 배포한 과제 원본과 템플릿 코드는 원저작자에게 권리가 있으며 이 라이선스의 대상이 아닙니다.

---
