# [Spring 2026] Computer Linear Algebra

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Computer-Linear-Algebra)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Computer-Linear-Algebra)

This repository archives bilingual lecture notes and coursework for the Computer Linear Algebra course.

*Author: Cheolwon Park (Korea University Sejong, CSE) - Year 3 (Junior) as of 2026*
<br><br>

## 📑 Table of Contents

- [About This Repository](#about-this-repository)
- [Course Information](#course-information)
- [Prerequisites](#prerequisites)
- [Repository Structure](#repository-structure)
- [License](#license)

---


<br><a name="about-this-repository"></a>
## 📝 About This Repository

This repository contains the bilingual lecture notes of a university-level Computer Linear Algebra course that follows Gilbert Strang's *Introduction to Linear Algebra* and covers vectors, linear equations, vector spaces, orthogonality, determinants, eigenvalues, SVD, and linear transformations.

- `Lecture-Notes/` holds lecture notes on the concepts of every lecture topic. The lecture deck numbers (N1-N8) match the textbook chapters, so the units are also numbered by textbook chapter (`Ch##`), covering Chapters 1-8.
- Every lecture note comes as an English file (`.md`) and a Korean file (`.ko.md`) with the same structure.
- The solutions to the weekly homework assignments are not uploaded, so there is no `Coursework/` folder.

> **🤖 AI Agent Usage**
> This course encourages the use of AI agents.
> [Claude Code](https://claude.ai/download) and [Gemini CLI](https://github.com/google-gemini/gemini-cli) were used as coding assistants throughout the course, and Claude Code was also used to organize and translate the lecture notes and to restructure the repository.

<br><a name="course-information"></a>
## 📚 Course Information

- **Semester:** Spring 2026 (March - June)
- **Affiliation:** Korea University Sejong

| Course&nbsp;Code| Course            | Type          | Instructor      | Department                              |
|:----------:|:------------------|:-------------:|:---------------:|:----------------------------------------|
|`DCSS321-01`|COMPUTER LINEAR ALGEBRA (English)|Major Elective|Prof. Shinhoo&nbsp;Kang|Department of Computer Science and Software Engineering|

### Course Overview

An English-taught course on the fundamental concepts of linear algebra, including vector spaces, linear transformations, matrices, systems of linear equations, and eigenvalues. Following Gilbert Strang's textbook, it covers vectors and matrices, solving linear equations, vector spaces and subspaces, orthogonality, determinants, eigenvalues and eigenvectors, the singular value decomposition (SVD), and linear transformations in turn, and aims to develop the ability to bridge abstract theory with computational techniques such as matrix factorizations and eigenvalue decomposition to solve problems in science, engineering, and data analysis. The instructor, schedule and class format, assessment, course policies, weekly plan and progress, and references are in the [course information document](Course-Information.md).

<br><a name="prerequisites"></a>
## ✅ Prerequisites

- An understanding of Discrete Mathematics is required.
- A Python interpreter must be installed.
- Familiarity with scientific computing libraries is expected.

- **💻 Development Environment**

| Tool | Company |  OS  | Notes |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|Numerical computing|
|Matplotlib|NumFOCUS|macOS|Visualization|
|JupyterLab|Project Jupyter|macOS|Interactive notebooks|

<br><a name="repository-structure"></a>
## 🗂 Repository Structure

```plaintext
26Spring_Computer-Linear-Algebra
├── Lecture-Notes               Lecture notes by textbook chapter
├── Course-Information.ko.md    Course operations
├── Course-Information.md
├── LICENSE                     Documentation license (CC BY-NC-SA 4.0)
├── LICENSE-CODE                Code license (PolyForm Noncommercial 1.0.0)
├── README.ko.md
└── README.md
```

The units are listed in the [`Lecture-Notes` index](Lecture-Notes/README.md).

<br><a name="license"></a>
## 🤝 License

The documents in this repository (lecture notes, explanations, and summary sheets) are licensed under [CC BY-NC-SA 4.0](LICENSE), and the code under [PolyForm Noncommercial 1.0.0](LICENSE-CODE). Figures and example code taken from lecture slides and textbooks, as well as assignment handouts and template code distributed by the course, belong to their original authors and are not covered by these licenses.

---
