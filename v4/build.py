#!/usr/bin/env python3
"""Markdown → LaTeX 結合スクリプト。本書専用の最小 MD パーサ。"""
import re
from pathlib import Path

HERE = Path(__file__).parent
CHAPTERS = [
    "00_to_reader.md",
    "00_preliminaries.md",
    "01_markowitz_basics.md",
    "02_efficient_frontier.md",
    "03_two_fund_separation.md",
    "04_capm.md",
    "05_apt_factor_models.md",
    "06_utility_theory.md",
    "07_stochastic_dominance.md",
    "08_risk_measures.md",
    "09_black_litterman.md",
    "10_continuous_time_merton.md",
    "11_performance_evaluation.md",
    "12_shrinkage_rmt.md",
    "13_robust_optimization.md",
    "14_risk_parity_hrp.md",
    "15_factor_zoo_ml.md",
    "16_backtest_overfit.md",
    "17_exercises.md",
]

LATEX_ESCAPES = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def escape_text(s: str) -> str:
    out = []
    for ch in s:
        out.append(LATEX_ESCAPES.get(ch, ch))
    return "".join(out)


def protect_math(line: str):
    """Replace $...$ and $$...$$ regions with placeholders so escaping skips them."""
    placeholders = []

    def repl_display(m):
        placeholders.append(("display", m.group(1)))
        return f"\x00MATHD{len(placeholders)-1}\x00"

    def repl_inline(m):
        placeholders.append(("inline", m.group(1)))
        return f"\x00MATHI{len(placeholders)-1}\x00"

    line = re.sub(r"\$\$(.+?)\$\$", repl_display, line, flags=re.S)
    line = re.sub(r"\$([^\$\n]+?)\$", repl_inline, line)
    return line, placeholders


def restore_math(line: str, placeholders):
    def repl(m):
        idx = int(m.group(2))
        kind, content = placeholders[idx]
        if kind == "display":
            return "\\[" + content + "\\]"
        return "$" + content + "$"

    return re.sub(r"\x00MATH(I|D)(\d+)\x00", repl, line, flags=re.S)


def process_inline(s: str) -> str:
    """Handle inline formatting AFTER math protection: bold/italic/code/links."""
    # ``code``
    s = re.sub(r"`([^`]+)`", lambda m: r"\texttt{" + escape_text(m.group(1)) + "}", s)
    # **bold** / __bold__
    s = re.sub(r"\*\*(.+?)\*\*", lambda m: r"\textbf{" + m.group(1) + "}", s)
    # *italic*
    s = re.sub(r"(?<!\*)\*([^*\n]+?)\*(?!\*)", lambda m: r"\emph{" + m.group(1) + "}", s)
    # [text](url) -> \href{url}{text}
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", lambda m: r"\href{" + m.group(2) + r"}{" + m.group(1) + "}", s)
    return s


def convert_paragraph(text: str) -> str:
    text, ph = protect_math(text)
    text = escape_text(text)
    text = process_inline(text)
    text = restore_math(text, ph)
    return text


def convert_table(lines):
    """Convert a markdown table (list of lines including header & separator) to LaTeX tabular."""
    rows = [l.strip() for l in lines if l.strip()]
    # Drop the separator row (---)
    cells = []
    for r in rows:
        if re.match(r"^\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)+\|?$", r):
            continue
        # Split by |, trim outer
        parts = r.strip("|").split("|")
        cells.append([c.strip() for c in parts])
    ncol = max(len(r) for r in cells)
    col_spec = "|" + "l|" * ncol
    out = [r"\begin{center}", r"\begin{tabular}{" + col_spec + "}", r"\hline"]
    for i, row in enumerate(cells):
        row = row + [""] * (ncol - len(row))
        out.append(" & ".join(convert_paragraph(c) for c in row) + r" \\ \hline")
    out += [r"\end{tabular}", r"\end{center}"]
    return "\n".join(out)


IMAGE_RE = re.compile(r"^\s*!\[([^\]]*)\]\(([^)]+)\)\s*$")


def convert_md(md: str, chapter_label: str) -> str:
    lines = md.splitlines()
    i = 0
    out = []
    in_list = False
    list_type = None  # 'itemize' or 'enumerate'
    in_quote = False

    def close_list():
        nonlocal in_list, list_type
        if in_list:
            out.append("\\end{" + list_type + "}")
            in_list = False
            list_type = None

    def close_quote():
        nonlocal in_quote
        if in_quote:
            out.append("\\end{quote}")
            in_quote = False

    while i < len(lines):
        line = lines[i]
        # Skip nav links lines like "[← ...] | [次 → ...]"
        if re.match(r"^\s*\[←.*\]\(.*\)\s*｜\s*\[次.*\]\(.*\)\s*$", line) or re.match(
            r"^\s*\[次章.*\]\(.*\)\s*$", line
        ) or re.match(r"^\s*\[←.*\]\(.*\)\s*$", line) or re.match(r"^\s*\[.*トップ.*\]\(.*\)\s*$", line):
            i += 1
            continue
        # horizontal rule
        if re.match(r"^---+\s*$", line):
            close_list()
            close_quote()
            out.append(r"\noindent\hrulefill")
            i += 1
            continue
        # standalone image ![caption](path)
        m_img = IMAGE_RE.match(line)
        if m_img:
            close_list()
            close_quote()
            caption = m_img.group(1).strip()
            path = m_img.group(2).strip()
            out.append(r"\begin{figure}[H]")
            out.append(r"\centering")
            out.append(r"\includegraphics[width=0.85\linewidth]{" + path + "}")
            if caption:
                out.append(r"\caption{" + convert_paragraph(caption) + "}")
            out.append(r"\end{figure}")
            i += 1
            continue
        # headers
        h = re.match(r"^(#{1,6})\s+(.*)$", line)
        if h:
            close_list()
            close_quote()
            level = len(h.group(1))
            title = convert_paragraph(h.group(2))
            if level == 1:
                # Each file's first '#' becomes a chapter
                out.append("\\chapter{" + title + "}")
            elif level == 2:
                out.append("\\section{" + title + "}")
            elif level == 3:
                out.append("\\subsection{" + title + "}")
            else:
                out.append("\\subsubsection{" + title + "}")
            i += 1
            continue
        # blockquote
        if line.startswith(">"):
            close_list()
            if not in_quote:
                out.append(r"\begin{quote}")
                in_quote = True
            content = line.lstrip(">").strip()
            out.append(convert_paragraph(content))
            i += 1
            continue
        else:
            close_quote()

        # Table detection: a row starting with | followed by a separator row
        if line.strip().startswith("|") and i + 1 < len(lines) and re.match(
            r"^\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)+\|?$", lines[i + 1].strip()
        ):
            close_list()
            tbl = [line]
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith("|"):
                tbl.append(lines[j])
                j += 1
            out.append(convert_table(tbl))
            i = j
            continue

        # Unordered list
        ul = re.match(r"^(\s*)-\s+(.*)$", line)
        ol = re.match(r"^(\s*)(\d+)\.\s+(.*)$", line)
        if ul:
            content = convert_paragraph(ul.group(2))
            if not in_list:
                out.append(r"\begin{itemize}")
                in_list = True
                list_type = "itemize"
            elif list_type != "itemize":
                out.append("\\end{" + list_type + "}")
                out.append(r"\begin{itemize}")
                list_type = "itemize"
            out.append(r"\item " + content)
            i += 1
            continue
        if ol:
            content = convert_paragraph(ol.group(3))
            if not in_list:
                out.append(r"\begin{enumerate}")
                in_list = True
                list_type = "enumerate"
            elif list_type != "enumerate":
                out.append("\\end{" + list_type + "}")
                out.append(r"\begin{enumerate}")
                list_type = "enumerate"
            out.append(r"\item " + content)
            i += 1
            continue

        close_list()

        # Blank line
        if line.strip() == "":
            out.append("")
            i += 1
            continue

        # Plain paragraph (consume until blank or block boundary)
        para_lines = [line]
        j = i + 1
        while j < len(lines):
            nxt = lines[j]
            if (
                nxt.strip() == ""
                or re.match(r"^#{1,6}\s", nxt)
                or re.match(r"^(\s*)-\s+", nxt)
                or re.match(r"^(\s*)\d+\.\s+", nxt)
                or nxt.startswith(">")
                or re.match(r"^---+\s*$", nxt)
                or (nxt.strip().startswith("|") and j + 1 < len(lines)
                    and re.match(r"^\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)+\|?$", lines[j + 1].strip()))
            ):
                break
            para_lines.append(nxt)
            j += 1
        para = " ".join(l.strip() for l in para_lines)
        out.append(convert_paragraph(para))
        i = j

    close_list()
    close_quote()
    return "\n".join(out)


PREAMBLE = r"""\documentclass[11pt,a4paper]{report}
\usepackage{luatexja-preset}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage{bm}
\usepackage{graphicx}
\usepackage{float}
\usepackage[margin=25mm]{geometry}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{array}
\hypersetup{colorlinks=true,linkcolor=blue!60!black,urlcolor=blue!50!black,citecolor=blue!50!black}
\setcounter{tocdepth}{1}
\setcounter{secnumdepth}{3}

\newcommand{\E}{\mathbb{E}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Var}{\mathrm{Var}}
\newcommand{\Cov}{\mathrm{Cov}}
\newcommand{\one}{\mathbf{1}}

\title{Investment Science \\ \large マルコヴィッツから連続時間最適化まで}
\author{}
\date{Version 4.0 — 学部 4 年向け版 \\ \today}

\begin{document}
\maketitle
\tableofcontents
\clearpage
"""

POSTAMBLE = r"""
\end{document}
"""


def main():
    parts = [PREAMBLE]
    for c in CHAPTERS:
        md = (HERE / c).read_text(encoding="utf-8")
        parts.append("% ===== " + c + " =====")
        parts.append(convert_md(md, c))
        parts.append("")
    parts.append(POSTAMBLE)
    out_path = HERE / "book.tex"
    out_path.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
