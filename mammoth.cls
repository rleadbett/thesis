\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mammoth}[The cozy mammoth]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% KOMA %

% This template is based of KOMA-scrbook
% Change all KOMA settings here

% Page setup settings
\PassOptionsToPackage{
    DIV  = 12,
    BCOR = 10mm,
}{typearea}

% Load KOMA
\LoadClass[
    a4paper,
    twoside,
    12pt,
    cleardoublepage=empty,
    numbers=noenddot,
    titlepage,
    toc=bibliography,
    toc=index,
]{scrbook}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compiler and Font 

% Automatically determines compiler and uses correct fonter
\RequirePackage{iftex}
\newif\ifxetexorluatex\xetexorluatexfalse
\ifxetex
    \xetexorluatextrue
\fi
\ifluatex
    \xetexorluatextrue
\fi

\ifxetexorluatex
    \RequirePackage{fontspec}
    \RequirePackage{unicode-math}
    \RequirePackage{microtype}
    \setmainfont{Libertinus Serif}
    \setmathfont{Libertinus Math}
    % Load some missing symbols from another font.
    \setmathfont{STIX Two Math}[%
        range = {
                \checkmark
            }
    ]
    \setsansfont[Scale=MatchLowercase]{Libertinus Sans}
    \setmonofont[Scale=MatchLowercase]{Roboto Mono}

\else
    \RequirePackage[T1]{fontenc}
    \RequirePackage[utf8]{inputenc}
    \RequirePackage{amssymb}
\fi

\renewcommand{\baselinestretch}{1.1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Referencing

\RequirePackage{csquotes}
\RequirePackage[english]{babel}
\RequirePackage[style=apa, natbib=true]{biblatex}
\AtEveryBibitem{%
    \clearfield{note}%
    \clearfield{url}%
}

%\renewcommand{\cite}[1]{\textcite{#1}}
%\newcommand{\citep}[1]{\parencite{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Colours

\RequirePackage[dvipsnames]{xcolor}
\RequirePackage[table]{xcolor}
\definecolor{accent}{HTML}{B48C09}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Graphics

\RequirePackage{graphicx}
\RequirePackage{caption}
\RequirePackage{subcaption}
\RequirePackage{tikz}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Tables

\RequirePackage{booktabs}
\RequirePackage{makecell}
\RequirePackage{multirow}
\RequirePackage{longtable}
\RequirePackage{multicol}
\RequirePackage{tabularx}
\RequirePackage{array}

\newcolumntype{Y}{>{\raggedright\arraybackslash}X}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hyperlinks & Bookmarks
\RequirePackage[%
    colorlinks = true,
    citecolor  = accent,
    linkcolor  = accent,
    urlcolor   = accent,
    unicode,
]{hyperref}
\RequirePackage{bookmark}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mathematics

\RequirePackage{amsmath}
\RequirePackage{amsthm}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Algorithms

\RequirePackage{amsmath}
\RequirePackage{algorithm}
\RequirePackage{algpseudocode}
%\DontPrintSemicolon
%\SetKwProg{Function}{function}{:}{}
%\SetKw{Continue}{continue}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Headers

\RequirePackage{scrlayer-scrpage}
\pagestyle{scrheadings}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Utils

\RequirePackage{enumitem}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Footnotes

%% Make bar full length
\renewcommand{\footnoterule}{\noindent\smash{\rule[3pt]{\textwidth}{0.4pt}}}

%% Macro for chapter with footnote
\newcommand{\chapterfootnote}[2]{\chapter[#1]{#1\footnote{#2}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Typefaces for parts, chapters, and sections

\renewcommand{\partformat}{\huge\partname~\thepart\autodot}
\renewcommand{\raggedpart}{\flushleft}

\setkomafont{part}{\normalfont\huge\scshape}

\setkomafont{sectioning}{\normalfont\scshape}
\setkomafont{descriptionlabel}{\normalfont\bfseries}

\setkomafont{caption}{\small}
\setkomafont{captionlabel}{\usekomafont{caption}}

% Large number for chapter
\renewcommand*{\chapterformat}{%
    \color{accent}\fontsize{50}{55}\selectfont\thechapter\autodot\enskip
}

% color the chapter numbers
\let\oldsectionformat\sectionformat
\renewcommand*{\sectionformat}{\color{accent}\oldsectionformat}
\let\oldsubsectionformat\subsectionformat
\renewcommand*{\subsectionformat}{\color{accent}\oldsubsectionformat}
\let\oldsubsubsectionformat\subsubsectionformat
\renewcommand*{\subsubsectionformat}{\color{accent}\oldsubsubsectionformat}

% Ad subsubsections
\setcounter{secnumdepth}{\subsubsectionnumdepth}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Environments

%% SC title
\newcommand{\titlesc}[1]{\textsc{\Huge #1}}

%% Prepages
\newenvironment{prepage}[1]{
    \begin{center}
        \titlesc{#1}
    \end{center}
    \quotation
    \noindent\ignorespaces
}{
    \endquotation
    \clearpage
}
