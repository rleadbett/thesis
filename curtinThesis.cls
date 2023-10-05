%%
%% --------------------------------------------------------------------------------
%%
%% 
%%   This work is licensed under the conditions of the Creative
%%   Commons Attribution 4.0 International (CC BY 4.0) license:
%%
%%		<https://creativecommons.org/licenses/by/4.0/>
%% 
%%   This work may be remixed, tweaked, and build upon, even
%%	 commercially, as long as the original Author is credited
%%   for the creation of the original work. 
%% 
%%   This code is offered as-is without any warranty either expressed
%%   or implied; User assumes all risk.
%%
%% --------------------------------------------------------------------------------
%% 
%%	Functionality
%%  -------------
%%  This class facilitates a Curtin thesis style template format
%%  which should serve as an aid to generating your Curtin Thesis in
%%	LaTeX.
%%
%%	Several class options are available to change the template based
%%  on the requirements. These options are:
%%
%%		-\DeclareOption{12pt} 		-- Set the font-size to 12pt - loaded by default
%%		-\DeclareOption{colorlinks}	-- Use the colorlinks options in xcolor
%%		-\DeclareOption{emptypage}	-- Use the emptypage package to remove page numbering on empty pages
%%
%%		-\DeclareOption{apa}		-- Use the APA bibliography style
%%		-\DeclareOption{chicago}	-- Use the Chicago bibliography style
%%		-\DeclareOption{vancouver}	-- Use the Vancouver bibliography style
%%		-\DeclareOption{harvard}	-- Use the Harvard bibliography style
%%
%%		Furthermore, all book class options can be used.
%%
%%	Packages used in this class are:
%%		
%%		standard graphics, math and biliographic packages:
%%		-amsfonts,amsmath,amssymb,graphicx,xcolor,pdfpages,[optional] emptypage,[either] apacite or natbib
%%		
%%		get appendices to work nicely
%%		-appendix
%%
%%		get bibliography to appear in TOC
%%		-tocbibind
%%
%%		packages for customising margins and spacings
%%		-setspace,xspace
%%
%%		table, array, float manipulation packages
%%		-array,multicol,rotfloat
%%
%%		code formatting package
%%		-listings
%%
%%		pdf configuration and link formatting packages
%%		-hyperref
%%
%%	This class also uses several custom functions to generate the thesis based on the template.
%%	The functions defined are:
%%
%%		globalvariables
%%		-\thesistitle			-- global variable to define thesis title
%%		-\authorname			-- global variable to define author name
%%		-\auxassoc				-- [optional] global variable to define auxiliary association
%%		-\graduationdate		-- global variable to define graduation date
%%		-\thesisquote			-- [optional] global variable to define thesis quote
%%
%%		markup variables
%%		-\linkcolor				-- [optional] markup variable to define the link color used
%%		-\citecolor				-- [optional] markup variable to define the cite color 
%%		
%%		markup functions
%%		-\setmetadata			-- function to generate the pdf metadata based on your global variables
%%		-\createtitlepage		-- function to generate the title page based on your global variables
%%		-\plagiarismstatement	-- function to generate a standard plagiarism statement required in the thesis
%%		-\showthesisquote		-- function to generate the thesis quote if \thesisquote has been used
%%		-\ackstatement			-- function to generate a standard acknowledgement regarding acknowledging authors in references
%%		-\includetoc			-- function to include a table of contents to the thesis
%%		-\includelof			-- function to include a list of figures to the thesis
%%		-\addapptotoc			-- function to include the appendices to the table of contents
%%
%%		formatting functions
%%		-\arabicnumbering{n}	-- changes the page numbering to arabic starting at number 'n'
%%		-\alphnumbering{n}		-- changes the page numbering to alphabetic starting at number 'n'
%%
%%		environments
%%		-\begin{customappendix}	-- starts the custom appendix environment - used to properly markup the appendices
%%		-\end{customappendix}	-- ends the custom appendix environment
%%
%% --------------------------------------------------------------------------------
%%		

\NeedsTeXFormat{LaTeX2e}[2017/04/15]
\ProvidesClass{curtinThesis}[2018/11/07 Curtin Thesis Class]

%% class markup declaration options %%
\DeclareOption{12pt}{
	\PassOptionsToClass{\CurrentOption}{book}
}
\newif\if@colorlinks\@colorlinksfalse
\DeclareOption{colorlinks}{
	\@colorlinkstrue
}
\newif\if@emptypage\@emptypagefalse
\DeclareOption{emptypage}{
	\@emptypagetrue
}

%% class bibliography definition options %%
\newif\if@APA\@APAfalse
\DeclareOption{apa}{
	\@APAtrue	
}
\newif\if@Chicago\@Chicagofalse
\DeclareOption{chicago}{
	\@Chicagotrue
}
\newif\if@Vancouver\@Vancouverfalse
\DeclareOption{vancouver}{
	\@Vancouvertrue
}
\newif\if@Harvard\@Harvardfalse
\DeclareOption{harvard}{
	\@Harvardtrue
}

%% class parse undefined options %%
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}

%% execute default options %%
\ExecuteOptions{12pt}

%% process given options %%
\ProcessOptions\relax

%% load the base class definition %%
\LoadClass[a4paper]{book}

%% load the packages required to properly style the thesis %%

% standard graphics, math and bibliographic packages
\RequirePackage{amsfonts}
\RequirePackage{amsmath}
\RequirePackage{amssymb}
\RequirePackage{graphicx}
\RequirePackage[table]{xcolor}
\RequirePackage{pdfpages}

\if@emptypage
	\RequirePackage{emptypage}
\fi

\if@APA
\RequirePackage[natbibapa]{apacite} 
\bibliographystyle{apacite}
\fi
\if@Vancouver
	\RequirePackage[numbers]{natbib}
	\bibliographystyle{unsrtnat}
\else
	\RequirePackage{natbib}
	\if@Chicago
		\bibliographystyle{chicago}
	\fi
	\if@Harvard
		\bibliographystyle{harvard}
	\fi
\fi


% get appendices to work nicely
\RequirePackage[toc,page]{appendix}

% get bibliography to appear in TOC
\RequirePackage[notindex,nottoc,notlof]{tocbibind} 

% packages for customising margins and spacings
\RequirePackage{setspace}
\RequirePackage{xspace}

% table, array, float manipulation packages
\RequirePackage{array}
\RequirePackage{multicol}
\RequirePackage{rotfloat}

% code formatting package
\RequirePackage{listings}
\lstset{frame=tb,
	aboveskip=3mm,
	belowskip=3mm,
	basicstyle={\small\ttfamily}
}

% highlighting
\usepackage{xcolor}
\usepackage{soul}

% load hyperref last to ensure it overrides commands necessary
\usepackage{hyperref} 

%% custom functions %%

% global variables
\def\thesistitle#1{\def\@thesistitle{#1}}
\def\authorname#1{\def\@authorname{#1}}
\def\auxassoc#1{\def\@auxassoc{#1}}\auxassoc{}
\def\graduationdate#1{\def\@graduationdate{#1}}
\def\thesisquote#1#2{\def\@quote{#1}\def\@quoteauthor{#2}}\thesisquote{}{}

% markup variables
\if@colorlinks\def\@colorlinks{true}\else\def\@colorlinks{false}\fi
\def\linkcolor#1{\def\@linkcolor{#1}}\linkcolor{}
\def\citecolor#1{\def\@citecolor{#1}}\citecolor{}

% metadata
\newcommand{\setmetadata}{%
	\ifx\@linkcolor\@empty\linkcolor{cyan}\fi
	\ifx\@citecolor\@empty\citecolor{cyan}\fi
	\makeatletter
	\hypersetup{
		pdfsubject=PhD Thesis,
		pdftitle={\@thesistitle},
		pdfauthor={\@authorname},
		pdfcreator={\@authorname},
		unicode=true,
		colorlinks={\@colorlinks},
		linkcolor={\@linkcolor},
		citecolor={\@citecolor},
		linktoc=page}
	\makeatother
}%

%% page styling %%
\pagestyle{plain}
\pagenumbering{roman}

% fix margins 
\renewcommand{\baselinestretch}{1.5}
\setlength{\voffset}{-1in}
\setlength{\hoffset}{-1in}
\setlength{\oddsidemargin}{4cm}
\setlength{\evensidemargin}{2.5cm}
\setlength{\textwidth}{14.5cm}
\setlength{\textheight}{22.5cm}
\setlength{\topmargin}{2.5cm}

%% define how to punctuate citations %%
\bibpunct{(}{)}{;}{a}{,}{,}

%% define titlepage %%
\newcommand{\createtitlepage}{%
	\newif\if@assoc\@assocfalse
	\ifx\@auxassoc\@empty\@assocfalse\else\@assoctrue\fi
	\begin{titlepage}
		\begin{center}
			\vspace{\stretch{0.3}}
			%%% School %%%
			{%
				\fontsize{12}{14.4} 
				\textbf{School of Electrical Engineering, Computing and\\ Mathematical Sciences}\\%
				\vspace{\stretch{1}}
			}%
			\if@assoc
			%%% Auxiliary association %%%
			{%
			\fontsize{12}{14.4} 
			\textbf{\@auxassoc}\\%
			\vspace{\stretch{1}}
			}%
			\fi
			%%% Thesis title %%%
			{%
				\fontsize{14}{16.4} 
				\textbf{\@thesistitle}\\%
				\vspace{\stretch{1.25}}
			}%
			%%% Author name %%%
			{%
				\fontsize{12}{14.4}
				\textbf{\@authorname}\\%
				\vspace{\stretch{0.2}}
			}%
			%%% Declaration %%%
			\vspace{\stretch{1.25}}
			{%
				\setstretch{1}
				\fontsize{10}{12} 
				\textbf{This thesis is presented for the Degree of}\\%
				\textbf{Doctor of Philosophy}\\%
				\textbf{of}\\%
				\textbf{Curtin University of Technology}\\%
				\vspace{\stretch{0.7}}
			}%
			%%% Submission date %%%
			{%
				\fontsize{12}{14.4} 
				\textbf{\@graduationdate}
			}%
		\end{center}
	\end{titlepage}
}%

\newcommand{\plagiarismstatement}{%
	\setcounter{secnumdepth}{3} % allow lower level numbering
	\setcounter{tocdepth}{4} % allow lower level numbering in table of contents
	\setcounter{page}{2}
	\chapter*{}
	%%% Plagiarism statement %%%
	\begin{center}
		\vspace*{\fill}
		To the best of my knowledge and belief this thesis contains no material previously published by any other person except where due acknowledgement has been made. This thesis contains no material which has been accepted for the award of any other degree or diploma in any university.\\%
		\vspace{\stretch{0.2}}
		\begin{flushright}
			\@authorname
		\end{flushright}
		\vspace*{\fill}
	\end{center}	
}%

\newcommand{\showthesisquote}{%
	\newif\if@quotepresent
	\ifx\@quote\@empty\@quotepresentfalse\else\ifx\@quoteauthor\@empty\@quotepresentfalse\else\@quotepresenttrue\fi\fi
	\if@quotepresent
		\chapter*{}
		%%% Thesis Quote %%%
		\begin{center}    
			\vspace*{\fill}
			\begin{flushleft}
				``\textit{\@quote}''
			\end{flushleft}
			\begin{flushright}
				--- \@quoteauthor
			\end{flushright}
			\vspace*{\fill}
		\end{center}
	\fi
}%

\newcommand{\includetoc}{%
	\addtolength{\baselineskip}{3pt}
	\tableofcontents
	\addtolength{\baselineskip}{-3pt}
	\vspace*{\fill}
}%

\newcommand{\includelof}{%
	\addtolength{\baselineskip}{1pt}
	\listoffigures
	\addtolength{\baselineskip}{-1pt}
	\vspace*{\fill}
	\newpage{\pagestyle{empty}\cleardoublepage}
}%

\def\arabicnumbering#1{\pagebreak\pagenumbering{arabic}\setcounter{page}{#1}}
\def\alphnumbering#1{\pagenumbering{Alph}\setcounter{page}{#1}}
\def\addapptotoc{\addtocontents{toc}{\protect\setcounter{tocdepth}{1}}}

\newenvironment{customappendix}[1]{
	\newpage{\pagestyle{empty}\cleardoublepage}
	#1
	\begin{appendices}
	\addapptotoc % add the appendices to the table of contents 
}{
	\end{appendices}
}

%% define acknowledgement statement as part of references %%
\def\ackstatement{\vspace*{\fill}%
	Every reasonable effort has been made to acknowledge the owners of copyright material.
	I would be pleased to hear from any copyright owner who has been omitted or incorrectly acknowledged.}%

\endinput
