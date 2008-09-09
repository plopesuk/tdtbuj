% thesis.cls -- A LaTeX PhD Thesis class by Alin M Elena
% created on 27th of March 2008
%modified on 28th of August 2008
\ProvidesClass{aprilthesis}[2008/08/28]
\typeout{Based on standard book class, with customisations}
 \usepackage[pdftex]{graphicx}
 \usepackage{setspace}
 \usepackage{extramath}
 \usepackage{fancyhdr} % funky headers
 \usepackage[T1]{fontenc}
 \usepackage{times}
% 
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
%
\ProcessOptions\relax
\LoadClass{book}
% load calligraphic font
\DeclareFontFamily{OT1}{pzc}{}
\DeclareFontShape{OT1}{pzc}{m}{it}{<-> s * [1.700] pzcmi7t}{}
\DeclareMathAlphabet{\mathscr}{OT1}{pzc}{m}{it}
%
%------------------------------------------------------------------
% Use just over 1.5x line spacing (assuming 12pt font)
 \onehalfspacing
%------------------------------------------------------------------
% Make pages a bit larger, given we are using 1.5 spacing
%%% all these things were computed, well approximated
%% and guessed for a4paper,12pt options of the book package
%horizontal part
   \setlength{\paperwidth}{210mm}   
   \setlength{\hoffset}{0.0mm}   
   \setlength{\oddsidemargin}{14.6mm}
   \setlength{\evensidemargin}{0.9mm}
   \setlength{\textwidth}{145.0mm}
%vertical part
   \setlength{\paperheight}{297mm}
   \setlength{\voffset}{-20.4mm}
   \setlength{\headheight}{10.0mm}
   \setlength{\headsep}{5.0mm}
   \setlength{\textheight}{242.0mm}
   \setlength{\footskip}{10.0mm}
%------------------------------------------------------------------
% Increase float fraction and discourage float-only pages
 \renewcommand\topfraction{.9}
 \renewcommand\bottomfraction{.5}
 \renewcommand\dbltopfraction{.9}
 \renewcommand\textfraction{.1}
 \renewcommand\floatpagefraction{0.8}
 \renewcommand\dblfloatpagefraction{0.8}
%------------------------------------------------------------------
% Title page stuff
 \newcommand*{\affiliation}[1]{\gdef\@affiliation{#1}}
 \newcommand*{\acadtitle}[1]{\gdef\@acadtitle{#1}}
 \newcommand*{\tdesc}[1]{\gdef\@tdesc{#1}}
 \newlength{\logo@width}
 \newcommand*{\logo}[2][0pt]{\gdef\@logo{#2}\setlength{\logo@width}{#1}}

 \renewcommand\maketitle{\begin{titlepage}%
     \let\footnotesize\small
     \let\footnoterule\relax
     \let \footnote \thanks
     \null\vfil
% %    \vskip 60\p@
     \begin{center}%
      {\huge \sc \@affiliation \par}
      \ifx\@logo\@undefined
         \vskip 10em
         \else%
         \ifdim\logo@width>0\p@\relax
         \begin{figure*}[!ht]
            \centering
             \includegraphics[width=\logo@width]{\@logo}%
          \end{figure*}
         \else
          \begin{figure*}[!ht]
             \centering
              \includegraphics{\@logo}%
          \end{figure*}
         \fi
       \fi
       {\Huge\bf \@title \par}%
       \vskip 5em%
       {\Large\bf \sc
       \lineskip .75em%
       \begin{tabular}[t]{c}%
         \@author
       \end{tabular}\par}%
       {\Large \@acadtitle}
       \vskip 5em%
       {\large \@tdesc\par}%
       \vskip 3em%
       {\large \@date \par}%       % Set date in \large size.

     \end{center}\par
     \vfil%
   \null
   \end{titlepage}%
   \setcounter{footnote}{0}%
   \global\let\thanks\relax
   \global\let\maketitle\relax
   \global\let\@thanks\@empty
   \global\let\@author\@empty
   \global\let\@date\@empty
   \global\let\@title\@empty
   \global\let\title\relax
   \global\let\author\relax
   \global\let\date\relax
   %\global\let\and\relax
 }
%
 \newcommand{\makededication}[2]{
 \cleardoublepage
 \thispagestyle{empty}
  \null\vfil%
   \begin{center}   
     \textbf{#1}\\
       {#2}    
    \end{center}
   \null\vfil% 
 }

 \newcommand{\mainquote}[2]{%
   \cleardoublepage
   \thispagestyle{empty}
   \null\vfil%
       \begin{center}
       \begin{spacing}{1.25}
  	     #1\\
 	\quad\quad \emph{#2}
      \vspace{12pt}
      \end{spacing}
       \end{center}%
    \null\vfil% 
 }



% %------------------------------------------------------------------
% % Chapter headings from:
% % http://zoonek.free.fr/LaTeX/LaTeX_samples_chapter/0.html
% % Amended to know about not numbering outside the mainmatter
%

\def\thickhrulefill{\leavevmode \leaders \hrule height 1ex \hfill
 \kern \z@}
%
 \def\@makechapterhead#1{%
   %\vspace*{50\p@}%
   \vspace*{10\p@}%
   {\parindent \z@ \raggedright \reset@font
     \if@mainmatter
        $\mathscr{{\@chapapp{}\; \thechapter}}$%
     \fi
         \par\nobreak
         \vspace*{1\p@}%
         \interlinepenalty\@M
         \hrule
         \vspace*{10\p@}%
         {\raggedleft \Huge \bfseries #1\par\nobreak}
         \par
         %\vspace*{10\p@}%
         %\hrule
     \vskip 40\p@
     %\vskip 100\p@
   }}

 \def\@makeschapterhead#1{%
   %\vspace*{50\p@}%
   \vspace*{10\p@}%
   {\parindent \z@ \raggedleft \reset@font
         \vspace*{1\p@}%
         \interlinepenalty\@M
         \hrule
         \vspace*{10\p@}%
         \Huge \bfseries #1\par\nobreak
         \par

     \vskip 40\p@
     %\vskip 100\p@
   }}
\pagestyle{fancy}
%customize headers
\renewcommand{\chaptermark}[1]{\markboth{}{#1}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}
\fancyhf{}
\fancyhead[LE,RO]{\bfseries \thepage}
\fancyhead[LO]{\nouppercase{\emph{\rightmark}}}
\fancyhead[RE]{\nouppercase{\emph{\rightmark}}}
%\lhead{\nouppercase{\emph{\tmark}}}
%\rhead{\nouppercase{\emph{\leftmark}}}
\cfoot{}
% %------------------------------------------------------------------
% % Fix table of contents to use chapter{} not chapter* so that running
% % head is set correctly. It will still appear un-numbered as it
% % is used only in 'frontmatter'
%
 \renewcommand\tableofcontents{%
  \if@twocolumn
   \@restonecoltrue\onecolumn \else \@restonecolfalse \fi
   \chapter{\contentsname}
    \markboth{}{\contentsname}
   \@starttoc{toc}%
   \if@restonecol\twocolumn\fi
 }
 \renewcommand\listoffigures{%
  \if@twocolumn
   \@restonecoltrue\onecolumn \else \@restonecolfalse \fi
   \chapter{\listfigurename}
    \markboth{}{\listfigurename}
   \@starttoc{lof}%
   \if@restonecol\twocolumn\fi
 }
 \renewcommand\listoftables{%
  \if@twocolumn
   \@restonecoltrue\onecolumn \else \@restonecolfalse \fi
   \chapter{\listtablename}
    \markboth{}{\listtablename}
   \@starttoc{lot}%
   \if@restonecol\twocolumn\fi
 }


 %
 %------------------------------------------------------------------
 % chapter quote.
 % it gets right aligned
 % two parameters first the quote next one the author of the quote
 \newcommand{\chapterquote}[2]{
 \begin{spacing}{1.0}
 \begin{quote}
 \begin{flushright}
  #1\\
 \hfill \emph{#2}
 \end{flushright}
 \end{quote}
 \end{spacing}
 \vskip 5pt
 }

     %%% Clear Header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Clear Header Style on the Last Empty Odd pages
     \makeatletter
     \def\cleardoublepage{\clearpage\if@twoside \ifodd\c@page\else%
         \hbox{}%
         \thispagestyle{empty}%              % Empty header styles
         \newpage%
         \if@twocolumn\hbox{}\newpage\fi\fi\fi}
     \makeatother

 %Fix hyperref labels in PDFTeX
% 	% \@ifpackageloaded{hyperref}{%
% 	%   \pdfstringdefDisableCommands{%
% 	%       \def\pi{pi}%
% 	%       \def\gamma{gamma}%
% 	%       \def\alpha{alpha}%
% 	%   }%
% 	% }{}