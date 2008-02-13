\ProvidesClass{thesis}
%% For normal draft builds
%\LoadClass[hyperpdf,nobind,draft,sftitles,twoside]{hepthesis}
%% For short draft builds (breaks citations by necessity)
%\LoadClass[hyperpdf,nobind,draft,hidefrontback]{hepthesis}
%% For soft-bound version
%\LoadClass[hyperpdf,bindnopdf]{hepthesis}
%% For hard-bound version (must be one-sided)
\LoadClass[bind,hyperpdf,a4paper,booktabs,index]{hepthesis}

%% -------------------------------------
%% Standard packages
%% -------------------------------------
\RequirePackage{cite}
\RequirePackage{mathrsfs} % script font
\RequirePackage{afterpage}
\RequirePackage{glossaries}

%% Using Babel allows other languages to be used and mixed-in easily
\RequirePackage[english]{babel}
\selectlanguage{english}

%% Quotes in different languages are different.
%% Using \enquote{bla} automatically uses the correct quotes
%% according to the Babel language that is set
\RequirePackage[babel]{csquotes}
\RequirePackage{subfig}