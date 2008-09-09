
PROJECT(thesis NONE)
INCLUDE(UseLATEX.cmake)
cmake_minimum_required(VERSION 2.4.8)
if(COMMAND cmake_policy)
  cmake_policy(VERSION 2.6)
endif(COMMAND cmake_policy)


SET(LATEX_OUTPUT_PATH bin)
SET(CMAKE_CURRENT_SOURCE_DIR src)

ADD_LATEX_DOCUMENT(thesis.tex
  INPUTS
  abstract.tex glossary.tex acknowledgement.tex
  introduction.tex chapter1.tex chapter2.tex chapter3.tex chapter4.tex
  conclusions.tex
  appendixA.tex appendixB.tex appendixC.tex appendixD.tex
  bfield.tex
  xkeyval.tex xkvtxhdr.tex  pst-xkey.tex  keyval.tex
  BIBFILES mythesis.bib
  IMAGE_DIRS figures
  USE_GLOSSARY
  USE_INDEX
  DEFAULT_PDF
  )
