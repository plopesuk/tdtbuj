
PROJECT(thesis NONE)
INCLUDE(UseLATEX.cmake)
cmake_minimum_required(VERSION 2.4.8)
if(COMMAND cmake_policy)
  cmake_policy(VERSION 2.6)
endif(COMMAND cmake_policy)


SET(LATEX_OUTPUT_PATH bin)
SET(CMAKE_CURRENT_SOURCE_DIR src)

ADD_LATEX_DOCUMENT(main.tex
  INPUTS chap1.tex chap2.tex chap3.tex 
  appendices.tex SphericalHarmonics.tex
  backmatter.tex frontmatter.tex glossary.tex
  xkeyval.tex xkvtxhdr.tex  pst-xkey.tex  keyval.tex
  BIBFILES mythesis.bib
  IMAGE_DIRS figures
  USE_INDEX  
  USE_GLOSSARIES
  DEFAULT_PDF
  )
