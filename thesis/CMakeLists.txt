
PROJECT(thesis NONE)
INCLUDE(UseLATEX.cmake)


SET(LATEX_OUTPUT_PATH bin)
SET(CMAKE_CURRENT_SOURCE_DIR src)

ADD_LATEX_DOCUMENT(main.tex
  INPUTS chap1.tex chap2.tex chap3.tex appendices.tex backmatter.tex frontmatter.tex
  BIBFILES mythesis.bib
  IMAGE_DIRS figures
  USE_INDEX  
  USE_GLOSSARY
  DEFAULT_PDF
  )
