# Specify a minimum required version
cmake_minimum_required(VERSION 2.8)

# Specify the tex file name
set(TEX_NAME hartree-fock)


# Define a LaTeX project, and specify where to find the UseLATEX commands
project(hartree-fock NONE)			# The NONE option has to be included, since we are not using the langue c, cxx, or Fortran
include(cmake/UseLATEX.cmake)		# I like anti-cluttering by moving the UseLATEX.cmake file into a separate directory

# BUILD: Add a file to the specified project
add_latex_document( ${TEX_NAME}.tex
                    BIBFILES ../../../../Documents/Archief/Bibliotheek/research_bib.bib)


# POST-BUILD: move the created .pdf-file to the top directory (copy and delete)
add_custom_command(TARGET "${TEX_NAME}_pdf" POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy "${TEX_NAME}.pdf" ..
COMMAND ${CMAKE_COMMAND} -E remove "${TEX_NAME}.pdf"
COMMENT "Moving PDF-file to top directory.")
