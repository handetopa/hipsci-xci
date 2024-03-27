# -----------------------------------------------------------
# Script Name: extract_text_before_first_dot.R
# Author: Hande Topa
# Date: 2024-03-20
# -----------------------------------------------------------

extract_text_before_first_dot <- function(name) {
  name=gsub("\\..*", "", name)
  return(name)
}