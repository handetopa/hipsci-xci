extract_text_before_first_dot <- function(name) {
  
  name=gsub("\\..*", "", name)
  return(name)
  
  
}