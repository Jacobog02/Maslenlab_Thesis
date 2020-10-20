## Jacob Gutierrez gutierja@ohsu.edu
## 2/4/20 
## Pretty Table Printing Helper Functions

## Libraries
require(dplyr)
require(knitr)
require(kableExtra)

## This function takes a df and creates pretty html output for knitting. 
jg_pretty_print <- function(df, cap = NA){
  
  df %>%  kable(caption = cap) %>% 
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive") , full_width = F, position = "left") %>%  
    scroll_box(width = "100%", height = "300px")
  
} 