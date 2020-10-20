## ----information,eval = F------------------------------------------------
## # Jacob Gutierrez gutierja@ohsu.edu 7/29/19
## # This file is a wrapper to dynamically render the .Rmd File into the target working directory. 
## Overall, this script catches bash variables to direct the proper input, output and parameters. 
## See this post for inspiration: https://stackoverflow.com/questions/32479130/passing-parameters-to-r-markdown

## Catch Arguements
args <- commandArgs(trailingOnly = TRUE)
#args <- commandArgs()

print(args)
if (length(args) == 3) {
  expname <- args[1]
  wrkdir <- args[2]
  script_dir <- args[3]
  set_title <- sprintf('%s: Probe Validation Visualization',expname)
  #notation= 'chr'
  outfile <- sprintf('%s/%s_visuals.html',wrkdir,expname)
  print(outfile)
  #print('HERE')
  #print(set_title)
}  else{
  stop("1) Experiment Name, 2) Working Directory, and 3) Script Directory required to generate html output", call.=FALSE)
}


## Build Parameter object
parameters <-  list(expname = expname,
                    wrkdir = wrkdir,
                    set_title = set_title)


## Render Document
template <- "updated_probe_visuals.Rmd"
render_path <- sprintf("%s/%s",script_dir,template)
#print(render_path)

rmarkdown::render(input=render_path,"html_document", output_file = outfile, params = parameters)

