## JG SVA functions 
## 1/20/20 gutierja@ohsu.edu
## JG 10/20/20: I did a major refactoring to make PCA and PCPR2 functionalized so you only need to store plotting data and not the large plots. 
## JG EDA function


### Libraries
require(sva) # For vanilla SVA
#require(isva) # nm.sv helper function
#require(SmartSVA) # For smartSVA function
require(broom) # for tidy function
require(reshape) # for correlation function data wrangling
require(gplots) # textplot
require(corrplot) # Correlation Plots
require(RColorBrewer) ## JG_custom_mds
require(ggfortify) ## PCA/MSA biplots

## Correct SVA function
svaBatchCor <- function(mat, sva_out){
  ## This function removes all of the surrogate variables from a given matrix 
  ## This action is not reccomended for any downstream analysis
  ## THIS IS A VISUALIZATION TECHNIQUE ONLY >:|
  
  # Do something to the matrix
  Y <- t(mat)
  
  # Grab surrogates
  W <- sva_out$sv
  # Solve linear algebra based on the link above
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  # Remove the correction from the original matrix (i think...) and store it back in object.
  sva_out$corrected <- t(Y - W %*% alpha)
  
  
  # Return to the overworld 
  return(sva_out)
  
}

## JG updated 2/1/2020 
## Adding smartSVA to be extra sure for cell type composition adjusmtent. 
## Changing functition to accept new argument algorithm Leek, smartSVA


## SVA Asssessment Function
## This was done because ~50% of the varaibility of the data cannot be explained by the batch variables. ## THIS WAS ABOUT THE MICRO ARRAY DATA
## algorith choices are: "Leek" and "smart"
## surrogate variable estimation choices are: "be", "leek", "RMT"
sva_assessment <- function(mat,pd,batch_factors,main_factor,bio_covariates,title,discover = FALSE, algorithm = "Leek", sv_est_alg = "be"){
  
  ## Remove NA's matrix
  mat <- remove_invariance(mat)
  
  # I should make it check for places with no variablity here but I do it outside. 
  
  # Make model matries! 
  ## If in discover mode then only include the main factor in the model
  ## If in adjust mode then include the main facotr only in mod and not in the rest of the model. 
  if(discover){
    ## Change the title
    title <- sprintf('Discover: ',title)
    ## Of interest
    interest <- as.formula(sprintf('~ %s', main_factor))
    mod <- model.matrix(interest,data = pd)
    null <- as.formula('~ 1')
    mod0 <- model.matrix(null,pd)
  } else{
    ## Not the main factor included ONLY adjusting for biological covariates. 
    adjust <- bio_covariates[ !(bio_covariates %in% main_factor)]
    
    
    interest <- as.formula(sprintf('~  %s', paste(c(main_factor,adjust),collapse = ' + ')))
    mod <- model.matrix(interest,data = pd)
    
    
    null <- as.formula(sprintf('~ %s', paste(adjust,collapse = ' + ')))
    mod0 <- model.matrix(null,data = pd)
    
    
  }
  
  
  ## Compute number of surrogates
  trad_sva = c("be", "leek")
  ## use the two surrogate variable estimation algorithms
  if (class(sv_est_alg) == "numeric"){
    #print("Detected manually set number of SVs")
    warning("Detected manually set number of SVs to save computational time.\nPlease ensure this estimate is accurate.")
    n.sv <- sv_est_alg
    
      
  }else if(sv_est_alg %in% trad_sva){
    ## Be method is permutation based
    ## leek method is asymptotic approach
    n.sv <- num.sv(mat, mod, method = sv_est_alg, seed = 2020)
  }else if (sv_est_alg == "RMT"){
    ## Random Matrix Theory Estimation : Better for many samples. 
    df <- data.frame(pred=pd[,main_factor])
    ## Determine the number of SVs
    mat.r <- t(resid(lm(t(mat) ~ pred, data=df)))
    
    ## Add one extra dimension to compensate potential loss of 1 degree of freedom
    ##  in confounded scenarios (very important)
    n.sv <- EstDimRMT(mat.r, FALSE)$dim + 1
  }else{
    stop("Undefined Surrogate Estimation Algorithm Selected")
  }
  print(sprintf("Number of significant surrogate variables is: %i",n.sv ))
  
  
  
  ## Check algorithm 
  ## Leek algorithm is the original SVA implementation CITE HERE
  if(algorithm == "Leek"){
    ## Run SVA
    sva_mod <- sva(mat,mod,mod0,n.sv = n.sv)
  }
  ## smartSVA algorithm was made to optimize cell type composition adjustment in EWAS
  else if(algorithm == "smart"){

    
    ## FAILS BC INVARIATE GENES. << LIES I guess RMT overestimated the number, we cant get effecient sv extraction bc limited samples.  
    sva_mod <- smartsva.cpp(mat, mod, mod0,n.sv= n.sv)
  }
  ## Raise Error
  else{
    stop("Undefined SVA Algorithm Selected")
  }
  
  ## Linear Model the technical 
  n <- sva_mod$n.sv
  sv <- sva_mod$sv
  
  # Loop through each SV
  for (i in seq(n)) {
    # Make a buffer df that can be added into a simple model col1 = SV col2-coln = Batch factors
    linear_frame <- cbind(SV = sv[,i],pd[,batch_factors])
    # make the formula model (I make SV always the single surrogate variable to simply the formula expression)
    mod_l <- as.formula(sprintf('SV ~ %s', paste(batch_factors,collapse = ' + ')))
    # Run the model 
    linear_mod <- lm(mod_l, data = linear_frame)
    # Use tidy to make into structured df
    lin_results <- tidy(linear_mod)
    # Label surrogate variable based on it number
    lin_results$var_num <- i
    
    # Start table if its first otherwise append table. 
    if (i == 1) {
      out_tbl <- lin_results
    } else {
      out_tbl <- rbind(out_tbl, lin_results)
    }
  }
  
  
  ## PLOTTING!
  # Now I need to make a plot of the heatmap using a similar design as below. 
  
  ## Convert pvalues into stars
  plotting <- out_tbl %>% mutate(stars = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")))
  ## Make the surrogates factors
  plotting <- plotting %>% mutate(var_num = as.factor(var_num))
  ## Make heat map plotting the x= SV y=batch effect text + color = signficance (stars)
  heatmap <- ggplot(plotting,aes(var_num,term,fill=stars)) + geom_tile() + theme_bw() + 
    geom_text(aes(label=stars)) +
    ggtitle(sprintf('Heat Plot of Batch Effect Correlations to Surrogates (Linear)\n%s', title))
  #print(heatmap)
  
  # Return the sig components, and the full model
  sva_mod$full_lm <- out_tbl
  sva_mod$sig_lm <- out_tbl %>% filter(p.value < .05 )
  sva_mod$batch_factors <- batch_factors
  sva_mod$heatplot <- heatmap
  
  ## Adding the title into the object to be used downstream!!! 1/24/20
  sva_mod$title <- title
  
  ## Put model info in 
  sva_mod$model <- mod
  sva_mod$null_model <- mod0
  
  # Correct data FOR VISUALIZATIONS ONLY!!! >:|
  sva_mod <- svaBatchCor(mat,sva_mod)
  
  ## PLOT THE SV BAR CHART
  to_plot <- as.data.frame(sva_mod$sv) 
  colnames(to_plot) <- paste("SV", seq(ncol(to_plot)), sep = '_')
  
  to_plot <- pivot_longer(to_plot, everything(), names_to = 'SV', values_to = 'value') %>% 
    mutate(SV = as.factor(SV))
  
  bar <- to_plot %>% ggplot(aes(x=SV,y=value,fill=SV)) + geom_boxplot() + ggtitle(sprintf('Boxplot of Surrogates Values \n%s', title))
  sva_mod$barplot <- bar
  
  return(sva_mod)
}

## PC-PR2 Function
## ADAPTED FROM the R code available in Fages et al.’s supplementary material
## Changed to accept a PCA object 
PC_PR2 <- function(pca, pd, batch_factors,title,thrs= .8){
  
  ## Extract the phenodata
  ## Subsetting to only of interest
  Z_InteresFactors <- pd[,batch_factors]
  Z_RowN <- nrow(Z_InteresFactors)
  Z_ColN <- ncol(Z_InteresFactors)
  ColNames <- names(Z_InteresFactors)
  
  # So we need a dataframe of just the batch data for each sample in addition to a list of the factors to analyze.
  #mat <-  na.omit(mat)
  #mat <- t(mat)
  # CALCUATE COVARIANCE MATRIX I FOUND THIS IN THIS PACKAGE
  # https://cran.r-project.org/web/packages/coop/vignettes/coop.pdf
  #cov_mat <- 1/(NROW(mat) -1) * crossprod(mat)
  #cov_mat <- coop::tcovar(mat) ## THESE DONT WORK :)))))
  
  ## NOw trying big cov function 
  #cov_mat <- bigcor(x = t(mat),fun = 'cor')
  
  
  #mat2 <- mat %*% t(mat) # THis estimates the covariance matrix! idk why they did it this way it runs out of memory
  #cov_mat <- cov(t(mat)) # This is memory effecient way to do the computation above. 
  
  
  ## I will use PCA function to do the deal. 
  #pca = prcomp(t(mat))
  #pca_sum <- summary(pca)
  
  
  
  # Compute eigenrepresentation 
  #eigenData <- eigen(cov_mat) # Cant compute full cov matrix
  #eigenValues <- eigenData$values
  
  # Using prcomp
  eigenValues <- (pca$sdev)^2
  
  
  ev_n <- length(eigenValues)   
  
  #eigenVectorMatrix <- eigenData$vectors  #CANNOT COMPUTE full covmat
  #eigenVectorMatrix <- pca$rotation
  eigenVectorMatrix <- pca$x
  eigenValuesSum <- sum(eigenValues)
  percentPCs <- eigenValues/eigenValuesSum
  
  my_counter_2 = 0
  my_sum_2 = 1
  
  for(i in ev_n:1){
    my_sum_2 = my_sum_2 - percentPCs[i]
    if((my_sum_2) <= thrs){
      my_counter_2 = my_counter_2 + 1
    }
    
  }
  
  if(my_counter_2 < 3){
    pc_n = 3
  }else{
    pc_n = my_counter_2
  }
  
  pc_data_matrix <- matrix(data=0, nrow = (Z_RowN * pc_n) ,ncol = 1 )
  
  my_counter = 0
  
  for(i in 1:pc_n){
    for(j in 1:Z_RowN){
      my_counter <- my_counter + 1
      pc_data_matrix[my_counter,1] = eigenVectorMatrix[j,i]
    }
  } # End of pc_data_matrix form
  
  AAA <- Z_InteresFactors[rep(1:Z_RowN,pc_n),]
  Data <- cbind( data_matrix = pc_data_matrix,AAA)
  # Be sure all the facotrs and numeric data is proper at this point.
  
  DataCol <- ncol(Data)
  typeIIIMatrix <- matrix(data=0, nrow = pc_n ,ncol = DataCol)
  ST_ResidualR2 <- matrix(data=0, nrow = pc_n, ncol = 2)
  
  
  ## First column is the data and we model against it
  ## Now making a formula to run linear model.
  data_names <- colnames(Data)
  explain <- data_names[1]
  other <- paste(data_names[2:length(data_names)], collapse = ' + ')
  form <- paste(explain,other,sep= ' ~ ')
  form <- as.formula(form)
  
  
  
  for (i in 1:pc_n){
    
    y = (((i-1)*Z_RowN) + 1 )
    TotSumSq <- var(Data[y:(((i - 1)*Z_RowN)+ Z_RowN),1]) * (Z_RowN-1)
    mod_dat <- Data[y:(((i-1)*Z_RowN) + Z_RowN),]
    Model <- lm(form, mod_dat)
    
    AnalysisVariance <- anova(Model)
    SumSq <- AnalysisVariance[2]
    #Residuals <- SumSq[DataCol,]
    ## JG 10/15/20 making explicit Residuals calling. 
    res_idx <- rownames(SumSq) == "Residuals"
    Residuals <- SumSq[res_idx,]
    RR <- Residuals/TotSumSq
    R2 <- 1 - RR
    ST_ResidualR2[i,] <- c(R2,RR)
    ST_ResidualR2_Names <- c('ST_R2','ST_Residuals')
    colnames(ST_ResidualR2) <- ST_ResidualR2_Names
    for (j in 1:DataCol){
      typeIIIMatrix[i,j] <- as.numeric(SumSq[j,1])
      typeIIIMatrixNames <- c(ColNames,'SumSqResiduals')
      colnames(typeIIIMatrix) <- typeIIIMatrixNames
    }  
    
  }# End of the typeIIIMatrix computation
  
  partialR2Matrix <- matrix(data=0, nrow = pc_n, ncol = DataCol -1)
  
  for (i in 1:pc_n){
    for (j in 1:DataCol-1){
      partialR2Matrix[i,j] <- typeIIIMatrix[i,j] / (typeIIIMatrix[i,DataCol] + typeIIIMatrix[i,j])
    }
  }
  
  partialR2MatrixWtProp <- matrix(data=0, nrow = pc_n, ncol = DataCol)
  
  for (i in 1:pc_n){
    weight = eigenValues[i]/sum(eigenValues[1:pc_n])
    for (j in 1:DataCol-1){
      partialR2MatrixWtProp[i,j] <- partialR2Matrix[i,j] * weight
      partialR2MatrixWtProp[i,DataCol] = ST_ResidualR2[i,1] * weight
    }
  }
  
  pR2Sums <- colSums(partialR2MatrixWtProp)*100
  plotnames <- c(ColNames,'R2')
  
  bp <- barplot(pR2Sums, xlab = 'Factors', ylab='Weighted Rpartial2', ylim =c(0,100.1),col = c('red'), las =2,main = title)
  axis(1,at = bp, labels = plotnames, xlab = 'Factors',cex.axis =0.5,las =2)
  values <- pR2Sums
  new_values = round(values,3)
  text(bp,pR2Sums,labels=new_values,pos=3,cex=0.8)
  plot <- recordPlot()
  
  out_mat <- partialR2MatrixWtProp
  colnames(out_mat) <- plotnames
  rownames(out_mat) <- paste("PC",seq(pc_n), sep='')
  
  ## Making Heat Map of R2 Values
  plotting <- out_mat %>% as.data.frame()%>% rownames_to_column("PC") %>% pivot_longer(-PC,names_to = "variable", values_to = "pr2")
  ## Adding Clean labels
  plotting <- plotting %>% mutate(labs = sprintf("%g%%",round(pr2*100, 1)))
  ## Removing R2 since it becomes confusing and isnt a variable. 
  plotting <- plotting %>% filter(variable != "R2")

  
  ## JG 10/20/20: Adding prop_var data to make a pretty X lab. 
  prop_var <- summary(pca)$importance[2,1:pc_n]
  buffer <- reshape::melt(prop_var)
  buffer$PC <- rownames(buffer) ## Manually setting this in THIS case is ok bc it is within pca function
  buffer <- buffer %>% as.data.frame() %>% dplyr::rename(perc = value)
  plotting <- left_join(plotting,buffer) ## add to plotting
  
  plotting <- plotting %>% mutate(PC_label = sprintf('%s\n(%g%%)', PC, round(perc *100, 0) ))
  ## JG 10/15/20: Adding levels ot newPC so the x axis is stable
  plotting <- plotting %>% mutate(PC_label = factor(PC_label, levels = unique(PC_label)))
  ## Please see functional approach to this plot below.   
#   heatmap <- ggplot(plotting,aes(PC,variable,fill=pr2)) + geom_tile() + 
#     scale_fill_gradient(low = 'white', high="red", breaks= c(0, 1,5, 10,20,30,40))+ 
#     ggtitle(sprintf('heatplot of partial R2 for each batch variable
# for all %i PCs', pc_n)) + geom_text(aes(label=labs)) + theme(legend.position = "left")
  
  #print(heatmap)
  
  
  out <- list()
  #out$plot <- plot
  out$r2wtprop <- out_mat
  #out$pcheatmap <- heatmap
  out$plotting <- plotting
  out$pcpr2_heatmap <- pcpr2_heatmap
  return(out)
  
}

## JG 10/20/20: PCPR Heatmap Function (STOLEN from the correlation heatmap below)
pcpr2_heatmap <- function(plot_df, X = "PC_label", Y = "variable", fill_col = "pr2", lab_col = "labs", title = ""){
 
  ## MUST USE QUASIQUOTATION TO WORK WITH VARIABLES in ggplot :)
  #https://rlang.r-lib.org/reference/nse-force.html
  
  ## Color system is for percent variance explained. 
  heatmap <- ggplot(plot_df,aes(!!sym(X),!!sym(Y),fill=!!sym(fill_col))) + geom_tile() + 
    scale_fill_gradient(low = 'white', high="red", breaks= c(0, 1,5, 10,20,30,40))+ 
    ggtitle(sprintf('Heatplot of Partial R2 for each Variable')) + geom_text(aes(label=!!sym(lab_col))) + theme(legend.position = "left")
  
  #print(heatmap)
  return(heatmap)
  
} # End of PCA assessment plotting function






## JG 10/15/20: This function does not accept tibbles they have specialized behavior that makes this script sad.
#assess_pca <- function(mat, pd,batch_factors,title){
## JG 10/20/20: I am modifying this script to offically retire the batch_factors argument. Instead it will be a complete pheno matrics
## Pheno (or pd) is a df with a column corresponding to the batch_factors and can be derived from the pd matrix. 
assess_pca <- function(mat, pd,title, total = NA, heatmap = F){
  
  ## NA's mess up PCA so I will remove them here
  mat <- na.omit(mat)
  
  if (is.na(total)){
    total <- ncol(mat)
  } else if (total < 3){
    total <-  3
  }
  #total <- 10
  one.pca <- prcomp(t(mat), rank. = total)
  one.sum <- summary(one.pca)
  prop_var <- one.sum$importance[2,1:total] # second row is proportion explained
  
  ## JG 10/20/20: Perform correlation and return tidy data. (does not have importance data)
  plotting <- pca_tidycors(one.pca$x, pd,x = "PC",y="Batch")
  
  ## JG 10/20/20: Adding prop_var data to make a pretty X lab. 
  buffer <- reshape::melt(prop_var)
  buffer$PC <- rownames(buffer) ## Manually setting this in THIS case is ok bc it is within pca function
  buffer <- buffer %>% as.data.frame() %>% dplyr::rename(perc = value)
  plotting <- left_join(plotting,buffer) ## add to plotting
  
  plotting <- plotting %>% mutate(PC_label = sprintf('%s\n(%g%%)', PC, round(perc *100, 0) ))
  ## JG 10/15/20: Adding levels ot newPC so the x axis is stable
  plotting <- plotting %>% mutate(PC_label = factor(PC_label, levels = unique(PC_label)))
  
  ## To plot 
  #jg_plot_pca(one.pca) ## This can be done with any  PCA. 
  #cor_heatmap(plot_df = plotting,X = "PC_label", Y = "Batch", fill_col = "cor", lab_col = "stars", title = "")
  
  
  out <- list()
  out$pca <- one.pca
  #out$heatmap <- heatmap
  ## JG 10/15/20: Adding the plotting data
  out$plot_df <- plotting
  ## JG 10/20/20: Adding the text plot
  #out$txt_plot <- txt_plot
  out$jg_plot_pca <- jg_plot_pca
  ## JG 10/22/20: Adding biplot function to output. 
  out$jg_autoplot <- jg_autoplot
  ## JG 10/23/20: ADding cor_heatmap 
  out$cor_heatmap <- cor_heatmap
  
  ## CONSIDER CAPTURING THE PLOTS TO BE USED IN THE MARKDOWN... 
  return(out)
} # End of PCA assessment plotting function
  
## JG 10/20/20: Plot PCA Summary
## Takes in a PCA object performs summary and then returns the plot.
jg_plot_pca <- function(one.pca, title = ""){
  ## Compute Summary
  one.sum <- summary(one.pca)
  prop_var <- one.sum$importance[2,] # second row is proportion explained
  

  # Skee plot
  plot(prop_var, type='b')
  title(main=sprintf('Skree Plot: Amount of Variance Explained in each PC\n%s',title))
  
  # Cum Sum plot
  plot(cumsum(prop_var), type = 'b')
  title(main=sprintf('Cumulative Sum of Variance Explained\n%s',title))
  
  # Text Plot of the summary 
  txt_plot <- textplot(capture.output(one.sum),valign='top', mar =c(5,4.5,4,2))
  #return(txt_plot)
  print(txt_plot)
}## End of jg_plot_pca


## JG 10/20/20: Pairwise Correlation df
## I need a function that takes in two dfs and creates a pairwise correlations which is returned in tidy format. 
## x corresponds to the df name of df1 and y is df2. 
pca_tidycors <- function(df1, df2,x = "PC",y="Batch"){
  ## In this case df1 is the PCA and df2 is pheno. 
  ## All columns in df1 will be correlated to the df2 matrix. 
  df1 <- as_numeric_mat(df1) ## check if numeric if not make it so
  df1_n <- nrow(df1) ## should match number of rows (samples) in df2
  df1_k <- ncol(df1) ## Should be number of principle components or "total"
  total <- df1_k ## Total Number of PCs or number of x observations
  df1_cols <- colnames(df1)
  
  
  ## Coherce into numeric column. 
  df2 <- as_numeric_mat(df2)
  df2_n <- nrow(df2) ## Number of samples should match the 
  df2_k <- ncol(df2) ## Number of Factors or number of y observations
  #n <- length(batch_factors)
  df2_cols <- colnames(df2)
  n <- df2_k ## Number of df2 columns
  
  ## Create output names
  out_names <- list(df2_cols,df1_cols)
  
  pval <- matrix(0,nrow=n, ncol=total, dimnames = out_names)
  cors <- matrix(0,nrow=n, ncol=total,dimnames = out_names)
  
  for (k in seq(n)){
    ## Select a column from df2  
    a_df2_col <- df2[,k]
    
    for (i in seq(total)){
      ## Select a column from df1 
      a_df1_col <- df1[,i]
      
      ## Perform correlation and save cor & pval
      cor.out <- cor.test(x =a_df1_col,y =a_df2_col, method= c('pearson'))
      pval[k,i] <- cor.out$p.value
      cors[k,i] <- cor.out$estimate
    }
  }
  
  
  

  # Combine the two long data frames
  plotting <- reshape::melt(pval)
  colnames(plotting) <- c(y,x,'p.value')
  
  
  buffer <- reshape::melt(cors)
  colnames(buffer) <- c(y,x,'cor')
  
  ## More search space use when the two dfs are unordered
  #plotting <- left_join(plotting,buffer)
  plotting$cor <- buffer$cor ## Add the cor column to the df. 
  

  #buffer <- reshape::melt(one.sum$importance[2,][1:total])
  #buffer$PC <- rownames(buffer)
  #buffer <- buffer %>% as.data.frame() %>% dplyr::rename(perc = value)
  
  #plotting <- left_join(plotting,buffer)
  
  # Make the labels
  plotting <- plotting %>% mutate(lab = round(cor,3))
  plotting <- plotting %>% mutate(stars = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")))
  
  # Make X axis 
  #plotting <- plotting %>% mutate(PC = gsub('PC','',PC), PC = sprintf('%s (%g%%)', PC, round(perc *100, 0) ))
  #plotting <- plotting %>% mutate(PC = gsub('PC','',PC), newPC = sprintf('%s\n(%g%%)', PC,round(perc *100, 0) ))
  #plotting <- plotting %>% mutate(newPC = sprintf('%s\n(%g%%)', PC, round(perc *100, 0) ))
  ## JG 10/15/20: Adding levels ot newPC so the x axis is stable
  #plotting <- plotting %>% mutate(newPC = factor(newPC, levels = unique(newPC)))
  return(plotting)
}


## JG 10/20/20: Correlation Heatmap Function
cor_heatmap <- function(plot_df, X = "PC_label", Y = "Batch", fill_col = "cor", lab_col = "stars", title = ""){
  ## JG 10/15/20: I dont want to colors to reflect the pvalues like that I want it to be colored by significance
  #heatmap <- ggplot(plotting,aes(PC,batch,fill=p.value)) + geom_tile() + 
  #  scale_fill_gradient(low = 'red', high="white", breaks= c(-Inf,0.001, 0.01, 0.05,Inf))+ 
  #  geom_text(aes(label=stars)) +
  #  ggtitle(sprintf('Heat Plot of Batch Effect Correlations (pearsons)\n%s', title))
  
  ## MUST USE QUASIQUOTATION TO WORK WITH VARIABLES :)
  #https://rlang.r-lib.org/reference/nse-force.html
  
  ## I will add the color based on the correlation 
  heatmap <- ggplot(plot_df,aes(!!sym(X),!!sym(Y),fill=!!sym(fill_col))) + geom_tile() + 
    scale_fill_gradient2(low = scales::muted('blue'), mid = "white", high=scales::muted("red"), breaks= c(-1, -.5, 0, .5,1))+ 
    geom_text(aes(label=!!sym(lab_col))) +
    ggtitle(sprintf('Heat Plot of Batch Effect Correlations (pearsons)\n%s', title))
  
  #print(heatmap)
  return(heatmap)
  
} # End of PCA assessment plotting function





## JG 10/16/20: Correlation Distance Function
## This function takes in a matrix and computes a correlation matrix to be used in distace. 
dist.cor <- function(x, method = "pearson", abs=TRUE, diag=FALSE, upper=FALSE){
  
  xcor <- cor(x, method = method)
  
  if(abs){
    xcor = 1-abs(xcor)
  } else{
    xcor = 1-xcor
  }
  
  d <- xcor[lower.tri(xcor,diag=FALSE)]
  
  
  attr(d, "Size") <- nrow(x)
  attr(d, "Labels") <- dimnames(x)[[1L]]
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- method
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  
  return(d)
}## End of dist.cor

## JG 10/21/20: See function below; I also made it for the mds function. 
JG_dist <- function(mat, dist.method = "euclidean"){
  ## Check to see if it is cor 
  if (dist.method == "cor"){
      dist_mat <- dist.cor(mat)
    } else{ ## Overwise use the standard method. 
      dist_mat <- dist(scale(t(mat)), method = dist.method)
    }
  return(dist_mat)
} ## End of JG_dist


## JG 10/21/20: Plot hclust
jg_plot_hc <- function(hc_out, pd,batch_factors){
  ## Plot the dendrogram
  plot(hc_out)
  
  ## Make serial plots of the 
  for (b in batch_factors){
    #print(b)
    new_hc <- hc_out
    #print(c_mappings[b])
    ## CHECK HERE IF FUNCTION FAILS! CONSIDER ADDING MORE CHECKS FOR INPUT. 
    new_hc$labels <- pd[b] %>% pull()
    plot(new_hc, main = b)
    
  }
}

## JG 10/16/20 Making automated hclust function to return fully analyzed data of matrix. 
## NOTE: PLEASSE t() THE DATA BEFORE PASSING INTO THE FUNCTION
## See ?dist for the full list of methods (distance measure to be used):  This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
## I will add a custom if statement for cor which will use my custom dist.cor function
## I also allowed the function to pass along a clustering method too: This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
## I am not sure about ALL the options but its better to be flexible. 
assess_hcclust <- function(mat,pd,batch_factors,dist.method = "euclidean", hclust.method = "complete"){
  
  ## JG 10/21/20: Made distance a callable function
  dist_mat <- JG_dist(mat,dist.method)
  ## Check to see if method == "cor"
  # if (dist.method == "cor"){
  #   dist_mat <- dist.cor(mat)
  # } else{
  #   dist_mat <- dist(scale(t(mat)), method = dist.method)
  # }
  #dist_mat <- dist(scale(t(clean_mat)))
  #cor_dist <- dist.cor(t(clean_mat))
  
  ## Compute Hierarcical Clustering 
  #hc_out <- hclust(cor_dist,method = "complete")
  #hc_out <- hclust(cor_dist,method = "ward.D")
  hc_out <- hclust(dist_mat,method = hclust.method)
  
  out <- list()
  out$hc <- hc_out
  out$plothc <- jg_plot_hc
  return(out)
  #jg_plothc
  
  
}## End of assess_hcclust

## JG 10/20/20: Make Numeric Matrix
## I found out I need to do this numeric matrix function in other parts of this script so here I go!
as_numeric_mat <- function(df){
  ## Make sure is numeric and make it numeric if it isnt.  
  rows <- rownames(df)
  cols <- colnames(df)
  not.numeric <- function(x) !is.numeric(x)
  if(any(mapply(df,FUN=not.numeric))){
    df <- mapply(df,FUN=as.numeric)
  }
  
  rownames(df) <- rows
  colnames(df) <- cols
  return(df)
}

## JG 10/20/20: Corrplot Function
## THis function will accept a matrix, assumed to be the pheno matrix with all columns being coherced into numeric and plotted
## https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
df_corrplot <- function(df, siglvl = .05){
  ## Make sure is numeric and make it numeric if it isnt. 
  df <- as_numeric_mat(df)
  
  corM <- cor(df)
  resM <- cor.mtest(df)
  ## Position Title: https://stackoverflow.com/questions/40509217/how-to-have-r-corrplot-title-position-correct
  corrplot(corM, p.mat = resM$p, sig.level = siglvl, title = "Correlation Plot",mar=c(0,0,1,0))
  plt <- recordPlot()
  return(plt)
}


##JG 10/21/20: Factorized MDS function this accepts the first 3 parameters of the assess function and gives back a fit object. 
## Consider making multiple MDS models available with parameter like assess_hclust. 
JG_MDS <- function(dat,numPositions= 1000, dist.method = "euclidean", mds.method = "classic"){
  
  # Compute eucledian distance 
  Mpca <- dat
  o <- order(rowVars(Mpca), decreasing =T)[seq(numPositions)]
  dist_in <- Mpca[o, ]
  #d <- dist(t(Mpca[o, ]))
  d <- JG_dist(dist_in, dist.method) ## Compute distance with method of choice. 
  
  ## Fit Multidimensional scaling model
  ## 10/21/20: Just classic for now.
  if (mds.method == "classic"){
    fit <- cmdscale(d, eig = TRUE) ## 10/22/20 Added list output to make autoplot. 
  }## Can insert more function calls here! 
  class(fit) <- "mds"
  
  return(fit)
}## end of JG_MDS


## JG 10/21/20: Now Factorizing code. I will modify the JG_mdsplot function to accept a fit MDS plot object and parameters to make a plot. 
assess_MDS <- function(dat,numPositions= 1000, dist.method = "euclidean"){
  
  ## Compute MDS with input
  fit <- JG_MDS(dat,numPositions,dist.method)
  
  out <- list()
  out$mds <- fit
  out$jg_autoplot <- jg_autoplot
  
  return(out)
  
}## End of assess_MDS 


## JG 10/20/20: IT HAPPENED I took some of my old code and it should work out of the box. 
## Consider refactoring code to return mds and harmonizing pca (biplot) and mds plotting function. 
## THIS IS NOW FACOTRIZED. Modifying so that it only accepts the fit object and that is it. 
## Jacob Gutierrez 6/10/19 gutierj@ohsu.edu
## This is a custom plot function that performs and MDS plot and allows for extra meta data to be plotted along side
#require(RColorBrewer)
JG_mdsPlot <- function(mdsfit, sampNames = NULL, 
                      sampGroups = NULL, xlim, ylim, pch = 21, 
                      pal = brewer.pal(8,'Dark2'), legendPos = 'bottomleft', 
                      main = NULL,rev = F,shapes = F , 
                      xyleg =list(one = c(25,20), two=c(25,30)), legtitles = list(one = 'one', two= 'two') ){
  
  
  fit <- mdsfit
  
  ## Starting Plot!!! -
  # calculate ranges now that data is generated
  xlim <- range(fit[, 1]) * 1.5
  ylim <- range(fit[, 2]) * 1.5
  
  
  # find sample groups
  if (is.null(sampGroups)) sampGroups <- rep(1, numPositions)
  sampGroups <- as.factor(sampGroups)
  
  
  # extract number of colors to make
  colourCount <- length(levels(sampGroups))
  
  # Check for reverse
  if(rev){
    sampGroups <- fct_rev(sampGroups)
  }
  
  
  ## Generate color pallet if the pallet is large enough then we good otherwise make gradient
  if(colourCount > length(pal)){
    #print('EXPAND')
    getcol <- colorRampPalette(pal)
    palx <- getcol(colourCount)
    col <- palx[sampGroups]
    
  } else{
    palx <- pal
    col <- palx[sampGroups]
    
  }
  
  
  
  
  # Plot with text. Maybe make it to plot with shapes
  # xpd allows for plotting outside the area
  if (!is.null(xyleg)){
    par(xpd = T, mar = par()$mar + c(0,0,0,5))
  }
  
  if (is.null(sampNames) | shapes) {
    if (!is.null(sampNames)){
      sampNames <- as.factor(sampNames)
      pch = c(21,22,23,24,25)[sampNames]
    }
    plot(
      x = fit[, 1],
      y = fit[, 2],
      pch = pch,
      bg = col,
      lwd = 1,
      xlim = xlim,
      ylim = ylim,
      xlab = "",
      ylab = "",
      main = main)
  } else {
    plot(
      x = 0,
      y = 0,
      type = "n",
      xlim = xlim,
      ylim = ylim,
      xlab = "",
      ylab = "",
      main = main)
    text(x = fit[, 1], y = fit[, 2], sampNames, col = col)
  }
  
  legendNCol <- colourCount
  
  if (!is.null(xyleg)){ # dont plot legends otherwise
    if (colourCount >= 1) {
      if (legendNCol > 10) legendNCol = as.integer(legendNCol/2)
      legend(xyleg$one[1], xyleg$one[2],
             #x = legendPos,
             title = legtitles$one,
             legend = levels(sampGroups),
             ncol = 1,
             text.col = palx[seq_len(colourCount)],
             cex = 1.0,
             pt.cex = 1.5)
      #xpd = TRUE)
      #title.cex = .75)
    }
    
    if (length(pch) > 1){
      #print('conditon?')
      shape_nrow = length(levels(sampNames))
      legend(xyleg$two[1], xyleg$two[2],
             #x = 'topleft',
             title = legtitles$two,
             legend = levels(sampNames),
             ncol = 1,
             #text.col = palx[seq_len(shape_nrow)],
             pch = unique(pch),
             #col = palx[seq_len(shape_nrow)],
             cex = 1.0,
             pt.cex = 1)
      #xpd = TRUE)
      #title.cex = .75)
    }
  }
  
  if (!is.null(xyleg)){
    par(mar=c(10, 10, 5, 5))
  }
  
}## END of MDS function

## Making a simple tidy command to return long MDS data. 
## DOESNT HELP GGBIPLOT
jg_mds_tidy <- function(mdsobj){
  
  ## catch data from mds object and name columns
  dat <- mdsobj$points %>%as.data.frame()
  colnames(dat) <- paste(seq(ncol(dat)))
  
  ## Add Rownames
  dat <- dat %>% rownames_to_column(var = "row") %>% relocate #%>% select(2,3,1)
  
  ## make longer and put the x and y names in the first column. 
  #dat <- pivot_longer(dat, cols= -row, names_to = "mds", values_to = "value") %>% select(2,3,1)
  
  return(dat)
  
}


## jgfortify ## This takes matrix (or df) input, adds the rownames on. 
## Then sorts the two matrices by rowname and merges them together. 
jg_fortify <- function(dat,pd){
  
  ## Check if it is a matrix otherwise continue and assume is dataframe
  if (is.matrix(dat)) dat <- dat %>% as.data.frame()
  
  ## add rownames to both
  dat <- dat %>% rownames_to_column(var = "row") %>% relocate(row, .after = last_col()) #%>% select(2,3,1)
  
  pd <- pd  %>% rownames_to_column(var = "row")
  
  
  out_df <- left_join(dat,pd,by="row")
  rownames(out_df) <- out_df$row
  
  return(out_df)
  
}

## JG 10/22/20: PCA & MDS Autoplot Function
## I am now scrapping useful pieces of code from the 
jg_autoplot <- function(obj, data = NA , x = 1, y = 2, scale = 1, color = NULL, label = NA, label.size = 3, shape = TRUE){
  
  ## Prepare plotting data
  if ( is(obj, "mds")){
    #print("mds")
    
    mds_data <- obj$points
    colnames(mds_data) <- paste("MDS",seq(ncol(mds_data)), sep = "_")
    #colnames(mds_data) <- c("1","2")
    ## Add phenotype data onto mds points
    
    ## JG 10/23/20: Adding data == NA conditional to make this script work. 
    if (is.na(data)){
      plot.data <- mds_data %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "row") %>% 
        relocate(row, .after = last_col()) #%>% select(2,3,1)
      
    } else{ ## Plot data exists. 
    plot.data <- jg_fortify(mds_data,data)
    }
    
    ## Add labels
    xlab = "MDS1"
    ylab = "MDS2"
    
    #plt <- ggplot(plot.data,aes(x = !!sym(xlab), y = !!sym(ylab))) + geom_point(aes(color =!!sym(color)))
    
  } else if (is(obj,"prcomp")){
    #print("prcomp")
    
    ## Make Plot data and accessory data to improve later. 
    ## Directly taken from ggfortify
    #plot.data <- tidy(obj) ## the ggplot2::fortify command is being defunct just use broom package
    
    ## Do the same as fortify
    ve <- obj$sdev^2 / sum(obj$sdev^2)
    PC <- paste0("PC", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    loadings.column <- 'rotation'
    
    lam <- obj$sdev[c(x, y)]
    #lam <- lam * sqrt(nrow(plot.data))
    lam <- lam * sqrt(length(ve))
    
    
    ## Select the two PCs to use.
    pca_data <-obj$x
    pca_names <- colnames(pca_data) %>% sprintf('%s\n(%g%%)', ., round(ve *100, 0) )
    #colnames(pca_data) <- pca_names 
    pca_data <- pca_data[,c(x,y)]
    
    ## Scaling data
    # scaling
    if (scale != 0) {
       lam <- lam ^ scale
       pca_data <- t(t(pca_data) / lam)
     }
    
    xlab = pca_names[x]
    ylab = pca_names[y]
    #pca_labels <- pca_names[c(x,y)]
    
    ## JG 10/23/20: Adding data == NA conditional to make this script work. 
    if (is.na(data)){
      plot.data <- pca_data %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "row") %>% 
        relocate(row, .after = last_col()) #%>% select(2,3,1)
      
    } else{ ## Plot data exists. 
      plot.data <- jg_fortify(pca_data,data)
    }
    
    
  } else{
    print("Unrecognized input")
    ## Add stop function here. 
  }
  
  
  ## Unpack Label Function
  ## JG 10/23/20: adding label = T conditiona; 
  if (is.na(label)){
    label.bool = FALSE
    label.col = "row"
  } else if( is.character(label)){
    label.bool = TRUE
    label.col = label
  } else if (isTRUE(label)){
    ## Exactly as is.na but I want to see the labels. 
    label.bool = TRUE
    label.col = "row"
  }

  
  ## My idea is to make the ggbiplot function to think all the data is acceptable prcomp input although it may be masked MDS data.
  #plt <- autoplot(object = obj,data = data, colour = color, label = label, label.size = label.size )
  #plt <- ggbiplot(plot.data = plot.data)
  #ggplot(plot.data,aes(x = !!sym(xlab), y = !!sym(ylab))) + geom_point(aes(color =!!sym(color)))
  #print(plot.data)
  plt <- ggbiplot(plot.data, colour = color,
                  label = label.bool, label.label = label.col, 
                  shape = shape,
                  xlab = xlab, ylab = ylab)
  
  return(plt)
  
  
}## End of jg_autoplot()


## END OF EDA_analysis_tools.R
