## JG SVA functions 
## 1/20/20 gutierja@ohsu.edu

### Libraries
require(sva) # For vanilla SVA
#require(isva) # nm.sv helper function
#require(SmartSVA) # For smartSVA function
require(broom) # for tidy function
require(gplots) # textplot

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
  
  heatmap <- ggplot(plotting,aes(PC,variable,fill=pr2)) + geom_tile() + 
    scale_fill_gradient(low = 'white', high="red", breaks= c(0, 1,5, 10,20,30,40))+ 
    ggtitle(sprintf('heatplot of partial R2 for each batch variable
for all %i PCs', pc_n)) + geom_text(aes(label=labs)) + theme(legend.position = "left")
  
  print(heatmap)
  
  
  out <- list()
  out$plot <- plot
  out$r2wtprop <- out_mat
  out$pcheatmap <- heatmap
  return(out)
  
}



## Custom Row Invariance removal function 
remove_invariance <- function(mat){
  # compute row variance
  all_rows <- matrixStats::rowVars(mat, na.rm = TRUE)
  
  
  ## Find invariant posoitons where == 0
  invariant_idxs <- all_rows == 0   
  
  # Now say some things. 
  nrows <-  nrow(mat)
  to_remove <- sum(invariant_idxs)
  perc <- round(to_remove/nrows, 3) * 100 
  str1 <- sprintf("Total Number of CpGs to compute SVA with: %i", nrows)
  str2 <- sprintf("Number of CpGs that are invariant across all samples (removed): %i",to_remove)
  str3 <- sprintf("As percentage of total complete NA's: %f", perc)
  print(str1)
  print(str2)
  print(str3)
  
  new_mat <- mat[!invariant_idxs,]
  return(new_mat)
}


## Custom function to make pretty output for analysis
assess_diffmethobj <- function(obj, total_cpgs, diff = 10, thres = .05){
  
  
  print(sprintf("Subsetting to Signficaint sites with methylation difference greater than %i and adjusted pvalues less than %.1f", diff, thres))
  
  print(summary(obj))
  
  myDiff10p.hyper=getMethylDiff(obj,difference=diff,qvalue=thres,type="hyper")
  myDiff10p.hypo=getMethylDiff(obj,difference=diff,qvalue=thres,type="hypo")
  myDiff10p=getMethylDiff(obj,difference=diff,qvalue=thres)
  
  #total_cpgs <- dim(meth)[1]
  
  total_size <- dim(myDiff10p)[1]
  hyper <- dim(myDiff10p.hyper)[1]
  hypo <- dim(myDiff10p.hypo)[1]
  not_sig <- total_cpgs - total_size
  
  
  display <- data.frame(bav_tav_sig = c(hyper, not_sig, hypo))
  rownames(display) <-  c("hyper","not_sig","hypo")
  print(display)
  
}


## JG 10/15/20: This function does not accept tibbles they have specialized behavior that makes this script sad.
assess_pca <- function(mat, pd,batch_factors,title){
  
  ## NA's mess up PCA so I will remove them here
  mat <- na.omit(mat)
  
  ## ------------------------------------------------------------------------
  total <- 10
  one.pca <- prcomp(t(mat))
  one.sum <- summary(one.pca)
  prop_var <- one.sum$importance[2,1:total] # second row is proportion explained
  
  
  ## ------------------------------------------------------------------------
  # Skee plot
  plot(prop_var, type='b')
  title(main=sprintf('Skree Plot: Amount of Variance Explained in each PC\n%s',title))
  
  # Cum Sum plot
  plot(cumsum(prop_var), type = 'b')
  title(main=sprintf('Cumulative Sum of Variance Explained\n%s',title))
  
  
  ## ------------------------------------------------------------------------
  #require(gplots) # textplot
  txt_plot <- textplot(capture.output(one.sum),valign='top', mar =c(5,4.5,4,2))
  txt_plot
  
  ## ----eval = F------------------------------------------------------------
  ## 
  ## autoplot(one.pca,x=2,y=3)
  
  
  ## ------------------------------------------------------------------------
  n <- length(batch_factors)
  pca.pval <- matrix(0,nrow=n, ncol=total)
  pca.cors <- matrix(0,nrow=n, ncol=total)
  
  for (k in seq(n)){
    # select factor of interest
    
    ftr.k <- batch_factors[k] 
    #print(k)
    
    for (i in seq(total)){
      # Loop through PCs computing the correlation
      #print(i)
      pc.i <- one.pca$x[,i]
      #print(ftr.k)
      a_batch <- pd[,ftr.k] 
      
      if (is.factor(a_batch)) {
        a_batch <- as.numeric(a_batch)
      }
      
      ## Print PC
      #print(i)
      #print(length(pc.i))
      #print(length(a_batch))
      cor.out <- cor.test(x =pc.i,y =a_batch, method= c('pearson'))
      pca.pval[k,i] <- cor.out$p.value
      pca.cors[k,i] <- cor.out$estimate
    }
  }
  
  rownames(pca.pval) <- batch_factors
  colnames(pca.pval) <- paste('PC', seq(total), sep = '')
  #pca.pval
  
  rownames(pca.cors) <- batch_factors
  colnames(pca.cors) <- paste('PC', seq(total), sep = '')
  
  
  ## ------------------------------------------------------------------------
  # Combine the two long data frames
  plotting <- reshape::melt(pca.pval)
  colnames(plotting) <- c('batch','PC','p.value')
  
  buffer <- reshape::melt(pca.cors)
  colnames(buffer) <- c('batch','PC','cor')
  
  plotting <- left_join(plotting,buffer)
  
  buffer <- reshape::melt(one.sum$importance[2,][1:total])
  buffer$PC <- rownames(buffer)
  buffer <- buffer %>% as.data.frame() %>% dplyr::rename(perc = value)
  
  plotting <- left_join(plotting,buffer)
  
  # Make the labels
  plotting <- plotting %>% mutate(lab = round(cor,3))
  plotting <- plotting %>% mutate(stars = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")))
  
  # Make X axis 
  #plotting <- plotting %>% mutate(PC = gsub('PC','',PC), PC = sprintf('%s (%g%%)', PC, round(perc *100, 0) ))
  #plotting <- plotting %>% mutate(PC = gsub('PC','',PC), newPC = sprintf('%s\n(%g%%)', PC,round(perc *100, 0) ))
  plotting <- plotting %>% mutate(newPC = sprintf('%s\n(%g%%)', PC, round(perc *100, 0) ))
  ## JG 10/15/20: Adding levels ot newPC so the x axis is stable
  plotting <- plotting %>% mutate(newPC = factor(newPC, levels = unique(newPC)))
  
  ## JG 10/15/20: I dont want to colors to reflect the pvalues like that I want it to be colored by significance
  #heatmap <- ggplot(plotting,aes(PC,batch,fill=p.value)) + geom_tile() + 
  #  scale_fill_gradient(low = 'red', high="white", breaks= c(-Inf,0.001, 0.01, 0.05,Inf))+ 
  #  geom_text(aes(label=stars)) +
  #  ggtitle(sprintf('Heat Plot of Batch Effect Correlations (pearsons)\n%s', title))
  
  ## I will add the color based on the correlation 
  heatmap <- ggplot(plotting,aes(newPC,batch,fill=cor)) + geom_tile() + 
    scale_fill_gradient(low = 'lightblue', high="red", breaks= c(-1, -.5, 0, .5,1))+ 
    geom_text(aes(label=stars)) +
    ggtitle(sprintf('Heat Plot of Batch Effect Correlations (pearsons)\n%s', title))
  
  print(heatmap)
  
  #+ scale_fill_gradient(low="#D7191C",  high="white")+ geom_text(aes(label=stars))
  out <- list()
  out$pca <- one.pca
  out$heatmap <- heatmap
  ## JG 10/15/20: Adding the plotting data
  out$plot_df <- plotting
  
  
  ## CONSIDER CAPTURING THE PLOTS TO BE USED IN THE MARKDOWN... 
  return(out)
} # End of PCA assessment plotting function
