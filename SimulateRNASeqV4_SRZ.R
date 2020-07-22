#########################################################################:####################################
#############################################################################################################
#############################################################################################################
###
### Simulate RNA-Seq data
### 1/2/2016 Qike Li
### revised on 2/9/2017 Samir Rachid Zaim
###
### the following code is modified based on the
### ~/Dropbox/Qike/Adaptive cutoff/DEseq related/suppl_ii.Rnw line 1094
###
#############################################################################################################w
#############################################################################################################
#############################################################################################################

library(data.table)

#############################################################################################################
#############################################################################################################
####
####  rGeneCount 
####  - input
####    @ parameters: tuple(mean, delta)
####
####  - output
####    @ rNegBinomial(1, mu, delta) V rpois(mu) 
####      --> depends on dispersion parameter
####
####
####  rGeneCountDysReg
####  
####  - input (same as rGeneCount except includes dysregulation parameter)
####    @ parameters: tuple(mean, delta)
####    @ Up_Reg: numeric [0,1] (percent up regulated)
#### 
####  - output
####    @ rNegBinomial(1, mu, delta) V rpois(mu) 
####      --> depends on dispersion parameter
####

#### rGeneCount

rGeneCount <- function(parameters){
  require(MASS)
  
  mu= parameters[1]
  delta=parameters[2]
  
  ## all zero counts handling 
  if ( mu == 0 ) {
    return(0)
  } else {
    ## Poisson case, underdispersed or mean = variance
    if ( delta <= 0 ) { 
      return(rpois(1, mu))
    } else {
      ## mu > 0 and delta > 0
      return(rnegbin(1, mu = mu, theta = 1/delta))
    }
  }
}

#### rGeneCountDysReg

rGeneCountDysReg <- function(parameters){
  require(MASS)
  mu= parameters[1]
  delta=parameters[2]
  Up_REG = parameters[3]
  Fold_Change_Dysregulation = T
  
#   print(mu)
#   print(delta)
#   print(Up_REG)
#   
  #stopifnot(0 <= Up_REG  & Up_REG<= 1)
  
  
  #### SAMPLES FROM RUNIF(3,5), RUNIF(-3,-5) WITH PROB = C(UP_REG, 1-UP_REG)
  fc_dysreg <- sample(c(runif(1, min = 3, max=5), 
                      -1*runif(1, min=3, max=5)), 
                      size=1, 
                      prob = c(Up_REG, 1-Up_REG))
  
  ## all zero counts handling 
  if ( mu == 0 ) {
    return(0)
  } else {
    ## Poisson case, underdispersed or mean = variance
    if ( delta <= 0 ) {
      
      #### UP/DOWN MULTIPLIER FOR DYSREGULATION 
      ## take max(new_mu, 0) so that dysreg_mu always >=0
      
      # dysreg_mu = max(mu + fc_dysreg * sqrt(mu),0)
    
      if(fc_dysreg<0) {
        dysreg_mu1 <- mu /abs(fc_dysreg)
        dysreg_mu2 = max(mu + fc_dysreg * sqrt(mu),0)
      } else if(fc_dysreg>0){
        dysreg_mu1 = mu * fc_dysreg
        dysreg_mu2 = max(mu + fc_dysreg * sqrt(mu),0)
      }
      
      if(Fold_Change_Dysregulation){
        dysreg_mu = dysreg_mu1
      } else {
        dysreg_mu = dysreg_mu2
      }
      
      return(rpois(1, dysreg_mu))
      
    } else {
      
      #### UP/DOWN MULTIPLIER FOR DYSREGULATION 
      ## fc * standard deviation = mu* for poisson
      ## take max(new_mu, 0) so that dysreg_mu always >=0
      
      theta = 1/delta
      
      if(fc_dysreg<0) {
        dysreg_mu1 <- mu /abs(fc_dysreg)
        dysreg_mu2 = max(mu + fc_dysreg * sqrt(mu + mu^2/theta),0)
      }else if(fc_dysreg>0){
        dysreg_mu1 = mu * fc_dysreg
        dysreg_mu2 =max(mu + fc_dysreg * sqrt(mu + mu^2/theta),0)
      }
      
      if(Fold_Change_Dysregulation){
        dysreg_mu = dysreg_mu1
      } else {
        dysreg_mu = dysreg_mu2
      }
        
      return(rnegbin(1, mu = dysreg_mu, theta = 1/delta))
    }
  }
}

#############################################################################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################
####
#### simulatePathway
####  - input
####    @ geneset: character vector | GO_List geneset
####    @ neg.bin.param: data.table | NB/Poisson empirically estimated parameters (entire set)
####    @ DysReg: boolean | simulate dysregulated (T) or regulated (F) pathway 
####
####  - output
####    @ pathway with simulated Poisson/NB counts based on empirical parameters

#### Simulate Pathway 
simulatePathway = function(geneset, neg.bin.param, DysReg = F, Up_REG=.1) {
  
  ## simulate a pathway under independence
  if(DysReg){
    sim.pathways <- unlist(apply(neg.bin.param[geneset,2:4, with=F], 1, rGeneCountDysReg))
    names(sim.pathways) <- neg.bin.param[geneset]$Gene
    
  }else{
  sim.pathways <- unlist(apply(neg.bin.param[geneset,2:3, with=F], 1, rGeneCount))
  }
  return(sim.pathways)
}

#############################################################################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################
####
#### Generate Dysregulated Pathway
####  - input
####    @ geneset: character vector | GO_List geneset
####    @ neg.bin.param: data.table | NB/Poisson empirically estimated parameters (entire set)
####    @ DEG_pct: numeric [0,1] | proportion of genes dysregulated in pathway
####    @ Up_DEG: numeric [0,1]  | proportion of dysregulated genes that are Up Regulated
#### 
####  - output 
####    @ pairedTranscriptomes: 2 x N matrix | normal, tumor samples for each gene in geneset
####    * each iteration samples a different set of genes in pathway to dysregulate to create
####      desired effect of pathway-specific (rather than biomarker specific) signal
####

generateDysregulatedPathway <- function(geneset, neg.bin.param, DEG_pct=1, Up_REG=.1, verbose=F){
  require(plyr)
  options(scipen=999)
  
  neg.bin.param$Up_Reg <- rep(Up_REG, length(neg.bin.param$Gene))
  
  pairedTranscriptomes <- data.frame(normal = (simulatePathway(geneset, neg.bin.param)), 
                     tumor   = (simulatePathway(geneset, neg.bin.param)), row.names = geneset)
  
  ngenes = length(geneset)
  true_isDE_ind <- sample(1: ngenes, size = ngenes*DEG_pct)
  DEG_status <- rep(FALSE,ngenes)
  DEG_status[true_isDE_ind] <- TRUE
  
  dysreg_genes <- geneset[DEG_status]
  pairedTranscriptomes$tumor[DEG_status] <-  simulatePathway(dysreg_genes,neg.bin.param = neg.bin.param, 
                                                             DysReg = T, Up_REG =Up_REG )
  
  pairedTranscriptomes$DEG <- DEG_status
  pairedTranscriptomes$Gene <- row.names(pairedTranscriptomes)
  
  if(verbose){
    pairedTranscriptomes <- join(pairedTranscriptomes, neg.bin.param, by='Gene')
  }
  return(pairedTranscriptomes)
}

#############################################################################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################
####
#### Generate Dysregulated Pathway
####  - input
####    @ geneset: character vector | GO_List geneset
####    @ neg.bin.param: data.table | NB/Poisson empirically estimated parameters (entire set)
####    @ DEG_pct: numeric [0,1] | proportion of genes dysregulated in pathway
####    @ Up_DEG: numeric [0,1]  | proportion of dysregulated genes that are Up Regulated
#### 
####  - output 
####    @ pairedTranscriptomes: 2 x N matrix | normal, tumor samples for each gene in geneset
####    * each iteration samples a different set of genes in pathway to dysregulate to create
####      desired effect of pathway-specific (rather than biomarker specific) signal
####

generatePairedTranscriptomes <- function(geneset, neg.bin.param, DEG_pct, Up_REG, verbose=F){
  require(plyr)
  options(scipen=999)
  
  ## generate entire transcriptome paire
  pairedTranscriptomes <- data.frame(normal = (simulatePathway(neg.bin.param$Gene, neg.bin.param)), 
                                     tumor   = (simulatePathway(neg.bin.param$Gene, neg.bin.param)), 
                                     row.names = neg.bin.param$Gene)
  
  pairedTranscriptomes$DEG<- rep(FALSE, length(neg.bin.param$Gene))
  
  ## generate dysregulated simulated pathway pair
  pathway_pairedTranscriptomes <- generateDysregulatedPathway(geneset , neg.bin.param=neg.bin.param , 
                                                              DEG_pct = DEG_pct, Up_REG = Up_REG, verbose = verbose)
  
  ### Replace values in paired transcriptome with dysregulated pathway
  idx <- which(row.names(pairedTranscriptomes) %in% geneset)
  
  pairedTranscriptomes[idx, 'normal'] <- pathway_pairedTranscriptomes$normal
  pairedTranscriptomes[idx, 'tumor'] <- pathway_pairedTranscriptomes$tumor
  pairedTranscriptomes[idx, 'DEG'] <- pathway_pairedTranscriptomes$DEG
  
  return(pairedTranscriptomes)
}


#############################################################################################################
#############################################################################################################
####
#### Generate cohort 
####  - input
####    @ neg.bin.param: data.table | NB/Poisson empirically estimated parameters (entire set)
####    @ DEG_pct: numeric [0,1] | proportion of genes dysregulated in pathway
####    @ Up_DEG: numeric [0,1]  | proportion of dysregulated genes that are Up Regulated
####    @ geneset_size: numeric {40,200} | 
####    @ go_list_40_200: list of 2 lists | list(names of pathways of size 40, names of pathways of size 200)
#### 
####  - output 
####    @ mat: 2*N x M matrix | normal, tumor samples for each gene in geneset for each N patient
####    * each iteration samples a different set of genes in pathway to dysregulate to create
####      desired effect of pathway-specific (rather than biomarker specific) signal
####

# simulate N patients with same condition, 
generate_cohort <- function(neg.bin.param=neg.bin.param_brain , 
                            DEG_pct = .2, 
                            Up_REG = 1, 
                            N=10, 
                            geneset_size=40, 
                            go_list_40_200,
                            verbose = F){

  #require(snowfall)
  require(parallel)
 
  if(geneset_size == 40){
    geneset_name <- sample(go_list_40_200[["go_list_size40"]], size = 1)
    geneset = GO_list[[geneset_name]]
    
  }else if(geneset_size==200){
    geneset_name<- sample(go_list_40_200[["go_list_size200"]], size = 1)  
    geneset = GO_list[[geneset_name]]
    
  }
  
  new_mat <- lapply(1:N, function(x) generatePairedTranscriptomes(neg.bin.param=neg.bin.param,
							 	geneset=geneset,
								DEG_pct=DEG_pct,
								Up_REG=Up_REG))
  mat = do.call(cbind,new_mat)

  return(list(simulated_matrix =mat,geneset_name = geneset_name ))
}

#############################################################################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################
####
####  get_t.test_pvalues 
####  - input
####    @ simulated_matrix| 2*NxM matrix | output from generate cohort
#### 
####  - output 
####    @ DEG <- vector of booleans | pvalue < alpha for each gene
####
####
####  get_kmeans_results 
####  - input
####    @ simulated_matrix| 2*NxM matrix | output from generate cohort
#### 
####  - output 
####    @ DEG <- vector of boolean predictions | kMEn idenfitied DEF gene
####
####

get_t.test_pvalues <- function(simulated_matrix, alpha=.05){
  
 
  #require(parallel)
  cat('\n getting t.test preds \n')
  #require(snowfall)
  #sfInit(parallel = TRUE, cpus = detectCores(), type = "SOCK")
  sfExport("simulated_matrix")
  # sfExport(simulated_matrix) #export appropriate objects that will be needed inside a function, if applicable
  # sfLibrary() #call to any special library
  get_pvals <- function(i) {
    return(t.test(x=log2(simulated_matrix[i, seq(from=1, to=length(simulated_matrix),3)]+1),
                  y=log2(simulated_matrix[i, seq(from=2, to=length(simulated_matrix),3)]+1),
		  na.action='na.exclude')$p.value)
  }
  
  pvals <- sapply(1:length(simulated_matrix[,1]),  get_pvals)
  #sfStop()

  DEG =pvals<alpha
  return(DEG)
}

get_kmeans_results <- function(simulated_matrix, geneset_name){
  cat('\n getting kmeans preds\n')
  #require(snowfall)
  #sfInit(parallel = TRUE, cpus = detectCores(), type = "SOCK", nostart=T)
  #sfExport("simulated_matrix","kMEn","f.absLogFC","neg.bin.param","GO_list","geneset_name","makeContingencyTable","cal.OR")
  
  # sfExport(simulated_matrix) #export appropriate objects that will be needed inside a function, if applicable
  # sfLibrary() #call to any special library
  get_predictions <- function(i) {
	pred_vec = kMEn(f.absLogFC(simulated_matrix[,i], simulated_matrix[,i+1]), gene.symbol = neg.bin.param$Gene, 
				GeneSet.list = GO_list[geneset_name], alternative='greater')
	return(pred_vec)
  }
  
  preds <- sfLapply( seq(from=1, to=length(simulated_matrix),3),  get_predictions)
  #sfStop()
  # preds = as.data.frame(preds); names(preds) = paste('kMEn',1:length(simulated_matrix)/3,sep='_')
  return(preds)
}

#############################################################################################################
#############################################################################################################
####
#### Match_FullList
####  - input
####    @ GO_list: list of lists | GO_bp_list
####    @ pathways_of_interest: character vector | subset of pathway to check if we have full match with tissue specific dataset
####    @ tissue_specific_params | matrix: 3xn matrix containing the estimated NB parameters from real dataset
####
####  - output
####    @ new_list: list | list of GO_bp pathways with full match in the the tissue-specific-set
####

#### eliminate all reduced lists who don't have complete genes
Match_FullList <- function(GO_list, pathways_of_interest, tissue_specific_params){
  
  new_list <- character()
  for(go_list in pathways_of_interest){
    go_list.genes <- GO_list[[go_list]]
    if(length(intersect(go_list.genes, tissue_specific_params$Gene)) == length(go_list.genes)){
      new_list = append(new_list, go_list)
    }
  }
  return(new_list)
}

#############################################################################################################
#############################################################################################################

