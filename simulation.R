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
#############################################################################################################
#############################################################################################################
#############################################################################################################
start <- Sys.time()

source('SimulateRNASeqV4_SRZ.R')
source('kMEn.R')
load('../data/GO_list.RData')

library(data.table)
library(stringr)
library(parallel)
library(doParallel)
library(snow)
library(snowfall)
#library(foreach)

#no_cores = detectCores()
#registerDoParallel(no_cores)

##################################################
###
### Parameters of Interest
### - Geneset: 40|200
### - % responsive: 5%|10%|25%
### - FC of Responsive Transcripts: ~N(mu=4, sigma^2)
### - % up v. down regulated: 100%|50%|25%
### - Num Patients: 10|20|30
###
##################################################

### Load config file & set params
go.bp <- read.delim2("../data/go_bp_filtered15-500.txt", stringsAsFactors = F)
neg.bin.param_brain <- fread('../data/negative_binomial_parameters_brain.csv',header=T); names(neg.bin.param_brain)[1] <- 'Gene'
setkey(neg.bin.param_brain, Gene)

#### restrict gene ontologies for which we have all parameter values estimated
#### from brain data 

GO_list_lengths = unlist(lapply(GO_list,length))

go_list_size40_reduced <- Match_FullList(GO_list = GO_list, 
                                         pathways_of_interest =  names(which(GO_list_lengths==40)), 
                                         tissue_specific_params = neg.bin.param_brain)

go_list_size200_reduced <- Match_FullList(GO_list = GO_list, 
                                          pathways_of_interest =  names(which(GO_list_lengths>150 & GO_list_lengths < 250)), 
                                          tissue_specific_params = neg.bin.param_brain)

go_list_40_200 <- list(go_list_size40 = go_list_size40_reduced, 
                       go_list_size200=go_list_size200_reduced)

rm(go_list_size40_reduced, go_list_size200_reduced)



## test all conditions 
# simulate N patients with same condition,

##############################
### need to create a function 
### that will take a set of 
### parameters, and run the 
### simulation m=2000 times 
### and save the

get_t.test_kmeans <- function(i){
	
			cohort <- generate_cohort(neg.bin.param= neg.bin.param , 
                        	            DEG_pct = deg, 
                                	    Up_REG = up, 
                                	    N=N, 
                                	    geneset_size=size,
					    go_list_40_200 = go_list_40_200, 
                                	    verbose = F)
			
			geneset_name = cohort$geneset_name
			pvals <- get_t.test_pvalues(simulated_matrix = cohort$simulated_matrix )
			preds <- get_kmeans_results(simulated_matrix = cohort$simulated_matrix, geneset_name=geneset_name )
	
			return(list(T_test =pvals, KMeans = preds, Cohort=cohort))
}

args = commandArgs(trailingOnly=TRUE)
print(args)

reps=28
size = as.numeric(args[1])
N    = as.numeric(args[2])
deg  = as.numeric(args[3])
up   = as.numeric(args[4])
file_name = as.character(args[5])

print(c(size,N,deg,up,reps))
neg.bin.param=neg.bin.param_brain

sfInit(parallel=T, cpus=detectCores(), type='SOCK')
sfLibrary(snowfall)
sfSource('kMEn.R')
sfSource('SimulateRNASeqV4_SRZ.R')
sfExport('go_list_40_200','neg.bin.param','size','deg','up','N','GO_list','reps','get_t.test_kmeans')

#cl <- makeCluster(no_cores, type='FORK')

full_list = sfLapply( 1:reps, function(x) get_t.test_kmeans(x))
sfStop()

setwd('/rsgrps/yves/samir/simulation_study2017/results/output')
save(full_list, file=paste(file_name,'.R',sep=''))

print(paste('saved ', file_name,'.R in output', sep=''))

print(start - Sys.time())
setwd('/rsgrps/yves/samir/simulation_study2017/code')
