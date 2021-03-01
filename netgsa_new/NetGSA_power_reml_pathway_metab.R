library(netgsa, lib.loc = "")

setwd("")

##********************************
##****   Main function   *******
##********************************
calcPower <- function(muval, iter, dysreg, dataset, dataset_name, penalize_diag){
  p            <- dim(dataset$dataset[[1]])[1]
  n            <- sapply(dataset$dataset, ncol)
  dysreg_genes <- dataset[[dysreg]]
  ncond        <- length(dataset$dataset)
  stopifnot(all(dysreg_genes %in% rownames(dataset$dataset[[1]])) & !is.null(dysreg_genes))
  mu <- list(rep(0,p), 
             muval * (match(rownames(dataset$dataset[[1]]), dysreg_genes, nomatch = 0)>0))

  #Set-up data
  set.seed(iter) #Should reproduce exactly the epislonVar
  dysreg_dat <- vector("list", ncond)
  for (i in 1:ncond) {
    ## Add some random noise to the data
    epsilonVar <- matrix(rnorm(p * n[i]), p, n[i])
    dysreg_dat[[i]] <- dataset$dataset[[i]] + matrix(rep(mu[[i]], n[i]), p, n[i]) + 0.1*epsilonVar
  }
  dysreg_dat_full <- do.call(cbind, dysreg_dat) #Center, scale + mu + noise
  ##*******************************  
  ## Estimating network
  ##*******************************
  
  cat(paste0("..... ", "Iteration: ", iter, " Estimating Adj with Cluster = FALSE ....\n"))
  start_adj_cluster_F <- Sys.time()
      adj_info_cluster_F <- prepareAdjMat(dysreg_dat_full, group = dataset$gene_groups, databases = dataset$database_edges, cluster = FALSE, lambda_c = dataset$lambda_noncluster, penalize_diag = penalize_diag) #not soooo much faster b/c reduces to one big 2400 group
    stop_adj_cluster_F  <- Sys.time()
  time_adj_cluster_F <- stop_adj_cluster_F - start_adj_cluster_F
  
  cat(paste0("..... ", "Iteration: ", iter, " Estimating Network with Cluster = FALSE ....\n"))
  start_net_cluster_F <- Sys.time()
      out_cluster_F  <- NetGSA(A = adj_info_cluster_F[["Adj"]],  x = dysreg_dat_full,  group=dataset$gene_groups, pathways = dataset$pathway_matrix, lklMethod = "REML")
  stop_net_cluster_F <- Sys.time()
  time_net_cluster_F <- stop_net_cluster_F - start_net_cluster_F
  
  timing_df <- data.frame(seed = iter, dysregulation_framework = dysreg, mu = muval, dataset_name = dataset_name, penalize_diag = penalize_diag, time_adj_cluster_F = time_adj_cluster_F, time_net_cluster_F = time_net_cluster_F)
  
  #Return
  #list(p_vals = data.frame(pathways = , cluster_true_p_vals = , cluster_false_p_vals = , dysreg = , iter = , mu = ),
  #     timing = data.frame(adj_cluster_T_timing = , adj_cluster_F_timing = , netgsa_cluster_T_timing = , netgsa_cluster_F_timing = ))
  return(list(results_noncluster = out_cluster_F[c("results", "beta", "s2.epsilon", "s2.gamma")], 
              timing = timing_df,
              seed = iter,
              dysregulation_framework = dysreg,
              mu = muval,
              dataset_name = dataset_name,
              penalize_diag = penalize_diag,
              sessionInfo = sessionInfo()))
}

## Read in arguments
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]])) #Should load "file" into R environment
  
  #Now parse through file. Can do this more flexibly, but requires funky looking code
  #If do this in a function make sure you eval() in global environment. See ?eval
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    eval(parse(text=line))
  }
  close(con)
  
}


curr_data <- (if (dataset_name == "Metabolites"){readRDS("./Inputs/metabolites.rds")}
              else {stop("This is metabolites Rscript. Change R script in execute_one_pwr_reml_pathway.sh")})
#Run power calculation
power_one_iteration <- calcPower(muval=mu, iter=iter, dysreg=dysreg, dataset = curr_data, dataset_name = dataset_name, penalize_diag=penalizediag)

#Save
file_save_name <- paste0(dataset_name, "_", mu, "_", penalizediag, "_", dysreg, "_", iter)
save_path <- "./Power results reml pathway/"
saveRDS(power_one_iteration, paste0(save_path, file_save_name, "_reml_pway.rds"))
