library(netgsa, lib.loc = "")
library(Hotelling, lib.loc = "")
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
  ##************************************
  ## Estimating network for each pway
  ##************************************
  pway_results <- lapply(1:nrow(dataset$pathway_matrix), function(pway){
    
    #Pathway matrix subset
    pway_pathway_matrix <- dataset$pathway_matrix[pway, dataset$pathway_matrix[pway,] == 1, drop = FALSE]
    pway_genes <- colnames(pway_pathway_matrix)
    
    #Data subset
    pway_dysreg_dat <- dysreg_dat_full[rownames(dysreg_dat_full) %in% pway_genes, ]

    #Edgelist subset
    pway_database_edges <- data.table::copy(dataset$database_edges)
    pway_database_edges$edgelist <- pway_database_edges$edgelist[(paste0(base_id_src, ":", base_gene_src) %in% pway_genes) & (paste0(base_id_dest, ":", base_gene_dest) %in% pway_genes)]
    pway_genes_not_in_dbs <- paste0(names(pway_database_edges$genes_not_in_dbs), ":", pway_database_edges$genes_not_in_dbs) %in% pway_genes
    pway_database_edges$genes_not_in_dbs <- if(any(pway_genes_not_in_dbs)){pway_database_edges$genes_not_in_dbs[pway_genes_not_in_dbs] } else{ NULL}
    
    
    cat(paste0("..... ", "Evaluating pathway: ", pway," with ", length(pway_genes), " genes \n"))
    start_adj_cluster_F <- Sys.time()
      adj_info_cluster_F <- prepareAdjMat(pway_dysreg_dat, group = dataset$gene_groups, databases = pway_database_edges, cluster = FALSE, lambda_c = dataset$lambda_noncluster, penalize_diag = penalize_diag)
    stop_adj_cluster_F  <- Sys.time()
    time_adj_cluster_F <- stop_adj_cluster_F - start_adj_cluster_F
    
    adj_mats_zero <- sapply(adj_info_cluster_F[["Adj"]][1:2], function(A) {sapply(A, function(a) all(a == 0))})
    if( all(adj_mats_zero)){
      used_hotelling = TRUE
      start_net_cluster_F <- Sys.time()
        ds1 <- pway_dysreg_dat[,dataset$gene_groups == unique(dataset$gene_groups)[1]]
        ds2 <- pway_dysreg_dat[,dataset$gene_groups == unique(dataset$gene_groups)[2]]
        hotel <- hotelling.test(x = t(ds1), y = t(ds2))
        results <- data.frame(pathway = rownames(pway_pathway_matrix), pSize = length(pway_genes), pval = hotel$pval, pFdr = hotel$pval, teststat = hotel$stats$statistic*hotel$stats$m)
        betas <- list(matrix(rowMeans(ds1), ncol = 1), matrix(rowMeans(ds2), ncol = 1))
        out_cluster_F <- list(results = results, beta = betas, "s2.epsilon" = NA, "s2.gamma" = NA)
      stop_net_cluster_F <- Sys.time()
    } else {
    used_hotelling = FALSE
    start_net_cluster_F <- Sys.time()
      out_cluster_F  <- NetGSA(A = adj_info_cluster_F[["Adj"]],  x = pway_dysreg_dat,  group=dataset$gene_groups, pathways = pway_pathway_matrix, lklMethod = "REML")
    stop_net_cluster_F <- Sys.time()
    
    
    }
    
    time_net_cluster_F <- stop_net_cluster_F - start_net_cluster_F
    timing_df <- data.frame(pathway = rownames(pway_pathway_matrix), seed = iter, dysregulation_framework = dysreg, mu = muval, dataset_name = dataset_name, penalize_diag = penalize_diag,
                            time_adj_cluster_F = time_adj_cluster_F, time_net_cluster_F = time_net_cluster_F, used_hotelling = used_hotelling)
    
    
    return(list(res = out_cluster_F[c("results", "beta", "s2.epsilon", "s2.gamma")], time = timing_df))
  })
  
  #Compile results from p-way level analysis to be consistent with overall
  out_results <- do.call(rbind, lapply(pway_results, function(x) x[["res"]][["results"]]))
  ##FDR adjust p-values...
  out_results[["pFdr"]] <-  p.adjust(out_results[["pval"]],"BH")
  
  out_beta <- lapply(1:ncond, function(i) do.call(rbind, lapply(pway_results, function(x) x[["res"]][["beta"]][[i]])))
  out_s2_epsilon <- do.call(rbind, lapply(pway_results, function(x) x[["res"]][["s2.epsilon"]]))
  out_s2_gamma <- do.call(rbind, lapply(pway_results, function(x) x[["res"]][["s2.gamma"]]))
  
  results_noncluster <- list(results = out_results, beta = out_beta, "s2.epsilon" = out_s2_epsilon, "s2.gamma" = out_s2_gamma)
  
  #Compile timing df
  timing_df <- do.call(rbind, lapply(pway_results, function(x) x[["time"]]))
  
  
  
  #Return
  #list(p_vals = data.frame(pathways = , cluster_true_p_vals = , cluster_false_p_vals = , dysreg = , iter = , mu = ),
  #     timing = data.frame(adj_cluster_T_timing = , adj_cluster_F_timing = , netgsa_cluster_T_timing = , netgsa_cluster_F_timing = ))
  return(list(results_noncluster = results_noncluster,
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

options(warn = 1)
curr_data <- (if (dataset_name == "BreastCancer"){readRDS("./Inputs/breastcancer2012.rds")}
              else if(dataset_name == "ProstateCancer"){readRDS("./Inputs/prostatecancer2015.rds")}
              else {stop("This is ProstateCancer and BreastCancer script. Change R script in execute_one_pwr_reml_pathway.sh")})
#Run power calculation
power_one_iteration <- calcPower(muval=mu, iter=iter, dysreg=dysreg, dataset = curr_data, dataset_name = dataset_name, penalize_diag=penalizediag)

#Save
file_save_name <- paste0(dataset_name, "_", mu, "_", penalizediag, "_", dysreg, "_", iter)
save_path <- "./Power results reml pathway/"
saveRDS(power_one_iteration, paste0(save_path, file_save_name, "_reml_pway.rds"))
