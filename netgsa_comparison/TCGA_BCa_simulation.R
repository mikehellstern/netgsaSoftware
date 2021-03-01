############################################
## simulation following framework of NetGSA review paper
## to compare nine gene network enrichment analysis approaches
##
##
## Kun Yue
## yuek@uw.edu
##
## edit: 11/14/2020 (set glasso threshold 1e-9)
## edit: 01/08/2021 (set glasso threhsold 1e-4)
## edit: 01/13/2021 (confirm using scaled and centered breast cancer data)
## edit: 02/07/2021 (fix colnames for topoGSA method)
############################################
## note: except NetGSA, all other methods analyze each pathway seperately. 
## TCGA-BCa data set




## Packages for different methods

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SPIA")
# BiocManager::install("ROntoTools")
# BiocManager::install("ToPASeq")
# install.packages('CePa', repos='http://cran.us.r-project.org')
# BiocManager::install("qpgraph")
# install.packages('topologyGSA', repos='http://cran.us.r-project.org')
# BiocManager::install("DEGraph")
# BiocManager::install("limma")
# BiocManager::install("PathNet")

library(devtools)
library(graphite)
library(org.Hs.eg.db)
library(SPIA)
library(ROntoTools)
library(ToPASeq)
library(CePa)
library(topologyGSA)
library(netgsa)
library(DEGraph)
library(glasso)
library(limma)
library(PathNet)
library(Matrix)
library(Hotelling)
## load data set, which contains required data formats for each method
filepath = ''
filepath = ''

setwd(filepath)



source("code/preprocess_lib.r")

set.seed(2020)

OutFolder <- "results/bca/"
filename <- "breastcancer2012_ready.rda"
newFileName <- gsub('.rda','',filename)
newFileName <- gsub('_ready','', newFileName)

# ## this load the centered and scaled breast cancer data (processed by Mike); the useful object is dataset
# processed_data = readRDS(paste0(OutFolder, newFileName, '.rds'))
# processed_data$dataset

## this load the old breast cancer dataset for the pathways and the genes2affect variables
load(paste0(OutFolder, filename))

## general parameters
today <- '20201130'
ncond <- 2
ncores <- 1


inputs = commandArgs(T)
DEBUG <- inputs[1]
muvals <- as.numeric(inputs[2])
jid <- as.integer(inputs[3])
perm <- inputs[4]

cat('settings: ', DEBUG, muvals, jid, perm, '\n')


newFileName <- paste0(newFileName, '_', DEBUG, '_perm_', perm, '_mu', muvals*10, '_ncores', ncores, '_', today, '_', jid)
check_file= tryCatch(load(paste0(OutFolder, newFileName, '.rda')), error=function(e)e)
if(! ('simpleError' %in% class(check_file))) stop('results already exist')

writeLines(capture.output(devtools::session_info()), paste0(OutFolder,newFileName,"_sessionInfo.txt"))


p <- nrow(dataset$dat)
npath <- length(pathSelected)
path.names <- names(pathSelected)
mu <- vector("list", ncond)
mu[[1]] <- rep(0, p)
mu[[2]] <- muvals * (match(dataset$gene_info$EntrezID, base::eval(as.name(paste0("genes2affect_",DEBUG))), nomatch = 0)>0)

##********************************
##****   Main function   *******
##********************************

## prepare the input dataset, change the gene IDs to make sure it can be searched in the database
entrezid_rownms <- dataset$gene_info$EntrezID[match(rownames(dataset$dat), dataset$gene_info$symbol, nomatch = NA)]
# sum(is.na(match(rownames(dataset$dat), dataset$gene_info$symbol, nomatch = NA))) # there is no missing symbols
# tmp = match(rownames(dataset$dat), dataset$gene_info$symbol, nomatch = NA)
# all(tmp == order(tmp)) # TRUE: the sequence of genes in the data is the same as in the gene_info
rownames(dataset$dat) <- entrezid_rownms

## prepare the edge lists provided to network estimation, based on searching database
try_load = tryCatch({
  load(paste0(OutFolder,today, 'edgeList_database.rda'))
}, error = function(e)e)

if('simpleError' %in% class(try_load)){
  database_search <- obtainEdgeList(rownames(dataset$dat), 
                                    c("biocarta", "kegg", "nci", "panther", "pathbank", "pharmgkb","reactome", "smpdb")) #"humancyc"
  save(list='database_search', file = paste0(OutFolder,today, 'edgeList_database.rda'))
}

## prepare the lambda sequence for the NetGSA network estimation algorithm (it seems only lambda_for_cluster is used)

## bca data
lambda_for_cluster = rev(c(2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4))
lambda_for_noncluster =  rev(c(49.5, 50, 50.5, 51, 51.5, 52, 52.5, 53))


dat <- vector("list", ncond)
dat[[1]] <- dataset$dat[,(dataset$sample_info$er_status_by_ihc=="Positive")]## p by n1
dat[[2]] <- dataset$dat[,(dataset$sample_info$er_status_by_ihc=="Negative")]## p by n2
dat <- lapply(dat, function(d) t(scale(t(d))))
n <- sapply(dat, ncol)
dat_combined <- cbind(dat[[1]], dat[[2]])



iters = jid

one_iteration = function(iters){
  set.seed(iters)
  
  new_dat <- vector("list", ncond)
  
  if(perm == 'T'){
    ## Permute sample labels
    ind.1 <- sample(ncol(dat_combined), n[1])
    new_dat[[1]] <- dat_combined[,ind.1]
    new_dat[[2]] <- dat_combined[,-ind.1] 
  }else{
    ## Or use the original sample labels
    new_dat[[1]] <- dat_combined[,seq(1,n[1])]
    new_dat[[2]] <- dat_combined[,-seq(1,n[1])]
  }

  
  dat <- new_dat
  
  ## output
  col.names <- c("PE.noCut", 
                 "PE.Cut", 
                 "SPIA", 
                 "PRS",
                 "cepa.ORA.equal.weight", 
                 "cepa.ORA.in.degree",
                 "cepa.ORA.out.degree",
                 "cepa.ORA.betweenness",
                 "cepa.ORA.in.reach",
                 "cepa.ORA.out.reach",
                 "cepa.GSA.equal.weight", 
                 "cepa.GSA.in.degree",
                 "cepa.GSA.out.degree",
                 "cepa.GSA.betweenness",
                 "cepa.GSA.in.reach",
                 "cepa.GSA.out.reach",
                 "NetGSA2_perpath_REHE", # NetGSA2 for estimating group specific network
                 "DEGraph", 
                 "topologyGSA", 
                 "topologyGSA.var", 
                 "CAMERA", 
                 "PathNet",
                 "NetGSA_perpath_REHE",  # NetGSA for estimating shared network among groups
                 'NetGSA2_perpath_REML', # perpath for analyze each pathway (and estimate nework) seperately by pathway
                 'NetGSA_perpath_REML',
                 'NetGSA2_allpath_cl_REHE', # allpath for analyze all pathways together (estimate the network accounting for all pathways)
                 'NetGSA_allpath_cl_REHE',  # cl for using clustering algorithm, ncl for not using clustering algorithm
                 'NetGSA2_allpath_ncl_REHE', # reREHE for sampling REHE, and REHE for normal REHE without sampling
                 'NetGSA_allpath_ncl_REHE',
				 'NetGSA2_allpath_cl_reREHE',
				 'NetGSA2_allpath_ncl_reREHE'
                 )
  sigInd <- matrix(NA, npath, length(col.names))
  rownames(sigInd) <- path.names
  colnames(sigInd) <- col.names
  
  ## store any warnings/comments for the pathway analysis (mainly for netgsa)
  notes = matrix('', nrow=npath+1, ncol=1, dimnames = list(c(path.names, 'all_pathways'), NULL))
  
  
  out.obj <- list()
  out.obj$pe <- 0
  out.obj$spia <- 0
  out.obj$prs <- 0
  out.obj$cepa <- 0
  
  ## store the computation time
  col.names_time <- c("PE.noCut", 
                      "PE.Cut", 
                      "SPIA", 
                      "PRS",
                      "cepa.ORA.equal.weight", 
                      "cepa.ORA.in.degree",
                      "cepa.ORA.out.degree",
                      "cepa.ORA.betweenness",
                      "cepa.ORA.in.reach",
                      "cepa.ORA.out.reach",
                      "cepa.GSA.equal.weight", 
                      "cepa.GSA.in.degree",
                      "cepa.GSA.out.degree",
                      "cepa.GSA.betweenness",
                      "cepa.GSA.in.reach",
                      "cepa.GSA.out.reach",
                      "NetGSA2_perpath_REHE_net", # NetGSA2 for estimating group specific network
                      "DEGraph", 
                      "topologyGSA", 
                      "topologyGSA.var", 
                      "CAMERA", 
                      "PathNet",
                      "NetGSA_perpath_REHE_net",  # NetGSA for estimating shared network among groups
                      'NetGSA2_perpath_REML_net', # perpath for analyze each pathway (and estimate nework) seperately by pathway
                      'NetGSA_perpath_REML_net',
                      'NetGSA2_allpath_cl_REHE_net', # allpath for analyze all pathways together (estimate the network accounting for all pathways)
                      'NetGSA_allpath_cl_REHE_net',  # cl for using clustering algorithm, ncl for not using clustering algorithm
                      'NetGSA2_allpath_ncl_REHE_net',
                      'NetGSA_allpath_ncl_REHE_net',
                      ## store separately network estimation time and analysis time for netgsa
                      "NetGSA2_perpath_REHE_analy",
                      "NetGSA_perpath_REHE_analy",  # NetGSA for estimating shared network among groups
                      'NetGSA2_perpath_REML_analy', # perpath for analyze each pathway (and estimate nework) seperately by pathway
                      'NetGSA_perpath_REML_analy',
                      'NetGSA2_allpath_cl_REHE_analy', # allpath for analyze all pathways together (estimate the network accounting for all pathways)
                      'NetGSA_allpath_cl_REHE_analy',  # cl for using clustering algorithm, ncl for not using clustering algorithm
                      'NetGSA2_allpath_ncl_REHE_analy',
                      'NetGSA_allpath_ncl_REHE_analy',
					  'NetGSA2_allpath_cl_reREHE_analy',
					  'NetGSA2_allpath_ncl_reREHE_analy'
  )
  out.obj$comp_time <- matrix(NA, nrow=1, ncol = length(col.names_time), dimnames = list(NULL, col.names_time))
  
  
  cat('... Setting up data ....\n')
  
  for (i in 1:ncond) {
    ## Add some random noise to the data
    epsilonVar <- matrix(rnorm(p * n[i]), p, n[i])
    dat[[i]] <- dat[[i]] + matrix(rep(mu[[i]], n[i]), p, n[i]) + 0.1*epsilonVar
  }
  
  ## ******************************
  ## (1) SPIA, Pathway-Express----
  ## ******************************
  cat('... Running PE ....\n')
  set.seed(iters)
  
  ## Two-sample t-test with FDR correction to first determine the DE genes;
  time = sum(system.time({
    pvals4genes <- sapply(1:p, function(i)  t.test(dat[[2]][i,], dat[[1]][i,], var.equal = FALSE)$p.value)
    
    pvals4genes.fdr <- p.adjust(pvals4genes, "BH")
    logFC.ALL <- rowMeans(dat[[2]]) - rowMeans(dat[[1]])
    names(logFC.ALL) <- dataset$gene_info$EntrezID
    
    res.pe.all <- pe(x = logFC.ALL, graphs = pe_paths, ref = dataset$gene_info$EntrezID,  nboot = 2000, verbose = FALSE)
    res.pe.all <- summary(res.pe.all)
    sigInd[,1] <- res.pe.all$pPert[match(path.names,rownames(res.pe.all))]
  })[1:2])
  
  out.obj$comp_time[,1] <- time 
  
  logFC.th <- logFC.ALL[which(pvals4genes.fdr<0.05)]
  setwd(OutFolder)
  if (length(logFC.th)>0){
    time = sum(system.time({
      # Use DE genes
      out.obj$spia <- 1
      out.obj$pe <- 1
      
      res.pe <- pe(x = logFC.th, graphs = pe_paths, ref = dataset$gene_info$EntrezID,  nboot = 2000, verbose = FALSE)
      res.pe <- summary(res.pe)
      sigInd[,2] <- res.pe$pComb[match(path.names,rownames(res.pe))]
      
    })[1:2])
    out.obj$comp_time[,2] <- time 
    
    time = sum(system.time({
      cat('... Running SPIA ....\n')
      res.spia <- runSPIA(de=logFC.th, all=dataset$gene_info$EntrezID, "BCA_KEGG")
      sigInd[,3] <- res.spia$pG[match(path.names,res.spia$Name)]
      
    })[1:2])
    out.obj$comp_time[,3] <- time 
    
  }
  setwd(filepath)
  
  ## ******************************
  ## (2) PRS----
  ## ******************************
  cat('... Running PRS ....\n')
  set.seed(iters)
  
  if (length(logFC.th)>0){
    time = sum(system.time({
      out.obj$prs <- 1
      res.prs <- prs(logFC.th, all=dataset$gene_info$EntrezID, pwys = pathSelected, nperm=2000)
      sigInd[,4] <- res.prs$p.value[match(path.names,rownames(res.prs))]
      
    })[1:2])
    out.obj$comp_time[,4] <- time 
    
  }
  
  ## ******************************
  ## (3) CePa----
  ## ******************************
  cat('... Running CePa ORA...\n')
  set.seed(iters)
  
  # Note cepa DE gene list is named by gene symbols, and the background genes are also named by gene symbols
  names(logFC.ALL) <- dataset$gene_info$symbol
  logFC.th <- logFC.ALL[which(pvals4genes.fdr<0.05)]
  if (length(logFC.th)>0){
    time = sum(system.time({
      out.obj$cepa <- 1
      res.cepa.ORA <- cepa.all(dif=names(logFC.th), bk=dataset$gene_info$symbol, pc = cepa_paths$pathways)
      
    })[1:2]) 
    out.obj$comp_time[,c(5,6,7,8,9,10)] <- time 
    
    sigInd[,5] <- sapply(res.cepa.ORA, function(a) a$equal.weight$p.value)
    sigInd[,6] <- sapply(res.cepa.ORA, function(a) a$in.degree$p.value)
    sigInd[,7] <- sapply(res.cepa.ORA, function(a) a$out.degree$p.value)
    sigInd[,8] <- sapply(res.cepa.ORA, function(a) a$betweenness$p.value)
    sigInd[,9] <- sapply(res.cepa.ORA, function(a) a$in.reach$p.value)
    sigInd[,10] <- sapply(res.cepa.ORA, function(a) a$out.reach$p.value)
  }
  
  cat('... Running CePa GSA ... \n')
  classlabels <- list()
  classlabels$label <- c(rep("Pos", n[1]), rep("Neg", n[2]))
  classlabels$treatment <- "Pos"
  classlabels$control <- "Neg"
  
  time = sum(system.time({
    res.cepa.GSA <- cepa.all(mat = cbind(dat[[1]], dat[[2]]), label = classlabels, pc = cepa_paths$pathways)
    
  })[1:2])    
  out.obj$comp_time[,c(11, 12, 13, 14, 15, 16)] <- time 
  
  sigInd[,11] <- sapply(res.cepa.GSA, function(a) a$equal.weight$p.value)
  sigInd[,12] <- sapply(res.cepa.GSA, function(a) a$in.degree$p.value)
  sigInd[,13] <- sapply(res.cepa.GSA, function(a) a$out.degree$p.value)
  sigInd[,14] <- sapply(res.cepa.GSA, function(a) a$betweenness$p.value)
  sigInd[,15] <- sapply(res.cepa.GSA, function(a) a$in.reach$p.value)
  sigInd[,16] <- sapply(res.cepa.GSA, function(a) a$out.reach$p.value)
  
  ##*******************************  
  ## (4) NetGSA, DEGraph, topologyGSA----
  ##*******************************
  cat('... Running NetGSA, DEGraph and topologyGSA ....\n')
  set.seed(iters)
  
  classlabels <- c(rep(1, n[1]), rep(2, n[2]))
  
  cat('... per pathway analysis...\n')
  
  time_est <- matrix(NA, nrow=npath, ncol=4, dimnames = list(NULL, c('2REHE', '2REML', 'REHE', 'REML')))
  time_analy <- matrix(NA, nrow=npath, ncol=6, dimnames = list(NULL, c('2REHE', '2REML', 'REHE', 'REML', 'DEGraph', 'topo')))
  
  for (loop_path in seq(1,npath)){
    set.seed(iters)
    cat("current path ...", loop_path, '/', npath, '...\n')
    index = camera_paths$membership[[loop_path]]
    pp <- length(index)
    current_data <- lapply(1:ncond, function(a) dat[[a]][index,])
    oneMat <- netgsa_paths$adjacency[[loop_path]] # 0-1 adjacency matrix, equivalent to provide the database edgelists 
    colnames(oneMat) <- rownames(oneMat) <- dataset$gene_info$EntrezID[match(rownames(oneMat), dataset$gene_info$symbol, nomatch = NA)]
    
    
    ## for each pathway, need to subset the edgelist database
    pway_database_edges <- data.table::copy(database_search)
    pway_genes <- rownames(current_data[[1]])
    pway_database_edges$edgelist <- pway_database_edges$edgelist[(paste0(base_id_src, ":", base_gene_src) %in% pway_genes) 
                                                                 & (paste0(base_id_dest, ":", base_gene_dest) %in% pway_genes)]
    pway_genes_not_in_dbs <- paste0(names(pway_database_edges$genes_not_in_dbs), ":", pway_database_edges$genes_not_in_dbs) %in% pway_genes
    pway_database_edges$genes_not_in_dbs <- if(any(pway_genes_not_in_dbs)){pway_database_edges$genes_not_in_dbs[pway_genes_not_in_dbs] } else{ NULL}
    
    ## initialize with NA, if the pathway fails to analyze will return results NA
    sigInd[loop_path,c("NetGSA2_perpath_REHE", "NetGSA2_perpath_REML",
                       "NetGSA_perpath_REHE", "NetGSA_perpath_REML")] <- NA
    
    ## use the netgsa package function prepareAdjMat to compute the partial correlation matrix
    time_est[loop_path,c('2REHE', '2REML')] = sum(system.time({
      network_info_path = prepareAdjMat(x= do.call(cbind,current_data), 
                                        group=classlabels,
                                        databases = pway_database_edges, 
                                        cluster = F,
                                        lambda_c = lambda_for_noncluster,
                                        penalize_diag = FALSE) # supply the lambda sequence, predetermined
      
    })[1:2])    
    
      
    ## for this provided lambda sequence, the estimated network information could be entirely zero; use hotelling test if this happens
    if(all(sapply(1:ncond, function(k)sum(abs(bdiag(network_info_path$Adj[[k]])))==0))){
      notes[loop_path,] <- paste(notes[loop_path,], 'NetGSA: seperate network estimation yielded zero network; use Hotelling.')
      time = sum(system.time({
        ds1 <- current_data[[1]]
        ds2 <- current_data[[2]]
        hotel <- Hotelling::hotelling.test(x = t(ds1),
                                           y = t(ds2))
        results <- data.frame(pathway = path.names[loop_path], pSize = length(pway_genes), pval = hotel$pval, pFdr = hotel$pval, teststat = hotel$stats$statistic*hotel$stats$m)
        # betas <- list(matrix(rowMeans(ds1), ncol = 1), matrix(rowMeans(ds2), ncol = 1))
        # out_cluster_F <- list(results = results, beta = betas, "s2.epsilon" = NA, "s2.gamma" = NA)
      })[1:2])

      time_analy[loop_path,c('2REHE', '2REML')] <- time
      sigInd[loop_path, c('NetGSA2_perpath_REHE', 'NetGSA2_perpath_REML')] <- results$pval
        
    }else{
      ## Run NetGSA, for each pathway
      B <- matrix(rep(1,pp), nrow=1)
      colnames(B) <- rownames(oneMat)
      rownames(B) <- path.names[loop_path]
      
      
      ## prepare the data matrix and the pathway matrix to match the order of Network matrix (to avoid error in NetGSA)
      A_mat_name <- lapply(network_info_path[["Adj"]][1:ncond], function(a) {
        do.call(c, lapply(a, rownames))
      })
      
      netgsa_data = do.call(cbind,current_data)
      netgsa_data = netgsa_data[A_mat_name[[1]],]
      B = B[,A_mat_name[[1]], drop=F]
      
      
      ## run NetGSA, based on separate group networks
      ## run for both REML and REHE
      time = sum(system.time({
        out <- NetGSA(network_info_path[["Adj"]], 
                      x=netgsa_data, 
                      group=classlabels, 
                      pathways = B,
                      lklMethod = "REHE", sampling = F, sample_n = 0.1, sample_p = 0.1 )
        
      })[1:2])
      cat('REHE s2e, s2g:', out$s2.e, out$s2.g, '\n')
      time_analy[loop_path, '2REHE'] <- time
      sigInd[loop_path,"NetGSA2_perpath_REHE"] <- out$results$pval
      
      time = sum(system.time({
        out <- NetGSA(network_info_path[["Adj"]], 
                      x=netgsa_data, 
                      group=classlabels, 
                      pathways = B,
                      lklMethod = "REML")
        
      })[1:2])
      cat('REML s2e, s2g:', out$s2.e, out$s2.g, '\n')
      
      time_analy[loop_path, '2REML'] <- time
      sigInd[loop_path,'NetGSA2_perpath_REML'] <- out$results$pval
      
      
    }

      

    if(F){
      
    # Also run NetGSA assuming the network is shared between the two conditions
    time_est[loop_path, c('REHE', 'REML')] = sum(system.time({
      network_info_path_share = prepareAdjMat(x= do.call(cbind,current_data), 
                                              group=rep(1, sum(n)),
                                              databases =  pway_database_edges, 
                                              cluster = F,
                                              lambda_c = lambda_for_noncluster,
                                              penalize_diag = FALSE)
      
    })[1:2])
    
    ## if the network is zero, use Hotelling test
    if(sum(abs(bdiag(network_info_path_share$Adj[[1]])))==0){
      notes[loop_path,] <- paste(notes[loop_path,], 'NetGSA: shared network estimation yielded zero network; use Hotelling.')
      time = sum(system.time({
        ds1 <- current_data[[1]]
        ds2 <- current_data[[2]]
        hotel <- Hotelling::hotelling.test(x = t(ds1),
                                           y = t(ds2))
        results <- data.frame(pathway = path.names[loop_path], pSize = length(pway_genes), pval = hotel$pval, pFdr = hotel$pval, teststat = hotel$stats$statistic*hotel$stats$m)
        # betas <- list(matrix(rowMeans(ds1), ncol = 1), matrix(rowMeans(ds2), ncol = 1))
        # out_cluster_F <- list(results = results, beta = betas, "s2.epsilon" = NA, "s2.gamma" = NA)
      })[1:2])
      
      time_analy[loop_path,c('REHE', 'REML')] <- time
      sigInd[loop_path, c('NetGSA_perpath_REHE', 'NetGSA_perpath_REML')] <- results$pval
      
    }else{
      
      A_mat_name <- lapply(network_info_path_share[["Adj"]][1], function(a) {
        do.call(c, lapply(a, rownames))
      })
      
      netgsa_data = do.call(cbind,current_data)
      netgsa_data = netgsa_data[A_mat_name[[1]],]
      B = B[,A_mat_name[[1]], drop=F]
      
      # make sure the structure of A is correct here
      time = sum(system.time({
        out <- NetGSA(A = list(network_info_path_share[["Adj"]][[1]],
                               network_info_path_share[["Adj"]][[1]]), 
                      x = netgsa_data, 
                      group= classlabels, 
                      pathways = B, 
                      lklMethod = "REHE", sampling = F, sample_n = 0.1, sample_p = 0.1 )
        
      })[1:2])
      time_analy[loop_path, 'REHE'] <- time
      sigInd[loop_path,"NetGSA_perpath_REHE"] <- out$results$pval
      
      time = sum(system.time({
        out <- NetGSA(A = list(network_info_path_share[["Adj"]][[1]],
                               network_info_path_share[["Adj"]][[1]]), 
                      x = netgsa_data, 
                      group= classlabels, 
                      pathways = B, 
                      lklMethod = "REML", sampling=F)
        
      })[1:2])
      time_analy[loop_path, 'REML'] <- time
      sigInd[loop_path,"NetGSA_perpath_REML"] <- out$results$pval
    }
    
    }
    
    
    time = sum(system.time({
      sigInd[loop_path,18] <- graph.T2.test(X1=t(current_data[[1]]), 
                                            X2=t(current_data[[2]]), 
                                            G=deGraph_paths$pathways[[loop_path]])$p.value
    })[1:2])
    time_analy[loop_path, 'DEGraph'] <- time
    
    ## topoGSA
    topoGSA_dat1 <- t(current_data[[1]])
    colnames(topoGSA_dat1) <- dataset$gene_info$symbol[match(colnames(topoGSA_dat1),dataset$gene_info$EntrezID, nomatch = 0)]
    topoGSA_dat2 <- t(current_data[[2]])
    colnames(topoGSA_dat2) <- dataset$gene_info$symbol[match(colnames(topoGSA_dat2),dataset$gene_info$EntrezID, nomatch = 0)]
    time = sum(system.time({
      ## when the graph is not DAG, the result is null
      test <- tryCatch(pathway.mean.test(y1=topoGSA_dat1, 
                                         y2=topoGSA_dat2, 
                                         dag = topoGSA_paths$pathways[[loop_path]], 
                                         alpha = 0.05, 
                                         perm.num = 1000), error=function(e) NULL)

    })[1:2])
    time_analy[loop_path, 'topo'] <- time
    
    if (is.null(test)){
      sigInd[loop_path,19] <- NA; sigInd[loop_path,20]  <- NA
    } else {
      sigInd[loop_path,19] <- test$p.value
      sigInd[loop_path,20] <- test$p.value.var
    }
  }
  
  tmp1 = colSums(time_est)
  
  tmp2 = colSums(time_analy)
  
  out.obj$comp_time[,c("NetGSA2_perpath_REHE_net",
                       "NetGSA2_perpath_REML_net",
                       "NetGSA_perpath_REHE_net",
                       "NetGSA_perpath_REML_net")] <- tmp1
  
  out.obj$comp_time[,c("NetGSA2_perpath_REHE_analy",
                       "NetGSA2_perpath_REML_analy",
                       "NetGSA_perpath_REHE_analy",
                       "NetGSA_perpath_REML_analy",
                       "DEGraph",
                       "topologyGSA")] <- tmp2 
  
  
  
  
  
  cat('... all pathway analysis...\n')
  set.seed(iters)
  
  ## we run version of cluster and noncluster, for REHE only; run for separate group network and shared group network
  
  ## separate group, noncluster
  time = sum(system.time({
    ## use the netgsa package function prepareAdjMat to compute the partial correlation matrix
    network_info_all_ncl = prepareAdjMat(x= do.call(cbind,dat), 
                                         group=classlabels,
                                         databases = database_search, 
                                         cluster = F,
                                         lambda_c = lambda_for_noncluster,
                                         penalize_diag = FALSE) # supply the lambda sequence, predetermined
    
  })[1:2])
  
  out.obj$comp_time[,'NetGSA2_allpath_ncl_REHE_net'] <-time
  
  if(all(sapply(1:ncond, function(k)sum(abs(bdiag(network_info_all_ncl$Adj[[k]])))==0))){
    notes[npath+1,] <- paste(notes[npath+1,], 'NetGSA: separate network nonclustering estimation yielded zero network; use Hotelling.')
    time = sum(system.time({
      ds1 <- dat[[1]]
      ds2 <- dat[[2]]
      hotel <- Hotelling::hotelling.test(x = t(ds1),
                                         y = t(ds2))
      results <- data.frame(pathway = path.names, pval = hotel$pval)
    })[1:2])
    
    out.obj$comp_time[,'NetGSA2_allpath_ncl_REHE_analy'] <- time
    sigInd[, 'NetGSA2_allpath_ncl_REHE'] <- results$pval
    
  }else{
    
    ## Run NetGSA, for all pathways
    B <- matrix(0, nrow=npath, ncol = nrow(dat[[1]]))
    rownames(B) <- names(camera_paths$membership)
    colnames(B) <- rownames(dat[[1]])
    for(i in 1:length(camera_paths$membership)) B[i,camera_paths$membership[[i]]] <- 1
    
    
    ## run NetGSA, based on separate group networks
    ## run for REHE only; run both sampling and nonsampling version
    time = sum(system.time({
      out <- NetGSA(network_info_all_ncl[["Adj"]], 
                    x= do.call(cbind,dat), 
                    group=classlabels, 
                    pathways = B,
                    lklMethod = "REHE", sampling = F, sample_n = 0.1, sample_p = 0.1 )
    })[1:2])
    
    out.obj$comp_time[,'NetGSA2_allpath_ncl_REHE_analy'] <- time
    
    ######
    sigInd[,"NetGSA2_allpath_ncl_REHE"] <- out$results[match(rownames(sigInd),out$results$pathway),"pval"]
    
    #--
    time = sum(system.time({
      out <- NetGSA(network_info_all_ncl[["Adj"]], 
                    x= do.call(cbind,dat), 
                    group=classlabels, 
                    pathways = B,
                    lklMethod = "REHE", sampling = T, sample_n = 0.1, sample_p = 0.1 )
    })[1:2])
    
    out.obj$comp_time[,'NetGSA2_allpath_ncl_reREHE_analy'] <- time
    
    sigInd[,"NetGSA2_allpath_ncl_reREHE"] <- out$results[match(rownames(sigInd),out$results$pathway),"pval"]
    
    
  }
  
  ## separate group, cluster
  out.obj$comp_time[,'NetGSA2_allpath_cl_REHE_net'] = sum(system.time({
    network_info_all_cl = prepareAdjMat(x= do.call(cbind,dat), 
                                        group=classlabels,
                                        databases = database_search, 
                                        cluster = T,
                                        lambda_c = lambda_for_cluster,
                                        penalize_diag = FALSE) # supply the lambda sequence, predetermined
    
  })[1:2])
    
  if(all(sapply(1:ncond, function(k)sum(abs(bdiag(network_info_all_cl$Adj[[k]])))==0))){
    notes[npath+1,] <- paste(notes[npath+1,], 'NetGSA: separate network clustering estimation yielded zero network; use Hotelling.')
    time = sum(system.time({
      ds1 <- dat[[1]]
      ds2 <- dat[[2]]
      hotel <- Hotelling::hotelling.test(x = t(ds1),
                                         y = t(ds2))
      results <- data.frame(pathway = path.names, pval = hotel$pval)
    })[1:2])
    
    out.obj$comp_time[,'NetGSA2_allpath_cl_REHE_analy'] <- time
    sigInd[, 'NetGSA2_allpath_cl_REHE'] <- results$pval
  }else{
    ## Run NetGSA, for all pathways
    B <- matrix(0, nrow=npath, ncol = nrow(dat[[1]]))
    rownames(B) <- names(camera_paths$membership)
    colnames(B) <- rownames(dat[[1]])
    for(i in 1:length(camera_paths$membership)) B[i,camera_paths$membership[[i]]] <- 1
    
    
    ## run NetGSA, based on seperate group networks
    ## run for REHE only
    out.obj$comp_time[,'NetGSA2_allpath_cl_REHE_analy'] <- 
      sum(system.time({
      out <- NetGSA(network_info_all_cl[["Adj"]], 
                    x= do.call(cbind,dat), 
                    group=classlabels, 
                    pathways = B,
                    lklMethod = "REHE", sampling = F, sample_n = 0.1, sample_p = 0.1 )
    })[1:2])
    
    sigInd[,"NetGSA2_allpath_cl_REHE"] <- out$results[match(rownames(sigInd),out$results$pathway),"pval"]
    
    #-- reREHE for sampling as well
    out.obj$comp_time[,'NetGSA2_allpath_cl_reREHE_analy'] <- 
      sum(system.time({
        out <- NetGSA(network_info_all_cl[["Adj"]], 
                      x= do.call(cbind,dat), 
                      group=classlabels, 
                      pathways = B,
                      lklMethod = "REHE", sampling = T, sample_n = 0.1, sample_p = 0.1 )
      })[1:2])
    
    sigInd[,"NetGSA2_allpath_cl_reREHE"] <- out$results[match(rownames(sigInd),out$results$pathway),"pval"]
    

  }

  if(F){
    
  ## shared group, noncluster
  out.obj$comp_time[,'NetGSA_allpath_ncl_REHE_net'] <- sum(system.time({
    ## use the netgsa package function prepareAdjMat to compute the partial correlation matrix
    network_info_all_ncl = prepareAdjMat(x= do.call(cbind,dat), 
                                         group= rep(1, sum(n)),
                                         databases = database_search, 
                                         cluster = F,
                                         lambda_c = lambda_for_noncluster,
                                         penalize_diag = FALSE) # supply the lambda sequence, predetermined
    
  })[1:2])
   
  if(sum(abs(bdiag(network_info_all_ncl$Adj[[1]])))==0){
    notes[npath+1,] <- paste(notes[npath+1,], 'NetGSA: shared network nonclustering estimation yielded zero network; use Hotelling.')
    time = sum(system.time({
      ds1 <- dat[[1]]
      ds2 <- dat[[2]]
      hotel <- Hotelling::hotelling.test(x = t(ds1),
                                         y = t(ds2))
      results <- data.frame(pathway = path.names, pval = hotel$pval)
    })[1:2])
    
    out.obj$comp_time[,'NetGSA_allpath_ncl_REHE_analy'] <- time
    sigInd[, 'NetGSA_allpath_ncl_REHE'] <- results$pval
  }else{
    ## Run NetGSA, for all pathways
    B <- matrix(0, nrow=npath, ncol = nrow(dat[[1]]))
    rownames(B) <- names(camera_paths$membership)
    colnames(B) <- rownames(dat[[1]])
    for(i in 1:length(camera_paths$membership)) B[i,camera_paths$membership[[i]]] <- 1
    
    
    ## run NetGSA, based on shared group networks
    ## run for REHE only
    out.obj$comp_time[,'NetGSA_allpath_ncl_REHE_analy'] <-sum(system.time({
      out <- NetGSA(list(network_info_all_ncl[["Adj"]][[1]], 
                         network_info_all_ncl[["Adj"]][[1]]),
                    x= do.call(cbind,dat), 
                    group=classlabels, 
                    pathways = B,
                    lklMethod = "REHE", sampling = T, sample_n = 0.1, sample_p = 0.1 )
    })[1:2])
    
    sigInd[,"NetGSA_allpath_ncl_REHE"] <- out$results$pval
    
    
  }
    
  
  
  ## shared group, cluster
  out.obj$comp_time[,'NetGSA_allpath_cl_REHE_net'] <- sum(system.time({
    ## use the netgsa package function prepareAdjMat to compute the partial correlation matrix
    network_info_all_cl = prepareAdjMat(x= do.call(cbind,dat), 
                                        group= rep(1, sum(n)),
                                        databases = database_search, 
                                        cluster = T,
                                        lambda_c = lambda_for_cluster,
                                        penalize_diag = FALSE) # supply the lambda sequence, predetermined
    
    
  })[1:2])

  if(sum(abs(bdiag(network_info_all_cl$Adj[[1]])))==0){
    notes[npath+1,] <- paste(notes[npath+1,], 'NetGSA: shared network clustering estimation yielded zero network; use Hotelling.')
    time = sum(system.time({
      ds1 <- dat[[1]]
      ds2 <- dat[[2]]
      hotel <- Hotelling::hotelling.test(x = t(ds1),
                                         y = t(ds2))
      results <- data.frame(pathway = path.names, pval = hotel$pval)
    })[1:2])
    
    out.obj$comp_time[,'NetGSA_allpath_cl_REHE_analy'] <- time
    sigInd[, 'NetGSA_allpath_cl_REHE'] <- results$pval
  }else{
    ## Run NetGSA, for all pathways
    B <- matrix(0, nrow=npath, ncol = nrow(dat[[1]]))
    rownames(B) <- names(camera_paths$membership)
    colnames(B) <- rownames(dat[[1]])
    for(i in 1:length(camera_paths$membership)) B[i,camera_paths$membership[[i]]] <- 1
    
    
    ## run NetGSA, based on separate group networks
    ## run for REHE only
    out.obj$comp_time[,'NetGSA_allpath_cl_REHE_analy'] <- sum(system.time({
      out <- NetGSA(list(network_info_all_cl[["Adj"]][[1]], 
                         network_info_all_cl[["Adj"]][[1]]),
                    x= do.call(cbind,dat), 
                    group=classlabels, 
                    pathways = B,
                    lklMethod = "REHE", sampling = T, sample_n = 0.1, sample_p = 0.1 )
      
    })[1:2])
    
    sigInd[,"NetGSA_allpath_cl_REHE"] <- out$results$pval
    
  }

  }
  
  
  ##*******************************  
  ## (5) camera----
  ##*******************************
  cat('... Running CAMERA ....\n')
  set.seed(iters)
  
  classBca2012 <- c(rep(1, n[1]), rep(0, n[2]))
  names(classBca2012) <- dataset$sample_info$GSMID
  expBca2012 <- cbind(dat[[1]], dat[[2]])
  rownames(expBca2012) <- dataset$gene_info$EntrezID
  design <- cbind(Intercept = 1, Group = classBca2012)
  
  out.obj$comp_time[,21] <- sum(system.time({
    res.camera <- limma::camera(expBca2012, camera_paths$membership, design) 
    
  })[1:2])
  sigInd[,21] <- res.camera$PValue[match(path.names,rownames(res.camera))]
  
  ## ******************************
  ## (6) PathNet----
  ## ******************************
  cat('... Running PathNet ....\n')
  set.seed(iters)
  
  out.obj$comp_time[,22] <- sum(system.time({
    logPval <- -log10(pvals4genes)
    dat4PathNet <- cbind(as.numeric(gsub("ENTREZID:","",dataset$gene_info$EntrezID)),logPval)
    colnames(dat4PathNet) <- c("Gene.ID", "ERneg")
    res.PathNet <- PathNet(Enrichment_Analysis = TRUE, Contextual_Analysis = FALSE, 
                           DirectEvidence_info = dat4PathNet, Column_DirectEvidence = 2, 
                           Adjacency = pathNet_paths$adjacency, pathway = pathNet_paths$pathways, 
                           n_perm = 2000, threshold = 0.05, use_sig_pathways  = FALSE)
    
    res.PathNet <- res.PathNet$enrichment_results
    dict.pathnet <- data.frame("full.name" = sort(path.names), 
                               'path.name' = sort(as.character(res.PathNet$Name)), stringsAsFactors = FALSE)
    res.PathNet$Name <- dict.pathnet$full.name[match(as.character(res.PathNet$Name),dict.pathnet$path.name)]
    
  })[1:2])
  sigInd[,22] <- res.PathNet$p_PathNet[match(path.names, res.PathNet$Name)]
  
  out.obj$sigInd <- sigInd
  out.obj$notes <- notes
  
  return(out.obj)
}


res = one_iteration(iters)

save('res', file=paste0(OutFolder, newFileName, '.rda'))
