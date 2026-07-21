conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
if (dir.exists(conda_lib)) .libPaths(conda_lib)
invisible(suppressPackageStartupMessages(lapply(
  c("quadprog", "reshape", "e1071", "Seurat", "ROCR", "varhandle", "MAST"),
  library,
  character.only = TRUE
)))

#Deconvolution Functions

#trim bulk and single-cell data to contain the same genes
trimData<-function(Signature,bulkData){
	Genes<-intersect(rownames(Signature),names(bulkData))
  	B<-bulkData[Genes]
  	S<-Signature[Genes,]
  	return(list("sig"=S,"bulk"=B))
}


#solve using OLS, constrained such that cell type numbers>0
solveOLS<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D, "2")
  solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#return cell number, not proportion
#do not print output
solveOLSInternal<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D, "2")
  solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS<-function(S,B){
  #first solve OLS, use this solution to find a starting point for the weights
  solution<-solveOLSInternal(S,B)
  #now use dampened WLS, iterate weights until convergence
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-findDampeningConstant(S,B,solution)
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    #decrease step size for convergence
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#solve WLS given a dampening constant
solveDampenedWLSj<-function(S,B,goldStandard,j){
  multiplier<-1*2^(j-1)
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsDampened<-wsScaled
  wsDampened[which(wsScaled>multiplier)]<-multiplier
  W<-diag(wsDampened)
  D<-t(S)%*%W%*%S
  d<- t(S)%*%W%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#find a dampening constant for the weights using cross-validation
findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-cbind(solutions,sol)
    }
    solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j<-which.min(colMeans(solutionsSd^2))
  return(j)
}

solveSVR<-function(S,B){
  #scaling
  ub=max(c(as.vector(S),B)) #upper bound
  lb=min(c(as.vector(S),B)) #lower bound
  Bs=((B-lb)/ub)*2-1
  Ss=((S-lb)/ub)*2-1
  
  #perform SVR
  model<-svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",kernel ="linear",cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef<0)]<-0
  coef<-as.vector(coef)
  names(coef)<-colnames(S)
  print(round(coef/sum(coef),5))
  return(coef/sum(coef))
}

#perform DE analysis using Seurat
DEAnalysis<-function(scdata,id,path){
  exprObj<-CreateSeuratObject(raw.data=as.data.frame(scdata), project = "DE")
  exprObj2<-SetIdent(exprObj,ident.use=as.vector(id))
  print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- FindMarkers(object=exprObj2, ident.1 = i, ident.2 = NULL, 
                             only.pos = TRUE, test.use = "bimod")
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
  }
}

#build signature matrix using genes identified by DEAnalysis()
buildSignatureMatrixUsingSeurat<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){
  
  #perform differential expression analysis
  DEAnalysis(scdata,id,path)
  
  numberofGenes<-c()
  for (i in unique(id)){
    load(file=paste(path,"/de_",i,".RData",sep=""))
    DEGenes<-rownames(de_group)[intersect(which(de_group$p_val_adj<pval.cutoff),which(de_group$avg_logFC>diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)
    assign(paste("cluster_lrTest.table.",i,sep=""),de_group[which(rownames(de_group)%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }
  
  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 2:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
        Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1) #G is optimal gene number
  #
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
      Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)
  save(Sig,file=paste(path,"/Sig.RData",sep=""))
  return(Sig)
}

##alternative differential expression method using MAST

#functions for DE

Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}

stat.log2=function(data.m, group.v, pseudo.count){
  #data.m=data.used.log2
  log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = varhandle::unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}

v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
m.auc=function(data.m,group.v) {
  AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
  AUC[is.na(AUC)]=0.5
  return(AUC)
  
}  

#perform DE analysis using MAST	    
DEAnalysisMAST<-function(scdata,id,path){
  
  pseudo.count = 0.1
  data.used.log2   <- log2(scdata+pseudo.count)
  colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
  diff.cutoff=0.5
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (i in unique(id)){
    cells.symbol.list2     = colnames(data.used.log2)[which(id==i)]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))                          
    cells.symbol.list1     = colnames(data.used.log2)[which(id != i)]
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))   
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    DE <- bigtable[bigtable$log2_fc >diff.cutoff,] 
    dim(DE)
    if(dim(DE)[1]>1){
      data.1                 = data.used.log2[,cells.coord.list1]
      data.2                 = data.used.log2[,cells.coord.list2]
      genes.list = rownames(DE)
      log2fold_change        = cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) = c("gene.name", "log2fold_change")
      counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
      groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(i, length(cells.coord.list2) ) )
      groups  = as.character(groups)
      data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
      colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                            primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                            geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                            phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      # .3 MAST 
      head(colData(vbeta.1))
      zlm.output <- zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC=TRUE)
      zlm.lr <- lrTest(zlm.output, 'Population')
      zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
      
      
      
      lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",i,sep=""), "log2fold_change", "Auc")
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]
      
      #. 4 save results
      write.csv(cluster_lrTest.table, file=paste(path,"/",i,"_lrTest.csv", sep=""))
      save(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.RData", sep=""))
    }
  }
}

#build signature matrix using genes identified by DEAnalysisMAST()
buildSignatureMatrixMAST<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){
  #compute differentially expressed genes for each cell type
  DEAnalysisMAST(scdata,id,path)
  
  #for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the estimated fold-change
  #is greater than 0.5
  numberofGenes<-c()
  for (i in unique(id)){
    if(file.exists(paste(path,"/",i,"_MIST.RData", sep=""))){
      load(file=paste(path,"/",i,"_MIST.RData", sep=""))
      pvalue_adjusted<-p.adjust(cluster_lrTest.table[,3], method = "fdr", n = length(cluster_lrTest.table[,3]))
      cluster_lrTest.table<-cbind(cluster_lrTest.table,pvalue_adjusted)
      DEGenes<-cluster_lrTest.table$Gene[intersect(which(pvalue_adjusted<pval.cutoff),which(cluster_lrTest.table$log2fold_change>diff.cutoff))]
      nonMir = grep("MIR|Mir", DEGenes, invert = T)  # because Mir gene is usually not accurate 
      assign(paste("cluster_lrTest.table.",i,sep=""),cluster_lrTest.table[which(cluster_lrTest.table$Gene%in%DEGenes[nonMir]),])
      numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
    }
  }
  
  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #for each, iterate and choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 2:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
      Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1)
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
    temp<-paste("cluster_lrTest.table.",i,sep="")
    temp<-as.name(temp)
    temp<-eval(parse(text = temp))
    temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
    Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)
  save(Sig,file=paste(path,"/Sig.RData",sep=""))
  return(Sig)
}

# 加载必要库 ----------------------------------------------------------------
library(doParallel)
library(foreach)
library(ROCR)  # 假设prediction函数来自ROCR包
library(reshape)

# # 加载自定义函数（需替换为实际路径）-----------------------------------------
# source("/PATH/TO/YOUR/buildSignatureMatrixMAST.R")
# source("/PATH/TO/YOUR/solveDampenedWLS.R")

# 串行初始化配置 -----------------------------------------------------------
root_path <- Sys.getenv(
    "ROOT_PATH",
    unset = ""
)
if (!nzchar(root_path)) stop("Set ROOT_PATH to a dataset directory containing tissue folders.")
dwls_recursive <- !identical(tolower(Sys.getenv("DWLS_RECURSIVE", unset = "1")), "0")
required_files <- c("bulk.txt", "singleCellExpr.txt", "singleCellLabels.txt", "signature.txt")
folders <- list.dirs(root_path, full.names = TRUE, recursive = dwls_recursive)
folders <- folders[sapply(folders, function(folder) {
  all(file.exists(file.path(folder, required_files)))
})]
folders <- sort(folders)

if (length(folders) == 0) {
  cat("\n==== 处理统计 ====\n")
  print(table(character(0)))
  writeLines(character(0), file.path(root_path, "processing_log.txt"))
  quit(save = "no", status = 0)
}

# 串行处理主循环 -----------------------------------------------------------
result_log <- foreach(folder = folders, 
                     .packages = c("ROCR", "utils", 'quadprog', 'reshape', 'e1071', 'Seurat', 'varhandle', 'MAST', 'foreach'),  # 确保加载依赖包
                     .combine = rbind,
                     .errorhandling = "pass") %do% {

  data_dir <- folder
  if (!dir.exists(data_dir)) {
    return(paste("[Skipped]", folder, "缺少data目录"))
  }

  # 检查结果文件
  result_file <- file.path(folder, "result", "DWLS.txt")
  if (file.exists(result_file)) {
    print(paste("[Skipped]", folder, "DWLS结果已存在"))
    return(paste("[Skipped]", folder, "DWLS结果已存在"))
  }
  
  # 主处理流程
  tryCatch({
    t_folder_start <- Sys.time()
    print("\n==== DWLS START ====")
    print(paste("[folder]", folder))
    print(paste("[time] start", format(t_folder_start, "%Y-%m-%d %H:%M:%S")))
    print(paste("[bulk_path]", file.path(data_dir, "bulk.txt")))
    print(paste("[sc_expr_path]", file.path(data_dir, "singleCellExpr.txt")))
    print(paste("[sc_label_path]", file.path(data_dir, "singleCellLabels.txt")))

    # 加载数据（避免使用setwd）
    t_io_start <- Sys.time()
    bulk_path <- file.path(data_dir, "bulk.txt")
    sc_expr_path <- file.path(data_dir, "singleCellExpr.txt")
    sc_label_path <- file.path(data_dir, "singleCellLabels.txt")
    
    if (!file.exists(bulk_path) || !file.exists(sc_expr_path) || !file.exists(sc_label_path)) {
      print("[Skipped] 缺少必要的输入文件")
      return(paste("[Skipped]", folder, "缺少必要的输入文件"))
    }
    
    bulk <- read.delim(bulk_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    sc_m <- read.delim(sc_expr_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    sc_ct_raw <- read.delim(sc_label_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    t_io_end <- Sys.time()
    print(paste("[time] IO", round(as.numeric(difftime(t_io_end, t_io_start, units = "secs")), 2), "s"))

    print(paste("[bulk dim (gene x sample)]", paste(dim(bulk), collapse = " x ")))
    print(paste("[bulk colnames head]", paste(head(colnames(bulk)), collapse = ", ")))
    print(paste("[bulk rownames head]", paste(head(rownames(bulk)), collapse = ", ")))
    print(paste("[sc_m dim (gene x cell)]", paste(dim(sc_m), collapse = " x ")))
    print(paste("[sc_m colnames head]", paste(head(colnames(sc_m)), collapse = ", ")))
    print(paste("[sc_m rownames head]", paste(head(rownames(sc_m)), collapse = ", ")))
    print(paste("[sc_ct_raw dim]", paste(dim(sc_ct_raw), collapse = " x ")))
    print(paste("[sc_ct_raw colnames]", paste(colnames(sc_ct_raw), collapse = ", ")))

    # 标签文件为单列，顺序与 singleCellExpr 列顺序一致
    if (nrow(sc_ct_raw) != ncol(sc_m)) {
      print("[Failed] label 与表达矩阵不匹配")
      return(paste("[Failed]", folder, "label 与表达矩阵不匹配"))
    }
    label_col <- colnames(sc_ct_raw)[1]  # 兼容不同的列名
    ct <- gsub(" ", "_", sc_ct_raw[[label_col]])
    print(paste("[labels aligned] n =", length(ct), "; unique =", length(unique(ct))))
    print(paste("[labels table top]", paste(head(sort(table(ct), decreasing = TRUE)), collapse = ", ")))
    
    # 创建中间目录
    mid_dir <- file.path(data_dir, "mid_data")
    dir.create(mid_dir, showWarnings = FALSE)
    
    # 读取已有 signature 或计算
    t_sig_start <- Sys.time()
    signature_file <- file.path(data_dir, "signature.txt")
    if (file.exists(signature_file)) {
      print(paste("[Signature] loading", signature_file))
      sig_df <- read.delim(signature_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
      Signature <- as.matrix(sig_df)
      storage.mode(Signature) <- "numeric"
    } else {
      print("[Signature] building with MAST ...")
      Signature <- buildSignatureMatrixMAST(scdata = sc_m, id = ct,
                                           path = mid_dir,
                                           diff.cutoff = 0.5, pval.cutoff = 0.05)
    }
    t_sig_end <- Sys.time()
    print(paste("[time] Signature", round(as.numeric(difftime(t_sig_end, t_sig_start, units = "secs")), 2), "s"))
    print(paste("[Signature dim]", paste(dim(Signature), collapse = " x ")))
    if (!is.null(Signature) && length(Signature) > 0) {
      print(paste("[Signature rownames head]", paste(head(rownames(Signature)), collapse = ", ")))
      print(paste("[Signature colnames]", paste(colnames(Signature), collapse = ", ")))
    }

    # 若无重叠基因则跳过
    overlap_genes <- length(intersect(rownames(Signature), rownames(bulk)))
    print(paste("[overlap genes]", overlap_genes))
    if (overlap_genes == 0) {
      print("[Skipped] 无重叠基因，无法计算")
      return(paste("[Skipped]", folder, "无重叠基因，无法计算"))
    }
    
    # 处理样本
    t_deconv_start <- Sys.time()
    all_results_matrix <- foreach(sample_id = 1:ncol(bulk), .combine = rbind) %do% {
      t_sample_start <- Sys.time()
      print(paste("[Sample]", sample_id, "/", ncol(bulk)))
      bulk_single <- setNames(as.numeric(bulk[, sample_id]), rownames(bulk))
      
      tr <- trimData(Signature, bulk_single)
      print(paste("  [trim] sig dim =", paste(dim(tr$sig), collapse = " x "), "; bulk len =", length(tr$bulk)))
      print(paste("  [trim genes head]", paste(head(rownames(tr$sig)), collapse = ", ")))
      solDWLS <- solveDampenedWLS(tr$sig, tr$bulk)
      print(paste("  [sol] sum =", round(sum(solDWLS), 4)))
      print(paste("  [time] sample", sample_id, round(as.numeric(difftime(Sys.time(), t_sample_start, units = "secs")), 2), "s"))
      data.frame(t(solDWLS)) |> 
        `colnames<-`(colnames(Signature)) |> 
        `rownames<-`(paste0("Bulk_", sample_id))
    }
    t_deconv_end <- Sys.time()
    print(paste("[time] Deconvolution", round(as.numeric(difftime(t_deconv_end, t_deconv_start, units = "secs")), 2), "s"))
    
    # 保存结果
    t_save_start <- Sys.time()
    result_folder <- file.path(folder, "result")
    dir.create(result_folder, showWarnings = FALSE)
    write.table(all_results_matrix, file.path(result_folder, "DWLS.txt"), sep = "\t")
    t_save_end <- Sys.time()
    print(paste("[time] Save", round(as.numeric(difftime(t_save_end, t_save_start, units = "secs")), 2), "s"))
    print(paste("[Processed]", folder))
    print(paste("[time] total", round(as.numeric(difftime(Sys.time(), t_folder_start, units = "secs")), 2), "s"))
    
    return(paste("[Processed]", folder))
    
  }, error = function(e) {
    print(paste("[Failed]", folder, "Error:", e$message))
    return(paste("[Failed]", folder, "Error:", e$message))
  })
}

# 打印结果统计 ------------------------------------------------------------
cat("\n==== 处理统计 ====\n")
print(table(gsub("\\[.*?\\]", "", result_log)))

# 保存详细日志
writeLines(result_log, file.path(root_path, "processing_log.txt"))
