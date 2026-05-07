library(BoolNet)
library(doParallel)
library(foreach)
set.seed(42)
# ── 1. 함수 정의 (기존과 동일) ──────────────────────────────────
Large_Attr = function(Network, startStates = 100000){
  attr = getAttractors(Network, method = "random", startStates = startStates)
  N_Num = length(Network$genes)
  matrix_elements = c()
  attr_size = c()
  for(i in 1:length(attr$attractors)){
    attr_vector = rep(0, N_Num)
    for(j in 1:dim(attr$attractors[[i]]$involvedStates)[2]){
      attr_partial_vector = c()
      for(k in 1:dim(attr$attractors[[i]]$involvedStates)[1]){
        node_num = min(32, N_Num-32*(k-1))
        attr_partial_vector = c(attr_partial_vector, rev(Binary_Vector(2^32+attr$attractors[[i]]$involvedStates[k,j],node_num)))
      }
      attr_vector = attr_vector + attr_partial_vector
    }
    attr_vector = attr_vector / dim(attr$attractors[[i]]$involvedStates)[2]
    #matrix_elements = c(matrix_elements, 1)
    for(j in 1:length(attr_vector)){
      matrix_elements = c(matrix_elements, attr_vector[j])
    }
    attr_size = c(attr_size, attr$attractors[[i]]$basinSize)
  }
  
  result = list()
  result[["Attr_Mat"]] = matrix(matrix_elements, ncol = length(attr$attractors))
  result[["Attr_Size"]] = attr_size
  return(result)
}
Binary_Vector = function(Number, Length){
  Binary = c()
  if((Number+1) > 2^Length){
    Number = Number %% (2^Length)
  }
  for(i in 1:Length){
    Quotient = Number %/% (2^(Length-i))
    Binary = c(Binary, Quotient)
    if(Quotient == 1){
      Number = Number - (2^(Length-i))
    }
  }
  return(Binary)
}
cal_node_activity = function(k, index){
  na = sum(k$Attr_Mat[index,] * k$Attr_Size) /sum(k$Attr_Size)
  #print(paste(gene_list[index],"'s Node_Activity: " ,na))
  return(na)
}


# ── 2. 병렬 환경 설정 ──────────────────────────────────────────────
# 사용할 코어 수 설정 
num_cores <- 30
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# 작업 경로 및 설정
bnet_path <- "/home/yeohs0212/MM/IQcell/PDAC_boolean_output/pdac_boolean_network.bnet"
out_dir <- "/home/yeohs0212/MM/IQcell/PDAC_boolean_output/perturbation_results"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Nominal 정보 저장
bnet = loadNetwork(bnet_path)
attr= Large_Attr(bnet,  startStates = 100000)
for (i in 1:length(bnet$genes)){
  ca = cal_node_activity(attr,i)
  cat(bnet$genes[i], ca, '\n')
}

normal_attr_df = data.frame(t(attr$Attr_Mat))
colnames(normal_attr_df) = bnet$genes
normal_attr_df['Attr_Size'] = attr$Attr_Size
write.csv(normal_attr_df,'/home/yeohs0212/MM/IQcell/PDAC_boolean_output/normal_attr_df.csv', row.names = FALSE)



# ── 3. 모든 Perturbation 조합 생성 ──────────────────────────────────
bnet <- loadNetwork(bnet_path)
gene_names <- bnet$genes
# 모든 조합 생성: (유전자이름, 고정값)
tasks <- expand.grid(gene = gene_names, val = c(0, 1))

# ── 4. 병렬 처리 루프 ─────────────────────────────────────────────────
# .packages 옵션으로 worker들이 BoolNet을 사용할 수 있게 함
foreach(i = 1:nrow(tasks), .packages = "BoolNet") %dopar% {
  gene <- tasks$gene[i]
  val <- tasks$val[i]
  
  # 개별 worker에서 네트워크 로드 (메모리 안정성)
  net <- loadNetwork(bnet_path)
  
  # Perturbation 수행
  result <- tryCatch({
    perturb_net <- fixGenes(net, gene, val)
    # initial state 개수 1,000,000으로 설정
    Large_Attr(perturb_net, startStates = 100000)
  }, error = function(e) {
    return(NULL)
  })
  
  # 결과 저장
  if (!is.null(result)) {
    df <- data.frame(t(result$Attr_Mat))
    colnames(df) <- net$genes
    df$Attr_Size <- result$Attr_Size
    
    save_path <- file.path(out_dir, paste0("perturb_", gene, "_fix_", val, ".csv"))
    write.csv(df, save_path, row.names = FALSE)
  }
}

# 병렬 클러스터 종료
stopCluster(cl)
cat("병렬 계산이 완료되었습니다.\n")
