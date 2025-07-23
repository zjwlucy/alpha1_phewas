# ------------------------------------------------------------
#  Single‑variant, single‑phecode association with GENESIS
#  • handles missing data
#  • uses an external kinship matrix
# ------------------------------------------------------------
# packageVersion("SeqArray")     # ‘1.48.0’
# packageVersion("SeqVarTools")  # ‘1.46.0’
# packageVersion("GENESIS")      # ‘2.38.0’
# packageVersion("Matrix")       # ‘1.7.3’
# packageVersion("Biobase")      # ‘2.68.0’
# packageVersion("data.table")   # ‘1.17.0’
# packageVersion("SPAtest")      # ‘3.1.2’


# ===============================
# (0). command‑line args 
# ===============================
  args <- commandArgs(trailingOnly = TRUE)
  phecode    <- args[1]            # e.g. "phecode_270.34"
  variant_id <- args[2]            # e.g. "ZZ_MM"
  input_dir  <- args[3]
  main_file  <- args[4] 
  gds_file   <- args[5]            # path to .gds created by seqSNP2GDS
  out_prefix <- args[6]
  out_dir    <- args[7]            # output directory 
  
  suppressPackageStartupMessages({
    library(SeqArray)     # open the GDS
    library(SeqVarTools)  # SeqVarData helper
    library(GENESIS)      # fitNullModel + assocTestSingle :contentReference[oaicite:0]{index=0}
    library(Matrix)       # sparse matrices for kinship
    library(Biobase)      # AnnotatedDataFrame
    library(data.table)  
    library(SPAtest)
  })
  

# ==============================================================  
#                        prepare  
# ============================================================== 
# (1) genetic file 
  gds     <- seqOpen(paste0(input_dir, gds_file), allow.duplicate=T)
  gds.ids <- seqGetData(gds, "sample.id")
  seqGetData(gds, "variant.id")

# (2) main file 
  covar_cols    <- c("age_at_cdr", "genetic_sex", paste0("PC", 1:10)) # basic covariates
  input_file    <- paste0(input_dir, main_file)
  dat           <- fread(input_file, data.table = FALSE, select = c("person_id", covar_cols, phecode))
  dat$sample.id <- dat$person_id <- as.character(dat$person_id) # SeqVarData() need sample id
  dat           <- dat[match(gds.ids, dat$person_id), ]         # puts rows in GDS order # rows with no match become NA
  rownames(dat) <- dat$sample.id                                # Biobase needs row‑names   
  sampleADF     <- AnnotatedDataFrame(dat)
#
  #identical(rownames(dat), dat$sample.id) # TRUE
  #dim(dat)   # 317207  15
  #table(dat[[phecode]]) # 296570 + 139 = 296709  
  #dat[1,]

# ---------------
  svd  <- SeqVarData(gds, sampleData = sampleADF) 
  
# keep only individuals with no missing values in the covariates (but may still have missing in SNPs)
  keep     <- complete.cases(pData(sampleData(svd))[, c(covar_cols, phecode) ])  # logical vector
  keep_ids <- rownames(sampleData(svd))[keep];  length(keep_ids)    # 296709
  seqSetFilter(svd, sample.id = keep_ids, verbose = FALSE)
# ---------------  
  
  
# (3) build a sparse kinship matrix
  rel_df <- fread("relatedness.tsv", data.table=F) # columns: i.s, j.s, kin
  K      <- makeSparseMatrix(rel_df, sample.include = keep_ids, diag.value = 1) # GENESIS default for GRM       
  dim(K)  # 296709 296709
  
  
 
# ==============================================================  
#                           Fit models 
# ==============================================================  
# Null model
  nullmod <- fitNullModel( dat[keep_ids, ],      
                           outcome     = phecode,         # binary outcome column
                           covars      = covar_cols,
                           cov.mat     = list(K),         # kinship as random effect
                           family      = "binomial"       # logistic mixed model :contentReference[oaicite:1]{index=1}
  )
  
# Association test for the single variant 
  iter  <- SeqVarBlockIterator(svd, variantBlock = 3)
  assoc <- assocTestSingle(iter, null.model = nullmod, test = "Score.SPA", BPPARAM = BiocParallel::SerialParam())  
  
  
  
########################  
# save results
  #out_file <- paste0(out_dir, out_prefix, variant_id, "_", phecode, "_assoc.RDS")
  out_file <- paste0(out_dir, out_prefix, "3snps_", phecode, "_assoc.RDS")
  assoc    <- list(null=nullmod, assoc=assoc)
  saveRDS(assoc, file = out_file)
  
  seqClose(gds)
  cat("Finished. Results written to", out_file, "\n") 