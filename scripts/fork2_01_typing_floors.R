# Fork-2 front end: harmonize -> ONE cross-cohort cell-typing -> Gate B canary -> Gate A floors.
# STOPS before Gate C (no proximity/distance feature here).
#
# RUNTIME REQUIREMENT (NOT satisfiable on the local Windows box, 2026-06-16):
#   needs FlowSOM (Bioconductor). bioconductor.org is UNREACHABLE from local (curl times out;
#   CRAN/GitHub/Zenodo work). Rphenograph (GitHub) was the design's other option but needs a C++
#   compiler and Rtools is ABSENT locally. Per the project's Claude-Code/CPU vs Colab/heavy split,
#   RUN THIS ON COLAB (network-open + installs the IMC stack + handles ~1.2M cells). Data objects are
#   on Drive (data/meyer/sce_ALL_sub.rds, data/metabric/SingleCells.fst, IMCClinical.fst).
#   See DECISIONS.md 2026-06-16 (PART 3 blocker).
# Inputs (gitignored): data/meyer/sce_ALL_sub.rds, data/metabric/SingleCells.fst, data/metabric/IMCClinical.fst
# Outputs: results/fork2_typing.md, results/fork2_floors.md, results/fork2_celltypes_meyer.csv/_metabric.csv
suppressMessages({library(SingleCellExperiment); library(fst); library(FlowSOM); library(survival)})
set.seed(1)
RES <- "C:/penumbra/results"

## --- 17 shared analysis markers: canonical <- (Meyer rowname, METABRIC column) ---
mk <- data.frame(
 canon=c("CD3","CD4","CD8","CD11c","CD15","CD20","CD68","HLA_DR","SMA","CD31_vWF","CK5","CK8_18","panCK","ER","HER2","Ki67","c_Casp3"),
 meyer=c("CD3","CD4","CD8a","CD11c","CD15","CD20","CD68","HLADR","SMA","CD31_vWF","CK5","CK8_18","panCK","ER","HER2","Ki67","c_Cas3_PARP"),
 metab=c("CD3","CD4","CD8","CD11c","CD15","CD20","CD68","HLA-DR","SMA","CD31-vWF","CK5","CK8-18","panCK","ER","HER2 (3B5)","Ki-67","c-Caspase3"),
 stringsAsFactors=FALSE)

## --- load + harmonize. Meyer 'exprs' = arcsinh(counts,cofactor=1); METABRIC raw -> arcsinh(cofactor=1) ---
meyer <- readRDS("C:/penumbra/data/meyer/sce_ALL_sub.rds")
mey_e <- as.matrix(assay(meyer,"exprs"))[mk$meyer,]            # 17 x cells
rownames(mey_e) <- mk$canon
mey_pid <- as.character(colData(meyer)$PID)

mbcols <- c("metabric_id", mk$metab)
mb <- read_fst("C:/penumbra/data/metabric/SingleCells.fst", columns=mbcols)
mb_e <- t(asinh(as.matrix(mb[, mk$metab])))                   # 17 x cells, cofactor 1
rownames(mb_e) <- mk$canon
mb_pid <- as.character(mb$metabric_id)

## --- per-cohort per-marker 99th-pct censor + 0-1 scale (aligns platform scales for joint typing) ---
scale01 <- function(M){ t(apply(M,1,function(v){q<-as.numeric(quantile(v,.99)); v<-pmin(v,q); if(q>0) v/q else v})) }
mey_s <- scale01(mey_e); mb_s <- scale01(mb_e)

## --- pool + ONE FlowSOM run (consistent scheme across both cohorts) ---
pooled <- t(cbind(mey_s, mb_s))                               # cells x 17
cohort <- c(rep("Meyer",ncol(mey_s)), rep("METABRIC",ncol(mb_s)))
ff <- flowCore::flowFrame(pooled)
fs <- FlowSOM(ff, colsToUse=mk$canon, scale=FALSE, nClus=12, seed=1)
meta <- as.integer(as.character(GetMetaclusters(fs)))

## --- annotate metaclusters by marker medians (transparent rule) ---
med <- t(sapply(sort(unique(meta)), function(k) apply(pooled[meta==k,,drop=FALSE],2,median)))
rownames(med) <- paste0("mc",sort(unique(meta)))
lineage <- function(p){
  if (max(p[c("panCK","CK5","CK8_18")])>=0.30) return("Tumor_epithelial")
  if (p["CD3"]>=0.20 && p["CD8"]>=p["CD4"] && p["CD8"]>=0.20) return("CD8_T")
  if (p["CD3"]>=0.20 && p["CD4"]>0.20) return("CD4_T")
  if (p["CD20"]>=0.25) return("B_cell")
  if (max(p[c("CD68","CD11c","CD15")])>=0.25) return("Myeloid")
  if (p["SMA"]>=0.25 || p["CD31_vWF"]>=0.25) return("Stroma_endothelial")
  return("Other")
}
mc_type <- setNames(apply(med,1,lineage), rownames(med))
cell_type <- mc_type[paste0("mc",meta)]

## --- per-cohort fractions ---
frac <- function(coh){ t <- table(cell_type[cohort==coh]); round(100*t/sum(t),2) }
fr_m <- frac("Meyer"); fr_b <- frac("METABRIC")
types <- sort(union(names(fr_m),names(fr_b)))

## --- GATE B canary: biologically-stable quantity (per-type marker profiles) must be concordant
##     cross-cohort (batch must NOT be the dominant axis). ---
prof_cohort <- function(coh){ sapply(types, function(ty){ idx<-which(cohort==coh & cell_type==ty)
   if(length(idx)<10) rep(NA,nrow(mk)) else apply(pooled[idx,,drop=FALSE],2,median) }) }
pm <- prof_cohort("Meyer"); pb <- prof_cohort("METABRIC")
corr_by_type <- sapply(types, function(ty){ a<-pm[,ty]; b<-pb[,ty]; if(any(is.na(a))||any(is.na(b))) NA else round(cor(a,b),3) })

## --- write typing report ---
con <- file(file.path(RES,"fork2_typing.md")); sink(con)
cat("# Fork-2 cell-typing + Gate B batch canary (PENDING REVIEW)\n\n")
cat("Harmonized to 17 shared analysis markers. Meyer 'exprs' (arcsinh cofactor 1); METABRIC arcsinh\n")
cat("(cofactor 1) of raw intensity; both per-cohort 99th-pct censored + 0-1 scaled. ONE FlowSOM run\n")
cat("(nClus=12, seed=1) on pooled cells; metaclusters annotated by marker-median rule.\n\n")
cat("## Metacluster -> type (median scaled marker profile)\n\n")
print(round(med,2)); cat("\nassignments:\n"); print(mc_type)
cat("\n## Per-cohort cell-type fractions (%)\n\n")
ftab <- data.frame(type=types, Meyer=as.numeric(fr_m[types]), METABRIC=as.numeric(fr_b[types]))
print(ftab, row.names=FALSE)
cat("\n## GATE B canary: cross-cohort per-type marker-profile correlation\n")
cat("(stable biology preserved -> high corr -> batch NOT dominant; low corr -> batch ~= signal risk)\n\n")
print(data.frame(type=types, profile_corr=as.numeric(corr_by_type[types])), row.names=FALSE)
med_corr <- median(corr_by_type, na.rm=TRUE)
cat(sprintf("\nmedian cross-cohort profile corr = %.3f\n", med_corr))
cat(sprintf("VERDICT: %s\n", ifelse(med_corr>=0.7,
   "PASS (biology stable cross-cohort; batch not dominant) -- proceed",
   "FLAG kill criterion (i): batch ~= signal -- profiles diverge cross-cohort")))
sink(); close(con)

## --- GATE A floors (FIT ON MEYER; become fixed scalars for METABRIC) ---
cd <- as.data.frame(colData(meyer))
# per-patient CD8 fraction (composition floor) + SoC covariates
cd$ct <- cell_type[cohort=="Meyer"]
pat <- data.frame(PID=mey_pid, ct=cd$ct)
cd8frac <- tapply(pat$ct=="CD8_T", pat$PID, mean)
ptv <- cd[!duplicated(cd$PID), c("PID","grade","pT_simple","pN_simple","status_DFS","DFS_months")]
ptv$cd8frac <- as.numeric(cd8frac[as.character(ptv$PID)])
ptv <- ptv[is.finite(ptv$DFS_months) & ptv$DFS_months>0 & !is.na(ptv$status_DFS), ]
ptv$grade <- as.numeric(ptv$grade)

con <- file(file.path(RES,"fork2_floors.md")); sink(con)
cat("# Fork-2 Gate A floors (fit on Meyer DFS; fixed scalars for METABRIC) (PENDING REVIEW)\n\n")
cat("Endpoint: Surv(DFS_months, status_DFS). n patients used:", nrow(ptv), " events:", sum(ptv$status_DFS), "\n\n")
cat("## Composition floor: Cox ~ CD8 fraction\n")
f1 <- tryCatch(coxph(Surv(DFS_months,status_DFS) ~ cd8frac, data=ptv), error=function(e)e)
print(summary(f1)$coefficients)
cat("\n## Clinical SoC composite: Cox ~ grade + pT_simple + pN_simple (NO age -> transferable to METABRIC)\n")
f2 <- tryCatch(coxph(Surv(DFS_months,status_DFS) ~ grade + factor(pT_simple) + factor(pN_simple), data=ptv), error=function(e)e)
print(summary(f2)$coefficients)
cat("\nNOTE: these Meyer-fit coefficients are LOCKED and applied as fixed linear-predictor scalars to\n")
cat("METABRIC (grade + stage available there; CD8 fraction from the same typing). Gate C/D not run.\n")
sink(); close(con)

## save per-cell types for downstream (gitignored)
write.csv(data.frame(PID=mey_pid, cell_type=cell_type[cohort=="Meyer"]), file.path("C:/penumbra/data/meyer","celltypes.csv"), row.names=FALSE)
write.csv(data.frame(metabric_id=mb_pid, cell_type=cell_type[cohort=="METABRIC"]), file.path("C:/penumbra/data/metabric","celltypes.csv"), row.names=FALSE)
cat("DONE. typing+canary -> results/fork2_typing.md ; floors -> results/fork2_floors.md\n")
