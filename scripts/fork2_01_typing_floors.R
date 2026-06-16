# Fork-2 front end: harmonize -> ONE cross-cohort cell-typing -> Gate B canary -> Gate A floors.
# Runs locally AND on Colab (paths from env vars; defaults to repo layout). Logic lives HERE, not in
# the notebook. STOPS after Gate A floors -- does NOT compute the Gate-C proximity feature, does NOT
# hand-roll clustering (uses SingleCellExperiment + FlowSOM).
#
# WHY COLAB: needs FlowSOM (Bioconductor). bioconductor.org is unreachable from the local Windows box
# (CRAN/GitHub/Zenodo work); Rphenograph needs a compiler and Rtools is absent locally. So this runs
# on Colab. See DECISIONS.md 2026-06-16.
#
# Env vars (notebook sets these to the Drive paths; locally they default to the repo layout):
#   MEYER_SCE             Meyer SingleCellExperiment .rds
#   METABRIC_SC           METABRIC SingleCells.fst
#   PENUMBRA_RESULTS_DIR  output dir for R1-R4 artifacts (also printed to stdout)
# (METABRIC IMCClinical.fst is NOT needed here -- it is for Gate D, not run.)

suppressMessages({library(SingleCellExperiment); library(fst); library(FlowSOM); library(survival)})
set.seed(1)

## ---- resolve paths: env var > repo-relative default ----
.args <- commandArgs(FALSE); .fa <- grep("^--file=", .args, value=TRUE)
REPO <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=","",.fa[1])))) else getwd()
getp <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
MEYER_SCE  <- getp("MEYER_SCE",   file.path(REPO,"data","meyer","sce_ALL_sub.rds"))
METABRIC_SC<- getp("METABRIC_SC", file.path(REPO,"data","metabric","SingleCells.fst"))
RES        <- getp("PENUMBRA_RESULTS_DIR", file.path(REPO,"results"))
dir.create(RES, showWarnings=FALSE, recursive=TRUE)
emit <- function(txt, file) { cat(txt); writeLines(txt, file.path(RES, file)) }
cat(sprintf("[paths] MEYER_SCE=%s\n        METABRIC_SC=%s\n        RES=%s\n", MEYER_SCE, METABRIC_SC, RES))

## ---- the LOCKED 20-marker feature space, with cohort-specific channel names ----
## Combined channels: CK8 & CK18 share ONE channel (CK8-18); cleaved-Casp3 & cleaved-PARP share ONE
## (Meyer c_Cas3_PARP / METABRIC c-Caspase3); vWF is measured by the combined CD31-vWF channel.
locked20 <- data.frame(stringsAsFactors=FALSE,
 marker     = c("CD3","CD4","CD8","CD11c","CD15","CD20","CD68","HLA-DR","SMA","vWF",
                "CK5","CK8","CK18","panCK","ER","HER2","Ki-67","Histone H3","cleaved-Casp3","cleaved-PARP"),
 compartment= c(rep("immune",8),"stromal","stromal",rep("epithelial-tumor",7),"other","other","other"),
 canon      = c("CD3","CD4","CD8","CD11c","CD15","CD20","CD68","HLA_DR","SMA","CD31_vWF",
                "CK5","CK8_18","CK8_18","panCK","ER","HER2","Ki67","HistoneH3","c_Casp3","c_Casp3"),
 meyer_ch   = c("CD3","CD4","CD8a","CD11c","CD15","CD20","CD68","HLADR","SMA","CD31_vWF",
                "CK5","CK8_18","CK8_18","panCK","ER","HER2","Ki67",NA,"c_Cas3_PARP","c_Cas3_PARP"),
 metab_ch   = c("CD3","CD4","CD8","CD11c","CD15","CD20","CD68","HLA-DR","SMA","CD31-vWF",
                "CK5","CK8-18","CK8-18","panCK","ER","HER2 (3B5)","Ki-67","Histone H3","c-Caspase3","c-Caspase3"))

## ---- load Meyer (exprs = arcsinh cofactor 1) ----
meyer <- readRDS(MEYER_SCE)
mey_rows <- rownames(meyer)
## ---- load METABRIC columns we need (raw intensity) ----
mb_have <- fst::metadata_fst(METABRIC_SC)$columnNames

## ---- R1: MARKER STATUS TABLE (resolve 20 -> 19 -> 17) ----
locked20$in_Meyer    <- !is.na(locked20$meyer_ch) & locked20$meyer_ch %in% mey_rows
locked20$in_METABRIC <- locked20$metab_ch %in% mb_have
# a canonical channel is usable iff present in BOTH cohorts (combined channels handled via canon)
canon_ok <- tapply(seq_len(nrow(locked20)), locked20$canon, function(ix)
                   all(locked20$in_Meyer[ix]) && all(locked20$in_METABRIC[ix]))
locked20$harmonized <- ifelse(locked20$canon=="HistoneH3" & !locked20$in_Meyer, "—",
                       ifelse(canon_ok[locked20$canon], locked20$canon, "—"))
# harmonized channel set (unique canon usable in both)
harm_canon <- unique(locked20$canon[locked20$harmonized!="—"])
mk <- unique(locked20[locked20$harmonized!="—", c("canon","meyer_ch","metab_ch")])
mk <- mk[!duplicated(mk$canon), ]

## ---- harmonize matrices (Meyer exprs; METABRIC arcsinh cofactor 1) ----
mey_e <- as.matrix(assay(meyer,"exprs"))[mk$meyer_ch,,drop=FALSE]; rownames(mey_e) <- mk$canon
mb <- read_fst(METABRIC_SC, columns=c("metabric_id", mk$metab_ch))
mb_e <- t(asinh(as.matrix(mb[, mk$metab_ch, drop=FALSE]))); rownames(mb_e) <- mk$canon
mb_pid <- as.character(mb$metabric_id); mey_pid <- as.character(colData(meyer)$PID)

## zero-variance guard (drop a channel if flat in either cohort) + record reason
zv <- sapply(mk$canon, function(g) var(mey_e[g,])==0 || var(mb_e[g,])==0)
if (any(zv)) { mk <- mk[!zv,]; mey_e <- mey_e[mk$canon,]; mb_e <- mb_e[mk$canon,]
               harm_canon <- setdiff(harm_canon, names(zv)[zv]) }
status <- function(r){
  if (r$canon=="HistoneH3" && !r$in_Meyer) return("DROPPED: absent in Meyer analysis SCE (structural)")
  if (!(r$canon %in% harm_canon)) return(if (zv[r$canon] %in% TRUE) "DROPPED: zero-variance" else "DROPPED: not in both")
  dup <- sum(locked20$canon==r$canon)
  if (dup>1) return(sprintf("KEPT (combined channel '%s' covers %d markers)", r$canon, dup))
  "KEPT"
}
locked20$status <- vapply(seq_len(nrow(locked20)), function(i) status(locked20[i,]), character(1))
r1 <- locked20[,c("marker","compartment","in_Meyer","in_METABRIC","harmonized","status")]
r1txt <- paste0(
 "# R1 MARKER STATUS (locked 20 -> harmonized analysis channels)\n\n",
 paste(capture.output(print(r1, row.names=FALSE)), collapse="\n"),
 sprintf("\n\n20 markers -> %d present-in-both concepts -> %d unique analysis channels (combined: CK8-18, c-Casp3/PARP; vWF via CD31-vWF).\n",
         sum(r1$harmonized!="—"), length(unique(mk$canon))))
emit(r1txt, "fork2_R1_marker_status.md"); cat("\n")

## ---- HARD GUARD: primary feature = CD8<->CK+ proximity needs CD8 AND >=1 cytokeratin ----
has_cd8 <- "CD8" %in% mk$canon
has_ck  <- any(c("CK5","CK8_18","panCK") %in% mk$canon)
if (!has_cd8 || !has_ck)
  stop(sprintf("HALT: primary feature (CD8<->CK+ proximity) impossible -- CD8 present=%s, cytokeratin present=%s. Fix harmonization before typing.", has_cd8, has_ck))
cat(sprintf("[guard] CD8 present=%s ; cytokeratin present=%s -> OK to type\n", has_cd8, has_ck))

## ---- per-cohort 99th-pct censor + 0-1 scale; pool ----
scale01 <- function(M) t(apply(M,1,function(v){q<-as.numeric(quantile(v,.99)); v<-pmin(v,q); if(q>0) v/q else v}))
mey_s <- scale01(mey_e); mb_s <- scale01(mb_e)
pooled <- t(cbind(mey_s, mb_s)); colnames(pooled) <- mk$canon
cohort <- c(rep("Meyer",ncol(mey_s)), rep("METABRIC",ncol(mb_s)))

## ---- ONE FlowSOM: train SOM on balanced subsample (RAM-safe), map ALL cells ----
N_SUB <- 200000
im <- which(cohort=="Meyer"); ib <- which(cohort=="METABRIC")
sub <- c(sample(im, min(N_SUB/2, length(im))), sample(ib, min(N_SUB/2, length(ib))))
fsom <- FlowSOM::ReadInput(flowCore::flowFrame(pooled[sub,,drop=FALSE]), transform=FALSE, scale=FALSE)
fsom <- FlowSOM::BuildSOM(fsom, colsToUse=mk$canon)
fsom <- FlowSOM::BuildMST(fsom)
mc   <- as.integer(FlowSOM::metaClustering_consensus(fsom$map$codes, k=12, seed=1))
nodes_all <- FlowSOM:::MapDataToCodes(fsom$map$codes, pooled)[,1]
meta <- mc[nodes_all]

## ---- annotate metaclusters by marker medians (transparent rule) ----
med <- t(sapply(sort(unique(meta)), function(k) apply(pooled[meta==k,,drop=FALSE],2,median)))
rownames(med) <- paste0("mc",sort(unique(meta)))
g <- function(p,m) if (m %in% names(p)) p[m] else 0
lineage <- function(p){
  if (max(g(p,"panCK"),g(p,"CK5"),g(p,"CK8_18"))>=0.30) return("Tumor_epithelial")
  if (g(p,"CD3")>=0.20 && g(p,"CD8")>=g(p,"CD4") && g(p,"CD8")>=0.20) return("CD8_T")
  if (g(p,"CD3")>=0.20 && g(p,"CD4")>0.20) return("CD4_T")
  if (g(p,"CD20")>=0.25) return("B_cell")
  if (max(g(p,"CD68"),g(p,"CD11c"),g(p,"CD15"))>=0.25) return("Myeloid")
  if (g(p,"SMA")>=0.25 || g(p,"CD31_vWF")>=0.25) return("Stroma_endothelial")
  "Other"
}
mc_type <- setNames(apply(med,1,lineage), rownames(med))
cell_type <- mc_type[paste0("mc",meta)]

## ---- R2: per-cohort fractions side by side ----
frac <- function(coh){ t<-table(factor(cell_type[cohort==coh])); round(100*as.numeric(t)/sum(t),2) }
types <- sort(unique(cell_type))
fm <- setNames(frac("Meyer"), sort(unique(cell_type[cohort=="Meyer"])))
fb <- setNames(frac("METABRIC"), sort(unique(cell_type[cohort=="METABRIC"])))
r2 <- data.frame(cell_type=types,
                 Meyer_pct=as.numeric(fm[types]), METABRIC_pct=as.numeric(fb[types]))
r2[is.na(r2)] <- 0
r2txt <- paste0("# R2 CELL-TYPE ANNOTATION + per-cohort fractions\n\n",
  "## metacluster -> type (median scaled marker profile)\n",
  paste(capture.output(print(round(med,2))), collapse="\n"), "\n\n",
  paste(capture.output(print(mc_type)), collapse="\n"), "\n\n## fractions (%)\n",
  paste(capture.output(print(r2, row.names=FALSE)), collapse="\n"), "\n")
emit(r2txt, "fork2_R2_celltypes.md"); cat("\n")

## ---- R3: GATE B canary -- numbers + margin ----
## stable quantity: per-type lineage-marker median profile, cross-cohort divergence (batch magnitude)
prof <- function(coh) sapply(types, function(ty){ ix<-which(cohort==coh & cell_type==ty)
   if(length(ix)<20) rep(NA,length(mk$canon)) else apply(pooled[ix,,drop=FALSE],2,median) })
pm <- prof("Meyer"); pb <- prof("METABRIC")
corr <- sapply(types, function(ty){ a<-pm[,ty]; b<-pb[,ty]; if(any(is.na(a))||any(is.na(b))) NA else cor(a,b) })
prof_div <- 1 - median(corr, na.rm=TRUE)                    # stable-quantity cross-cohort divergence
## survival-signal axis: per-patient CD8 fraction; batch shift vs biological cross-patient range
cd8_m <- tapply(cell_type[cohort=="Meyer"]=="CD8_T", mey_pid, mean)
cd8_b <- tapply(cell_type[cohort=="METABRIC"]=="CD8_T", mb_pid, mean)
batch_shift <- abs(mean(cd8_m,na.rm=TRUE) - mean(cd8_b,na.rm=TRUE))      # cross-cohort mean shift
bio_range   <- mean(c(sd(cd8_m,na.rm=TRUE), sd(cd8_b,na.rm=TRUE)))       # within-cohort cross-patient SD
margin <- bio_range / batch_shift
verdict <- if (median(corr,na.rm=TRUE)>=0.7 && margin>=1) "PASS (biology stable; batch < signal axis)" else
           "FLAG kill (i): batch ~= signal"
r3txt <- paste0("# R3 GATE B BATCH CANARY (numbers, not a token)\n\n",
 "Stable quantity = per-cell-type lineage-marker median profile (should be biologically invariant).\n",
 "Signal axis    = per-patient CD8-fraction spread (the composition axis survival is read along).\n\n",
 paste(capture.output(print(data.frame(cell_type=types, profile_corr=round(as.numeric(corr[types]),3)), row.names=FALSE)), collapse="\n"),
 sprintf("\n\nstable-quantity cross-cohort divergence (1 - median profile corr) = %.3f\n", prof_div),
 sprintf("signal axis: CD8-frac cross-patient SD (mean of cohorts) = %.4f\n", bio_range),
 sprintf("batch: CD8-frac cross-cohort mean shift                  = %.4f\n", batch_shift),
 sprintf("MARGIN (signal SD / batch shift) = %.2f   (>=1 => batch below signal axis)\n", margin),
 sprintf("VERDICT: %s\n", verdict))
emit(r3txt, "fork2_R3_canary.md"); cat("\n")

## ---- R4: GATE A floors, fit on Meyer DFS (FIXED scalars for METABRIC) ----
cd <- as.data.frame(colData(meyer))
cd8frac_m <- tapply(cell_type[cohort=="Meyer"]=="CD8_T", mey_pid, mean)
pt <- cd[!duplicated(cd$PID), c("PID","grade","pT_simple","pN_simple","status_DFS","DFS_months")]
pt$cd8frac <- as.numeric(cd8frac_m[as.character(pt$PID)])
pt <- pt[is.finite(pt$DFS_months) & pt$DFS_months>0 & !is.na(pt$status_DFS) & !is.na(pt$cd8frac), ]
pt$grade <- suppressWarnings(as.numeric(pt$grade))
f_comp <- coxph(Surv(DFS_months,status_DFS) ~ cd8frac, data=pt)
f_soc  <- coxph(Surv(DFS_months,status_DFS) ~ grade + factor(pT_simple) + factor(pN_simple), data=pt)
r4txt <- paste0("# R4 GATE A FLOORS (Meyer-fit Cox; fixed scalars for METABRIC)\n\n",
 sprintf("Endpoint Surv(DFS_months, status_DFS); n=%d patients, %d events.\n\n", nrow(pt), sum(pt$status_DFS)),
 "## Composition floor: ~ CD8 fraction\n",
 paste(capture.output(print(round(summary(f_comp)$coefficients,4))), collapse="\n"),
 "\n\n## Clinical SoC composite: ~ grade + stage (pT_simple + pN_simple)  [NO age -> transferable to METABRIC]\n",
 paste(capture.output(print(round(summary(f_soc)$coefficients,4))), collapse="\n"),
 "\n\nThese Meyer-fit coefficients are LOCKED. The eventual METABRIC claim is calibrated to ",
 "'additive over grade+stage' (NOT 'over SoC broadly'), because METABRIC IMCClinical has no age.\n")
emit(r4txt, "fork2_R4_floors.md"); cat("\n")

## ---- persist per-cell types (gitignored / Drive) ----
write.csv(data.frame(PID=mey_pid, cell_type=cell_type[cohort=="Meyer"]),
          file.path(RES,"fork2_celltypes_meyer.csv"), row.names=FALSE)
write.csv(data.frame(metabric_id=mb_pid, cell_type=cell_type[cohort=="METABRIC"]),
          file.path(RES,"fork2_celltypes_metabric.csv"), row.names=FALSE)
cat("DONE (Gate A). Artifacts R1-R4 in", RES, "-- HARD STOP before Gate C (no proximity feature).\n")
