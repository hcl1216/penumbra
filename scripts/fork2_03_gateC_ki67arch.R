# Fork-2 GATE C / Secondary-2 (LAST feature): Ki67 tumour proliferative ARCHITECTURE, Meyer ONLY.
# Independent of the killed CD8<->tumour axis. METABRIC untouched. Spatial/graph/permutation via imcRtools
# (NO hand-rolled distance/permutation). Concordance via the fixed survival::concordance accessor.
# PRE-REGISTERED (DECISIONS.md 2026-06-16, commit 5e52b24) -- spec frozen, NOT tuned vs survival.
#
# PRIMARY = Ki67 assortativity ENRICHMENT among tumour cells:
#   Ki67+ tumour cell = tumour cell with Ki67(exprs) > Otsu(pooled Meyer tumour Ki67); if that pooled
#     distribution is unimodal (dip test p>=0.05) -> fall back to cohort MEDIAN + flag.
#   kNN k=6 graph on TUMOUR CELLS ONLY, per core (imcRtools buildSpatialGraph).
#   For each Ki67+ tumour cell: fraction of its 6 tumour-neighbours that are Ki67+ (aggregateNeighbors),
#     / patient baseline Ki67+ rate among graph-eligible tumour cells (=1 random, >1 clustered).
#   Per-patient scalar = mean enrichment over the patient's Ki67+ tumour cells. Require >=10 Ki67+ tumour
#   cells/patient (graph-eligible) else excluded+counted.
# FLOOR = Ki67+ tumour fraction per patient (proliferation index; likely NON-flat).
# Controls: type-preserving within-core permutation of Ki67+ labels among tumour cells (N=499);
#   density/confound Spearman. HARD STOP after Meyer fit + permutation; not locked for Gate D.

suppressMessages({library(SingleCellExperiment); library(imcRtools); library(survival)})
set.seed(1)

## ---- paths (env > repo default) ----
.args <- commandArgs(FALSE); .fa <- grep("^--file=", .args, value=TRUE)
REPO <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=","",.fa[1])))) else getwd()
getp <- function(v,d){x<-Sys.getenv(v); if(nzchar(x)) x else d}
MEYER_SCE <- getp("MEYER_SCE", file.path(REPO,"data","meyer","sce_ALL_sub.rds"))
RES       <- getp("PENUMBRA_RESULTS_DIR", file.path(REPO,"results"))
CELLTYPES <- getp("MEYER_CELLTYPES", file.path(RES,"fork2_celltypes_meyer.csv"))
N_PERM    <- as.integer(getp("PERM_N","499"))
K_NN      <- 6L; MIN_KI67 <- 10L
dir.create(RES, showWarnings=FALSE, recursive=TRUE)
emit <- function(txt,f){cat(txt); writeLines(txt, file.path(RES,f))}
cidx <- function(m) tryCatch(as.numeric(survival::concordance(m)$concordance), error=function(e) NA_real_)
for (p in c("EBImage","diptest")) if (!requireNamespace(p, quietly=TRUE))
  install.packages(p, repos="https://cloud.r-project.org", quiet=TRUE)
cat(sprintf("[paths] MEYER_SCE=%s\n        RES=%s  N_PERM=%d  k=%d  min_Ki67=%d\n",
            MEYER_SCE, RES, N_PERM, K_NN, MIN_KI67))

## ---- load + attach committed typing; pre-sort by core for stable imcRtools alignment ----
sce <- readRDS(MEYER_SCE)
ct  <- read.csv(CELLTYPES, stringsAsFactors=FALSE)
stopifnot(nrow(ct)==ncol(sce), all(as.character(ct$PID)==as.character(colData(sce)$PID)))
sce$cell_type <- ct$cell_type; sce$PID <- as.character(colData(sce)$PID)
sce <- sce[, order(as.character(colData(sce)$sample_id))]
is_tum <- sce$cell_type=="Tumor_epithelial"
ki67 <- as.numeric(assay(sce,"exprs")["Ki67",])
cat(sprintf("[cells] n=%d tumour=%d cores=%d patients=%d\n",
            ncol(sce), sum(is_tum), length(unique(colData(sce)$sample_id)), length(unique(sce$PID))))

## ---- Otsu threshold on pooled tumour Ki67 (survival-blind); median fallback if unimodal ----
ki67_t <- ki67[is_tum]
otsu <- EBImage::otsu(ki67_t, range=range(ki67_t), levels=256L)
dip  <- diptest::dip.test(ki67_t)
unimodal <- dip$p.value >= 0.05
thr <- if (unimodal) median(ki67_t) else otsu
thr_note <- sprintf("Otsu=%.4f ; dip p=%.4g (%s) ; THRESHOLD USED = %.4f (%s)",
                    otsu, dip$p.value, ifelse(unimodal,"UNIMODAL","bimodal"), thr,
                    ifelse(unimodal,"median fallback","Otsu"))
is_ki67pos <- is_tum & (ki67 > thr)

## ---- composition floor: Ki67+ tumour fraction per patient (all tumour cells) ----
pid_all <- sce$PID
ki67frac <- tapply(is_ki67pos[is_tum], pid_all[is_tum], mean)

## ---- tumour-only SCE; keep cores with >= k+1 tumour cells (kNN k=6 needs >=7) ----
sct <- sce[, is_tum]
core_t <- as.character(colData(sct)$sample_id)
big <- names(table(core_t))[table(core_t) >= (K_NN+1)]
sct <- sct[, core_t %in% big]
sct$ki67_status <- factor(ifelse(as.numeric(assay(sct,"exprs")["Ki67",]) > thr, "pos","neg"),
                          levels=c("neg","pos"))
pid_t <- sct$PID
n_dropped_smallcore <- sum(is_tum) - ncol(sct)

## ---- imcRtools kNN graph (k=6) on tumour cells, per core ----
sct <- buildSpatialGraph(sct, img_id="sample_id", type="knn", k=K_NN, coords=c("Pos_X","Pos_Y"), name="knn")

## ---- enrichment feature builder from a Ki67-status labeling (uses FIXED graph) ----
baseline <- tapply(sct$ki67_status=="pos", pid_t, mean)          # patient Ki67+ rate (graph-eligible) -- fixed under within-core shuffle
assort_enrichment <- function(status){                          # status: factor neg/pos over sct cells
  s <- sct; s$lab_tmp <- factor(status, levels=c("neg","pos"))
  s <- aggregateNeighbors(s, colPairName="knn", aggregate_by="metadata", count_by="lab_tmp", name="aggN")
  prop <- colData(s)[["aggN"]]
  fpos <- if (!is.null(prop) && "pos" %in% colnames(prop)) prop[,"pos"] else rep(0, ncol(s))
  enr  <- fpos / baseline[pid_t]                                # enrichment over patient baseline
  ispos <- status=="pos"
  tapply(enr[ispos], pid_t[ispos], mean)                        # per-patient mean enrichment over Ki67+ tumour cells
}

## ---- REAL feature + min-cell filter ----
primary <- assort_enrichment(sct$ki67_status)
nki67   <- tapply(sct$ki67_status=="pos", pid_t, sum)
excl    <- names(nki67)[nki67 < MIN_KI67]
primary[excl] <- NA

## ---- assemble patients + survival ----
clin <- as.data.frame(colData(sce))[!duplicated(pid_all), c("PID","status_DFS","DFS_months")]
P <- data.frame(PID=names(primary), primary=as.numeric(primary), stringsAsFactors=FALSE)
P$ki67frac <- as.numeric(ki67frac[P$PID]); P$n_ki67 <- as.integer(nki67[P$PID])
P <- merge(P, clin, by="PID"); P$status_DFS <- as.integer(P$status_DFS)
ok <- is.finite(P$DFS_months) & P$DFS_months>0 & !is.na(P$status_DFS)
fit_all  <- P[ok & is.finite(P$ki67frac), ]                     # floor available
fit_prim <- P[ok & is.finite(P$primary), ]                      # primary computable (>=10 Ki67+)
write.csv(P, file.path(RES,"fork2_gateC_ki67_features.csv"), row.names=FALSE)
n_excl_min <- sum(ok & !is.finite(P$primary)); ev_excl <- sum(P$status_DFS[ok & !is.finite(P$primary)]==1)

## ---- Cox: floor first, then primary (same patient set for fair dC); HR per 1 SD; two-sided ----
z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
fp <- fit_prim; fp$primary_z<-z(fp$primary); fp$ki67frac_z<-z(fp$ki67frac)
m_floor <- coxph(Surv(DFS_months,status_DFS) ~ ki67frac_z, data=fp)
m_prim  <- coxph(Surv(DFS_months,status_DFS) ~ primary_z,  data=fp)
m_both  <- coxph(Surv(DFS_months,status_DFS) ~ ki67frac_z + primary_z, data=fp)
sm <- function(m,term){ s<-summary(m); c(HR=s$coef[term,"exp(coef)"], lo=s$conf.int[term,"lower .95"],
                                         hi=s$conf.int[term,"upper .95"], p=s$coef[term,"Pr(>|z|)"]) }
C_floor<-cidx(m_floor); C_prim<-cidx(m_prim); C_both<-cidx(m_both)

## ---- type-preserving within-core permutation of Ki67+ labels among tumour cells ----
img_t <- as.character(colData(sct)$sample_id)
perm_status <- function(){ out<-as.character(sct$ki67_status)
  for(im in unique(img_t)){ ix<-which(img_t==im); out[ix]<-sample(out[ix]) }
  factor(out, levels=c("neg","pos")) }
real_z <- summary(m_prim)$coef["primary_z","z"]
null_z <- numeric(N_PERM); null_C <- numeric(N_PERM)
for(b in seq_len(N_PERM)){
  pe <- assort_enrichment(perm_status())
  dd <- fp; dd$pz <- z(as.numeric(pe[dd$PID]))
  mb <- tryCatch(coxph(Surv(DFS_months,status_DFS) ~ pz, data=dd), error=function(e)NULL)
  if(is.null(mb)){ null_z[b]<-NA; null_C[b]<-NA } else { null_z[b]<-summary(mb)$coef["pz","z"]; null_C[b]<-cidx(mb) }
}
null_z<-null_z[is.finite(null_z)]; null_C<-null_C[is.finite(null_C)]
perm_p <- (1+sum(abs(null_z)>=abs(real_z)))/(1+length(null_z))

## ---- density / confound ----
area_pp <- tapply(as.numeric(colData(sct)$area), pid_t, mean)
ntum_pp <- tapply(rep(1,ncol(sct)), pid_t, sum)
dC <- fp; dC$area<-as.numeric(area_pp[dC$PID]); dC$ntum<-as.numeric(ntum_pp[dC$PID])
sp <- function(a,b) suppressWarnings(cor.test(a,b,method="spearman"))
s_ntum<-sp(dC$primary,dC$ntum); s_area<-sp(dC$primary,dC$area); s_base<-sp(dC$primary,dC$ki67frac)

## ---- verdict ----
beats_floor <- (C_prim - C_floor) > 0 && !(sm(m_both,"primary_z")["lo"]<=1 && sm(m_both,"primary_z")["hi"]>=1)
spatial_ok  <- perm_p < 0.05 && C_prim > mean(null_C)
verdict <- if (beats_floor && spatial_ok) "PASS (additive over Ki67-fraction floor AND exceeds type-preserving null)" else
           "FAIL -> KILL: not additive over floor and/or does not exceed the spatial null"

out <- paste0(
"# R6 GATE C Secondary-2: Ki67 tumour proliferative architecture (Meyer only; PRE-REGISTERED; PENDING REVIEW)\n\n",
"Ki67+ threshold: ", thr_note, "\n",
sprintf("Tumour cells in cores with <%d tumour cells (excluded from kNN graph): %d\n", K_NN+1, n_dropped_smallcore),
sprintf("Patients: primary computable (>=%d Ki67+ tumour cells) + DFS = %d (of 215; %d excluded by min-cell, %d events among excluded).\n",
        MIN_KI67, nrow(fit_prim), n_excl_min, ev_excl),
sprintf("Events in analysis set: %d (vs 48 in full Meyer DFS).\n\n", sum(fit_prim$status_DFS)),
"## Distributions\n",
sprintf("primary assortativity-enrichment: median %.3f IQR [%.3f, %.3f] range [%.3f, %.3f]  (1=random)\n",
        median(fit_prim$primary), quantile(fit_prim$primary,.25), quantile(fit_prim$primary,.75),
        min(fit_prim$primary), max(fit_prim$primary)),
sprintf("Ki67+ tumour fraction (floor): median %.3f IQR [%.3f, %.3f]\n\n",
        median(fit_prim$ki67frac), quantile(fit_prim$ki67frac,.25), quantile(fit_prim$ki67frac,.75)),
"## Cox (Meyer DFS; HR per 1 SD; two-sided; same patient set)\n",
sprintf("FLOOR Ki67+ fraction:   HR %.3f [%.3f-%.3f] p=%.3f  C=%.3f\n",
        sm(m_floor,"ki67frac_z")["HR"],sm(m_floor,"ki67frac_z")["lo"],sm(m_floor,"ki67frac_z")["hi"],sm(m_floor,"ki67frac_z")["p"],C_floor),
sprintf("PRIMARY assortativity:  HR %.3f [%.3f-%.3f] p=%.3f  C=%.3f   dC vs floor = %+.3f\n",
        sm(m_prim,"primary_z")["HR"],sm(m_prim,"primary_z")["lo"],sm(m_prim,"primary_z")["hi"],sm(m_prim,"primary_z")["p"],C_prim, C_prim-C_floor),
sprintf("JOINT (floor+primary):  primary HR %.3f [%.3f-%.3f] p=%.3f ; C=%.3f (dC over floor = %+.3f)\n",
        sm(m_both,"primary_z")["HR"],sm(m_both,"primary_z")["lo"],sm(m_both,"primary_z")["hi"],sm(m_both,"primary_z")["p"],C_both, C_both-C_floor),
ifelse(sm(m_prim,"primary_z")["HR"]>1,"direction: more-compartmentalised -> WORSE DFS (as hypothesised)\n",
       "direction: more-compartmentalised -> BETTER DFS (SURPRISING -- flag)\n"),
"\n## Type-preserving permutation null (Ki67 labels shuffled within core among tumour cells; n=", length(null_z), ")\n",
sprintf("real |z| = %.2f ; null |z| mean %.2f, 95th %.2f ; permutation p = %.4f\n",
        abs(real_z), mean(abs(null_z)), quantile(abs(null_z),.95), perm_p),
sprintf("real C = %.3f ; permuted-feature C mean %.3f (should collapse toward floor C=%.3f)\n", C_prim, mean(null_C), C_floor),
"\n## Density / confound (Spearman vs primary)\n",
sprintf("vs tumour-cell count: rho %+.3f p=%.3g\n", s_ntum$estimate, s_ntum$p.value),
sprintf("vs core area:         rho %+.3f p=%.3g\n", s_area$estimate, s_area$p.value),
sprintf("vs baseline Ki67+ rate: rho %+.3f p=%.3g\n", s_base$estimate, s_base$p.value),
"\n## VERDICT (pre-registered kill criteria)\n", verdict, "\n",
"\nHARD STOP: Meyer-only, LAST pre-registered feature. NOT locked for Gate D (reviewer call). METABRIC untouched.\n")
emit(out, "fork2_R6_gateC_ki67.md")
cat("\nDONE R6. Artifacts: results/fork2_R6_gateC_ki67.md + fork2_gateC_ki67_features.csv\n")
