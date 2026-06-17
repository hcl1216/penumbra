# Fork-2 GATE C: primary spatial feature = CD8<->CK+ tumor proximity, on MEYER ONLY.
# METABRIC is NOT touched (Gate D is separate). Distances/neighbors via imcRtools (NOT hand-rolled).
# Runs on Colab via the same env-var pattern; artifacts to stdout + Drive results/.
#
# LOCKED, PRE-REGISTERED spec (do NOT tune against survival):
#   PRIMARY   = infiltrated-fraction: proportion of a patient's CD8 T cells within 20 um of a CK+ tumor
#               cell. 20 um = contact/juxtacrine range, biology-locked (one frozen radius, not swept).
#   SECONDARY = median CD8 -> nearest CK+ tumor distance per patient (exploratory, multiplicity-FLAGGED).
#   Cell types come from the committed cross-cohort typing (results/fork2_celltypes_meyer.csv).
#   Distances computed PER IMAGE/core (img_id=sample_id); cells never pooled across cores in coordinate
#   space -- only the per-patient SUMMARY pools a patient's CD8 cells across their cores.
#   Coordinate units assumed micrometres (steinbock Pos_X/Pos_Y, 1 px = 1 um at Hyperion 1 um res).
#
# HARD STOP after the Meyer fit + type-preserving permutation control. Does NOT lock the feature for
# Gate D and does NOT read METABRIC.

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
RADIUS    <- 20   # um -- LOCKED
dir.create(RES, showWarnings=FALSE, recursive=TRUE)
emit <- function(txt,f){cat(txt); writeLines(txt, file.path(RES,f))}
cat(sprintf("[paths] MEYER_SCE=%s\n        CELLTYPES=%s\n        RES=%s  N_PERM=%d  RADIUS=%dum\n",
            MEYER_SCE, CELLTYPES, RES, N_PERM, RADIUS))

## ---- load SCE + attach committed cell types (assert alignment) ----
sce <- readRDS(MEYER_SCE)
ct  <- read.csv(CELLTYPES, stringsAsFactors=FALSE)
stopifnot(nrow(ct)==ncol(sce))
stopifnot(all(as.character(ct$PID)==as.character(colData(sce)$PID)))   # same SCE column order
sce$cell_type <- ct$cell_type
sce$PID <- as.character(colData(sce)$PID)
## order cells by core up front, so imcRtools graph/aggregation (which return objects ordered by
## sample_id) never reorder us out of alignment with the label/PID/distance vectors.
sce <- sce[, order(as.character(colData(sce)$sample_id))]
img <- as.character(colData(sce)$sample_id)
is_cd8 <- sce$cell_type=="CD8_T"; is_tum <- sce$cell_type=="Tumor_epithelial"
cat(sprintf("[cells] n=%d  CD8=%d  Tumor=%d  cores=%d  patients=%d\n",
            ncol(sce), sum(is_cd8), sum(is_tum), length(unique(img)), length(unique(sce$PID))))

## ---- nearest CK+ tumour distance per cell, PER IMAGE ----
## (imcRtools::minDistToCells is not exported in imcRtools 1.18.1, so use RANN::nn2 -- the ANN
##  nearest-neighbour library imcRtools itself uses; a library NN query per image, NOT a hand-rolled
##  distance matrix. Used only for the SECONDARY + reporting; the PRIMARY uses the imcRtools graph below.)
if (!requireNamespace("RANN", quietly=TRUE))
  install.packages("RANN", repos="https://cloud.r-project.org", quiet=TRUE)
stopifnot(requireNamespace("RANN", quietly=TRUE))
coords_m <- as.matrix(as.data.frame(colData(sce))[, c("Pos_X","Pos_Y")])
d <- rep(Inf, ncol(sce))                          # Inf where a core has no tumour cell
for (im in unique(img)) {
  ix <- which(img==im); ti <- ix[is_tum[ix]]
  if (length(ti)==0) next
  nn <- RANN::nn2(data=coords_m[ti,,drop=FALSE], query=coords_m[ix,,drop=FALSE], k=1)
  d[ix] <- nn$nn.dists[,1]
}

## ---- imcRtools 20um expansion graph (built ONCE; positions fixed across permutations) ----
sce <- buildSpatialGraph(sce, img_id="sample_id", type="expansion", threshold=RADIUS,
                         coords=c("Pos_X","Pos_Y"), name="exp20")

## per-patient feature builder from a cell_type labeling (uses the FIXED graph for primary)
infiltrated_fraction <- function(labels){
  s <- sce; s$lab_tmp <- factor(labels)
  s <- aggregateNeighbors(s, colPairName="exp20", aggregate_by="metadata",
                          count_by="lab_tmp", name="aggN")   # uses the FIXED graph; sce pre-sorted -> order kept
  prop <- colData(s)[["aggN"]]                       # per-cell neighbour-type proportions (from RETURNED object)
  tum_col <- if (!is.null(prop) && "Tumor_epithelial" %in% colnames(prop)) prop[,"Tumor_epithelial"] else rep(0, ncol(s))
  cd8 <- labels=="CD8_T"
  has_tum_nb <- cd8 & !is.na(tum_col) & tum_col > 0
  tapply(has_tum_nb[cd8], sce$PID[cd8], mean)        # fraction of patient's CD8 with a tumour neighbour <=20um
}

## ---- REAL per-patient features ----
pid <- sce$PID
primary  <- infiltrated_fraction(sce$cell_type)
finite_d <- is.finite(d)
secondary<- tapply((d)[is_cd8 & finite_d], pid[is_cd8 & finite_d], median)   # median nearest-tumour dist (um)
cd8frac  <- tapply(is_cd8, pid, mean)                                        # composition floor
ncd8     <- tapply(is_cd8, pid, sum)
# CD8 cells sitting in tumour-free cores (no finite distance) -- reported, NOT silently dropped
cd8_inf  <- tapply(is_cd8 & !finite_d, pid, sum)

## ---- assemble patient table + survival ----
clin <- as.data.frame(colData(sce))[!duplicated(pid), c("PID","status_DFS","DFS_months")]
P <- data.frame(PID=names(primary), primary=as.numeric(primary), stringsAsFactors=FALSE)
P$secondary <- as.numeric(secondary[P$PID]); P$cd8frac <- as.numeric(cd8frac[P$PID])
P$n_cd8 <- as.integer(ncd8[P$PID]); P$cd8_in_tumourfree_core <- as.integer(cd8_inf[P$PID])
P <- merge(P, clin, by="PID")
P$status_DFS <- as.integer(P$status_DFS)
n_zero_cd8 <- sum(is.na(P$primary) | P$n_cd8==0)
fit_set <- P[is.finite(P$primary) & P$n_cd8>0 & is.finite(P$DFS_months) & P$DFS_months>0 & !is.na(P$status_DFS), ]
write.csv(P, file.path(RES,"fork2_gateC_features.csv"), row.names=FALSE)

## standardize predictors (HR per 1 SD; C-index is scale-invariant)
z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
fs <- fit_set; fs$primary_z<-z(fs$primary); fs$secondary_z<-z(fs$secondary); fs$cd8frac_z<-z(fs$cd8frac)
cox1 <- function(f,dat) coxph(as.formula(f), data=dat)
cidx <- function(m) summary(m)$concordance[1]
m_floor <- cox1("Surv(DFS_months,status_DFS) ~ cd8frac_z", fs)
m_prim  <- cox1("Surv(DFS_months,status_DFS) ~ primary_z", fs)
m_sec   <- cox1("Surv(DFS_months,status_DFS) ~ secondary_z", fs[is.finite(fs$secondary_z),])
sm <- function(m){ s<-summary(m); c(HR=s$coef[1,"exp(coef)"], lo=s$conf.int[1,"lower .95"],
                                    hi=s$conf.int[1,"upper .95"], p=s$coef[1,"Pr(>|z|)"], C=cidx(m)) }
sp <- sm(m_prim); ss <- sm(m_sec); sfl <- sm(m_floor)
dC_prim <- sp["C"]-sfl["C"]; dC_sec <- ss["C"]-sfl["C"]

## ---- TYPE-PRESERVING within-core permutation null for the PRIMARY ----
## shuffle cell_type WITHIN each core (preserves #CD8/#tumour and the position set; breaks type<->position).
## positions fixed -> reuse the prebuilt 20um graph; recompute infiltrated-fraction + Cox each iter.
perm_labels <- function(){ out<-sce$cell_type
  for(im in unique(img)){ ix<-which(img==im); out[ix]<-sample(out[ix]) }; out }
real_z <- summary(m_prim)$coef[1,"z"]; real_C <- sp["C"]
null_z <- numeric(N_PERM); null_C <- numeric(N_PERM)
for(b in seq_len(N_PERM)){
  pf <- infiltrated_fraction(perm_labels())
  dd <- fs; dd$pf <- z(as.numeric(pf[dd$PID]))
  mb <- tryCatch(coxph(Surv(DFS_months,status_DFS) ~ pf, data=dd), error=function(e)NULL)
  if(is.null(mb)){ null_z[b]<-NA; null_C[b]<-NA } else { null_z[b]<-summary(mb)$coef[1,"z"]; null_C[b]<-cidx(mb) }
}
null_z <- null_z[is.finite(null_z)]; null_C <- null_C[is.finite(null_C)]
perm_p <- (1+sum(abs(null_z)>=abs(real_z)))/(1+length(null_z))     # two-sided, |z|

## ---- DENSITY / CORE confound check ----
area_pp  <- tapply(as.numeric(colData(sce)$area), pid, mean)
ncell_pp <- tapply(rep(1,ncol(sce)), pid, sum)
dens <- data.frame(PID=fs$PID, primary=fs$primary, n_cd8=fs$n_cd8,
                   area=as.numeric(area_pp[fs$PID]), ncell=as.numeric(ncell_pp[fs$PID]))
sp_cd8  <- suppressWarnings(cor.test(dens$primary, dens$n_cd8, method="spearman"))
sp_area <- suppressWarnings(cor.test(dens$primary, dens$area, method="spearman"))
sp_ncell<- suppressWarnings(cor.test(dens$primary, dens$ncell, method="spearman"))

## convergent vs divergent (more infiltration should track shorter distance -> features anti-correlated)
conv <- suppressWarnings(cor(fs$primary, fs$secondary, use="complete.obs", method="spearman"))

## ---- report ----
out <- paste0(
"# R5 GATE C: CD8<->CK+ tumour proximity (Meyer only; PRE-REGISTERED; PENDING REVIEW)\n\n",
sprintf("Feature set: %d patients with computable primary + DFS (of %d patients; %d excluded for zero CD8).\n",
        nrow(fs), length(unique(pid)), n_zero_cd8),
sprintf("CD8 cells in tumour-free cores (no finite distance) total = %d; these CD8 count as NOT infiltrating\n",
        sum(P$cd8_in_tumourfree_core, na.rm=TRUE)),
"(numerator) and remain in the patient denominator; they are excluded from the SECONDARY median only.\n\n",
"## Distributions\n",
sprintf("primary infiltrated-fraction (<=%dum): median %.3f  IQR [%.3f, %.3f]  range [%.3f, %.3f]\n",
        RADIUS, median(fs$primary), quantile(fs$primary,.25), quantile(fs$primary,.75), min(fs$primary), max(fs$primary)),
sprintf("secondary median dist (um):           median %.1f  IQR [%.1f, %.1f]\n",
        median(fs$secondary,na.rm=TRUE), quantile(fs$secondary,.25,na.rm=TRUE), quantile(fs$secondary,.75,na.rm=TRUE)),
"\n## Cox (Meyer DFS; HR per 1 SD; C-index)\n",
sprintf("composition floor (CD8 frac): HR %.3f [%.3f-%.3f] p=%.3f  C=%.3f\n", sfl["HR"],sfl["lo"],sfl["hi"],sfl["p"],sfl["C"]),
sprintf("PRIMARY (infiltr. frac):      HR %.3f [%.3f-%.3f] p=%.3f  C=%.3f   dC vs floor = %+.3f\n",
        sp["HR"],sp["lo"],sp["hi"],sp["p"],sp["C"],dC_prim),
sprintf("secondary (median dist) [EXPLORATORY]: HR %.3f [%.3f-%.3f] p=%.3f  C=%.3f  dC=%+.3f\n",
        ss["HR"],ss["lo"],ss["hi"],ss["p"],ss["C"],dC_sec),
"\n## Type-preserving within-core permutation null for PRIMARY (n=", length(null_z), " valid)\n",
sprintf("real |z| = %.2f ; null |z|: mean %.2f, 95th pctile %.2f ; permutation p = %.4f\n",
        abs(real_z), mean(abs(null_z)), quantile(abs(null_z),.95), perm_p),
sprintf("real C = %.3f ; permuted-feature C: mean %.3f (should collapse toward floor C=%.3f)\n",
        real_C, mean(null_C), sfl["C"]),
sprintf("=> spatial-ness: %s\n", ifelse(perm_p<0.05 && real_C > mean(null_C),
        "real signal EXCEEDS type-preserving null (spatial, not composition)",
        "real signal NOT clearly above null -- treat as composition/inconclusive")),
"\n## Density / core confound (Spearman vs primary)\n",
sprintf("vs n_CD8:      rho %+.3f  p=%.3g\n", sp_cd8$estimate, sp_cd8$p.value),
sprintf("vs core area:  rho %+.3f  p=%.3g\n", sp_area$estimate, sp_area$p.value),
sprintf("vs core ncell: rho %+.3f  p=%.3g\n", sp_ncell$estimate, sp_ncell$p.value),
"(fraction formulation should largely control density; large rho => covariate caveat for Gate D)\n",
"\n## Primary vs secondary agreement\n",
sprintf("Spearman(primary, secondary) = %+.3f  -> %s\n", conv,
        ifelse(conv <= -0.4, "CONVERGENT (more infiltration ~ shorter distance) = architecture signal",
        ifelse(conv >= 0.4, "unexpected same-sign -- inspect", "weak/divergent = signal may be in the infiltrating subset specifically"))),
"\nHARD STOP: Meyer-only. Feature NOT locked for Gate D (reviewer call). METABRIC untouched.\n")
emit(out, "fork2_R5_gateC.md")
cat("\nDONE Gate C (Meyer). Artifacts: results/fork2_R5_gateC.md + fork2_gateC_features.csv\n")
