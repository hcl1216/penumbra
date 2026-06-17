# Fork-2 COHORT-VALIDITY / POSITIVE-CONTROL check (Meyer ONLY). Gate before any negative writeup.
# Tests whether Meyer can detect KNOWN prognostic signal at all -> makes the R5/R6 spatial nulls
# interpretable. NOT a feature hunt, NOT salvage. Killed spatial features NOT re-run. METABRIC untouched.
# FULL cohort (no spatial / min-cell / zero-CD8 filters). Both endpoints (DFS + OS). Pre-registered
# (DECISIONS.md 2026-06-16, commit ca3c2f6). Concordance via the fixed survival::concordance accessor.

suppressMessages({library(SingleCellExperiment); library(survival)})

## ---- paths ----
.args <- commandArgs(FALSE); .fa <- grep("^--file=", .args, value=TRUE)
REPO <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=","",.fa[1])))) else getwd()
getp <- function(v,d){x<-Sys.getenv(v); if(nzchar(x)) x else d}
MEYER_SCE <- getp("MEYER_SCE", file.path(REPO,"data","meyer","sce_ALL_sub.rds"))
RES       <- getp("PENUMBRA_RESULTS_DIR", file.path(REPO,"results"))
CELLTYPES <- getp("MEYER_CELLTYPES", file.path(RES,"fork2_celltypes_meyer.csv"))
dir.create(RES, showWarnings=FALSE, recursive=TRUE)
emit <- function(txt,f){cat(txt); writeLines(txt, file.path(RES,f))}
cidx <- function(m) tryCatch(as.numeric(survival::concordance(m)$concordance), error=function(e) NA_real_)

## ---- load SCE + committed typing; build per-patient table (full cohort, no filters) ----
sce <- readRDS(MEYER_SCE)
ct  <- read.csv(CELLTYPES, stringsAsFactors=FALSE)
stopifnot(nrow(ct)==ncol(sce), all(as.character(ct$PID)==as.character(colData(sce)$PID)))
cell_type <- ct$cell_type; pid <- as.character(colData(sce)$PID)
lymph <- tapply(cell_type %in% c("CD8_T","CD4_T","B_cell"), pid, mean)            # stromal-TIL analog
immune<- tapply(cell_type %in% c("CD8_T","CD4_T","B_cell","Myeloid"), pid, mean)  # +myeloid
cd <- as.data.frame(colData(sce)); cd <- cd[!duplicated(cd$PID), ]
P <- data.frame(PID=cd$PID, grade=suppressWarnings(as.numeric(cd$grade)),
                pT=cd$pT_simple, pN=cd$pN_simple,
                status_DFS=as.integer(cd$status_DFS), DFS_months=as.numeric(cd$DFS_months),
                status_OS =as.integer(cd$status_OS),  OS_months =as.numeric(cd$OS_months),
                stringsAsFactors=FALSE)
P$lymph <- as.numeric(lymph[P$PID]); P$immune <- as.numeric(immune[P$PID])
ord_num <- function(x){ n<-suppressWarnings(as.numeric(as.character(x))); if(all(is.na(n))) as.integer(factor(x)) else n }
P$pN_ord <- ord_num(P$pN); P$pT_ord <- ord_num(P$pT)
z <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
P$lymph_z <- z(P$lymph); P$immune_z <- z(P$immune)

## ---- endpoint definitions ----
EP <- list(DFS=c(t="DFS_months", s="status_DFS"), OS=c(t="OS_months", s="status_OS"))
surv_of <- function(P,ep){ ok <- is.finite(P[[ep["t"]]]) & P[[ep["t"]]]>0 & !is.na(P[[ep["s"]]]); ok }

## ---- per-endpoint follow-up / censoring summary ----
ep_summary <- function(ep){
  ok <- surv_of(P,ep); t <- P[[ep["t"]]][ok]; s <- P[[ep["s"]]][ok]
  rk <- tryCatch(summary(survfit(Surv(t, 1-s)~1))$table["median"], error=function(e) NA)  # reverse-KM median FU
  sprintf("n=%d  events=%d  censored=%d | median follow-up (rev-KM)=%.1f mo | time range [%.1f, %.1f] | median event time=%.1f",
          sum(ok), sum(s==1), sum(s==0), as.numeric(rk),
          min(t), max(t), suppressWarnings(median(t[s==1])))
}

## ---- Cox helper: HR/CI/p/C for one term on one endpoint ----
cox_term <- function(formula_rhs, term, ep){
  ok <- surv_of(P,ep); d <- P[ok,]
  f <- as.formula(sprintf("Surv(%s,%s) ~ %s", ep["t"], ep["s"], formula_rhs))
  m <- tryCatch(coxph(f, data=d), error=function(e) NULL)
  if (is.null(m)) return(c(HR=NA,lo=NA,hi=NA,p=NA,C=NA))
  s <- summary(m); tm <- if (term %in% rownames(s$coef)) term else rownames(s$coef)[1]
  c(HR=s$coef[tm,"exp(coef)"], lo=s$conf.int[tm,"lower .95"], hi=s$conf.int[tm,"upper .95"],
    p=s$coef[tm,"Pr(>|z|)"], C=cidx(m))
}
fmt <- function(v) if (any(is.na(v))) "   NA" else sprintf("HR %.2f [%.2f-%.2f] p=%.3f C=%.3f", v["HR"],v["lo"],v["hi"],v["p"],v["C"])

## ---- controls (each on BOTH endpoints) ----
controls <- list(
  "pN ordinal trend"        = list(rhs="pN_ord",   term="pN_ord"),
  "pT ordinal trend"        = list(rhs="pT_ord",   term="pT_ord"),
  "grade (per unit)"        = list(rhs="grade",    term="grade"),
  "Total lymphocyte frac /SD"= list(rhs="lymph_z", term="lymph_z"),
  "Total immune frac /SD"   = list(rhs="immune_z", term="immune_z"))
rows <- lapply(names(controls), function(nm){
  c(name=nm, DFS=fmt(cox_term(controls[[nm]]$rhs, controls[[nm]]$term, EP$DFS)),
              OS =fmt(cox_term(controls[[nm]]$rhs, controls[[nm]]$term, EP$OS)))})

## ---- pN per-level Cox + KM log-rank (both endpoints); strata n/events ----
pn_detail <- function(ep){
  ok <- surv_of(P,ep); d <- P[ok,]; d$pNf <- factor(d$pN)
  tab <- sapply(levels(d$pNf), function(L) c(n=sum(d$pNf==L), ev=sum(d[[ep["s"]]][d$pNf==L]==1)))
  lr <- tryCatch({ sd<-survdiff(as.formula(sprintf("Surv(%s,%s)~pNf",ep["t"],ep["s"])), data=d)
                   1-pchisq(sd$chisq, length(sd$n)-1)}, error=function(e) NA)
  m <- tryCatch(summary(coxph(as.formula(sprintf("Surv(%s,%s)~pNf",ep["t"],ep["s"])), data=d))$coef,
                error=function(e) NULL)
  list(strata=tab, logrank_p=lr, coef=m)
}
pnD <- pn_detail(EP$DFS); pnO <- pn_detail(EP$OS)

## ---- verdict logic (pre-specified) ----
Cval <- function(rhs,term,ep) cox_term(rhs,term,ep)["C"]
pval <- function(rhs,term,ep) cox_term(rhs,term,ep)["p"]
t1_C_DFS <- max(Cval("pN_ord","pN_ord",EP$DFS), Cval("pT_ord","pT_ord",EP$DFS), na.rm=TRUE)
t1_C_OS  <- max(Cval("pN_ord","pN_ord",EP$OS),  Cval("pT_ord","pT_ord",EP$OS),  na.rm=TRUE)
t1_sig_DFS <- min(pval("pN_ord","pN_ord",EP$DFS), pval("pT_ord","pT_ord",EP$DFS), na.rm=TRUE) < 0.05
t1_sig_OS  <- min(pval("pN_ord","pN_ord",EP$OS),  pval("pT_ord","pT_ord",EP$OS),  na.rm=TRUE) < 0.05
ok_DFS <- (t1_C_DFS >= 0.60) && t1_sig_DFS
ok_OS  <- (t1_C_OS  >= 0.60) && t1_sig_OS
verdict <- if (ok_DFS) "TIER-1 DETECTS under DFS -> cohort can see signal -> R5/R6 DFS nulls INTERPRETABLE -> publishable negative (DFS)." else
           if (ok_OS)  "TIER-1 FLAT under DFS but DETECTS under OS -> WRONG ENDPOINT: R5/R6 (DFS) nulls may be artifacts; flag OS re-analysis (decide separately)." else
           "TIER-1 FLAT under BOTH endpoints -> cohort UNDERPOWERED / cannot support the question -> Meyer result is INCONCLUSIVE, not a biological negative."
til_DFS <- (Cval("lymph_z","lymph_z",EP$DFS)>=0.58) && (pval("lymph_z","lymph_z",EP$DFS)<0.05)
til_OS  <- (Cval("lymph_z","lymph_z",EP$OS) >=0.58) && (pval("lymph_z","lymph_z",EP$OS) <0.05)
tier2 <- if (til_DFS || til_OS) "Tier-2 total-TIL IS prognostic but spatial features were not -> the negative is specifically SPATIAL-BEYOND-COMPOSITION (clean claim)." else
         "Tier-2 total-TIL is ALSO flat -> the immune compartment is quiet in this cohort (broader caveat on the negative)."

## ---- report ----
out <- paste0(
"# R7 COHORT VALIDITY / POSITIVE CONTROLS (Meyer, full cohort, both endpoints; PENDING REVIEW)\n\n",
"## Endpoint follow-up / censoring\n",
"DFS: ", ep_summary(EP$DFS), "\n",
"OS : ", ep_summary(EP$OS),  "\n\n",
"## Positive controls -- Cox (HR per level/unit/SD), 95% CI, C-index\n",
sprintf("%-26s | %-34s | %-34s\n", "control", "DFS", "OS"),
paste(sapply(rows, function(r) sprintf("%-26s | %-34s | %-34s", r["name"], r["DFS"], r["OS"])), collapse="\n"),
"\n\n## Nodal pN detail\n",
"DFS strata (n / events):\n", paste(capture.output(print(pnD$strata)), collapse="\n"),
sprintf("\nDFS KM log-rank p = %.4g\n", pnD$logrank_p),
"DFS per-level Cox:\n", paste(capture.output(print(round(pnD$coef,3))), collapse="\n"),
"\nOS strata (n / events):\n", paste(capture.output(print(pnO$strata)), collapse="\n"),
sprintf("\nOS KM log-rank p = %.4g\n", pnO$logrank_p),
"OS per-level Cox:\n", paste(capture.output(print(round(pnO$coef,3))), collapse="\n"),
"\n\n## Notes\n",
"- Grade is grade-3-dominant in Meyer -> low variance; a flat grade HR is a variance issue, NOT cohort failure.\n",
sprintf("- Tier-1 best C: DFS=%.3f (sig=%s), OS=%.3f (sig=%s). Threshold for 'detects' = C>=0.60 & p<0.05.\n",
        t1_C_DFS, t1_sig_DFS, t1_C_OS, t1_sig_OS),
"\n## VERDICT (pre-specified)\n", verdict, "\n", tier2, "\n",
"\nHARD STOP: Meyer-only. Do NOT write the negative / re-run spatial under OS / touch METABRIC (reviewer call).\n")
emit(out, "fork2_R7_cohort_validity.md")
cat("\nDONE R7. Artifact: results/fork2_R7_cohort_validity.md\n")
