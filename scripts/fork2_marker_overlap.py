"""
Fork 2 -- PROTEIN-MARKER overlap, compartment-split, for each gate-clean survivor cohort.

NOT the Step-0 gene-level logic. Fork 2 has no RNA. This intersects ANTIBODY TARGETS measured
in BOTH the discovery cohort AND METABRIC-IMC, canonicalized at the marker/epitope level via the
explicit override table below (CD8a=CD8; ER=estrogen receptor; HER2=c-erbB-2; panCK=Pan-Keratin;
combined channels expand: CK8-18->{CK8,CK18}, CD31-vWF->{CD31,vWF}, c-Caspase3c-PARP->{CASP3,PARP}).
Then it splits the SHARED markers into IMMUNE / STROMAL / EPITHELIAL-TUMOR / OTHER via a documented
compartment table. NO recommendation -- the thesis axis is open; this surfaces the lists so the
biology choice can be read off them.

All panels are parsed from real source files in data/panels/ (see PROVENANCE.md):
  metabric_markers.txt                  METABRIC-IMC 37 (validation target)
  basel_zuri_stainingpanel_RAW.csv      Jackson/Basel IMC (Zenodo 3518284)
  keren_mibi_panel_RAW.txt              Keren MIBI-TNBC (cellData.csv header)
  engelhardt_cycif_proteins_RAW.csv     Engelhardt/Chang CycIF (GitHub engjen/cycIF_TMAs)
  meyer2025_imc_panel_RAW.tsv           Meyer 2025 IMC (Zenodo 10890543 raw headers; CK clone-walled)

HARD RULE: marker facts from these real files only; CK identities that are clone-truncated in the
Meyer deposit are NOT counted as strict matches (reported provisional). Output: PENDING REVIEW.
"""

from __future__ import annotations
import sys, csv
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PATHS  # noqa: E402

# --- Canonicalization: raw antibody label -> set of canonical marker tokens. ---
# Only entries that can match METABRIC matter; anything unmapped becomes its own token and
# simply won't intersect. Combined channels expand to multiple tokens.
OVERRIDE = {
    # combined / multi
    "ck8-18": {"CK8", "CK18"}, "cytokeratin 8/18": {"CK8", "CK18"},
    "cd31-vwf": {"CD31", "VWF"},
    "c-caspase3c-parp": {"CASP3", "PARP"},
    "cleaved caspase3 / cleaved parp (apoptosis)": {"CASP3", "PARP"},
    # cytokeratins (identities known -> specific token; only METABRIC's CK5/CK8/CK18/panCK match)
    "ck5": {"CK5"}, "cytokeratin 5": {"CK5"}, "ck8": {"CK8"}, "ck18": {"CK18"},
    "ck7": {"CK7"}, "cytokeratin 7": {"CK7"}, "ck19": {"CK19"}, "cytokeratin 19": {"CK19"},
    "ck14": {"CK14"}, "keratin 14": {"CK14"}, "keratin 14 (krt14)": {"CK14"},
    "ck17": {"CK17"}, "keratin17": {"CK17"}, "keratin6": {"CK6"},
    "panck": {"PANCK"}, "pan cytokeratin": {"PANCK"}, "pan-keratin": {"PANCK"},
    "keratin epithelial": {"PANCK"},  # AE3 component
    # receptors / tumor
    "er": {"ER"}, "estrogen receptor (er)": {"ER"}, "estrogen receptor": {"ER"},
    "her2": {"HER2"}, "her2 (3b5)": {"HER2"}, "her2 (d8f12)": {"HER2"},
    "c-erbb-2 - her2": {"HER2"}, "her2 (c-erbb-2)": {"HER2"},
    "ki-67": {"KI67"}, "ki67": {"KI67"},
    # immune
    "cd8": {"CD8"}, "cd8a": {"CD8"}, "foxp3": {"FOXP3"}, "pd1": {"PD1"}, "pd-1": {"PD1"},
    "cd45ro": {"CD45RO"}, "cd45ra": {"CD45RA"}, "hla-dr": {"HLA_DR"},
    "hla_class_1": {"HLA_ABC"}, "hla-abc": {"HLA_ABC"}, "cd11c": {"CD11C"},
    "sma": {"SMA"}, "asma": {"SMA"},
    "podoplanin": {"PDPN"}, "pdpn": {"PDPN"},
    "histone h3": {"HISTONE_H3"},
    "cleaved parp": {"PARP"}, "cleaved caspase3": {"CASP3"}, "cleaved casp3": {"CASP3"},
}
# tokens treated as pan-CD45 (NOT the same epitope as METABRIC's CD45RA/CD45RO isoforms -> excluded)
PAN_CD45 = {"cd45", "cd45 (pan)"}

# Non-marker columns in the Keren cellData.csv header
KEREN_NONMARKER = {"sampleid","celllabelinimage","cellsize","c","na","si","p","ca","fe","background",
                   "ta","au","tumoryn","tumorcluster","group","immunecluster","immunegroup","dsdna"}

# --- Compartment table for METABRIC tokens (shared markers are a subset of these). ---
COMPARTMENT = {
    # immune
    "CD3":"immune","CD4":"immune","CD8":"immune","CD20":"immune","CD68":"immune","CD163":"immune",
    "CD11C":"immune","CD16":"immune","CD15":"immune","CD38":"immune","CD45RA":"immune","CD45RO":"immune",
    "FOXP3":"immune","PD1":"immune","ICOS":"immune","OX40":"immune","GITR":"immune","HLA_DR":"immune",
    "HLA_ABC":"immune","B2M":"immune","CD57":"immune",
    # stromal
    "SMA":"stromal","FSP1":"stromal","PDPN":"stromal","PDGFRB":"stromal","CAV1":"stromal",
    "CD31":"stromal","VWF":"stromal","CXCL12":"stromal",
    # epithelial-tumor
    "CK5":"epithelial-tumor","CK8":"epithelial-tumor","CK18":"epithelial-tumor","PANCK":"epithelial-tumor",
    "HER2":"epithelial-tumor","ER":"epithelial-tumor","KI67":"epithelial-tumor",
    # other (structural / apoptosis)
    "HISTONE_H3":"other","CASP3":"other","PARP":"other",
}

def canon(raw: str) -> set[str]:
    k = raw.strip().lower()
    if k in PAN_CD45:
        return {"CD45_PAN"}              # deliberately non-matching to RA/RO isoforms
    if k in OVERRIDE:
        return set(OVERRIDE[k])
    if "clone-walled" in k:              # ONLY Meyer's clone-truncated CKs are unresolved
        return {"CK_UNRESOLVED"}         # -> provisional, non-strict
    return {raw.strip().upper().replace(" ", "_")}

def tokset(raws):
    out = set()
    for r in raws:
        out |= canon(r)
    return out

def load_metabric():
    raws = [l.strip() for l in (PATHS["panels"]/"metabric_markers.txt").read_text(encoding="utf-8").splitlines()
            if l.strip() and not l.startswith("#") and l.strip() not in ("DNA1","DNA2")]
    return raws, tokset(raws)

def load_basel():
    raws = []
    with (PATHS["panels"]/"basel_zuri_stainingpanel_RAW.csv").open(encoding="utf-8") as f:
        for row in csv.reader(f):
            if len(row) < 3: continue
            t = row[2].strip()
            if t in ("Target","",) or t.lower() in ("argondimers","rutheniumtetroxide","rabbit igg h l",
                    "undefined","dna1","dna2"):
                continue
            raws.append(t)
    return raws, tokset(raws)

def load_keren():
    line = [l for l in (PATHS["panels"]/"keren_mibi_panel_RAW.txt").read_text(encoding="utf-8").splitlines()
            if l.strip() and not l.startswith("#")][0]
    raws = [c for c in line.split(",") if c.strip().lower() not in KEREN_NONMARKER]
    return raws, tokset(raws)

def load_engelhardt():
    raws = []
    txt = (PATHS["panels"]/"engelhardt_cycif_proteins_RAW.csv").read_text(encoding="utf-8")
    for cell in txt.replace("\n", ",").split(","):
        cell = cell.strip()
        if not cell: continue
        parts = cell.split()           # entries look like "0 ER" / "21 CK14"
        name = parts[-1] if parts else ""
        if name and not name.isdigit():
            raws.append(name)
    return raws, tokset(raws)

def load_meyer():
    raws = []
    with (PATHS["panels"]/"meyer2025_imc_panel_RAW.tsv").open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("metal\t"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3: continue
            t = cols[2].strip()
            if t in ("DNA1","DNA2"): continue
            raws.append(t)
    return raws, tokset(raws)

def main():
    results = PATHS["results_dir"]; results.mkdir(parents=True, exist_ok=True)
    mb_raw, mb = load_metabric()
    cohorts = {
        "Jackson/Basel IMC": (load_basel(), "IMC", "281 (+72)", "OS+DFS", "same as METABRIC"),
        "Keren MIBI-TNBC":   (load_keren(), "MIBI", "41 (TNBC)", "OS+recurrence", "cross-platform"),
        "Engelhardt CycIF":  (load_engelhardt(), "CycIF", "102", "OS+RFS", "cross-platform"),
        "Meyer 2025 IMC":    (load_meyer(), "IMC", "215 (TNBC)", "OS+DFS/recurrence", "same as METABRIC"),
    }
    KEY = ["CD8","FOXP3","ER","HER2","CD3","PANCK","CD31","HLA_DR"]
    out = ["# Fork-2 PROTEIN-MARKER overlap with METABRIC-IMC -- compartment-split (PENDING REVIEW)\n",
           f"METABRIC-IMC validation panel = {len(mb_raw)} markers -> {len(mb)} canonical tokens.",
           "Shared = antibody targets measured in BOTH (canonicalized; clone-walled Meyer CKs = provisional).",
           "Pan-CD45 is NOT matched to METABRIC's CD45RA/CD45RO isoforms (different epitope).\n"]
    summary = []
    for name,((raws,toks),plat,n,endpoint,platnote) in cohorts.items():
        shared = sorted(toks & mb)
        prov_ck = "CK_UNRESOLVED" in toks   # Meyer clone-walled CKs (provisional epithelial overlap)
        comp = {"immune":[],"stromal":[],"epithelial-tumor":[],"other":[]}
        for t in shared:
            comp[COMPARTMENT.get(t,"other")].append(t)
        out.append(f"## {name}  ({plat}, n={n}, endpoint={endpoint}; platform {platnote})")
        out.append(f"- panel size: {len(toks)} canonical tokens; **shared with METABRIC: {len(shared)}**")
        for c in ("immune","stromal","epithelial-tumor","other"):
            out.append(f"  - {c} ({len(comp[c])}): {', '.join(comp[c]) if comp[c] else '-'}")
        if prov_ck:
            out.append("  - PROVISIONAL: Meyer has >=1 clone-walled cytokeratin channel; if any is CK5/CK8-18/"
                       "panCK (likely by metal layout, UNCONFIRMED) add to epithelial-tumor.")
        disc = {k:("YES" if (canon_tok:=k) in toks else "no") for k in KEY}
        out.append(f"  - key discriminators: " + ", ".join(f"{k}={disc[k]}" for k in KEY))
        out.append("")
        summary.append((name,plat,len(shared),len(comp["immune"]),len(comp["stromal"]),
                        len(comp["epithelial-tumor"]),disc["CD8"],disc["FOXP3"],disc["ER"],disc["HER2"]))
    out.append("## Side-by-side (shared counts by compartment + key discriminators)\n")
    out.append("| cohort | platform | shared | immune | stromal | epi-tumor | CD8 | FOXP3 | ER | HER2 |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in summary:
        out.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(*r))
    (results/"fork2_marker_overlap.md").write_text("\n".join(out)+"\n", encoding="utf-8")
    print("\n".join(out[-(len(summary)+4):]))
    print("\nwritten:", results/"fork2_marker_overlap.md")

if __name__ == "__main__":
    main()
