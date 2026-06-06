"""
Build filter1_age_fixed experiment CSVs.

Fixed features (always in model):
  ecog_z  — ECOG performance status, z-normalized (closest to sequencing day 0)
  age_z   — age at sequencing, z-normalized
  SAMPLE_TYPE_Local_Recurrence, SAMPLE_TYPE_Metastasis
  GENE_PANEL_IMPACT341, GENE_PANEL_IMPACT410, GENE_PANEL_IMPACT505
  CRC only: site_rectal, site_colorectal_NOS  (colon = reference)

Removed from candidate pool vs original filter1:
  FGA_z, ANEUPLOIDY_SCORE_z  (would absorb CNA signal → over-adjustment)
  TMB_log_z       (reserved for sensitivity analysis, not in main model)
  TUMOR_PURITY_z  (technical covariate, not in main model)

Missing ECOG: impute with training set median per cancer.
Missing age:  impute with training set median per cancer.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import os

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR     = ROOT / "model_input"
FILTER1_DIR  = DATA_DIR / "filter1_data"
ID_MAP_DIR   = DATA_DIR / "sample_id_mapping"
CLINICAL_DIR = ROOT / "raw_data" / "all_timeline_csv"
OUT_DIR      = DATA_DIR / "filter1_age_fixed"

os.makedirs(OUT_DIR, exist_ok=True)

CANCER_CFG = {
    "BRCA": {"begin": 18, "end": 8},
    "LUAD": {"begin": 28, "end": 16},
    "PAAD": {"begin": 26, "end": 15},
    "PRAD": {"begin": 20, "end": 8},
    "CRC":  {"begin": 28, "end": 14},
}


def load_clinical_base():
    patient = pd.read_csv(
        f"{CLINICAL_DIR}/data_clinical_patient.csv",
        usecols=["PATIENT_ID", "CURRENT_AGE_DEID", "GENDER"]
    ).rename(columns={"CURRENT_AGE_DEID": "age"})
    sample = pd.read_csv(
        f"{CLINICAL_DIR}/data_clinical_sample.csv",
        usecols=["SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE_DETAILED"]
    )
    perf = pd.read_csv(
        f"{CLINICAL_DIR}/data_timeline_performance_status.txt", sep="\t",
        usecols=["PATIENT_ID", "START_DATE", "ECOG"]
    )
    return patient, sample, perf


def get_closest_ecog(perf_df, patient_ids):
    """
    For each patient, take the ECOG record closest to day 0 (sequencing).
    No window restriction — use any available record.
    """
    sub = perf_df[
        perf_df["PATIENT_ID"].isin(patient_ids) & perf_df["ECOG"].notna()
    ].copy()
    sub["abs_date"] = sub["START_DATE"].abs()
    closest = (sub.sort_values("abs_date")
                  .groupby("PATIENT_ID")["ECOG"]
                  .first()
                  .reset_index())
    return closest  # PATIENT_ID, ECOG


def get_patient_ids(cancer, sample_df):
    """Map sample IDs (train+test) to patient IDs."""
    train_ids = pd.read_csv(f"{ID_MAP_DIR}/{cancer}_train_ids.csv")["SAMPLE_ID"].values
    test_ids  = pd.read_csv(f"{ID_MAP_DIR}/{cancer}_test_ids.csv")["SAMPLE_ID"].values
    id_to_pat = sample_df.set_index("SAMPLE_ID")["PATIENT_ID"].to_dict()
    return train_ids, test_ids, id_to_pat


CRC_SITE_MAP = {
    "Colon Adenocarcinoma":      "colon",
    "Rectal Adenocarcinoma":     "rectal",
    "Colorectal Adenocarcinoma": "colorectal_NOS",
}


def attach_clinical(df, ids, id_to_pat, patient_df, ecog_df, sample_df=None, cancer=None):
    """Attach age, ECOG (and CRC site_group) by row-order SAMPLE_IDs."""
    n = min(len(df), len(ids))
    df = df.iloc[:n].copy()
    sample_ids = ids[:n]

    pat_ids = [id_to_pat.get(s, None) for s in sample_ids]
    df["_PAT"] = pat_ids
    df["_SID"] = sample_ids

    age_map  = patient_df.set_index("PATIENT_ID")["age"].to_dict()
    ecog_map = ecog_df.set_index("PATIENT_ID")["ECOG"].to_dict()

    df["age"]  = df["_PAT"].map(age_map)
    df["ecog"] = df["_PAT"].map(ecog_map)

    if cancer == "CRC" and sample_df is not None:
        ctd_map = sample_df.set_index("SAMPLE_ID")["CANCER_TYPE_DETAILED"].to_dict()
        site_raw = df["_SID"].map(ctd_map).map(CRC_SITE_MAP).fillna("colorectal_NOS")
        df["site_rectal"]         = (site_raw == "rectal").astype(int)
        df["site_colorectal_NOS"] = (site_raw == "colorectal_NOS").astype(int)

    df.drop(columns=["_PAT", "_SID"], inplace=True)
    return df


def build_csv(cancer, patient_df, sample_df, perf_df):
    print(f"\n{'='*55}\n  {cancer}\n{'='*55}")

    train_ids, test_ids, id_to_pat = get_patient_ids(cancer, sample_df)
    all_patient_ids = list({id_to_pat[s] for s in list(train_ids) + list(test_ids)
                            if s in id_to_pat})

    # Closest ECOG per patient
    ecog_df = get_closest_ecog(perf_df, all_patient_ids)

    # Load filter1 CSVs
    train = pd.read_csv(f"{FILTER1_DIR}/{cancer}_3flag_enh_filter_train.csv")
    test  = pd.read_csv(f"{FILTER1_DIR}/{cancer}_3flag_enh_filter_test.csv")

    # Attach clinical (age, ecog, CRC site) by row order
    train = attach_clinical(train, train_ids, id_to_pat, patient_df, ecog_df, sample_df, cancer)
    test  = attach_clinical(test,  test_ids,  id_to_pat, patient_df, ecog_df, sample_df, cancer)

    # Z-normalize using TRAIN stats only
    for col in ["age", "ecog"]:
        mu  = train[col].mean(skipna=True)
        std = train[col].std(skipna=True)
        std = std if std > 0 else 1.0
        median_train = train[col].median(skipna=True)

        train[f"{col}_z"] = (train[col] - mu) / std
        test[f"{col}_z"]  = (test[col]  - mu) / std

        # Impute missing with 0 (= train mean after z-norm)
        n_missing_tr = train[f"{col}_z"].isna().sum()
        n_missing_te = test[f"{col}_z"].isna().sum()
        train[f"{col}_z"] = train[f"{col}_z"].fillna(0)
        test[f"{col}_z"]  = test[f"{col}_z"].fillna(0)

        cov = (len(train) - n_missing_tr) / len(train) * 100
        print(f"  {col}: train coverage={cov:.1f}%, missing imputed={n_missing_tr} "
              f"(median_raw={median_train:.1f})")

        train.drop(columns=[col], inplace=True)
        test.drop(columns=[col], inplace=True)

    # Define column groups
    fixed_orig = [c for c in train.columns
                  if c.startswith("SAMPLE_TYPE_") or c.startswith("GENE_PANEL_")]
    site_fixed = ["site_rectal", "site_colorectal_NOS"] if cancer == "CRC" else []
    remove     = {"FGA_z", "ANEUPLOIDY_SCORE_z", "TMB_log_z", "TUMOR_PURITY_z"}
    survival   = ["Event_OS", "OS_MONTHS"]
    new_fixed  = ["ecog_z", "age_z"]
    all_fixed  = new_fixed + fixed_orig + site_fixed
    skip       = set(fixed_orig) | set(site_fixed) | remove | set(survival) | set(new_fixed)
    candidates = [c for c in train.columns if c not in skip]

    print(f"  Fixed ({len(all_fixed)}): {all_fixed}")
    print(f"  Removed: FGA_z, ANEUPLOIDY_SCORE_z, TMB_log_z, TUMOR_PURITY_z")
    print(f"  Candidates: {len(candidates)}")

    col_order = new_fixed + fixed_orig + site_fixed + candidates + survival
    train_out = train[col_order]
    test_out  = test[col_order]

    train_out.to_csv(f"{OUT_DIR}/{cancer}_train.csv", index=False)
    test_out.to_csv( f"{OUT_DIR}/{cancer}_test.csv",  index=False)
    print(f"  Saved: {cancer}_train.csv {train_out.shape}, {cancer}_test.csv {test_out.shape}")

    return fixed_orig, site_fixed, candidates


def write_cfg(cancer, fixed_orig, site_fixed, candidates):
    cfg = CANCER_CFG[cancer]
    new_fixed = ["ecog_z", "age_z"]
    feature_line = ",".join(new_fixed + fixed_orig + site_fixed)

    cfg_text = f"""[IBCGA_Config]
data = {cancer}_train.csv

begin = {cfg['begin']}
end = {cfg['end']}

cv = 5
iga_type = 3
estimator = CoxPH
fitness = cindex
crossover = 0.8
mutation = 0.05
generation = 100
first_generation = 300
population = 50
conv_gen = 35
first_conv_gen = 100
tol = 0.001
iscg = True
job = 8
mix_cross = True
plot_record = True
pop_record = True
random = -1
resume_path = None
sample = 1
transformer = None
med = True
feature = {feature_line}
"""
    with open(f"{OUT_DIR}/{cancer}_train.cfg", "w") as f:
        f.write(cfg_text)
    print(f"  CFG feature= {feature_line}")
    print(f"  begin={cfg['begin']}, end={cfg['end']}")


def main():
    print("Loading clinical data...")
    patient_df, sample_df, perf_df = load_clinical_base()

    summary = []
    for cancer in ["BRCA", "LUAD", "PAAD", "PRAD", "CRC"]:
        fixed_orig, site_fixed, candidates = build_csv(cancer, patient_df, sample_df, perf_df)
        write_cfg(cancer, fixed_orig, site_fixed, candidates)
        cfg = CANCER_CFG[cancer]
        summary.append(dict(cancer=cancer,
                            n_fixed=2 + len(fixed_orig) + len(site_fixed),
                            n_cand=len(candidates),
                            **cfg))

    print(f"\n\n{'='*65}")
    print("  SUMMARY — filter1_age_fixed (ecog_z + age_z as fixed)")
    print(f"{'='*65}")
    print(f"  {'Cancer':<8} {'Fixed':>6} {'Cand':>6} {'begin':>6} {'end':>5}")
    print(f"  {'-'*34}")
    for r in summary:
        print(f"  {r['cancer']:<8} {r['n_fixed']:>6} {r['n_cand']:>6} "
              f"{r['begin']:>6} {r['end']:>5}")
    print(f"\n  Output: {OUT_DIR}/")


if __name__ == "__main__":
    main()
