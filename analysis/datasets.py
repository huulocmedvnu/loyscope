# the canonical set of datasets: each cohort counted once (OneK1K through the author count matrices; the library
# regrouping of the retina atlas is a robustness analysis, not a dataset). Used by every summary and figure script.
import os
EXCLUDE = {"onek1k", "retina_lib"}
NAMES = {"aida": "AIDA", "onek1k_geo": "OneK1K", "iaging": "Immune Aging", "blood": "Indonesia blood", "covid": "COVID-19 blood",
         "lupus": "Lupus blood", "bone_marrow": "Bone marrow", "kidney": "Kidney KPMP cells", "ckd_sn": "Kidney KPMP nuclei",
         "gb2024": "Kidney, published study", "heart": "Heart", "dcm": "Heart DCM", "encode_heart": "Heart ENCODE", "hlca": "Lung HLCA",
         "ild": "Lung ILD", "muscle": "Skeletal muscle", "pancrea": "Pancreas", "adipose": "White adipose", "skin": "Skin",
         "retina": "Retina/choroid", "ad_mssm": "Brain MSSM", "ad_radc": "Brain RADC", "ad_microglia": "Microglia",
         "seaad_mtg": "Brain MTG", "tsapiens": "Tabula Sapiens"}
def canonical(extract_dir):
    """dataset names in NAMES order that have a fit in extract_dir"""
    return [k for k in NAMES if os.path.exists(os.path.join(extract_dir, f"{k}.fitamb.tsv"))]
