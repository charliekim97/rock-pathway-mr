import pandas as pd, numpy as np
snps=[l.strip() for l in open("instrument_snps_for_outcome.txt") if l.strip()]
lead="rs12468344"; r201=["rs12468344","rs13034888","rs7574150"]
r21=pd.read_csv("blood_clumped_independent.txt",sep="\t")["SNP"].tolist()
muscle=["rs4335920","rs12996712"]
tiers={"Blood lead (Primary)":[lead],"Blood r2<.01 (Support IVW)":r201,
       "Blood r2<.1 (LD-corr sens)":r21,"Muscle 2-SNP":muscle}

found={}
# API outcomes
api=pd.read_csv("api_outcomes_long.tsv",sep="\t")
for o in api.outcome.unique():
    found[o]=set(api[(api.outcome==o)&(api.beta.notna())].SNP)
# hepatic fat (all 17 present)
lf=pd.read_csv("liverfat_harmonization_preview.txt",sep="\t")
found["HepaticFat"]=set(lf.SNP)
# BMI pulit
bmi=pd.read_csv("bmi_rock2_pulit.tsv",sep="\t")
found["BMI_Pulit"]=set(bmi["SNP"].str.replace(r":.*","",regex=True))
# HOMA-IR
try:
    hm=pd.read_csv("homair_rock2.tsv",sep="\t");found["HOMA-IR_Dupuis"]=set(hm.iloc[:,0])
except: found["HOMA-IR_Dupuis"]=set()

order=["HepaticFat","FI_EUR","FG_EUR","HbA1c_EUR","2hGlu_EUR","T2D_Xue2018","BMI_Pulit","HOMA-IR_Dupuis"]
print("=== OUTCOME × TIER coverage matrix ===\n")
hdr="%-26s"%"Tier"+"".join("%-14s"%o[:13] for o in order)
print(hdr); print("-"*len(hdr))
for tn,ts in tiers.items():
    row="%-26s"%tn
    for o in order:
        c=len(set(ts)&found.get(o,set())); row+="%-14s"%f"{c}/{len(set(ts))}"
    print(row)
tot="%-26s"%"TOTAL unique /17"
for o in order: tot+="%-14s"%f"{len(found.get(o,set()))}/17"
print("-"*len(hdr)); print(tot)
