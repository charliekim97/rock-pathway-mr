import pandas as pd, numpy as np
api=pd.read_csv("api_outcomes_long.tsv",sep="\t")
lf=pd.read_csv("liverfat_harmonization_preview.txt",sep="\t")
bmi=pd.read_csv("bmi_rock2_pulit.tsv",sep="\t"); bmi["rs"]=bmi["SNP"].str.replace(r":.*","",regex=True)
conv=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")

print("=== PALINDROMIC AF alignment ===")
print("rs4335920: exposure EA=A freq 0.333(GTEx)/0.353(EUR)  [muscle]")
print("rs7598570: exposure EA=%s freq %.3f  [blood, eQTLGen]"%(
   conv[conv.SNP=='rs7598570'].AssessedAllele.iloc[0], conv[conv.SNP=='rs7598570'].eaf.iloc[0]))
def show(rs):
    print(f"\n  {rs}:")
    for o in ["FI_EUR","FG_EUR","HbA1c_EUR","2hGlu_EUR","T2D_Xue2018"]:
        r=api[(api.outcome==o)&(api.SNP==rs)]
        if len(r) and pd.notna(r.EA.iloc[0]):
            print(f"    {o:12s} EA={r.EA.iloc[0]} freq={r.eaf.iloc[0]}")
    r=lf[lf.SNP==rs]
    if len(r): print(f"    {'HepaticFat':12s} EA={r.out_EA.iloc[0]} freq={r.out_eaf.iloc[0]:.3f}")
    r=bmi[bmi.rs==rs]
    if len(r): print(f"    {'BMI_Pulit':12s} EA(A1)={r.Allele1.iloc[0].upper()} freq={r.Freq1.iloc[0]}")
show("rs4335920"); show("rs7598570")

print("\n=== SIGNAL PREVIEW: lead + muscle single-SNP p across outcomes ===")
for rs,lab in [("rs12468344","blood lead"),("rs4335920","muscle"),("rs12996712","muscle")]:
    line=f"  {rs} ({lab}): "
    for o in ["FI_EUR","FG_EUR","HbA1c_EUR","2hGlu_EUR","T2D_Xue2018"]:
        r=api[(api.outcome==o)&(api.SNP==rs)]
        if len(r) and pd.notna(r.p.iloc[0]): line+=f"{o.split('_')[0]}={r.p.iloc[0]:.2g}  "
    r=lf[lf.SNP==rs]; 
    if len(r): line+=f"HepFat={r.out_p.iloc[0]:.2g}  "
    print(line)
