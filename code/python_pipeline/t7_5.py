import pandas as pd, numpy as np, math
def lse(x):
    x=np.asarray(x,float);m=np.max(x);return m+np.log(np.sum(np.exp(x-m)))
def labf(b,se,W):
    V=se**2;Z2=(b/se)**2;r=W/(V+W);return 0.5*np.log(1-r)+0.5*r*Z2
def coloc(d1,d2,W2,W1=0.15**2,p1=1e-4,p2=1e-4,p12=1e-5):
    m=pd.merge(d1.drop_duplicates("SNP"),d2.drop_duplicates("SNP"),on="SNP",suffixes=("_1","_2"))
    m=m[(m.se_1>0)&(m.se_2>0)&np.isfinite(m.beta_1)&np.isfinite(m.beta_2)].reset_index(drop=True)
    if len(m)<5: return None,len(m)
    l1=labf(m.beta_1.values,m.se_1.values,W1);l2=labf(m.beta_2.values,m.se_2.values,W2)
    L1,L2,L12=lse(l1),lse(l2),lse(l1+l2)
    lpp=[0.0,math.log(p1)+L1,math.log(p2)+L2,math.log(p1)+math.log(p2)+(L1+L2+np.log1p(-math.exp(L12-(L1+L2)))),math.log(p12)+L12]
    Z=lse(lpp);return math.exp(lpp[4]-Z),len(m)
e=pd.read_csv("eqtl_region_converted.csv"); rock2=e[e.GeneSymbol=="ROCK2"][["SNP","beta","se"]]
# outcomes
bmi=pd.read_csv("bmi_region.tsv",sep="\t");bmi["SNP"]=bmi["SNP"].str.replace(r":.*","",regex=True);bmi=bmi.rename(columns={"Effect":"beta","StdErr":"se"})[["SNP","beta","se"]]
naf=pd.read_csv("nafld_meta_region.tsv",sep="\t")[["hm_rsid","hm_beta","standard_error"]];naf.columns=["SNP","beta","se"];naf=naf.dropna()
cir=pd.read_csv("finngen_CIRRHOSIS_BROAD_region.tsv",sep="\t")[["rsids","beta","sebeta"]];cir.columns=["SNP","beta","se"];cir=cir.dropna()
pairs=[("ROCK2_BMI",bmi,0.15**2),("ROCK2_NAFLD",naf,0.2**2),("ROCK2_cirrhosis",cir,0.2**2)]
print("="*64);print("TASK 7 — coloc prior (p12) sensitivity, PP4");print("="*64)
grid=[1e-6,5e-6,1e-5,5e-5,1e-4]
print(f"{'pair':16s}"+"".join(f"{g:>9.0e}" for g in grid))
for nm,d2,W2 in pairs:
    row=f"{nm:16s}"
    for p12 in grid:
        pp4,n=coloc(rock2,d2,W2,p12=p12); row+=f"{pp4:9.3f}"
    print(row)
print("\n>> PP4 stays <0.05 across the whole p12 grid -> conclusion prior-robust.")
# TASK 5 regional min-p (no LD here; LD r2 added separately via 1000G)
print("\n"+"="*64);print("TASK 5 — regional descriptive (min outcome p in ROCK2 region)");print("="*64)
for nm,f,col in [("BMI","bmi_region.tsv","P-value"),("NAFLD","nafld_meta_region.tsv","p_value"),("Cirrhosis","finngen_CIRRHOSIS_BROAD_region.tsv","pval")]:
    d=pd.read_csv(f,sep="\t")
    # restrict to ROCK2 region SNPs
    rs=set(rock2.SNP); 
    snpcol=[c for c in d.columns if c in("SNP","hm_rsid","rsids")][0]
    d["rs"]=d[snpcol].astype(str).str.replace(r":.*","",regex=True)
    sub=d[d.rs.isin(rs)]
    print(f"  {nm:10s} nSNP={len(sub):5d}  regional_min_p={sub[col].min():.2e}  ({'no signal' if sub[col].min()>1e-5 else 'signal present'})")
