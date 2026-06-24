import pandas as pd, numpy as np, math
def logsumexp(x):
    x=np.asarray(x); m=np.max(x); return m+np.log(np.sum(np.exp(x-m)))
def labf(beta,se,W):
    V=se**2; Z2=(beta/se)**2; r=W/(V+W)
    return 0.5*np.log(1-r)+0.5*r*Z2
def coloc_abf(d1,d2,W1=0.15**2,W2=0.15**2,p1=1e-4,p2=1e-4,p12=1e-5):
    d1=d1.drop_duplicates("SNP"); d2=d2.drop_duplicates("SNP")
    m=d1.merge(d2,on="SNP",suffixes=("_1","_2")).reset_index(drop=True)
    m=m[(m.se_1>0)&(m.se_2>0)&np.isfinite(m.beta_1)&np.isfinite(m.beta_2)&np.isfinite(m.se_1)&np.isfinite(m.se_2)]
    if len(m)<5: return None
    l1i=labf(m.beta_1.values,m.se_1.values,W1)
    l2i=labf(m.beta_2.values,m.se_2.values,W2)
    L1=logsumexp(l1i); L2=logsumexp(l2i); L12=logsumexp(l1i+l2i)
    # H3 = sum_{i!=j} = exp(L1+L2)-exp(L12)
    lpp=[0.0, math.log(p1)+L1, math.log(p2)+L2,
         math.log(p1)+math.log(p2)+ (L1+L2 + np.log1p(-np.exp(L12-(L1+L2)))),
         math.log(p12)+L12]
    Z=logsumexp(lpp); pp=[math.exp(x-Z) for x in lpp]
    return dict(nsnp=len(m),PP0=pp[0],PP1=pp[1],PP2=pp[2],PP3=pp[3],PP4=pp[4],
                top=m.loc[(l1i+l2i).argmax(),"SNP"])
# eQTL ROCK2
e=pd.read_csv("eqtl_region_converted.csv")
rock2=e[e.GeneSymbol=="ROCK2"][["SNP","beta","se"]].copy()
# BMI
bmi=pd.read_csv("bmi_region.tsv",sep="\t"); bmi["SNP"]=bmi["SNP"].str.replace(r":.*","",regex=True)
bmi=bmi.rename(columns={"Effect":"beta","StdErr":"se"})[["SNP","beta","se"]]
# HepFat
hf=pd.read_csv("hepfat_region.tsv",sep="\t")[["hm_rsid","hm_beta","standard_error"]]; hf.columns=["SNP","beta","se"]
hf=hf.dropna()
for name,d2,W2 in [("BMI",bmi,0.15**2),("HepaticFat",hf,0.15**2)]:
    r=coloc_abf(rock2,d2,W2=W2)
    print(f"ROCK2(blood) x {name}: nsnp={r['nsnp']} PP3={r['PP3']:.3f} PP4={r['PP4']:.3f} (top {r['top']})")
