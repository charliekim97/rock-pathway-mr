import pandas as pd, numpy as np, math
def lse(x):
    x=np.asarray(x,float);m=np.max(x);return m+np.log(np.sum(np.exp(x-m)))
def labf(b,se,W):
    V=se**2;Z2=(b/se)**2;r=W/(V+W);return 0.5*np.log(1-r)+0.5*r*Z2
def coloc(d1,d2,W1=0.15**2,W2=0.15**2,p1=1e-4,p2=1e-4,p12=1e-5):
    d1=d1.drop_duplicates("SNP");d2=d2.drop_duplicates("SNP")
    m=d1.merge(d2,on="SNP",suffixes=("_1","_2")).reset_index(drop=True)
    m=m[(m.se_1>0)&(m.se_2>0)&np.isfinite(m.beta_1)&np.isfinite(m.beta_2)&np.isfinite(m.se_1)&np.isfinite(m.se_2)].reset_index(drop=True)
    if len(m)<5:return None
    l1=labf(m.beta_1.values,m.se_1.values,W1);l2=labf(m.beta_2.values,m.se_2.values,W2)
    L1,L2,L12=lse(l1),lse(l2),lse(l1+l2)
    lpp=[0.0,math.log(p1)+L1,math.log(p2)+L2,
         math.log(p1)+math.log(p2)+(L1+L2+np.log1p(-math.exp(L12-(L1+L2)))),math.log(p12)+L12]
    Z=lse(lpp);pp=[math.exp(x-Z) for x in lpp]
    return dict(n=len(m),PP3=pp[3],PP4=pp[4],top=m.loc[int((l1+l2).argmax()),"SNP"])
e=pd.read_csv("rho_eqtl_region_conv.csv")
def eq(g):return e[e.GeneSymbol==g][["SNP","beta","se"]]
bmi=pd.read_csv("rho_bmi_region.tsv",sep="\t");bmi["SNP"]=bmi["SNP"].str.replace(r":.*","",regex=True)
bmi=bmi[["SNP","Effect","StdErr"]];bmi.columns=["SNP","beta","se"]
hf=pd.read_csv("rho_hf_region.tsv",sep="\t")[["hm_rsid","hm_beta","standard_error"]];hf.columns=["SNP","beta","se"]
print(f"{'eQTL':6s} {'outcome':10s} {'n':>5s} {'PP3':>6s} {'PP4':>6s}  top")
for g in ["RHOA","ROCK1"]:
    for oc,d2 in [("HepaticFat",hf),("BMI",bmi)]:
        r=coloc(eq(g),d2)
        if r: print(f"{g:6s} {oc:10s} {r['n']:5d} {r['PP3']:6.3f} {r['PP4']:6.3f}  {r['top']}")
