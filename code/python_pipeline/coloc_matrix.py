import pandas as pd, numpy as np, math, os
def lse(x):
    x=np.asarray(x,float);m=np.max(x);return m+np.log(np.sum(np.exp(x-m)))
def labf(b,se,W):
    V=se**2;Z2=(b/se)**2;r=W/(V+W);return 0.5*np.log(1-r)+0.5*r*Z2
def coloc(d1,d2,W2,W1=0.15**2,p1=1e-4,p2=1e-4,p12=1e-5):
    d1=d1.drop_duplicates("SNP");d2=d2.drop_duplicates("SNP")
    m=d1.merge(d2,on="SNP",suffixes=("_1","_2")).reset_index(drop=True)
    m=m[(m.se_1>0)&(m.se_2>0)&np.isfinite(m.beta_1)&np.isfinite(m.beta_2)&np.isfinite(m.se_1)&np.isfinite(m.se_2)].reset_index(drop=True)
    if len(m)<5:return None
    l1=labf(m.beta_1.values,m.se_1.values,W1);l2=labf(m.beta_2.values,m.se_2.values,W2)
    L1,L2,L12=lse(l1),lse(l2),lse(l1+l2)
    lpp=[0.0,math.log(p1)+L1,math.log(p2)+L2,math.log(p1)+math.log(p2)+(L1+L2+np.log1p(-math.exp(L12-(L1+L2)))),math.log(p12)+L12]
    Z=lse(lpp);pp=[math.exp(x-Z) for x in lpp]
    return dict(n=len(m),PP3=pp[3],PP4=pp[4])
e=pd.read_csv("rho_eqtl_region_conv.csv")
def eq(g):return e[e.GeneSymbol==g][["SNP","beta","se"]]
def apifile(gene,oc):
    L={"ROCK1":"21000000","RHOA":"48900000"}[gene]
    for cand in [f"out_{gene}_{oc}_{L}.csv", f"out_{gene}_{oc}.csv"]:
        if os.path.exists(cand):
            d=pd.read_csv(cand); return d[["SNP","beta","se"]]
    return None
bmi=pd.read_csv("rho_bmi_region.tsv",sep="\t");bmi["SNP"]=bmi["SNP"].str.replace(r":.*","",regex=True);bmi=bmi[["SNP","Effect","StdErr"]];bmi.columns=["SNP","beta","se"]
hf=pd.read_csv("rho_hf_region.tsv",sep="\t")[["hm_rsid","hm_beta","standard_error"]];hf.columns=["SNP","beta","se"]
Q=0.15**2;CC=0.2**2
outs=[("HepaticFat",hf,Q),("FI",None,Q),("FG",None,Q),("HbA1c",None,Q),("2hGlu",None,Q),("T2D",None,CC),("BMI",bmi,Q)]
res={}
for g in ["RHOA","ROCK1"]:
    res[g]={}
    for oc,d2,W2 in outs:
        if d2 is None: d2=apifile(g,oc)
        r=coloc(eq(g),d2,W2) if d2 is not None else None
        res[g][oc]=r
        print(f"{g:6s} {oc:11s} -> {('PP4=%.3f PP3=%.3f n=%d'%(r['PP4'],r['PP3'],r['n'])) if r else 'NA'}")
# ROCK2 from block4
r2=pd.read_csv("coloc_results.csv")
r2map={"BMI":"BMI","HepaticFat":"HepaticFat","FI":"FI","FG":"FG","HbA1c":"HbA1c","2hGlu":"2hGlu","T2D":"T2D"}
print("\n===== FULL PP4 MATRIX (3 gene x 7 outcome) =====")
ocs=["HepaticFat","FI","FG","HbA1c","2hGlu","T2D","BMI"]
print(f"{'gene':7s}"+"".join(f"{o:>11s}" for o in ocs))
# ROCK2 row
row="ROCK2  "
for o in ocs:
    v=r2[(r2.eQTL=="ROCK2")&(r2.outcome==o)]
    row+=f"{(v.PP4.iloc[0] if len(v) else float('nan')):>11.3f}"
print(row)
for g in ["ROCK1","RHOA"]:
    row=f"{g:7s}"
    for o in ocs:
        r=res[g].get(o); row+=f"{(r['PP4'] if r else float('nan')):>11.3f}"
    print(row)
