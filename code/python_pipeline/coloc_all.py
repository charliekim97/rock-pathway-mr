import pandas as pd, numpy as np, math
def logsumexp(x):
    x=np.asarray(x,float); m=np.max(x); return m+np.log(np.sum(np.exp(x-m)))
def labf(b,se,W):
    V=se**2; Z2=(b/se)**2; r=W/(V+W); return 0.5*np.log(1-r)+0.5*r*Z2
def coloc(d1,d2,W1=0.15**2,W2=0.15**2,p1=1e-4,p2=1e-4,p12=1e-5):
    d1=d1.drop_duplicates("SNP"); d2=d2.drop_duplicates("SNP")
    m=d1.merge(d2,on="SNP",suffixes=("_1","_2")).reset_index(drop=True)
    m=m[(m.se_1>0)&(m.se_2>0)&np.isfinite(m.beta_1)&np.isfinite(m.beta_2)&np.isfinite(m.se_1)&np.isfinite(m.se_2)].reset_index(drop=True)
    if len(m)<5: return None
    l1=labf(m.beta_1.values,m.se_1.values,W1); l2=labf(m.beta_2.values,m.se_2.values,W2)
    L1,L2,L12=logsumexp(l1),logsumexp(l2),logsumexp(l1+l2)
    lpp=[0.0, math.log(p1)+L1, math.log(p2)+L2,
         math.log(p1)+math.log(p2)+(L1+L2+np.log1p(-math.exp(L12-(L1+L2)))),
         math.log(p12)+L12]
    Z=logsumexp(lpp); pp=[math.exp(x-Z) for x in lpp]
    return dict(nsnp=len(m),PP0=pp[0],PP1=pp[1],PP2=pp[2],PP3=pp[3],PP4=pp[4],
                topSNP=m.loc[int((l1+l2).argmax()),"SNP"])
e=pd.read_csv("eqtl_region_converted.csv")
def eqtl(g): return e[e.GeneSymbol==g][["SNP","beta","se"]]
def api(name):
    d=pd.read_csv(f"out_{name}.csv"); return d.rename(columns={"beta":"beta","se":"se"})[["SNP","beta","se"]]
def localf(path,snpcol,bcol,scol,strip=False):
    d=pd.read_csv(path,sep="\t")
    if strip: d[snpcol]=d[snpcol].str.replace(r":.*","",regex=True)
    d=d[[snpcol,bcol,scol]]; d.columns=["SNP","beta","se"]; return d
bmi=localf("bmi_region.tsv","SNP","Effect","StdErr",strip=True)
hf=pd.read_csv("hepfat_region.tsv",sep="\t")[["hm_rsid","hm_beta","standard_error"]]; hf.columns=["SNP","beta","se"]
Q=0.15**2; CC=0.2**2
jobs=[
 ("ROCK2","BMI",bmi,Q),("ROCK2","HepaticFat",hf,Q),
 ("ROCK2","FI",api("ROCK2_FI"),Q),("ROCK2","FG",api("ROCK2_FG"),Q),
 ("ROCK2","HbA1c",api("ROCK2_HbA1c"),Q),("ROCK2","2hGlu",api("ROCK2_2hGlu"),Q),
 ("ROCK2","T2D",api("ROCK2_T2D"),CC),
 ("KCNJ11","T2D[posctrl]",api("KCNJ11_T2D"),CC),
 ("GIPR","T2D[posctrl]",api("GIPR_T2D"),CC),
 ("SORT1","LDL[posctrl]",api("SORT1_LDL"),Q),
]
rows=[]
print(f"{'eQTL':8s} {'outcome':16s} {'nSNP':>5s} {'PP3':>6s} {'PP4':>6s} {'PP4/(PP3+PP4)':>13s}  topSNP")
for g,oc,d2,W2 in jobs:
    r=coloc(eqtl(g),d2,W2=W2)
    if r:
        cond=r['PP4']/(r['PP3']+r['PP4']) if (r['PP3']+r['PP4'])>0 else float('nan')
        print(f"{g:8s} {oc:16s} {r['nsnp']:5d} {r['PP3']:6.3f} {r['PP4']:6.3f} {cond:13.3f}  {r['topSNP']}")
        rows.append([g,oc,r['nsnp'],round(r['PP0'],3),round(r['PP1'],3),round(r['PP2'],3),round(r['PP3'],3),round(r['PP4'],3),r['topSNP']])
pd.DataFrame(rows,columns=["eQTL","outcome","nSNP","PP0","PP1","PP2","PP3","PP4","topSNP"]).to_csv("coloc_results.csv",index=False)
