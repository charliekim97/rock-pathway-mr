import pandas as pd, numpy as np, math
from statistics import NormalDist
nd=NormalDist()
def z2p(z): return 2*(1-nd.cdf(abs(z)))
def lse(x):
    x=np.asarray(x,float);m=np.max(x);return m+np.log(np.sum(np.exp(x-m)))
def labf(b,se,W):
    V=se**2;Z2=(b/se)**2;r=W/(V+W);return 0.5*np.log(1-r)+0.5*r*Z2
def coloc(d1,d2,W2,W1=0.15**2,p1=1e-4,p2=1e-4,p12=1e-5):
    d1=d1.drop_duplicates("SNP");d2=d2.drop_duplicates("SNP")
    m=d1.merge(d2,on="SNP",suffixes=("_1","_2")).reset_index(drop=True)
    m=m[(m.se_1>0)&(m.se_2>0)&np.isfinite(m.beta_1)&np.isfinite(m.beta_2)&np.isfinite(m.se_1)&np.isfinite(m.se_2)].reset_index(drop=True)
    if len(m)<5: return None
    l1=labf(m.beta_1.values,m.se_1.values,W1);l2=labf(m.beta_2.values,m.se_2.values,W2)
    L1,L2,L12=lse(l1),lse(l2),lse(l1+l2)
    lpp=[0.0,math.log(p1)+L1,math.log(p2)+L2,math.log(p1)+math.log(p2)+(L1+L2+np.log1p(-math.exp(L12-(L1+L2)))),math.log(p12)+L12]
    Z=lse(lpp);pp=[math.exp(x-Z) for x in lpp]
    return dict(n=len(m),PP3=pp[3],PP4=pp[4])
# eQTL regions
e2=pd.read_csv("eqtl_region_converted.csv"); rock2=e2[e2.GeneSymbol=="ROCK2"][["SNP","beta","se"]]
erho=pd.read_csv("rho_eqtl_region_conv.csv")
def eqrho(g): return erho[erho.GeneSymbol==g][["SNP","beta","se"]]
def fg(f):
    d=pd.read_csv(f,sep="\t").rename(columns={"#chrom":"c"})
    d=d[["rsids","beta","sebeta"]].dropna(); d.columns=["SNP","beta","se"]; return d
def nf(f):
    d=pd.read_csv(f,sep="\t")[["hm_rsid","hm_beta","standard_error"]].dropna(); d.columns=["SNP","beta","se"]; return d
outs={"Cirrhosis":fg("finngen_CIRRHOSIS_BROAD_region.tsv"),
      "NAFLD_meta":nf("nafld_meta_region.tsv"),
      "FibrChirK74":fg("finngen_K11_FIBROCHIRLIV_region.tsv"),
      "NASH":fg("finngen_NASH_region.tsv")}
CC=0.2**2
print("=== ROCK2 (blood eQTL) x fibrosis-axis coloc ===")
for nm,d2 in outs.items():
    r=coloc(rock2,d2,CC); print(f"  ROCK2 x {nm:12s}: PP3={r['PP3']:.3f} PP4={r['PP4']:.3f} (n={r['n']})" if r else f"  ROCK2 x {nm}: NA")
print("\n=== ROCK1 / RHOA (blood eQTL) x Cirrhosis ===")
for g in ["ROCK1","RHOA"]:
    r=coloc(eqrho(g),outs["Cirrhosis"],CC); print(f"  {g} x Cirrhosis: PP3={r['PP3']:.3f} PP4={r['PP4']:.3f} (n={r['n']})" if r else f"  {g}: NA")
# ROCK1/RHOA cirrhosis MR (lead Wald)
rho=pd.read_csv("rho_converted.csv"); exp_rho={r.SNP:dict(EA=r.AssessedAllele,OA=r.OtherAllele,beta=r.beta,se=r.se) for r in rho.itertuples()}
cir=pd.read_csv("finngen_CIRRHOSIS_BROAD_region.tsv",sep="\t").rename(columns={"#chrom":"c"})
codd={r.rsids:(r.alt,r.ref,r.beta,r.sebeta,r.pval) for r in cir.itertuples() if isinstance(r.rsids,str)}
def wald(snp):
    e=exp_rho.get(snp); o=codd.get(snp)
    if not e or not o: return None
    EA,OA,bo,so,po=o
    if EA==e["OA"] and OA==e["EA"]: bo=-bo
    elif not(EA==e["EA"] and OA==e["OA"]): return None
    b=bo/e["beta"]; se=abs(so/e["beta"]); return b,se,z2p(b/se)
print("\n=== ROCK1/RHOA x Cirrhosis MR (lead Wald) ===")
for g,snp in [("ROCK1","rs142716063"),("RHOA","rs115431681")]:
    w=wald(snp); print(f"  {g} ({snp}):", f"{w[0]:+.3f} [{w[0]-1.96*w[1]:+.3f},{w[0]+1.96*w[1]:+.3f}] p={w[2]:.2g}" if w else "NA")
