import urllib.request,json,pandas as pd,numpy as np
from statistics import NormalDist
nd=NormalDist()
def z2p(z): return 2*(1-nd.cdf(abs(z)))
import math
def zfp(p):
    if p<=0: p=1e-300
    arg=1-p/2
    if arg>=1.0:
        L=-math.log(p); return math.sqrt(max(2*L-math.log(2*math.pi*L),1e-6))
    return nd.inv_cdf(arg)
def se_from(b,p):
    if b==0 or p>=1: return float('nan')
    return abs(b)/zfp(p)
def fetch(rs,g):
    try:
        with urllib.request.urlopen(f"https://www.ebi.ac.uk/gwas/summary-statistics/api/associations/{rs}?study_accession={g}",timeout=20) as r:
            d=json.load(r)
        return d if "beta" in d else None
    except: return None
# exposures (eQTLGen lead) : snp, EA, OA, beta, se
ctrl=pd.read_csv("ctrl_converted.csv")
def lead(gene):
    s=ctrl[ctrl.GeneSymbol==gene].sort_values("Pvalue").iloc[0]
    return dict(SNP=s.SNP,EA=s.AssessedAllele,OA=s.OtherAllele,beta=s.beta,se=s.se,F=s.F)
s1=pd.read_csv("sort1_eqtl.txt",sep="\t")
# SORT1 lead recompute eaf/beta
import numpy as np
af=pd.read_csv("sort1_af.txt",sep="\t")[["SNP","AlleleA","AlleleB","AlleleB_all"]]
m=s1.merge(af,on="SNP");m["eaf"]=[r.AlleleB_all if r.AssessedAllele==r.AlleleB else 1-r.AlleleB_all for r in m.itertuples()]
Z=m.Zscore;p=m.eaf;N=m.NrSamples;den=np.sqrt(2*p*(1-p)*(N+Z**2));m["beta"]=Z/den;m["se"]=1/den
sr=m.sort_values("Pvalue").iloc[0]; sort1=dict(SNP=sr.SNP,EA=sr.AssessedAllele,OA=sr.OtherAllele,beta=sr.beta,se=sr.se)

def wald(exp,out_d):
    bo=out_d['beta']
    if out_d['effect_allele']==exp['EA']: pass
    elif out_d['effect_allele']==exp['OA']: bo=-bo
    else: return None
    so=se_from(out_d['beta'],out_d['p_value'])
    b=bo/exp['beta']; se=abs(so/exp['beta']); return b,se,z2p(b/se),out_d['p_value']

print("=== POSITIVE CONTROLS (eQTLGen blood cis-eQTL -> outcome, lead-SNP Wald) ===\n")
# KCNJ11 -> T2D
k=lead("KCNJ11"); d=fetch(k['SNP'],"GCST006867")
r=wald(k,d); print(f"KCNJ11 -> T2D  | eQTL {k['SNP']} F={k['F']:.0f} | Wald logOR={r[0]:+.3f} [{r[0]-1.96*r[1]:+.3f},{r[0]+1.96*r[1]:+.3f}] p={r[2]:.2g} (outcome SNP p={r[3]:.1e})")
# GIPR -> T2D (find SNP present in Xue)
for cand in "rs12972094,rs4994276,rs12972282,rs4380143,rs12974310,rs2238689,rs2302382,rs35845603".split(","):
    row=ctrl[(ctrl.GeneSymbol=="GIPR")&(ctrl.SNP==cand)]
    if not len(row): continue
    d=fetch(cand,"GCST006867")
    if d:
        g=dict(EA=row.iloc[0].AssessedAllele,OA=row.iloc[0].OtherAllele,beta=row.iloc[0].beta,se=row.iloc[0].se,F=row.iloc[0].F)
        r=wald(g,d); print(f"GIPR   -> T2D  | eQTL {cand} F={g['F']:.0f} | Wald logOR={r[0]:+.3f} [{r[0]-1.96*r[1]:+.3f},{r[0]+1.96*r[1]:+.3f}] p={r[2]:.2g} (outcome SNP p={r[3]:.1e})")
        break
# SORT1 -> LDL (Willer2013 GCST002222)
d=fetch(sort1['SNP'],"GCST002222")
r=wald(sort1,d); print(f"SORT1  -> LDL  | eQTL {sort1['SNP']} | Wald betaSD={r[0]:+.3f} [{r[0]-1.96*r[1]:+.3f},{r[0]+1.96*r[1]:+.3f}] p={r[2]:.2g} (outcome SNP p={r[3]:.1e})")

print()
gl=ctrl[ctrl.GeneSymbol=="GIPR"].sort_values("Pvalue").iloc[0]
gexp=dict(EA=gl.AssessedAllele,OA=gl.OtherAllele,beta=gl.beta,se=gl.se,F=gl.F)
for trait,acc in [("2hGlu","GCST90002227"),("FG","GCST90002232"),("T2D-Xue","GCST006867")]:
    d=fetch(gl.SNP,acc)
    if d:
        r=wald(gexp,d); print(f"GIPR   -> {trait:8s}| eQTL {gl.SNP} F={gexp['F']:.0f} | Wald={r[0]:+.3f} [{r[0]-1.96*r[1]:+.3f},{r[0]+1.96*r[1]:+.3f}] p={r[2]:.2g} (SNP p={r[3]:.1e})")
    else: print(f"GIPR   -> {trait:8s}| {gl.SNP} not in study")
