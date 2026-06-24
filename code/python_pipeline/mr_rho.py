import urllib.request,json,pandas as pd,numpy as np,math
from statistics import NormalDist
nd=NormalDist()
def z2p(z):return 2*(1-nd.cdf(abs(z)))
def zfp(p):
    if p<=0:p=1e-300
    a=1-p/2
    return math.sqrt(max(2*(-math.log(p))-math.log(2*math.pi*(-math.log(p))),1e-6)) if a>=1 else nd.inv_cdf(a)
m=pd.read_csv("rho_converted.csv")
exp={r.SNP:dict(EA=r.AssessedAllele,OA=r.OtherAllele,beta=r.beta,se=r.se,gene=r.GeneSymbol) for r in m.itertuples()}
RHOA=["rs115431681"]; ROCK1_lead="rs142716063"; ROCK1_2=["rs142716063","rs72881399"]
def fetch(rs,g):
    try:
        with urllib.request.urlopen(f"https://www.ebi.ac.uk/gwas/summary-statistics/api/associations/{rs}?study_accession={g}",timeout=20) as r:
            d=json.load(r)
        return d if "beta" in d else None
    except:return None
accs={"FI":"GCST90002238","FG":"GCST90002232","HbA1c":"GCST90002244","2hGlu":"GCST90002227","T2D":"GCST006867"}
# local
bmi=pd.read_csv("rho_bmi.tsv",sep="\t");bmi["rs"]=bmi["SNP"].str.replace(r":.*","",regex=True)
hf=pd.read_csv("rho_hf.tsv",sep="\t")
def out_get(rs,oc):
    if oc=="BMI":
        r=bmi[bmi.rs==rs]
        if len(r):return dict(EA=r.Allele1.iloc[0].upper(),OA=r.Allele2.iloc[0].upper(),beta=float(r.Effect.iloc[0]),se=float(r.StdErr.iloc[0]))
    elif oc=="HepFat":
        r=hf[hf.hm_rsid==rs]
        if len(r):return dict(EA=r.hm_effect_allele.iloc[0],OA=r.hm_other_allele.iloc[0],beta=float(r.hm_beta.iloc[0]),se=float(r.standard_error.iloc[0]))
    else:
        d=fetch(rs,accs[oc])
        if d:
            b=d['beta'];p=d['p_value'];se=abs(b)/zfp(p) if (b!=0 and 0<p<1) else float('nan')
            return dict(EA=d['effect_allele'],OA=d['other_allele'],beta=b,se=se)
    return None
def harm(rs,oc):
    e=exp[rs];o=out_get(rs,oc)
    if not o or not np.isfinite(o['se']) or o['se']<=0:return None
    bo=o['beta']
    if o['EA']==e['EA'] and o['OA']==e['OA']:pass
    elif o['EA']==e['OA'] and o['OA']==e['EA']:bo=-bo
    else:return None
    return dict(bx=e['beta'],sx=e['se'],by=bo,sy=o['se'])
def wald(rs,oc):
    h=harm(rs,oc)
    if not h:return None
    b=h['by']/h['bx'];se=abs(h['sy']/h['bx']);return b,se,z2p(b/se)
def ivw(snps,oc):
    d=[harm(s,oc) for s in snps];d=[x for x in d if x]
    if len(d)<2:return None
    bx=np.array([x['bx'] for x in d]);by=np.array([x['by'] for x in d]);sy=np.array([x['sy'] for x in d])
    w=bx**2/sy**2;b=np.sum(bx*by/sy**2)/np.sum(w);se=np.sqrt(1/np.sum(w));return b,se,z2p(b/se),len(d)
ocs=["HepFat","FI","FG","HbA1c","2hGlu","T2D","BMI"]
def f(r):return f"{r[0]:+.3f}[{r[0]-1.96*r[1]:+.3f},{r[0]+1.96*r[1]:+.3f}]p={r[2]:.2g}"
print("=== RHOA (1-SNP lead Wald, rs115431681 F=27) ===")
for oc in ocs:
    w=wald("rs115431681",oc); print(f"  {oc:8s}: {f(w) if w else 'NA'}")
print("\n=== ROCK1 lead Wald (rs142716063 F=86) ===")
for oc in ocs:
    w=wald(ROCK1_lead,oc); print(f"  {oc:8s}: {f(w) if w else 'NA'}")
print("\n=== ROCK1 r2<.01 IVW (2 SNP) ===")
for oc in ocs:
    iv=ivw(ROCK1_2,oc); print(f"  {oc:8s}: {f(iv) if iv else 'NA'}")
