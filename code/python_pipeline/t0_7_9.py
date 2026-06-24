import pandas as pd, numpy as np, math, urllib.request, json
from statistics import NormalDist
nd=NormalDist()
def p_from_z(z): return 2*(1-nd.cdf(abs(z)))
print("="*64); print("TASK 0a — Z->beta/se internal round-trip (eQTLGen ROCK2)"); print("="*64)
e=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")
err=(e.beta/e.se - e.Zscore).abs().max()
print(f"max |beta/se - Zscore| = {err:.2e}  -> {'PASS' if err<1e-6 else 'FAIL'}")
bre=e.Zscore/np.sqrt(2*e.eaf*(1-e.eaf)*(e.n+e.Zscore**2))
print(f"max |beta - recompute(Z,eaf,n)| = {(e.beta-bre).abs().max():.2e}  (should ~0)")
print(f"EAF in (0,1): {bool(((e.eaf>0)&(e.eaf<1)).all())} | alleles non-missing: {bool((e.AssessedAllele.notna()&e.OtherAllele.notna()).all())}")

print("\n"+"="*64); print("TASK 0a-2 — direction concordance blood vs GTEx muscle (eQTL Catalogue)"); print("="*64)
def ec_dataset(tissue):
    u="https://www.ebi.ac.uk/eqtl/api/v2/datasets/?study_label=GTEx&quant_method=ge&size=1000"
    d=json.load(urllib.request.urlopen(u,timeout=30))
    for r in d:
        if tissue.lower() in (str(r.get('sample_group',''))+str(r.get('tissue_label',''))).lower(): return r['dataset_id']
    return None
def ec_gene(ds,ensg,maxp=20):
    rows=[];start=0
    for _ in range(maxp):
        u=f"https://www.ebi.ac.uk/eqtl/api/v2/datasets/{ds}/associations?molecular_trait_id={ensg}&size=1000&start={start}"
        try: pg=json.load(urllib.request.urlopen(u,timeout=30))
        except Exception as ex: break
        if not pg: break
        rows+=pg
        if len(pg)<1000: break
        start+=1000
    return pd.DataFrame(rows)
ds=ec_dataset("muscle"); print("GTEx muscle dataset_id:",ds)
g=ec_gene(ds,"ENSG00000134318") if ds else pd.DataFrame()
if len(g):
    g=g[['rsid','ref','alt','beta']].dropna()
    m=e[['SNP','AssessedAllele','OtherAllele','beta']].merge(g,left_on='SNP',right_on='rsid')
    def al(r): return r['beta_y'] if r['alt']==r['AssessedAllele'] else (-r['beta_y'] if r['ref']==r['AssessedAllele'] else np.nan)
    m['slope']=m.apply(al,axis=1); m=m.dropna(subset=['slope'])
    print(f"common SNP: {len(m)} | sign concordance: {(np.sign(m.beta_x)==np.sign(m.slope)).mean():.3f}  (magnitude NOT a pass/fail)")
else: print("muscle fetch empty — skip")

print("\n"+"="*64); print("TASK 9 — power / minimum detectable effect"); print("="*64)
inst=pd.read_csv("/sessions/determined-eager-davinci/mnt/후속연구/O2_scripts/data/instruments_frozen.csv")
r2=inst[(inst.gene=="ROCK2")&(inst.tier.isin(["ROCK2_blood_lead","ROCK2_blood_r2.01"]))]
R2x=float((2*r2.eaf*(1-r2.eaf)*r2.beta**2).sum())
za,zb=nd.inv_cdf(0.975),nd.inv_cdf(0.8)
def mde_cont(N): return (za+zb)/math.sqrt(N*R2x)
def mde_bin(nc,nk): Neff=4/(1/nc+1/nk); return (za+zb)/math.sqrt(Neff*R2x)
print(f"ROCK2 instrument R2 (r2<.01 set) = {R2x:.4f}")
rows=[("BMI cont N=806834",mde_cont(806834),"SD",""),
      ("HepaticFat cont N=32974",mde_cont(32974),"SD",""),
      ("NAFLD 8434/770180",mde_bin(8434,770180),"logOR",""),
      ("Cirrhosis 5545/494803",mde_bin(5545,494803),"logOR",""),
      ("T2D DIAMANTE 80154/853816",mde_bin(80154,853816),"logOR","")]
for n,v,sc,_ in rows:
    extra=f"  OR={math.exp(v):.2f}" if sc=="logOR" else ""
    print(f"  {n:30s} MDE={v:.3f} {sc}{extra}")
