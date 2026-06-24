import pandas as pd, numpy as np, math
from statistics import NormalDist
nd=NormalDist()
def z2p(z): return 2*(1-nd.cdf(abs(z)))
# ---- exposures ----
conv=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")
exp_blood={r.SNP:dict(EA=r.AssessedAllele,OA=r.OtherAllele,beta=r.beta,se=r.se,eaf=r.eaf) for r in conv.itertuples()}
exp_musc={"rs4335920":dict(EA="A",OA="T",beta=-0.187116,se=0.023170,eaf=0.333),
          "rs12996712":dict(EA="T",OA="G",beta=0.133808,se=0.021600,eaf=0.483)}
rho=pd.read_csv("rho_converted.csv")
exp_rho={r.SNP:dict(EA=r.AssessedAllele,OA=r.OtherAllele,beta=r.beta,se=r.se) for r in rho.itertuples()}
lead="rs12468344"; r201=["rs12468344","rs13034888","rs7574150"]; muscle=["rs4335920","rs12996712"]
r21=pd.read_csv("blood_clumped_independent.txt",sep="\t")["SNP"].tolist()
snps17=[l.strip() for l in open("instrument_snps_for_outcome.txt")]
# ---- outcome loaders -> dict SNP-> (EA,OA,beta,se,p) ----
def load_finngen(f):
    d=pd.read_csv(f,sep="\t")
    d=d.rename(columns={"#chrom":"chrom"})
    out={}
    for r in d.itertuples():
        rs=r.rsids
        if isinstance(rs,str) and rs.startswith("rs") and np.isfinite(r.beta) and r.sebeta>0:
            out[rs]=(r.alt,r.ref,float(r.beta),float(r.sebeta),float(r.pval))
    return out
def load_nafld(f):
    d=pd.read_csv(f,sep="\t")
    out={}
    for r in d.itertuples():
        rs=r.hm_rsid
        if isinstance(rs,str) and rs.startswith("rs") and pd.notna(r.hm_beta) and pd.notna(r.standard_error) and r.standard_error>0:
            out[rs]=(r.hm_effect_allele,r.hm_other_allele,float(r.hm_beta),float(r.standard_error),float(r.p_value))
    return out
OUT={
 "Cirrhosis(5545)":load_finngen("finngen_CIRRHOSIS_BROAD_region.tsv"),
 "NAFLD_meta(8434)":load_nafld("nafld_meta_region.tsv"),
 "FibrChirK74(2653)":load_finngen("finngen_K11_FIBROCHIRLIV_region.tsv"),
 "NASH(254)":load_finngen("finngen_NASH_region.tsv"),
}
def harm(expd,snp,od):
    e=expd.get(snp); o=od.get(snp)
    if not e or not o: return None
    EA,OA,bo,so,po=o
    if EA==e["EA"] and OA==e["OA"]: pass
    elif EA==e["OA"] and OA==e["EA"]: bo=-bo
    else: return None
    if not(np.isfinite(bo) and so>0): return None
    return dict(bx=e["beta"],by=bo,sy=so,po=po)
def wald(expd,snp,od):
    h=harm(expd,snp,od)
    if not h or h["bx"]==0: return None
    b=h["by"]/h["bx"]; se=abs(h["sy"]/h["bx"]); return b,se,z2p(b/se)
def ivw(expd,snps,od):
    H=[harm(expd,s,od) for s in snps]; H=[x for x in H if x]
    if len(H)<2: return None
    bx=np.array([x["bx"] for x in H]);by=np.array([x["by"] for x in H]);sy=np.array([x["sy"] for x in H])
    w=bx**2/sy**2;b=np.sum(bx*by/sy**2)/np.sum(w);se=np.sqrt(1/np.sum(w));return b,se,z2p(b/se),len(H)
def expall(snp):
    return exp_musc.get(snp) or exp_blood.get(snp)
def cov(snps,od): return sum(1 for s in snps if s in od)
print("="*70); print("ROCK2 — fibrosis axis"); print("="*70)
for name,od in OUT.items():
    print(f"\n### {name}  (tier coverage: lead {cov([lead],od)}/1, r01 {cov(r201,od)}/3, r2.1 {cov(r21,od)}/13, muscle {cov(muscle,od)}/2)")
    # single-SNP preview
    for s,lab in [(lead,"blood lead"),("rs4335920","muscle"),("rs12996712","muscle")]:
        h=harm(expall(s),s,od) if False else None
    w=wald(exp_blood,lead,od)
    iv=ivw({**{k:exp_musc[k] for k in exp_musc}},muscle,od)
    def f(r): return f"{r[0]:+.3f} [{r[0]-1.96*r[1]:+.3f},{r[0]+1.96*r[1]:+.3f}] p={r[2]:.2g}" if r else "NA"
    print("  Blood lead Wald   :", f(w))
    print("  Muscle 2-SNP IVW  :", f(iv[:3]) if iv else "NA")
    iv2=ivw(exp_blood,r201,od)
    print("  Blood r2<.01 IVW  :", f(iv2[:3]) if iv2 else "NA")
