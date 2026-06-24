import pandas as pd, numpy as np
from statistics import NormalDist
_nd=NormalDist()
def z2p(z): return 2*(1-_nd.cdf(abs(z)))
import math
def chi2sf(Q,df):
    if df==1: return 2*(1-_nd.cdf(math.sqrt(max(Q,0))))
    if df==2: return math.exp(-Q/2)
    # general (regularized upper incomplete gamma, series) for small df
    from math import exp,lgamma
    a=df/2.0; x=Q/2.0
    if x<=0: return 1.0
    if x< a+1:
        term=1.0/a; ssum=term; n=a
        for _ in range(200):
            n+=1; term*=x/n; ssum+=term
            if term<ssum*1e-12: break
        return 1- ssum*exp(-x+a*math.log(x)-lgamma(a))
    else:
        b=x+1-a; c=1e300; d=1/b; h=d
        for i in range(1,200):
            an=-i*(i-a); b+=2; d=an*d+b; d=1/d if d else 1e-300; c=b+an/c; del_=d*c; h*=del_
            if abs(del_-1)<1e-12: break
        return exp(-x+a*math.log(x)-lgamma(a))*h

conv=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")
exp_blood={r.SNP:dict(EA=r.AssessedAllele,OA=r.OtherAllele,beta=r.beta,se=r.se,eaf=r.eaf) for r in conv.itertuples()}
exp_musc={"rs4335920":dict(EA="A",OA="T",beta=-0.187116,se=0.023170,eaf=0.333),
          "rs12996712":dict(EA="T",OA="G",beta=0.133808,se=0.021600,eaf=0.483)}
lead="rs12468344"; r201=["rs12468344","rs13034888","rs7574150"]
r21=pd.read_csv("blood_clumped_independent.txt",sep="\t")["SNP"].tolist()
muscle=["rs4335920","rs12996712"]

rows=[]
api=pd.read_csv("api_outcomes_long.tsv",sep="\t")
for r in api.itertuples():
    if pd.notna(r.beta) and pd.notna(r.se): rows.append([r.outcome,r.SNP,r.EA,r.OA,float(r.beta),float(r.se),r.eaf])
lf=pd.read_csv("liverfat_harmonization_preview.txt",sep="\t")
for r in lf.itertuples(): rows.append(["HepaticFat",r.SNP,r.out_EA,r.out_OA,float(r.out_beta),float(r.out_se),float(r.out_eaf)])
bmi=pd.read_csv("bmi_rock2_pulit.tsv",sep="\t"); bmi["rs"]=bmi["SNP"].str.replace(r":.*","",regex=True)
for r in bmi.itertuples(): rows.append(["BMI",r.rs,r.Allele1.upper(),r.Allele2.upper(),float(r.Effect),float(r.StdErr),float(r.Freq1)])
out=pd.DataFrame(rows,columns=["outcome","SNP","EA","OA","beta","se","eaf"])

def get_exp(s): return exp_musc[s] if s in exp_musc else exp_blood.get(s)
def harm(snp,o):
    e=get_exp(snp); row=out[(out.outcome==o)&(out.SNP==snp)]
    if e is None or len(row)==0: return None
    row=row.iloc[0]; bo,eaf_o=row.beta,row.eaf
    if row.EA==e["EA"] and row.OA==e["OA"]: pass
    elif row.EA==e["OA"] and row.OA==e["EA"]: bo=-bo; eaf_o=1-eaf_o
    else: return None
    if not (np.isfinite(bo) and np.isfinite(row.se) and row.se>0): return None
    return dict(bx=e["beta"],sx=e["se"],by=bo,sy=row.se,eaf_x=e["eaf"],eaf_o=eaf_o,bx_se=e["se"])
def wald(snp,o):
    h=harm(snp,o); 
    if not h or h["bx"]==0: return None
    b=h["by"]/h["bx"]; se=abs(h["sy"]/h["bx"]); return b,se,z2p(b/se)
def ivw(snps,o):
    d=[harm(s,o) for s in snps]; d=[x for x in d if x]
    if len(d)<2: return None
    bx=np.array([x["bx"] for x in d]);by=np.array([x["by"] for x in d]);sy=np.array([x["sy"] for x in d])
    w=bx**2/sy**2; b=np.sum(bx*by/sy**2)/np.sum(w); se=np.sqrt(1/np.sum(w))
    Q=float(np.sum((by-b*bx)**2/sy**2)); df=len(d)-1; Qp=chi2sf(Q,df)
    return b,se,z2p(b/se),Q,Qp,len(d)
def ivw_corr(snps,o,R,Rsnps):
    d={}
    for s in snps:
        h=harm(s,o)
        if h: d[s]=h
    order=[s for s in Rsnps if s in d]
    if len(order)<2: return None
    bx=np.array([d[s]["bx"] for s in order]);by=np.array([d[s]["by"] for s in order]);sy=np.array([d[s]["sy"] for s in order])
    idx=[Rsnps.index(s) for s in order]; Rs=R[np.ix_(idx,idx)].copy()
    Rs=(Rs+Rs.T)/2; np.fill_diagonal(Rs,1.0)
    Sig=np.outer(sy,sy)*Rs + np.eye(len(order))*1e-12
    try: Si=np.linalg.solve(Sig,np.eye(len(order)))
    except: Si=np.linalg.pinv(Sig,rcond=1e-10)
    A=bx@Si@bx; b=(bx@Si@by)/A; se=np.sqrt(1/A)
    return b,se,z2p(b/se),len(order)
def steiger(snp,o):
    h=harm(snp,o)
    if not h: return None
    rx2=2*h["eaf_x"]*(1-h["eaf_x"])*h["bx"]**2
    ro2=2*h["eaf_o"]*(1-h["eaf_o"])*h["by"]**2
    return rx2,ro2,("exp->out" if rx2>ro2 else "AMBIG/rev")

Rdf=pd.read_csv("ld_matrix_13.csv",index_col=0); R=Rdf.values; Rsnps=list(Rdf.index)
outcomes=["HepaticFat","FI_EUR","FG_EUR","HbA1c_EUR","2hGlu_EUR","T2D_Xue2018","BMI"]
def fmt(b,se,p):
    return f"{b:+.3f} [{b-1.96*se:+.3f},{b+1.96*se:+.3f}] p={p:.2g}"
res=[]
print("="*92); print("ROCK2 MR — per-SD-ROCK2-expression effect (T2D = log-OR)"); print("="*92)
for o in outcomes:
    print(f"\n### {o}")
    w=wald(lead,o);  print(f"  Blood lead Wald        : {fmt(*w)}") if w else None
    if w: res.append([o,"blood_lead_Wald",*w])
    iv=ivw(muscle,o)
    if iv: print(f"  Muscle 2-SNP IVW       : {fmt(iv[0],iv[1],iv[2])} | Q={iv[3]:.2f} Qp={iv[4]:.2f}"); res.append([o,"muscle_IVW",iv[0],iv[1],iv[2]])
    for s in muscle:
        ws=wald(s,o)
        if ws: print(f"    {s} Wald        : {fmt(*ws)}")
    iv2=ivw(r201,o)
    if iv2: print(f"  Blood r2<.01 IVW (n={iv2[5]})  : {fmt(iv2[0],iv2[1],iv2[2])} | Q={iv2[3]:.2f} Qp={iv2[4]:.2f}"); res.append([o,"blood_r2.01_IVW",iv2[0],iv2[1],iv2[2]])
    ic=ivw_corr(r21,o,R,Rsnps)
    if ic: print(f"  Blood r2<.1 LDcorr IVW(n={ic[3]}): {fmt(ic[0],ic[1],ic[2])}"); res.append([o,"blood_r2.1_LDcorr_IVW",ic[0],ic[1],ic[2]])
    # tissue concordance
    if w and iv: print(f"  >> tissue concordance  : blood {w[0]:+.3f} vs muscle {iv[0]:+.3f}  ({'CONCORDANT sign' if np.sign(w[0])==np.sign(iv[0]) else 'discordant'})")
    st=steiger(lead,o)
    if st: print(f"  >> Steiger (lead)      : r2_exp={st[0]:.3f} r2_out={st[1]:.2e} -> {st[2]}")
pd.DataFrame(res,columns=["outcome","method","beta","se","p"]).to_csv("rock2_mr_results.csv",index=False)
print("\nsaved rock2_mr_results.csv")
