import pysam, numpy as np, pandas as pd
conv=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")
bypos={int(r.SNPPos):r for r in conv.itertuples()}; positions=set(bypos)
eur=set(l.strip() for l in open("eur.ids"))
url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
tb=pysam.TabixFile(url); samples=list(tb.header)[-1].split("\t")[9:]
eur_idx=[i for i,s in enumerate(samples) if s in eur]
def dose(g): a=g[:3]; return (a[0]=='1')+(a[2]=='1') if (a[0] in '01' and a[2] in '01') else np.nan
geno={};meta={}
for line in tb.fetch("2",11100000,11500000):
    f=line.split("\t");pos=int(f[1])
    if pos not in positions:continue
    ref,alt=f[3],f[4]
    if ","in alt:continue
    row=bypos[pos]
    if {ref,alt}!={row.AssessedAllele,row.OtherAllele}:continue
    geno[row.SNP]=np.array([dose(f[9:][i]) for i in eur_idx],float)
    meta[row.SNP]=(pos,row.beta,row.se,row.F,row.Pvalue)
matched=list(geno);M=np.vstack([geno[s] for s in matched]);idx={s:i for i,s in enumerate(matched)}
def r2(a,b):
    x=M[idx[a]];y=M[idx[b]];m=~(np.isnan(x)|np.isnan(y))
    if m.sum()<50:return 0.0
    c=np.corrcoef(x[m],y[m])[0,1];return 0.0 if np.isnan(c) else c*c
def clump(thr):
    order=sorted(matched,key=lambda s:meta[s][4]);kept=[];rem=set()
    for s in order:
        if s in rem:continue
        kept.append(s)
        for t in order:
            if t!=s and t not in rem and r2(s,t)>=thr:rem.add(t)
    return kept
for thr in (0.1,0.01,0.001):
    print("r2<%-5s -> %d independent signals"%(thr,len(clump(thr))))
