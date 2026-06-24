import pysam, numpy as np, pandas as pd
conv=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")
bypos={int(r.SNPPos):r for r in conv.itertuples()};positions=set(bypos)
eur=set(l.strip() for l in open("eur.ids"))
url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
tb=pysam.TabixFile(url);samples=list(tb.header)[-1].split("\t")[9:]
eur_idx=[i for i,s in enumerate(samples) if s in eur]
def dose(g):a=g[:3];return (a[0]=='1')+(a[2]=='1') if (a[0] in '01' and a[2] in '01') else np.nan
geno={};meta={}
musc_pos={11415777:"rs4335920",11261678:"rs12996712"}
for line in tb.fetch("2",11100000,11500000):
    f=line.split("\t");pos=int(f[1]);ref,alt=f[3],f[4]
    if ","in alt:continue
    if pos in positions:
        row=bypos[pos]
        if {ref,alt}=={row.AssessedAllele,row.OtherAllele}:
            geno[row.SNP]=np.array([dose(f[9:][i]) for i in eur_idx],float)
            meta[row.SNP]=(pos,row.beta,row.se,row.F,row.Pvalue)
    if pos in musc_pos and musc_pos[pos] not in geno:
        geno[musc_pos[pos]]=np.array([dose(f[9:][i]) for i in eur_idx],float)
matched=[s for s in geno if s in meta];M=np.vstack([geno[s] for s in geno]);idx={s:i for i,s in enumerate(geno)}
def r2(a,b):
    x=M[idx[a]];y=M[idx[b]];m=~(np.isnan(x)|np.isnan(y))
    if m.sum()<50 or a not in idx or b not in idx:return np.nan
    c=np.corrcoef(x[m],y[m])[0,1];return np.nan if np.isnan(c) else c*c
def clump(thr):
    order=sorted(matched,key=lambda s:meta[s][4]);kept=[];rem=set()
    for s in order:
        if s in rem:continue
        kept.append(s)
        for t in order:
            if t!=s and t not in rem and r2(s,t)>=thr:rem.add(t)
    return kept
s10=clump(0.1); s01=clump(0.01)
print("=== Blood supporting set r2<0.01 (3 near-independent) ===")
for s in sorted(s01,key=lambda x:meta[x][4]):
    pos,b,se,F,p=meta[s];print(f"  {s:12s} pos{pos} beta={b:+.4f} F={F:.1f} p={p:.1e}")
print("\n=== Muscle mutual r2 (EUR) ===")
print("  rs4335920 - rs12996712 r2 = %.4f"%r2("rs4335920","rs12996712"))
print("  in 1000G EUR: rs4335920=%s rs12996712=%s"%("rs4335920" in geno,"rs12996712" in geno))
# outcome-extraction union list = 13 blood (r2<0.1) + 2 muscle
union=sorted(set(s10)|{"rs4335920","rs12996712"})
pd.Series(union).to_csv("instrument_snps_for_outcome.txt",index=False,header=False)
print("\n=== FROZEN outcome-extraction SNP list: %d unique ==="%len(union))
print(", ".join(union))
