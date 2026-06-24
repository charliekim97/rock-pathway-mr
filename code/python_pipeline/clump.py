import pysam, numpy as np, pandas as pd

conv = pd.read_csv("rock2_eqtlgen_converted.txt", sep="\t")
# map by hg19 pos
bypos = {int(r.SNPPos): r for r in conv.itertuples()}
positions = set(bypos)
eur = set(l.strip() for l in open("eur.ids"))

url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
tb = pysam.TabixFile(url)
samples = list(tb.header)[-1].split("\t")[9:]
eur_idx = [i for i,s in enumerate(samples) if s in eur]
print("EUR cols matched:", len(eur_idx))

geno={}; meta={}
def dose(gt):
    a=gt[:3]; 
    return (a[0]=='1')+(a[2]=='1') if (a[0] in '01' and a[2] in '01') else np.nan
for line in tb.fetch("2",11100000,11500000):
    f=line.split("\t")
    pos=int(f[1])
    if pos not in positions: continue
    ref,alt=f[3],f[4]
    if "," in alt: continue  # skip multiallelic
    row=bypos[pos]
    aset={row.AssessedAllele, row.OtherAllele}
    if {ref,alt}!=aset: continue   # allele match
    fmt=f[9:]
    d=np.array([dose(fmt[i]) for i in eur_idx],dtype=float)
    geno[row.SNP]=d
    meta[row.SNP]=(pos,ref,alt,row.beta,row.se,row.F,row.Pvalue)

matched=list(geno)
print("SNPs matched to 1000G EUR:", len(matched), "/ 621")
unmatched=[r.SNP for r in conv.itertuples() if r.SNP not in geno]
print("unmatched (kept aside):", len(unmatched))

# genotype matrix
M=np.vstack([geno[s] for s in matched])
# greedy clump by p ascending, r2<0.1, whole-region window
order=sorted(matched, key=lambda s: meta[s][6])
idx={s:i for i,s in enumerate(matched)}
def r2(a,b):
    x=M[idx[a]];y=M[idx[b]]
    m=~(np.isnan(x)|np.isnan(y))
    if m.sum()<50: return 0.0
    c=np.corrcoef(x[m],y[m])[0,1]
    return 0.0 if np.isnan(c) else c*c
kept=[]; removed=set()
for s in order:
    if s in removed: continue
    kept.append(s); 
    for t in order:
        if t!=s and t not in removed and r2(s,t)>=0.1:
            removed.add(t)
print("\n=== INDEPENDENT SIGNALS (r2<0.1): %d ===" % len(kept))
rows=[]
for s in kept:
    pos,ref,alt,b,se,F,p=meta[s]
    rows.append((s,pos,ref,alt,round(b,4),round(se,5),round(F,1),p))
res=pd.DataFrame(rows,columns=["SNP","pos_hg19","ref","alt","beta","se","F","pval"]).sort_values("pval")
print(res.to_string(index=False))
res.to_csv("blood_clumped_independent.txt",sep="\t",index=False)
print("\nlead SNP:", res.iloc[0].SNP, "F=",res.iloc[0].F)
