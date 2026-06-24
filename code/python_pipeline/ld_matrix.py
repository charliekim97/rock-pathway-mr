import pysam, numpy as np, pandas as pd
conv=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")
b13=pd.read_csv("blood_clumped_independent.txt",sep="\t")["SNP"].tolist()
sub=conv[conv.SNP.isin(b13)][["SNP","SNPPos","AssessedAllele","OtherAllele"]].copy()
eur=set(l.strip() for l in open("eur.ids"))
url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
tb=pysam.TabixFile(url);samples=list(tb.header)[-1].split("\t")[9:]
ei=[i for i,s in enumerate(samples) if s in eur]
def dose_alt(g):a=g[:3];return (a[0]=='1')+(a[2]=='1') if(a[0]in'01'and a[2]in'01') else np.nan
geno={}
pos2snp={int(r.SNPPos):r.SNP for r in sub.itertuples()}
ea={r.SNP:r.AssessedAllele for r in sub.itertuples()}
for line in tb.fetch("2",11100000,11500000):
    f=line.split("\t");pos=int(f[1])
    if pos in pos2snp:
        ref,alt=f[3],f[4]
        if ","in alt: continue
        d=np.array([dose_alt(f[9:][i]) for i in ei],float)
        # code dosage as count of exposure effect allele
        snp=pos2snp[pos]
        if alt==ea[snp]: geno[snp]=d
        elif ref==ea[snp]: geno[snp]=2-d
        else: geno[snp]=d  # fallback
snps=[s for s in b13 if s in geno]
M=np.vstack([geno[s] for s in snps])
# signed Pearson r matrix (pairwise complete)
n=len(snps); R=np.eye(n)
for i in range(n):
    for j in range(i+1,n):
        x,y=M[i],M[j];m=~(np.isnan(x)|np.isnan(y))
        r=np.corrcoef(x[m],y[m])[0,1]; R[i,j]=R[j,i]=0 if np.isnan(r) else r
pd.DataFrame(R,index=snps,columns=snps).to_csv("ld_matrix_13.csv")
print("LD matrix built for %d SNPs (coded to exposure effect allele)"%n)
print("max |off-diag r|: %.3f  (r2: %.3f)"%(np.max(np.abs(R-np.eye(n))), np.max((R-np.eye(n))**2)))
