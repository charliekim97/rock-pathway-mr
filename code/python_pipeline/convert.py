import pandas as pd, numpy as np
e = pd.read_csv("rock2_eqtlgen.txt", sep="\t")
af = pd.read_csv("rock2_af.txt", sep="\t")[["SNP","AlleleA","AlleleB","AlleleB_all"]]
m = e.merge(af, on="SNP", how="inner")
# align EAF to AssessedAllele
def eaf(r):
    if r.AssessedAllele == r.AlleleB: return r.AlleleB_all
    if r.AssessedAllele == r.AlleleA: return 1 - r.AlleleB_all
    return np.nan
m["eaf"] = m.apply(eaf, axis=1)
bad = m["eaf"].isna().sum()
m = m.dropna(subset=["eaf"]).copy()
m["n"] = m["NrSamples"]
Z = m["Zscore"]; p = m["eaf"]; N = m["n"]
denom = np.sqrt(2*p*(1-p)*(N + Z**2))
m["beta"] = Z/denom
m["se"]   = 1/denom
m["F"]    = (m["beta"]/m["se"])**2   # == Z^2
m = m.sort_values("F", ascending=False)
m.to_csv("rock2_eqtlgen_converted.txt", sep="\t", index=False)
print("converted SNPs:", len(m), "| allele-unmatched dropped:", bad)
print("F range: %.1f – %.1f | all F>10: %s" % (m["F"].min(), m["F"].max(), bool((m["F"]>10).all())))
cols=["SNP","SNPPos","AssessedAllele","OtherAllele","eaf","n","Zscore","beta","se","F","Pvalue"]
print("\n=== TOP 5 by F ===")
print(m[cols].head(5).to_string(index=False))
print("\n=== rs7581184 (packet blood instrument) ===")
sub=m[m.SNP=="rs7581184"]
print(sub[cols].to_string(index=False) if len(sub) else "NOT in significant set")
