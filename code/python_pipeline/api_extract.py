import urllib.request, json, math
from statistics import NormalDist
nd=NormalDist()
snps=[l.strip() for l in open("instrument_snps_for_outcome.txt") if l.strip()]
outcomes={
 "FG_EUR":"GCST90002232","FI_EUR":"GCST90002238","HbA1c_EUR":"GCST90002244",
 "2hGlu_EUR":"GCST90002227","T2D_Xue2018":"GCST006867",
}
def fetch(rs,g):
    url=f"https://www.ebi.ac.uk/gwas/summary-statistics/api/associations/{rs}?study_accession={g}"
    try:
        with urllib.request.urlopen(url,timeout=20) as r:
            d=json.load(r)
        if "error" in d: return None
        return d
    except Exception:
        return None
def se_from(beta,p):
    try:
        if p<=0 or p>=1 or beta==0: return float('nan')
        z=nd.inv_cdf(1-p/2); return abs(beta)/z
    except: return float('nan')
rows=[]
cov={o:0 for o in outcomes}
for o,g in outcomes.items():
    for rs in snps:
        d=fetch(rs,g)
        if d:
            cov[o]+=1
            b=d.get('beta'); p=d.get('p_value')
            rows.append([o,rs,d.get('effect_allele'),d.get('other_allele'),
                         b,se_from(b,p) if b is not None and p else None,
                         d.get('effect_allele_frequency'),p])
        else:
            rows.append([o,rs,None,None,None,None,None,None])
import csv
with open("api_outcomes_long.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t");w.writerow(["outcome","SNP","EA","OA","beta","se","eaf","p"]);w.writerows(rows)
print("=== API coverage (of 17) ===")
for o in outcomes: print(f"  {o:14s} {cov[o]}/17  ({outcomes[o]})")
