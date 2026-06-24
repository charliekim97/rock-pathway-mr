import urllib.request,json,csv,sys,math,time
from statistics import NormalDist
nd=NormalDist()
def zfp(p):
    if p<=0:p=1e-300
    a=1-p/2
    return math.sqrt(max(2*(-math.log(p))-math.log(2*math.pi*(-math.log(p))),1e-6)) if a>=1 else nd.inv_cdf(a)
name,chrom,L,U,acc=sys.argv[1:6]
csvf=f"out_{name}.csv"
url=f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations?study_accession={acc}&bp_lower={L}&bp_upper={U}&size=1000"
rows={}; t0=time.time()
while url and time.time()-t0<38:
    d=None
    for _ in range(6):
        try:
            with urllib.request.urlopen(url,timeout=25) as r: d=json.load(r); break
        except Exception: time.sleep(1.5)
    if d is None: print(f"{name}: FAILED (api)"); sys.exit(1)
    for k,v in d.get('_embedded',{}).get('associations',{}).items():
        rs,b,p=v.get('variant_id'),v.get('beta'),v.get('p_value')
        if rs and b is not None and p is not None:
            se=abs(b)/zfp(p) if (b!=0 and 0<p<1) else float('nan')
            rows[rs]=[rs,v.get('effect_allele'),v.get('other_allele'),b,se,p]
    url=d.get('_links',{}).get('next',{}).get('href')
with open(csvf,"w",newline="") as f:
    csv.writer(f).writerow(["SNP","EA","OA","beta","se","p"]); csv.writer(f).writerows(rows.values())
print(f"{name}: {len(rows)} SNPs, {'DONE' if not url else 'PARTIAL(rerun)'}")
