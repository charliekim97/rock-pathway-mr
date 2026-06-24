import urllib.request,json,csv,sys,math,time,os
from statistics import NormalDist
nd=NormalDist()
def zfp(p):
    if p<=0:p=1e-300
    a=1-p/2
    return math.sqrt(max(2*(-math.log(p))-math.log(2*math.pi*(-math.log(p))),1e-6)) if a>=1 else nd.inv_cdf(a)
name,chrom,L,U,acc=sys.argv[1:6]
tag=f"{name}_{L}"; csvf=f"out_{tag}.csv"; statef=f"out_{tag}.next"
url=None
if os.path.exists(statef):
    c=open(statef).read().strip()
    if c=="DONE": print(f"{name}: DONE"); sys.exit()
    if c.startswith("http"): url=c
if url is None:
    url=f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations?study_accession={acc}&bp_lower={L}&bp_upper={U}&size=1000"
    with open(csvf,"w",newline="") as f: csv.writer(f).writerow(["SNP","EA","OA","beta","se","p"])
t0=time.time()
while url and time.time()-t0<36:
    d=None
    for _ in range(6):
        try:
            with urllib.request.urlopen(url,timeout=25) as r: d=json.load(r); break
        except Exception: time.sleep(1.5)
    if d is None: open(statef,"w").write(url); print(f"{name}: api-stall, saved, rerun"); sys.exit(1)
    rows=[]
    for k,v in d.get('_embedded',{}).get('associations',{}).items():
        rs,b,p=v.get('variant_id'),v.get('beta'),v.get('p_value')
        if rs and b is not None and p is not None:
            se=abs(b)/zfp(p) if (b!=0 and 0<p<1) else float('nan')
            rows.append([rs,v.get('effect_allele'),v.get('other_allele'),b,se,p])
    with open(csvf,"a",newline="") as f: csv.writer(f).writerows(rows)
    url=d.get('_links',{}).get('next',{}).get('href')
    open(statef,"w").write(url if url else "DONE")

n=sum(1 for _ in open(csvf))-1
print(f"{name}: {n} SNPs, {'DONE' if not url else 'more(rerun)'}")
