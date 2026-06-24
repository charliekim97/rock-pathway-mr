import urllib.request,json,csv,sys,math,os,time
from statistics import NormalDist
nd=NormalDist()
def zfp(p):
    if p<=0:p=1e-300
    a=1-p/2
    if a>=1.0:
        L=-math.log(p);return math.sqrt(max(2*L-math.log(2*math.pi*L),1e-6))
    return nd.inv_cdf(a)
name,chrom,L,U,acc=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
csvf=f"out_{name}.csv"; statef=f"out_{name}.next"
if os.path.exists(statef):
    url=open(statef).read().strip()
    if url=="DONE": print(name,"already DONE"); sys.exit()
else:
    url=f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations?study_accession={acc}&bp_lower={L}&bp_upper={U}&size=1000"
    with open(csvf,"w",newline="") as f: csv.writer(f).writerow(["SNP","EA","OA","beta","se","p"])
t0=time.time();pages=0
while url and time.time()-t0<35:
    d=None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url,timeout=30) as r: d=json.load(r); break
        except Exception as e:
            time.sleep(2)
    if d is None: print("err after retries");break
    a=d.get('_embedded',{}).get('associations',{}); rows=[]
    for k,v in a.items():
        rs,b,p=v.get('variant_id'),v.get('beta'),v.get('p_value')
        if rs and b is not None and p is not None:
            se=abs(b)/zfp(p) if (b!=0 and 0<p<1) else float('nan')
            rows.append([rs,v.get('effect_allele'),v.get('other_allele'),b,se,p])
    with open(csvf,"a",newline="") as f: csv.writer(f).writerows(rows)
    pages+=1
    url=d.get('_links',{}).get('next',{}).get('href')
open(statef,"w").write(url if url else "DONE")
n=sum(1 for _ in open(csvf))-1
print(f"{name}: +{pages}pages, total {n} SNPs, {'DONE' if not url else 'more'}")
