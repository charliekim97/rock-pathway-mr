import urllib.request,json,csv,sys,math
from statistics import NormalDist
nd=NormalDist()
def zfp(p):
    if p<=0: p=1e-300
    arg=1-p/2
    if arg>=1.0: 
        L=-math.log(p); return math.sqrt(max(2*L-math.log(2*math.pi*L),1e-6))
    return nd.inv_cdf(arg)
def fetch_region(chrom,L,U,acc,out):
    base=f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations?study_accession={acc}&bp_lower={L}&bp_upper={U}&size=1000"
    url=base; rows={}; pages=0
    while url and pages<60:
        try:
            with urllib.request.urlopen(url,timeout=30) as r: d=json.load(r)
        except Exception as ex:
            print("  err",ex); break
        a=d.get('_embedded',{}).get('associations',{})
        for k,v in a.items():
            rs=v.get('variant_id'); b=v.get('beta'); p=v.get('p_value')
            if rs and b is not None and p is not None:
                se=abs(b)/zfp(p) if (b!=0 and 0<p<1) else float('nan')
                rows[rs]=(rs,v.get('effect_allele'),v.get('other_allele'),b,se,p,v.get('effect_allele_frequency'))
        pages+=1
        url=d.get('_links',{}).get('next',{}).get('href')
    with open(out,"w",newline="") as f:
        w=csv.writer(f);w.writerow(["SNP","EA","OA","beta","se","p","eaf"]);w.writerows(rows.values())
    print(f"  {out}: {len(rows)} SNPs ({pages} pages)")
    return len(rows)

jobs=eval(sys.argv[1])
for name,chrom,L,U,acc in jobs:
    print(name); fetch_region(chrom,L,U,acc,f"out_{name}.csv")
