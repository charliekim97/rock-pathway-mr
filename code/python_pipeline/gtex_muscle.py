import math
from statistics import NormalDist
nd = NormalDist()
def z_from_p(p):  # two-sided
    return nd.inv_cdf(1 - p/2)

# GTEx v8 Muscle_Skeletal independent cis-eQTL (official endpoint values)
musc = [
 # snp, rsid, variantId, nes(slope), p_indep, maf, dyn_nes, dyn_se(error), dyn_t
 ("rs4335920","chr2_11415777_T_A_b38",-0.187116,6.65406e-16,0.332861,-0.17783529,0.02321254,-7.661174),
 ("rs12996712","chr2_11261678_G_T_b38",0.133808,5.81572e-10,0.483003, 0.12443030,0.02231840, 5.575),
]
rows=[]
for snp,vid,nes,p,maf,dnes,dse,dt in musc:
    Z=z_from_p(p)
    se_off=abs(nes)/Z          # SE consistent with official slope+p
    F_off=Z**2
    F_dyn=(dnes/dse)**2
    rows.append((snp,nes,se_off,maf,p,F_off,F_dyn))
    print(f"{snp:11s} beta(nes)={nes:+.4f}  se={se_off:.5f}  maf={maf:.3f}  p={p:.2e}  F_official={F_off:.1f}  F_dyneqtl={F_dyn:.1f}")

# write muscle instrument file
with open("rock2_gtex_muscle_instrument.txt","w") as f:
    f.write("SNP\tbeta\tse\teaf\tpval\tFstat\ttissue\n")
    for snp,nes,se,maf,p,Fo,Fd in rows:
        f.write(f"{snp}\t{nes}\t{se:.6f}\t{maf}\t{p}\t{Fo:.2f}\tMuscle_Skeletal\n")
print("\nwritten: rock2_gtex_muscle_instrument.txt")
