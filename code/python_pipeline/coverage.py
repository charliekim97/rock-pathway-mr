import pandas as pd, numpy as np
o=pd.read_csv("liverfat_rock2.tsv",sep="\t")
o=o.drop_duplicates("hm_rsid")
# exposure
b=pd.read_csv("rock2_eqtlgen_converted.txt",sep="\t")[["SNP","AssessedAllele","OtherAllele","beta","eaf"]]
b.columns=["SNP","exp_EA","exp_OA","exp_beta","exp_eaf"]
musc=pd.DataFrame({"SNP":["rs4335920","rs12996712"],"exp_EA":["A","T"],"exp_OA":["T","G"],
                   "exp_beta":[-0.187116,0.133808],"exp_eaf":[0.333,0.483]})  # GTEx ALT effect
# tiers
blood13=pd.read_csv("blood_clumped_independent.txt",sep="\t")["SNP"].tolist()
blood01=["rs12468344","rs13034888","rs7574150"]
lead="rs12468344"
def tier(s):
    t=[]
    if s==lead:t.append("BloodLead")
    if s in blood01:t.append("Blood_r2<.01")
    if s in blood13:t.append("Blood_r2<.1")
    if s in ["rs4335920","rs12996712"]:t.append("Muscle")
    return ",".join(t)
exp=pd.concat([b[b.SNP.isin(blood13+blood01)],musc],ignore_index=True).drop_duplicates("SNP")
m=exp.merge(o[["hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","hm_effect_allele_frequency","p_value"]],
            left_on="SNP",right_on="hm_rsid",how="left")
m["tier"]=m.SNP.map(tier)
def palin(ea,oa):
    return {ea,oa} in ({"A","T"},{"C","G"})
m["palindromic"]=[palin(r.exp_EA,r.exp_OA) for r in m.itertuples()]
m=m.sort_values("p_value")
cols=["SNP","tier","exp_EA","exp_OA","out_EA","out_OA","out_beta","out_se","out_eaf","out_p","palindromic"]
m=m.rename(columns={"hm_effect_allele":"out_EA","hm_other_allele":"out_OA","hm_beta":"out_beta",
                    "standard_error":"out_se","hm_effect_allele_frequency":"out_eaf","p_value":"out_p"})
pd.set_option("display.width",200,"display.max_columns",30)
print("=== OUTCOME (liver fat GCST90029073) coverage: %d/17 SNPs ==="%m.hm_rsid.notna().sum())
print(m[cols].to_string(index=False))
m[cols].to_csv("liverfat_harmonization_preview.txt",sep="\t",index=False)
print("\n=== TIER COVERAGE ===")
for name,snps in [("Blood lead (Primary)",[lead]),("Blood r2<0.01 (Supporting IVW)",blood01),
                  ("Blood r2<0.1 (LD-corr sensitivity)",blood13),("Muscle (2-SNP)",["rs4335920","rs12996712"])]:
    have=m[m.SNP.isin(snps)].hm_rsid.notna().sum()
    print(f"  {name}: {have}/{len(snps)}")
print("\n=== rs4335920 PALINDROMIC strand check ===")
r=m[m.SNP=="rs4335920"].iloc[0]
print(f"  exposure(GTEx muscle): EA=A(ALT) freq~0.353(1000G EUR)/0.333(GTEx)")
print(f"  outcome(liverfat): EA={r.out_EA} OA={r.out_OA} EAF={r.out_eaf:.3f} beta={r.out_beta} p={r.out_p}")
print(f"  -> AF concordant if outcome A-allele freq ~0.35 (far from 0.5 => strand resolvable)")
