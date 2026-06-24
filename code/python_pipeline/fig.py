import os
os.environ["MPLCONFIGDIR"]="/tmp/mplcfg"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none",
    "font.family":"DejaVu Sans", "font.size":9,
    "axes.linewidth":0.8, "axes.edgecolor":"#333333",
    "savefig.dpi":300, "figure.dpi":150,
})
OUT="/sessions/determined-eager-davinci/mnt/outputs/rock2_mr/figs"
# Okabe-Ito colorblind-safe
C=dict(blue="#0072B2",orange="#E69F00",green="#009E73",purple="#CC79A7",
       vermillion="#D55E00",sky="#56B4E9",yellow="#F0E442",grey="#999999",black="#000000")
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.png",bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight")
    plt.close(fig)

OUTCOMES=["HepaticFat","FI","FG","HbA1c","2hGlu","T2D","BMI","Cirrhosis","NAFLDmeta","FibrK74"]
FIBROSIS={"Cirrhosis","NAFLDmeta","FibrK74"}

# ---------------- FIG 1 : pipeline + F-stat ----------------
fig=plt.figure(figsize=(11.5,5.0))
gs=fig.add_gridspec(1,2,width_ratios=[1.15,1.0],wspace=0.30)
axp=fig.add_subplot(gs[0]); axp.axis("off"); axp.set_xlim(0,10); axp.set_ylim(-0.45,10.1)
def box(yc,t,fc,h=1.35):
    w,x=9.4,0.3
    axp.add_patch(FancyBboxPatch((x,yc-h/2),w,h,boxstyle="round,pad=0.05,rounding_size=0.12",
        fc=fc,ec="#333333",lw=0.9))
    axp.text(x+w/2,yc,t,ha="center",va="center",fontsize=7.6)
stages=[
 (9.25,"Exposure: cis-eQTL in 5 tissues — eQTLGen blood (N≤31,684);\nGTEx muscle, liver, subcutaneous & visceral adipose","#DDEBF7"),
 (7.55,"Targets: ROCK2 / ROCK1 / RHOA  (RhoA–ROCK pathway)","#E2F0D9"),
 (5.95,"Instruments: Z→β/SE (Zhu) · F-statistic · LD clump (1000G EUR); tiered","#FFF2CC"),
 (4.10,"10 primary outcomes (7 cardiometabolic + 3 fibrosis/liver);\nT2D = DIAMANTE EUR / Finnish / trans-ancestry; Chen NAFLD replication","#FFE3D5"),
 (2.30,"MR (Wald / IVW / LD-corrected) + Steiger  ·  coloc.abf + coloc.susie","#FCE4D6"),
 (0.60,"Robustness: multi-tissue · prior/power · Genebass pLoF ·\npositive controls (KCNJ11→T2D, SORT1→LDL liver)","#EDEDED"),
]
for yc,t,fc in stages: box(yc,t,fc)
for i in range(len(stages)-1):
    axp.add_patch(FancyArrowPatch((5,stages[i][0]-0.70),(5,stages[i+1][0]+0.70),
        arrowstyle="-|>",mutation_scale=11,lw=1.2,color="#555555"))
axp.set_title("a  Study design & pipeline",loc="left",fontsize=11,fontweight="bold")
# F-stat panel
axf=fig.add_subplot(gs[1])
inst=[("ROCK2 blood (rs12468344)",409,C["blue"]),
      ("ROCK2 muscle (rs4335920)",65,C["orange"]),
      ("ROCK2 muscle (rs12996712)",38,C["orange"]),
      ("ROCK1 blood (rs142716063)",86,C["green"]),
      ("RHOA blood (rs115431681)",27,C["purple"])]
y=np.arange(len(inst))[::-1]
for yi,(lab,F,col) in zip(y,inst):
    axf.barh(yi,F,color=col,height=0.46,ec="#333",lw=0.6)
    axf.text(4,yi+0.30,lab,va="bottom",ha="left",fontsize=7.4)   # name above bar (inside panel)
    axf.text(F+8,yi,f"F={F}",va="center",fontsize=8)
axf.axvline(10,ls="--",color=C["vermillion"],lw=1.2)
axf.text(10,len(inst)-0.35,"F=10",color=C["vermillion"],fontsize=7.5,ha="center",va="bottom")
axf.set_yticks([]); axf.set_ylim(-0.6,len(inst)-0.2)
axf.set_xlim(0,480)
axf.set_xlabel("Instrument strength (F-statistic)")
axf.set_title("b  Instrument strength",loc="left",fontsize=11,fontweight="bold")
for s in ["top","right","left"]: axf.spines[s].set_visible(False)
save(fig,"Fig1_design_instruments")
print("Fig1 done")

# ---------------- FIG 2 : MR forest ----------------
# data: gene -> outcome -> list of (tier,label,color, beta,lo,hi)
TIER={"lead":("Blood lead (Wald)",C["blue"]),"muscle":("Muscle IVW",C["orange"]),
      "r01":("Blood r²<0.01 IVW",C["green"]),"r1":("Blood r²<0.1 LD-corr",C["purple"]),
      "ivw":("Blood r²<0.01 IVW",C["green"]),"wald":("Blood Wald (1 SNP)",C["blue"])}
ROCK2={
 "HepaticFat":[("lead",-0.029,-0.124,0.065),("muscle",-0.030,-0.098,0.038),("r01",-0.008,-0.096,0.080),("r1",0.004,-0.078,0.085)],
 "FI":[("lead",-0.016,-0.037,0.005),("muscle",-0.005,-0.020,0.011),("r01",-0.011,-0.032,0.009),("r1",-0.009,-0.027,0.009)],
 "FG":[("lead",-0.009,-0.027,0.009),("r01",-0.011,-0.026,0.005),("r1",-0.004,-0.016,0.008)],
 "HbA1c":[("lead",0.008,-0.004,0.020),("muscle",0.004,-0.008,0.017),("r01",0.007,-0.003,0.018),("r1",0.012,0.001,0.023)],
 "2hGlu":[("lead",0.010,-0.140,0.160),("muscle",0.025,-0.176,0.226),("r01",-0.007,-0.144,0.131),("r1",0.018,-0.120,0.155)],
 "T2D":[("lead",0.021,-0.058,0.099),("muscle",-0.012,-0.068,0.044),("r01",0.046,-0.026,0.119)],
 "BMI":[("lead",0.037,0.014,0.060),("muscle",-0.000,-0.017,0.017),("r01",0.039,0.017,0.060),("r1",0.025,0.006,0.045)],
 "Cirrhosis":[("lead",-0.099,-0.336,0.139),("muscle",-0.182,-0.351,-0.013),("r01",-0.132,-0.352,0.088)],
 "NAFLDmeta":[("lead",-0.091,-0.291,0.109),("muscle",-0.031,-0.176,0.114),("r01",-0.109,-0.294,0.077)],
 "FibrK74":[("lead",0.053,-0.288,0.395),("muscle",-0.009,-0.253,0.234),("r01",0.026,-0.290,0.342)],
 "NASH":[("lead",1.123,0.021,2.224),("muscle",0.192,-0.592,0.976),("r01",0.676,-0.343,1.694)],
}
ROCK1={
 "HepaticFat":[("lead",-0.131,-0.259,-0.003),("ivw",-0.143,-0.265,-0.021)],
 "FI":[("lead",0.038,-0.031,0.106),("ivw",0.009,-0.043,0.061)],
 "FG":[("lead",0.008,-0.042,0.058),("ivw",-0.015,-0.055,0.024)],
 "HbA1c":[("lead",-0.022,-0.065,0.020),("ivw",-0.032,-0.070,0.005)],
 "2hGlu":[("lead",0.013,-0.113,0.140),("ivw",0.007,-0.117,0.130)],
 "T2D":[("lead",0.012,-0.112,0.135),("ivw",0.009,-0.106,0.125)],
 "BMI":[("lead",-0.005,-0.036,0.027),("ivw",-0.005,-0.035,0.025)],
 "Cirrhosis":[("lead",-0.107,-0.743,0.529)],"NAFLDmeta":[],"FibrK74":[],"NASH":[],
}
RHOA={
 "HepaticFat":[("wald",-0.092,-0.229,0.046)],"FI":[("wald",-0.014,-0.052,0.025)],
 "FG":[("wald",0.040,-0.014,0.094)],"HbA1c":[("wald",0.008,-0.166,0.181)],
 "2hGlu":[("wald",0.211,0.002,0.420)],"T2D":[("wald",-0.044,-0.16,0.072)],"BMI":[("wald",0.012,-0.023,0.047)],
 "Cirrhosis":[("wald",0.159,-0.257,0.576)],"NAFLDmeta":[],"FibrK74":[],"NASH":[],
}
NOMINAL={("ROCK2","BMI"),("ROCK1","HepaticFat"),("ROCK2","Cirrhosis")}
genes=[("ROCK2",ROCK2),("ROCK1",ROCK1)]
fig,axes=plt.subplots(1,2,figsize=(10.5,7.2),sharey=True)
for ax,(gname,gd) in zip(axes,genes):
    ax.axvline(0,color="#888",lw=1.0,zorder=1)
    seen=set()
    for oi,oc in enumerate(OUTCOMES):
        base=len(OUTCOMES)-1-oi
        tiers=gd[oc]
        n=len(tiers)
        offs=np.linspace(0.26,-0.26,n) if n>1 else [0.0]
        for off,(tier,b,lo,hi) in zip(offs,tiers):
            lab,col=TIER[tier]
            ax.errorbar(b,base+off,xerr=[[b-lo],[hi-b]],fmt="o",ms=5.6,color=col,
                ecolor=col,elinewidth=1.7,capsize=2.8,zorder=3,
                label=lab if lab not in seen else None)
            seen.add(lab)
        if (gname,oc) in NOMINAL:
            ax.text(0.965,base,"★",color=C["vermillion"],fontsize=12,va="center",ha="right",
                    transform=ax.get_yaxis_transform())
    ax.axhline(2.5,color="#cccccc",lw=0.8,ls=(0,(4,3)))   # divider: cardiometabolic | fibrosis
    ax.set_yticks(range(len(OUTCOMES))); ax.set_yticklabels(OUTCOMES[::-1],fontsize=10.5)
    ax.set_xlim(-0.7,0.7); ax.set_xlabel("MR estimate (per-SD expression)\nβ (continuous) or log-OR (binary) [95% CI]",fontsize=10)
    ax.set_title(gname,fontweight="bold",fontsize=13)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.set_ylim(-0.6,len(OUTCOMES)-0.4)
    if gname=="ROCK2":
        ax.text(-0.69,6.55,"cardiometabolic",fontsize=6.5,color="#999",style="italic")
        ax.text(-0.69,2.62,"fibrosis / liver",fontsize=6.5,color="#999",style="italic")
handles,labels=[],[]
for ax in axes:
    for h,l in zip(*ax.get_legend_handles_labels()):
        if l not in labels: handles.append(h); labels.append(l)
fig.legend(handles,labels,fontsize=8,frameon=False,ncol=4,loc="lower center",bbox_to_anchor=(0.5,0.02))
fig.text(0.5,-0.02,"Continuous traits are per-SD expression (β); binary traits (T2D, fibrosis/liver) are log-OR. T2D = DIAMANTE-EUR (80,154 cases); replicated in FinnGen (Finnish) and DIAMANTE trans-ancestry. Exploratory RHOA shown in Supplementary Fig. S2.",
         fontsize=7,color="#666",ha="center")
fig.text(0.5,-0.05,"★ nominal p<0.05; none supported by colocalization. ROCK2 muscle–cirrhosis is tissue-discordant (blood null) and negative on muscle colocalization (PP4=0.032).",
         fontsize=7,color="#666",ha="center")
fig.suptitle("cis-eQTL MR of ROCK2 and ROCK1 expression vs 10 primary cardiometabolic and fibrosis/liver outcomes",
             fontsize=12,fontweight="bold",y=1.0)
fig.subplots_adjust(bottom=0.16)
save(fig,"Fig2_MR_forest")
print("Fig2 done")

# ---------------- FIG 3 : coloc heatmap + positive control ----------------
nan=np.nan
M=np.array([
 [0.005,0.003,0.001,0.001,0.005,0.003,0.045, 0.010,0.008,0.010],
 [0.070,0.018,0.013,0.007,0.012,0.016,0.004, 0.055,nan,  nan ],
])
genenames=["ROCK2","ROCK1"]
import matplotlib.cm as cm
cmap=cm.get_cmap("magma").copy(); cmap.set_bad("#DInvalid".replace("DInvalid","E8E8E8"))
fig=plt.figure(figsize=(13.5,5.6))
gs=fig.add_gridspec(1,2,width_ratios=[3.3,0.95],wspace=0.22)
ax=fig.add_subplot(gs[0])
Mm=np.ma.masked_invalid(M)
im=ax.imshow(Mm,cmap=cmap,vmin=0,vmax=1,aspect="auto")
ax.set_xticks(range(10)); ax.set_xticklabels(OUTCOMES,rotation=40,ha="right",fontsize=10)
ax.set_yticks(range(2)); ax.set_yticklabels(genenames,fontsize=12)
ax.axvline(6.5,color="white",lw=2.0)  # divider cardiometabolic | fibrosis
for i in range(2):
    for j in range(10):
        v=M[i,j]
        if np.isnan(v): ax.text(j,i,"n.t.",ha="center",va="center",color="#999",fontsize=9)
        else: ax.text(j,i,f"{v:.3f}",ha="center",va="center",color="white" if v<0.55 else "black",fontsize=9.5)
ax.set_title("a  Colocalization PP4 (shared causal variant) — ROCK1 and ROCK2 × 10 primary outcomes",loc="left",fontweight="bold",fontsize=10)
cb=fig.colorbar(im,ax=ax,fraction=0.030,pad=0.02); cb.set_label("PP4")
# positive control panel
axp=fig.add_subplot(gs[1])
pc=[("KCNJ11→T2D\n(blood)",0.888,C["green"]),("SORT1→LDL\n(liver)",1.000,C["green"]),("SORT1→LDL\n(blood)",0.000,C["grey"])]
xs=np.arange(len(pc))
axp.bar(xs,[p[1] for p in pc],color=[p[2] for p in pc],ec="#333",lw=0.6,width=0.6)
for x,(lab,v,c) in zip(xs,pc): axp.text(x,v+0.03,f"{v:.3f}",ha="center",fontsize=7.5)
axp.axhline(0.75,ls="--",color=C["vermillion"],lw=1.0); axp.text(2.0,0.78,"PP4=0.75",color=C["vermillion"],fontsize=7)
axp.set_xticks(xs); axp.set_xticklabels([p[0] for p in pc],fontsize=7.5)
axp.set_ylim(0,1.0); axp.set_ylabel("PP4")
axp.set_title("b  Positive controls",loc="left",fontweight="bold",fontsize=10.5)
for s in ["top","right"]: axp.spines[s].set_visible(False)
fig.text(0.5,0.00,"Every PP4 < 0.11 across all ROCK1/ROCK2 gene×trait colocalizations.  T2D = DIAMANTE-EUR (80,154 cases).  n.t. = not tested (ROCK1 fibrosis beyond cirrhosis).  "
         "Exploratory single-instrument RHOA is shown in Supplementary Fig. S2.  (b) Positive controls require the mechanistically correct tissue: SORT1→LDL colocalizes in liver (1.00) not blood (0.00); KCNJ11→T2D in blood (0.89).",
         fontsize=6.4,color="#555",ha="center")
fig.subplots_adjust(bottom=0.30)
save(fig,"Fig3_coloc_heatmap")
print("Fig3 done")

# ---------------- FIG 4 : nominal hits PP3 vs PP4 ----------------
hits=[("ROCK2 × BMI\n(blood IVW)",0.605,0.045),("ROCK1 × HepaticFat",0.033,0.070),("RHOA × 2hGlu",0.058,0.061),("ROCK2 × Cirrhosis\n(muscle IVW)",0.042,0.032)]
fig,ax=plt.subplots(figsize=(8.4,4.4))
x=np.arange(len(hits)); w=0.34
ax.bar(x-w/2,[h[1] for h in hits],w,label="PP3 (distinct variants)",color=C["sky"],ec="#333",lw=0.6)
ax.bar(x+w/2,[h[2] for h in hits],w,label="PP4 (shared causal variant)",color=C["vermillion"],ec="#333",lw=0.6)
for xi,h in zip(x,hits):
    ax.text(xi-w/2,h[1]+0.015,f"{h[1]:.3f}",ha="center",fontsize=7.5)
    ax.text(xi+w/2,h[2]+0.015,f"{h[2]:.3f}",ha="center",fontsize=7.5)
ax.axhline(0.75,ls="--",color="#333",lw=0.9)
ax.text(x[-1]+0.5,0.755,"coloc threshold 0.75",fontsize=7,ha="right",va="bottom")
ax.set_xticks(x); ax.set_xticklabels([h[0] for h in hits],fontsize=8.5)
ax.set_xlim(-0.6,len(hits)-0.4); ax.set_ylim(0,0.82); ax.set_ylabel("Posterior probability")
ax.set_title("Nominal MR associations lack evidence of a shared causal signal",fontweight="bold",fontsize=10.5,pad=10)
ax.legend(fontsize=8,frameon=False,ncol=2,loc="lower center",bbox_to_anchor=(0.5,-0.28))
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.subplots_adjust(bottom=0.24)
save(fig,"Fig4_nominal_hits_rejected")
print("Fig4 done")

# ---------------- FIG 6 : ROCK2 tissue x outcome coloc ----------------
TIS=["Blood (eQTLGen)","Muscle (GTEx)","Liver (GTEx)","Adipose-SC (GTEx)","Adipose-Visc (GTEx)"]
OC6=["Hepatic fat","NAFLD","Cirrhosis","BMI"]
T=np.array([
 [0.005,0.008,0.010,0.045],
 [0.008,0.009,0.032,0.010],
 [0.006,0.004,0.005,0.005],
 [0.008,0.007,0.013,0.079],
 [0.005,0.004,0.007,0.010],
])
fig,ax=plt.subplots(figsize=(5.6,3.6))
im=ax.imshow(T,cmap="magma",vmin=0,vmax=1,aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(OC6,rotation=25,ha="right",fontsize=8.5)
ax.set_yticks(range(5)); ax.set_yticklabels(TIS,fontsize=8.5)
for i in range(5):
    for j in range(4):
        ax.text(j,i,f"{T[i,j]:.3f}",ha="center",va="center",color="white",fontsize=8)
ax.set_title("ROCK2 colocalization (PP4) across exposure tissues",fontweight="bold",fontsize=10.5,pad=8)
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.set_label("PP4")
fig.text(0.5,-0.04,"ROCK2 cis-eQTL has usable instruments in all five tissues (F = 13.5–409); no colocalization with any metabolic/liver outcome (all PP4 ≤ 0.08).",
         fontsize=7,color="#555",ha="center")
fig.subplots_adjust(bottom=0.30)
save(fig,"Fig6_tissue_coloc")
print("Fig6 done")

# ---------------- FIG S3 : RHOA exploratory (single weak instrument) ----------------
rh=[("Hepatic fat",-0.092,-0.229,0.046),("Fasting insulin",-0.014,-0.052,0.025),
    ("Fasting glucose",0.040,-0.014,0.094),("HbA1c",0.008,-0.166,0.181),
    ("2-h glucose",0.211,0.002,0.420),("T2D",-0.044,-0.160,0.072),
    ("BMI",0.012,-0.023,0.047),("Cirrhosis",0.159,-0.257,0.576)]
fig,ax=plt.subplots(figsize=(6.4,4.2))
ax.axvline(0,color="#888",lw=1.0)
y=np.arange(len(rh))[::-1]
for yi,(lab,b,lo,hi) in zip(y,rh):
    ax.errorbar(b,yi,xerr=[[b-lo],[hi-b]],fmt="o",ms=4.6,color=C["purple"],ecolor=C["purple"],
                elinewidth=1.3,capsize=2.5)
ax.set_yticks(y[::-1]); ax.set_yticklabels([r[0] for r in rh][::-1],fontsize=8.5)
ax.set_xlim(-0.7,0.7); ax.set_xlabel("RHOA MR estimate (per-SD expression) [95% CI]")
ax.set_title("Supplementary Fig. S2  RHOA: exploratory single-instrument analysis",fontweight="bold",fontsize=10,pad=8)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.text(0.5,-0.04,"RHOA had a single genome-wide-significant blood cis-eQTL (rs115431681, F=27); estimates are Wald ratios and are exploratory. "
         "No outcome reached the colocalization threshold. In DIAMANTE-EUR, the RHOA cis-window contained a genome-wide T2D signal "
         "(regional p=2×10⁻⁹) that did NOT colocalize with the RHOA eQTL (PP4=0.02, PP3=0.52) — a positive demonstration that the pipeline "
         "returns a distinct-variant verdict where regional outcome signal is genuinely present.",
         fontsize=6.6,color="#555",ha="center",wrap=True)
fig.subplots_adjust(bottom=0.26)
save(fig,"FigS3_RHOA")
print("FigS3 done")
print("ALL FIGURES SAVED to",OUT)
