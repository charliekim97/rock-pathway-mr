import os
os.environ["MPLCONFIGDIR"]="/tmp/mplcfg"
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"DejaVu Sans","savefig.dpi":300})
OUT="/sessions/determined-eager-davinci/mnt/outputs/rock2_mr/figs"
C=dict(meas="#56B4E9", drug="#E69F00", kin="#009E73", dis="#CC79A7", concl="#D55E00", grey="#888888")

fig=plt.figure(figsize=(11,5.4)); ax=fig.add_subplot(111); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,10)
def box(x,y,w,h,t,fc,fs=8.4):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.06,rounding_size=0.12",fc=fc,ec="#333",lw=0.9))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs)
def arr(x1,y1,x2,y2,style="-|>",col="#555",lw=1.3):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=12,lw=lw,color=col))

ax.text(0.1,9.7,"Conceptual limits of expression-based genetic validation for intracellular kinases",
        fontsize=12,fontweight="bold")

# LEFT branch — what we measured (genetics)
box(0.3,7.7,4.3,1.25,"A  cis-eQTL  →  ROCK2 transcript abundance\n(lifelong, small-magnitude; what this study proxies)",C["meas"])
arr(2.45,7.7,2.45,6.95)
box(0.3,5.5,4.3,1.25,"Germline variant effect on\nROCK2 mRNA level  (blood, muscle, liver, adipose)",C["meas"],8.0)

# RIGHT branch — drug action (pharmacology)
box(5.4,7.7,4.3,1.25,"B  Drug exposure  →  intracellular free [inhibitor]\n(acute, dose-dependent; what ROCK2 drugs do)",C["drug"])
arr(7.55,7.7,7.55,6.95)
box(5.4,5.5,4.3,1.25,"C  ATP-competitive inhibition of catalytic activity\n(modulated by ATP 1–10 mM, enzyme abundance, substrate, context)",C["drug"],7.8)

# converge to kinase output
arr(2.45,5.5,4.6,4.55); arr(7.55,5.5,5.4,4.55)
box(2.7,3.5,4.6,1.0,"D  ROCK2 kinase output (substrate phosphorylation)\n→ disease-relevant signalling",C["kin"],8.4)
arr(5.0,3.5,5.0,2.75)
box(2.7,1.7,4.6,1.0,"Disease (cardiometabolic / fibrotic) phenotype",C["dis"],8.6)

# bottom conclusion bar
box(0.3,0.2,9.4,1.0,
    "E  Expression and activity are linked but non-identical perturbations: a null for genetically proxied expression\n"
    "does not, by itself, refute pharmacologic inhibition of kinase activity.",C["concl"],8.6)

# annotations
ax.text(0.3,6.6,"measured here",fontsize=7,color=C["grey"],style="italic")
ax.text(9.7,6.6,"not measured",fontsize=7,color=C["grey"],style="italic",ha="right")
fig.savefig(f"{OUT}/Fig5_conceptual_limits.png",bbox_inches="tight")
fig.savefig(f"{OUT}/Fig5_conceptual_limits.pdf",bbox_inches="tight")
print("Fig5 done")
