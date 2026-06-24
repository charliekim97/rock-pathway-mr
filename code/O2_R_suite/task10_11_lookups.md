# Tasks 10 & 11 — web lookups (no script; 5–10 min each)

## Task 10 — Genebass rare-variant pLoF burden (orthogonal perturbation)
Go to **https://genebass.org** (UKB exome, ~400k). For each gene (ROCK2, ROCK1, RHOA):
1. Search gene → open the gene page.
2. Annotation = **pLoF** (and optionally pLoF+missense). Burden test.
3. Record association for: BMI, T2D, ALT, AST, "non-alcoholic fatty liver disease", cirrhosis/fibrosis (if present).
4. Tabulate: gene × phenotype → beta/OR, p, N carriers (cases).

Interpretation: **null** → one Supplementary line ("rare-variant pLoF burden underpowered / no signal"); do not over-push.
A signal would be interesting (note carrier N — usually tiny). Output: one table `task10_genebass.csv` (hand-entered).

## Task 11 — ROCK2 sQTL feasibility (activity-adjacent proxy)
eQTL Catalogue (https://www.ebi.ac.uk/eqtl/), splicing QTL (quant_method = "leafcutter" / "txrevise"):
1. Query ROCK2 (ENSG00000134318) for **sQTL** in liver, adipose, muscle.
2. Note any strong sQTL (p<5e-8) and whether the spliced event hits the **kinase domain, coiled-coil, or Rho-binding domain**.
3. If a strong, domain-relevant sQTL exists → flag as a possible "activity-adjacent" instrument for a future targeted MR (decide then).
   If none → one "future direction" sentence.

Output: `task11_sqtl_feasibility.csv` (presence/absence + domain note). **Do not start a new MR axis pre-emptively (HARD STOP).**
