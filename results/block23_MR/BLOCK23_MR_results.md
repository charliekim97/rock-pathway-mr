# Block 2/3 — ROCK2 drug-target MR + positive controls

효과 단위: **per-SD ROCK2 발현** → outcome SD (T2D = log-OR). 노출: eQTLGen 혈액(Z→β/SE) + GTEx 근육(nes, SE=|β|/Z(p)). 결과 SE: API outcome=|β|/Z(p), flat-file(HepFat·BMI)=exact. LD-보정 IVW = 1000G EUR signed LD matrix(Burgess correlated-IVW).

## ROCK2 — 추정치 (β [95% CI] p)

| Outcome | 혈액 lead Wald | 근육 2-SNP IVW (Q p) | 혈액 r²<.01 IVW | 혈액 r²<.1 LDcorr (sens) | tissue concord. |
|---|---|---|---|---|---|
| **HepaticFat** | −0.029 [−0.124,0.065] 0.54 | −0.030 [−0.098,0.038] 0.39 (Qp .54) | −0.008 0.86 | +0.004 0.93 | concordant null |
| **FI** (IR) | −0.016 [−0.037,0.005] 0.13 | −0.005 [−0.020,0.011] 0.53 (Qp .26) | −0.011 0.28 | −0.009 0.33 | concordant null |
| **FG** | −0.009 [−0.027,0.009] 0.34 | (lead 결손) | −0.011 0.18 | −0.004 0.51 | — |
| **HbA1c** | +0.008 [−0.004,0.020] 0.19 | +0.004 0.50 (Qp .35) | +0.007 0.17 | +0.012 [0.001,0.023] **0.036** | concordant |
| **2hGlu** | +0.010 0.90 | +0.025 0.81 (Qp .89) | −0.007 0.92 | +0.018 0.80 | concordant null |
| **T2D** | −0.013 [−0.111,0.085] 0.79 | −0.031 [−0.100,0.038] 0.38 (Qp .43) | (1 SNP만) | +0.005 0.92 | concordant null |
| **BMI** | **+0.037 [0.014,0.060] 0.0019** | −0.000 [−0.017,0.017] 0.99 (Qp **.02**) | **+0.039 [0.017,0.060] 4.2e-4** | +0.025 [0.006,0.045] 0.012 | **DISCORDANT** |

Steiger: 모든 outcome에서 r²_exp(0.013) ≫ r²_out(~1e-5) → 방향 exp→out 확정(역인과 배제). 단 outcome 신호가 ~0이라 trivial.

### 읽을 거리 (해석은 다음 턴)
- **6/7 outcome 전방위 null.** HepFat·FI·FG·2hGlu·T2D 전부 lead·근육·IVW 일치 null. HbA1c는 LD-corr 13-SNP만 nominal(p=0.036) — lead/근육 null이라 robust 아님, 다중검정(≈28테스트, Bonferroni~0.0018) 통과 못함.
- **BMI만 혈액에서 신호:** lead Wald p=0.0019, r²<.01 IVW p=4.2e-4 (Bonferroni 통과). **그러나 근육은 완전 null(p=0.99)이고 tissue DISCORDANT**, 근육 2-SNP 간 이질성(Q p=0.02: rs4335920 −0.015 vs rs12996712 +0.026). → 조직-일치하는 깨끗한 ROCK2-발현→BMI 인과로 보기 어려움. 혈액-특이 신호인지, cis 영역의 다른 유전자/변이 pleiotropy인지 **coloc(블록 4)로 판별 필요.**

## Positive controls — 같은 파이프라인(eQTLGen 혈액 cis-eQTL → outcome, lead Wald)

| Control | eQTL lead (F) | 추정치 | p |
|---|---|---|---|
| **KCNJ11 → T2D** (Xue) | rs2074310 (F=166) | log-OR **+0.409** [0.310,0.508] | **6.7e-16** |
| **SORT1 → LDL** (Willer2013) | rs657420 (F=105) | βSD **−0.766** [−0.858,−0.673] | **7e-60** (SNP) |
| GIPR → 2hGlu/FG | rs12972094 (F=57) | −0.218 / −0.032 | 0.12 / 0.20 (ambiguous*) |

\* GIPR lead 혈액 eQTL이 glycemic-causal 변이와 정렬 안 됨(GIPR locus는 변이별 방향 복잡). KCNJ11·SORT1 두 강신호로 검증 충분.

**결론(방법 검증):** 동일 파이프라인이 진짜 신호를 강하게 잡는다(KCNJ11→T2D 6.7e-16, SORT1→LDL 7e-60, 둘 다 알려진 방향). → **ROCK2 전방위 null은 방법 아티팩트가 아니라 진짜 무신호.** PCSK9는 혈액 미발현(eQTLGen 부재)이라 SORT1로 대체.

## 소스
- 노출: eQTLGen blood cis-eQTL; GTEx v8 Muscle_Skeletal.
- Outcome: HepFat GCST90029073; FI GCST90002238; FG GCST90002232; HbA1c GCST90002244; 2hGlu GCST90002227; T2D GCST006867(Xue2018); BMI GCST009004(Pulit2019); LDL GCST002222(Willer2013).
