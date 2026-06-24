# 완성된 Coloc 매트릭스 — 3 유전자 × 7 형질 (PP4)

coloc.abf, eQTLGen 혈액 cis-eQTL × 각 형질 GWAS 영역. PP4 = 공유 인과변이 사후확률. **colocalization 통상 임계 PP4≥0.75.**

| 유전자(혈액 eQTL) | HepFat | FI | FG | HbA1c | 2hGlu | T2D | BMI |
|---|---|---|---|---|---|---|---|
| **ROCK2** (F=409) | 0.005 | 0.003 | 0.001 | 0.001 | 0.005 | 0.004 | 0.045 |
| **ROCK1** (F=86) | 0.070 | 0.018 | 0.013 | 0.007 | 0.012 | 0.206 | 0.004 |
| **RHOA** (F=27, 1 SNP) | 0.020 | 0.005 | 0.399* | 0.012 | 0.061 | 0.249* | 0.309* |

**전 셀 PP4 < 0.40 — 어떤 유전자×형질도 colocalization 임계(0.75) 근처도 못 감. 경로 전체 × 표현형 전체에서 공유 인과변이 0건.**

## 해석

1. **ROCK2·ROCK1 = 깨끗한 null.** 강한 eQTL(F=409, 86) → coloc 검정력 충분 → 전 형질 PP4≤0.07. 신뢰도 높은 null.

2. **RHOA(*) 중간값은 신호 아니라 weak-instrument 아티팩트.** RHOA는 혈액 FDR-유의 cis-eQTL 1개(F=27)뿐 → eQTL 신호 약함 → coloc가 H0를 못 좁혀 사후확률이 H3/H4로 분산 → FG 0.399, BMI 0.309, T2D 0.249처럼 "따뜻해" 보임. 그러나 (a) 전부 임계 0.75 미달, (b) RHOA MR 전부 null(FG p=0.14, BMI 0.51), (c) PP4/(PP3+PP4)도 0.5 안팎. → **RHOA coloc는 underpowered, 인과 근거 아님.** 논문엔 weak-instrument 한계로 명시하고 ROCK1/ROCK2를 robust 결론의 축으로.

3. **ROCK1×T2D 0.206**: 중간값이나 임계 미달 + ROCK1 T2D MR은 lead가 Xue 부재로 NA. 영역 coloc도 비유의. 신호 아님.

## 명목 MR hit 3개 — 전부 coloc 기각 (Fig4 완성)
| 명목 hit | MR p | coloc PP4 | 판정 |
|---|---|---|---|
| ROCK2 × BMI | 4e-4 (IVW) | 0.045 | 기각(다른 변이) |
| ROCK1 × HepaticFat | 0.022 (IVW) | 0.070 | 기각 |
| RHOA × 2hGlu | 0.048 (Wald) | 0.061 | 기각 |

세 명목 연관 전부 colocalization 없음 → LD/pleiotropy. 인과 아님.

## 결론 (분석 종료)
Rho–ROCK 축(RHOA·ROCK1·ROCK2) × 7 심대사형질 = **MR null + coloc no-coloc**, positive control(KCNJ11 PP4=0.89)로 방법 검증. **pathway-level null 확정.** 분석 블록 0–5 + coloc 매트릭스 완료. 다음 = figure 생성 → draft.
