# Block 5 — Pathway-wide (RhoA / ROCK1), 동일 파이프라인

목적: Rho–ROCK 축 상류(RHOA)·동류(ROCK1)도 사람 대사형질에 인과 없으면 → **pathway-level null로 격상**.

노출: eQTLGen 혈액 cis-eQTL (Z→β/SE). 효과 단위 = per-SD 발현. coloc = full 영역(사용자 로컬 추출 `rho_region_full_eqtl.txt`: RHOA 4,335 / ROCK1 2,230 cis SNP).

## Instrument
- **RHOA**: eQTLGen 혈액 FDR-유의 cis-eQTL **단 1개** (rs115431681, F=27). → 단일-SNP Wald만. (RHOA는 혈액에서 약한 eQTL 유전자.)
- **ROCK1**: 229개 유의, r²<.01 독립 **2개** (lead rs142716063 F=86, rs72881399). → lead Wald + 2-SNP IVW.

## MR (β [95%CI] p)

| Outcome | RHOA Wald(1SNP) | ROCK1 lead Wald | ROCK1 r²<.01 IVW |
|---|---|---|---|
| HepaticFat | −0.092 [−0.23,0.05] 0.19 | −0.131 [−0.26,−0.003] **0.045** | −0.143 [−0.27,−0.02] **0.022** |
| FI | −0.014 0.49 | +0.038 0.28 | +0.009 0.74 |
| FG | +0.040 0.14 | +0.008 0.75 | −0.015 0.45 |
| HbA1c | +0.008 0.93 | −0.022 0.31 | −0.032 0.094 |
| 2hGlu | +0.211 [0.002,0.42] **0.048** | +0.013 0.84 | +0.007 0.92 |
| T2D | NA | NA | NA |
| BMI | +0.012 0.51 | −0.005 0.77 | −0.005 0.75 |

T2D = NA: RHOA·ROCK1 instrument SNP이 Xue2018 T2D sparse 커버리지에 부재(proxy 필요, 미적용). 명목 hit 2개(ROCK1→HepFat, RHOA→2hGlu)는 다중검정(≈18테스트, Bonferroni~0.0028) 통과 못함.

## Coloc (PP4 = 공유 인과변이)

| eQTL | outcome | nSNP | PP3 | PP4 |
|---|---|---|---|---|
| **ROCK1** | **HepaticFat** | 1387 | 0.033 | **0.070** |
| ROCK1 | BMI | 2059 | 0.046 | 0.004 |
| RHOA | HepaticFat | 2529 | 0.183 | 0.020 |
| RHOA | BMI | 3931 | 0.382 | 0.309* |

## 판정

**1. ROCK1→HepaticFat 명목신호 = 가짜 (확인사살).** MR nominal(IVW p=0.022)였으나 **coloc PP4=0.070** → ROCK1 발현 변이와 간지방 신호는 다른/무신호. thesis(ROCK1-NAFLD)의 가장 그럴듯한 지점인데도 인과 colocalization 없음. ROCK2 BMI 때와 동일 패턴 — MR 명목 → coloc로 기각.

**2. RHOA×BMI PP4=0.309(*)는 신호 아님.** MR null(p=0.51)인데 PP4 중간값 — **RHOA가 혈액 약한 eQTL(유의 SNP 1개)이라 coloc 검정력 낮아 PP가 가설 간 분산**된 것. PP4<0.7 + MR null → 인과 근거 아님. weak-instrument 한계로 명시.

**3. 나머지 전부 null.** MR·coloc 양쪽.

## 종합 — pathway-level null 격상
Rho–ROCK 축 3개 노드(상류 RHOA / ROCK1 / ROCK2) **모두** genetically-proxied 발현이 사람 대사형질에 인과 없음. 단일 타겟(ROCK2) null → **경로-수준 null**로 격상됨. single-target보다 publishable(경로 전체의 translatability gap).

## 한계 / 미완 (확인사살은 끝, 보강 선택)
- 글리세믹(FI/FG/HbA1c/2hGlu)·T2D × RHOA/ROCK1 coloc 미실행(MR null이라 confirmatory). 원하면 API 영역추출로 채움.
- T2D MR: RHOA/ROCK1 lead가 Xue sparse에 부재 → proxy(r²>0.8) 또는 dense T2D(Mahajan) 필요.
- RHOA 단일-SNP/약한 eQTL → RHOA 결론은 ROCK1/ROCK2보다 약함(명시).

## 소스
eQTLGen cis-eQTLs_full (RHOA ENSG00000067560 / ROCK1 ENSG00000067900); HepFat GCST90029073; BMI GCST009004.
