# Block 1 (batch) — Outcome × Instrument coverage matrix

동결 instrument 17개(혈액 13 r²<0.1 ∪ r²<0.01 3 ∪ 근육 2)를 6개 metabolic outcome에서 추출. 전부 EUR, rsID 매칭.

## Outcome 소스 (확정)

| Outcome | 소스 | N | accession | 경로 | 비고 |
|---|---|---|---|---|---|
| Hepatic fat | Haas 2021 Cell Genom | 32,974 | GCST90029073 | EBI flat (GRCh37) | thesis 직결 |
| Fasting insulin (FI, BMI-adj) | MAGIC Chen 2021 | 151,013 | GCST90002238 | GWAS Catalog API | **IR 주력 readout** |
| Fasting glucose (FG) | MAGIC Chen 2021 | 200,622 | GCST90002232 | API | |
| HbA1c | MAGIC Chen 2021 | 146,806 | GCST90002244 | API | |
| 2h-glucose | MAGIC Chen 2021 | 63,396 | GCST90002227 | API | bonus |
| T2D | **Xue 2018** Nat Commun | 62k ca/596k co (EUR) | GCST006867 | API | Mahajan 대체(아래) |
| BMI | **Pulit 2019** GIANT+UKB | 806,834 | GCST009004 | EBI flat | Yengo 대체(아래) |
| HOMA-IR | Dupuis 2010 | ≤46,186 | — | MAGIC flat | 너무 sparse, 폐기 권고 |

## Coverage 매트릭스 (tier별)

| Tier | HepFat | FI | FG | HbA1c | 2hGlu | T2D | BMI | HOMA-IR |
|---|---|---|---|---|---|---|---|---|
| 혈액 lead (Primary Wald) | 1/1 | 1/1 | 1/1 | 1/1 | 1/1 | **1/1** | 1/1 | 0/1 |
| 혈액 r²<.01 (Support IVW) | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 | 1/3 | 3/3 | 0/3 |
| 혈액 r²<.1 (LD-corr sens) | 13/13 | 13/13 | 13/13 | 13/13 | 13/13 | 4/13 | 13/13 | 3/13 |
| 근육 2-SNP | 2/2 | 2/2 | 2/2 | 2/2 | 2/2 | **2/2** | 2/2 | 0/2 |
| **TOTAL /17** | **17** | **17** | **17** | **17** | **17** | **6** | **17** | **3** |

→ **6개 핵심 outcome(HepFat·FI·FG·HbA1c·2hGlu·BMI)는 primary+supporting+sensitivity 전부 완전.** proxy 불필요.
→ T2D는 **primary tier(lead+근육 2개) 완전**, supporting/sensitivity만 결손.
→ HOMA-IR 폐기 — FI(151k, 17/17)가 동일 IR 축을 6× N으로 대체.

## Palindromic — 둘 다 전 outcome에서 깨끗하게 해소

- **rs4335920 (T/A, 근육 primary):** EA=A, freq 0.34–0.35가 모든 outcome에서 일치(exposure 0.333/0.353). 0.5에서 멀어 strand 무모호. 근육 primary 그대로 사용.
- **rs7598570 (C/G, sensitivity tier):** 모든 outcome out EA=C freq~0.66 = exposure G freq 0.31의 반대 allele. G-freq ~0.34로 정합. 해소됨.

## 소스 대체 — 판정 필요 3건

1. **T2D: Mahajan 2022 → Xue 2018.** Mahajan 2018/2022 둘 다 GWAS Catalog harmonised/API 미제공(DIAMANTE 데이터접근 동의서 필요, 스크립트 불가). Xue 2018(EUR, dense)은 API 제공이고 **primary tier(lead+근육) 완전**. supporting/sensitivity 결손은 secondary. T2D에 신호 뜨면 그때 Mahajan 2018 EUR을 동의서로 받아 full-tier 보강 가능.
2. **IR: HOMA-IR → FI.** HOMA-IR(Dupuis 2010)은 HapMap-era라 lead·근육 0개. FI(Chen 2021)가 우월. HOMA-IR 폐기.
3. **BMI: Yengo 2018 → Pulit 2019.** Yengo는 HapMap2 sparse(3/17, lead 결손). Pulit(GIANT+UKB, 807k, 1000G-dense)로 17/17.

## SE 처리 메모
API outcome(FI/FG/HbA1c/2hGlu/T2D)은 SE 미제공 → SE=|β|/Z(p)로 유도(정규근사, 표준). flat-file outcome(HepFat·BMI)은 exact SE 보유. block 3 primary 추정에서 특정 SNP이 결정적이면 해당 flat file로 exact SE 정제 가능(현재 전부 null이라 실익 없음).

## ⚠ 신호 preview — 6 outcome 전부 NULL (단일-SNP)
| SNP | FI | FG | HbA1c | 2hGlu | T2D | HepFat |
|---|---|---|---|---|---|---|
| rs12468344 (lead) | 0.13 | 0.34 | 0.19 | 0.90 | 0.79 | 0.71 |
| rs4335920 (근육) | 0.66 | 1.0 | 0.30 | 0.98 | 0.24 | 0.32 |
| rs12996712 (근육) | 0.22 | 0.12 | 0.60 | 0.78 | 0.90 | 0.79 |

어떤 instrument도, 어떤 outcome에서도 연관 없음(전 p>0.12). **정식 MR은 전방위 null로 갈 가능성 매우 높다.** = 패킷의 null 분기(translatability 갭). 정식 추정·CI·tissue concordance는 블록 2/3.
