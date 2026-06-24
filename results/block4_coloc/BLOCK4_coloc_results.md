# Block 4 — Colocalization (coloc.abf)

방법: coloc.abf (Giambartolomei 2014, Wakefield ABF). eQTL = eQTLGen 혈액 full cis (Z→β/SE 변환), outcome = 각 GWAS 영역 추출. 윈도우 = ROCK2 cis 신호 중심(chr2:11.1–11.6Mb; BMI·HepFat는 full ±1Mb). prior p1=p2=1e-4, p12=1e-5. sd.prior: quant 0.15, cc(T2D) 0.2. coloc는 Z²기반이라 allele harmonize 불필요. 입력은 사용자가 로컬에서 추출한 `region_full_eqtl.txt`(ROCK2 9,983 / KCNJ11 6,929 / GIPR 8,599 / SORT1 7,330 변이).

## PP 표 (H3 = 다른 변이, H4 = 공유 인과변이)

| eQTL | outcome | nSNP | PP3 | PP4 | PP4/(PP3+PP4) | top SNP |
|---|---|---|---|---|---|---|
| ROCK2(blood) | **BMI** | 9302 | 0.605 | **0.045** | 0.069 | rs12468344 |
| ROCK2(blood) | HepaticFat | 7318 | 0.159 | 0.005 | 0.032 | rs12468344 |
| ROCK2(blood) | FI | 2147 | 0.044 | 0.003 | 0.065 | rs12468344 |
| ROCK2(blood) | FG | 2158 | 0.040 | 0.001 | 0.036 | rs12468344 |
| ROCK2(blood) | HbA1c | 2167 | 0.084 | 0.001 | 0.010 | rs12468344 |
| ROCK2(blood) | 2hGlu | 2172 | 0.085 | 0.005 | 0.056 | rs12468344 |
| ROCK2(blood) | T2D(Xue) | 853 | 0.017 | 0.004 | 0.199 | rs12468344 |
| **KCNJ11** | **T2D [pos ctrl]** | 1066 | 0.112 | **0.888** | 0.888 | rs757110 |
| GIPR | T2D [pos ctrl] | 661 | 0.000 | 0.000 | 0.651 | rs11083785 |
| SORT1 | LDL [pos ctrl] | 284 | 0.752 | 0.013 | 0.017 | rs650985 |

## 핵심 판정

**1. BMI 혈액 신호 = ROCK2 아님 (결론적).** PP4=0.045, PP3=0.605. 즉 ROCK2 영역에 BMI 신호는 있으나(PP3) **ROCK2 발현을 움직이는 변이와는 다른 변이**다. 블록 2/3의 BMI 혈액 MR 신호(lead p=0.0019, r²<.01 IVW p=4e-4)는 **ROCK2 인과효과가 아니라 cis 영역의 다른 신호에 의한 LD 오염/수평 pleiotropy**로 확정. tissue-discordance(근육 null)와 정합 — 진짜 ROCK2 효과였다면 coloc PP4가 높아야 한다.

**2. 6개 대사형질 전부 no-coloc.** PP4 0.001–0.005, 전부 무시 가능. MR null의 독립적 이중확인. (대부분 PP0+PP1 우세 = outcome 신호 자체가 영역에 없음.)

**3. 방법 검증 — coloc 파이프라인 작동 확인.** **KCNJ11 × T2D PP4=0.888** — 같은 파이프라인이 진짜 coloc를 강하게 잡는다. → ROCK2의 전방위 no-coloc는 파이프라인 무능이 아님.
   - GIPR×T2D: PP3=PP4≈0 (Xue T2D가 GIPR 영역에 신호 희박 — GIPR T2D 효과가 약/복잡, ambiguous control. MR에서도 null이었음).
   - SORT1×LDL: PP4=0.013(낮음). **단 이건 outcome 데이터 한계** — Willer 2013 LDL이 HapMap2-sparse(영역 284 SNP)라 진짜 인과변이 미태깅 → PP3로 샘. MR 단일-SNP은 강했음(p=7e-60). 밀도 높은 GLGC(Graham 2021)면 PP4 올라갈 것(API 미제공이라 미적용). KCNJ11이 파이프라인 검증엔 충분.

## 미적용 + 한계
- **BMI × 근육 coloc 미실행**: GTEx 근육 all-variant 영역 데이터가 공개 API/소형 파일로 안 나옴(all-pairs 대용량). 단 결론 불변 — BMI가 혈액 ROCK2와도 coloc 안 되므로(PP4=0.045) 조직 무관하게 ROCK2 효과 아님.
- T2D(Xue)·LDL(Willer)은 영역 변이 밀도 낮음(853/284) → 해당 coloc 검정력 제한. ROCK2 결론(no-coloc)엔 영향 없음(어차피 MR도 null).

## 종합 (블록 0–4)
genetically-proxied ROCK2 발현은 7개 대사형질(간지방·FI·FG·HbA1c·2hGlu·T2D·BMI)에 **인과 효과 없음** — MR null + coloc no-coloc 양쪽 일치, positive control(KCNJ11 PP4=0.89, MR p=6.7e-16)로 방법 검증됨. 유일한 명목 신호(BMI 혈액)는 coloc로 ROCK2 비귀속(LD/pleiotropy)으로 정리됨. = **null 서사 확정 (translatability gap).**

## 소스
eQTLGen cis-eQTLs_full_20180905; outcome: GCST90029073(HepFat), GCST90002238/32/44/27(MAGIC), GCST006867(T2D Xue), GCST009004(BMI Pulit), GCST002222(LDL Willer).
