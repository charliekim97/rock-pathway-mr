# Fibrosis/cirrhosis 축 — ROCK2 (+ROCK1/RHOA), MR + coloc

YBK 타겟 직격 축. 동결 ROCK2 17 instrument + ROCK1/RHOA(cirrhosis만). FinnGen=hg38, NAFLD meta=GRCh37, rsID 매칭.

## Outcome (case 수)
| Outcome | 소스 | cases | 검정력 |
|---|---|---|---|
| Cirrhosis (broad) | FinnGen R12 CIRRHOSIS_BROAD | 5,545 | 양호 |
| NAFLD meta | Ghodsian 2021, GCST90091033 | 8,434 | 최고 |
| Fibrosis+cirrhosis K74 | FinnGen K11_FIBROCHIRLIV | 2,653 | 보통 |
| NASH | FinnGen NASH | 254 | ⚠ 매우 약함 |

ROCK2 coverage: lead 1/1, r²<.01 3/3, r²<.1 12/13, 근육 2/2 (모든 outcome 공통, ~16/17).

## MR (per-SD ROCK2 발현 → log-OR)

| Outcome | 혈액 lead Wald | 근육 2-SNP IVW | 혈액 r²<.01 IVW |
|---|---|---|---|
| **Cirrhosis** | −0.099 p=0.41 | **−0.182 [−0.351,−0.013] p=0.035** | −0.132 p=0.24 |
| NAFLD meta | −0.091 p=0.37 | −0.031 p=0.67 | −0.109 p=0.25 |
| FibrChir K74 | +0.053 p=0.76 | −0.009 p=0.94 | +0.026 p=0.87 |
| NASH | +1.123 p=0.046† | +0.192 p=0.63 | +0.676 p=0.19 |

† NASH 254 cases — log-OR +1.12(OR≈3)에 CI [0.02,2.22], 사실상 노이즈. 검정력 한계로 해석 불가(명시).

## Coloc (PP4)

| eQTL × outcome | PP3 | PP4 | n |
|---|---|---|---|
| **ROCK2(blood) × Cirrhosis** | 0.301 | **0.010** | 8199 |
| ROCK2 × NAFLD meta | 0.186 | 0.008 | 6426 |
| ROCK2 × FibrChir K74 | 0.294 | 0.010 | 8199 |
| ROCK2 × NASH | 0.441 | 0.093 | 8199 |
| ROCK1 × Cirrhosis | 0.095 | 0.055 | 1662 |
| RHOA × Cirrhosis | 0.069 | 0.023 | 3182 |

ROCK1/RHOA × Cirrhosis MR(lead Wald): ROCK1 p=0.74, RHOA p=0.45 (둘 다 null).

## 판정 — fibrosis 축도 NULL (프레임 유지)

**핵심 질문(cirrhosis × ROCK2 coloc PP4) 답: 0.010 = no-coloc.** 혈액 MR도 null(p=0.41). NAFLD 메타(최고 검정력) 전부 null + no-coloc. K74도 null. → **fibrosis 축에서 ROCK2 발현↔질환 인과 없음.** pathway-level null이 가장 thesis-그럴듯한 축에서도 확인됨.

### 단 하나의 flag: ROCK2 **근육** → Cirrhosis IVW p=0.035 (protective 방향)
- 근육-특이 nominal. **혈액은 MR·coloc 모두 null** → tissue-discordant, 비-robust.
- 다중검정(이제 11개+ outcome) 통과 못함.
- 방향이 protective(발현↑→cirrhosis↓) = thesis 가설(ROCK2 억제가 항섬유화)과 **반대** — 억제가 오히려 해로울 수 있다는 쪽. 단 비-robust라 과대해석 금지.
- **근육 eQTL coloc은 GTEx 근육 all-variant 필요(현 미보유) → O2에서 확인사살 예정.** 혈액 no-coloc + 다중검정 + 단일조직으로 이미 정황상 기각.

## 논문 프레임 결론
fibrosis 축 추가가 **null 서사를 강화**한다 — "가장 그럴듯한 fibrosis/NAFLD/cirrhosis에서도, 최고 검정력 메타에서도 신호 없음". Discussion에 fibrosis 축을 별도 문단으로, 근육-cirrhosis nominal은 limitation/future(O2 muscle coloc)로 처리. 프레임은 pure pathway-level null 유지.

## 미완 (O2)
- ROCK2 근육 × cirrhosis coloc (GTEx 근육 all-variant) — 근육 nominal 확인사살.
- (기존) BMI × 근육, SORT1 dense LDL.

## 소스
FinnGen R12 (CIRRHOSIS_BROAD, K11_FIBROCHIRLIV, NASH); Ghodsian 2021 NAFLD meta GCST90091033; eQTLGen blood / GTEx muscle eQTL.
