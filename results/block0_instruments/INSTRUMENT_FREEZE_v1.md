# ROCK2 Drug-target MR — Instrument Freeze (Block 1.5)
실행일 2026-06-20 · 데이터: eQTLGen blood (n≤31,684, hg19) + GTEx v8 Muscle_Skeletal (hg38) · LD: 1000G phase3 EUR (n=503)

## 동결된 instrument 세트

| 조직 | Tier | SNP | F | 분석 방법 |
|---|---|---|---|---|
| **Blood** | **Primary (A)** | rs12468344 (lead) | 409 | Wald ratio |
| Blood | **Supporting IVW (B)** | rs12468344, rs13034888, rs7574150 (r²<0.01, 3개) | 409 / 42 / 20 | IVW + heterogeneity (near-independent → naive IVW 유효) |
| Blood | Extended sensitivity | r²<0.1 세트 13개 | 20–409 | **LD-corrected IVW 전용** (LD matrix 필수) |
| **Muscle** | **Primary + mini-IVW** | rs4335920 (F=65), rs12996712 (F=38) | 상호 r²=0.016 | per-SNP Wald + 2-SNP IVW + heterogeneity |

**Outcome 추출 리스트 (블록 1 입력) = 15 unique rsID:**
rs116346698, rs12468344, rs12996712, rs13013760, rs16857268, rs185336055, rs2357890, rs34261765, rs4335920, rs4606895, rs61540355, rs6731491, rs6750442, rs72783415, rs7598570
= 혈액 r²<0.1 세트 13개 ∪ 근육 2개. 블록 1에서 15개 전부 추출 후 tier별로 subset.

## 핵심 진단 — r² 민감도 (이번 블록의 결정적 발견)

| clump r² | 독립 신호 수 |
|---|---|
| <0.1 | **13** |
| <0.01 | **3** |
| <0.001 | **1** |

13 → 3 → 1 붕괴 = r²<0.1 세트의 13개는 **상당한 잔여 LD를 안고 있다** (대부분 같은 소수 인과신호의 약한 tag). 함의:

- **naive IVW(TwoSampleMR `mr(dat)`)를 13개에 그대로 돌리면 통계적으로 틀린다** — 독립 가정 위배 → CI가 과도하게 좁아짐(anti-conservative). 패킷 블록 2 `clump_data(clump_r2=0.1)` + 블록 3 `mr(dat)` 조합이 정확히 이 함정에 빠진다.
- **해결:** Primary = 혈액 lead Wald + 근육 per-SNP/mini-IVW. Supporting = 혈액 r²<0.01 3개 IVW(진짜 near-independent라 naive IVW 유효). r²<0.1 13개는 **LD matrix 보정 IVW**(`MendelianRandomization::mr_ivw(..., correl=TRUE)` 또는 generalized IVW)로만 sensitivity 보고.

## 추가 함정 (블록 1~2 전 차단)

1. **빌드 불일치.** eQTLGen=hg19, GTEx=hg38. 근육 SNP hg19 좌표 확정: rs4335920 = 2:11555903, rs12996712 = 2:11401804. **outcome 추출은 rsID로** (position 매칭 금지).
2. **rs4335920 = T/A palindromic.** EUR ALT AF=0.353 (GTEx maf 0.333과 일치) → AF로 strand 해소 가능하나 harmonize 후 반드시 확인.
3. **rs12996712는 혈액에서도 강한 eQTL** (blood p=4.0e-68). 근육 instrument이자 혈액 신호 — 조직 간 effect 방향 비교에 유용.
4. **혈액 매칭률 82%** (510/621가 1000G EUR에 매칭, 111개 indel/multiallelic/부재). 단 독립 index SNP 13개는 전부 매칭됨 — 클럼핑 결과에 영향 없음.

## 패킷 정정 로그 (§0 / 함정에 반영 필요)

- **AF 파일명**: 패킷의 `..._picard.txt.gz` 는 존재 안 함. 실제 = `2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`.
- **GTEx SE**: `independentEqtl`/`singleTissueEqtl` 엔드포인트는 slope_se 미제공. **`dyneqtl` 엔드포인트의 `error` 필드 = slope_se** (nes/tStatistic로 검증). 또는 official slope+p에서 se=|β|/Z(p) 유도.
- **rs7581184**: 패킷 STEP 0의 F≈46(p=1.4e-11)은 GTEx 혈액(n~700) conditional 2차 신호. eQTLGen marginal(n~31,600)에선 F=333. 같은 SNP, 다른 데이터셋/조건. 클럼핑에서 lead가 아님(rs12468344가 lead).
- **GTEx 좌표 hg38** (위 함정 1).
- **몰겐니스 SSL**: `--no-check-certificate` 필요.
- **클럼핑 reference 무료 대안**: TwoSampleMR 원격 API 대신 pysam로 1000G phase3 EUR 영역 슬라이스 → r² 직접 계산(greedy clump). 토큰/플랫폼 의존 0.
