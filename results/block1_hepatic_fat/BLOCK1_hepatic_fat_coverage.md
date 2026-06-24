# Block 1 — Hepatic fat outcome coverage

**Outcome GWAS:** Haas et al. 2021, *Cell Genomics* — UKB MRI-derived liver fat (ML-quantified). GWAS Catalog **GCST90029073**, N=32,974 European, build GRCh37, harmonised full sumstats. PMID 34957434.

## Coverage — 17/17 instrument SNP 전부 존재, proxy 불필요

| Tier | SNP | outcome coverage |
|---|---|---|
| Blood lead (Primary, Wald) | rs12468344 | **1/1** |
| Blood r²<0.01 (Supporting IVW) | rs12468344, rs13034888, rs7574150 | **3/3** |
| Blood r²<0.1 (LD-corr sensitivity) | 13개 | **13/13** |
| Muscle (2-SNP IVW) | rs4335920, rs12996712 | **2/2** |

동결 리스트 17개로 갱신(혈액 13 ∪ r²<0.01 3 ∪ 근육 2). r²<0.01 supporting 세트의 rs13034888·rs7574150는 r²<0.1 13-세트에 없어(greedy clump이 임계값마다 다른 index SNP 선택) 추가 추출함.

## Palindromic 체크

- **rs4335920 (T/A, muscle primary) — 해소됨.** exposure EA=A freq 0.353(1000G EUR)/0.333(GTEx), outcome EA=A freq 0.341. 같은 allele, AF 일치, 0.5에서 멀음 → strand 명확. **근육 fallback(rs12996712 단독) 불필요.**
- rs7598570 (C/G, sensitivity tier만) — exposure/outcome effect allele 표기 반대 + palindromic. r²<0.1 보조 tier 전용이고 outcome p=0.55(null)이라 영향 미미하나, harmonize 시 AF 기반 정렬 명시 필요.
- 나머지 비-palindromic allele flip(rs6731491, rs4606895 등)은 harmonize 자동 처리.

## ⚠ 단일-SNP preview — 간지방은 전 instrument NULL

| SNP | tier | outcome beta | se | p |
|---|---|---|---|---|
| rs12468344 (blood lead) | Primary | −0.0047 | 0.0077 | **0.71** |
| rs4335920 (muscle) | Muscle | +0.0086 | 0.0081 | **0.32** |
| rs12996712 (muscle) | Muscle | −0.0001 | 0.0078 | **0.79** |
| rs34261765 (최소 p) | sens. | +0.022 | 0.012 | 0.081 |

17개 중 최소 p=0.081. 어떤 instrument도 간지방과 연관 없음. **정식 MR(블록 2/3)은 hepatic fat에서 거의 확실히 null이 나온다.** 패킷이 예상한 null 분기 — 단, 가장 thesis-그럴듯한 outcome에서의 null이라 정보가 높은 negative(translatability 갭 서사). 정식 추정·CI는 블록 2/3에서 확정.
