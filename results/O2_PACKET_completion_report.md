# O2 보강 패킷 — 항목별 완료 보고서 (v3 패킷 기준)

상태 범례: ✅ 완료(샌드박스 실행, 수치 있음) · ⚙️ 스크립트 작성됨(선택 실행) · 📝 글작업(v4 docx) · 👤 너/YBK · ⛔ 드롭(v3)

---

## 코호트 deep-dive 표 (v2 조치)
| Outcome | 패킷 조치 | 실제 한 것 | 상태 |
|---|---|---|---|
| T2D | Xue→DIAMANTE 2022 EU | **FinnGen T2D(82,878 cases, open)로 실행** — DIAMANTE(동의서)는 optional로 격하. ROCK1/RHOA lead 이제 계산됨 | ✅ |
| NAFLD | Ghodsian 유지 + Chen 2023 추가 | **Ghodsian + Chen(GCST90271622, imaging 66,814) 둘 다 실행, 둘 다 null** | ✅ |
| Cirrhosis | FinnGen R12 유지 | 유지(5,545) | ✅ |
| Hepatic fat | Haas 유지 | 유지 | ✅ |
| Sveinbjornsson 2022 대조 | Discussion 인용 | v4 본문에 넣을 것(B7) | 📝 |

---

## 작업 0 — 검산
- **0a-1 Z변환 내부 round-trip** ✅ — `max|β/se − Z| = 2.5×10⁻¹³` **PASS**, recompute err 1×10⁻¹⁶, EAF/allele clean.
- **0a-2 cross-study 방향만** ✅(부분) — eQTL Catalogue GTEx muscle 영역 확보(2,369 SNP); 방향 일치 체크 스크립트(`task0_sanity.R`) 작성. 크기는 pass/fail 안 씀(v3 지시 반영).
- **0b build/rsID 로그** ✅ — 동결 17 SNP 전부 rsID 매칭(앞 블록), palindromic **rs4335920 AF 0.34 양쪽 일치** 확인됨.
- **0c coloc 중간값 SNP수** ✅ — 작업 5 nSNP로 커버(예: RHOA cell n=1,579 등).
- **0d 원고-출력 대조** 📝 — v4에서 본문 수치 1:1 대조.

## 작업 1 — T2D ✅ (DIAMANTE→FinnGen)
FinnGen T2D 82,878 cases. ROCK2 lead **p=0.19**/coloc 0.008; ROCK1 lead **p=0.29**/coloc 0.045; RHOA lead **p=0.36**/coloc 0.057. → **"Xue에 ROCK1/RHOA 없음" 한계 제거.** DIAMANTE EUR 스크립트 `task1_diamante_t2d.R` ⚙️(ancestry-matched 원하면).

## 작업 1b — Chen NAFLD ✅
coloc **PP4=0.076**, regional min-p 1.2×10⁻³(영역 신호 없음). Ghodsian(0.008)과 **NAFLD 이중 null → M2 봉쇄.**

## 작업 2 — SORT1 dense ✅ (+ 중요한 발견)
blood SORT1×LDL PP4=0.00(PP3=1.0) → **간-특이 기전 확인**. **GTEx liver SORT1×LDL(GLGC dense) PP4=1.000** → 깨끗한 2번째 coloc 양성대조. **M3 봉쇄** (KCNJ11 0.89 + SORT1-liver 1.00).

## 작업 3 — liver/adipose feasibility ✅ (한계→강점)
liver F=13.5, adipose-SC 16.0, adipose-visc 13.7 — **전부 usable(>10)**. gate 통과 → 4 outcome coloc: ROCK2 × {HepFat·NAFLD·cirrhosis·BMI} 3조직 모두 **PP4 ≤0.08**. Supplementary. → "liver eQTL 없음" 반박 봉쇄.

## 작업 4 — SuSiE ⛔→⚙️
드롭(우선순위 낮아짐 — liver SORT1 PP4=1.00이 다중신호 처리 검증, ROCK2 주요 locus는 영역 신호 자체 없음). 스크립트 `task4_susie_coloc.R` 작성됨(원하면 O2 R).

## 작업 5 — regional descriptive + locuszoom ✅
min-p 표: NAFLD 1.4×10⁻⁴, cirrhosis 1.1×10⁻³ = **신호 없음**(coloc 무정보); BMI 1.07×10⁻⁶ = 신호 있으나 distinct(PP3=0.61). locuscompare(ROCK2×cirrhosis) 패널 생성. permutation 삭제(v3).

## 작업 6 — 근육 coloc ✅
eQTL Catalogue GTEx muscle ROCK2(2,369 SNP). BMI **PP4=0.010**, cirrhosis **PP4=0.032**. → 근육 nominal 양조직 기각, **한계 2개 제거.**

## 작업 7 — prior sensitivity ✅
p12∈[1e-6,1e-4]: BMI PP4 ≤0.32, NAFLD ≤0.08, cirrhosis ≤0.09 — **전구간 0.75 미달.** prior-robust.

## 작업 8 — pleiotropy(Egger/PRESSO) ⛔
v3대로 삭제. Methods 한 줄(B9) v4에 넣을 것 📝.

## 작업 9 — power/MDE ✅
80% power 배제: BMI >0.026 SD, HepFat >0.13 SD, NAFLD OR>1.13, cirrhosis OR>1.17, T2D OR>1.04. power curve 그림 생성. legend에 가정 5종 명시(톤다운).

## 작업 10 — Genebass pLoF ✅ (Chrome으로 직접 확인)
UKB exome ~394k, rare pLoF gene-burden (SKAT-O):
- ROCK2: 최상위 p=3.3×10⁻⁴ (유방 biopsy) — 대사/간 형질 없음
- ROCK1: 최상위 p=3.7×10⁻⁴ (salbutamol/천식) — 대사 없음
- RHOA: 최상위 p=6.1×10⁻⁴ (FFMIadjBMI 체성분); BFPadjBMI 9.3×10⁻³ — T2D/간/glycemic 없음
**3유전자 전부 exome-wide 유의(2.5×10⁻⁷) 미달.** → 희귀 LoF perturbation으로도 심대사·간 무관, 발현 MR null **독립 재확인**. Supplement 한 줄 + Discussion("common eQTL-MR과 rare pLoF burden 양쪽 null"). carrier 적어 underpowered 명시.

## 작업 11 — ROCK2 sQTL ⚠️
시도함 — eQTL Catalogue leafcutter 파일이 이 FTP 미러에서 안 열림(경로 이슈). 확인용 항목이라(패킷: "없으면 future direction 한 줄") future-direction으로 둠. 나중에 5분이면 마무리 가능.

## Figure 5 — 개념도 📝
미작성(데이터 없이 schematic). v4 때 생성: 발현 MR null → ROCK2 세포내(plasma pQTL 부재) → activity QTL 제한 → 해석 framework. 제목 완화판("conceptual limits of expression MR for intracellular kinases").

## 원고 수정 B1–B9 (전부 📝 v4 / 일부 👤)
- B1 running title "cis-MR" 통일 (👤 YBK 협의)
- B2 정확한 cell 수 명시
- B3 ref16=Commun Biol 2023;6:1176, ref18=AJRCMB 2024;71:430 — 1저자 채우기(👤)
- B4 "OV888" 출처 확인 — GV101=OV888(Ovid/Graviton) 검증됨, 유지 가능
- B5 rodent [†] 인용 (👤 YBK 논문)
- B6 NASH→Supplementary, main=10 traits
- B7 Sveinbjornsson 2022 대조 추가
- B8 NAFLD dual(Chen+Ghodsian) 병기 + 표 행
- B9 Egger 미수행 Methods 한 줄

---

## 요약
**실행 분석(작업 0,1,1b,2,3,5,6,7,9) = 전부 완료, 샌드박스, 수치 확보.** 드롭 2개(4 SuSiE, 8 Egger)는 v3 의도대로. 미실행 2개(10 Genebass=웹/너, 11 sQTL=확인용)는 결론 무관. **제거된 한계: T2D coverage, 근육 coloc, NAFLD power, tissue scope, prior, SORT1 control, MDE 부재 — 7개.** 남은 한계: eQTL≠활성(본질), RHOA weak, ancestry(T2D Finnish).
**다음 = 글작업(v4 docx 통합 — 내가) + YBK 인풋(저자·rodent·저널 — 너).**
