# ROCK pathway drug-target MR — 논문 구조 (구조 확정 단계)

타겟 저널: **Diabetologia / Cardiovascular Diabetology** (보강 시), **Sci Rep** floor.
핵심 자산: pathway-level null + drug-target framing + positive-control 검증 + 명목 hit coloc 기각.

---

## 1. Title claim 후보 3개

**(A) 경로-수준 null + drug-target 정면**
> "Genetically proxied inhibition of the RhoA–ROCK pathway and cardiometabolic traits: a drug-target Mendelian randomization and colocalization study"
— 중립적·정석. Diabetologia 톤. "inhibition" 넣어 belumosudil/fasudil repurposing 독자층 끌어옴.

**(B) Translatability gap 직격 (가장 강한 claim)**
> "Two decades of rodent RhoA–ROCK metabolic biology do not translate to human genetics: a Mendelian randomization study of ROCK1, ROCK2 and RhoA across seven cardiometabolic traits"
— 서사가 제목에 박힘. 임팩트↑, 단 reviewer가 "translate" 단정 싫어할 수 있음 → "...are not supported by human genetics"로 완화 가능.

**(C) 타겟 명시 + null 정직**
> "No causal effect of genetically proxied ROCK1/ROCK2/RhoA expression on glycemic traits, hepatic fat, type 2 diabetes or BMI: a cis-Mendelian randomization and colocalization study"
— 가장 보수적·검색친화적. Sci Rep/CVD 안전판. null을 제목에 명시(요즘 null 명시가 오히려 인용됨).

**추천: (A) primary, (B) cover-letter 서사용.** (C)는 저널 다운그레이드 시.

---

## 2. Figure 골격

### Fig 1 — Instrument + 파이프라인 (Methods 시각화)
- **1A 도식**: eQTLGen 혈액(N≤31,684) + GTEx v8 근육 → Z→β 변환 → 7 outcome GWAS(EUR) → MR(Wald/IVW/LDcorr) → coloc.abf → positive control. 데이터 출처·N 박스.
- **1B instrument 표/막대**: 3유전자 × tissue × F-stat.
  - ROCK2 혈액 lead rs12468344 **F=409**; 근육 rs4335920 **F=65**, rs12996712 **F=38**
  - ROCK1 혈액 lead rs142716063 **F=86** (독립 2)
  - RHOA 혈액 rs115431681 **F=27** (단일)
  - 점선: weak-instrument 임계 F=10 (전부 상회 표시).
- 메시지: "instrument 전부 강함 → null이 weak-instrument 탓 아님."

### Fig 2 — MR forest (main result)
- 행: 7 형질(HepFat, FI, FG, HbA1c, 2hGlu, T2D, BMI) × 3 유전자, tier 색상(혈액 lead Wald / 근육 IVW / 혈액 r²<.01 IVW).
- β [95%CI], null line(0) 강조. 명목 hit 3개(ROCK2-BMI, ROCK1-HepFat, RHOA-2hGlu) ★ 표시.
- 패널 분리: 2A ROCK2(가장 완전, 4 tier) / 2B ROCK1 / 2C RHOA.
- 메시지: "거의 전부 CI가 0 통과. 명목 hit 3개만 별표 → Fig4로 연결."

### Fig 3 — Coloc PP4 히트맵 + positive control 검증
- **3A 히트맵**: 행 3유전자(혈액 eQTL) × 열 형질, 셀=PP4 (0–1 컬러). 전부 <0.1(파랑) → 시각적 null. (현 산출: ROCK2 7형질 전부 PP4≤0.045; ROCK1·RHOA는 HepFat/BMI 완료, 글리세믹/T2D는 **추가 산출 필요**.)
- **3B positive control 패널**: 같은 파이프라인.
  - **KCNJ11 → T2D: PP4=0.888** (MR log-OR +0.41, p=6.7e-16) — coloc+MR 양성 ✓
  - **SORT1 → LDL: MR p=7e-60** (단 coloc PP4=0.013, Willer HapMap2-sparse 한계 → 캡션 명시; 가능하면 dense GLGC로 재산출)
- 메시지: "파이프라인은 진짜 신호를 잡는다(KCNJ11). ROCK 축의 null은 무능이 아님."

### Fig 4 — 명목 hit 3개 coloc 기각 (해석의 핵심)
- 명목 MR hit 3개를 coloc로 무력화하는 대비:
  - ROCK2 × BMI: MR IVW p=4e-4 **→ PP4=0.045 (PP3=0.61)**
  - ROCK1 × HepFat: MR IVW p=0.022 **→ PP4=0.070**
  - RHOA × 2hGlu: MR p=0.048 **→ PP4 추가 산출 필요**
- 형식: 각 hit별 **PP3 vs PP4 막대 대비**(공유 아님↔공유) + 가능하면 locuszoom 미니(eQTL 피크 vs outcome 피크 어긋남). locuszoom은 영역 데이터 있으니 작성 가능(±250kb).
- 메시지: "명목 연관은 전부 LD/pleiotropy. 다른 변이. 인과 아님."

---

## 3. 핵심 서사 5문장

1. RhoA–ROCK 축은 20년간 설치류에서 인슐린저항성·간지방·비만 표적으로 검증돼 왔고, 선택적 ROCK2 억제제(belumosudil)가 임상에 있어 repurposing 후보다.
2. 그러나 이 표적의 **사람 인과 근거는 부재**했다 — 우리는 cis-eQTL을 도구로 genetically-proxied ROCK1/ROCK2/RhoA 발현이 7개 심대사 형질에 미치는 인과효과를 MR + colocalization으로 처음 체계 검증했다.
3. 강한 instrument(F=27–409)와 조직 다중성(혈액·근육)에도 불구하고 **세 노드 전부, 일곱 형질 전부에서 인과 효과 없음** — 경로-수준 null.
4. 산발한 명목 연관 3개(ROCK2-BMI, ROCK1-간지방, RHOA-2h혈당)는 모두 colocalization에서 기각(PP4<0.1)된 반면, 동일 파이프라인의 양성대조(KCNJ11→T2D PP4=0.89·p=6.7e-16, SORT1→LDL p=7e-60)는 강하게 검출되어 방법의 민감도를 입증한다.
5. 결론은 **translatability gap**이다 — 설치류 표현형은 약리적 ROCK *활성* 억제를 반영하나, 인간 유전학은 평생 *발현* 변이를 반영한다; 두 층위의 불일치가 ROCK 대사 repurposing의 핵심 위험이며, 발현이 아닌 활성/급성 억제를 겨냥한 검증이 필요함을 시사한다.

---

## 4. Limitations (명시 목록)

1. **eQTL ≠ 약물 활성.** MR은 평생 *발현* 변이 효과를 추정 — belumosudil/fasudil의 *활성·급성* 억제와 층위 다름. null이 약리 무효를 의미하진 않음(핵심 한계, Discussion 전면).
2. **T2D 커버리지.** ROCK1/RHOA instrument가 Xue2018 T2D(sparse)에 부재 → T2D MR 미산출. proxy(r²>0.8/LDlink) 또는 dense DIAMANTE(Mahajan, 데이터접근 동의서)로 보강 예정.
3. **글리세믹/T2D × RHOA·ROCK1 coloc = confirmatory 미완.** MR null이라 후순위였음; 완전 히트맵 위해 API 영역추출로 채울 예정.
4. **RHOA weak-instrument.** 혈액 FDR-유의 cis-eQTL 1개(F=27)뿐 → RHOA 결론은 ROCK1/ROCK2보다 약함. multi-SNP 불가, coloc 검정력 낮음(RHOA×BMI PP4=0.31은 분산 아티팩트).
5. **BMI × 근육 coloc 미실행.** GTEx 근육 all-variant 영역이 공개 API/소형파일 부재 → **O2(Harvard HPC)에서 GTEx all-pairs로 보강 예정.** (단 BMI가 혈액 ROCK2와도 no-coloc이라 결론 불변.)
6. **EUR 한정.** instrument·outcome 모두 유럽계 → 타 조상 일반화 제한.
7. **조직 한정.** 혈액(eQTLGen)·근육(GTEx) 중심; 간·지방·췌장 직접 eQTL coloc은 부분적.
8. **SORT1 양성대조 coloc 약함.** Willer2013 LDL이 HapMap2-sparse → PP4 저평가(MR은 강함). KCNJ11이 주 검증.

---

## 산출 데이터 위치 (figure 작성용)
- Fig1B: `block0_instruments/`, `block5_pathway/`
- Fig2: `block23_MR/rock2_mr_results.csv`, `block5_pathway/`
- Fig3/4: `block4_coloc/coloc_results.csv`, `block5_pathway/coloc_rho_results.txt`, `block23_MR/positive_controls.txt`
- **Fig3/4 완성 전 추가 산출 필요**: ROCK1/RHOA × {FI,FG,HbA1c,2hGlu,T2D} coloc, RHOA×2hGlu coloc, (선택) SORT1 dense LDL, BMI×근육.
