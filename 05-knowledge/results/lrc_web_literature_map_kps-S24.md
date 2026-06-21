# LRC(14) web-literature map: the repo's machinery IS the established frontier (kps-S24)

Web trawl of the current LRC literature, mapped onto the repo's sector route.

## The state of the art and WHY 14 is open
- **arXiv:2604.23906 (Rosenfeld-Trakulthongchai et al., Apr 2026)** proves up to 13 runners. **k=13 (14 runners)
  EXPLICITLY OPEN.** Polynomial method (Prop 4.4) certifies the canonical tuple (1,..,k) but needs **k+1 PRIME**
  (Fermat in the field Z_{k+1}); for k+1=14=2*7 composite it FAILS, and the c=14 lift is intractable (NOT
  CRT-factorable). => the repo's apex-prime-7 harmonic route is THE tool for the composite case (HYP-2758).

## The repo's machinery = the established methods, pushed for n=14
- **LP bounds via trig polynomials (arXiv:2010.02271, Brazilian Bull.)** = the repo's MOMENT DUAL (THM-534,
  meas(S7) <= L_y = sum_r y_r S_r). The kps-S24 decorrelated closed form `sum_t P_t^{(r)} p_t(B)` is exactly
  this LP, with the far-element surjection probabilities P_t^{(r)} as the explicit weights.
- **Tao's blog / "Some remarks" (arXiv:1701.02048):** the identified path is "variants of the
  INCLUSION-EXCLUSION formula combined with bounds on the size of HIGHER-RANK BOHR SETS." This IS the repo's
  approach: inclusion-exclusion over missed sectors (the p_t/W_a decomposition) + the RELATION LATTICE
  (Lambda(E) = the higher-rank Bohr sets = the resonances/commensurable far pairs, HYP-2757). Tao: only "slight
  improvements" known generally -- the repo pushes it for the SPECIFIC composite n=14 via apex-prime-7.
- **Riesz products (arXiv:2511.16636, Nov 2025):** constructs a product measure adapting to the speed-relation
  lattice -- the ASYMPTOTIC analog of the repo's DECORRELATED PRODUCT (the far elements as independent factors,
  Sigma P_t p_t). Asymptotic (not explicit n=14), but confirms the product/dissociation structure.
- **"Correlation among runners" (arXiv:1407.3381):** the runner-correlation = the repo's resonance/joint
  discrepancy (the L7 commensurable-curve correction to the decorrelated product).

## Net
The repo's sector route is NOT ad hoc: meas(S7)<=cap = the LP/moment-dual bound; the resonance lattice = Tao's
higher-rank Bohr sets; the decorrelated product = a finite-n Riesz product. The repo's novelty is doing all
three EXPLICITLY for the composite n=14 that the polynomial method cannot reach. The decorrelated main term is
the LP bound (closed form, <cap); the resonance correction is the Bohr-set/Riesz piece (HYP-2757 atlas). 
-> HYP-2758, HYP-2757, THM-534, OPEN-Q-108.
