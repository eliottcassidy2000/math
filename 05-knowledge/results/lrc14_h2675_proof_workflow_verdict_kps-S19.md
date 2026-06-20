# HYP-2675 proof workflow — final adversarial verdict (kps-S19)

4-angle workflow (cross-scale-decorrelation, plateau-recursion, coverage-deficit-direct, dichotomy-finite) +
independent verifiers + completeness critic, on the SOLE LRC(14) residual `span(E)>B ⟹ p0(E) ≤ cap_k`.
All four returned REDUCTION; 3 of 4 verifiers returned holds=false on the strong forms. EXACT-Fraction throughout.

## REFUTED (over-claims caught by adversarial verification)
- **(W*) wide ⟹ p0 ≤ Q(k−1):** REFUTED. Witness `E=[0,19,20,21,22,23,24,25]` (k=8): p0 = 9524621/47108600
  = 0.20218 > Q(7) = 289/1470 = 0.19660 (still < cap_8). Q(k−1) is the decorrelated LIMIT, NOT a finite upper
  bound. (= kps's MISTAKE-080; the "= Q(k-1)" crux of lrc14_h2675_decorrelation_foundation was over-claimed.)
- **Plateau-recursion induction premise (wide E' ⇒ p0 ≤ Q(k−2)):** REFUTED. `[0,12,18,24,30,36,42]`:
  p0 = 0.0456 > Q(6) = 11/294 = 0.0374.
- **"Plat(E')≤Q(k−1) PROVED via THM-535":** MISLABEL — this is HYP-2603 (AP/consec-extremality of the plateau),
  a CONJECTURE (exhaustively VERIFIED at small span, NOT proved). THM-535 proves only cap_k ≥ (k−6)/7 + the
  specific consec meas value.

## PROVED (genuine advances, re-verified exact)
- **L1 cardinality lemma:** at any x, |E| points hit ≤ |E| sectors; a cluster of size s≤5 covers all 6 inner
  sectors on measure EXACTLY 0 (q(s)=0, s≤5; q(6)≤0.0526, q(7)≤0.2624). Per-cluster coverage deficit δ_0 ≥ 0.7376.
- **Comb bound (THM-546/547, SHARPER than the σ-bound):** for E=E'∪{w}, w=max, |Δ_w| ≤ 2·c1(E')/(7w) = (6/49)V(E')/w,
  c1(E') = #interval-components of the miss-exactly-one region of E'. Tooth-counting proof; 0 violations over
  thousands of exact peels (loose ~5–11×). This is the rigorous one-far decorrelation RATE.
- caps (THM-535), margins μ_k = cap_k − Q(k−1) > 0 (μ_min = 129643/980980 ≈ 0.132 at k=9), scale-invariance
  p0(λE)=p0(E) (THM-531-B): PROVED/VERIFIED.

## The CORRECTED residual (CONJECTURE, target the CAP level not the Q level)
> **wide ⟹ p0(E) ≤ cap_k** — 0 counterexamples in ~10⁴–10⁵ adversarial wide primitive sets; max wide p0 ≈
> 0.073/0.144/0.251/0.352/0.448 (k=8..12), **margins ≥ 0.30 to cap** (vs thin ≤0.06 to Q — the Q-level is the
> wrong, refuted target). Reduces to the CAP-level joint decorrelation bound
> `p0(E) ≤ p0_inf(shapes) + Σ_far (6/49)V_i/g_i` with an explicit JOINT multi-dim Erdős–Turán–Koksma constant
> for the carrier vector (frac(c_i·x))_i on the r-torus. The r=2/one-far case IS the proved comb bound (THM-546);
> only the **multi-cluster constant** + the balanced-wide branch (c1 grows with span) remain OPEN.

## Highest-value next step (synthesis)
Prove the CAP-level joint bound directly (margin ≥0.30 ⟹ even a LOSSY constant suffices — the tractable
product/IE problem the C(k) workflow predicted, NOT signed cancellation): (a) extend THM-546's tooth-counting
to two far clusters simultaneously (2-D comb bound); (b) bound c1 against the cap margin (5× the Q budget, so
its span-growth no longer breaks the cutoff); (c) run the small finite check on p0_inf < cap over the ≥2-cluster
shape family (already exact: max p0_inf = Q(k−1), margin ≥0.132). → HYP-2675, HYP-2603, THM-546/547, MISTAKE-080.

## CANON marks
HYP-2675 (wide⟹p0≤cap): CONJECTURE. (W*) wide⟹p0≤Q(k−1): REFUTED. THM-546/547 comb bound: PROVED.
L1 cardinality: PROVED. HYP-2603 (consec maximizes Plat): CONJECTURE (verified small). LRC(14): NOT PROVED.
