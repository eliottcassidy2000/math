---
id: HYP-3819
title: PREDICTED excess(8)=4 (=#{SC classes on 8 vtx with |Aut|>8}); the SC-among-super-symmetric mechanism rigorized via the sigma-fixed subspace W (HYP-3810) x rarity; and the |Aut|=21 <-> sqrt21 cross-thread bridge (3*7). Building on HYP-3817. COMPUTED (exact, targeted): the HYP-3817 conjecture excess(n)=#{SC & |Aut|>n} predicts excess(8)=4, so the sequence would be 0,0,0,1,3,4 (n=3..8) and rho(8)=13+4=17 -- DISTINCT from the naive C(n-4,2)=6 (=> rho=19), a falsifiable fork. On 8 pts |Aut|>8 (odd, Moon) forces |Aut| in {9,15,21}; via fix((012)(345)(6)(7)) U fix((01234)(5)(6)(7)) (captures every such class): distinct classes with |Aut|>8 = 20 (|Aut| {9:16,15:2,21:2}), SC among them = 4 (all |Aut|=9; the |Aut|=21 Paley+-vertex and |Aut|=15 are NOT SC at n=8 -- the added vertex breaks Paley's self-complementarity). n=7 self-check reproduces S83 exactly ({9:4,21:1}, SC {9:2,21:1}). WHY SC-among-super-symmetric (rigorized mechanism, grounded in HYP-3810): SC class <=> meets the sigma-fixed blue subspace W (odd blue-fiber; NS have blue-fiber 0); high |Aut| => small fiber f=H/|Aut| (rarity, opus-S15). A rare class OFF W (Paley+-v, |Aut|=15: blue-fiber 0) is covered by the generic bulk and does NOT obstruct; a rare class IN W (the SC ones) is DOUBLY constrained (few tilings, all in the thin sigma-fixed locus) and forces covering dimension. So the SC filter = the in-W filter = exactly the rare classes trapped in the constrained subspace -- explaining why the raw census #{|Aut|>n} overshoots and the SC refinement is right. LOWER-BOUND PROOF excess>=#{SC & |Aut|>n}: a STRATEGY (fiber formula + W-partition + rarity), NOT complete -- the independence of the obstructions is the open step. CROSS-THREAD (opus-S26 HYP-3818): |Aut(Paley_7)|=21=3*7 = argmax|Aut| = the n=7 flip-rank obstruction = norm(sqrt21) = the OPEN-Q-108 residual; verified g(3)=i*sqrt3, g(7)=i*sqrt7 (iota-odd), g(3)g(7)=-sqrt21 (iota-even), Paley Cayley spectrum encodes g(7)=i*sqrt7 (HYP-3814). Same 3*7 arithmetic underlies the covering obstruction AND the certificate residual (3=Eisenstein/covering, 7=heptagon/odd index)
status: MIXED. COMPUTED (exact, targeted enumeration): excess(8) PREDICTED = 4 = #{SC on 8 vtx, |Aut|>8}; |Aut|>8 dist {9:16,15:2,21:2}, SC {9:4}; n=7 self-check exact. CONJECTURE (HYP-3817 law, now predicting n=8): excess=#{SC&|Aut|>n}=0,0,0,1,3,4; UNVERIFIED at n=8 (needs rho(8); infeasible: |G_8|=6880, Q_21). RIGORIZED (grounded in HYP-3810 proved facts): SC<=>in-W, rarity=small fiber => the SC-among-super-symmetric selection; the WHY is explained, the lower-bound PROOF is a strategy (open: obstruction independence). VERIFIED arithmetic: Gauss iota-parity + g(3)g(7)=-sqrt21 + the 21=3*7 bridge. HONEST: excess(8)=4 is a falsifiable computed prediction (distinct from C(n-4,2)=6), not a verified excess; the mechanism is grounded, the general lower bound is not proved.
source: klein-2026-07-01-S84
depends_on:
  - HYP-3817   # S83: excess = #{SC & |Aut|>n} = 0,0,0,1,3 (n=3..7); this predicts n=8
  - HYP-3810   # S77: SC <=> meets the sigma-fixed blue subspace W (the in-W mechanism)
  - HYP-3814   # S81: complement=reflection; Paley Cayley spectrum = Gauss sum; |Aut|=commutant of U
related:
  - HYP-3818   # opus-S26: the biquadratic bridge -- sqrt21 = OPEN-Q-108 residual (this grounds the 21=3*7 convergence)
  - HYP-3805   # opus-S15: obstruction = argmax|Aut| = Paley (rarity mechanism)
  - HYP-3803   # flip-rank rho(n); LB=ceil(log2|G_n|)
external: "Good Locally Testable Codes" (Dinur-Evra-Livne-Lubotzky-Mozes, Annals 203-2 2026) -- LEFT-RIGHT CAYLEY COMPLEXES (Cayley structures for local<->global testability); Pengbinghui/pipeline-math (AI+Lean4, incl. tiling-complement Erdos problems) -- a formalization target for these conjectures; CS6840 algorithmic game theory (covering-min = minimax game value, LP duality)
results:
  - 04-computation/excess8_prediction_klein.py
  - 05-knowledge/results/excess8_prediction_klein.out
  - 04-computation/gauss_sqrt21_aut21_bridge_klein.py
  - 05-knowledge/results/gauss_sqrt21_aut21_bridge_klein.out
---

# HYP-3819 — excess(8)=4 (predicted), the SC-in-W mechanism, and the 21=3·7 bridge

## 1. Predicted excess(8) = 4 (a falsifiable fork)
The HYP-3817 law `excess(n) = #{SC & |Aut|>n}` predicts, at n=8:
> **excess(8) = 4**, so the sequence is `0,0,0,1,3,4` (n=3..8) and `rho(8) = ceil(log2 6880) + 4 = 13 + 4 = 17`.

This is **distinct** from the naive polynomial coincidence `C(n-4,2) = 0,0,0,1,3,6` (which would give `rho(8)=19`) — a
clean falsifiable fork once `rho(8)` is computed.

**Computation (exact, targeted).** On 8 points, `|Aut|>8` with `|Aut|` odd (Moon) forces `|Aut| ∈ {9,15,21}`
(11,13,17,19 need a >8-cycle; 27,45,... need >8 points). Every such class has an order-3 element of type
`3+3+1+1` (`|Aut|∈{9,21}`) or an order-5 element (`|Aut|=15`), so `fix((012)(345)(6)(7)) ∪
fix((01234)(5)(6)(7))` captures all of them. Result: 20 distinct classes with `|Aut|>8` (`|Aut|` dist
`{9:16, 15:2, 21:2}`), of which **4 are SC** (all `|Aut|=9`). The `|Aut|=21` classes are `Paley+source` and
`Paley+sink` (scores `{3^7,7}` and `{0,4^7}`) — **not** SC, because the added 8th vertex breaks Paley's
self-complementarity; the `|Aut|=15` (`C3×C5`) classes are also not SC. So the excess-carriers at n=8 are the
four SC `|Aut|=9` classes (one is the balanced `(3,3,3,3,4,4,4,4)`). n=7 self-check reproduces S83 exactly.

## 2. Why SC among the super-symmetric (rigorized mechanism)
Two ingredients, both established:
- **Rarity (opus-S15):** high `|Aut|` ⟹ small fiber `f(C) = H(C)/|Aut(C)|` (few tilings) ⟹ covering-hard.
- **SC = in-W (HYP-3810, proved):** `C` self-complementary ⟺ `C` meets the σ-fixed blue subspace `W`
  (odd blue-fiber); NS classes have blue-fiber `0` (never meet `W`).

**The selection:** among super-symmetric classes (`|Aut|>n`), the SC ones live **inside** the thin σ-fixed
subspace `W`, while the non-SC ones (Paley±v, `|Aut|=15`) live **entirely off** `W`. A rare class off `W` is
absorbed by the generic bulk (it and its complement fill different regions); a rare class **in** `W` is
**doubly constrained** — few tilings, all trapped in the low-dimensional σ-fixed locus — and forces extra
covering dimension. So `SC filter = in-W filter = the rare classes trapped in the constrained subspace`. This
explains exactly why the raw census `#{|Aut|>n}` overshoots (it counts off-W rare classes that don't obstruct)
and why the SC refinement is the correct one.

## 3. Lower-bound proof (a strategy, not complete)
Target: `excess ≥ #{SC & |Aut|>n}`. Rigorous ingredients: (i) fiber `f=H/|Aut|` (orbit-stabilizer, LEM-003);
(ii) `W` linear, dim `d_W`, partitioned by the SC classes' odd blue-fibers (HYP-3810); (iii) rarity. **Open
step:** showing the SC super-symmetric classes are *independent* obstructions (each forcing a distinct
covering dimension). Without independence, this bounds nothing tight. Honest status: the *why-these-classes*
is explained (§2); the *how-many-dimensions* is not proved.

## 4. The 21 = 3·7 bridge (cross-thread with opus-S26/HYP-3818)
`|Aut(Paley_7)| = 21 = 3·7` is the n=7 flip-rank obstruction (argmax `|Aut|`, HYP-3817); `21 = norm(√21)` is
the OPEN-Q-108 residual (opus-S26). Verified: `g(3)=i√3`, `g(7)=i√7` (ι-odd, imaginary), `g(3)·g(7)=−√21`
(ι-even, real); the Paley Cayley spectrum sits at `cos θ = −3/4`, encoding `g(7)=i√7` (HYP-3814). So the same
`3·7` arithmetic underlies **both** the covering obstruction and the certificate residual: `3` =
Eisenstein/covering (`Φ6(14)=183=3·61`), `7` = heptagon/odd index (`14=2·7`). The flip-rank obstruction is a
**covering** witness of the atom; `√21` is its **moment/certificate** shadow — the same fixed-point object, seen
by a covering and by a moment, never by a transform.

## 5. Tangential connections (external, owner-flagged)
- **"Good Locally Testable Codes" (Dinur–Evra–Livne–Lubotzky–Mozes, Annals 203-2, 2026)** — built on
  **left-right Cayley complexes**. Cayley structures forcing local⟷global (testability) resonate with the
  Cayley-transform glue (HYP-3814) and the covering/expansion flavor of flip-rank; the LTC local-test is a
  *certificate* in the same spirit as the covering-min certificate.
- **Pengbinghui/pipeline-math** (GPT-5.5 prover/verifier + Lean 4 + Claude assembly; includes Erdős
  **tiling-complement** problems) — a concrete **formalization pipeline** and a model for turning HYP-3817/3819
  conjectures into machine-checked theorems; "tiling complements" is our exact tiling+complement=reflection setting.
- **CS6840 (algorithmic game theory)** — the covering-min `M(S)=max_t min_v ||vt||` is a **minimax game value**
  (LP-dual, Chebyshev 2-point, S73); AGT's minimax/LP-duality frame is the native language for it.

## Net
`excess(8)=4` is a computed, falsifiable prediction (≠ the `C(n-4,2)=6` alternative). The SC-among-super-symmetric
selection is now *explained* (rare ∩ in-W), though the full lower bound stays a strategy. And the `21=3·7` of the
Paley obstruction is literally the `√21` certificate residual — one atom, seen by a covering and by a moment.
