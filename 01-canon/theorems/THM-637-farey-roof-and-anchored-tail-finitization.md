# THM-637 — The Farey roof (AP max-gap mechanism) and the anchored-tail finitization of μ_θ

**Source:** opus-2026-07-07-S134; **proofs completed opus-2026-07-07-S135** (self-contained, no
external citations; every lemma machine-verified in exact rational arithmetic, 0 violations).
**Status:** parts (a),(b),(c),(d) **PROVED** (complete proofs below); part (e) VERIFIED-computational
(candidate lemma, OPEN).
**Load-bearing for:** kps-S59's diameter floor (k=13 leg of `hlarge` for primitive diameter ≤ 75 —
the "modulo roof" caveat there is discharged by this file), monad-S2's exact crossing, klein-S154's
LOO composite.
**Scripts:** `lrc_farey_roof_apmaxgap_opus_S134.py`, `lrc_roof_proof_verification_opus_S135.py`
(step-by-step proof verification, exact rationals, k ≤ 40), plus the S134 anchored-tail suite.

Throughout: `config(E, x) = {frac(e·x) : e ∈ E}` on the circle ℝ/ℤ, `maxgap` = largest circular gap,
`gap∋a` = the circular gap containing the point `a`, `μ_θ(E) = Leb{x ∈ [0,1] : maxgap > θ}`.
`F_k` = the Farey fractions of order k in [0,1]. For an **open Farey-k cell** `(p/q, p′/q′)`
(consecutive in `F_k`): `qp′ − pq′ = 1`, `q + q′ > k`, and no fraction with denominator ≤ k lies
strictly inside. Write `ε = x − p/q`, `ε′ = p′/q′ − x` (so `ε + ε′ = 1/(qq′)`, `ε, ε′ > 0`).
Note: for `x` in an open cell and `1 ≤ i ≤ k`, `ix ∉ ℤ` (else `x = a/i ∈ F_k`), so `frac(ix) ∈ (0,1)`
and `frac(−ix) = 1 − frac(ix)`.

---

## Lemma A (first return from above)

**For `x` in the open Farey-k cell `(p/q, p′/q′)` and every `1 ≤ i ≤ k`:
`frac(ix) ≥ frac(qx) = qε`, with equality iff `q | i` and `i = q`** (multiples `i = tq` give `tqε`,
minimized at `t = 1`; non-multiples are strictly larger).

*Proof.* Let `a = ⌊ix⌋`, so `frac(ix) = ix − a` and `a/i ≤ x`. The fraction `a/i` has denominator
`i ≤ k`, so it cannot lie strictly inside the cell; since `a/i ≤ x < p′/q′`, we get **`a/i ≤ p/q`**.

*Case `a/i = p/q`.* Lowest terms force `q | i`, say `i = tq`, `t ≥ 1`; then `a = tp` and
`frac(ix) = t(qx − p) = t·qε ≥ qε`, equality iff `t = 1`.

*Case `a/i < p/q`.* Let `m := ip − aq ≥ 1`. Then `frac(ix) = ix − a = iε + (ip − aq)/q = iε + m/q`.
Suppose for contradiction `frac(ix) < qε` (or `= qε`; the argument below excludes both). Then
`m/q ≤ (q − i)ε`, which forces `i < q`, and since `ε < 1/(qq′)` (open cell):
`m/q < (q − i)/(qq′)`, i.e. **`mq′ + i < q`**, and trivially `mq′ + i ≥ 1·1 + 1 > 0`.
Now the determinant identity `qp′ = pq′ + 1` gives
> `q·(ip′ − aq′) = i(pq′ + 1) − aqq′ = q′(ip − aq) + i = mq′ + i`,

so `q` divides `mq′ + i` — impossible for `0 < mq′ + i < q`. ∎
*(Boundary cells `p/q = 0/1` or `p′/q′ = 1/1`: the case split degenerates but the statement is
immediate — e.g. on `(0, 1/k)`, `a = 0` for all `i` and `frac(ix) = ix ≥ x = qε` with `q = 1`.)*

## Lemma B (first return from below)

**For `x` in the open cell: `1 − frac(ix) ≥ q′ε′` for every `1 ≤ i ≤ k`, with the minimum
`q′ε′` attained at `i = q′`.** Equivalently `max_i frac(ix) = frac(q′x) = 1 − q′ε′`.

*Proof.* The reflection `x ↦ 1 − x` is an order-reversing bijection of `F_k` preserving
consecutivity, and maps the cell to the open Farey-k cell `((q′−p′)/q′, (q−p)/q)` whose **left**
denominator is `q′`, with `(1−x) − (q′−p′)/q′ = ε′`. Since `ix ∉ ℤ`,
`frac(i(1−x)) = frac(−ix) = 1 − frac(ix)`. Apply Lemma A on the reflected cell:
`1 − frac(ix) = frac(i(1−x)) ≥ q′ε′`, equality at `i = q′`. ∎

## Lemma C (gap bound — no three-distance fine structure needed)

**Every circular gap of `config({1..k}, x)` has length ≤ `qε + q′ε′`.**

*Proof.* Fix a config point `u = frac(ax)`, `1 ≤ a ≤ k`; we exhibit a config point at circular
distance ≤ `qε + q′ε′` above `u` (this bounds the gap above `u`, and every gap lies above some
config point):
- if `a + q ≤ k`: `frac((a+q)x) = u + frac(qx) (mod 1)` sits at distance `qε`;
- else if `a − q′ ≥ 1`: `frac((a−q′)x) = u + (1 − frac(q′x)) (mod 1)` sits at distance `q′ε′`;
- else `k − q < a ≤ q′` (the two failures are compatible — the window `(k−q, q′]` is nonempty
  precisely because `q + q′ > k`). Then `b := a + q − q′` satisfies `1 ≤ b ≤ k`
  (from `a > k − q`: `b > k − q′ ≥ 0`, so `b ≥ 1`; from `a ≤ q′`: `b ≤ q ≤ k`), and
  `frac(bx) = u + qε + q′ε′ (mod 1)` sits at distance exactly `qε + q′ε′ < 1`. ∎

## Theorem (the Farey roof)

**For `x` in the open Farey-k cell `(p/q, p′/q′)`:**
> **`maxgap(config({1..k}, x)) = gap∋0 = qε + q′ε′ = q(x − p/q) + q′(p′/q′ − x)`** —
> linear on the cell, interpolating the node values `1/q → 1/q′`.

*Proof.* By Lemmas A and B, the config's minimum is `qε` (at `frac(qx)`) and its maximum is
`1 − q′ε′` (at `frac(q′x)`), so the gap containing the point 0 has length exactly
`qε + q′ε′`. By Lemma C no gap is longer. ∎

**Machine verification (S135):** exact `Fraction` arithmetic at random rational `x` with
denominators ~10⁶ inside random cells — Lemma A + the divisibility identity: 0 violations
(k = 5, 8, 13, 21, 40; 300 cells each); roof = gap∋0 = maxgap and Lemma-C case-split coverage:
0 failures (k = 3..13, 400 cells each). The theorem also gives, on each cell, the identity
`maxgap = m⁺ + m⁻` (two-sided closest approach to 0 — the observer-lens form; klein-S153's
a.e. numeric, now an identity on every open cell).

## Corollaries (all exact; machine-verified)

1. **Mean:** `E[maxgap(AP_k)] = Σ_cells ∫ roof = Σ 1/(q·q′²) = (1/2)Σ(q+q′)/(q²q′²)` over
   consecutive Farey-k pairs. Values: `1, 3/4, 7/12, 1/2, 5/12, 23/60, 1/3, 43/140, 47/168,
   19/72, 151/630, 796/3465, 93/440` (k = 1..13). `E[min_i frac(ix)] = half` (reflection).
2. **Tail:** `μ_θ(AP_k)` = the roof's superlevel measure — per-cell linear tail, exact rationals.
   Reproduces the canon `691/735, 247/294, 38/49, 1381/2205, 13823/24255, 477/1078` (θ = 1/7,
   k = 8..13). The good set of the AP is **exactly the q ≤ ⌊1/θ⌋ Farey neighborhoods**: pure
   resonance-window measure, zero bulk. At θ = 1/7, `q = 7` nodes are exactly marginal — **the
   apex prime 7 is invisible to the density floor.**
3. **Diameter floor (with kps-S59's subset lemma — now fully proved end-to-end):** for any
   integer set `E` with primitive diameter `D`: `μ_θ(E) ≥ μ_θ(AP_{D+1})` and
   `E[maxgap(E)] ≥ A(D+1)`. Exact crossings (recomputed independently, S135, matching kps-S59 &
   monad-S2): `μ_{1/7}(AP_n) ≥ m_P = 14249/252252` through **n = 76**
   (`μ(AP₇₆) = 2314528732/40290957525 ≈ 0.057445`), `< m_P` from n = 77; `A(n) > 1/7` through
   **n = 22** (`A(22) = 3029671/21162960`), `≤ 1/7` from n = 23. Hence **every 13-element E with
   primitive diameter ≤ 75 has `μ_{1/7}(E) ≥ m_P`** — the k=13 leg of `hlarge` on bounded
   diameter, now proved with no outstanding caveats.

## (d) The anchored-tail finitization (exact, every E) — PROVED

`F₇` has maximal spacing `1/7` (attained at `(0, 1/7)`), so any circular gap of length `> 1/7`
contains a point of `F₇`; conversely anchored gaps are gaps. Hence for every finite integer set:
> `{x : max_{a∈F₇} gap∋a > 1/7} = {x : maxgap > 1/7}` up to boundary, so
> **`μ_{1/7}(E) = P_x(max_{a∈F₇} gap∋a > 1/7)` — an 18-anchor closest-approach statistic.**

## (e) OPEN candidate (A″-F₆) — unchanged from S134

`t_{F₆}(E) ≥ μ_{1/7}(AP_k)` for affine-normalized k-sets, equality iff AP; implies (A′). Survived
normalized corpus, descents (k = 13, 12, 10), exhaustive 1-swap and sampled 2-swap scans; the
window-exactness mechanism (each inter-cluster gap contains its flanking `j/q`) and the escape
confinement `(0,1/6) ∪ (5/6,1)` are the structural reasons. Not proved; see HYP-4782.

## Formalization status (updated opus-S135)

**The pointwise roof is GREEN:** `TournamentH7/LRCFareyRoof.lean` (kernel-pure, axioms
`[propext, Classical.choice, Quot.sound]`, no sorry/native_decide, in the root manifest) —
`no_middle_fraction`, `lemmaA`, `lemmaB` (the divisibility-contradiction proofs), `lemmaC`
(three exhibited indices), `zero_gap_empty` (the 0-gap of length roof is config-free),
`fract_q_mul`. Hypotheses are the cleared forms (`p < q·x`, `q′·x < p′`, `q·p′ − p·q′ = 1`,
`k < q + q′`); all Farey facts derived, no continued fractions, no measure theory.

**Remaining seam to a fully GREEN k=13 bounded-diameter floor:** mac-mini-S42's
`LRCTailDiameter.lean` (GREEN) carries the diameter chain conditional on one certificate Prop —
the **AP₇₆ Farey ledger** `μ_{1/7}(AP₇₆) ≥ m_P`. That certificate is a finite rational
computation consuming these pointwise theorems: per Farey-76 cell, the superlevel of the linear
roof is an interval with rational endpoints (`fract_q_mul` + `zero_gap_empty` + `lemmaC` pin the
maxgap function; the cell sum is exact arithmetic). Discharging it closes the loop.
