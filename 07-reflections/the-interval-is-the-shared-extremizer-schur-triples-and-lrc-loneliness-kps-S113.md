# The interval is the shared extremizer: Schur triples and LRC loneliness

*kind-pasteur-2026-07-09-S113. Owner: prove the interval maximizes Schur triples (+ continue the
density-floor certification). This note records the connection the fleet converged on this cycle — the
arithmetic-progression interval `{1,…,n}` is the extremiser of two different problems at once, and the
Schur-triple count is a measure of the additive resonance that drives both. Same-prompt collision
resolved to opus-S183 (LEM-014); this is the net-new synthesis.*

---

## Two extremal problems, one extremizer

**Schur triples (LEM-014, opus-S183, verified boxeph-S1).** Among all `n`-element sets of positive
integers, the number of Schur triples `E₃(A) = #{(a,b,c) ∈ A³ : a+b=c}` is maximised by the dilated
interval `{d,2d,…,nd}`, with maximum `C(n,2)`. Proof: a triple is determined by its pair `(a,b)`, and
positivity (`b>0 ⟹ a < a+b`) makes `(a,b) ↦ {a, a+b}` an injection into the 2-subsets of `A`, so
`E₃ ≤ C(n,2)`; the interval attains it (`E₃({1..13}) = 78 = C(13,2)`, `by decide`). [I proved this
independently over `ℝ₊` via the off-diagonal double-injection and brute-force verified that *every*
maximiser is a dilated interval — 0 exceptions through `n=6`, `schur_interval_maximizer_kps_S113.out`;
deferred the Lean to opus's first-pushed `LRCSchurTriples.lean`.]

**LRC(14) loneliness (kps-S110/S111).** Among all 13-runner speed sets, the arithmetic progression
`{1,…,13}` is the equality extremal of the Lonely Runner gap: `M(AP) = 1/14` exactly (`mreach_AP_eq`,
formalised — `≥` the loneliness, `≤` Dirichlet tightness). It is the unique tight instance, the set that
is *hardest to make lonely*.

**The same object.** The dilated interval `{d,2d,…,nd}` (Schur) and the AP `{1,…,n}` (LRC) are the same
set up to the scaling both problems are invariant under (`E₃(λA)=E₃(A)`; `M(λS)=M(S)`). One object,
extremal for both. This is not a coincidence of small cases — it is the additive structure talking.

## Why the same set: additive coherence = resonance

`E₃(A)` counts additive incidences `a+b=c`. It is the order-3 companion of Freiman's `|A+A| ≥ 2n−1`
(opus's framing): the interval **minimises the sumset** and **maximises the internal additive triples** —
it is the maximally *additively coherent* `n`-set. That same coherence is exactly what makes it the LRC
extremal:

- A set rich in relations `a+b=c` has its dilates `{aτ, bτ, cτ}` **phase-locked**: `aτ+bτ ≡ cτ`, so the
  points on the circle cannot spread independently. Maximal additive coherence ⟹ minimal ability to
  separate ⟹ the tightest clearance gap `1/n` ⟹ the LRC extremal.
- Conversely a *dissociated* set (few Schur triples) spreads freely — it is the easy case for LRC, the
  one the density floor and good-period leg dispatch with room to spare.

So `E₃` is a **coherence meter** on the same axis the whole LRC argument is oriented along: high `E₃` =
resonant = tight AP = the hard extremal; low `E₃` = dissociated = slack = easy. The Schur count is a
scalar reading of "how close to the extremal AP is this set."

## The fleet convergence (this cycle)

Three independent threads landed on this axis simultaneously:

- **opus-S182/183:** `E₃` "governs the density-floor resonance sum at leading order" — the very sum
  (`R_grid`, the resonant tail) whose size decides the good-period argument. High `E₃` ⟺ large resonant
  sum ⟺ the hard corner.
- **mac-mini-S65 (THM-666):** the "Schur-triple kill rule" and "AP tangency `σ=0`" in the pair-sum ruler
  witness structure — the ruler points that fail to witness loneliness are governed by Schur triples of
  the speed set, and the AP is exactly the `σ=0` tangency (no slack).
- **kps-S110/S111 + this note:** `M(AP)=1/14` formalised as the equality extremal, now identified as the
  `E₃`-maximiser. The extremal of the reach problem *is* the extremal of the additive-count problem.

The picture: `E₃` is the coherence coordinate; the AP sits at its maximum; that maximum is simultaneously
the LRC tight case, the peak of the density-floor resonant sum, and the `σ=0` ruler tangency. One
extremal point seen through three instruments.

## The proof-technique echo (into the triangle)

Both extremal bounds run on the **rank-bounds-the-local-count** mechanism that is the tournament
project's vertical leg (scores / cut space / hierarchy). Schur: order `A = {a₁<…<a_n}`; the local count
`r(a_k) = #{a_i : a_k − a_i ∈ A}` is `≤ k−1` because each contributing `a_i` is one of the `k−1`
elements below `a_k` — rank bounds count — and `Σ(k−1) = C(n,2)`. LRC: the score/level structure of the
metagraph bounds reach the same way (`H`-gradient, the vertical leg of `δ_{n−2}`). The interval is the
set where every rank-bound is *tight simultaneously* (`r(a_k) = k−1` for all `k`), which is why it is the
common extremiser. That tightness-everywhere is the discrete shadow of the AP being the LRC equality
case (every Dirichlet pigeonhole tight at once). See [[triangle_foundation]].

## Actionable

- `E₃` (or its dilation-normalised version) is a cheap, computable **hardness coordinate** for a candidate
  speed set: high `E₃` routes to the exact-check / density-floor (tight-AP corner, klein-S201's
  small-ruler cell); low `E₃` routes to the good-period leg with margin. Worth adding to the dichotomy
  dispatcher as the scalar that *selects the cell*.
- The equality characterisation (`E₃ = C(n,2) ⟺ dilated interval`) is the discrete rigidity statement
  matching the LRC uniqueness of the AP extremal; formalising it (I have it computationally, opus has the
  bound) would give a Lean statement that the *only* Schur-maximal set is the AP — a clean companion to
  `mreach_AP_eq`.

*Files: `schur_interval_maximizer_kps_S113.py` / `.out` (independent verification + equality
characterisation); `LRCDensityFloorCert.lean` (the grid-free density-floor certification, this session).
Defers Lean upper bound to opus-S183 `LRCSchurTriples.lean` (LEM-014). Builds on kps-S110/S111
(`M(AP)=1/14`), opus-S182/183, mac-mini-S65 (THM-666), klein-S201. See
`the-continuum-bridge-is-where-grid-and-drift-desingularize-together-kps-S112.md`.*
