---
id: THM-1271
title: Gate completeness, three-quarters proved — the RUNG-FLOOR LEMMA (closed-form witness a = D·p⁻¹, unique danger slot = the deleted N−1, so M(F_D(N)) ≥ D/((N+1)D−1) whenever gcd(2D−1, Q) = 1), the SMALL-MODULI SEAL (every candidate at q′ ≤ 2N has value ≤ 1/(N+1), N odd), and the PACKING SEAL (any candidate with e ≥ c has value ≤ 1/(N+1)) — jointly: the (N−1)-slot e-channel at pair-sum moduli > 2N is the SOLE escape, collapsing exact M(F_D(N)) to a ~10⁵-op computation; the pruned proof-backed evaluator reproduces the full table and decides the D=7 rung at N=2311.
status: >
  PROVED — Lemmas 1, 2, 3 (elementary; full proofs in this file), and the resulting
  reduction of M(F_D(N)) to the e-channel maximum (given THM-1002). VERIFIED-EXACT —
  the closed-form witness on every (D,N) in {3..7}×[3D−1,100] with gcd(p,Q)=1 (zero
  failures); L2 exhaustively at N=31/61/91 (zero violations); the pruned evaluator
  against the full evaluator on every odd-N table row (gate). The N=2311 verdict:
  see §5 (filled from the run). REMAINING for full gate completeness: the e-channel
  branch analysis in general (the per-N scan IS a finite certificate; the general-N
  law stays mechanism-derived).
source: death-star-2026-07-19-S59d (HYP-7905; owner: prove gate completeness and test the D=7 rung at N=2311)
depends_on:
  - THM-1002  # maximizer at pair-sum moduli (the candidate-completeness input)
  - THM-1286  # the gate law being completed; the tower context
related: [THM-1285, HYP-7900, HYP-4516, THM-996]
scripts:
  - 04-computation/lrc_gate_completeness_and_N2311_deathstar_S59d.py -> 05-knowledge/results/lrc_gate_completeness_and_N2311_deathstar_S59d.out
---

# THM-1271 — the rung floor, the two seals, and the e-channel reduction

Throughout: `D ≥ 3`, `p = 2D−1`, `N > 3D−2`, `Q = (N+1)D−1`, `x = D(N−1)`,
`F = F_D(N) = {1..N}∖{N−1} ∪ {x}`, base `B = {1..N}∖{N−1}`, and `|y|_q` the distance
of `y mod q` to `0`. A *candidate* is a pair `(q′, a)`; its *value* is
`min_{v∈F} |va|_{q′} / q′`. `M(F) = max_t min_v ‖vt‖`.

## 1. Lemma 1 (rung floor; closed-form witness)

> If `gcd(p, Q) = 1`, then at `a ≡ D·p⁻¹ (mod Q)` every `v ∈ F` has `|va|_Q ≥ D`,
> with equality at `v = p` and `v = x`. Hence `M(F_D(N)) ≥ D/Q`.

*Proof.* `(N+1)D = Q+1 ≡ 1 (mod Q)`, so `D⁻¹ ≡ N+1`, and
`p(N+1) = (2D−1)(N+1) = 2D(N+1) − (N+1) ≡ 2 − (N+1) = 1−N (mod Q)`.
Suppose `v·a ≡ r` with `|r| ≤ D−1`. Then
`v ≡ r·a⁻¹ ≡ r·p·D⁻¹ ≡ r·p(N+1) ≡ −r(N−1) (mod Q)`.
- `r = 0`: `Q | v`, impossible (`0 < v ≤ x = Q−p < Q`).
- `r ∈ [1, D−1]`: the representative is `Q − r(N−1) = N(D−r) + D+r−1 > N`, and it is
  `< x` (their difference is `N − 3D + 2 + (r−1)(N−1)·…` ≥ `N−3D+2 > 0`), and `= x`
  would force `r(N−1) = p`, i.e. `r | p` with `r ≤ D−1 < p`; then `r = 1`, `N = 2D`,
  excluded. So `v ∉ F`.
- `r ∈ [−(D−1), −1]`: the representative is `|r|(N−1) ≤ (D−1)(N−1) < Q`; it lies in
  `F` only if `|r| = 1` (giving `v = N−1`, the DELETED element — not in `F`) since
  for `|r| ≥ 2` it exceeds `N` and equals `x` only at `|r| = D`, excluded.
So no `v ∈ F` is within `D−1` of `0`. Equality: `p·a ≡ D` and
`x·a = (Q−p)a ≡ −D`. ∎

Two structural corollaries. (i) **The canonical defect position is forced**: the
unique danger slot of this witness is `u ≡ −(−1)(N−1) = N−1` — the family works
*because* it deletes exactly the element the arithmetic condemns. (ii) The floor
holds at EVERY `N` with `gcd(p,Q) = 1` — e.g. `M({1..11,13,48}) ≥ 4/55` at N=13 —
the gate governs only whether competitors push `M` above the window, never whether
the rung value is available from below.

## 2. Lemma 2 (small-moduli seal; N odd)

> For N odd, every candidate with `q′ ≤ 2N` has value `≤ 1/(N+1)` (below the window).

*Proof.* `q′ ≤ N`: if `q′ ≠ N−1` then `q′ ∈ B` gives distance `0`; if `q′ = N−1`
then `x = D(N−1) ≡ 0`. Value `0`. Now `N < q′ ≤ 2N`. By Dirichlet (boxes
`0..⌈q′/2⌉`), there is `1 ≤ u₀ ≤ ⌈q′/2⌉ ≤ N` with `|u₀a|_{q′} ≤ q′/(⌈q′/2⌉+1) < 2`,
so `≤ 1`. If `u₀ ≠ N−1`: `u₀ ∈ B`, value `≤ 1/q′ ≤ 1/(N+1)`. If `u₀ = N−1` with
`|u₀a| = 0`: `x ≡ 0`, value 0. If `u₀ = N−1` with `|(N−1)a| = 1`: then
`gcd(N−1, q′) = 1`; since `N−1` is even (N odd), `q′` must be odd. Set
`u′ = q′ − (N−1) ∈ [2, N+1]`; `u′ = N+1` forces `q′ = 2N` (even, excluded) and
`u′ = N−1` forces `q′ = 2N−2` (even, excluded), so `u′ ∈ B`, and
`|u′a| = |−(N−1)a| = 1`. ∎

## 3. Lemma 3 (packing seal) and the e-channel reduction

> If a candidate `(q′, a)` has `e := |(N−1)a|_{q′} ≥ c := min_{u∈B} |ua|_{q′}`,
> then `(N+1)c ≤ q′`; its value is `≤ c/q′ ≤ 1/(N+1)`.

*Proof.* Consider the `N+1` points `{u·a mod q′ : u = 0, 1, …, N}` (indices
`{0} ∪ B ∪ {N−1}` = `{0,…,N}`). Every pairwise index difference `d ∈ [1, N]` is
either in `B` (`|da| ≥ c`) or `d = N−1` (`|da| = e ≥ c`), so the points are pairwise
`≥ min(c,e) = c` apart on a circle of circumference `q′`: `(N+1)c ≤ q′`. The
family's value is `≤ c/q′` (some base element realizes `c`). ∎

**Reduction.** `M(F) ≥ D/Q > 1/(N+1)` (Lemma 1). By THM-1002 the maximizer sits at
`(q*, a*)` with `q*` dividing a pair sum of `F`. By Lemma 2, `q* > 2N`; by Lemma 3,
`e* < c*`; then the value is `min(c*, |xa*|)/q*` with `|xa*| = D·e*` (valid since
`D·e* ≤ D·q*/N < q*/2`), and the `N`-point packing on `{0} ∪ B` (differences again
exhaust `[1,N]`) gives `e* ≤ q*/N`. Writing the same time non-reduced at the full
pair sum `S` (`q* | S`, `S ∈ {x+u : u ∈ B} ∪ {2x}` — base sums are `≤ 2N−1`,
sealed), the congruence `(N−1)a ≡ ±e (mod S)` with `gcd(N−1,S) | e ≤ S/N`
enumerates every surviving candidate. Hence

```text
M(F_D(N)) = max( D/Q ,  max over that finite e-channel of min(c, De)/S ),
```

a `~10⁵`-operation exact computation even at `N = 2311`. At `S = 2x`,
`gcd(N−1, 2x) = N−1 > S/N` — the channel is empty and `2x` is sealed outright; the
rung itself reappears as the `e = 1` solution at `S = Q = x+p` (its `a` is Lemma 1's
witness), a consistency check the code confirms.

## 4. What this proves, and what remains

- Every candidate reaching the window travels the `(N−1)`-slot channel: the ONLY way
  `F_D(N)` can be beaten below-or-above is through small `|(N−1)a|` — the deleted
  element's ghost. Gate arithmetic is exactly the bookkeeping of which moduli admit
  a profitable ghost.
- Per-`N` gate completeness is now a FINITE CERTIFICATE (the e-channel scan); the
  pruned evaluator built from Lemmas 1–3 reproduces the full evaluator on every
  odd-`N` table row (script gate), and its member outputs are proof-backed exact
  values, upgrading `4/127, 4/247, 4/367, 6/1271` from evaluator-exact to
  lemma-plus-finite-certificate exact.
- Remaining OPEN for the general-`N` gate law: a uniform description of the
  e-channel maxima (the branch analyses of THM-1286 §3, now localized to the ghost
  channel). Even-`N` needs a variant of Lemma 2 (the `q′ = 2N` corner genuinely
  hosts the `b=2`-type loose configuration — that is WHY even `N` fails the gate,
  so the odd-`N` restriction is not a loss for members).

## 5. The D=7 rung at N=2311

`F_7(2311) = {1..2309, 2311} ∪ {16170}`, `p = 13`, `Q = 16183`, `gcd(13, Q) = 1`,
window `W_2311 = (1/2312, 2/4623)` of width `≈ 9.4·10⁻⁸`. Gate check:
`2311 ≡ 1 (mod 2310 = L_7)` and `2311 ≡ 10 ≢ 1 (mod 13)` — predicted OPEN with
`M = 7/16183`. Result:

```text
M(F_7(2311)) = 7/16183  EXACTLY — rung ATTAINED, IN-WINDOW.
Witness: a = 12449 = 7·13⁻¹ mod 16183 (Lemma 1's closed form; min distance 7).
Computed in 3 seconds by the pruned evaluator (naive cost ~10¹⁰ ops);
gate: 142/142 table rows reproduced exactly, 0 mismatches.
```

**The third consecutive out-of-sample confirmation of the tower** (after 4/247 at
N=61 and 6/1271 at N=211), this one at a 2311-speed scale inside a window one part
in 10⁷ wide. Confirmed rungs now ride the primorials 6, 30, 210, 2310 (D = 3, 4, 6,
7). The next rung is D=9 (binder 17; D=8's binder 15 = 3·5 composite, skipped like
D=5): predicted first opening at N ≡ 1 (mod 30030), N ≢ 1 (mod 17) — N = 30031 =
59·509; the same machinery computes it (the e-channel scan scales linearly).
