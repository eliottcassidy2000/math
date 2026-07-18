# THM-999 — The bounded-denominator lemma; the residual R is proven for the tight locus (death-star-2026-07-17-S56)

**Status:**
- **Lemma A (PROVED, general):** every loneliness time of a tight family has denominator `≤ 2·Vmax`;
  in fact its reduced denominator divides `v₊ + v₋` for two active runners on opposite sides.
- **Corollary (PROVED):** the residual R ("no loneliness of denominator `> n`") is a FINITE, decidable
  check for any tight family — only the multiples of `n` in `(n, 2·Vmax]` need ruling out.
- **Theorem R for the tight locus (PROVED):** R holds for `AP` (no check needed) and for `GW`
  (only denom 28, 42 to rule out; both absent). Hence R holds for the entire known primitive tight
  locus `{AP, GW}`.
- **General primitive R:** decidable per family; a *uniform* proof is equivalent to bounding `Vmax`
  over primitive tight families — a facet of the tight-locus classification (THM-995 equality horn).

Answers the S56 residual of THM-997 ("primitive tight ⟹ no denominator `> n` loneliness"). Source
HYP-7305. Script `04-computation/lrc_bounded_denominator_deathstar_S56.py` (+ `.out`).
Notation: `Vmax = max_i |v_i|`, `f(t) = min_v ‖vt‖`, tight `⟺ M(V) = sup_t f(t) = 1/n`. WLOG speeds
positive (`‖vt‖` depends only on `|v|`).

---

## Lemma A — loneliness times of a tight family have bounded denominator

**Statement.** Let `V` be tight. Every loneliness time `t*` (i.e. `f(t*) = 1/n`) satisfies: there exist
two runners `v₊, v₋ ∈ V` with `v₊ t* ≡ +1/n` and `v₋ t* ≡ −1/n (mod 1)`, whence
`(v₊ + v₋) t* ≡ 0 (mod 1)` and `t* = k/(v₊ + v₋)`. In particular the reduced denominator of `t*`
divides `v₊ + v₋`, so `denom(t*) ≤ 2·Vmax`.

**Proof.**
1. **The maxima are isolated.** `f` is continuous and 1-periodic; `M(V) = 1/n` is its maximum. If `f`
   were flat at `1/n` on an interval, the safe set `{t : f(t) ≥ 1/n}` would have positive measure,
   contradicting tightness. So each maximizer `t*` is a *strict* local maximum.
2. **`f` is piecewise linear with nonzero slopes.** On a maximal interval where the argmin runner `v`
   and its nearest integer are constant, `f(t) = ±(vt − integer)` has slope `±v ≠ 0`. So `f` has no
   interior local max on a linear piece; every max is at a breakpoint.
3. **Two active runners, opposite sides.** For a strict local max, `f(t* + δ) < 1/n` for all small
   `δ ≠ 0`. Write each active runner's position as `v t* ≡ ε_v/n (mod 1)`, `ε_v ∈ {+1, −1}`
   (active `⟺ ‖vt*‖ = 1/n`). For small `δ`, `‖vt* + vδ‖ = |ε_v/n + vδ| = 1/n + ε_v v δ`. So an
   `ε = +1` runner makes `f` drop on the **left** (`δ<0`), an `ε = −1` runner on the **right**
   (`δ>0`). A strict local max needs `f` to drop on both sides, so there is an active `v₊` with
   `ε = +1` (`v₊ t* ≡ +1/n`) and an active `v₋` with `ε = −1` (`v₋ t* ≡ −1/n`).
4. **The sum vanishes.** `v₊ t* + v₋ t* ≡ 1/n − 1/n ≡ 0 (mod 1)`, so `(v₊ + v₋) t* ∈ ℤ`, i.e.
   `t* = k/(v₊ + v₋)`. The reduced denominator divides `v₊ + v₋ ≤ 2·Vmax`. ∎

**Verified** (script): every loneliness time of `AP` and `GW` with denominator `≤ 120` has denominator
exactly 14 (`≤ 2·Vmax = 26` and `48` resp.); the active pair at `t = 1/14` is `{1, 13}`, sum `= 14`.

*(This composes with THM-996-I: loneliness ⟹ `n | denom` and now `denom ≤ 2·Vmax`. Together the
loneliness locus of any tight family lives on the finite set of `t = k/(nm')`, `1 ≤ m' ≤ 2Vmax/n`.)*

---

## Corollary — R is a finite check; and the active-runner lemma

By THM-996-I, `denom(t*)` is a multiple of `n`; by Lemma A it is `≤ 2·Vmax`. So a tight family can be
lonely at denominator `> n` **only** at `nm'` with `2 ≤ m' ≤ ⌊2·Vmax/n⌋`. Ruling these out is a finite
computation. Moreover (THM-997, restated): at a denom-`nm'` loneliness time every active runner is
divisible by `m'` (`va ≡ ±m' (mod nm')` with `gcd(a,nm') = 1` forces `m' | v`), so the pair `v₊, v₋`
of Lemma A are both multiples of `m'` and `nm' | (v₊ + v₋)`.

---

## Theorem R for the tight locus `{AP, GW}`

**AP `= {1,…,13}` (`Vmax = 13`).** `2·Vmax = 26 < 28 = 2n`, so there is **no** multiple of 14 in
`(14, 26]`: by Lemma A no loneliness time can have denominator `> 14`. R holds with **no computation**.

**GW `= {1,…,11,13,24}` (`Vmax = 24`).** `2·Vmax = 48`; the multiples of 14 in `(14, 48]` are `28` and
`42`. A genuine denom-28 (resp. 42) loneliness time would make `liveCount(28) > 6` (resp.
`liveCount(42) > 6`), since the six base-resonance times of denom 14 already contribute 6 at every
resonance. But `liveCount(28) = liveCount(42) = 6` (verified). So GW has **no** loneliness time of
denominator 28 or 42, hence none `> 14`. R holds by a two-value check.

**Conclusion.** R is PROVED for both primitive tight families. Since the primitive tight locus is
`{AP, GW}` (the classification, THM-995 VII), R holds throughout it — and, crucially, the proof is now
a bounded finite check, not an unbounded verification.

---

## The general primitive case, honestly

For a hypothetical primitive tight family `V ∉ {AP, GW}`, Lemma A still makes R a finite check, but
`Vmax` is a priori unbounded, so the check is not uniform. The obstruction is genuinely global:

- A putative denom-`nm'` loneliness time (`m' ≥ 2`) needs two active runners `v₊, v₋` divisible by
  `m'` with `v₊ + v₋ = nm'·j ≤ 2·Vmax`; primitivity (`gcd V = 1`) forces some non-active runner `w`
  with `m' ∤ w`, sitting strictly inside the safe band. **No local contradiction arises** — the
  arrangement vertex and the two active constraints are consistent with such a `w`.
- Only *global* tightness (no better time anywhere) rules it out — which is exactly the tight-locus
  rigidity. Empirically no counterexample exists: `≈ 87k` primitive candidates at denom 28 and `≈ 73k`
  at denom 42 were tested (all pruned by `M > 1/14`), none tight; GW clean to `m = 150`.

**Verdict.** R is proven exactly as far as the tight locus is known: unconditionally for the
difference-closed primitive families (THM-997, via `m' | gcd`), and for `{AP, GW}` (Lemma A + finite
check). A uniform proof for *all* primitive tight families is equivalent to bounding `Vmax` over the
primitive tight locus, i.e. to the classification/rigidity core itself — so R beyond `{AP, GW}` is not
an independent gap but the same hard core as THM-995's equality horn. Consistent with THM-996-III: the
loneliness census separates the tight locus no further than the classification does.

→ THM-996, THM-997, THM-995 (VII), HYP-7305, MISTAKE-100.
