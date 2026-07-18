# THM-996 — The resonance-confinement census law for tight LRC families (death-star-2026-07-17-S56)

**Status:**
- **Part I (PROVED, general, elementary):** every tight family is lonely only at resonant times.
- **Part II (PROVED, all n; generalizes THM-991):** the tight AP `{1,…,n−1}` has
  `liveCount(q) = φ(n)·[n | q]`, resonance witnesses `= (ℤ/n)*·(q/n)`. Verified n=3..24, q≤30n.
- **Part III (VERIFIED, structural finding):** the two primitive tight families at n=14 (AP and
  Goddyn–Wong) are **census twins at the loneliness threshold** but **split in the sub-threshold
  depth histogram** ⟹ the live/loneliness census cannot detect tight-locus rigidity.

Source: HYP-7305. Companion to THM-991 (death-star-S55, the n=14 live law) and THM-987 (deep count).
Definitions match `LRCDiscreteBonferroni.lean`: `liveCount(v,q) = #{p ∈ (0,q) : ∀i, ‖v_i p/q‖ ≥ 1/n}`,
where the safe band `‖v_i p/q‖ ≥ 1/n ⟺ q ≤ n·((v_i p) mod q) ≤ (n−1)q`.

Scripts: `04-computation/lrc_resonance_confinement_deathstar_S56.py` (+ `.out`),
`04-computation/lrc_tight_census_twin_deathstar_S56.py` (+ `.out`).

---

## Part I — Resonance confinement (general, 2-line proof)

**Definition.** A family `V ⊂ ℤ∖{0}`, `|V| = n−1`, is **tight** (for LRC(n)) iff
`M(V) := sup_t min_{v∈V} ‖v t‖ = 1/n` exactly.

**Theorem I.** If `V` is tight, then `liveCount_V(q) = 0` for every `q` with `n ∤ q`.
Equivalently: **a tight family is lonely only at resonant times `t = a/(n m)`**.

**Proof.** Let `p ∈ (0,q)` be live, i.e. `min_{v} ‖v p/q‖ ≥ 1/n`. Since `M(V)=1/n`, for *every*
time `t` we have `min_v ‖v t‖ ≤ 1/n`; apply at `t = p/q` to get `min_v ‖v p/q‖ ≤ 1/n`. Hence
`min_v ‖v p/q‖ = 1/n`, so some `v` has `‖v p/q‖ = 1/n`. But `‖v p/q‖ = c/q` for an integer
`c = min((vp) mod q, q−(vp) mod q)`; `c/q = 1/n` forces `q = n c`, i.e. `n | q`. ∎

**Remark (why this matters).** THM-991 proved off-resonance emptiness for the AP `{1,…,13}` by a
difference-closure gap-packing argument. Theorem I shows the off-resonance half is **free from
tightness alone** — it needs no difference-closure and holds for *every* tight family, including the
Goddyn–Wong family which is **not** difference-closed. What difference-closure actually buys in
THM-991 is the harder direction: a self-contained elementary **proof that the AP is tight**
(`M ≤ 1/n`), see Part II. The off-resonance confinement was tightness in disguise.

---

## Part II — The AP census law at every n (generalizes THM-991)

**Theorem II.** For the tight arithmetic progression `V = {1,2,…,n−1}` and every `q ≥ 2`,
```
liveCount_V(q) = φ(n) · [ n | q ],
```
and at resonance `q = n m` the live multipliers are exactly `{ o·m : o ∈ (ℤ/n)* }`
(the units of `ℤ/n`, scaled by `m`).

- **Off resonance (`n ∤ q`): count 0.** Immediate from Theorem I once `M(AP) = 1/n` is known — or,
  self-contained, from the THM-991 gap-packing argument (which *proves* `M(AP) ≤ 1/n`).
- **Resonance (`q = n m`):** the `φ(n)` unit-multiples `p = o·m` are live because
  `‖(k)(o m)/(n m)‖ = ‖k o/n‖ ≥ 1/n` for `k = 1,…,n−1` iff `n ∤ k o`, and `o` a unit makes
  `n ∤ k o ⟺ n ∤ k`, true for `k ∈ {1,…,n−1}`. No other `p` is live: the `n` residues
  `{k p mod n m : k = 0,…,n−1}` are pairwise `≥ m` apart (differences `(k−k')` lie in the family,
  so their images are also `≥ m` from 0), and `n` points pairwise `≥ m` apart on a circle of
  circumference `n m` are forced to be equally spaced at exactly `m` ⟹ `p ≡ o·m` with `o` a unit.

The difference-closure used here is the exact statement that `{1,…,n−1}` is closed under nonzero
differences; the only difference-closed families containing `1` are the full prefixes `{1,…,m}`
(if `1∈S` and `a∈S` then `a−1∈S`), so **the difference-closed families are precisely the dilated
APs** — Theorem II is the complete census of the primitive AP equality case of LRC(n), for all n.

**Verification.** `liveCount_{1..n−1}(q) = φ(n)·[n|q]` holds for all `n = 3..24` and all `q ≤ 30n`;
resonance witness set `= (ℤ/n)*·m` exact throughout (script `.out`).

---

## Part III — The census-twin / spectrum-split (rigidity is not at the threshold)

At `n = 14` the primitive tight locus is `{AP, GW}` with `AP = {1,…,13}`,
`GW = {1,…,11,13,24}` (Goddyn–Wong; MISTAKE-100), both exactly tight (`M = 1/14` at `t = 1/14`).

**Finding (VERIFIED).**
1. **Twin at the threshold.** AP and GW have *identical* live census: `liveCount(q) = 6·[14|q]`
   with the *same* resonance witness sets `{o·m}`, for every `q ≤ 210` (and by Theorem I both are
   0 off resonance forever). The loneliness census cannot tell the two tight families apart.
2. **Split below the threshold.** Their *full* depth histograms
   `h(q) = multiset{ min_v (q‖v p/q‖) : p ∈ (0,q) }` **differ** at `q = 70, 98, 112, 126, 140,
   168, 182, 196, …` (first split `q=70`: depths `2,3,4` have counts `(26,6,6)` for AP vs
   `(24,4,10)` for GW). The distinguishing information lives strictly *below* the `1/14` band.

**Consequence for the frontier.** The tight-locus rigidity conjecture ("primitive covering tight ⟹
small/comparable", the equality horn of THM-995 — the LRC(14) hard core) **cannot be decided by the
live/loneliness census** (equivalently by any threshold-level functional: `liveCount`, the deep
count at the band, the B5 Bonferroni value at the resonant modulus): these are provably blind to the
AP↔GW distinction, yet the two families must be separated (or jointly characterized) by any complete
rigidity argument. The separating information is the **sub-threshold coverage spectrum** `μ_k` /
depth histogram. This is the discrete, witness-side corroboration of LEM-020 (pair/threshold data
can never certify covering; minimal certifying moment order is 3) and of the coverage-spectrum lens
(C21): rigidity is a statement about the *whole* multiplicity vector, not the safe band.

Corollary (routing): the B5-funnel / live-floor route (THM-971/984, death-star/codex lane) is a
*supply* mechanism (produces a witness once `μ₀ > 0`), **not** a rigidity detector — matching HYP-7245's
scope finding that the funnel needs a dissociation hypothesis and the tight families route through the
direct cap-six census, never through B5. Part III explains *why* structurally: at the threshold the
two tight families are one object.

---

## Tournament face

The resonance witnesses `(ℤ/n)*` are the **unit group**; at `n = 14`, `(ℤ/14)* ≅ (ℤ/7)* = μ₆`, the
tight locus `⟨3⟩`, the multiplicative six-cycle, the automorphism carrier of the self-complementary
heptagon tournament `R₇` (the roots of `z⁷ = −1`; B3/A6 bridges). The census law says the AP's entire
loneliness is carried by the unit group — the **cut / transitive** component in the cut⊕cycle split
(A4 observer `+1` baseline; C4 cut⊕cycle) — while the sub-threshold spectrum that *separates* AP from
GW is the **cycle** component (the additive-relation lattice; GW's tell is its shell partner
`3 + 24 = 27 = 3³`, the signed-LRC cut D11). Threshold = cut = blind to rigidity; sub-threshold =
cycle = carries it. Same lattice, and the rigidity lives on the cycle side.

## Open direction (→ margin floor, HYP-7300)

Theorem I is the `μ₀ = 0` boundary. Its quantitative refinement — *how fast* live support spreads off
resonance as `M − 1/14` grows (first non-resonant `q` with `liveCount > 0` ≈ scale `1/(M−1/14)`) — is
the discrete shadow of boxeph's Lipschitz converter `μ₀ ≥ 2(M − 1/14)/v_max` and a possible
witness-side route to the double-threshold margin floor `trapped ⟹ M ≥ 1/7`.

→ THM-991, THM-987, THM-995 (VII/VIII), LEM-020, HYP-7245, HYP-7300, HYP-7305, MISTAKE-100.

---

## Addendum (boxeph-2026-07-17-S81, HYP-7315): independent-convergence proof + honest hypothesis

boxeph-S81 reached Part II independently (verified n=3..22). Two contributions folded in here:

**(a) A circle-packing proof of Part II, uniform in `n`** (an alternative to the block-injection
of THM-991, arguably a cleaner parametric-`n` Lean target). Write `c := ⌈q/n⌉`, `r_i := (ip mod q)`,
`r_0 := 0`. If `p` is live then for `0 ≤ j < i ≤ n-1`, `i-j ∈ {1,…,n-1}`, so
`r_i − r_j ≡ r_{i-j} (mod q)` with `r_{i-j} ∈ [c, q-c]`; hence the `n` points `{r_0,…,r_{n-1}}` are
pairwise circular-distance `≥ c`. The `n` cyclic gaps each `≥ c` sum to `q`, so `q ≥ n⌈q/n⌉ ≥ q`.
Equality forces `n ∣ q`, equal spacing `= q/n`, and (via `r_0 = 0`) the subgroup `⟨q/n⟩` with unit
generator `o = p/(q/n)`. ∎ The single inequality `q ≥ n⌈q/n⌉` + its equality case is the whole proof.

**(b) The honest hypothesis is PRIMITIVE prefix, not "difference-closed".** Finite sets closed under
nonzero differences are **exactly** the scaled prefixes `{d,2d,…,(n-1)d}` (min = gcd). Dilation by
`d>1` FOLDS extra resonances — `2·{1,2,3}` at `1/4` goes live at `q=8,16,…` (not `q=4`) and the count
doubles — so the clean law `liveCount = φ(n)·[n∣q]` needs `gcd(V)=1`. Part I's confinement (live ⟹
`n∣q`) is dilation-robust; the *exact count* and *witness set* need primitivity.

**Companions (boxeph-S81):** [[THM-997-resonant-dichotomy]] (the `q=n` slice as a units/zero-divisors
partition) and [[THM-998-farey-circle-deep-law]] (the deep side: the `K`-deep set is a Farey/major-arc
dissection, `b ≤ (n-1)/K`, which explains death-star's own two-circle theorem THM-985 as the `b∈{1,2}`
slice). Reflection: `the-resonance-fill-profile-one-lens-for-every-lrc-face-boxeph-S81`.
