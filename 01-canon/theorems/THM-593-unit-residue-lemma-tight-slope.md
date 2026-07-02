# THM-593: The unit-residue lemma for tight sets + the exact tight-set slope formula

**Status:** PROVED (elementary; exact-rational verification on all known tight sets)
**Author:** mac-mini-2026-07-01-S93
**Verification:** `04-computation/lrc_radius_slope_frame_macmini_S93.py`, `04-computation/lrc_tight13_locus_slope_floor_macmini_S93.py` (+ .out in `05-knowledge/results/`).

---

## Setting

`S ⊂ Z_{>0}` finite, `q ≥ 2`. Call `S` **q-tight** if `M(S) = max_t min_{v∈S} ||vt|| = 1/q` exactly. (For LRC(N) with `|S| = N−1` and `q = N`, tight sets are the conjectured extremals — the equality locus.) Assume `S` contains no multiple of `q` (all known/classified tight sets: residue systems in `{1,…,q−1}`, HYP-3750).

## Part A — the unit-residue lemma

**Lemma.** If `M(S) = 1/q` and no `v ∈ S` is `≡ 0 (mod q)`, then for **every** `a` with `gcd(a,q) = 1` the residue set `R_a = {av mod q : v ∈ S}` contains **both** `+1` and `−1`.

**Proof.** Fix a unit `a` and suppose `1 ∉ R_a`. Take `t_ε = a/q − ε` with `0 < ε < 1/(q·max S)`. For `v` with residue `r_v ∈ {2,…,q−2}`: `||v t_ε|| ≥ min(r_v/q − vε, 1 − r_v/q + vε) ≥ 2/q − vε > 1/q`. For `r_v = q−1`: `||v t_ε|| = 1/q + vε > 1/q`. No residue is `0` or `1`, so `min_v ||v t_ε|| > 1/q = M(S)`, contradicting the definition of `M` as the supremum. Hence `1 ∈ R_a`; applying this at the unit `−a` gives `1 ∈ R_{−a} = −R_a`, i.e. `−1 ∈ R_a`. ∎

**Corollary A1 (all units represented).** Taking `a` over all units: every unit residue `u ∈ (Z/q)^*` is represented in `S mod q`. Hence `|S| ≥ φ(q)`, and a tight `(q−1)`-set has exactly `q−1−φ(q)` "free" (non-unit) slots.

**Corollary A2 (every unit point is a witness).** For tight `S`, at every unit `a` the point `a/q` is a **witness**: `min_v ||v·a/q|| = 1/q` exactly (residue 1 is present by the lemma and no residue is 0). So tight sets have at least `φ(q)` witnesses, at the unit fractions.

**Corollary A3 (drops are non-units).** In the duplication+drop classification of non-difference-closed tight sets (klein-S48, HYP-3750: residues = `{1..q−1}` with one residue dropped, one duplicated), the **dropped residue is always a non-unit**.

**Corollary A4 (prime modulus ⟹ only permutation type).** If `q` is prime, all nonzero residues are units, so no residue can be dropped: every tight set's residues are a permutation of `{1,…,q−1}`. This **explains the census zeros** of klein-S48 at `n = 4, 6` (`q = 5, 7` prime) with no enumeration. (At `q = 9` the census zero is a second-order fact — drops of `{3,6}` are lemma-legal but no completion covers; verified by sweep.)

## Part B — the exact slope of a tight set

By THM-592(ii) and Corollary A2, near `r = 1/q` the lonely set of a tight `S` is a union of one shrinking component per unit witness `a/q`, bounded on each side by the runners achieving residues `±1` at `a`; when a residue `u` is carried by several runners the **binding endpoint owner is the fastest**, `v_max(u) = max{v ∈ S : v ≡ u (mod q)}`. Each witness component has width `(1/q − r)(1/v_max(a^{-1}) + 1/v_max(−a^{-1}))`. Summing over units (as `a` runs over units so does `a^{-1}`):

```
m_S(r) = c_S · (1 − q r)   on the last linearity cell below 1/q,  with

c_S = (2/q) · Σ_{u ∈ (Z/q)^*}  1 / v_max(u).          (exact)
```

**Corollary B1 (AP value & the S92 constant).** For `S = {1,…,q−1}`: `v_max(u) = u`, so
`c_AP(q) = (2/q) Σ_{u∈(Z/q)^*} 1/u`. At `q = 14`: `c_AP = 2(1/13 + 1/45 + 1/33) = 1666/6435 = 0.258897…` — the exact value of the empirical "0.26" slope of HYP-3824/S92.

**Corollary B2 (slope rigidity iff lifts avoid unit residues).** `c_S = c_AP(q)` iff every unit residue is carried by its minimal representative (`v_max(u) = u`), i.e. all lifted elements sit on **non-unit** residues. Verified rigid: AP; GW `{1,2,3,4,5,7,12}` (q=8, lift on residue 4); GW `{1,…,11,13,24}` (q=14, lift on residue 10); cross `{1,3,4,5,9}` (q=6, lift on residue 3).

**Corollary B3 (rigidity FAILS in general — beaters exist).** The q=8 tight set `{1,4,5,6,7,11,13}` (drop 2, dup 5 via 13, relift 3→11 — the relift is a forced covering patch) carries units 3 and 5 on fast runners 11 and 13:
`c = (2/8)(1/1 + 1/7 + 1/11 + 1/13) = 328/1001 = 0.3277 < c_AP(8) = 44/105 = 0.4190.`
So the AP is **not** universally the slope minimizer over the tight locus; the slope floor is a classification-dependent quantity `(2/q)·min_{tight S} Σ_u 1/v_max(u)`.

**Corollary B4 (q = 14).** Exhaustive sweep over the classified families (single lifts j ≤ 12, double dup+drop/relifts j ≤ 4): the only non-AP tight 13-set is the GW set `{1,…,11,13,24}`, which is rigid. Hence at `N = 14` the slope floor over the enumerated tight locus is exactly
`c_min(14) = c_AP(14) = 1666/6435`,
confirming and sharpening HYP-3824's `inf_S |L_S(r)| ≈ 0.26·(1−14r)`: the constant is `1666/6435` and it is attained by BOTH extremals (AP and GW). [Sweep: `lrc_tight13_locus_slope_floor_macmini_S93.out`. Caveat: exhaustiveness is relative to the HYP-3750 residue classification and the lift bounds; higher lifts verified to break tightness upward (M > 1/q, LRC-safe) in all sampled families.]

## Interpretation

Tightness is a **critical covering**: danger arcs at radius exactly `1/q` cover the circle, and every unit fraction is pinned by a `±1` residue pair. Modifying an element opens a hole elsewhere (verified: plain drop-2-dup-5 at q=8 has a hole at `t = 9/20`, `M = 3/20`), and further lifts must patch exactly the hole they open — the domino constraint that makes the tight locus finite and sporadic, and the quantitative content of "the atom is a fixed point with no measure of its own" (S92 reflection): at `r = 1/q` the extremal lonely set is `φ(q)` isolated points; the measure lives entirely in the slope, and the slope is an explicit unit-harmonic sum.

-> THM-592, HYP-3840, HYP-3750, HYP-3824, HYP-2893 (GW/Jacobsthal), OPEN-Q-108.
