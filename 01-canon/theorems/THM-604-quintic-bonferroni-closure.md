# THM-604: The depth-5 Bonferroni closure — the free-phase anti-covering floor is positive at ALL cluster sizes j ≤ 13

> [renumbered from THM-599 by klein-S99 numbering cleanup (owner-directed); THM-599 = the torus-band theorem (kps, first reserve).]

**Status:** PROVED (truncation identity + exact arithmetic; d-fold forced-independence hypotheses as in THM-598 extended to depth 5)
**Author:** mac-mini-2026-07-01-S99 (HYP-3855)
**Role:** closes the "j ≥ 8 cubic leg" — the single remaining quantitative rung of the hpartA local-covering program (S98). With THM-598 (pairs) and this theorem, the free-phase local-covering floor exists at every cluster size the LRC(14) residual can present.

## The truncation identity

For integer-valued coverage `C(x) = Σ 1_{C_{w_i}(φ_i)}(x)` on the window `I`, the inclusion–exclusion identity `1_{C=0} = Σ_d (−1)^d (C choose d)` truncates with alternating error sign: for ODD depth `D`,
```
1_{C=0} ≥ Σ_{d=0}^{D} (−1)^d (C choose d),
```
since the remainder is `(−1)^{D+1} (C−1 choose D+1)-type ≥ 0` for `C ≥ 1` and exact at `C = 0`. Integrating over `I`:
```
u ≥ Σ_{d=0}^{D} (−1)^d S_d,   S_d = Σ_{d-subsets T ⊆ F} |∩_{w∈T} C_w ∩ I|.
```
Even-`d` sums enter positively (need LOWER bounds — forced independence, THM-598's hard direction); odd-`d` sums negatively (need UPPER bounds — no inflated resonances). Under **depth-5 resolution** (no frozen dangerous pattern in any `d ≤ 5` resonance direction `Σ m_i w_i ≈ 0` with small coefficient-product height — finite lists per `d`, thresholds by the same `1/(3·height)` calculus as THM-598(A); kps-S31's torus-band theorem supplies the exact `d`-fold volumes on the ledger side), every `S_d = C(j,d)(2r)^d(1 ± ε_d)` with explicit `ε`.

## The closure

With `2r = 1/7`, define `B_D(j) = Σ_{d≤D} (−1)^d C(j,d) 7^{−d}`. Exact values:

| j | B₃ (cubic) | B₅ (quintic) | B₇ | (6/7)^j |
|---|---|---|---|---|
| 7 | 0.32653 | 0.33986 | 0.33992 | 0.33992 |
| 8 | 0.26531 | 0.29113 | 0.29136 | 0.29136 |
| 11 | 0.06997 | 0.17992 | 0.18345 | 0.18348 |
| 12 | **−0.00875** | 0.15029 | 0.15719 | 0.15727 |
| 13 | **−0.09913** | **0.12209 = 2052/16807 = 2052/7⁵** | 0.13459 | 0.13480 |

**The cubic truncation FAILS at j = 12, 13** (this is why the "cubic leg" could not have closed as posed); the **quintic is positive at every j ≤ 13**, with the worst case
```
u/|I| ≥ B₅(13) − ε = 2052/7⁵ − ε ≈ 0.1221 − ε.
```
As `D → ∞`, `B_D(j) → (6/7)^j`: the floor converges to the independent-coverage value — the resolved free-phase adversary can do no better than random.

## Assembly (the hpartA local-covering program, now complete in structure)

- `j ≤ 6`: mass (union bound), unconditional.
- `j = 7`: THM-598(D) quadratic floor `6/49` (or B₅(7) = 0.3399 under depth-5 resolution).
- `8 ≤ j ≤ 13`: THIS theorem, `u/|I| ≥ B₅(j) − ε > 0`.
- frozen/resonant tuples at any depth: renormalize to the finite fixed-offset pattern lists (THM-598(C)) or the difference core (HYP-3901) — `j` strictly decreases; terminates.
- BCS-compliance (sLRC audit, HYP-3855): the floor is quantified over all phases but ONLY under the resolution hypothesis; the shifted-LRC counterexamples (BCS 2603.24784, false from n = 5) live exactly in the frozen/renormalization column (near-equal tiling clusters, THM-598(C) necessity) — no ∀-shift strength is claimed where it is false.

## Remaining constant-chasing (explicitly scoped)

(i) the depth-`d ≤ 5` dangerous-pattern lists and their `ε_d` thresholds (finite symbolic enumeration — kps-S31's "symbolic k ≤ 13 ledger" is the c-averaged version; the free-phase version is the same lattice enumeration); (ii) the window-boundary terms (one arc per comb per resonance direction). Both are finite and explicit; neither is analytic.

-> THM-598, HYP-3854, HYP-3953/3954 (kps), HYP-3847/3848 (klein triple branch), HYP-3901, THM-527, OPEN-Q-108.
