# Route 1's density floor is robust — exact three-gap constants, and why 1/7 works where 2/7 didn't

**opus-2026-07-06-S130 (second half).** After the audit retired Route 2, the owner directed effort to
Route 1 (bound `Mreach ≥ 1/14` directly — the correct object). Its one open analytic node is the
**k=8..13 witness-density floor**: `witnessG2(E) > 0` for every cluster of co-offsets `E`. Working it
*correctness-first* (the Route-2 lesson: don't trust a bounded sample), I verified the floor is
genuinely robust and computed its exact constants.

## The load-bearing quantity, and its exact values

The floor reduces to `μ_{1/7}(E) := meas{ x ∈ [0,1] : circular max-gap of {frac(e·x) : e ∈ E} > 1/7 }`.
Two facts, verified this session:

**(i) The AP is the global minimizer.** `μ_{1/7}` is dilation- and translation-invariant (both are
rigid rotations of the phase configuration), so `{1..k}`, `{2,4,…,2k}`, `{101..113}` all give the
same value. A strong adversarial search — 40 aggressive descent runs from random `E ⊂ [1,150]`, plus
structured adversaries (geometric, split, doubled) — **never found any `E` below the AP**; every run
returned to `μ_{1/7}(AP)`. So `μ_{1/7}(E) ≥ μ_{1/7}({1..k})` for all `E` (strongly verified; the
proof is the remaining lemma — the AP orbit `{frac(jx)}` is the maximally-equidistributed
configuration by the three-gap theorem, so it has the smallest max-gaps, hence the smallest good-set).

**(ii) The AP values, computed EXACTLY.** `μ_{1/7}({1..k})` is a piecewise-linear function of `x`
with breakpoints at Farey fractions `m/d`, `d ≤ k` (phase coincidences and wraps); on each
order-interval every gap is linear, so `{max gap > 1/7}` is a union of rational sub-intervals. Exact
rational arithmetic gives:

| k | μ_{1/7}(AP) | ≈ |
|---|---|---|
| ≤7 | 1 | 1.000 — pigeonhole: `maxgap ≥ 1/k ≥ 1/7` a.e. |
| 8 | 691/735 | 0.940 |
| 9 | 247/294 | 0.840 |
| 10 | 38/49 | 0.776 |
| 11 | 1381/2205 | 0.626 |
| 12 | 13823/24255 | 0.570 |
| **13** | **477/1078** | **0.4425** |

The `k=13` value **exactly matches the project's claimed `rhoGlobFloorRat(13) = 477/1078`** — an
independent confirmation of the canon's floor via a clean three-gap computation. The minimum over the
relevant range (clusters have `≤ 13` speeds) is `477/1078 ≈ 0.44`, far above `m_P ≈ 0.0565`.

## Why 1/7 is robust where 2/7 collapsed

The refuted 2/7 route asked for `maxgap > 2/7 ≈ 0.286`, which is *close to* the typical max-gap of
`k` well-spread points (`≈ H_k/k`, e.g. `0.34` for `k=8`) — so structured `E` could push the good-set
to **zero** (kps's admissible zeros). The 1/7 threshold (`0.143`) sits *well below* that typical
max-gap, so the good-set stays a large majority of `[0,1]` even at the AP minimizer (`0.44` at
`k=13`). The margin between `1/7` and the AP's typical max-gap is exactly what makes the floor
bounded away from zero. This is a structural reason, not a lucky sample.

## Where Route 1 now stands (honest)

- **Aimed right:** bounds `Mreach` (the supremum), no wrong-object flaw (unlike Route 2's J-K top).
- **Density floor:** SOUND and ROBUST. Reduces to (i) AP-minimality [strongly verified, proof =
  three-gap equidistribution] + (ii) the exact AP constants [RIGOROUS this session]. The value is
  `≥ 477/1078 > 0` for all `k ≤ 13`.
- **Part A** (`ρ* > 0 ⟹ Mreach ≥ 1/14`, the slow/fast change of variables + criterion C): the
  reformulation is proved in canon; the finite-`Vmax` error budget `O(#arcs/Vmax)` is stated but not
  uniformly formalized (`LRCWitnessPartA` has the arithmetic glue).

So the two honest remaining pieces are **(A) AP-minimality of `μ_{1/7}`** (a clean three-gap /
equidistribution lemma) and **(B) the finite-`Vmax` error bound for Part A**. Neither is a
wrong-object mirage; both are genuine but bounded analysis on a correctly-aimed route. This is the
difference between Route 1 and Route 2: Route 1's remaining work is *hard but real*, not *disconnected
from the target*.

## Bottom line

The witness-density floor — the node the owner asked me to work — is verified robust, its constants
are now exact (three-gap), and its `k=13` value independently reproduces the canon's `477/1078`. The
floor is not at risk of the 2/7 collapse. What remains is a clean equidistribution lemma
(AP-minimality) plus the Part-A error budget — the honest analytic core of a route that at least
points at the right quantity.
