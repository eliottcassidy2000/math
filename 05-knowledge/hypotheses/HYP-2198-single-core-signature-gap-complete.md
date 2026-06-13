# HYP-2198 — Single-core signature gap, complete for all lengths

**Status:** CONFIRMED
**Source:** monad-compute-2026-06-04-S2
**Resolves:** OPEN-Q-055 sub-question "Prove or refute the single-core signature gap: r_core(s) never equals 3 or 10."

## Setup

A *single-core complete-Ω* tournament is a transitive order `v_0 > v_1 > ... > v_{m-1}`
(`v_0` the source) together with one extra "core" vertex whose arc pattern is a bit
string `s ∈ {0,1}^m`:

- `s[i] = 1` : core → v_i  (core beats v_i)
- `s[i] = 0` : v_i → core

The number of **directed odd cycles through the core** is

```
r(s) = Σ_{i<j, s_i=1, s_j=0} f(j - i - 1),   f(0)=1,  f(t)=2^{t-1} for t≥1
```

(an odd cycle is `core → v_i → (increasing transitive chain) → v_j → core`; the interior
of the chain between positions `i` and `j` has `t = j-i-1` vertices, of which an
even-sized subset must be used, and `#even-subsets(t) = f(t)`).

When Ω is the **complete** conflict graph on these `r` cycles,
`H = I(K_r, 2) = 1 + 2r`. So:

| r  | H  | role |
|----|----|------|
| 3  | 7  | the K3 obstruction (THM-029) |
| 10 | 21 | the K10 obstruction (THM-079) |
| 31 | 63 | the THM-344 n=8 unlock |
| 94 | 189| the THM-025 n=9 count |

## The length bound (makes it a finite theorem)

Stripping **leading 0s** (a 0 with no 1 to its left contributes nothing) and
**trailing 1s** (a 1 with no 0 to its right contributes nothing) does not change `r`.
A nonzero value therefore has a **canonical** witness `s'` that starts with `1` and ends
with `0`. If `|s'| = L ≥ 3`, the (first-1, last-0) pair alone contributes
`f(L-2) = 2^{L-3}`, hence `r(s') ≥ 2^{L-3}`; if `L = 2`, `r ≥ 1`.

**Consequence.** Every achievable `r ∈ (0, R]` has a canonical witness of length
`L ≤ 3 + ⌊log₂ R⌋`. Enumerating *all* canonical strings up to that length therefore finds
*every* achievable value `≤ R` — so any value in `[1,R]` not found is provably
**unachievable at any length**, not merely "absent so far".

## Result (complete to R = 2^17 = 131072, length cap L = 20)

- **r = 3 (H=7): PERMANENT GAP.** Unreachable by any single-core signature of any length.
  (Its witness would have length ≤ 3+⌊log₂3⌋ = 4, all enumerated — fully general.)
- **r = 10 (H=21): PERMANENT GAP.** Witness would be length ≤ 6, all enumerated.
- **r = 31 (H=63): REACHABLE**, first at length 7 — matches THM-344's cores
  `1001100`, `1100110`.
- **r = 94 (H=189): GAP** (single-core; the actual THM-025 example is not single-core).
- The single-core gap set is **dense**: ~50% of `[1, 2^17]`. Smallest gaps
  `{3, 6, 10, 14, 17, 20, 21, 24, 27, 28, 29, 33, …}`.

This upgrades S11/S12's "no r=3 or r=10 through length 16 / 40" to a finite, all-length theorem.

## Interpretation

Single-core complete-Ω is a **strict sub-construction** of all tournaments (its gap set is
huge, whereas only `{7,21}` are globally forbidden H). The result therefore:

- **Explains** why H=63 unlocks via a single-core complete-Ω construction (r=31 reachable)
  while H=7 (r=3) and H=21 (r=10) cannot be produced this way.
- Does **NOT** by itself prove H=21 is globally forbidden — that is HYP-1753 / THM-079, which
  must rule out *all* Ω shapes (connected non-complete, multi-core, etc.), not just the
  single-core complete one. Cf. MISTAKE-024/050: H=63 *is* achievable; only `{7,21}` are
  permanent.

## Cross-checks (all pass)

- THM-344 single-core signatures `1001100 / 1100110 / 10011001 / 11001101` all give `r=31`.
- Canonicalisation r-invariance on assorted strings.
- `1·0^k → 2^{k-1}`.
- S11 convention (all `2^m` strings, `m=2..16`): r=3 absent, r=10 absent, r=31 first at m=7.

## Files

- `04-computation/single_core_signature_complete_monad_s2.py`
- `05-knowledge/results/single_core_signature_complete_monad_s2.out`
