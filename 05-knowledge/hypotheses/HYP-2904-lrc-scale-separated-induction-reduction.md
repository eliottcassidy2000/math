# HYP-2904: LRC scale-separated induction needs a component budget

**Status:** PROOF FRAGMENT / induction guardrail
**Owner:** codex-2026-06-22
**Tags:** LRC, induction, finite-comb, equidistribution, scale separation, Node 3
**Scripts:** `04-computation/lrc_scale_separated_induction_codex.py`
**Output:** `05-knowledge/results/lrc_scale_separated_induction_codex.out`
**Depends on:** HYP-2900, HYP-2901, HYP-2902, HYP-2903, HYP-2895, THM-523, THM-565, THM-566, KPS-S31v, KPS-S31w

## Statement

Use the workspace convention: `LRC(n)` means `n-1` speeds and threshold
`1/n`.  Let `B` be a seed speed set and let

`A = Safe_n(B) = {t in [0,1): ||b t|| >= 1/n for every b in B}`.

If `A` is a finite union of `c` intervals and has measure `mu>0`, then for
every added integer speed `v`,

`meas(A cap Safe_n({v})) >= (1-2/n) mu - 2c/v`.

Consequently the extension `B union {v}` is LRC-safe at threshold `1/n` for
every

`v > 2c / ((1-2/n) mu)`.

This is a genuine inductive reduction for the large-speed / Node-3 branch:
once the smaller seed has positive safe measure at the next threshold, every
sufficiently large new speed can be added by an exact finite-comb estimate.

## Proof

The unsafe set for the new speed is the periodic comb

`U_v = {t : ||v t|| < 1/n}`.

In each period of length `1/v`, its measure is exactly `2/(nv)`, so a whole
period contributes the fixed fraction `2/n`.  On a seed-safe interval `I`, all
whole periods inside `I` therefore contribute exactly `(2/n)` times their
length.  Only the two boundary partial periods of `I` can deviate from that
fraction, and their total length is at most `2/v`.  Thus

`meas(I cap U_v) <= (2/n) meas(I) + 2/v`.

Summing over the `c` components of `A` gives

`meas(A cap U_v) <= (2/n) mu + 2c/v`,

and subtracting from `mu` gives the displayed lower bound for
`A cap Safe_n({v})`.

If `LRC(n-1)` is already known for the seed size, then `A` has positive measure
at threshold `1/n`: the `LRC(n-1)` witness is at distance at least
`1/(n-1)>1/n` from every seed runner, so continuity gives a nonempty open
interval inside `A`.  Thus every fixed seed has a finite large-speed threshold.

## Exact audit

The script checks the bound with exact rational interval arithmetic.  For the
S46 AP-core seed

`B={1,2,3,4,5,6,7,8,9,10,11,13}`, `n=14`,

the seed safe set has `4` components and measure `426/35035`.  The comb bound
certifies every added speed `v>=768`; in particular it certifies the committed
radical/lcm speeds:

| added speed | exact after adding | comb lower bound |
|---:|---:|---:|
| `30030` | `313/30030` | `7472/735735` |
| `60060` | `337/32340` | `1514/147147` |
| `510510` | `191/18326` | `26032/2501499` |

Sub-14 training rows show the same mechanism.  For seed `{1,...,8}` at
`n=10`, the least comb-certified speed is `305`: `v=210` is safe in fact but
not by the crude boundary bound, while `v=420` is certified.  For seed
`{1,...,6}` at `n=8`, the least comb-certified speed is `183`, and `v=840`
is certified.

## Relation to KPS-S31v

Incoming KPS-S31v proves the natural multi-large extension for the Node-3
branch.  If `G_C` is a bounded core lonely set and `U_v` is one danger comb,
then the same boundary-period argument gives

`meas(G_C cap U_v) <= (1/7) meas(G_C) + arcCount(G_C)/(7v)`

at threshold `1/14`.  A union bound yields

`meas(G_C minus union_i U_{v_i}) >= (1-r/7)meas(G_C) - (arcCount(G_C)/7) sum_i 1/v_i`.

For `r<=6`, the leading coefficient remains positive, so sufficiently large
speeds cannot cover the bounded core.  Thus HYP-2904 is the one-speed
induction atom; KPS-S31v is its Bonferroni-1 / multi-large package.  Both
carry the same invariant: a positive seed floor plus a component or arc-count
budget.  For `r>=7`, the union bound becomes vacuous and the remaining target
is the second-moment / resonant-pair defect ledger.

## Relation to KPS-S31w

Incoming KPS-S31w organizes the same atom into a global reduction tree:

- **R1 remove-large:** if one speed is scale-separated from the rest, peel it
  by HYP-2904 / KPS-S31v and reduce the size by one while multiplying safe
  measure by approximately `6/7`.
- **R2 omit-prime:** if the set misses all multiples of some prime `p<=13`,
  then `t=1/p` is already a direct `>1/14` witness.
- **R3 dilation:** normalize by `M(dS)=M(S)`.

Iterating R1 peels a scale hierarchy down to a bounded core and uses the
already-proven smaller LRC cases for the seed floors.  Thus the reduction
tree is:

- `unbounded or scale-separated -> smaller LRC seed`;
- `missing small prime -> direct witness`;
- `bounded covering core -> finite Node-2 extremality`.

HYP-2904 supplies the rigorous local inequality needed for R1.  KPS-S31w
explains the largest possible payoff of size induction: it kills the
unbounded/non-covering cases, but it cannot by itself resolve the bounded
covering core where all speeds are comparable.  That remaining base is the
three-gap / AP-hull / Legendre-Venn finite extremality problem.

## S117 local one-interval sharpening

HYP-2906 gives the existence-first version of this induction atom.  If the
seed has one witness margin `alpha>1/n` and max speed `m`, then the seed is
safe on a connected interval of length `2(alpha-1/n)/m`.  The added speed's
unsafe comb has tooth length `2/(nv)`, so `v>m/(n(alpha-1/n))` is enough.

Using only `LRC(n-1)` gives `alpha=1/(n-1)` and therefore the sharp symbolic
gate `v>(n-1)m`.  In the LRC14 row this is:

`v_max > 13 v_second  ->  safe by LRC13`.

This is stronger than the global component-budget threshold for one-speed
existence.  It is weaker in a different direction: it gives a witness, not a
positive measure floor, and it does not handle several large speeds whose
individual ratios stay below the gate.  So the proof order should be:

- first peel a locally large top speed by HYP-2906;
- then use HYP-2904/KPS-S31v for density or multi-large budgets;
- finally route top-balanced bounded cores to Node-2 extremality.

## Guardrail

This is not a pure runner-count induction.  Dilation preserves the seed safe
measure but multiplies its component count.  For the same AP-core seed at
`n=14`, scaling by `q` keeps measure `426/35035` while the component counts and
least certified speeds grow linearly:

| scale | components | least certified `v` |
|---:|---:|---:|
| `1` | `4` | `768` |
| `2` | `8` | `1536` |
| `5` | `20` | `3838` |
| `10` | `40` | `7676` |
| `25` | `100` | `19190` |
| `50` | `200` | `38380` |

So `LRC(n-1)` alone gives positivity for each fixed seed, but it does not give
a uniform threshold depending only on `n`.  Any real induction must carry one
of the following extra invariants:

- a scale-normalized component/arc budget;
- a bounded-core Node-2 atlas before the large-speed handoff;
- an exact-period analytic estimate replacing the crude `2c/v` boundary loss.

## Tournament Analysis

Vertices are proof carriers, not runners.  The useful pairwise observable is:
which carrier preserves the LRC predicate while retaining component, scale, and
threshold information.

The audited tournament ranks

`finite_comb_budget > scale_normalized_seed > node3_large_speed > node2_bounded_core > pure_size_induction > raw_runner_vertices`.

This is transitive in the script (`directed_3cycles=0`).  The challenged
assumption is that induction should quotient only by number of runners.  The
quotient must remember the topology of the safe set: interval components are
the exact place where the large-speed Weyl/equidistribution error is paid.

## Consequence for LRC14

HYP-2904 sharpens the S46 Node-3 induction skeleton.  The large committed
speed branch is analytically reducible to smaller-size seed facts, but the
reduction is scale-separated and topology-aware.  The remaining proof
obligation is to connect arbitrary LRC14 rows to either:

1. a bounded/scale-normalized Node-2 atlas with positive margin; or
2. a Node-3 committed-speed branch where the seed component budget is small
   enough, or the exact-period Weyl error is sharper than `2c/v`.

This fits the HYP-2901/HYP-2902 split: finite AP/Legendre-Venn structure closes
bounded seeds, while exact-period equidistribution closes committed speeds
beyond the lcm denominator wall.
