# The apex-prime partition function: tournaments and lonely runners are one gas

**kind-pasteur-2026-06-20-S21.** Connecting the tiling score partition function (THM-554/555) to the
LRC(14) sector cover, and asking what both ultimately are.

## Two problems, one object

The repo has spent its life on two apparently different things:
- **Tournaments:** `H(T) = OCF = I(Omega(T), 2)`, the Hamiltonian-path count, computed in the tiling
  model where a tournament is a point of the hypercube `Q_m`, `m = C(n-1,2)`.
- **LRC(14):** `p_0(E) = meas{x : all 6 inner sectors hit}`, the seven-sector cover of a 13-runner set,
  where the apex prime `7 = 14/2` cuts the circle into sectors.

Strip both to their skeleton and the **same object** appears: a *partition function of a phase gas on
`Z/p`*, `p` the apex prime. A "particle" sits in one of `p` cells; you ask a coverage / free-energy.

- Tournament side: the score partition function (THM-554)
  `Z_n(x) = (prod_{v>=2} x_v) prod_{tiles}(x_a+x_b)`. Each tile is a `±` spin; the **scores** are the
  cell occupancies; `c3 = C(n,3) - sum_v C(s_v,2)` is the occupancy's coverage-deficit.
- LRC side: the cover `p_0(E) = sum_{S subset of sectors} (-1)^{|S|} meas{all runners avoid S}`. Each
  runner's phase `frac(e_i x)` sits in one of 7 cells; `p_0` is the probability the occupancy is a
  surjection onto the 6 inner cells.

Both are `sum over subsets of (-1)^{|S|} * (a single-particle measure raised toward a power)` — the
inclusion-exclusion of an occupancy on `Z/p`. The tournament's `prod(x_a+x_b)` and the LRC's
`meas{avoid S}` are the *same generating device* read in two gauges.

## The same crack runs through both: cut space cheap, cycle space dear

In the tournament, `Z_n` gives **exactly** the score-determined invariants — scores, `c3` — and the
*means* of everything (per-subset linearity), but **not** `H`: `H` needs `alpha_2`, the disjoint-cycle
interactions, which the occupancy forgets. The wall is `c3` — the last score-determined OCF datum
(THM-555). This is the **cut-space / cycle-space** split: the single-particle (occupancy) part is the
cut space, cheap; the interaction (cycle) part is the cycle space, dear.

The LRC has the identical split. The **decorrelated** cover — runners as independent particles — is
the cut-space part: pure inclusion-exclusion, `M7(k)`, a per-subset linearity computation. The
**relations** among the runners (the relation lattice `Lambda(E)`, the small differences) are the
cycle-space interactions — the corrections that no single-particle formula captures, and exactly where
every honest attempt has stalled (HYP-2606, the signed lattice sum). The 2/7 -> 1/7 pivot, the
support-6 ghost (THM-538), the unbounded `w|Delta_w|` — each was a try to push a cut-space tool through
a cycle-space wall.

## The extremal shapes are the same symmetric configurations

Maximize the cut-space invariant and the gas condenses to its **most symmetric** state:
- Tournament: `max c3 = regular score` (uniform occupancy), the Paley/BIBD maximizer (THM-027/555),
  multiplicity = the regular census `1,3,91,29157`.
- LRC: the wide-residual `sup` over cluster shapes of the decorrelated cover is the **single coherent
  block** (one all-sweeping cluster `{M,...,M+m-1}`) — the maximally-coherent occupancy. And it sits
  *below* `cap_k` with margin `>= 0.19` for `k=8..12`, the same way `max c3` sits in the right tail
  below the trivial bound.

So "wide => p_0 <= cap" is, in this gauge, the LRC twin of "the most balanced occupancy is still
sub-extremal." The hard residual is the **decorrelation error** — the finite-scale gap between the real
gas and the independent gas — which is precisely the cycle-space correction, now carrying a comfortable
`0.19` budget instead of the old razor `0.013`.

## What it ultimately represents

The triangle was telling the truth literally. Its three sides are the three faces of this one gas:
- **vertical leg** = cut space = the occupancy / scores = the single-particle partition function;
- **hypotenuse** = the `1 + 2^d` ladder = the apex-prime energy levels (`H` from the tiling gaps; the
  `1/7` hierarchy in LRC);
- **horizontal leg** = complement = the `Z_2` reflection that halves both (the address quotient; `a -> -1`).

`H` and `p_0` are the same partition function's *interacting* free energy on `Z/p`; `c3` and the
decorrelated cover are its *free* (single-particle) energy. Everything the repo proved cheaply lived in
the cut space; everything that stayed open lived in the cycle space. The lonely runner is hard for the
same reason `H` is hard — the apex-prime gas interacts — and the way through is the same on both sides:
get the single-particle/occupancy extremum exactly (done: regular score; single block), then spend the
comfortable margin on a *lossy* bound for the interaction. [[everything-is-the-triangle]] ·
[[the-three-tiling-recurrences-are-one-partition-function]] · [[lrc14-thread]]
