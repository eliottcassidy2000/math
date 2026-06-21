---
id: THM-562
title: LRC14 L7 bounded-ratio handoff skeleton — generated atom compatibility plus Delsarte packets reduce the remaining balanced two-far window to a finite resonance atlas and one non-resonant 2D discrepancy lemma
status: OPEN proof skeleton with exact scout evidence; HYP-2730 now supplies the torus-discrepancy reduction and an elementary D<=14/p tail proof; not a proof of LRC14
source: codex-2026-06-21-S72
depends_on:
  - THM-534   # sector moment / Delsarte LP dual certificate
  - THM-561   # miss-zeta factorial atom inversion
  - HYP-2728  # generated tail45 strip separator
  - HYP-2726  # Delsarte/Krawtchouk formulation
related:
  - HYP-2729
  - HYP-2727
  - OPEN-Q-108
---

# THM-562 — LRC14 L7 Bounded-Ratio Handoff Skeleton

## Object

Let

`E = B union {f_1,f_2}`, with `B subset {0,...,14}`, `f_1 < f_2`, and
`1 < f_2/f_1 < 43/20`.

For `x in [0,1)`, let `M_E(x)` be the subset of inner sectors `{1,...,6}`
not hit by the sector labels `floor(7 e x) mod 7`, `e in E`.  Define

`p_t(E) = meas{x : |M_E(x)| = t}`.

The L7 row-level cover atom is `p_0(E)`.  The raw row diagnostics
`tail45(E)=p_5(E)+5p_6(E)` and `U4(E)=p_0(E)+tail45(E)` are useful, but
they are not the HYP-2728 normalized generated-atom strip.

## Conditional Statement

The remaining L7 bounded-ratio case follows from the following two premises.

**Finite Resonance Atlas.**  There is a finite set of rational ratio channels
`p/q in (1,43/20)` and a finite residue decomposition in the multiplier
`m` such that for every primitive row with `f_1=qm`, `f_2=pm` in these channels,
the exact cell law gives `p_0(E) <= cap_|E|`.

**Non-Resonant 2D Discrepancy Lemma.**  For ratios not in the finite atlas,
the map `x -> (f_1 x, f_2 x)` has a uniform two-dimensional discrepancy bound
strong enough to place the row-level Delsarte packet below `cap_|E|` after the
HYP-2726 Krawtchouk-positive dual is applied.

HYP-2730 sharpens this premise to the explicit torus-curve inequality

`|R(p/q)| <= D_{p,q} <= 14/p`,

where `D_{p,q}` is the L1 discrepancy of the cell law of
`v -> (qv,pv)` against the uniform `7x7` grid.  The `14/p` proof is
elementary: fix a first-coordinate sector, use `gcd(p,q)=1` to get equally
spaced starts for the second-coordinate arcs, apply Koksma to a trapezoid of
variation `2/7`, and sum the 49 cell errors.

If both premises hold, then KPS L7 closes: the near-merge side is already
routed to the single-cluster theorem, the ratio `>=43/20` side is routed to
the geometric comb chain, and the residual `1 < f_2/f_1 < 43/20` is covered by
the two premises above.

## Exact Scout Evidence

`04-computation/lrc14_l7_tail45_bounded_ratio_codex_s72.py` scanned `3111`
exact rows of this form, with `237` nonprimitive rows skipped and `0` cap
violations.  The named KPS row

`E=(0,2,4,6,8,10,12,25,28)`

has

`p_0=299/1050`, `cap_9-p_0=62911/300300`, `tail45=607/14700`.

The thinnest sampled margins by `k=8..12` are

`1471/5880`, `62911/300300`, `16/91`, `424/1911`, `4891/17640`.

Tournament Analysis on ratio channels is transitive under the risk key
`thin margin, then larger p_0`, with Hamiltonian path

`5/3 > 7/4 > 7/6 > 5/4 > 6/5 > 2/1 > 7/5 > 4/3 > 9/5 > 8/7 > 3/2 > 15/7 > 28/25 > 8/5`.

Raw row `tail45` disagrees with risk on `27/91` channel edges, confirming that
the generated `tail45` strip should not be scalarized directly into the L7 row
bound.  It is a compatibility gate before the Delsarte handoff.

Incoming HYP-2730 then supplies the expected resonance correction:

`R(p/q) = p0_inf(B,p/q) - P2(B)`, with `|R(p/q)| <= D_{p,q} <= 14/p`.

`04-computation/lrc_q108_L7_discrepancy_bound_codex_s72.py` records a cruder
but formalization-friendly backup route:

`D_{p,q} <= 24/(7q)`.

Exact checks through `q<=80` show zero failures for this crude bound and for
the sharper observed `D_{p,q} <= 12/(7q)`.  Since `24/(7*17)<0.21`, every
`q>=17` resonance tail is safe by the crude bound; exact `D` for `5<=q<=16` is
already below `0.21`, so only `2/1,3/2,4/3,5/3,5/4` need row-level atlas
checks.

## Formalization Frontier

The Lean-ready finite layer is the exact cell-law statement: for a finite row
`E`, the breakpoints are `a/(7e)`, the sector word is constant between adjacent
breakpoints, and the resulting `p_t(E)` is a rational sum of interval lengths.
This can be formalized independently of the analytic non-resonant lemma.

The non-finite formal target can now be the elementary seven-bin discrepancy
lemma: condition on `u={qv}`; the second coordinate is a shifted q-point lattice
in seven equal bins; every bin count differs from `q/7` by at most one; summing
the conditional L1 discrepancy gives `D_{p,q} <= 24/(7q)`.

This is weaker than KPS's proved `14/p` bound in some ranges and stronger in
others; both reduce the analytic tail to elementary one-dimensional
discrepancy plus finite checks.

## Honesty

This is not a proof of LRC14.  It is a formal proof skeleton plus exact
bounded-ratio evidence.  HYP-2730 and the S72 discrepancy probe reduce the
remaining mathematical work to wiring the finite row atlas, the finite `f1`
window, and the elementary torus-line discrepancy lemma into the L7 ledger.
HYP-2732 also blocks a tempting overclaim: the sector-cover cap does not
directly imply a positive lonely-measure lower bound.  The theorem remains a
sector-cover/L7 handoff inside OPEN-Q-108's reduction stack.
