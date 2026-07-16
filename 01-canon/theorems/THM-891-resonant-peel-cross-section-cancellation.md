---
id: THM-891
title: Exact cross-section cancellation for a resonant far peel
status: PROVED. The far-peel limit, exact miss-pattern formula, seven-residue reduction, and consecutive-core constant are rigorous. The diameter-20 extremal statement is finite-exact evidence, not a universal theorem.
source: codex-2026-07-16-S17
depends_on: [THM-727, THM-883, THM-884, THM-887-uniform-maxS-and-affine-witness-coordinate]
related: [THM-888, THM-889, HYP-7021, HYP-7024]
verification: 04-computation/lrc14_resonant_cross_section_cancellation_codex_S17.py -> 05-knowledge/results/lrc14_resonant_cross_section_cancellation_codex_S17.out
---

# THM-891 — exact cross-section cancellation for a resonant far peel

## 1. Statement

Let `E={0,e_1,...,e_5}` be a fixed six-offset core.  For almost every `x`, let
`M_E(x)` be the set of sectors in `Z/7Z` missed by the six points `{e x}`.  The
stationary offset forces `0 notin M_E(x)`.  Write

- `B_s = meas{x : M_E(x)={s}}`,
- `A_{s,c} = meas{x : M_E(x)={s,c}}` for `s<c`.

Append a far offset `t`, let `R_s(E,t)` be the set on which the seven-offset family
misses exactly sector `s`, and use the `THM-727` centered section function

`g_s(y)=1_[s/7,(s+1)/7)({y})-1/7`.

For a fixed positive integer `a`, define the two-scale error at the owner resonance

`Error_t(a)=sum_s integral_[R_s(E,t)] g_s(atx) dx`.

Then

`Error_t(a) = C_a(E) + O_E(1/t)`,

where

`C_a(E) = sum_s B_s sum_(c != s) J_a(c,s)
          + sum_(s<c) A_{s,c}(J_a(c,s)+J_a(s,c))`

and

`J_a(c,s)=integral_[c/7,(c+1)/7] g_s(ay) dy`.

Moreover, if `a=7q+r`, `0<=r<=6`, then

`a J_a(c,s) = (7 delta_r(c,s)-r)/49`,

where `delta_r(c,s)=1` exactly when `s` belongs to the cyclic interval
`{cr,cr+1,...,cr+r-1} mod 7`.  Consequently

> **`a C_a(E)` depends only on `a mod 7`, and it is zero when `7|a`.**

For the consecutive slow core `E={0,1,2,3,4,5}`, the six nonzero values of
`a C_a(E)` are

| `a mod 7` | exact `a C_a(E)` | decimal |
|---:|---:|---:|
| 1 | `-239/5145` | `-0.046452866861` |
| 2 | `209/20580` | `0.010155490768` |
| 3 | `-39/6860` | `-0.005685131195` |
| 4 | `47/20580` | `0.002283770651` |
| 5 | `-4/1715` | `-0.002332361516` |
| 6 | `83/6860` | `0.012099125364` |

Thus every owner resonance of this asymptotic family satisfies

`|C_a(E)| <= 239/(5145 a) <= 239/5145 < 0.097`.

This proves the signed cross-section mechanism behind the exact `THM-884` audit.  The
sectionwise absolute-value constant `0.8287` from the colliding `THM-887` is a valid
uniform envelope, but it discards precisely the signs retained above.

## 2. Conditioning identity

Let `c_t(x)=floor(7{tx})`.  Adding one point can delete at most one missed sector.
Therefore the seven-offset family misses exactly `s` in precisely two disjoint cases:

1. `M_E(x)={s}` and `c_t(x) != s`;
2. `M_E(x)={s,c}` and `c_t(x)=c`.

Ignoring the null boundary set gives the exact identity

`Error_t(a)
 = sum_s integral 1_[M={s}] 1_[c_t != s] g_s(atx) dx
 + sum_(s!=c) integral 1_[M={s,c}] 1_[c_t=c] g_s(atx) dx`.

For an unordered pair `{s,c}`, the two ordered terms are exactly the two summands
`J_a(c,s)` and `J_a(s,c)` in the theorem.

## 3. Averaging lemma and the limit

Each indicator `1_[M_E=T]` is a finite union of rational intervals, because its only
walls are `k/(7e_i)`.  If `F` is any such step function and `h` is a bounded
one-periodic step function, cell decomposition into the `t` intervals
`[j/t,(j+1)/t]` gives

`integral_0^1 F(x)h(tx) dx = (integral F)(integral h) + O_F(1/t)`.

Indeed, cells not meeting a wall of `F` contribute their constant `F` value times
`(integral h)/t` exactly; only `O_F(1)` wall cells remain, each of length `1/t`.

Apply this lemma with

- `h(y)=1_[floor(7{y}) != s] g_s(ay)`, whose mean is `sum_(c!=s)J_a(c,s)`;
- `h(y)=1_[floor(7{y}) = c] g_s(ay)`, whose mean is `J_a(c,s)`.

Substitution into the conditioning identity proves `Error_t(a)=C_a(E)+O_E(1/t)`.

## 4. Exact seven-microcell calculation

Under `z=ay`, the source interval `[c/7,(c+1)/7]` becomes `a` consecutive
microintervals `[k/7,(k+1)/7]`.  Let

`N_a(c,s)=#{k in {ca,...,(c+1)a-1}: k=s mod 7}`.

The target-sector part has measure `N_a(c,s)/(7a)`, while the centered `-1/7` part
integrates to `-1/49`.  Hence

`J_a(c,s)=N_a(c,s)/(7a)-1/49`.

If `a=7q+r`, every residue occurs `q` times and the cyclic block of `r` residues
starting at `cr` occurs once more.  Thus `N_a(c,s)=q+delta_r(c,s)` and

`aJ_a(c,s)=(q+delta_r)/7-(7q+r)/49=(7delta_r-r)/49`.

This proves the seven-residue law.

## 5. Exact consecutive-core calculation

The breakpoint sweep for `E={0,1,2,3,4,5}` gives, in sector order `s=1,...,6`,

`B=(1/35, 1/35, 5/84, 11/210, 1/28, 2/35)`.

The nonzero upper-triangular pair masses are

`A_12=1/35, A_14=13/420, A_16=1/28,
 A_24=2/105, A_25=1/21, A_26=11/210,
 A_35=3/140, A_36=5/84, A_46=1/35, A_56=1/35`.

All omitted `A_sc` are zero.  Substitution into the microcell formula produces the
six-entry table in the statement.  The verification script performs the same sweep on
the integer wall grid and checks every rational exactly.

As a direct finite-`t` referee at `a=1`, the exact errors are

| `t` | `Error_t(1)` |
|---:|---:|
| 25 | `-13/300` |
| 50 | `-689/14700` |
| 100 | `-169/3675` |
| 200 | `-449/9800` |
| 400 | `-391/8400` |

Their differences from `-239/5145` are `O(1/t)` exactly as the proof predicts.

## 6. Finite extremal evidence and the remaining crux

The same exact event sweep was run over all `15,246` primitive cores
`{0,e_1,...,e_5}` with `1<=e_1<...<e_5<=20`.  The unique largest value of
`max_r |r C_r(E)|` is

`16/343 = 0.046647230321`

at `E={0,1,2,3,4,6}`, in residue `r=1`.  This is below the propagation slack but is
not promoted to a theorem outside the finite box.  `HYP-7024` records the sharp
all-core target and two important negative findings: pair-collision mass is arithmetic,
not universally `1/7`, and residue `1` does not dominate pointwise on every core.

## 7. Synthesis with the concurrent routes

- `THM-883` finds the owner comb in each section.  This theorem composes the sections
  before taking absolute values and identifies their exact cancellation.
- `THM-884` observes ratios `0.022--0.060` in full finite audits.  The value
  `239/5145` explains the consecutive-family resonant limit within that band.
- The other `THM-887` proves the coarse `0.8287` sectionwise envelope.  The gap to
  `0.0465` is now an explicit signed miss-pattern identity, not unexplained slack.
- `THM-888` diagonalizes the quadratic comb weights by residue.  `THM-891` is the
  linear, signed cross-section counterpart.
- `HYP-7021` retains short arithmetic relations for balanced owner combs.  Those
  relations are also necessary here: the attempted fixed collision-moment shortcut is
  false precisely on relation-rich pairs such as `(1,8)`.
- The constant-propagation ledger closes the coarse tail at `3.4 diam`.  A universal
  form of `HYP-7024`, plus a finite-`t` remainder bound, would attack the factor-17.6
  loss and could contract that residual band toward the natural `w>=diam` boundary.

## 8. Tournament and assumption challenge

The computation treats the six nonzero residues as vertices.  Its pairwise observable
is `|rC_r|-|qC_q|`; majority over primitive cores supplies the binary relation, and
lexicographic order is the tie gauge.  Both diameter-14 and diameter-20 tournaments are
transitive with Hamiltonian path `(1,2,3,4,6,5)`, score histogram one vertex at every
score `0,...,5`, singleton SCCs, no directed triangles, one Hamiltonian path, and zero
edge flips between the two boxes.

This tournament is useful telemetry, not a proof engine.  Runners, swing arcs, fixed
sections, section boundaries, wall events, residues, Fourier modes, miss patterns, and
proof obligations were all considered as vertices.  The miss patterns are the faithful
carrier: they preserve the exact two-scale coefficient but destroy finite-`t` aliasing
and wall chronology.  Residues preserve only the final error ordering.  Treating runners
or arcs as vertices would lose the signed quotient that creates the cancellation.
