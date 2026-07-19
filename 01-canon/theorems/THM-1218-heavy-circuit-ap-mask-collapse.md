---
id: THM-1218
title: THE HEAVY-CIRCUIT AP-MASK COLLAPSE — mass above 60/637 forces the unique arithmetic-progression completion and at most three beat-mask classes
status: PROVED (THM-1203 analytic provider plus exact 99-pair rational core; deterministic referee; Lean kernel certificate for the strict cutoff, four-obligation rigidity, residue collapse, and three-mask consumer)
source: codex-2026-07-19-S82
depends_on: [THM-1203, THM-1216, THM-1217]
related: [THM-1172, THM-1174, THM-1177, THM-1181, THM-1192, HYP-7845]
script: 04-computation/lrc14_heavy_circuit_ap_mask_collapse_referee_codex_S82.py
output: 05-knowledge/results/lrc14_heavy_circuit_ap_mask_collapse_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCHeavyCircuitAPMaskCollapse.lean
script_sha256: 2a25d558243128fc581e8205295e72bca8a51cd2d9c7bdc7267c8fefc27a23fb
output_sha256: 319418a89f9dd6fbde0e787d1e99d570ed66799b098ebfca8d4c568ad776a53b
formalization_sha256: ad36a51d851a8bcbc458b183f9f748ecccd8a13dc0395e57c50465abd0db0336
---

# THM-1218 — the heavy-circuit AP-mask collapse

## Statement

Let

```text
b1<b2<b3<b4<b5<b6,                 q=b6-b5>0.             (1)
```

For `1<=i<j<=4`, let `BAD_ij` be the continuum four-comb bad event for the
quartet `(bi,bj,b5,b6)`, equivalently

```text
BAD_ij={u in R/Z : the largest cyclic gap of
                    {bi u,bj u,b5 u,b6 u} is at most 2/7}. (2)
```

Put

```text
H=60/637.                                                     (3)
```

If

```text
mu(BAD_ij)>H,                                                 (4)
```

then the quartet is forced to be the unique arithmetic-progression
completion of the defining pair:

```text
(bi,bj,b5,b6)=(b5-2q,b5-q,b5,b5+q).                          (5)
```

Consequently at most one of the six pairs among `b1,b2,b3,b4` is heavy in
the sense of (4).

At a beat numerator `p`, define the direct strict danger mask

```text
D_q(b)={p mod q : 14 min(bp mod q,-bp mod q)<q}.               (6)
```

The four masks belonging to (5) coincide.  In particular, after identifying
the defining masks of `b5,b6`, the usual five beat-mask obligations contain
at most three distinct masks.  This equality persists on the actual
mixed-period master clock

```text
d0=gcd(q,b1,...,b6),                 L=q/d0,                  (7)
```

and `L` is a divisor of `q`.

In the common reduced-period case, with

```text
A(Q)=2ceil(Q/14)-1,
```

the common-zero union estimate therefore improves from five mask classes to
three:

```text
|U| <= 1+3(A(Q)-1)=3A(Q)-2,
B3(Q)=|U|_cap+1=3A(Q)-1=6ceil(Q/14)-4 <= Q.                (8)
```

Thus any `B3(Q)` consecutive beat numerators contain a safe residue.  In the
mixed-period case the corresponding honest statement is: choose one
representative `N_s` from each distinct lifted-mask class, `r<=3`; then

```text
|U| <= 1+sum_(s=1)^r (|N_s|-1).                            (9)
```

As in THM-1217, using this bound to escape a slow gap still requires the
right number of consecutive beat points and a strict cap below `L` (or a
direct proper-union/run certificate).  The present theorem supplies the
class collapse, not those geometric hypotheses.

## The strict additive-triangle spectrum

Let

```text
A=[1/7,2/7] union [3/7,4/7] union [5/7,6/7]
J(p,r)=mu{u:{pu},{ru},{(p+r)u} in A}.                     (10)
```

THM-1203 supplies the following sharp implication:

```text
J(p,r)>60/637  ==>  p=r or p=2r or r=2p.                 (11)
```

Here the strict inequality is essential.  Divide by `gcd(p,r)` and, when
the two reduced entries are unequal, exchange them so `p<r`.  The sheared
section argument of THM-1203 gives

```text
J(p,r) <= 3/49+6/[7(p+r)] <= 60/637       when p+r>=26.  (12)
```

Equality in the second inequality occurs at `p+r=26`, so (12) excludes the
tail only under `J>H`, not under `J>=H`.  The remaining domain has exactly
`99` coprime pairs `1<=p<r`, `p+r<=25`.  Exact interval merging and the
independent carry formula agree on every pair.  The top three reduced
strata are

```text
J(1,2)=2/21,
J(1,6)=13/147,
J(3,10)=8/105.                                           (13)
```

The first is the unique value strictly above `H`; indeed

```text
2/21-H=2/1911,                 H-13/147=11/1911.         (14)
```

The equal-gap case also satisfies (11) directly.  This proves the complete,
unbounded spectral implication (11).  The exact finite core is a theorem
certificate; the referee's additional diameter-24 quartet census and
height-16 sextuple census are diagnostics, not substitutions for (12).

## Four deletion circuits force equal gaps

Fix a quartet in (2) and write its positive consecutive gaps as

```text
x=bj-bi,                  y=b5-bj,                  z=b6-b5=q. (15)
```

The K4 three-band implication in THM-1203 says that `BAD_ij` is contained in
each of the four deletion-triangle events.  Hence

```text
mu(BAD_ij) <= J(x,y),
mu(BAD_ij) <= J(y,z),
mu(BAD_ij) <= J(x,y+z),
mu(BAD_ij) <= J(x+y,z).                                  (16)
```

Under (4), (11) therefore puts all four pairs in the ratio set

```text
R={(s,t):s=t or s=2t or t=2s}.                            (17)
```

The four obligations are rigid.  Normalize `y=1` merely to inspect ratios.
The first two obligations put `(x/y,z/y)` in
`{1,2,1/2}^2`.  Of these nine possibilities, the two aggregate obligations

```text
(x,y+z) in R,                    (x+y,z) in R             (18)
```

hold simultaneously only at `(1,1)`.  Thus `x=y=z=q`, which is exactly
(5).  This is a symbolic integer lemma in Lean; the referee also checks all
`64^3` positive triples through `64`, obtaining precisely the `64` equal-gap
survivors.

Because (5) fixes both values below `b5`, strict ordering permits at most one
index pair `(i,j)`.  Equivalently, the heavy-pair graph on the four
complementary speeds is either empty or a single edge.

## Residue and quotient-clock collapse

Equation (5) gives

```text
bi == bj == b5 == b6  (mod q).                            (19)
```

The predicate in (6) depends only on this residue, proving equality of the
four direct masks pointwise.  For the actual mixed clock, set

```text
gi=gcd(bi,q),       Qi=q/gi.
```

THM-1217 proves `Qi|L`.  Restricting (6) to `p=0,...,L-1` is exactly the lift
of the reduced `Qi`-mask to `Z/LZ`, so the same pointwise equality proves the
master-mask collapse.  No claim that `q|L` is used: the actual direction is

```text
L|q.                                                       (20)
```

The referee realizes every divisor `L|q` for `q<=100` as an actual clock.
For `d=q/L`, it uses the ordered packet

```text
d(1,2,8L+1,9L+1,10L+1,11L+1),                            (21)
```

whose common gcd with `q` is `d`, whose last gap is `q`, and whose last four
entries form the required AP.  It checks the direct predicate against every
reduced quotient lift for all six rows.  It separately replays the AP
predicate on ranges of lengths `2q,3q,5q`; these supplementary pointwise
checks are deliberately **not** called master-clock tests.

Finally, the heavy AP identifies three of the five relevant labels
`b1,...,b5`.  There are at most two other classes, proving the first claim
behind (8)--(9).  Every strict danger mask contains zero, so adjoining the
distinct classes spends at least one overlap per new class.  This proves
(9), and equal size `A(Q)` gives (8).  THM-1216 supplies the elementary
inequality `B3(Q)<=Q` for every `Q>=2` and the consecutive-block consumer.

## Tournament and assumption audit

The first useful graph has the six candidate pairs among
`b1,b2,b3,b4` as vertices.  Its observable is the continuum mass in (2),
switched at `H`; the heavy set has cardinality at most one.  A lexicographic
tie orientation gives the transitive tournament with score histogram
`(0,1,2,3,4,5)`, no directed cycles, six singleton SCCs, and one Hamiltonian
path, but this orientation is only a gauge.

The information-preserving object is instead the four-deletion circuit
hypergraph (16), followed by the equality quotient of the five-mask nerve.
Runner vertices retain neither the Haar mass, the compatibility of the four
additive triangles, nor the resulting mask equality.  Alternative vertices
considered here are pair candidates, consecutive gaps, deletion obligations,
quotient residues, and mask classes.  This challenges the assumption that a
tournament associated to LRC must have runners or arcs as its vertices.

## Formalization and scope

`LRCHeavyCircuitAPMaskCollapse.lean` kernel-checks the exact rational margins,
the weak tail endpoint (hence the need for strict `>`), the denominator-
cleared 99-pair core, four-obligation rigidity from an explicit
`HeavyTriangleProvider`, AP uniqueness, direct and arbitrary-clock pointwise
mask equality, the three-class union ledger, and the threshold consumer.

The Lean provider boundary is intentional.  The identification of
`triangleMeasure` with Haar measure, gcd-dilation invariance, the
BAD-to-four-deletion inclusions (16), and the sheared analytic tail remain
THM-1203 prose providers.  They are represented in the formal theorem by a
provider structure and four domination hypotheses rather than being
silently assumed as kernel facts.

This theorem does not say that some pair must be heavy.  It says that any
pair above the exact `60/637` frontier becomes maximally rigid and saves two
mask classes.  Closing LRC(14) still requires a global argument ensuring
that either such a heavy circuit appears with enough beat points, or the
uniformly non-heavy branch has sufficient transverse overlap or survivor
curvature.
