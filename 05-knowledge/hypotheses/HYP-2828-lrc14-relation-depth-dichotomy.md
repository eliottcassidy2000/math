---
id: HYP-2828
title: LRC14 relation-depth dichotomy for the genuine-wide seam
status: OPEN theorem target; exact finite evidence only
source: codex-2026-06-22-S80
depends_on:
  - HYP-2822
  - HYP-2820
  - HYP-2817
  - HYP-2818
  - HYP-2819
  - HYP-2810
  - HYP-2807
  - HYP-2808
  - THM-564
  - THM-563
related:
  - HYP-2823
  - HYP-2816
  - HYP-2815
  - HYP-2814
  - HYP-2701
  - HYP-2637
  - HYP-2606
  - OPEN-Q-108
---

# HYP-2828: Relation-Depth Dichotomy For The Genuine-Wide Seam

## Claim

The failed one-peel dichotomy should be replaced by a relation-depth dichotomy.
For a normalized primitive genuine-wide LRC14 row `E`, define its peel depth as
the least number of deleted elements after which the affine-normalized row has
span at most `14`.  The new theorem target is:

```text
If p0(E) > Q(k-1), then E has peel depth <= 2.

More sharply: in the k=12 genuine-wide branch, every p0>Q(11) row belongs
to a two-peel bounded atlas, hence to the generalized-doublet/R-tail lane.
Rows of peel depth >= 3 should satisfy p0(E) <= Q(k-1), and therefore are
cap-safe with the old slack currency.
```

This is not a proof.  It is the first exact finite audit that turns the
mac-mini S23 covolume/dichotomy prompt into a concrete low-depth separator.
It also refines the older HYP-2606 covolume route: absolute saturated covolume
alone was too lossy, but low covolume after two deletions appears to be the
right structural coordinate.

## Exact S80 Evidence

Script:

```text
04-computation/lrc14_covolume_dichotomy_bridge_codex_s80.py
```

Stored output:

```text
05-knowledge/results/lrc14_covolume_dichotomy_bridge_codex_s80.out
```

The exact audit scans the same normalized primitive genuine-wide window used by
S78 for `span<=18`, keeping exact `p0` arithmetic.  It records full, one-peel,
and two-peel affine-normalized squared norms, plus the least peel depth needed
to reduce span to `<=14`.

```text
k=10, span<=18:
  total=48620, primitive=48619, genuine_wide=27484, over_Q=0.
  Depth histogram all: d1:1, d2:23353, d3:4075, d4:55.
  Pearson p0-Q against inverse sqrt norm:
    full +0.1236, one-peel +0.1646, two-peel +0.2736.

k=11, span<=18:
  total=43758, primitive=43758, genuine_wide=29724, over_Q=0.
  Depth histogram all: d2:23110, d3:6449, d4:165.
  Pearson p0-Q against inverse sqrt norm:
    full +0.1625, one-peel +0.2109, two-peel +0.3307.

k=12, span<=18:
  total=31824, primitive=31824, genuine_wide=24816, over_Q=4.
  Depth histogram all: d2:17097, d3:7389, d4:330.
  Depth histogram over_Q: d2:4.
  Over_Q two-peel span: le14:4.
  Pearson p0-Q against inverse sqrt norm:
    full +0.2044, one-peel +0.2639, two-peel +0.3630.
```

Thus every exact positive `p0-Q` row in the `span<=18` scan has peel depth
`2`, and after two deletions affine-reduces to span `<=14`.  No depth-3 or
depth-4 row beats the decorrelated floor in this window.

The same pattern holds on the seven positive k=12 witnesses from the S78
`span<=20` exact audit.  The S80 script includes a targeted witness-bank mode
so this can be checked without rerunning the larger span-20 scan:

```text
S78 span<=20 positive witness bank:
  depth histogram: d2:7.
  two-peel span: le14:7.
```

Representative witnesses:

```text
(0,2,4,6,8,9,10,11,12,14,16,18)
  p0=238949/388080, p0-Q=36613/2716560,
  depth=2, remove (9,11), two-peel span=9.

(0,2,3,4,6,8,9,10,12,14,15,18)
  p0=10697/17640, p0-Q=257/61740,
  depth=2, remove (0,18), two-peel span=13.

(0,2,4,6,7,8,10,11,12,14,18,20)
  p0=117647/194040, p0-Q=919/226380,
  depth=2, remove (7,11), two-peel span=10.
```

The best all-row covolume proxy is not the full normalized norm but the
two-peel norm.  Its correlation with `p0-Q` is the largest at every audited
`k`:

```text
k=10: +0.2736.
k=11: +0.3307.
k=12: +0.3630.
```

This is weak correlation, not a proof, but it is directional evidence that the
danger signal is "nearly bounded after two deletions" rather than "globally
small norm."

## Proposed Proof Split

The relation-depth dichotomy would let the current proof DAG become more
modular:

1. **Depth 1: single-far lane.**  If one deletion affine-reduces to bounded
   span, the row is in the THM-563/HYP-2820 endpoint-period world.  HYP-2822
   supplies the needed correction that small `f>=15` q6 contraction can be as
   large as `14/15`, so the proof must use a boundary envelope and exact margin
   ledger rather than a uniform `1/7` story.
2. **Depth 2: generalized-doublet lane.**  If two deletions affine-reduce to
   bounded span, the row should be routed to the generalized-doublet atlas:
   a bounded base plus a far pair `{M,M+g}`, frozen room, the
   Mordell-Tornheim `12*zeta(3)` R-tail, and a finite low-`M` window.
3. **Depth >= 3: separator lemma.**  Prove directly that depth-3+ genuine-wide
   rows satisfy `p0<=Q(k-1)`, or at least `p0<cap_k` with a uniform slack margin.
   This should be the easiest branch because the S80 audit shows many such
   rows but no floor beaters, and their two-peel span remains large.

The open separator can be phrased in Freiman language:

```text
High p0 above the decorrelated floor forces high additive energy
and a rank-low GAP-like relation lattice.  In this finite LRC14 address
space, "rank-low enough to beat Q" should imply two-peel reducibility.
```

Equivalently, depth-3+ rows are too genuinely spread.  Their relation lattice
has insufficient coherent same-sign packets to create the even-AP plus
odd-bridge resonance seen at k=12.

## Link To gK8 Concentration

HYP-2828 is compatible with the cleaner HYP-2812 gK8 route and the incoming
HYP-2823 variance reframe.  HYP-2823 says concentration extremality is a
second-moment problem for the miss count

```text
N = number of empty inner sectors,
Var(N) = S1 + 2*S2 - S1^2,
q0+q6 = P(N in {0,6}).
```

The later S23 exact-moment sharpening makes this more concrete:

```text
L_yK8 = 10 - 10*S1 + 10*S2 - 9*S3 + 6*S4.
```

Thus the gK8 route can be read as a degree-4 feasible-moment-region problem.
HYP-2828 is an exception atlas for that moment-region proof: if the global
degree-4 inequality needs to discharge resonant rows separately, the audit
says the only audited positive `p0-Q` rows are already depth-2 and therefore
belong to the generalized-doublet/R-tail lane.

In that language, HYP-2828 is not a competing proof of concentration
extremality.  It is the finite resonant taxonomy underneath it: positive
`p0-Q` rows survive only when two deletions reveal a bounded, highly correlated
core.  If the global moment-region proof of HYP-2823 succeeds, this
relation-depth theorem becomes optional for the cap inequality.  If the
moment-region proof has resonant exceptions, HYP-2828 says where to send them.

The relation-depth atlas is also compatible with the latest Leg-C closure:
incoming HYP-2817/claude-opus S3 exhaustively verifies the generalized-doublet
finite window over all bounded bases, gaps `g=1..4`, and `M in [15,50]` with
zero violations.  That makes the depth-2 lane much stronger computationally:
depth 2 can now be treated as the checked generalized-doublet/R-tail atlas,
while the remaining structural novelty in HYP-2828 is the depth-3+ separator.

Thus if gK8 concentration extremality is proved globally, this relation-depth
theorem is not needed for the final cap inequality.  It remains useful as the
finite resonant atlas underneath the concentration proof:

```text
depth 1 explains single-far q6 boundary packets;
depth 2 explains two-far odd-bridge packets;
depth >=3 is the decorrelated middle-mass branch.
```

In that reading, relation depth is an address quotient for the same smoothing
phenomenon seen by `L_yK8=10q0+q3+10q6`: high `q0=p0` can survive only when a
small relation lattice keeps the row close to a bounded core.  Otherwise
deconcentration suppresses the weighted extremes.

## Relation To The Witness-Floor Route

Post-fetch KPS HYP-2825/HYP-2826/HYP-2827 changes the global proof order.  The
sector-cap branch and the witness-floor branch are parallel, not the same
inequality in disguise.  HYP-2826 suggests the last LRC14 link may be the
1/7-scale witness floor

```text
G2(P,E) = meas{x in G_P : maxgap({e_i x}) > 1/7} >= m_P,
m_P = 14249/252252.
```

HYP-2828 should not be mistaken for that final global crux.  It lives inside
the sector/genuine-wide branch: it is useful for explaining why the exact
sector cap evidence localizes to low relation depth, and for routing any
remaining sector-cap exceptions into already-developed finite atlases.  If the
KPS witness-floor route is proved, HYP-2828 remains a proof-order diagnostic
and a guardrail for formalization rather than a necessary final theorem.

## Assumption Challenge

The script explicitly rejects runner-level Tournament Analysis for this task.
The useful vertices are proof buckets:

```text
(peel depth, two-peel span bucket, far count).
```

The observable is exact bucket maximum `p0`; the switch orients toward larger
bucket maximum; the Hamiltonian path is the displayed bucket order.  This
quotient preserves the predicate "can this profile beat Q or cap" and destroys
the exact phase, the exact far placement, and the internal endpoint ledger.
Therefore a positive HYP-2828 proof would still need to hand depth-1 and
depth-2 profiles back to the endpoint-period and R-tail machinery.

Alternate vertex sets considered:

- runners;
- far elements;
- remove-one witnesses;
- two-peel bounded cores;
- affine relation-lattice buckets;
- proof obligations.

The selected quotient is proof obligations.  Runner vertices would miss the
fact that the same two far speeds can be safe or unsafe depending on which two
bridges make the row affine-bounded after deletion.

## Status

HYP-2828 is open.  It contributes:

- an exact finite `span<=18` relation-depth audit at k=10,11,12;
- a targeted audit of the seven known S78 k=12 positive witnesses through
  `span<=20`;
- the observation that every positive `p0-Q` witness in those banks has peel
  depth `2` and two-peel span `<=14`;
- the sharper diagnostic that two-peel normalized norm tracks `p0-Q` better
  than full or one-peel norm in the full exact row banks;
- a new proof obligation: prove the depth-3+ separator and route depth-1/depth-2
  rows to already-developed endpoint-period and R-tail ledgers.

No LRC14 proof is claimed.
