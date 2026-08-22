---
id: THM-3457
title: "Four-point line-metric turnpike preorder atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Four ordered line
  points realize exactly 25 six-distance preorders, with census 10/10/4/1.
  Reversal and S4 give 14 classes; arbitrary labellings give 300.  XOR
  recovers K4 edge incidence, not order.  Fibonacci occupies one wall and
  still needs origin/Pell sidecars.  No LRC, Rule 30, or JC consequence.
source: codex-2026-08-15-four-point-turnpike-atlas
audit: >
  independent algebra, equality-wall, stratum, reversal, S4, XOR-incidence,
  Gram, Fibonacci-boundary, hostile, source-scope, and non-consequence audit;
  exact normal/optimized/stored replay, dependency/hash, AST/security, and
  documentation gates
depends_on:
  - THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order
related:
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
external:
  - Rybin--Zhang--Luo, "XX^t Can Be Faster", arXiv:2505.09814v2 (CITED algorithmic context only)
script: 04-computation/four_point_line_metric_turnpike_preorder_atlas_thm3457.py
output: 05-knowledge/results/four_point_line_metric_turnpike_preorder_atlas_thm3457.out
script_sha256: 1c0f33c69bc3cde0ef9be284b662e7d74e40cc9fbd82d5c35046eb1472e3ed21
output_sha256: 52e0a0000677b81a6a3a63078b36d3419f994a69d543f7219ffc48100b457008
semantic_sha256: a04174188df0edc25bef47c3f67dd5309f28e835519e4efbfa26e3df03c2cb0a
hash_basis: LF-normalized bytes
---

# THM-3457 -- four-point line-metric turnpike preorder atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is an elementary repository classification, not a literature-priority
claim and not a solution of the classical turnpike reconstruction problem.
The word *turnpike* here refers only to the six distances of four labelled
points on a line.

## 1. Typed setup and inheritance gate

Let

```text
x_0<x_1<x_2<x_3,
a=x_1-x_0,       b=x_2-x_1,       c=x_3-x_2,            (1)
```

so `a,b,c>0`.  Write the six labelled distances as

```text
A=01=a,          B=12=b,          C=23=c,
D=02=a+b,        E=13=b+c,        F=03=a+b+c.           (2)
```

Their increasing-value relation is a total preorder on the six edges of
`K4`.  Equal values remain tied.  Its strict part orients every unequal pair
from the smaller distance to the larger one.

The tournament validity card is therefore:

| field | declaration |
|---|---|
| four geometric vertices | the ordered line points `0<1<2<3` |
| six comparison vertices | the labelled edges `01,12,23,02,13,03` |
| pairwise observable | comparison of the two positive distances |
| orientation gauge | smaller distance points to larger distance |
| ties | missing strict edges, or bidirected edges in the weak relation |
| preserved target | the complete labelled distance preorder |
| lost | magnitudes, common scale, translation, reflection gauge, and recurrence origin |
| needed sidecars for Fibonacci recovery | one marked origin coordinate and the Cassini/Pell norm |

The closest proved mechanism is
[THM-3454](THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order.md),
which finds one additive wall of this atlas on a selected Berggren `U`-spine.
The canonical hostile is the gap triple `(1,3,4)`: it lies on the same
one-tie wall as the stable Fibonacci windows but need not be a Fibonacci
window.  [MISTAKE-402](../MISTAKES.md)
is the corrected near miss: spine index and actual rooted depth differ by one.
The least-used sidecars are the positional line gauge and the XOR-detectable
adjacent/opposite incidence of two `K4` edges.

## 2. Complete positional atlas

Positive scaling of `(a,b,c)` changes no signature.  The representatives in
the tables are primitive projective representatives.  Signatures are written
from least to greatest distance.

### 2.1 Ten strict signatures

| ID | sharp chamber | representative | signature | reversal |
|---|---|---:|---|---|
| `S1` | `a<b` and `a+b<c` | `(1,2,4)` | `01<12<02<23<13<03` | `S1r` |
| `S1r` | `c<b` and `b+c<a` | `(4,2,1)` | `23<12<13<01<02<03` | `S1` |
| `S2` | `a<b<c<a+b` | `(2,3,4)` | `01<12<23<02<13<03` | `S2r` |
| `S2r` | `c<b<a<b+c` | `(4,3,2)` | `23<12<01<13<02<03` | `S2` |
| `S3` | `a<c<b` | `(1,3,2)` | `01<23<12<02<13<03` | `S3r` |
| `S3r` | `c<a<b` | `(2,3,1)` | `23<01<12<13<02<03` | `S3` |
| `S4` | `b<a` and `a+b<c` | `(2,1,4)` | `12<01<02<23<13<03` | `S4r` |
| `S4r` | `b<c` and `b+c<a` | `(4,1,2)` | `12<23<13<01<02<03` | `S4` |
| `S5` | `b<a<c<a+b` | `(3,2,4)` | `12<01<23<02<13<03` | `S5r` |
| `S5r` | `b<c<a<b+c` | `(4,2,3)` | `12<23<01<13<02<03` | `S5` |

### 2.2 Ten one-tie signatures

| ID | sharp wall stratum | representative | signature | reversal |
|---|---|---:|---|---|
| `T1` | `c=a+b`, `a<b` | `(1,2,3)` | `01<12<23=02<13<03` | `T1r` |
| `T1r` | `a=b+c`, `c<b` | `(3,2,1)` | `23<12<01=13<02<03` | `T1` |
| `T2` | `c=a+b`, `b<a` | `(2,1,3)` | `12<01<23=02<13<03` | `T2r` |
| `T2r` | `a=b+c`, `b<c` | `(3,1,2)` | `12<23<01=13<02<03` | `T2` |
| `T3` | `a=b`, `c<a` | `(2,2,1)` | `23<01=12<13<02<03` | `T3r` |
| `T3r` | `b=c`, `a<c` | `(1,2,2)` | `01<12=23<02<13<03` | `T3` |
| `T4` | `a=b`, `a<c<2a` | `(2,2,3)` | `01=12<23<02<13<03` | `T4r` |
| `T4r` | `b=c`, `c<a<2c` | `(3,2,2)` | `12=23<01<13<02<03` | `T4` |
| `T5` | `a=b`, `c>2a` | `(1,1,3)` | `01=12<02<23<13<03` | `T5r` |
| `T5r` | `b=c`, `a>2c` | `(3,1,1)` | `12=23<13<01<02<03` | `T5` |

### 2.3 Four two-tie signatures and the maximal tie

| ID | sharp stratum | representative | signature | reversal |
|---|---|---:|---|---|
| `D1` | `a=b`, `c=2a` | `(1,1,2)` | `01=12<23=02<13<03` | `D1r` |
| `D1r` | `b=c`, `a=2c` | `(2,1,1)` | `12=23<01=13<02<03` | `D1` |
| `D2` | `a=c<b` | `(1,2,1)` | `01=23<12<02=13<03` | fixed |
| `D3` | `a=c>b` | `(2,1,2)` | `12<01=23<02=13<03` | fixed |
| `M` | `a=b=c` | `(1,1,1)` | `01=12=23<02=13<03` | fixed |

The last row has two nontrivial tie blocks, of sizes three and two.  It has
four tied unordered edge-pairs, not four tie blocks.

Consequently the positional census is exactly

```text
strict / one-tie / two-tie / maximal = 10 / 10 / 4 / 1. (3)
```

## 3. Why the atlas is complete

Positivity permanently gives

```text
01<02<03,       12<02<03,
12<13<03,       23<13<03.                              (4)
```

Inspecting the remaining pairs in `(2)`, equality can occur only on

```text
a=b,       b=c,       a=c,       a=b+c,       c=a+b.   (5)
```

More precisely,

```text
01=12 iff a=b,              12=23 iff b=c,
01=23 iff a=c iff 02=13,
01=13 iff a=b+c,            23=02 iff c=a+b.            (6)
```

Every other proposed equality forces a positive gap to vanish.  In
particular, the relation `b=a+c` creates no equality among the six distances.

Away from `(5)`, reverse the line if necessary so that `a<c`.  If `c<b`,
there is the single chamber `a<c<b`.  If `b<c`, then either `a<b` or `b<a`,
and independently either `c<a+b` or `c>a+b`.  This gives
`1+2*2=5` strict chambers with `a<c`; reversal gives the ten strict rows.

On `a=b`, the cuts `c=a` and `c=2a` leave three generic pieces.  The wall
`b=c` gives their reversals.  Each additive wall `a=b+c` and `c=a+b` has two
generic pieces according to the order of its two smaller gaps.  These are the
ten one-tie rows.

The wall `a=c` automatically creates both equalities in the middle line of
`(6)` and has the two generic pieces `b<a` and `b>a`.  The remaining positive
wall intersections are the rays `(1,1,2)`, `(2,1,1)`, and `(1,1,1)`.  They
give the last three rows.  Thus no unlisted stratum exists and `(3)` follows.

## 4. Reversal, abstract `S4`, and all labellings

Line reversal is

```text
(a,b,c) -> (c,b,a),
01 <-> 23,       02 <-> 13,       12,03 fixed.          (7)
```

It fixes exactly `D2,D3,M`.  Burnside's lemma therefore gives

```text
(25+3)/2=14.                                             (8)
```

reversal orbits: five strict, five one-tie, three two-tie, and one maximal.

The 25-row atlas uses the positional gauge `0<1<2<3`; it is not itself
closed under arbitrary vertex relabelling.  Nevertheless abstract `S4`
equivalence makes no further identifications.  Indeed, `03` is the unique
largest edge, so an isomorphism between two canonical signatures must
stabilize the unordered endpoint edge `{0,3}`.  Its stabilizer is

```text
{id, (0 3), (1 2), (0 3)(1 2)}.                         (9)
```

The transposition `(1 2)` sends the universal inequality `01<02` to the
impossible `02<01`.  The transposition `(0 3)` sends it to `13<23`, contrary
to the universal `23<13`.  Only the identity and full line reversal survive.
Hence `(8)` is also the number of abstract `S4` classes.

The eleven non-fixed classes have orbit size `24`; the three reversal-fixed
classes have stabilizer two and orbit size `12`.  Therefore the number of
distinct signatures after arbitrary labelling of the four vertices is

```text
11*24+3*12=300.                                          (10)
```

## 5. Missing edges, bidirected edges, and the XOR sidecar

The strict comparison graph on the six edge labels is a transitive `T6`
only in the ten strict chambers.  For a preorder with tie blocks of sizes
`m_1,...,m_s`, it has

```text
15-sum_i binom(m_i,2)                                   (11)
```

strict arcs and `prod_i m_i!` linear extensions.  Equivalently, the weak
comparison digraph bidirects each tied pair.  The four atlas types have:

| type | missing strict comparisons | strict arcs | linear extensions |
|---|---:|---:|---:|
| strict | `0` | `15` | `1` |
| one tied pair | `1` | `14` | `2` |
| two tied pairs | `2` | `13` | `4` |
| maximal `3+2+1` blocks | `4` | `11` | `12` |

There are two distinct objects here.  The geometric vertex order itself is a
transitive `T4`.  The distance comparison object lives on the six edges of
that `K4`.  Cardinality four does not turn the latter into the former.

XOR gives the exact incidence sidecar.  Let `e_0,...,e_3` be the standard
basis of `F_2^4`, and encode edge `ij` by

```text
v_ij=e_i+e_j.                                            (12)
```

For two distinct edges `ij,kl`,

```text
wt(v_ij XOR v_kl)=2  iff the edges share one endpoint,
wt(v_ij XOR v_kl)=4  iff the edges are opposite.        (13)
```

Thus the six labels are the weight-two layer of the Boolean four-cube.  XOR
recovers the octahedral line graph `L(K4)` and its three complementary
opposite pairs.  It is symmetric and supplies no orientation; the numerical
distance comparison is still required.

Equations `(6)` and `(13)` yield a sharp tie-incidence rule:

```text
every one-tie stratum ties adjacent K4 edges;
an opposite-edge tie never occurs alone;
01=23 iff 02=13.                                        (14)
```

This is the precise survivor behind an informal “XOR is a size-four
tournament” analogy: size four supplies the Boolean incidence geometry, but
the tournament or preorder lives on six edge labels and needs a separate
ordered observable.

## 6. Rank-one Gram realization and the `XX^t` boundary

Let `x=(x_0,x_1,x_2,x_3)^t`, let

```text
H=I-(1/4)11^t,       z=Hx,       G=zz^t.                (15)
```

Then `G` is a centred rank-one positive semidefinite `4 x 4` Gram matrix and

```text
(x_i-x_j)^2=G_ii+G_jj-2G_ij.                            (16)
```

Squaring is strictly increasing on positive distances, so `(16)` carries
exactly the same six-edge preorder.  Conversely every rank-one PSD Gram
matrix `zz^t` with four distinct coordinates gives a four-point line metric,
up to common translation and reflection.

Rybin--Zhang--Luo's cited RXTX algorithm uses a `4 x 4` block scheme for a
general product `XX^t`, with eight recursive symmetric calls and 26 general
products; see the maintained
[core-paper role map](../../05-knowledge/reference/CORE-PAPERS-COMPOSITIONAL-RELATIONS.md#rybin--zhang--luo----xxt-can-be-faster).
The common `4+6` indexing is genuine: a symmetric `4 x 4` output has four
diagonal entries and six off-diagonal `K4` edge entries.  But RXTX does not
impose rank one, the line identities `(2)`, or any distance orientation.
It supplies algorithmic span-and-cover context, not this atlas and not an
LRC certificate.

## 7. Fibonacci wall and sharp failure of the converse

For THM-3454's consecutive Fibonacci spine indices at `k>=4`, the three
successive gaps are

```text
(a,b,c)=(F_(k-2),F_(k-1),F_k),
c=a+b,             a<b.                                 (17)
```

Hence every stable Fibonacci window lies in `T1`:

```text
01<12<23=02<13<03.                                      (18)
```

At `k=3`, the gaps are `(1,1,2)` and the additional tie `01=12` gives `D1`.
At `k=2`, the four Fibonacci indices are not strictly increasing, so the
positive-gap atlas does not apply.

The wall `(18)` is necessary but far from sufficient for Fibonacci ancestry.
For example `(a,b,c)=(1,3,4)` has exactly the same signature.  If positioned
at spine index `x_0=1`, it gives `(1,2,5,9)`, for which the marked seam
`b=x_0` fails.  THM-3454 proves that the full recurrence lift requires

```text
b=x_0,       c=a+b,                                     (19)
```

and that the Cassini/Pell unit then selects the consecutive Fibonacci
windows.  In actual rooted-depth coordinates the first equation is affine,
as recorded by MISTAKE-402.  The preorder forgets all of `(19)` except the
single wall equality `c=a+b` and the comparison `a<b`.

## 8. Exact companion and controls

The companion uses only integer and rational arithmetic.  It certifies,
independently of a bounded grid:

- feasibility or infeasibility of all `3^5=243` sign strata of the five
  forms in `(5)`, after the exact normalization `a+b+c=1`: exactly `25`
  are feasible and `218` are infeasible;
- the complete `10/10/4/1` signature atlas and all frozen representatives;
- the three reversal-fixed rows, 14 reversal/`S4` classes, stabilizers, and
  all 300 arbitrary labelled signatures;
- the strict-arc and linear-extension invoice `(11)`;
- the XOR incidence law `(13)`, including the forced opposite-tie pair; and
- the Fibonacci wall, its `k=3` boundary, and the non-Fibonacci hostile
  `(1,3,4)`.

Normal, optimized, and stored replay agree byte for byte after LF
normalization.  The source, output, and semantic hashes are frozen in the
frontmatter.  The companion also freezes THM-3454's LF-normalized dependency
hash `2f551870...7d058`.

## 9. Scope and non-consequences

| field | exact boundary |
|---|---|
| source | four distinct labelled points on the real line |
| target | total preorder on their six labelled distances |
| quotient | common translation and positive scale; reflection only after orbiting |
| preserved | all distance comparisons, tie blocks, and `K4` incidence |
| destroyed | magnitudes, origin, Pell sign, ancestry word, owner, phase, and current |
| cheapest positive | `(a,b,c)=(1,2,3)` |
| cheapest recurrence hostile | `(a,b,c)=(1,3,4)` at `x_0=1` |

The theorem does not reconstruct arbitrary point sets from an unlabelled
distance multiset.  It does not identify XOR with an oriented tournament,
does not classify general rank-one Gram algorithms, and does not transport a
Rule 30 trace or a Berggren ternary branch.  It supplies no common physical
time, owner, phase, LRC current, spectral nonvanishing, Keller map, or
Jacobian degree exclusion.  LRC(14) and `JC(2)` remain open.
