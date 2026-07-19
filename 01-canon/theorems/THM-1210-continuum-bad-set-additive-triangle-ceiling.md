---
id: THM-1210
title: The continuum four-comb ceiling by the non-arithmetic additive-triangle quotient
status: PROVED, COMPUTER-ASSISTED.  The comb/max-gap identity, Hamilton-cycle identity, BAD-to-three-band K4 implication, non-AP triangle deletion, six-triangle torus geometry, and p+q>=26 sheared-grid tail are analytic.  The remaining 99 coprime pairs 1<=p<q, p+q<=25 are exhaustively certified by two exact rational evaluators; their unique reduced maximum is J(1,2)=2/21.  Thus mu(BAD)<=2/21, with equality exactly for four-term arithmetic-progression frequencies.  The finite arithmetic core and equality-obligation lemma also have a Lean kernel certificate; the Haar/alcove bridge remains analytic prose.  This theorem does NOT by itself supply the finite-comb/eroded-start glue needed for uniform r=5.
source: codex-2026-07-18-S77
related: [THM-527, THM-530, THM-1092, THM-1147, THM-1172, THM-1174, THM-1177, THM-1181, THM-1211, MISTAKE-181]
script: 04-computation/lrc14_continuum_bad_hamilton_cycle_referee_codex_S77.py
output: 05-knowledge/results/lrc14_continuum_bad_hamilton_cycle_referee_codex_S77.out
script_sha256: 6379d0b83f045688f7866766425cb1ed1c43b4a3b5fbb0966ed2781851041e7f
output_sha256: c6334b9a31accfb1f63120162191343a3610589dc513d54349acd23da653bc27
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCContinuumTriangleCeiling.lean (kernel-pure arithmetic cutoff, non-AP deletion, exact carry core, and rational ceiling/equality bridge; Haar/alcove connectors remain external)
lean_sha256: fafdd2de6a8b00c92e4cd423e1bae40baff230f1ed2532b9d29f03a95175ea42
---

# THM-1210 -- the additive-triangle ceiling

## 1. Statement

Let `0<d2<d3<d4` be integers, put

```text
P_d(u)={0,{d2 u},{d3 u},{d4 u}} subset R/Z,
```

and let `Delta_max(P_d(u))` be its largest cyclic gap.  In the continuum
four-comb model define `F_infty(u)` to be the longest component left in
`[0,1]` by the three teeth

```text
[(7/6){-di u}-1/6, (7/6){-di u}] intersect [0,1].
```

Then

```text
F_infty(u)=(7/6)(Delta_max(P_d(u))-1/7)_+,                 (1)
BAD_d={u:F_infty(u)<=1/6}
     ={u:Delta_max(P_d(u))<=2/7},                         (2)
mu(BAD_d)<=2/21.                                           (3)
```

The family `d=(m,2m,3m)` attains `2/21`, by common-dilation invariance and
the exact two runs for `(1,2,3)`.  Conversely, equality forces this form.
Thus `2/21` is the sharp **uniform value**, and the equality directions are
exactly the four-term arithmetic progressions `{0,m,2m,3m}`.

All interval boundaries below are finite rational sets and therefore Haar
null.  Closed intervals are used in (1)--(3); open intervals are used when
convenient in the correlation proof without changing any measure.

## 2. The exact comb/max-gap dictionary

Change coordinates from `x in [0,1]` to `y=6x/7 in [0,6/7]`.  A tooth with
right edge `g={-di u}` becomes

```text
[g-1/7,g] intersect [0,6/7].
```

Adjoin the fixed arc `[6/7,1]`, itself of length `1/7`.  We now have four
directed circle arcs of length `1/7`, with right endpoints

```text
{0,{-d2 u},{-d3 u},{-d4 u}}.
```

If their cyclic endpoint gaps are `Delta_0,...,Delta_3`, the complementary
gaps after the arcs have lengths `(Delta_j-1/7)_+`.  Scaling back by `7/6`
proves (1).  Negation preserves the cyclic gap multiset, so the endpoint set
may be written as `P_d(u)`.  Equation (2) follows immediately:

```text
(7/6)(Delta_max-1/7)_+ <= 1/6  iff  Delta_max<=2/7.
```

This is also the precise correction to the informal language that the bad
region consists of boxes or of a tube: in each fixed cyclic-order chamber it
is a simplex in the four gap variables.

## 3. The Hamilton-cycle precursor

The six cyclic orders of the four labeled points, based at label `0`, pair
under reversal into the three undirected Hamilton cycles of `K4`.  Fix one
order

```text
pi=(0,pi_1,pi_2,pi_3),   pi_4=pi_0,
h_i=d_{pi_(i+1)}-d_{pi_i},   d_0=0.
```

Then `sum_i h_i=0`.  Away from collision walls, the four gaps in this order
are `({h_i u})_i`.  For `f=1_[0,2/7]`, the order-cell bad indicator is

```text
product_i f({h_i u}).                                      (4)
```

The converse implicit in (4) is worth spelling out.  If all four fractions
lie in `(0,2/7]`, their sum is an integer because `sum h_i=0`; it lies in
`(0,8/7]`, so it equals `1`.  Their partial sums therefore reconstruct the
stated cyclic order and its four gaps.  Reflection `u -> -u` reverses the
order and preserves measure.  Choosing one orientation of each undirected
Hamilton cycle gives the exact identity

```text
mu(BAD_d)=2 sum_{C in HC(K4)/reversal}
               integral_T product_{e in C} f({h_e u}) du.  (5)
```

Identity (5) exposes the compatibility question, but bounding its three
terms separately loses too much.  The next quotient deletes to one additive
triangle and proves the ceiling directly.

## 4. The three-band edge predicate

Define the symmetric three-band set

```text
A=[1/7,2/7] union [3/7,4/7] union [5/7,6/7] subset T.       (6)
```

**K4 implication.**  If four circle points have every cyclic gap at most
`2/7`, every gap is also at least `1/7`: the other three gaps total at most
`6/7`.  An oriented difference of two points is the sum of one, two, or
three consecutive gaps.  Such sums lie respectively in

```text
[1/7,2/7], [3/7,4/7], [5/7,6/7].                           (7)
```

For the middle case, the two selected gaps give the upper bound `4/7` and
the complementary two gaps give the lower bound `1-4/7=3/7`.  The
three-gap case is the complement of one gap.  Consequently BAD forces all
six pair differences of `P_d(u)` into `A`.

The implication, rather than an equality, is all that is needed.  There are
some boundary configurations whose six differences lie in the closed `A`
but whose max-gap is `3/7`; they form a null wall set and are not silently
used in the argument.

## 5. Delete to a non-arithmetic triangle

Every four distinct integers contain a non-arithmetic three-point subset.
Indeed, for sorted `v0<v1<v2<v3`, use `(v0,v1,v2)` unless its two gaps are
equal; then use `(v1,v2,v3)` unless those gaps are equal.  If both triples
are arithmetic, all three consecutive gaps are equal, and `(v0,v1,v3)` has
successive gaps in ratio `1:2`.

Choose such a triple `a<b<c` and write

```text
p=b-a>0,   q=c-b>0,   p!=q.
```

Translation by `a` does not change pair differences.  Section 4 gives

```text
BAD_d subset {u:{pu},{qu},{(p+q)u} in A}.                  (8)
```

Set

```text
J(p,q)=mu{u:{pu},{qu},{(p+q)u} in A}.                      (9)
```

It remains to prove the following two-frequency lemma.

> **Additive-triangle lemma.**  If `p,q` are positive integers and `p!=q`,
> then `J(p,q)<=2/21`.

Common dilation is harmless: multiplication by `g=gcd(p,q)` preserves Haar
measure, so `J(p,q)=J(p/g,q/g)`.  Symmetry permits `1<=p<q` and
`gcd(p,q)=1`.

## 6. Six exact triangles on the two-torus

Let

```text
R={(x,y) in T^2:x in A, y in A, x+y in A}.                 (10)
```

Write `x=(i+alpha)/7`, `y=(j+beta)/7`, with
`i,j in {1,3,5}` and `0<alpha,beta<1`.  The carry in
`floor(7(x+y))` gives exactly six possible cells:

```text
alpha+beta>1 in (i,j)=(1,1),(1,3),(3,1),
alpha+beta<1 in (i,j)=(3,5),(5,3),(5,5).                  (11)
```

For `i+j=6` the parity after the carry is never admissible.  For `i+j<6`
one needs a carry; for `i+j>6` the subtraction of `7` reverses parity and
one needs no carry.  Thus `R` is the disjoint union, up to boundaries, of
six right triangles of leg `1/7`.  In particular

```text
area(R)=6*(1/2)*(1/7)^2=3/49,                              (12)
```

and every horizontal or vertical section of `R` is a union of at most two
circle intervals.  The sharper tail below uses the sheared sections, whose
support has measure only `3/7`.

This is the exact version of the observed six-lump geometry.  They are
triangular alcoves, not boxes, and their full inequalities matter.

## 7. The sheared analytic tail `p+q>=26`

Put `N=p+q` and shear (10) by

```text
(x,y) -> (x,z=x+y).                                        (13)
```

For fixed `z`, let `B_z={x:x in A,z-x in A}`.  Write
`z=(r+gamma)/7`, `0<gamma<1`.  The six triangles (11) give the following
complete list, where the table displays `7 B_z`:

```text
r=1: (3,3+gamma) union (5,5+gamma),
r=3: (1+gamma,2) union (5,5+gamma),
r=5: (1+gamma,2) union (3+gamma,4).                         (14)
```

The section is empty for every other `r`.  Thus `B_z` is supported on
`z in A` and is a union of at most two intervals there.  Its three bandwise
area contributions are respectively `1/49,1/49,1/49`, which also rederives
`area(R)=3/49`.

Now condition on `z={Nu}`.  Its inverse branches are

```text
u=(k+z)/N,   0<=k<N.
```

Because `gcd(p,N)=gcd(p,q)=1`, the values

```text
{p(k+z)/N},   0<=k<N,
```

are a shifted `N`-grid in the `x`-circle.  For one circle interval, its grid
frequency exceeds its length by at most `1/N`; for a union of at most two
intervals the excess is at most `2/N`.  Crucially, (14) says that this error
is integrated only over `A`, of measure `3/7`.  Hence

```text
J(p,q)<=3/49+(3/7)(2/N)=3/49+6/(7N).                       (15)
```

For `N>=26`,

```text
6/(7N)<=3/91<5/147=2/21-3/49.                              (16)
```

At the endpoint `N=26` the tail bound and its strict margin are exactly

```text
3/49+6/(7*26)=60/637,
2/21-60/637=2/1911>0.                                    (16a)
```

Therefore the additive-triangle lemma holds analytically and strictly for
every reduced pair with `p+q>=26`.  The Lean theorem `sum_tail_cutoff`
records the strict inequality (not merely a weak ceiling), and
`sum_tail_cutoff_endpoint` independently normalizes both rationals in
(16a).

## 8. The finite exact core

There is also a closed carry formula, independent of interval merging.  It
is useful both as a second referee and as the smaller formalization target.
Put `N=p+q`, and for `1<=k<N` define

```text
L_k = k/q       if k<=q,
      (N-k)/p   if k>q,
c_k = #{r in (2p-q,4p-q,2p-3q): k=r (mod 7)},             (17)
```

where the residue triple is a multiset.  Then exactly

```text
J(p,q)=2/[7N] sum_(k=1)^(N-1) c_k L_k.                    (18)
```

To derive (18), scale the two-torus coordinates by `7` and write `t=7u`.
The three upper triangles in (11) have integer cells
`S={(1,1),(1,3),(3,1)}`.  Set `s=Nt=n+theta`.  On the upper-carry part let
`P=floor(pt)` and `Q=floor(qt)`; then `P+Q=n-1`.  If
`k=pn mod N`, the allowed `theta`-length is the tent `L_k` in (17).  The
seven lifts `n -> n+lN` move `(P,Q)` by `l(p,q)` modulo `7`.  Their line
invariant is

```text
qP-pQ = p-k  (mod 7).                                     (19)
```

Evaluating (19) on the three points of `S` gives precisely the three
residues in (17).  The upper `t`-measure is `N^-1 sum c_k L_k`; scaling back
from `t` to `u` contributes `1/7`, and reflection supplies the factor `2`.

There remain exactly `99` coprime pairs

```text
1<=p<q,  p+q<=25.                                           (20)
```

This is a complete rational interval computation, not a float grid.  For
each positive integer `d`, the set `{u:{du} in A}` is the sorted union of
the `3d` intervals

```text
[(7k+j)/(7d),(7k+j+1)/(7d)],
0<=k<d,  j in {1,3,5}.                                     (21)
```

The referee intersects the lists for `d=p,q,p+q` by a linear merge and sums
their `Fraction` lengths.  It enumerates precisely the coprime domain (20)
and asserts

```text
number of rows=99,
number of distinct exact values=87,
three largest strata, each attained uniquely:
  J(1,2)=2/21,
  J(1,6)=13/147,
  J(3,10)=8/105.                                          (22)
```

In particular there are no violations and the unique reduced maximizer is
`(1,2)`.  The referee evaluates both (18) and the independent interval
merge (21) on all `99` rows and asserts equality row by row.  It also asserts
the exact tail endpoint and margin in (16a).  No chamber sampling or numerical
tolerance enters (16a)--(22).  A Lean file independently
evaluates the denominator-cleared form of (18) with ordinary kernel `decide`;
it checks the core cardinality, absence of violations, and unique equality
pair.  Its rational bridge proves that the cleared predicate is exactly
`triangleMeasure p q <= 2/21`.  The identification of `triangleMeasure` with
the Haar event (9) remains the analytic derivation above, rather than a
measure-theoretic Lean formalization.  Together, (15)--(22) prove the
additive-triangle lemma, and (8) proves the uniform ceiling (3).

The same referee independently audits all `C(30,3)=4,060` triples through
`d4<=30`, all `12,180` Hamilton-cycle correlations in (5), and every chosen
triangle implication.  It finds exact maximum `2/21` on the ten dilates
`(m,2m,3m)`, `1<=m<=10`, and second stratum `4/105` (19 triples).

## 9. Equality directions

Suppose `mu(BAD_d)=2/21`, and write the four sorted integer frequencies with
successive positive gaps `a,b,c`.  Apply the three-band implication to each
of the four three-point deletions.  Whenever its successive gaps `(r,s)` are
unequal, BAD is contained in its additive event `J(r,s)`.  Equality of the
outer measures and the additive ceiling force `J(r,s)=2/21`; strictness of
the tail and the unique reduced core equality pair then force `r:s` to be
`1:2` or `2:1`.  If `r=s`, the triple was already arithmetic.  Consequently
each of

```text
(a,b), (b,c), (a,b+c), (a+b,c)                             (23)
```

is either equal or has ratio `1:2`.

This tiny obligation system forces `a=b=c`.  From `(a,b)`, if `a=b`, then
`b+c>a` and `(a,b+c)` forces `b+c=2a`, hence `c=a`.  If `b=2a`, then
`b+c>2a` makes `(a,b+c)` impossible.  If `a=2b`, then `(a,b+c)` forces
`c=b` or `c=3b`; the first contradicts `(a+b,c)=(3b,b)`, and the second
contradicts `(b,c)=(b,3b)`.  Thus only `a=b=c` remains.  Conversely those
gaps give `{0,m,2m,3m}` and attain `2/21`.  The Lean certificate includes
this disjunctive lemma as `ratio12_four_triangles_rigid`.

## 10. Corrections to the preceding maximiser sketches

The superseded exact-balance draft identifies the centre of an alcove but
does not control the measure of the whole inequality region.  BAD does not
force four gaps to equal `1/8`; it only forces them into `[1/7,2/7]` with
sum `1`.  Positive BAD measure therefore does not require an edge ratio to
remain exactly `1:2:3` on an interval.

The superseded six-box draft's centre-hitting rigidity also drops the integer shifts in
its congruences.  A decisive counterexample is

```text
d=(1,6,7), u=3/4:
({-u},{-6u},{-7u})=(1/4,1/2,3/4),                          (24)
```

although `(1,6,7)` is not proportional to `(1,2,3)`.  Its exact bad measure
is `5/147`.  Thus neither "non-proportional directions miss the centres" nor
"non-proportional directions have zero sojourn" is a valid route.  The
triangle lemma proves the needed upper bound without either assertion.

## 11. Tournament and assumption challenge

There are three useful vertex choices, and only the last is needed for the
bound.

1. **Phase vertices.**  The four labeled phase points carry the pairwise
   observable "which comes first in the cyclic order based at `0`".  The
   switch is the order wall `(di-dj)u in Z`; ties use the numerical label
   path.  On every open chamber the resulting order tournament is transitive:
   score histogram `(0,1,2,3)`, four singleton SCCs, no directed cycle, and
   one Hamiltonian path.  This tournament remembers order and destroys gap
   size, hence cannot prove (3) alone.
2. **Hamilton-cycle obligations.**  Retaining the four signed increments of
   an order gives (5).  Reversal pairs orientations, but bounding the three
   correlations independently destroys their compatibility.
3. **Pair-difference obligations.**  Give an undirected edge `{i,j}` the bit
   `1_A({(di-dj)u})`.  BAD forces a `K4`; delete to a non-AP triangle whose
   edge labels are the additive circuit `(p,q,p+q)`.  This quotient preserves
   exactly the predicate needed for the upper bound and deliberately
   destroys the unused fourth vertex.  Its six alcoves and two-interval
   sections are the proof-bearing geometry.

The challenged assumption is that tournament vertices must be runners, or
that the proof must sum all three Hamilton cycles.  The useful vertices are
pair obligations; the useful circuit is a single non-arithmetic triangle.

## 12. Honest finite-comb boundary

This theorem closes the **uniform continuum measure inequality**.  It does
not by itself close uniform `r=5` in the original finite comb.  That passage
still needs, in the notation of THM-1162,

1. a uniform finite-to-continuum domination/error estimate when the offsets
   may scale with the carrier; and
2. a lower/address bound for the **eroded start complex** `E_k(P)`, not merely
   the raw safe-set measure `mu(S(P))`.

Those are now cleanly separated from the extremal continuum inequality.
