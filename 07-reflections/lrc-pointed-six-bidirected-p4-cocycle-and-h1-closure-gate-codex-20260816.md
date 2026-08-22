# The pointed six is the arc module of a bidirected Boolean path

**Status: FINITE-EXACT REPRESENTATION SIDECAR TO THE TWO-CURRENT-DIGIT
CANDIDATE; INDEPENDENT AUDIT OF THAT PARENT IS ACCEPTED; LRC(14) remains
OPEN.**  The six pointed carrier lines are not merely analogous to a
four-vertex tournament.  They have an exact intrinsic description as the six
directed arcs of

```text
0 <-> 1 <-> 3 <-> 2,                                 (1)
```

a four-vertex relation with three both-way pairs and three missing pairs.
Equivalently, they are the six edges of the alternating state/root incidence
tree `P7`.  This exact typing explains both why the six-line factorization is
natural and why it is not yet a flux: the static graph has `H1=0`.

The address-conditioned diagonal maps split into six `13 by 13` row-sum-one
kernels of total rank `68/78`.  Arc reversal descends at `143/169` digit
pairs and fails at exactly 26, all on the middle root.  Chamber reflection is
exact on all 1,014 scalar entries.  Thus the direction mark is mostly
redundant but globally load-bearing.

## Inheritance pass

The closest exact mechanism is the five-coordinate response tensor at commit
`b1baa781a`.  It retains

```text
(point=(state,u), r0, r1, s=u-q, relation),           (2)
```

has base/one-digit/two-digit/union carrier ranks all six, and gives unique
diagonal maps

```text
K_(r0,r1) in Mat_6(F),
sum_r1 K_(r0,r1)=I_6.                                (3)
```

The canonical hostile is its failed stationary ansatz: no `K_r1` independent
of `r0` exists.  The closest corrected near miss is the distinction between
amplitude rank `13/12` and carrier rank six.  The least-used sidecar is the
incidence relation between the four source subsets and their marked roots.

The present computation rebuilds the complete parent tensor from its
pre-integration workers.  It does not read a stored matrix bank.  The rebuilt
source, tensor, and diagonal-profile digests are exactly

```text
dd2f4837...d6e097,
f08ced17...b28986,
d1c7e561...7f595.                                    (4)
```

During this import, a stale documentation checksum was found.  The tracked
LF bytes of the parent script at `b1baa781a` have SHA-256

```text
9d1671e0...e3205f2,                                  (5)
```

not the previously printed `bc872773...c279`.  The semantic digest
`38725dc1...f204c` and stored output are unchanged; the index and parent
reflection are corrected together with this sidecar.

## 1. The exact four-state relation

Order the owner-visible states by their physical interval subsets:

| state | bits | source subset |
|---:|---:|---|
| `0` | `00` | `{0}` |
| `1` | `01` | `{0,6}` |
| `3` | `11` | `{6,12}` |
| `2` | `10` | `{12}` |

Two states are related when their source subsets meet.  Every nonempty
intersection is a singleton:

```text
0 intersect 1 = {0},
1 intersect 3 = {6},
3 intersect 2 = {12}.                                (6)
```

The other three unordered pairs have empty intersection.  Marking either
endpoint of a nonempty pair supplies both directions.  Therefore the pairwise
observable has census

```text
(both-way, one-way, missing)=(3,0,3).                 (7)
```

This is the precise “tournament with missing and/or both-way edges” object.
It is a bidirected `P4`, not a tournament and not a complete `K4`.  Choosing
one direction on each both-way pair is an orientation gauge; it does not fill
the three missing pairs.

The point-to-arc bijection is

```text
(0,0)  <->  0 -> 1,
(1,0)  <->  1 -> 0,
(1,6)  <->  1 -> 3,
(3,6)  <->  3 -> 1,
(3,12) <->  3 -> 2,
(2,12) <->  2 -> 3.                                  (8)
```

Thus the number six is structural: it counts the arcs of the three realized
both-way pairs.  Its equality with `binom(4,2)` is accidental.  The six
potential edges of `K4` would have different typing.

## 2. The alternating tree and the Boolean closure edge

Keeping each shared root as a vertex gives

```text
state 0 -- root 0 -- state 1 -- root 6
        -- state 3 -- root 12 -- state 2.             (9)
```

This is `P7`.  Its `7 by 6` boundary matrix has rank six, so

```text
dim H1(P7)=0.                                         (10)
```

Contracting the three root vertices gives the state path in (1).  Its
boundary matrix has rank three and again `H1(P4)=0`.

The state order `0,1,3,2` is a Hamilton path in the Boolean square `Q2`.
The one omitted square edge is

```text
2 -- 0.                                               (11)
```

Adjoining (11) changes `P4` to the cycle `C4`.  Orient the four edges as

```text
0->1, 1->3, 3->2, 2->0.
```

Then the boundary matrix has rank three and

```text
H1(C4;F)=F*(1,1,1,1).                                (12)
```

Equivalently, for a one-cochain `j=(j0,j1,j2,j3)`, the explicit cohomology
coordinate is

```text
Phi(j)=j0+j1+j2+j3,                                  (13)
```

because `ker Phi` is the space of vertex coboundaries.  Equations (12)--(13)
are the smallest local owner-state closure template.  On the present static
path every cochain is exact, so no nonzero `H1` flux can be inferred.

The most plausible semantic source for (11) is a later `U_clock` or complete
address return carrying state `2` back to state `0` without changing temporal
copy.  That edge is not present in the current tensor.

### Relation to the existing D5 theorems

This `C4` is **not** THM-3496's seven-chart cycle.  THM-3496 already proves,
after four explicit markings, the normalized coefficient-changing map

```text
H1_graph(C7;F13) -> H1_et(K((lambda));mu_13),
[g] |-> (sum_i g_i) kappa_lambda.                    (D5-1)
```

It also proves that this Kummer line does not map additively to the
characteristic-zero Hamiltonian response module.  THM-3450 supplies a
different characteristic-zero marked-amplitude isomorphism and proves that
the full Keller germ needs the two ANOVA margin sectors missing from a
doubly-centred source.  Neither theorem identifies the present four owner
states with the seven transported charts.

If a lawful clock eventually supplies (11) **and** an `F13`-valued cochain
`j` on that cycle, there is an explicit marked comparison cospan:

```text
H1(C4;F13) --s4--> F13 <--s7-- H1(C7;F13)
                              --Phi_lambda--> Kummer H1,
s4([j])=j0+j1+j2+j3.                                (D5-2)
```

Both seam maps are isomorphisms, so exponent-one normalization sends `[j]`
to `s4([j]) kappa_lambda`.  This is algebraically unique after the same
orientation and generator markings as THM-3496.  Three realization arrows
remain absent:

1. the physical `U_clock` edge producing `C4`;
2. a coefficient map from the present characteristic-zero/cyclotomic response
   amplitudes (represented here by good reduction in a large split field) to
   an `F13` word-current; and
3. a semantic map from this owner-state cycle to the seven-chart cycle.

Thus (D5-2) sharpens the interface without reversing HYP-9031's direct-map
no-go.  Even if all three arrows were supplied, THM-3496's additive-flux
hostile and THM-3450's order-eight/order-fourteen margin hostiles would still
have to be paid.

## 3. The diagonal bundle is six address kernels

Write

```text
K_(r0,r1)=diag(k_e(r0,r1))_(e in Arcs(P4)).           (14)
```

For each fixed arc `e`, assemble

```text
P_e=(k_e(r0,r1))_(r0,r1 in F_13).                    (15)
```

Equation (3) becomes the exact row partition

```text
P_e * 1 = 1                                           (16)
```

for every arc.  The six kernel records are:

| pointed arc | rank | right/left `1`-nullity | support SCCs | row supports |
|---|---:|---:|---:|---|
| `(0,0): 0->1` | 11 | `1/1` | 4 | twelve 10s, one 12 |
| `(1,0): 1->0` | 11 | `1/1` | 4 | twelve 10s, one 12 |
| `(1,6): 1->3` | 12 | `1/1` | 3 | eleven 10s, two 11s |
| `(3,6): 3->1` | 12 | `1/1` | 3 | eleven 10s, two 11s |
| `(3,12): 3->2` | 11 | `1/1` | 4 | twelve 10s, one 12 |
| `(2,12): 2->3` | 11 | `1/1` | 4 | twelve 10s, one 12 |

Their direct-sum rank is

```text
11+11+12+12+11+11=68 of 78.                          (17)
```

Every kernel has thirteen distinct rows.  No kernel has a column sum equal
to one.  Hence (16) is a one-sided cylinder partition, not a doubly stochastic
law.  The entries live in a finite split field and have no proved positivity.

There is nevertheless a useful finite-state formulation.  Algebraically,
each `P_e` is a weighted transition matrix on the 13 possible previous
digits.  The six together form a diagonal bundle over the order-two de
Bruijn address graph.  This is a static two-block code with memory `r0`, not
a chronological Markov process.  A genuine recurrence requires the next
digit to compose with these kernels in the actual three-digit tensor.

## 4. Orientation is a localized, load-bearing sidecar

Arc reversal pairs the six lines as

```text
(0,1), (2,3), (4,5).                                 (18)
```

A diagonal map in (14) descends to the three unoriented edges exactly when
opposite arc weights agree.  The exact census is

```text
962/1014 scalar equalities,
143/169 complete address-pair descents.               (19)
```

The outer-root pairs `(0,1)` and `(4,5)` agree everywhere.  Every failure in
(19) lies on the middle-root pair `(2,3)`.  The 26 failed digit pairs split
canonically:

```text
20 source-middle failures:
  r0 in {3,9}, r1 in F_13\{0,6,12};

6 endpoint residuals:
  (0,11),(0,12),(6,5),(6,7),(12,0),(12,1).            (20)
```

The first 20 coincide with the two middle pointed fibres in the parent's 60
pre-endpoint proportionality exceptions.  The last six survive even where
source proportionality itself is available, isolating an endpoint-level
orientation defect.

The distinct chamber reflection

```text
(point,r0,r1) -> (5-point,12-r0,12-r1)                (21)
```

holds on all 1,014 entries.  Thus orientation is not arbitrary noise: its
failure set is exactly reflection-stable.  Any D5 quotient that identifies
opposite arcs must either retain this 26-pair defect as a sidecar or prove
that the eventual closure/current annihilates it.

## 5. The boundary defect has a binomial, not Fibonacci, profile

Along the `P7` edge order, the six exceptional `r0` fibres are

```text
(12,12,9,3,0,0).                                     (22)
```

They satisfy reflection

```text
r_i+r_(5-i)=12.                                       (23)
```

Their consecutive drops are

```text
(0,3,6,3,0),                                         (24)
```

Arc reversal also splits the location section exactly as

```text
h_even=(12,12,6,6,0,0),
h_odd =(0,0,3,-3,0,0)=3*(e_(1->3)-e_(3->1)).         (24a)
```

The [independent typing audit](lrc-r5-tent-location-c4-cospan-hostile-audit-codex-20260816.md)
records the resulting formal closure cospan.  The normalized odd-basis
convention gives the conditional cycle `3*(1,1,1,1)` and seam `12=-1`, while
literal directed-arc chains introduce a factor two.  More importantly,
`h_odd` records exception locations, not response amplitudes; the missing
`2--0` edge and the location-to-current coefficient map both remain absent.

The drop generating polynomial is

```text
3*x*(1+x)^2.                                          (25)
```

This is a cumulative Pascal-row defect: the active drop kernel is
`3*(1,2,1)`.  It is an exact recurrence-shaped signal, but it is not a
Fibonacci law and not evidence for a ternary branch tree.  The cheapest
meaningful recurrence test is whether a third address digit transports or
convolves (25) while preserving the six-line carrier.

The three roots explain the alternating path and its two interval ends; they
do not by themselves make the ancestry tree ternary.  The retained address
digit is 13-ary.  These two arities belong to different coordinates.

## 6. Subsets of the harmonic series: what is retained and what is lost

Every address-support subset `S` of `{0,...,12}` can of course be represented
by the corresponding subseries

```text
h(S)=sum_(r in S) 1/(r+1).                            (26)
```

The five support subsets actually appearing in the source profiles have five
distinct values of `h`.  On this tiny realized family the scalar happens to
separate support.

It does not separate arbitrary subsets.  Among all `2^13=8192` subsets,
there are only

```text
3712 distinct reciprocal sums,
2944 collision values,
maximum collision multiplicity 3.                    (27)
```

The smallest hostile is

```text
{2} and {3,6}:       1/2 = 1/3+1/6.                  (28)
```

Thus “a subset of the naturals is a subset of the harmonic series” is faithful
when the indexed terms, equivalently the indicator word, are retained.  It is
false as an injective scalar encoding.  The lawful analogue in the LRC tensor
is the full cylinder indicator or its complete Fourier polynomial, not only
one reciprocal sum.

## Connection contract

| field | exact answer |
|---|---|
| source | rebuilt two-current-digit pointed/root-difference response tensor |
| vertices | four owner interval states `0,1,3,2` |
| pairwise observable | nonempty intersection, marked by its unique shared root |
| generalized tournament | bidirected `P4`: 3 both-way, 0 one-way, 3 missing pairs |
| six-object carrier | directed arcs of `P4`, equivalently edges of alternating `P7` |
| address action | six row-sum-one `13 by 13` finite-field kernels, direct-sum rank 68 |
| exact symmetry | chamber reflection on `1014/1014` entries |
| quotient loss | arc reversal fails at 26 middle-root address pairs |
| static homology | `H1(P7)=H1(P4)=0` |
| needed sidecar | lawful closure edge `2--0`, exact address/clock and temporal copy |
| cheapest D5 test | realize the closure over `F13`, then feed its seam through THM-3496's marked cospan; additive flux still faces the proved no-go |
| harmonic boundary | realized five supports separate; arbitrary subset sums collide |
| scope | finite representation sidecar only; no chronology, physical current, D5 bridge, row exclusion, or LRC(14) |

## Next decisive tests

1. **Completed:** the disjoint five-coordinate audit reconstructed the tensor
   and accepted all 169 diagonal maps; use its pinned artifacts as the fixed
   base for the remaining tests.
2. **Completed:** the
   [third-digit probe](lrc-r5-third-current-digit-p4-live-line-obstruction-codex-20260816.md)
   keeps all 2,197 children on the live arc lines, but full-map uniqueness is
   only `130/169`; the pair law is vacuous-only and every fixed child raises
   the 78-state union rank from 68 to 78.
3. Transport the four states to a lawful `U_clock` and ask whether state `2`
   returns to state `0` on the same temporal copy, supplying (11).
4. If the closure exists, project to antisymmetric edge cochains and test
   whether the 26-pair orientation defect is exact, annihilated, or remains as
   a second class.
5. Construct the missing coefficient and chart-realization arrows into
   (D5-2).  Its marked Kummer class is still not additive JC flux; matching
   dimensions or the count six is not a bridge.

## Reproduction

```text
python -B 04-computation/lrc_r5_pointed_bidirected_p4_cocycle_sidecar_20260816.py
python -B -O 04-computation/lrc_r5_pointed_bidirected_p4_cocycle_sidecar_20260816.py
```

Normal and optimized transcripts are byte-identical.  The semantic SHA-256 is
`b2ba313f88fbab0d36e95a63ade832492743ed6007fa0480434083d7dab0ecd3`.
