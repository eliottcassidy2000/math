# LRC(14): q=6 is the binary--ternary CRT fibre product

**Date:** 2026-08-14  
**Status:** FINITE-EXACT INDEPENDENT CONVERGENCE AUDIT; SUPERSEDED AS CURRENT
STRUCTURAL SOURCE BY THM-3395.  This note retains the CRT, tournament,
ancestry, and harmonic interpretations of the exact q=6 computation.  It
gives no refined-ledger decrement and does not prove LRC(14).

## 1. Inheritance board

The closest proved mechanisms are now sharply separated.

| lane | proved mechanism | hostile or corrected near miss | q=6 use |
|---|---|---|---|
| q=2 | THM-3387's edge law `u+v>7 gcd(u,v)` | MISTAKE-382's endpoint-only loss | parity-triple quotient |
| q=3 | THM-3388's affine triangle closure | `(1,4,7)` passes every pair test | antipodal-pair quotient |
| q=4 | THM-3389's typed complete cochain | `(2,7,11)` and `(1,3,11,5)` | template for mixed block sizes |
| ramified sheets | THM-3385's exact `gcd(u,q)` capacity | capacity equality is not cover | blocker species |
| ancestry | word tree versus exponent support versus harmonic mass | word collisions under commuting dilation | q=6 ternary orbit |

The least-used sidecar is the Chinese-remainder coordinate of a blocked sheet.
It turns the apparently heterogeneous q=6 blockers into partial assignments on
one binary and one ternary coordinate.

## 2. Six sheets form a `2 x 3` grid

Use the canonical Chinese-remainder chart

```text
Z/6Z  -->  Z/2Z x Z/3Z,
k      |-> (k mod 2, k mod 3).                         (1)
```

A transverse speed `u` blocks one coset of the kernel of multiplication by
`u` on `Z/6Z`.  Because `6` does not divide `u`, there are exactly three
species:

```text
gcd(u,6)=1:  one point of the 2 x 3 grid;
gcd(u,6)=2:  one whole binary fibre over a ternary point (size 2);
gcd(u,6)=3:  one whole ternary fibre over a binary point (size 3). (2)
```

Thus q=6 is not just six unrelated sheets.  It is the fibre product of the
binary and ternary sheet problems.  The binary coordinate can legitimately be
read as an XOR bit.  The ternary coordinate retains the phase triangle which
THM-3388 showed cannot be reconstructed from pairwise feasibility.

This also separates two superficially similar four-state objects.  A directed
pair with missing and bidirected possibilities has state set
`{none,left,right,both} ~= Z/2Z x Z/2Z`; the q=4 sheet circle is `Z/4Z`.
They have the same cardinality but different subgroup geometry.  Forgetting
that distinction destroys exactly the antipodal-pair information used by
THM-3389.

## 3. Why there are exactly 23 minimal sheet covers

In the `2 x 3` grid, a size-three blocker is a row, a size-two blocker is a
column, and a singleton is a cell.  An irredundant cover can choose `r=0,1,2`
rows and `c=0,1,2,3` columns.  Every cell outside the selected rows and columns
must then be supplied by its singleton; any other singleton would be
redundant.  If both rows or all three columns are chosen, minimality forbids
choosing anything else.  The possibilities are therefore

```text
six cells                                      1
one column + four cells                        3
two columns + two cells                        3
three columns                                  1
one row + three cells                          2
one row + one column + two cells               6
one row + two columns + one cell               6
two rows                                       1
                                               --
                                               23.      (3)
```

This is an analytic classification of the blocker-pattern layer.  It does not
say that arbitrary speeds of those types fire at one common source time.  That
remaining phase-gluing predicate is metric and is carried by an affine
cochain.

## 4. The complete affine cochain

Assign each selected blocker a speed `u_i` and a sheet representative `k_i`.
For a possible common firing phase write

```text
x_i=a_i/u_i-k_i/6,                  a_i in Z,
p_ij=6u_i u_j(x_i-x_j).                              (4)
```

Every oriented pair must satisfy

```text
p_ij == (k_j-k_i)u_i u_j  (mod 6 gcd(u_i,u_j)),
14|p_ij| < 6(u_i+u_j),
p_ji=-p_ij.                                           (5)
```

The normalized gaps `p_ij/(6u_i u_j)` come from vertex centres exactly when
their circulation vanishes on every triangle:

```text
u_h p_ij + u_i p_jh + u_j p_hi = 0.                  (6)
```

Necessity is immediate.  Conversely, `(6)` makes the complete rational
one-cochain a coboundary.  The congruences in `(5)` are precisely pairwise
compatibility of the common-shift congruences for the centres, hence the
generalized CRT supplies one common shift.  The strict inequalities give
pairwise-intersecting lifted real intervals.  There are at most six owners and
their total length is at most `6/7<1`, so the circular wraparound obstruction
is impossible and one source time lies in every interval.

The exact probe compares this criterion with a separate endpoint/mid-cell
event sweep.  They agree on every literal subset.  THM-3395 now proves the
general q<=7 criterion and independently replays the q=6 slice; this probe is
retained as a convergence artifact.

## 5. The quotient tower is exact, but incomplete on mixed patterns

The CRT chart makes two inherited strata literal.

- Two size-three blockers cover the two binary rows.  Writing their speeds as
  `3a,3b`, collapse the free ternary coordinate.  Common phase is exactly the
  q=2 law for `a,b`.  The probe checks `528` such pairs.
- Three size-two blockers cover the three ternary columns.  Writing their
  speeds as `2a,2b,2c`, collapse the free binary coordinate.  Common phase is
  exactly THM-3388's q=3 cochain for `a,b,c`.  The probe checks `1,140` such
  triples.  In particular `(2,8,14)` is the lift of hostile `(1,4,7)`: every
  pair gap exists, but no phase triangle closes.

These quotient maps preserve the full-cover predicate only on their pure
blocker strata.  A mixed pattern such as `3+2+1` uses both CRT coordinates at
once.  Separate q=2 and q=3 shadows no longer remember which row, column, and
cell belong to the same source phase.  The complete q=6 cochain is the needed
gluing sidecar.

This is the precise sense in which the obstruction is a pullback rather than
a product of two yes/no tests:

```text
q=6 cover witness
   --> q=2 binary shadow
   --> q=3 ternary shadow,

but compatible shadows + shared integral phase data --> q=6 witness. (7)
```

## 6. Literal clutter and atlas

For

```text
V={1,2,3,4,5,7,8,9,10,11,13,14},                     (8)
```

the independent event and cochain routes find the same `39`
inclusion-minimal cover edges:

```text
rank 3: 3,             rank 4: 29,             rank 5: 7. (9)
```

The three rank-three edges are exactly pure q=3 lifts:

```text
{2,8,10}, {2,10,14}, {4,10,14}.                       (10)
```

There is no rank-two edge in the literal pool: the only parity-triple speeds
are `3,9`, whose q=2 reduction `{1,3}` fails the gcd threshold.  The mixed
rank-four and rank-five edges are new q=6 information.

The independence profile is

```text
(I_0,...,I_12)=(1,12,66,217,441,515,304,76,5,0,0,0,0). (11)
```

A q=6 body row has either one or two core clocks, so the globally transverse-
safe count is

```text
2I_5+I_4=2*515+441=1471.                              (12)
```

Exactly seven additional rows are rescued by their core danger sets, giving
`1,478` pointwise-exact rows among `2,079` candidates, exactly the q=6 slice
of the THM-3387 atlas.  As in the q=2 endpoint repair, this seven-row difference
is why a global hypergraph is not by itself the pointwise theorem.

## 7. What survives of the tournament idea

There is an intrinsic complete graph in the phase criterion: its vertices are
selected blocker owners and its edges carry the signed integers `p_ij`.
Antisymmetry orients a chosen nonzero gap, but the orientation is only the sign
of the data.  The magnitude, affine residue class, owner type, and triangle
circulation remain indispensable.

Consequently a tournament can be a serialization or audit sidecar:

```text
vertex: selected blocker owner;
pair observable: sign(p_ij), with zero retained as a tie;
gauge: global reflection reverses every sign;
preserved: coarse cyclic ordering of centres;
lost: |p_ij|, congruence fibre, block size, and H^1 closure;
target: common source phase;
needed sidecar: full typed integral cochain.             (13)
```

The decisive hostile `(2,8,14)` has every pair feasible.  Any tournament made
by choosing pair orientations exists, while the q=6 cover witness does not.
This is a direct no-go for treating a tournament certificate as an iff.

## 8. Ternary ancestry, collisions, and subsets of the harmonic series

Common multiplication by any positive integer coprime to six permutes the
sheet labels and preserves blocker species and minimal covers.  Starting from
the edge `{2,8,10}` and the three commuting multipliers `{7,11,13}` therefore
gives a ternary ancestry tree.

At depth `d`:

```text
word addresses:                 3^d;
distinct exponent triples:      binom(d+2,2);
integer support:                {r 7^i 11^j 13^k:
                                 r in {2,8,10}, i+j+k=d}. (14)
```

The word tree is free, but integer multiplication commutes, so many addresses
collide onto one exponent triple.  The support-shell harmonic masses obey the
same elementary-symmetric recurrence as the q=3 and q=4 orbits:

```text
1001H_d=311H_(d-1)-31H_(d-2)+H_(d-3),  d>=3.          (15)
```

The whole orbit has finite harmonic mass

```text
(1/2+1/8+1/10)
 (1-1/7)^(-1)(1-1/11)^(-1)(1-1/13)^(-1)
=29029/28800.                                          (16)
```

This makes the user's “every subset of the naturals is a subset of the
harmonic series” idea precise on a controlled carrier.  Unique prime
factorization identifies the exponent lattice `N^3` with a countable subset
of the positive integers, and `N^3` is countably bijective with `N`.  Hence
every `A subset N` can be transported to a subfamily of this orbit.  Every
such subfamily has a finite weighted harmonic mass because it lies below
`29029/28800`.

The representation warning is essential: a subset of word addresses is not
the same as a subset of integer support, because commuting words collide.
One must either quotient by the Parikh/exponent map or retain multiplicity.
This is exactly the same missing-sidecar issue seen in pairwise tournaments.

The Fibonacci and Berggren trees live nearby but not identically.  Their
distinguished rays are controlled by additive or matrix recurrences, whereas
`(15)` is the characteristic recurrence of a commuting multiplicative monoid
with roots `1/7,1/11,1/13`.  The honest common grammar is

```text
free finite-word ancestry --> quotient support --> weighted boundary measure,
```

not equality of the numerical recurrences.

## 9. Cheapest next decisive tests

1. Independently reconstruct all q=6 edges from event geometry and audit the
   affine-cochain converse outside the literal pool, including strict endpoint
   hostiles.
2. Isolate the seven core rescues analytically as intersections between q=6
   cover cells and the one- or two-clock core danger union.
3. Complete q=5.  Prime five has only singleton blockers, so it is the clean
   rank-five control between typed q=4 and CRT-mixed q=6.
4. State one uniform `2<=q<=6` affine-coset clutter theorem, with q-dependent
   blocker partitions and the complete `K_r` cochain as the common sidecar.
5. Only after that theorem is audited, ask whether its quotient tower can
   transport to the physical-current/drift frontier.  The present probe has no
   such transport map.

Reproduce the finite probe with

```text
python 04-computation/lrc14_q6_typed_cover_clutter_probe_20260814.py
python -O 04-computation/lrc14_q6_typed_cover_clutter_probe_20260814.py
```

Normal and optimized runs agree exactly.  The semantic digest is
`c329b80567995d2ce6110b5e216695b99faa062d04cb70fdfa5bc0db7edc8019`.
