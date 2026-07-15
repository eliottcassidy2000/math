---
id: THM-810
title: Four-replacement oriented-deck classification and the order-three coset interface
status: PROVED (all-order scalar classification by symbolic attenuation/capacity reduction plus exact finite Cayley lemmas) + FINITE-EXACT (common-sheet overlap, four primitive base rows, and tournament audit)
source: codex-2026-07-15-S10 quartic continuation
depends_on:
  - THM-804
related: [THM-769, THM-770, THM-800, THM-806, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_four_oriented_deck_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_four_oriented_deck_codex_S10.out
---

# THM-810 — Four-replacement oriented-deck classification

Put `F=F_13^*`.  For labels `r,o in F` and a positive deck order
`D` not divisible by `13`, define the left-germ count and capacity

```text
N_D(r,o)=#{z in Z:-D<z<=D and z=D r o^(-1) (mod 13)},
C_D(r,o)=N_D(r,o)/D.                                      (1)
```

The label `r` is the replacement colour and `o` is the missing-owner family.
The half-open interval in (1) is the exact `(-1/13,1/13]` orientation from
THM-800/804.  In particular, `+1/13` covers a left germ and `-1/13` does not.

## 1. All-order scalar classification

### Theorem 1

Let `R={r_1,r_2,r_3,r_4}` be four distinct elements of `F`, and let
`D_1,...,D_4` be positive integers not divisible by `13`.  Suppose every
owner family is covered at scalar capacity:

```text
sum_(i=1)^4 C_(D_i)(r_i,o)>=1             for every o in R. (2)
```

Then exactly one of the following holds:

```text
(i)  D_1=D_2=D_3=D_4=1;

(ii) D_1=D_2=D_3=D_4=3 and
     R=aH for some a in F,
     H=<5>={1,5,8,12}.                                  (3)
```

Conversely, every row in (i) satisfies (2), and every order-three coset in
(ii) satisfies (2) with equality at all four owners.

Consequently, consider a full-residue four-replacement packet with positive
integers `c,w_r`

```text
A=(c[12] minus {cr:r in R}) union {w_r:r in R},
w_r=cr (mod 13),                    13 does not divide c,
D_r=c/gcd(c,w_r).                                      (4)
```

If `A` is tight at `1/13`, its four replacements cover all oriented
missing-owner germs, so (2) holds.  Thus a tight packet either descends to
common scale, or lies on the exceptional order-three coset interface (3).

The theorem classifies the necessary oriented sheet-capacity predicate.  It
does **not** assert that an order-three exceptional packet is tight.

## 2. Infinite orders reduce to twelve attenuated decks

Write

```text
D=13q+d,                         1<=d<=12.
```

The interval `(-D,D]` is `2q` complete residue blocks plus `(-d,d]`, hence

```text
N_D(r,o)=2q+N_d(r,o),
C_D(r,o)=2/13+(d/D)(C_d(r,o)-2/13).                     (5)
```

Thus every infinite deck is one of twelve residue-deck vectors attenuated
toward the constant vector `2/13`.  This identity is useful, but attenuation
cannot simply be discarded: it improves low entries while reducing high
entries.  The proof therefore uses two exact capacity cutoffs.

Let

```text
f(D)=max_o C_D(1,o),
S_4(D)=sum of the four largest entries of (C_D(1,o))_(o in F). (6)
```

Directly from (5),

```text
D>=3  => f(D)<=1/3, with equality iff D=3,
D>=4  => S_4(D)<=1,   with equality iff D=4.             (7)
```

More explicitly, for `d=1,...,12`,

```text
f(13q+d)  =(2q+m_d)/(13q+d),
S_4(13q+d)=(8q+T_d)/(13q+d),
m=(1,1,1,1,1,1,2,2,2,2,2,2),
T=(1,3,4,4,4,4,5,7,8,8,8,8).                          (7a)
```

These twelve short-deck values prove both sharp equality statements in (7),
rather than leaving them as sampled inequalities.

Also `f(D)<=2/13+1/D` and `S_4(D)<=8/13+4/D`.  Therefore
the two thresholds needed later are finite.  Exact evaluation of (5) gives

```text
L_(1/6)={D>=3:f(D)>=1/6}
 ={3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,20,21,22,
   23,24,27,28,29,30,33,34,35,36,40,41,42,46,47,48,
   53,54,59,60,66,72},                                  (8)

L_(2/3)={D>=4:S_4(D)>=2/3}
 ={4,5,6,7,8,9,10,11,12,15,16,17,18,21,22,23,24,28,
   29,30,34,35,36,42,48}.                               (9)
```

For `D>=79`, the displayed analytic bounds are already strict, so (8)--(9)
are infinite classifications, not arbitrary search cutoffs.

## 3. Proof of Theorem 1

The short decks have especially simple hit alphabets, stated in terms of the
owner/replacement ratio `o/r`:

```text
D=1: capacity 1 at o/r=1 and 0 elsewhere;
D=2: capacity 1/2 at o/r in {1,2,11} and 0 elsewhere;
D=3: capacity 1/3 at o/r in A_3={1,3,5,8,10} and 0 elsewhere. (10)
```

### 3.1 An order-one colour forces common scale

An order-one colour contributes zero at the other three owners.  Those three
owners must therefore be covered by their three remaining colours.  The
three-colour oriented capacity argument in THM-804 forces all three remaining
orders to be one.  This gives (i).

Assume henceforth that no order is one.

### 3.2 Order-two colours

Draw a directed cross-edge `r -> o` between order-two labels when
`o/r in {2,11}`.  Every edge ratio is a signed two.  A directed cycle of
length `m=2,3,4` would require

```text
(+/-2)^m=1 (mod 13),                                    (11)
```

but the possible products are respectively `{4,9}`, `{8,5}`, and `{3,10}`.

If all four orders are two, every owner needs a cross-edge in addition to its
own half-capacity, contradicting (11).  If exactly three orders are two, each
of their three owners still needs a cross-edge: the fourth colour contributes
at most `1/3`, so own `1/2` plus that colour is insufficient.  Again (11)
gives a directed two- or three-cycle.

If exactly one order is two, it hits at most three of the four selected
owners.  At a missed owner the other three capacities must sum to one.  By
(7), all three orders are then three.  The first finite Cayley lemma below
shows that no such labelled `2,3,3,3` row covers all four owners.

Finally suppose exactly two orders are two.  Every owner must receive at
least one order-two hit, because the other two colours total at most `2/3`.
The two order-two colours have at most six incidences on the four selected
owners, so at least two owners receive exactly one half-capacity.  At either
such owner the remaining pair must total at least `1/2`; since each member is
at most `1/3`, both have capacity at least `1/6`.  Their two orders therefore
belong to the finite set (8).  The second finite Cayley lemma rules out all
such rows.

### 3.3 No order-two colour

Let `k` be the number of order-three colours.  If `1<=k<=3`, every order-three
colour must hit every selected owner.  Indeed, if one misses an owner, the
`k-1` remaining order-three colours contribute at most `(k-1)/3`, while each
of the `4-k` other colours contributes strictly less than `1/3`; their total
is strictly less than one.

Two order-three colours hit one another in both directions only when their
ratio is `5` or `8`.  After normalizing them to `1` and `5`, their complete
hit neighborhoods have intersection

```text
A_3 intersect 5A_3={1,5}.                               (12)
```

They therefore cannot both hit two further owners.  Hence a mixed row has
exactly one order-three colour.  Normalize its label to one.  Its other three
owners must be chosen from `{3,5,8,10}`.  They must contribute total capacity
at least `8/3` across the four owner inequalities.  Since each colour has
owner-sum at most one, each of the three remaining colours has total capacity
at least `2/3` across the four owners, placing all three orders in (9).  The
third finite Cayley lemma closes these rows.

If `k=0`, sum the four owner inequalities.  Each colour contributes at most
`S_4(D)<=1`, so equality must hold everywhere.  Formula (7) forces every
order to be four and every selected ratio to be a maximal order-four entry.
The mutually maximal cross-ratio set is

```text
{3,4,9,10}.                                             (13)
```

After one label is normalized to one, its adjacency graph is the four-cycle

```text
3--4--10--9--3,
```

which has no triangle from which to choose the other three labels.  Thus this
case is impossible.

It remains that all four orders are three.  The last finite Cayley lemma says
that the unique normalized four-label set with at least three hits at every
owner is

```text
H={1,5,8,12}.                                           (14)
```

Multiplicative covariance gives every coset `aH`.  Within `H`, each owner is
hit by itself and by the two ratios `5,8`, and is missed only by its negative.
Every owner sum is therefore exactly one.  This proves (ii), its converse,
and Theorem 1.  ∎

### 3.4 Exact finite Cayley lemmas

All normalizations multiply every label by one unit and therefore lose no
rows.  The exact integer replay gives:

| case | normalized label/order bank | rows | largest minimum owner capacity | survivors |
|---|---:|---:|---:|---:|
| one `D=2`, three `D=3` | `C(11,3)` label triples | 165 | `2/3` | 0 |
| two `D=2` | 49 label configurations, two orders in `L_(1/6)` | 78,400 | `17/18` | 0 |
| one `D=3`, three `D>=4` | four label triples, three orders in `L_(2/3)` | 62,500 | `92/105` | 0 |
| four `D=3` | `C(11,3)` label triples | 165 | `1` | `{1,5,8,12}` only |

Every comparison is an integer cross-multiplication.  Floating point and a
sampled circle play no role.

## 4. The coset survivor passes exact sheet overlap

The scalar exception is not an overlap artefact.  Suppose `R=aH` and all
orders are three.  Then

```text
c=3g,                  w_r=g u_r,
gcd(u_r,3)=1,          u_r=3r (mod 13).                 (15)
```

Put `e_r=u_r mod 3 in {1,2}`.  The exact common-sheet cover exists if and only
if

```text
e_r=e_(-r) for the two negative pairs in aH.             (16)
```

Thus, in normalized label order `(1,5,8,12)`, the four and only four feasible
patterns are

```text
(1,1,1,1),    (1,2,2,1),    (2,1,1,2),    (2,2,2,2).   (17)
```

### Proof

Fix an owner `o`, and choose `a_o o=1 (mod 13)`.  For a colour `r` that hits
this owner, let `z in {-2,-1,1,2,3}` be the unique integer with

```text
z=3r/o (mod 13).
```

Its eligible deck sheet `ell in Z/3Z` is obtained by reducing
`u_r(a_o+13ell)=z (mod 39)` modulo three:

```text
ell=e_r z-a_o (mod 3).                                  (18)
```

For the three colours hitting `o`, the ratios `o/r` are `1,5,8`, and the
corresponding `z` residues modulo three are `0,1,-1`.  Their three sheet
classes are therefore, up to the common translation `-a_o`,

```text
0,       e_(8o),       -e_(5o).                         (19)
```

They are all distinct exactly when `e_(8o)=e_(5o)`.  The labels `8o` and
`5o` are negatives, so applying (19) at every owner is precisely (16).  Each
of the three classes repeats `g` times among the `c=3g` sheets; hence (16)
covers every left germ exactly once.  CRT realizes every sign pattern in
(17), proving necessity and sufficiency.  ∎

After dividing (4) by `g`, the exceptional packet is

```text
A/g=3([12] minus aH) union {u_r:r in aH}.                (20)
```

It is therefore an `s=3` deep packet with eight on-sheet speeds and four
off-sheet exceptions.  The first genuine Hamming-four obstruction to
common-scale descent is exactly a shallow/deep gluing interface, not an
isolated shallow anomaly.

For orientation only, take the least positive CRT representatives in (15)
when `a=1`.  Exact breakpoint evaluation gives:

| pattern | `(u_1,u_5,u_8,u_12)` | `M(A/g)` |
|---|---|---:|
| `(1,1,1,1)` | `(16,28,37,10)` | `6/43` |
| `(1,2,2,1)` | `(16,2,11,10)` | `2/17` |
| `(2,1,1,2)` | `(29,28,37,23)` | `6/43` |
| `(2,2,2,2)` | `(29,2,11,23)` | `1/8` |

All four base rows are loose.  This is not a uniform result for arbitrary
`u_r+39k_r`; the order-three collar remains a named structured residual.

### 4.1 The lift-invariant `q=39` equality clock

The exceptional family nevertheless has a uniform boundary skeleton.  For
each feasible parity pattern and every choice of nonnegative lift heights
`k_r`, put `u'_r=u_r+39k_r`.  Each of the three cosets has eight times

```text
coset aH         numerator clock modulo 39              active core pairs
H                2,10,11,16,23,28,29,37                {12,27},{18,21}
2H               1,5,8,14,25,31,34,38                  {3,36},{15,24}
4H               4,7,17,19,20,22,32,35                 {6,33},{9,30}. (20a)
```

For the row's coset, take `t` to be a listed numerator divided by 39.  The
complete packet in (20), with `u_r` replaced by `u'_r`, has margin exactly
`1/13`.  Reduction modulo 39 shows that the coset's same eight numerators work
for all four parity patterns.  At each time the two active speeds form one of
the listed core pairs and have signed residues `+3,-3`.  Thus the equality
point is a strict local cusp of the packet minimum: either one-sided
perturbation immediately lowers one active core constraint.

This does not prove the packet tight.  It proves that arbitrary lift height is
invisible on a common eight-point clock and that strict looseness cannot be
obtained by perturbing those clock witnesses.  A uniform closure must find a
different core-safe component or prove that the four lifted exception combs
cannot cover all such components.

## 5. Tournament Analysis and assumption challenge

On the four replacement labels use the antisymmetric cross-capacity difference

```text
C_3(r,s)-C_3(s,r)                                       (21)
```

as pairwise observable.  Switch from the left germ `(-1/13,1/13]` to the
right germ `[-1/13,1/13)`.  Every pair in `H` ties in both gauges, so declare
the numerical order `(1,5,8,12)` as the tie Hamiltonian path.

Both tournaments are the same transitive tournament: score sequence
`(3,2,1,0)`, no directed triangle, four singleton SCCs, and one Hamiltonian
path.  The gauge has zero tournament edge flips.  Nevertheless the exact
incidence changes on eight owner/colour entries: the left gauge hits each
owner by its own colour and its two `5,8` neighbors, whereas the right gauge
misses the own colour and admits the negative colour.  Moreover the four
feasible sheet tilings are distinguished only by the two parity bits in (16).

Thus the bare tournament destroys both the oriented incidence matrix and the
sheet allocation while reporting no change at all.  Runner, residue, gap, or
fixed-section vertices have the same defect.  The predicate-preserving local
carrier is the bipartite incidence between the twelve owner-sheet obligations
`(o,ell)` and the four replacement deck classes, decorated by the half-open
germ flag and `e_r`.  The tournament is useful telemetry only after that
carrier has been retained.

## Exact replay

`lrc13_hamming_four_oriented_deck_codex_S10.py` independently counts every
oriented interval `(-D,D]` through order 999 and verifies the residue-block
identity, derives the exact cutoff sets (8)--(9), runs all
141,230 finite Cayley rows in the proof, reconstructs all sixteen parity
patterns and the four exact tilings, verifies the lift-invariant clock (20a),
evaluates the four base-packet maxima, and records both tournament gauges.  The
infinite residue-block identity itself is proved by the two outer length-13
blocks in (5); the direct scan is an independent finite replay, not its proof.
The canonical output is stored beside the script.
