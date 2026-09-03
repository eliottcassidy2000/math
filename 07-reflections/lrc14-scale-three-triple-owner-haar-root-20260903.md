---
title: "LRC(14) scale-three owner permutation and triple-comb Haar scout"
status: >
  PROVED exact three-sheet owner criterion and body-specific Haar transfer;
  PROVED sharp 6/77 bound and unique (1,5,11) equality case for the complete
  signed-(1,2,1) resonant sector; VERIFIED-EXACT finite census through odd
  tail height 199. The nonresonant universal 6/77 bound remains OPEN, and
  LRC(14) remains open.
source: root / LRC14 continuation session, 2026-09-03
artifacts:
  - 04-computation/lrc14_scale_three_triple_owner_probe_root_20260903.py
  - 05-knowledge/results/lrc14_scale_three_triple_owner_probe_root_20260903.out
script_sha256: 72cd5142eb23570f94ccb909ec4270be05bb2c7c4d48c76adc5d470ff9a07e22
output_sha256: 637a7bafc4baa457d4c8cec959e6c97e06e5066096c62c9ae6ec767873c68226
hash_basis: raw LF bytes
related:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
  - THM-2060-crt-tail-coset-saturation
  - THM-2064-multitail-sheet-capacity-and-dyadic-seam
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
---

# Scale three: retain the owner permutation

## 1. Inheritance and the exact surviving seam

THM-2060/2064 gives the sharp common-fibre capacity

```text
sum_i ceil(q_i/7)/q_i < 1.
```

For three nonabsorbed tails, every odd relative order is at least three, so
the only equality left by this scalar bound has

```text
q_1=q_2=q_3=3.                                         (1)
```

Primitivity then reduces the common scale to three.  This is precisely the
prime-three boundary left visible in THM-4004.  Its control with owners
`(1,11,10)` proves that an arbitrary lower-dimensional witness cannot be
lifted blindly.  Replacing the even tail by `13` gives the all-odd local
hostile `(1,11,13)` at the same quotient phase `y=1/11`: parity alone does
not remove the three-colour obstruction.

The live object is therefore not the three capacities but their common owner
permutation.

## 2. Exact owner criterion

Let `C` be a finite positive core and let `a,b,c` be positive integers prime
to three.  For `y in G_C`, consider the three physical lifts

```text
x_j=(y+j)/3,                    j in Z/3Z.              (2)
```

For one tail `w`, it is dangerous on some lift exactly when

```text
||wy||<3/14.                                             (3)
```

When `(3)` holds, let `k_w(y)` be the unique nearest integer to `wy`.  The
unique dangerous lift is

```text
o_w(y)=-w^(-1) k_w(y) mod 3.                            (4)
```

Indeed, danger on lift `j` is equivalent to

```text
|wy+wj-3n|<3/14
```

for some integer `n`.  The integer within `3/14` of `wy` is unique, and its
congruence modulo three gives `(4)`.

It follows that all three lifts over `y` are spoiled by `{a,b,c}` if and
only if

```text
all three inequalities (3) hold, and
{o_a(y),o_b(y),o_c(y)}=Z/3Z.                            (5)
```

Thus the exact equality seam is a three-colour permutation, not a capacity
equality.

## 3. Six disjoint comb intersections and Haar transfer

On the physical `x`-circle put

```text
D_(w,j)={x:||w(x+j/3)||<1/14}.                          (6)
```

For `3` not dividing `w`, the three sets `D_(w,j)` are pairwise disjoint.
The complete three-tail failure comb is the disjoint union, up to endpoints,

```text
F_(a,b,c)=disjoint_union_(pi in S_3)
           D_(a,pi(0)) intersect D_(b,pi(1)) intersect D_(c,pi(2)). (7)
```

Translation by `1/3` cyclically rotates all three owner labels.  Hence the
six terms form two equal-measure orbits and

```text
mu(F_(a,b,c))
 =3 mu(D_(a,0) intersect D_(b,1) intersect D_(c,2))
 +3 mu(D_(a,0) intersect D_(b,2) intersect D_(c,1)).   (8)
```

Formula `(8)` is exact and computable by intersecting three sorted lists of
rational danger teeth; no sampling is required.

There is an immediate body-specific transfer.  Multiplication by three sends
the physical body-safe set for `3C` onto `G_C` and preserves Haar measure. If
`G_C` is nonempty and

```text
mu(G_C)>=mu(F_(a,b,c)),                                 (9)
```

then

```text
3C union {a,b,c}
```

has an LRC(14) witness.  For a strict inequality this is direct measure
comparison.  Equality also works: the body-safe preimage is nonempty compact,
whereas `F_(a,b,c)` is a proper open subset of the connected circle (`x=0`
is not in it).  A compact subset of a proper open set cannot have the same
Haar measure as that open set, by the same compact/open argument as THM-4150.

This is a proved transfer.  THM-4373 now supplies a universal `6/77` bound on
the full signed-`(1,2,1)` resonant sector, but the transfer remains tail-
specific outside that sector.  The nonemptiness hypothesis is automatic for
the ten-speed bodies in the proposed LRC(14) application by `LRCUpTo13`; it
is necessary in the arbitrary finite-core formulation when the failure comb
has zero measure.

## 4. Exact finite signal and the proved resonant sector

The companion performs two independent exact calculations of `(8)`:

1. direct intersection of the six rational tooth combs; and
2. a complete wall-cell classification on one `1/3`-period.

They agree on every primitive odd triple through `31`, under label
permutations and a common dilation by five.  The cached interval
implementation then checks all `47,499` primitive triples of distinct odd
nonmultiples of three through `199`.  The unique maximum is

```text
(a,b,c)=(1,5,11),              mu(F)=6/77,              (10)
```

and the next values are

```text
(1,11,23): 12/161,
(1,25,49): 24/343,
(5,7,17):  58/833.                                      (11)
```

The entire top thirty satisfies, after a permutation and choice of signs, a
relation with absolute coefficient pattern `(1,2,1)`.  Of the `47,499`
triples, `1,023` have such a relation.  Their maximum is `(10)`.  Among the
other `46,476`, the maximum drops sharply to

```text
(1,11,43):                 12/301.                      (11a)
```

THM-4373 proves that the resonant maximum is global within that sector.  If
`p<q` are the two coefficient-one speeds, then the coefficient-two speed is

```text
b=(p+q)/2 when q=p mod 12,
b=(q-p)/2 when q=-p mod 12.                             (12)
```

Primitivity is `gcd(p,q)=1`.  With

```text
A=3(q-p)/28,                       B=3(q+p)/28,
f(t)=2/(pq)((B-t)_+-(A-t)_+),
```

the exact resonant formula is

```text
mu(F_(p,b,q))=2 sum_(k>=1, 3 does not divide k) f(k).   (13)
```

The determinant index `2k=q n_p-p n_q` records a pair of nearest endpoint
integers.  Endpoint eligibility together with the even-determinant parity
condition makes the middle integer integral, and the middle error is their
half-sum or half-difference, so endpoint
eligibility already forces middle eligibility; the three owners permute
exactly when `3` does not divide `k`.

Define the period-three quadrature error

```text
E(t)=sum_(k>=1, 3 does not divide k)(t-k)_+-t^2/3.
```

It has range `[-1/3,0]`, and `(13)` becomes

```text
mu(F_(p,b,q))=3/49+4/(pq)(E(B)-E(A))
             <=3/49+4/(3pq).                            (14)
```

For `pq>=81`, the last expression is strictly below `6/77`.  The congruence
parametrization leaves exactly thirteen pairs with `pq<81`; direct use of
`(13)` has unique maximum `(p,q)=(1,11)`, hence `(p,b,q)=(1,5,11)`.

For `(1,5,11)`, the failure comb consists of six components, each of length
`1/77`; geometrically, six full teeth of the fastest tail are nested in the
two slower translated combs.

The remaining sharp conjecture is

```text
mu(F_(a,b,c))<=6/77                                    (15)
```

for every primitive triple of distinct positive odd integers prime to three,
with equality only at `(1,5,11)` up to permutation.  Before primitive
normalization, common dilation preserves the same measure because
multiplication by the dilation is Haar-preserving and, being prime to three,
permutes the three sheet labels.
THM-4373 proves this statement for every signed-`(1,2,1)` resonant triple.
The nonresonant statement is **OPEN**, not inferred from the finite census.
If completed, `(9)` would give a scale-three analogue of THM-4150: every
ten-speed body with safe-set measure at least `6/77` would accept every three
odd order-three tails.

## 5. Information ledger

```text
source:      all three labelled lift owners over the complete core-safe set
target:      the prime-three equality seam of common-fibre capacity
map:         nearest-integer residue -> owner colour in Z/3Z
preserved:   strict danger, physical sheet, tail label, common dilation
destroyed:   core-component address after Haar scalarization
sidecar:     the exact six-intersection comb F_(a,b,c)
hostile:     odd triple (1,11,13) realizes a local owner permutation
test:        exact rational interval intersection plus independent wall cells
```

The computation independently checks THM-4373's endpoint and quadrature
formulas on all `1,023` resonant triples through `199`, including the thirteen
small endpoint pairs.  It supplies finite evidence only on the nonresonant
sector, does not lower the general LRC(14) frontier, and does not replace the
physical entry conditions of THM-4004.
