---
id: THM-3542
title: "Fixed Keller depth-three newest-prime residue point/pair saturation"
status: >
  PROVED + VERIFIED-EXACT.  A good rational specialization of the depth-three
  predecessor cover has point-orbit degrees 1,2,6 and injective unordered-
  pair-sum orbit degrees 1,2,6,6,9,12.  Specialization orbits refine the
  generic newest-prime residue action, which in turn refines the marked
  ternary-tree stabilizer.  Equality of both orbit counts proves the full
  point-and-pair saturation gate at n=3 and reduces the raw carrier to nine
  valuation packets.  This does not identify the full decomposition group or
  prove saturation at n>=4.
source: codex/depth-three-residue-saturation/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
  - THM-3538-fixed-keller-newest-prime-prescribed-coordinate-index-criterion
  - THM-3539-fixed-keller-newest-prime-decomposition-centralizer-and-lca-packet-floor
related:
  - THM-3540-fixed-keller-depth-two-newest-prime-residue-saturation
script: 04-computation/keller_depth_three_newest_prime_residue_saturation_20260816.py
output: 05-knowledge/results/keller_depth_three_newest_prime_residue_saturation_20260816.out
script_sha256: 93ebfa6fdb0e176bb5c365fd6e5f7df33ba8ec6ac03e00f7fdff5b73612ae103
output_sha256: da0dba42b49fdb6a789e6ebbd99786cb4138b8c02328aa9c13ed4adb9eac2b8d
semantic_sha256: d79f24556807af9c7c335104564bb9f6baf047d42ca46a0123da669a228e6835
hash_basis: LF-normalized bytes
---

# THM-3542 -- one good fibre saturates the depth-three residue orbits

**PROVED + VERIFIED-EXACT.**

Retain the fixed sporadic Keller map `F`, its raw image-prime tower

```text
P0=L,             P1=H,             P2=J,
```

and THM-3539's image `bar D` of the newest-prime decomposition group on the
nine predecessor blocks over `kappa(J)`.  THM-3539 proves only

```text
bar D <= H2:=Stab_(W2)(00),                            (1)
```

where `W2=S3 wr S3` acts on the nine depth-two ternary words.  This theorem
proves that `bar D` and `H2` have identical orbits on both points and
unordered pairs.  It proves the saturation gate needed by THM-3539, without
claiming equality of the groups.

## 1. The depth-three residue field has a rational two-parameter chart

THM-2570 gives the normalization of `V(L)`:

```text
nu(tau,lambda)=
 (lambda^2(3-tau lambda)/27,
  lambda(4-tau lambda)/3,
  tau).                                                (2)
```

THM-3530 proves that each map

```text
V(Pj) -> V(P_(j+1))                                   (3)
```

has generic degree one.  Therefore `F^2 o nu` is dominant and birational
onto `V(J)`, and

```text
kappa(J)=Q(tau,lambda).                                (4)
```

Over the generic point of `J`, the cover `F^2` is finite etale of degree
nine by THM-3535.  The unique branch induced by `(2)` is the marked ancestry
block `00`; its splitting-field action on the nine blocks is exactly the
group `bar D` in `(1)`.

## 2. One exact specialization is a good nine-point fibre

Specialize `(2)` at

```text
(tau,lambda)=(0,1),       q0=(1/9,4/3,0) in V(L),
q1=F(q0),                  q2=F(q1).                    (5)
```

The first image and two useful old-`L` values are

```text
q1=(19840/2187,964/243,14/81),
L(q1)=7664/81,
L(q2)=-2129081689932006656253463984/
       984770902183611232881.                          (6)
```

Both values in `(6)` are nonzero.  The exact coordinates of `q2` are retained
in the companion's semantic ledger.

For a target `(a,b,c)`, THM-2473's inverse `x`-core is

```text
E_(a,b,c)(W)=L(a,b,c)W^3+(4-3bc)W-2c.                 (7)
```

Once `W` is a root, put

```text
D=(12a-b^2)W^2+bW+2,
Y=b-3aW((9ac-b)W+2)/D,
Z=(2W-c-3W^2Y)/W^3.                                  (8)
```

The companion applies `(7)`--`(8)` first to `q2`, then applies `(7)` to the
reconstructed middle point.  It checks exactly that:

1. the outer cubic is squarefree;
2. every denominator in `(8)` is a unit modulo that cubic;
3. substituting the reconstructed point into `F` gives `q2` modulo the
   outer cubic;
4. `L` of every middle point is a unit, so every child core remains cubic;
5. the second reconstruction denominator has nonzero resultant on every
   outer root; and
6. eliminating the middle `x` coordinate gives a degree-nine polynomial
   whose constant term is nonzero and whose derivative is coprime to it.

Thus the resultant is not carrying a root at infinity, a pole of the inverse
section, or a repeated coordinate.  Its primitive polynomial `P_s(X)` is the
actual `x`-coordinate eliminant of the nine reduced points in

```text
F^(-2)(q2).                                            (9)
```

The literal `x` coordinate separates all nine points in this fibre.

## 3. Exact point and pair orbit factorizations

Exact factorization over `Q` gives

```text
deg P_s=9,       irreducible-factor degrees 1,2,6.     (10)
```

The marked child cubic `E_(q1)` divides `P_s` and has factor degrees `1,2`;
its linear root is `X=1/9`, the marked point `q0`.  Hence `(10)` respects the
ancestry split rather than merely matching an unordered degree list.

Let the roots of the monic version of `P_s` be `r_1,...,r_9`.  Define

```text
D_s(T)=product_i (T-2r_i)=2^9 P_s(T/2),
R_s(T)=product_(i<j) (T-r_i-r_j).                     (11)
```

Up to a nonzero rational scalar, the symmetric resultant identity is

```text
Res_X(P_s(X),P_s(T-X))=D_s(T) R_s(T)^2.               (12)
```

The exact quotient in `(12)` is a square after monic normalization.  Its
square root `R_s` has degree `36`, is squarefree, and factors over `Q` with
degrees

```text
1,2,6,6,9,12.                                         (13)
```

Squarefreeness is load-bearing: it says all 36 unordered pair sums are
distinct, so `R_s` loses no pair information.  Exact characteristic-zero
factorization proves every factor in `(10)` and `(13)` irreducible.  As an
independent arithmetic fingerprint, the companion also intersects the
possible factor degrees from squarefree reductions.  For example, the
degree-nine factor is certified by the two patterns

```text
mod 13: 3+6,             mod 19: 1+4+4.               (14)
```

Their subset-degree intersection is `{0,9}`.  The degree-twelve factor stays
`6+6` in the recorded modular fingerprint; this is an imprimitive-cycle
shadow, not a rational factor, and is kept explicitly rather than being
misreported as a failed transitivity test.

Since `P_s` separates points, its three irreducible factors are precisely the
three orbits of the specialized absolute Galois group `G_s` on `(9)`.  Since
the pair-sum map is injective, the six irreducible factors of `R_s` are
precisely its six orbits on unordered pairs.

## 4. Specialization forces generic saturation

We use the elementary good-specialization lemma: if a finite etale cover of
an integral rational variety is specialized at a rational point where its
degree and discriminant survive, then, after compatible choices of geometric
points, the specialized Galois image is a subgroup of the generic Galois
image.  This follows by localizing an integral model until the Galois closure
is finite etale; a prime over the rational point has a decomposition group
inside the generic Galois group, and its residue action is the specialized
splitting action.

The gates in Section 2 place `(5)` in this good locus.  Therefore

```text
G_s <= bar D <= H2.                                   (15)
```

The companion independently enumerates all `144` elements of `H2`.  Its
point and pair orbit sizes are

```text
points:             1,2,6;
unordered pairs:    1,2,6,6,9,12.                     (16)
```

For subgroup inclusions, the smaller group's orbit partition refines the
larger group's partition.  Equations `(10)`, `(13)`, `(15)`, and `(16)` show
that the first and last partitions have the same number of cells: three on
points and six on pairs.  A proper refinement would increase that number.
Consequently every refinement in `(15)` is equality at the level of these
two orbit sets:

```text
Orb_(G_s)(B)=Orb_(bar D)(B)=Orb_(H2)(B),
Orb_(G_s)(binom(B,2))=Orb_(bar D)(binom(B,2))
                     =Orb_(H2)(binom(B,2)).            (17)
```

This proves THM-3539's gate `G_(2-orb)` at `n=3`.  Notice that `(17)` does
not imply `G_s=bar D=H2`; distinct subgroups can have the same point and pair
orbits.

## 5. The exact nine-packet descent

Use THM-3539's shell labels with predecessor depth `r=2`.  The two off-ray
shells have sizes

```text
s_0=6,                 s_1=2.                          (18)
```

For any one of THM-3538's integral observations `theta=y,z,u=1/x`, equation
`(17)` turns its raw 45-factor carrier into

```text
i_(theta,3)
 =v(A_*)+6v(A_0)+2v(A_1)
 +6v(R^*_0)+2v(R^*_1)
 +12v(R^<_(0,1))
 +9v(R^=_(0,0))+6v(R^=_(0,1))+v(R^=_(1,1)).           (19)
```

The multiplicities in `(19)` sum to

```text
internal: 1+6+2=9,       pairs: 6+2+12+9+6+1=36.     (20)
```

There are exactly

```text
3 internal packets + 6 pair packets = 9=3^2          (21)
```

as predicted by THM-3539's quadratic floor.  Since every valuation is
nonnegative, the carrier is a unit exactly when these nine representatives
are units.  THM-3538 already determines the prescribed-coordinate indices at
this depth; the new content is the lawful residue descent of the carrier, not
a new unit calculation.

## 6. Equality and failure boundary

THM-3542 proves no statement at `n>=4`.  In particular, it does not prove:

1. `bar D=H2`, full decomposition-centralizer equality, or the residue-field
   degree;
2. point/pair saturation at depth four or at all depths;
3. that a simple pair-sum remains injective in a larger specialization;
4. unitness of any previously open THM-3538 packet;
5. multiple atom decompositions for one fixed Keller map, an arbitrary-map
   classification, `JC(2)`, or a counterexample in dimension two; or
6. a current, spectral closure, row exclusion, or LRC(14).

The next economical target is not a full depth-four Galois group.  It is one
good 27-point specialization together with a separating observable on its
351 unordered pairs.  Matching the four point orbits and twelve pair orbits
of `H3` would prove the next saturation gate by the same refinement argument.

## 7. Exact companion

Reproduce with

```text
python -B 04-computation/keller_depth_three_newest_prime_residue_saturation_20260816.py
python -B -O 04-computation/keller_depth_three_newest_prime_residue_saturation_20260816.py
```

The normal and optimized transcripts match the stored output exactly.  The
companion contains no executable `assert`; it records LF-normalized file
hashes, semantic ledger hash, exact fibre gates, factor ledgers, modular
fingerprints, and the independent 144-element orbit census.

**QED.**
