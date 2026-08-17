---
id: THM-3535
title: "Fixed Keller full wreath monodromy and all-level constant-linear primitivity"
status: >
  PROVED + VERIFIED-EXACT.  For every n>=1, the geometric and arithmetic
  monodromy groups of the nth iterate of the fixed sporadic Keller map are
  the full iterated wreath products W_n=S3 wr ... wr S3 in their natural
  degree-3^n actions.  The only intermediate fields of the inverse tower are
  its ancestor fields.  Consequently every nonzero constant linear form
  alpha*x+beta*y+gamma*z, in particular each literal coordinate, is a
  primitive element at every depth.  The proof uses THM-3530's unique newest
  image prime to isolate one bottom-block transposition and does not use
  countable avoidance.
source: codex/fixed-linear-wreath/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3530-fixed-keller-all-level-image-prime-and-component-tower
related:
  - THM-3508-level-two-sporadic-keller-three-coordinate-primitive-discriminant-square-class
  - THM-3519-level-three-sporadic-keller-three-coordinate-primitivity-and-common-discriminant-class
  - THM-3531-fixed-keller-intrinsic-all-level-discriminant-square-class
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_fixed_linear_all_level_wreath_probe_20260816.py
output: 05-knowledge/results/keller_fixed_linear_all_level_wreath_probe_20260816.out
script_sha256: 54775a030eed27cf5a24b497d8271c18f5a58b95881b9b8f96c79ae1630c1ae3
output_sha256: 8f7fb51113a2e065467df91ee1314a91c101423038a52dba28cd97f44c4e3b23
semantic_sha256: 5810ca665b68ad0f8bdfc5a62445c7eb1b5d8d65ee6ae565db1e31bd1a31f9c9
hash_basis: LF-normalized bytes
---

# THM-3535 -- one fixed linear observation separates every inverse depth

**PROVED + VERIFIED-EXACT.**

Let `F:A^3->A^3` be the fixed sporadic Keller map of THM-2473.  For a
characteristic-zero constant field `k`, write an inverse chain as

```text
F(q_i)=q_(i-1),
K_i=k(q_i),
K_0 subset K_1 subset ... subset K_n,
[K_i:K_(i-1)]=3.                                      (1)
```

The coordinates of every `q_i` are algebraically independent over `k`.
Let `G_n` be the geometric monodromy group of `K_n/K_0` in its natural
action on the `3^n` leaves of the inverse tree, and put

```text
W_1=S_3,
W_n=S_3^(3^(n-1)) semidirect W_(n-1).                 (2)
```

## 1. The theorem

For every `n>=1`,

```text
G_n=W_n.                                               (3)
```

The arithmetic monodromy over `Q(a,b,c)` is also `W_n`: it contains the
geometric group and is already contained in the same compositional wreath
product.

The intermediate fields of `K_n/K_0` are exactly

```text
K_0,K_1,...,K_n.                                       (4)
```

Moreover, for every constant triple

```text
(alpha,beta,gamma) in k^3 \ {(0,0,0)},
theta_n=alpha*x_n+beta*y_n+gamma*z_n,                  (5)
```

one has

```text
K_0(theta_n)=K_n                                      (6)
```

for every `n>=1`.  In particular, one may choose one fixed rational form
once and for all; indeed every nonzero rational form works simultaneously.
The literal coordinates `x_n`, `y_n`, and `z_n` are primitive at all levels.

## 2. The newest prime isolates one bottom-block transposition

Use THM-3530's raw image primes `P_j`, with `P_0=L`.  It proves that they are
pairwise distinct and absolutely irreducible, that

```text
S_(F^r)=V(P_0...P_(r-1)),                              (7)
```

and that every restriction

```text
V(P_j) -> V(P_(j+1))                                  (8)
```

induced by `F` has generic degree one.

Fix `n>=2` and let `eta` be the generic point of `V(P_(n-1))`.  By (7),
`eta` lies outside `S_(F^(n-1))`, because `P_(n-1)` is distinct from all
older primes.  Therefore, on a neighbourhood of `eta`,

```text
F^(n-1):A^3 -> A^3                                    (9)
```

is proper and etale, hence finite etale of degree `3^(n-1)`.  Over a strict
henselization at `eta`, (9) splits into `3^(n-1)` unramified sections.

The composite of the degree-one maps (8) gives

```text
V(L) -> V(P_(n-1))                                    (10)
```

generic degree one.  Thus exactly one section in (9) meets `V(L)`.  The
meeting is simple: (9) is etale and both divisors are reduced.  On every
other section `L` is a unit.

For the final inverse step, THM-2473's first-coordinate core and its
reciprocal are

```text
E(w)=Lw^3+Tw-2c,
E(1/u)u^3=L+Tu^2-2cu^3,             T=4-3bc.           (11)
```

At the generic point of `V(L)`, the functions `c,T`, and
`S=27ac^2-9bc+8` are units.  For example, the point `(2/27,1,1)` on `V(L)`
has `(c,T,S)=(1,1,1)`, so none vanishes identically on the irreducible
divisor.  The exact discriminant identity is

```text
disc(E)=-4L(T^3+27c^2L)=-4LS^2.                       (12)
```

The reciprocal Newton polygon has points `(0,1),(2,0),(3,0)`.  Hence it has
one unit root and two roots of valuation `1/2`; equivalently, (11) has one
finite unramified branch and one tamely ramified quadratic pair at infinity.
Its inertia is one transposition.

The preceding cover (9) is unramified, only one of its sections has `L=0`,
and every other last-step cover is finite etale.  Consequently `G_n`
contains an inertia element which is a transposition on exactly one bottom
ternary block and fixes every other leaf.  This support statement, not only
the odd discriminant class, is the decisive sidecar.

## 3. One local transposition forces the full wreath product

We use the following elementary group lemma.

> Let `G` be a transitive subgroup of
> `S_3^m semidirect H`, preserving `m` blocks of size three and surjecting
> onto a transitive `H`.  If `G` contains a transposition supported on one
> block, then `G` contains the whole base group `S_3^m`.

Indeed, conjugation moves the supported transposition to every block.  For
one fixed block `B`, transitivity of `G` on all leaves implies that its
setwise stabilizer is transitive on `B`: an element carrying one point of
`B` to another must preserve `B` by the block property.  Its induced group
on `B` is therefore a transitive subgroup of `S_3` containing a
transposition, hence all of `S_3`.  Conjugating the supported transposition
inside that stabilizer produces a second transposition on `B`; the two
generate a copy of `S_3` supported only on `B`.  Conjugating between blocks
gives every factor of `S_3^m`.

THM-2473 proves `G_1=S_3`.  Suppose inductively that `G_(n-1)=W_(n-1)`.
The compositional block action gives

```text
G_n subset S_3^(3^(n-1)) semidirect W_(n-1),           (13)
```

and restriction to the preceding normal closure surjects onto
`W_(n-1)`.  The source of `F^n` is irreducible, so `G_n` is transitive on
its leaves.  Section 2 supplies the supported transposition.  The lemma
therefore gives the full base group in (13), and the surjective quotient
gives `G_n=W_n`.  This proves (3) by induction.  In particular, it upgrades
THM-2473's former finite-field evidence for `G_2=W_2` to a proof and proves
the same statement at every depth.

## 4. Full ternary wreath actions have only ancestor blocks

Identify the leaves of `W_n` with ternary words of length `n`.  Its bottom
blocks are the three children of each word of length `n-1`.

We prove by induction that a block containing a fixed leaf is exactly one of
its ancestor cylinders.  Let `B` be such a block.  If `B` lies inside one
bottom ternary block, primitivity of the local `S_3` action makes `B` either
a singleton or that whole ternary block.

Otherwise `B` meets at least two bottom blocks.  A local permutation
supported on either met block fixes a point of `B` in the other block, so
its translate intersects `B`; the block axiom forces equality.  Hence the
intersection with each met bottom block is invariant under its local
`S_3`, and is therefore the whole bottom block.  Thus `B` is a union of
bottom blocks, and its projection is a block for `W_(n-1)`.  The induction
hypothesis makes that projection an ancestor cylinder.  This proves the
classification.

For a separable nonnormal extension, subgroups between a leaf stabilizer
and the full monodromy group are in bijection with blocks containing that
leaf.  The ancestor cylinders correspond exactly to the fields in (1).
This proves (4).  There are no diagonal, cross-branch, or other hidden
intermediate fields.

## 5. No nonzero constant linear form descends one step

It remains to exclude membership of (5) in the immediate parent field.  The
following two exact fibres from THM-2473 do this uniformly.  For every
integer `m>=1`, over

```text
t_m=(-1/(4m^2),0,0)                                   (14)
```

the three points are

```text
p_m^0=(0,0,-1/(4m^2)),
p_m^+=(m,-3/(2m),13/(2m^2)),
p_m^-=(-m,3/(2m),13/(2m^2)).                           (15)
```

The target (14) is outside `V(L)`, since `L(t_m)=-4/m^2`, so this is a
reduced finite-etale fibre.  If a constant linear form descended through
`F`, it would be constant on every generic fibre.  The difference of its
two pullbacks is regular on the finite-etale fibre product over the
complement of `V(L)`; generic vanishing therefore extends to every fibre
there, including (14).

Use `p_1^+-p_1^-`, `p_2^+-p_2^-`, and `p_1^+-p_1^0`.  Their coordinate
rows are

```text
(2,-3,0),
(4,-3/2,0),
(1,-3/2,27/4),                                         (16)
```

whose determinant is

```text
243/4 != 0.                                            (17)
```

Thus the only constant linear form that is constant on both fibres is zero.
Consequently, for every nonzero triple in (5),

```text
theta_n notin K_(n-1).                                 (18)
```

If `K_0(theta_n)` were proper in `K_n`, (4) would put it inside some
`K_j` with `j<n`, hence inside `K_(n-1)`, contradicting (18).  This proves
(6).  Notice that the same three rows have full rank over every
characteristic-zero constant field, so the conclusion is uniform in the
constant direction, not a countable-choice statement.

## 6. The countable-avoidance trap is bypassed, not solved abstractly

At each fixed depth, primitive-element theory gives a nonempty Zariski-open
set of coefficient triples.  It would be invalid to infer a rational point
in the intersection of all those opens merely from countability.  Enumerate
`P^2(Q)={r_1,r_2,...}` and let `B_s` be the union of its first `s` points.
Every `P^2\B_s` is a nonempty rational Zariski open, while

```text
intersection_s (P^2(Q)\B_s)=empty.                    (19)
```

The proof above never invokes such an intersection.  It proves (18) for all
nonzero constant directions at once, and the wreath block theorem converts
that one-step statement to every depth.

## 7. Discriminant consequence and exact boundary

By THM-3531, the intrinsic trace-discriminant square class is

```text
d(K_n/K_0)=[(-1)^nP_(n-1)].                            (20)
```

Because every form (5) is primitive, its monic degree-`3^n` minimal
polynomial has exactly the class (20).  Any saturated full-degree eliminant
differs by a target-field scalar, whose discriminant exponent
`2*3^n-2` is even, so it has the same class.  This recovers all three
literal-coordinate statements at every level without expanding their
eliminants.

THM-3533 separately sharpens the normalization discriminant: its newest-prime
coefficient is one, while a primitive power order has exponent `1+2i`, with
`i` its local index length.  The present theorem proves primitivity but does
**not** prove `i=0` for any prescribed constant-linear form.  It also does not
determine exact old-prime multiplicities or discriminant-divisor
singularities, and does not say that every specialized finite fibre is
separated by every form.  It concerns constant linear forms for this one
fixed map; arbitrary polynomial observations, tame-family covariance, other
Keller maps, monoid-wide monodromy, `JC(2)`, `DC(2)`, and LRC(14) remain
outside its scope.  Full wreath monodromy also does not make `F^n` proper or
invertible: (7) still has `n` nonproperness components.

## 8. Exact companion

The companion independently checks:

1. (12), the reciprocal Newton slopes, and the one-transposition local type;
2. the three rows (16) and determinant (17);
3. minimal invariant blocks for the exact full wreath generators through
   depth four, with histograms
   `1`; `1,3`; `1,3,9`; `1,3,9,27`; `1,3,9,27,81` in the prescribed
   ancestor multiplicities;
4. the countable-avoidance hostile (19); and
5. an orthogonal finite-field probe over `F_41`.

For the last probe it enumerates all `41^3` source points.  The exact
one-level fibre histogram is

```text
1:35261, 3:11220,
```

and the two-level histogram is

```text
1:19926, 2:2851, 3:7919, 4:1782, 5:1393,
6:288, 7:469, 9:48.                                   (21)
```

Across those 48 split degree-nine fibres, every one of the
`41^2+41+1=1723` projective constant-linear directions separates all nine
points on at least one fibre.  The direction-to-first-witness ledger has
SHA-256

```text
3e291809e49a1c16652c8bbacfe04c746c7bcb35b0fe08c51968d3d661deb493.
```

This finite-field census is a hostile/supporting probe, not the all-level
proof; Sections 2--5 carry the theorem.

Reproduce with

```text
python -B 04-computation/keller_fixed_linear_all_level_wreath_probe_20260816.py
python -B -O 04-computation/keller_fixed_linear_all_level_wreath_probe_20260816.py
```

The ordinary and optimized transcripts match the stored output exactly.

**QED.**
