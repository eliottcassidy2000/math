# Prime defects, two triangle sectors, cube dissections and an elliptic tree

**Session begun 2026-09-21; closed 2026-09-22. PROVED scoped mechanisms;
CITED external classifications; FINITE-EXACT controls. Collatz OPEN.**
No new exceptional-prime discovery, full Collatz/Berggren isomorphism,
or new Lean formalization is claimed.

## Inheritance and research board

The closest proved mechanism is [THM-3756, odd-square ordinal Berggren
descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md),
extended by the [previous odd-square session](odd_square_20260921_synthesis.md).
The current hostile is not merely the existence of two u=17 cycles: it is
the mismatch between a sector-exchanging arithmetic action and a cube isometry.
The corrected near miss is treating a ternary symbolic tree as a positive
solution tree. The least-used sidecars are the coset character, prime-power
lift depth, torsion coordinate and sampling/metric information.

Anchor: how the two u=17 sectors differ and lift. Niche: harmonic and
binomial prime defects. Wildcard: exact dissections and group-law trees.

| Live concept | Exact operation | What the operation does not automatically retain |
|---|---|---|
| Triangle shell | Multiplication by2 modulo sign | A switch between different cosets |
| Prime defect | Valuation of a return minus1 | Other exceptional-prime predicates |
| Harmonic sums | Reciprocal pairing and binomial expansion | A forward Collatz orbit |
| Cube dissection | Isometries, affine transport, parity | Arithmetic multiplier action and all lengths |
| Elliptic tree | 3P+iH with a free coefficient | Positivity and the separate torsion coordinate |
| Collatz | Exact guarded affine edges | Global descent from any of these local models |

The research cards used were **Type every analogy and every implication**,
**Retain the local profile until its global counting weights are known**,
and **Separate unbounded local support from a height-bounded modular cover**
in [META-PATTERNS](../../00-navigation/META-PATTERNS.md). This session applies
them; it does not add a card based on one family of examples.

## 1. The two numerical sequences are one arithmetic invariant

The user's clarification identifies the primes as p=5,7,11,13,17,19,23,... .
In the localized rational ring where denominators are coprime to p, put

```
eta_p = H_(p-1)/p² mod p,
beta_p = (binom(2p-1,p-1)-1)/p³ mod p.
```

The [prime lane](prime_shells_20260921_primes.md) proves

```
beta_p = 2 eta_p mod p, for every prime p>=5.
```

Their initial rows are respectively `(3,6,6,7,10,14,18)` and
`(1,5,1,1,3,9,13)`. Thus the second sequence is a dependent observable,
not another independent source of irregularity. The first zero condition
is exactly the Wolstenholme extra divisibility. H4=25/12 and H6=49/20
are the initial cases of the underlying p² divisibility.

The claimed binomial is repaired as `binom(9,4)=binom(9,5)=126=5³+1`;
`binom(9,2)=36`. In fact n=5 uniquely solves
`binom(2n-1,n-1)=n³+1` in positive integers, by the exact ratio comparison
in the prime note. This is a genuine isolated equality with an elementary
mechanism, outside Catalan's perfect-power hypotheses.
The exact prime ranks are pi(16843)=1944 and
pi(2124679)=157504, so there are exactly155559 primes strictly between
these endpoints. This is a finite sieve count, not evidence of a universal
gap formula.

A complete finite scan of primes5<=p<=20000 gives

| eta_p | Primes in this bounded universe |
|---:|---|
| 0 | 16843 |
| 1 | 103 |
| 2 | 12101 |
| 3 | 5,2851 |

This recovers a mathematical role for12101 and2851 without assigning them
the Wolstenholme property. It is a plausible reconstruction of why those
numbers were supplied, not a claim that the user's intended indexing was
fully specified.

The standard exceptional-prime notions also separate. The prime note
records primary-source Bernoulli data, independently verified regularity
of1093 and2851, positive irregular witnesses for the larger inputs, and
a full irregular-prime census through300. In particular:

- 1093 is Wieferich and regular;
- 3511 is Wieferich and irregular;
- 12101 is irregular but not Wolstenholme;
- 16843 and2124679 are Wolstenholme, and therefore irregular.

These are new-to-session classifications and reproductions of known
arithmetic, not newly discovered primes.

There is also a precise shared lifting template: for the even binomial
multiplier Q=binom(2p-1,p-1), the affine map R(n)=Qn+(Q-1)/p fixes the
first r-1 source digits, where r=v_p(Q-1). Ordinary Wolstenholme
divisibility gives r>=3; a Wolstenholme prime gives at least one extra
fixed digit. The map retains the same LTE period mechanism as an inverse
Collatz braid. Its odd multiplier factor changes the target odd core,
however, so it is not an inverse fibre of the original map.

## 2. The u=17 separation has a character, a field and a lift law

The [shell lane](prime_shells_20260921_shells.md) changes the question from
coverage to the complete component invariant. In the group
`(Z/uZ)^*/{+1,-1}`, triangle doubling is multiplication by2. Its cycles
are cosets modulo `<2,-1>`. At17 that subgroup is exactly the squares:

```
+ sector: 1 ->15 ->13 ->9 ->1;
- sector: 3 ->11 ->5 ->7 ->3.
```

Their Legendre characters are +1 and -1. Multiplication by3 exchanges
them, but its square is D³: the full scalar action is C8, not C4 x C2.
A bit describing the sector alone discards this phase relation.

The same partition survives the exact map `v -> 2cos(2pi v/17)` to
period-four points of f(x)=x²-2. Their two quartic orbit polynomials are
irreducible over Q(sqrt17), exchanged by its quadratic conjugation, and
multiply to one irreducible degree-eight polynomial over Q. The orbit
sums solve `t²+t-4=0`. This is a concrete field-theoretic form of the
two-way split associated with the regular17-gon.

For any odd prime p, let d be the least signed return exponent of2, write
2^d=epsilon mod p, and r=v_p(2^d-epsilon). Then the shell at u=p^e has

```
period = d*p^max(0,e-r),
number of cycles = (p-1)/(2d) * p^min(e-1,r-1).
```

Every lift from p^e to p^(e+1) either splits each cycle into p cycles,
or lengthens it by p. The first happens while e<r; the second afterwards.
Since r=v_p(2^(p-1)-1), this makes the Wieferich defect geometric:

| Prime | Shell p | Shell p² |
|---|---|---|
| 17 | 2 cycles, length4 | 2 cycles, length68 |
| 1093 | 3 cycles, length182 | 3279 cycles, length182 |
| 3511 | 1 cycle, length1755 | 3511 cycles, length1755 |

The Chebyshev multiplier around each cycle is exactly epsilon*2^d.
Thus r also measures its p-adic closeness to1. At p these cyclotomic
points coalesce upon reduction, so this is not a claim of a distinct
periodic orbit surviving over the residue field.

The inherited inverse affine braid `R_p(n)=2^ord_p(2)n+(2^ord_p(2)-1)/p`
has the same valuation defect shifted by one level. This transports a
period/lift theorem between the two settings. It does not turn that
inverse-fibre recurrence into a forward Collatz trajectory.

For composite conductors the local return signs must also synchronize.
At63=9*7 the local period-three signs disagree, forcing a global period6;
at65=5*13 the local signs agree at the lcm6. Their equal period has
different local mechanisms. This retains information that a bare count
of threes and twos would conceal.

## 3. A cube realizes D, but cannot realize the full scalar action

The [dissection lane](prime_shells_20260921_dissections.md) recovers the
repo's five-piece cube construction with exact coverage:
the alternating-vertex regular tetrahedron has volume1/3 and the four
corner tetrahedra each have volume1/6. The corner pieces are half-volume,
not congruent to planar halves of the regular tetrahedron.

The two complementary central tetrahedra themselves overlap in an
octahedron and have union volume1/2; their vertices fill the cube's
vertex set, but their solids do not fill the cube.

The user's two-dimensional picture is exact in a rectangle: a central
isosceles triangle plus copies of its two altitude halves fill it. The
central triangle is equilateral at aspect sqrt3/2. An affine map turns
the rectangle into a square, retaining the dissection and area ratios
but changing the angles. The note proves that an uncut rigid
rearrangement of these equilateral pieces into a square is impossible
by their edge-length field Q(sqrt3).

The cube isometry `G(x,y,z)=(1-y,x,1-z)` has two four-cycles and preserves
vertex parity. An explicit table conjugates it with D17, identifying
quadratic character with parity. A vertical reflection supplies a
commuting sector-swapping involution. However no cube isometry has
order8, so multiplication by3 cannot be transported as an isometry
under any labelling. The finite dynamical model preserves D and the
sector label, while its extension of the action differs from arithmetic.

The three perfect matchings of a tetrahedron give the genuine quotient
S4/V4=S3. This is a finite action on three objects; no unique decreasing
parent or three-child generation rule follows from it.

## 4. A genuine elliptic ternary tree, with a sharp positivity boundary

The [elliptic lane](prime_shells_20260921_elliptic.md) uses the repaired
81/80/79-digit seed from the previous audit. The literal decimal strings
had an extra factor10 in the last two coordinates. The repaired triple
corresponds to H=9G, G=(-4,28), on

```
E: y²=x³+109x²+224x.
```

For any infinite-order rational H, the three fixed algebraic maps

```
P -> 3P-H, 3P, 3P+H
```

give an exact labelled ternary tree on `{nH:n>=1}`. Its integer
coefficient branches are n->3n-1,3n,3n+1; every n>=2 has unique parent
floor((n+1)/3). Match these three branches with the three Berggren
branches to obtain an explicit tree isomorphism through the common word.
This preserves root, adjacency, labels and depth, but not side lengths,
prime factors, odd-square shell or fruit positivity.

At H=9G all first children18G,27G,36G fail positivity. The positive point
17G is outside this subgroup half. More strongly, no single expanding
affine elliptic map P->mP+R with |m|>=2 preserves all positive rational
fruit solutions: its expanding circle action contradicts the proper
open real positivity region. Variable return rules are not ruled out.

Literal inverse tripling is a different structure. Keeping the torsion
coordinate proves the inverse tree of3^rG has1+3r vertices in <G,T>;
the root9G has7. The CITED full Mordell--Weil classification extends this
to all rational points. Dropping the torsion guard would incorrectly
replace its single continuing branch by a full infinite ternary tree.

The power-gap identity2^10=10³+24 gives the legitimate point(10,32) on
the separate curve y²=x³+24. It is not a gap-one instance of Catalan.
The two elliptic curves have different j invariants and different
good-reduction point counts at19, so they are neither geometrically
isomorphic nor Q-isogenous. The finite exponent scan is not an
all-exponent Diophantine classification.

## 5. What has moved, and what still needs proof

The session produces several exact maps rather than a shared numerical
motif: harmonic quotient to binomial quotient; shell to coset; shell
to cyclotomic cycle; prime valuation to cycle lift; shell permutation
to cube isometry; and Berggren branch word to a group-law elliptic tree.
Each map comes with its retained predicates and first lost coordinate.

For Collatz, the open obligation is still a rule that preserves actual
guarded edges and forces a well-founded decrease or universal root
arrival. The [earlier edge hostile](odd_square_20260921_edges.md) shows
why ordinary Berggren descent does not already do so. The new sector
and torsion obstructions suggest retaining a coset/phase or torsion
coordinate before attempting a quotient-based descent argument.

All four lanes have paper proofs and exact companions. Run
`python 04-computation/experiments/prime_shells_20260921_verify.py` for
fresh normal and optimized replays with frozen input hashes and output
sentinels. The [manifest](prime_shells_20260921_manifest.json) binds the
notes, programs, inherited inputs and certificates. External classification
claims retain their primary citations and are not inferred from a finite scan.
