# Odd-square triangle families, their cycles, and the legal Collatz subset

**2026-09-21. Status: PROVED scoped elementary results; FINITE-EXACT
controls; CITED inverse-orbit bound. Collatz remains OPEN.** No novelty
claim or new Lean-formalization claim. The user's annotation about integer
versus rational k is addressed explicitly below. The newly pasted density
argument is preserved as an [unaudited source](../reference/ODD-SQUARE-GLUED-LINE-2026-09-21-SOURCE.md).

## Inheritance and session board

Closest proved mechanism: [THM-3756, odd-square ordinal Berggren descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md)
already supplies both square roots, the totient fibre count and unique
tree descent. The [semicircle continuation](collatz_mod6_20260917_pythagorean_semicircle.md)
supplies shape coordinates and Gaussian doubling. The canonical hostile
is a primitive triangle whose legal Collatz interpretation disappears
after a legitimate triangle operation. The corrected near miss is to
forget the inner root or to treat an ambient spatial density as an orbit
law. The least-used sidecars are the sampling height and edge direction.

Anchor: the odd-square chart and its actual Collatz edge interpretation.
Niche: shape distributions under different arithmetic heights.
Wildcard: a finite triangle permutation conjugate to modular doubling,
with exact links to Fermat- and Mersenne-shaped numbers.

| Live concept | Genuine connection | Information requiring retention |
|---|---|---|
| Primitive triangles | Two coprime odd square roots | Both roots, not the outer root alone |
| Shape interval | Similar triangles / half-angle coordinates | Height and choice of probability measure |
| Fixed odd-square family | Folded doubling on units modulo u | Unit coset and possible separate cycles |
| Collatz edges | Coprime roots give a primitive triangle | Orientation, parameter and valuation guard |
| Inverse reachability | Actual exponent words and CRT | Height, sign, content and one fixed integer |
| Divisors / floors | Exponent boxes and boundary residues | Prime multiplicities and exceptional zero residues |

## 1. The odd square is exactly half an address

For every primitive triple with odd leg a, even leg b, and hypotenuse c,

```text
c+b=u^2, c-b=v^2,
u>v>0 odd, gcd(u,v)=1,
(a,b,c)=(uv,(u^2-v^2)/2,(u^2+v^2)/2).                 (S1)
```

The converse holds. This is the inherited THM-3756 chart, not a newly
discovered parametrization. The even leg is the one that works;
`c+odd_leg=2m^2` cannot be a square. At a fixed odd u>=3 there are exactly
`phi(u)/2` primitive triples. For instance u=5 gives both (5,12,13) and
(15,8,17). Thus indexing by odd squares is a useful stratification, but
does not uniquely index individual triangles.

**Annotation response.** The earlier statement about
`(k^2-1,2k,k^2+1)` was about similarity classes. Integer k gives a thin
subset. Even k gives a primitive triple directly; odd k gives common
factor two and must be reduced. After reduction its two families are
exactly `v=1` and `v=u-2`, with the 3-4-5 overlap counted once.
Allowing rational k>1 covers every primitive shape after clearing
denominators and removing the common factor. This agrees with the
odd-square chart: an integer address retains both roots, while their ratio
is the single shape parameter.
Details, exact count and equality cases:
[triangle lane](odd_square_20260921_triangles.md).

## 2. A genuine cycle operation on each family

Fix u and replace the inner odd root by

```text
D_u(v)=|u-2v|.                                        (S2)
```

This permutes the `phi(u)/2` allowed roots. It is multiplication by two
on `(Z/uZ)^*/{+1,-1}`, using the unique odd representative between 0
and u. Every cycle has the same length

```text
d_u=min{k>=1 : 2^k=+1 or -1 mod u}.                   (S3)
```

The side-length operation is especially concrete:

```text
(a,b,c) -> (|b+c-2a|, 2(a+b-c), 3c-b-2a).            (S4)
```

Its new even leg is four times the old inradius, and its new hypotenuse
plus even leg remains u^2. At u=7 this gives the genuine triangle cycle

```text
(7,24,25) -> (35,12,37) -> (21,20,29) -> (7,24,25).
```

The map `v -> 2cos(2pi v/u)` conjugates this finite permutation to
`z -> z^2-2` on its cyclotomic set. That angle is a modular coordinate,
not the actual triangle angle `2arctan(v/u)`.

There is a precise connection to the user's power-of-two examples:
`u=2^m+1` has d_u=m for every m>=1, since `2^m=-1 mod u` and smaller
positive powers cannot equal +/-1 modulo u by size. Likewise
`u=2^m-1` has d_u=m for m>=3. Primality is unnecessary. The number of
cycles is `phi(u)/(2m)` in either case. Hence u=63 has three six-cycles,
whereas u=65 has four six-cycles.

The integer-k boundary families lie in the orbit containing v=1, because
D_u(1)=u-2. That orbit first fails to cover a whole family at u=17:

```text
1 -> 15 -> 13 -> 9 -> 1;
3 -> 11 -> 5 -> 7 -> 3.
```

Thus a small generating pattern can be real and still miss a disjoint
component. The full proof and first-failure audit are in the triangle lane.

## 3. What the epsilon coordinates actually measure

For unit hypotenuse, with shorter projection e, longer projection d,
altitude l, and the user's smaller angle theta,

```text
e+d=1, l^2=ed,
r=e/d=tan^2(theta),
l=sqrt(r)/(1+r)=sin(2theta)/2.                         (S5)
```

These are three coordinates on one interval, not three independent
variables. The unit hypotenuse is the semicircle's diameter, so its radius
is one half. Simultaneous strict lower bounds by epsilon are equivalent,
for 0<epsilon<1/2, to

```text
theta>max(arctan(sqrt(epsilon)),arcsin(2epsilon)/2).    (S6)
```

There is a real crossover between which constraint binds, at
`epsilon(1+epsilon)^2=1`, namely epsilon=0.4655712319... . This crossover
is geometric; the coprimality constant enters a different calculation.

The relevant lattice densities are:

| Sampling universe | Coprime density / joint density |
|---|---:|
| All integer pairs | 6/pi^2 |
| Opposite-parity pairs, conditionally | 8/pi^2 |
| Odd-odd pairs, conditionally | 8/pi^2 |
| Coprime opposite-parity pairs among all pairs | 4/pi^2 |
| Coprime odd-odd pairs among all pairs | 2/pi^2 |

The factor at prime two is essential. In every fixed odd-square family,
coprimality has already been imposed, so all its triangles are primitive.

Two natural ways to sample all primitive triangles give different laws:

| Height / universe | Count | Limiting smaller-angle CDF |
|---|---|---|
| Hypotenuse c<=X | X/(2pi)+O(sqrt(X)log X) | 4theta/pi |
| Outer odd root u<=U | U^2/pi^2+O(U log U) | tan(theta/2)+1-tan(pi/4-theta/2) |
| One fixed odd-u family as u->infinity | phi(u)/2 | Same CDF as outer-root height |

The last limit holds through every sequence of odd u tending to infinity,
including highly composite u. The slope CDF differs from uniform by at
most `2^omega(u)/phi(u)<=sqrt(15)/(2sqrt(u))`.
The distinction comes from the radial weights: u<=U gives uniform v/u;
c<=X gives slope density `4/[pi(1+t^2)]`. Thus the specific probability
attached to an epsilon cutoff depends on a declared sampling rule.
All proofs and finite controls: [shape lane](odd_square_20260921_shapes.md).

## 4. The legal Collatz subset is primitive, thin, and angle-restricted

Every distinct positive odd U_(+/-1) edge x->y has gcd(x,y)=1. Therefore

```text
(xy, |x^2-y^2|/2, (x^2+y^2)/2)                        (S7)
```

is a primitive triple with outer and inner roots max(x,y),min(x,y).
Keep the direction: the exact remaining legality test is
`3x+b=2^k y` with actual k>=1.

This is a useful new-to-session bridge, but both natural triangle
operations fail to preserve it. The U_1 edge 7->11 gives (77,36,85);
its Berggren parent (21,20,29) has roots 7,3, which form no U_1 edge.
The U_(-1) edge 5->7 gives (35,12,37); applying (S2) gives the same
illegal parent-root pair 7,3. A descent theorem on all triangles therefore
does not descend automatically to legal Collatz edges.

There is also a sharp geometric exclusion:

```text
U_1 edges:    theta<=arctan(65/72)=42.0750... degrees,
             equality exactly at13->5, triangle(65,72,97);
U_(-1) edges: theta<arctan(48/55)=41.1121... degrees,
             supremum approached along x=3 mod16.
```

These bounds are proved by all four halving regimes k=1,2,3,>=4.
The near-isosceles end of the user's interval does not occur here.

There are only
`C_edge sqrt(X)+O(log X)` directed legal edges with encoded hypotenuse
at most X, where

```text
C_edge=sum_(k>=1) 1/sqrt(2(4^k+9))=0.5078191557... .
```

Their triangle image has density zero among all primitive triangles.
Uniform directed-edge sampling under this height has an atomic limiting
angle law, concentrated at the shapes of the branch ratios 3/2^k, rather
than the ambient continuous uniform law.
This is a proved sampling mismatch, not merely a warning that independence
might fail. [Edge proofs and exact certificate](odd_square_20260921_edges.md).

## 5. Inverse reachability and the minus-five extension

Reversing every arrow and asking whether 1 reaches every positive integer
is equivalent to the usual Collatz conjecture. Negating integers instead
changes the parameter sign while preserving forward arrows; doing both
operations is not the same as doing neither.

The three residue rows are useful exact local coordinates, but their
connectivity alone omits height and the successive valuation guards.
Krasikov--Lagarias explicitly retain a height-constrained counting
function as well as residue classes. **CITED:** their
[Theorem 6.1](https://arxiv.org/pdf/math/0205002) proves at least x^0.84
ancestors below x for each fixed positive target prime to three, for
sufficiently large x. Its counting universe is all positive starting
integers for the shortcut map T(n)=(3n+1)/2 on odds and n/2 on evens,
not an explicitly odd-only count. The primary text's definitions in section2 and
Theorem6.1/its proof were read. This historical bound does not assert all
integers are ancestors; no claim that it is the latest bound is made here.

For a general odd parameter b, every edge satisfies
`gcd(x,y)=gcd(x,b)`. For b=-5, content is either1 or5 and remains constant.
Content5 divides to the b=-1 dynamics and gives triangle content25.
Content1 gives primitive triangles when the endpoint magnitudes differ;
the portal1->-1 gives a zero leg and is excluded. Sign changes are retained.
The [inverse lane](odd_square_20260921_inverse.md) proves an exact
two-exponent completion theorem with simultaneous 2-, 3- and 5-adic
constraints, including arbitrary powers of five. Disjoint negative basins
can each realize every compatible finite pattern. This strengthens the
reason local coverage cannot prove that only one global basin exists.

## 6. Audit of the new pasted argument and the other repeated connections

| Pasted / repeated claim | Verdict and strongest survivor |
|---|---|
| Dividing once by2 leaves only odd states | False: n=1 goes to2. Use full odd acceleration consistently. |
| The two nonzero mod3 rows must alternate | False: U_1(11)=17 stays5mod6; U_1(1)=1 stays1mod6. Odd multiples of3 really are source-only under U_1. |
| The sign-reflected three cycles are another independent family | They are conjugate copies; known odd periods1,2,7. Completeness is open. Reflection does not reverse time or turn divergence into attraction. |
| F=S+U is a universal conservation law | False at N=4: F=1,S=1,U=1. Its complete N>=2 solution shapes are the inherited p,p^3,p^2qr. |
| Prime-factor means are assigned to log and loglog as entropy derivatives | No defined derivative or entropy is supplied; the prime-factor mean has loglog scale, while divisor counts have log scale under uniform integer sampling. These are separate means, not orbit laws. |
| Strong halving requires a smooth odd factor | False: 1077->101 has exponent5 and prime odd core101. Only the actual two-adic valuation controls that division. |
| Spatial squarefree density forces every orbit to descend | Unsupported; no trajectory measure or descent inequality is given. The exact edge measure above already differs from the ambient primitive-point law. |
| The three floor identities | Inherited all-modulus repair: first two hold iff squarefree; the product sum holds iff prime. Boundary-zero counts explain the difference. |
| p,p^3,p^2qr and almost-prime classes | The exponent profile, not Omega alone, controls divisor balance. p^3,p^2q,pqr share Omega3 but differ in balance. |
| {7,21}, Bott periodicity, octonion dimensions | Existing tournament exclusions and Fano/squareclass maps remain valid in their own domains. No defined map turns their shared counts into Collatz descent; this session adds none. |

The all-modulus floor formulas are in the
[floor reciprocity note](arithmetic_braids2_20260917_floor_reciprocity.md),
the divisor classification in the
[divisor note](arithmetic_braids_20260917_divisors.md), and the detailed
scope of the remaining cross-domain suggestions in
[the earlier map audit](collatz_mod6_20260917_wild_typing.md).
These routes retain the useful proved statements without treating a
similar small number as a transported theorem.

## Reproduction and next proof obligations

```text
python 04-computation/experiments/odd_square_20260921_verify.py
```

The [master verifier](../../04-computation/experiments/odd_square_20260921_verify.py)
replays all four lanes normally and with -O, rejects stale outputs, checks
embedded source hashes when supplied, freezes all notes and sources, and
pins the resulting [manifest](odd_square_20260921_manifest.json).
All infinite results are written proofs; finite certificates are controls
over their explicitly stated universes. New general Lean formalization
is not claimed.

The next structural target is precise: find an operation on the encoded
triangles that preserves `3x+b=2^k y` while decreasing a well-founded
quantity, or prove a controlled legal path to such a decrease. Berggren
descent and fixed-root doubling fail that preservation test. Another
target is a height-sensitive inverse certificate that covers each integer,
stronger than realizing each finite residue pattern somewhere. Those are
still open obligations, not consequences already contained in the density.
