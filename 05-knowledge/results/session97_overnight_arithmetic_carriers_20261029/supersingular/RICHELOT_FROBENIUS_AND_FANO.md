# The fifteen kernels, the Frobenius guard, and the thirty missing Fano carriers

**Session97, 2026-10-29. Status:** CITED standard arithmetic identifications; PROVED elementary finite-geometric maps and matching recurrence; FINITE-EXACT complete small-space certificates. No new moonshine theorem, Sha-finiteness theorem, or tournament Hamilton-spectrum theorem is claimed. Source snapshot: `5d6eb91dd76eee86b91e04d856daa32c32b3d61d`.

## Results first

1. The actual bridge behind one occurrence of **15** is
   six branch points → fifteen perfect matchings → fifteen Lagrangians in
   `J(C)[2]` → fifteen geometric Richelot kernels. Rational kernels are the
   **Frobenius-fixed matchings**, not all geometric matchings.
2. At the Ogg prime **11**, two explicit superspecial products have **15**
   versus **3** rational kernels. Supersingular coarse moduli being rational
   does not make their attached kernel data rational.
3. At eight branch points the matching recipe is incomplete:
   **135 = 105 matching + 30 extended-Hamming Lagrangians**.
   Puncturing the latter at a named point gives precisely the thirty labelled
   Fano planes. This is an explicit code/geometry bridge, not a shared count.
4. For Frobenius cycle type `(1,7)`, **zero** matching kernels but **two**
   nonmatching kernels survive. An arithmetic witness is the genus-three
   curve `y²=x⁷+x+1` over `F5`. Its two kernels correspond to the two cyclic
   Fano planes. Their cyclic-subgroup normalizers have order **21**, with
   multipliers `{1,2,4}`. This explains a genuine local 7/21 connection;
   it does not identify these integers with forbidden Hamilton-path counts.

## Inheritance and concept board

Closest mechanism: Session96 `f1/SYMPLECTIC_SKELETON_AND_GLUE.md`, particularly
`L2(2)=15=9+6` and its source-preserving graph-of-an-isomorphism proof.
Canonical hostile inherited from Session32: rational supersingular j-values
can have nonrational level-two edges. Its concrete four-prime certificate is
in `session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md`.
That note already uses Aricheta's modular rationality theorem; this session
must not rediscover it as new.
Least-used relevant sidecar: the **Frobenius action on marked two-torsion**,
not only its dimension or the field of the coarse j-invariant.

The code connection is inherited from Session85
`FANO_PALEY_INFORMATION.md`: its labelled Fano/Paley bridge needs a field
chart, dual line normals and a reversal origin. Here a new source for the
same thirty labelled Fano planes is supplied. No orientation is automatically
created by a Lagrangian. The live board is:
coarse moduli / symplectic kernels / Frobenius / quadratic refinements /
branch-point matchings / Fano codes.

## 1. Which fifteen, and which arithmetic bridge?

The fifteen Ogg primes are

`2,3,5,7,11,13,17,19,23,29,31,41,47,59,71`.

The relevant property is that **all supersingular j-invariants in characteristic
p** lie in `Fp`; equivalently `X0+(p)` has genus zero. These are also precisely
the prime divisors of the Monster order. This is not supersingularity of one
fixed elliptic curve. Aricheta extends the rationality question to level
structures and Monster centralizers; see Theorems 1.1 and 1.3 in
[Supersingular Elliptic Curves and Moonshine](https://arxiv.org/html/1809.07421v2).
There is no canonical map found from that set of primes to a fifteen-element
Lagrangian set. Their cardinalities alone cannot transport an action or pairing.

There is instead a direct arithmetic passage: products of supersingular
elliptic curves are superspecial abelian surfaces. Their principally polarized
2-isogeny kernels are Lagrangians. One must retain polarization, marked kernel,
Frobenius and source automorphism weights. [Katsura–Takashima](https://arxiv.org/abs/2003.00633)
classify decomposed Richelot isogenies, including automorphism identifications.
Fifteen geometric kernels do not mean fifteen nonisomorphic targets.

For a smooth genus-two curve with six distinct geometric branch points in
characteristic different from two, its two-torsion is

`V6 = {even subsets of the six points}/(S ~ S complement)`.

Addition is symmetric difference and the symplectic form is
`b([S],[T])=|S intersection T| mod2`. This is a four-dimensional space.
The Galois-equivariant subset description and pairing are standard; see
[Creutz–Viray, Two torsion in the Brauer group of a hyperelliptic curve](https://sites.math.washington.edu/~bviray/papers/CV_BrauerHyperellipticCurves.pdf), §2.
A matching `{P1,P2,P3}` maps to the span of its three pair classes. Their sum
is zero, any two are independent, and they pair trivially. Conversely the
three nonzero vectors of a Lagrangian have unique pair representatives;
orthogonality forces these pairs to be disjoint. This proves the bijection.

Richelot realizes these kernels arithmetically. A matching gives a splitting
`f=F1 F2 F3` into quadratics, allowed to be permuted by the ground-field
Galois group. The kernel descends exactly when the unordered splitting does.
The determinant of their coefficient rows distinguishes a Jacobian codomain
from a product codomain: zero gives the latter. Thus neither a chosen order
of factors nor nondegeneracy should be silently imposed on the fifteen.
See [Florit–Smith, An atlas of the Richelot isogeny graph](https://arxiv.org/pdf/2101.00917), §§1–2.

## 2. An all-size fixed-matching recurrence

Let a permutation have `m_l` cycles of length `l`. An invariant perfect
matching pairs two distinct cycles only when their lengths agree: the mate
map commutes with the permutation and is a bijection between those cycles.
There are exactly `l` compatible offsets. A cycle can be matched internally
only if `l` is even, by its unique antipodal pairing.

Consequently its number of invariant matchings is `product_l F_l(m_l)`, where

```
F_l(0)=1;
F_l(1)=0 for l odd, 1 for l even;
F_l(m)=l(m-1) F_l(m-2)                         (l odd),
F_l(m)=F_l(m-1)+l(m-1) F_l(m-2)               (l even).
```

Equivalently, odd l requires even m and gives
`m! l^(m/2)/(2^(m/2)(m/2)!)`; even l gives
`sum_(j=0)^floor(m/2) m! l^j/((m-2j)! 2^j j!)`.
The recurrence preserves cycle sources and phase offsets rather than replacing
them by counts of branch points. For an involution with t two-cycles and s
fixed points, multiply `(s-1)!!` (zero for odd s) by `f_t`, where
`f_t=f_(t-1)+2(t-1)f_(t-2)`. The values start `1,1,3,7,25`.

For six points this gives the complete rational-kernel table:

| Frobenius cycle type | Fixed kernels |
|---|---:|
| 1^6 | 15 |
| 2 1^4 | 3 |
| 3 1^3 | 0 |
| 2^2 1^2 | 3 |
| 4 1^2 | 1 |
| 3 2 1 | 0 |
| 5 1 | 0 |
| 2^3 | 7 |
| 4 2 | 1 |
| 3^2 | 3 |
| 6 | 1 |

In particular **all fifteen kernels rational iff all six branch points are
rational**, for genus-two curves over finite fields of odd characteristic.
This is about stable subgroups, not each element of each subgroup being
rational; the latter is stronger for an individual kernel.

The seven matchings fixed by `(12)(34)(56)` consist of the original matching
and two crossed replacements for each of its three choices of a retained
pair. Their shared-pair graph is three triangles with a common center:
seven vertices, nine edges, twelve nonedges among twenty-one unordered pairs.
This is an explicit 15→7 mechanism, but its undirected 21 pairs are not a
Hamilton-path count.

## 3. Product Frobenius intertwiners and the p=11 hostile

For `A=E×E'` with product polarization and both factors defined over the ground
field, write `U=E[2]`, `W=E'[2]`. Session96's nine separable Lagrangians are
products of lines. Its six mixed ones are graphs of symplectic isomorphisms
`a:U→W`; over F2 symplectic and antisymplectic agree.

If Frobenius acts by `u,v` respectively, the exact number of rational kernels is

`(# fixed lines of u)(# fixed lines of v) + #{a in Sp2(F2): va=au}`.

The second term vanishes unless u and v are conjugate, then equals the
centralizer order. `Sp2(F2)=S3` has types I,T,C (identity, order2, order3).
The table entries below are `(separable,mixed,total)`:

| u \\ v | I | T | C |
|---|---|---|---|
| I | (9,6,15) | (3,0,3) | (0,0,0) |
| T | (3,0,3) | (1,2,3) | (0,0,0) |
| C | (0,0,0) | (0,0,0) | (0,3,3) |

This is not a table for a product whose factors are exchanged by Frobenius;
that requires the exchange map as additional data.

Over F11, set `E1:y²=x³-x`, `E2:y²=x³+1`. Exhaustive point enumeration gives
12 points on each. The standard elliptic trace-zero criterion for p>3 makes
both supersingular. The first cubic splits completely; the second has exactly
one rational root and its residual quadratic is irreducible. Frobenius types
are I and T. Thus `E1²` has fifteen rational kernels but `E2²` has only three,
despite **11 being an Ogg prime and both surfaces being superspecial**.
This is the precise missing-sidecar hostile requested by the numerical
fifteen analogy. Over the quadratic extension both types become identity.

## 4. Eight branch points: the matching representation ceases to be complete

Put `V8=even subsets of eight points / complement`, dimension six. There are
`(2+1)(4+1)(8+1)=135` Lagrangians, whereas perfect matchings give `7!!=105`.
The matching map is still injective: the four pair vectors recover the matching.

**Dichotomy.** If a Lagrangian contains a pair vector e, pass to
`e perpendicular / <e>`. Deleting that pair identifies this quotient with
V6. Its induced Lagrangian is a matching plane by §1, and its full preimage
recovers a matching on all eight points. Consequently every other Lagrangian
contains no pair vector. Its seven nonzero vectors are represented by
four-subsets. The full preimage in the even-weight binary code of length eight
therefore has enumerator

`1 + 14 z^4 + z^8`.

It is a self-dual doubly-even `[8,4,4]` code, i.e. an extended Hamming code.
This description proves the split without treating it as a numerical deficit.

**Quadratic guard.** `q([S])=|S|/2 mod2` is well-defined on V8 because
complementing changes it by `4-|S|`, an even integer. Its polar form is b.
Exactly the thirty nonmatching Lagrangians are totally singular for q; each
of the other 105 has nonzero linear restriction of q. On V6 the analogous
formula is NOT well-defined: complementing changes it by `3-|S|`, odd.
A separate quadratic refinement is needed there. In particular the fifteen
symplectic planes of dimension four over F2 must not be confused with the
six totally singular planes for a split quadratic form, as in the distinct
BKLPR local model. No Selmer probability follows from our symplectic count.

Choose a named eighth point infinity. Among the fourteen weight-four words,
exactly seven contain infinity. Remove infinity to get seven triples. Any
pair of these triples meets in exactly one point: their corresponding words
intersect evenly, and cannot coincide or be complementary. The triples are
therefore a Fano plane on the remaining seven points. Conversely a Fano plane
extends to its `[8,4,4]` code by adjoining infinity to each line and taking
complements, zero and the all-one word. These maps are inverse. Thus the
thirty exotic Lagrangians map **bijectively and S7-equivariantly** to the thirty
labelled Fano planes once infinity is marked. Forgetting the mark loses the
chosen puncture map, not merely a count.

## 5. A genuine 7/21 arithmetic example

The polynomial `x^7+x+1` is irreducible over F5. The checker certifies
`x^(5^7)=x mod f` and `gcd(f,x^5-x)=1`; since seven is prime these are the full
irreducibility conditions. Hence the smooth hyperelliptic curve

`C/F5: y²=x⁷+x+1`

has Frobenius type `(1,7)` on its eight branch points (the fixed point is
infinity). No perfect matching is Frobenius-invariant by §2. Yet it has
exactly two stable Lagrangians, both exotic. Number the seven finite roots by
successive Frobenius conjugates. Their Fano planes are

`dev{0,1,3}` and `dev{0,1,5}` in Z/7.

The exact census independently verifies these two are all the fixed planes.
A group-theoretic check explains the count: a Fano plane's automorphism group
has order168 and eight Sylow7 groups; counting incidences with the120 cyclic
order7 subgroups of S7 gives `30*8/120=2` fixed planes. Inside the affine
normalizer `x→ax+b`, each plane is preserved exactly for `a∈{1,2,4}` and any
b, giving order `3*7=21`. Inversion swaps the two planes. The same chart and
multiplier group occur in the repo's explicitly oriented Paley7 construction.
The chart and the choice between the two planes remain necessary; no canonical
tournament orientation has been extracted from an unmarked curve.

These are rational subgroup schemes defining two polarized `(2,2,2)` isogenies.
We have not computed equations of the codomains, proved that they are
hyperelliptic Jacobians, or claimed supersingularity of this genus-three curve.
The example is a hostile to the matching-only carrier, not a new high-rank curve.

## 6. The other six-set and an outer-automorphism certificate

The full matching intersection graph has parameters `(15,6,1,3)`.
A **pentad** is a set of five matchings partitioning all fifteen edges of K6.
There are six pentads, and each matching belongs to two. Thus a matching maps
to an unordered pair in this *new* six-set. The map is bijective; matchings
share a pair iff their corresponding new duads are disjoint. This retains
incidence and gives a second six-point representation. The S6 action on
pentads is faithful and a transposition acts as three disjoint transpositions,
so it realizes an outer rather than inner automorphism after a labelling of
the new six-set. The small complete census certifies these assertions.
There is still no natural action of S6 on the fifteen Ogg primes supplied by
this construction.

## 7. Boundaries with previous six-root geometry

The prior `six_seven_gale_cube_20261025/BRANCH_CUBE_AND_DOUBLE_COVER.md`
constructs K3 surfaces from six **linear forms on P² and simultaneous square
roots**. Its branch curves can have genus five. They are not the six branch
points of a genus-two double cover of P¹. Consequently matching the number
six does not produce a Richelot map on those K3 surfaces. A valid future bridge
would require an exhibited genus-two Jacobian, Kummer presentation, or
fibration with an explicit compatible six-point branch divisor. None is
asserted here. The seven-root invariants `p_g=7`, `K²=16` do not by themselves
supply such data.

## Reproduction and audit scope

Run `python check_richelot_carrier.py` and also with `python -O`.
The output is identical and all checks use explicit exceptions, not assertions.
The universe includes all720 S6 permutations, all135 Lagrangians in V8,
all22 S8 cycle types, all six K6 one-factorizations, all30 Fano punctures,
the complete product-intertwiner table, and the two finite-field arithmetic
witnesses. The retained output lists the complete S8 table, including
`(2^4):25+14=39`, `(4²):5+6=11`, and `(1,7):0+2=2`.
The checker does not verify external theorems about Jacobian moduli or
superspecial isogeny graph connectedness; those are cited and distinguished.
