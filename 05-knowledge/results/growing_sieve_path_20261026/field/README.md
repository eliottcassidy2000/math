# The sixteen-state cover transports tournaments, but not arithmetic quality

**Status: PROVED elementary deductions and explicit counterexamples;
FINITE-EXACT finite-automaton and permutation controls; independently audited PASS.**
[Independent audit](../audit_field/AUDIT.md). This continues the
audited [paired-design packet](../../twin_pythagorean_124_20261026/finite_field/README.md).
Its determinant bit and 16 reachable point-and-bit states are inherited.
The new result identifies that observer with a homogeneous space of **actual
oriented tournaments**, gives its exact prime-wheel refinement, and exhibits
arithmetic information that does not descend through it. No general abc
bound, twin-prime infinitude, or global Hamilton-spectrum conclusion follows.

## 1. Inheritance and the six live objects

Closest proved mechanism: the paired `S(3,4,8)` designs are preserved by the
outer Berggren branches U,D and exchanged by the middle branch A. The
canonical hostile is UU versus A: equal conic endpoint, opposite determinant
bits. The corrected near miss is replacing PGL2(7) by PSL2(7) for all three
branches. The least-used sidecar now is the **oriented tournament on the
complement of a marked pole**.

The board is: branch words; the homogeneous 16-state cover; marked Paley
tournaments; prime-wheel residues; radical/height data; and pointed monoid
evaluations. The anchor is the exact finite observer; the niche is its
arithmetic failure boundary; the wildcard is a genuine F1-language diagram
of distinct monoid evaluations.

Recent inputs change the question rather than supply interchangeable bits:

- [Strong-atom phase and joint degree](../../strong_atom_continuation_vitali_20261026/FROBENIUS-ABC.md)
  retains the ancestral owner map and the joint field coordinates `(B,C)`.
  Actual path complement is not field Frobenius.
- [Fresh re-marking and divisor arrangements](../../fresh_remarking_quotient_and_divisor_arrangements_20261026.md)
  treats operation-stable refinements and weighted divisor cells. Its fresh
  re-marking minimality is OPEN beyond its finite universe; we do not use it
  as a theorem here. Its divisor weights and augmentation are distinguished
  explicitly in §5.
- [THM-4553, point-stripping ladder](../../../../01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md)
  already supplies the Paley order-21 pole stabilizer and its two arc orbitals.
  The following homogeneous-space and bracket formulas are an explicit
  continuation of those inherited facts, not a new Paley classification.

## 2. An arc-level realization of the two-sheet observer

Let `X=P1(F7)`, `G=PGL2(F7)`, and let `chi(g)` be the Legendre symbol of its
determinant. Write `chi(g)=(-1)^b(g)`, with b in F2. The common Euclid slope
is `r=n/m`; as in the prior packet we use the explicitly reciprocal auxiliary
coordinate `t=m/n=1/r` for formulas. A pole is a point of X; it is not
necessarily the integer tree root.

Take `v_t=(t,1)` and `v_infinity=(1,0)`, and write `[x,y]=det(v_x,v_y)`.
For three distinct points define

    K(p,u,v)=[p,u][u,v][v,p].                            (1)

For each `(p,e)` in `X x F2`, orient a tournament T_(p,e) on `X minus {p}` by

    u -> v  iff  Legendre(K(p,u,v))=(-1)^e.              (2)

**F1, the exact transport theorem.** These are 16 distinct marked-pole
tournaments, two converse orientations at each pole. The action is

    g : (p,e) -> (gp,e+b(g)),                            (3)

and it transports every arc in (2). This G-set is `G/H21`, where H21 is the
square-determinant stabilizer of infinity. The forgetful map to the eight
poles is precisely the degree-two homogeneous cover

    G/H21 -> G/B42,                                    (4)

where B42 is the full pole stabilizer.

*Proof.* Rescaling the representative at any one of p,u,v multiplies (1)
by a square. Interchanging u,v negates (1), and `Legendre(-1)=-1`, so (2)
is an honest tournament. Applying a matrix g multiplies each bracket by
det(g), before the individual representative rescalings. Thus the square
class of K is multiplied by det(g)^3, proving (3).

At p=infinity, (1) is `v-u`; the e=0 tournament is Paley7 and e=1 its
converse. The pole stabilizer consists of `t->a t+b`, a nonzero. Its kernel
under chi has a in `{1,2,4}`, so its orders are 42 and 21. Its two chi fibres
give the two orientations. Transitivity of G on poles now shows there are
exactly 16 marked-pole objects and that the stabilizer of `(infinity,0)` is
exactly H21. The coset map is `gH21 -> (g infinity,b(g))`. Its projection
forgets b and is (4). QED.

The inherited U matrix is `[2,-1;1,0]`, so U sends infinity to the actual
root coordinate t=2 (common slope r=4) and has square determinant. Therefore
the branch observer based at the root is `G/(U H21 U^-1)`. This conjugation
is required: the actual root and the distinguished pole have different
coordinates. In common-r coordinates the pole infinity in t is r=0.

The two outer branches generate PSL2(7); all three generate G, as proved in
the prior packet. Consequently the branch semigroup reaches all 16 objects.
The sheet e is exactly the parity of the number of A letters in the word.
This now transports literal tournament arcs, not just an abstract bit.

At the pole infinity, H21 acts regularly on the 21 arcs of Paley7: the unique
affine map taking `(0,1)` to an arc `(u,v)` has multiplier `v-u`, a square.
It also acts regularly on the 21 reversed ordered pairs. Hence its invariant
tournaments are exactly Paley7 and its converse. Every object in (2) has

    H=189,  |Aut|=21,  tc=H/|Aut|=9.                    (5)

The checker independently computes H for all 16 arc relations and exhausts
the 7! automorphisms at the base pole. Thus (1)--(4) genuinely connect the
number seven of vertices and the order-21 symmetry to the branch observer.
They do not supply a mechanism for the two forbidden Hamilton counts.

## 3. Why tournament phases do not become this determinant bit

Every finite tournament automorphism group has odd order. Otherwise Cauchy's
theorem supplies an involution; its vertex permutation contains a transposed
pair, reversing the unique arc between that pair, a contradiction. Therefore
every homomorphism from a tournament automorphism group to C2 is trivial.
In particular chi restricts trivially to H21. The nonsquare half of B42
contains anti-isomorphisms between the two orientations, not automorphisms
of one oriented tournament.

This gives an operation-level separation from the incoming strong-atom
joint-degree invariant. Under the Session69 marking, the pole-preserving
Frobenius is `t->2t`, whose determinant is a square and whose order is three.
For `(B,C)` in F8² the stabilizer among its powers has order three when both
coordinates lie in F2 and order one otherwise. For example `(0,0)` and
`(0,alpha)` have different joint degrees while the marginal B is identical.
All these Frobenius powers have determinant bit zero. The sign character
therefore cannot recover this joint degree or the corresponding residual
stabilizer. Retaining that incoming coordinate is necessary when its
predicate is the target; adding one bit does not universally close every
continuation interface.

Raw H and normalized tc retain their separate scopes. The all-order raw-H
exclusions are inherited from [THM-1370, H-spectrum omits 7 and 21](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The tc exclusions through eleven vertices are [THM-4790](../../../../01-canon/theorems/THM-4790-presentation-fibres-7-and-21-absent-through-eleven-vertices.md);
their unrestricted form remains [HYP-9352](../../../hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md).
Equation (5) is one transported family, not a proof of either exclusion.

## 4. Exact prime-wheel refinement and two arithmetic failures

For any positive modulus M, add the residue vector `(m,n) mod M` to the
point-and-bit observer. The state is

    (t mod7, e, (m,n) mod M),                            (6)

with the projective infinity convention for the first coordinate. The allowed
residue vectors lie in the ambient carrier `gcd(m,n,M)=1`; this condition
alone is not a characterization of reachability at arbitrary M. For example,
when 2 divides M the source's opposite-parity condition must also be retained.
Each branch applies its integer
matrix to the vector, its projectivity to t, and toggles e only for A.
Because every branch determinant is +/-1, this is a closed finite transition
system; forgetting any residue component commutes with every branch. When
7 divides M, retain the consistency of the two mod7 coordinates; arbitrary
cartesian states need not be realizable. This is a forward observer, not an
injective encoding of the integer tree.

On the D-ray, `(m,n)=(2+2k,1)` and e=0 for every k. For an odd squarefree
wheel `M=product q`, excluding divisibility of either candidate `m-1,m+1`
is exactly `gcd(m²-1,M)=1`, with `product(q-2)` admissible residues. The
exception that a candidate can equal a wheel prime must be retained when
using this as a primality test. Beyond the largest wheel prime the test is
necessary, and still not sufficient. The prior [CRT obstruction](../../twin_pythagorean_124_20261026/residue_obstruction/README.md)
already proves that a fixed wheel cannot certify twins; we do not repeat it.

**FINITE-EXACT refinement at M=15.** There are 192 primitive residue pairs,
and all `16*192=3072` joint states are reached. The D-ray has period 105 in
(6), with 21 wheel-surviving states; its allowed m residues are `{0,3,12}`
mod15. This last 21 is a wheel count, with no asserted bijection to the
21 arcs above.

For the arithmetic test use the actual primitive abc packet

    A=a², B=b², C=c²,  where (a,b,c)=Phi(m,n).             (7)

Its radical is `R=rad(abc)`. The statement that abc quality exceeds one is
exactly `C>R`, so the following checks use integer inequalities, not floating
point logarithms.

**Same full finite observer, opposite abc-hit predicates.**

| word | (m,n) | actual twins m−1,m+1 | (a,b,c) | R | C |
|---|---|---|---|---:|---:|
|D^5|(12,1)|(11,13)|(143,24,145)|124410|21025|
|D^215|(432,1)|(431,433)|(186623,864,186625)|8358844170|34828890625|

Both states have t=5, e=0, and `(m,n)=(12,1)` mod15, so they give the same
object in (6), survive its wheel, and even have the same **true** twin-prime
predicate. Yet the first has C<R and the second C>R. Thus this observer,
even supplemented by the twin predicate, does not determine the abc-hit
predicate. Its missing coordinates include valuations and Archimedean height.
This is a counterexample to determinacy, not a disproof of any abc inequality.

**Radical does not determine the determinant bit.** The exact words

    UUUDUDDDU  -> (187,100),  Phi=(24969,37400,44969),  e=0;
    DUUDDDDDUA -> (386,153),  Phi=(125587,118116,172405), e=1

have the identical radical

    R=2099693897070=2*3*5*7*11*17*29*41*193*233.           (8)

The checker records all term factorizations. No minimality of this witness
is claimed. The first failed implication is descent of e through the radical
observation. Conversely, even e alone permits both C<R and C>R: the words
empty, UU, A, ADU exhibit all four `(e, C>R)` possibilities. This rules out
the displayed determinacy claims, not every conceivable parity-dependent
bound on quality.

## 5. A literal monoid diagram for the F1 language

Here F1 language means explicit pointed monoids and their integral monoid
algebras; no ordinary field of cardinality one is assumed. The inherited
[marked Laurent-lift note](../../normseven_arithmetic_f1_20261026/f1/RESULT.md)
already warns that torus characters, positive monoids and affine boundaries
retain different information. Our maps are simpler and completely explicit.

Let W be the free branch-word monoid. Counting U,A,D gives its commutative
monoid quotient `W -> N³`. After adjoining zero and base extending to Z,
the contracted monoid algebra is `Z[x_U,x_A,x_D]` (the absorbing monoid zero
maps to ring zero). The branch determinant is exactly the
monomial evaluation

    word -> x_U^(#U) x_A^(#A) x_D^(#D) -> (-1)^(#A),      (9)

where the last map evaluates at `(1,-1,1)`. Ordinary augmentation evaluates
at `(1,1,1)` instead. It sends U and A to the same value 1, while (9) sends
them to opposite signs; hence the determinant character does not descend
through augmentation. Commutativization retains (9) but loses the conic
endpoint: UA and AU give `(8,3)` and `(8,5)`, with t=5 and t=3 modulo7.

For a squarefree integer `N=p1...pk`, the analogous proper-divisor polynomial
is

    P_N(x)=product_i(1+x_i)-product_i x_i.                (10)

Augmentation gives `2^k-1`; evaluation at the actual primes gives s(N);
evaluation at all −1 gives `(-1)^(k+1)`, for k>=1. These are different maps
from the same explicit algebra. This recovers exactly the count/weight/sign
distinction in the incoming divisor-arrangement note. No assignment of
branch letters to arithmetic primes is asserted by (9)--(10).

There is also a strict obstruction. Radical support is the idempotent monoid
of finite prime subsets under union: `S union S=S`. Any monoid homomorphism
from it to a group must be trivial, since its image y obeys y²=y. Thus a
nontrivial group character such as (9) cannot descend multiplicatively
through this radical-support monoid. The elementary sign `(-1)^|S|` is not
a character for arbitrary union: overlapping supports defeat it. It behaves
multiplicatively on disjoint supports, a different operation domain. This
abstract no-go and the actual primitive abc witness (8) cover distinct
claims: respectively homomorphic descent and arbitrary functional descent
of this particular bit.

## 6. Reproduction and scope

Run in this directory:

```sh
python3 -I -B check.py > replay.json
python3 -I -B -O check.py > replay.optimized.json
cmp replay.json normal.json
cmp replay.json replay.optimized.json
```

The checker enumerates every projective matrix action over F7, all 16 arc
relations and their 336 transports, the complete base 7! automorphism oracle,
all mod15 primitive-pair states and transitions, the full reachable observer,
and the stated integer factorizations. Source pins separate inherited proofs,
finite evidence and incoming open statements. All new mathematical claims are
elementary proofs or specified finite universes. The surviving positive
bridge is (1)--(4); the additional wheel is exact but loses the predicates
demonstrated above. Neither provides a reduction among the famous open
problems merely because the same small numbers occur.
