# Arithmetic seams: pointed products, Fermat path counts, and three-to-six lifts

**Date: 2026-09-21. Status: PROVED scoped statements; CITED tournament
identity; FINITE-EXACT controls; VERIFIED Lean arithmetic.**
Collatz convergence, a universal explanation of prime exceptions, and the
intended Anthropic reference remain OPEN or UNRESOLVED. No priority claim.

The productive common theme is a precise operation together with the
coordinate a quotient discards. The present session finds actual bridges:
Fermat numbers count Hamiltonian paths in an explicit tournament family;
Fermat primes control divisibility in the repeated-three family; a signed
lift of any rational three-point cycle carries `S_2 x S_3`; and adjoining a
distinguished root makes the proposed arithmetic product coherent.

## 1. Inheritance, portfolio, and proof boundaries

The closest proved mechanisms are the operation fibres in
[THM-2422, operation-fibres-summand-closure-and-twin-center-ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md),
the trace-one lift in
[THM-4146, rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md),
and the [inherited Mersenne/Fermat primitive-divisor comparison](collatz_mod6_20260917_zsigmondy_triad.md).
The canonical hostiles are a composite with entirely new prime factors,
two tournaments with equal scores but different path counts, and a
projective three-cycle whose vector representative has period six.
Corrected near misses concern suppressed roots, prime-power depth, and
identifying an incidence array with a planar embedding.

The underused coordinates are the first-hit offset of an affine residue
orbit, the overlaps of directed cycles, the central sign of a lift, and
the remainder omitted from a multiplication incidence. These formed the
live board throughout the session. The anchor was the user's arithmetic
and tournament synthesis; the niche was primitive-factor return ranks;
the wildcard was the reflection extending the rational-cycle lift. No
unrelated LRC or E8 result is claimed from matching small integers.

| Lane | Current proof and exact-control record |
|---|---|
| Repeated digits and prime divisors | [Prime-pattern note](arithmetic_seams_20260921_primes.md) |
| Tournament addition, multiplication, and path counts | [Tournament note](arithmetic_seams_20260921_tournaments.md) |
| Rational quadratic lifts and Chebyshev cycles | [Dynamics note](arithmetic_seams_20260921_dynamics.md) |
| Shifted arithmetic, graph complements, and plane addresses | [Operations note](arithmetic_seams_20260921_operations.md) |
| Machine-checked arithmetic | [Standalone Lean package](../../04-computation/lean/ArithmeticSeams/README.md) |

## 2. The tournament construction gives Fermat numbers directly

On vertices `0,...,n-1`, orient every edge from smaller to larger except
the edge `n-1 -> 0`. Call this tournament E_n, with n>=3. Then

```text
H(E_n)=1+2^(n-2).                                      (S1)
```

**PROVED by a bijection.** A Hamiltonian path avoiding the sole backward
edge must be the unique increasing listing. A path using that edge has
two increasing pieces: an arbitrary subset S of the interior vertices,
then n-1, then 0, then the complementary interior subset. There are
exactly `2^(n-2)` choices. This proof needs no external path-count theorem.

At `n=2^r+2`, (S1) is precisely `F_r=2^(2^r)+1`. Thus the tournaments on
3,4,6,10,18,34 vertices have 3,5,17,257,65537,4294967297 Hamiltonian
paths. The last count is composite. Nothing in the construction or
bijection changes there: the prime-to-composite transition is arithmetic
of the count, not a loss of the graph property.

The proposed ten-vertex array has 9 fixed path edges and 36 independent
chord directions. Relative to that path, these 36 chords index an
undirected fundamental-cycle basis; its cycle lengths are 3 through10,
with multiplicities 8,7,...,1. They are not 36 directed triangles.

Reverse every nonconsecutive chord instead, obtaining I_n. Its directed
simple cycles are exactly contiguous intervals, while E_n's cycles choose
interior subsets. Both have the same outdegree multiset and n-2 triangles.
Using the **CITED** odd-cycle collection identity of
[Grinberg--Stanley, Theorem1.39 and equation42](https://arxiv.org/pdf/2307.05569),
the interval family has a tribonacci Hamiltonian-path law. At order10,
I_n and E_n have193 and257 paths despite identical scores and eight
triangles. Cycle overlaps are the missing coordinate.

## 3. The repeated-three and Fermat families share a modular mechanism

Put `R_k=(10^(k+1)-7)/3`, where k counts the threes and R_0=1. The exact
first composite is

```text
R_8=333333331=17*19607843.
```

Both factors are prime and both are new to this sequence. More strongly:

**PROVED.** Every fifteen consecutive R_k are pairwise coprime, uniformly
over all starting indices. Fifteen is sharp: R_1 and R_16 share31.

The mechanism is the general append-digit identity
`gcd(A_m,A_(m+t))=gcd(A_m,(b^t-1)/(b-1))`. In the decimal family a common
prime at gap t must divide both `10^t-1` and `7^t-1`. For t=1,...,14 the
only possible primes are3,11,13,37. Three never divides R, and seven lies
outside the subgroup generated by ten modulo each other candidate. The
finite gap certificate therefore proves an all-index theorem. At gap15,
`ord_31(10)=15` and `10^2=7 mod31` give sharpness.

**PROVED conditional on the displayed number being prime.** For every
Fermat prime `p=2^(2^r)+1` with r>=2,

```text
ord_p(10)=p-1.                                         (S2)
```

Indeed p=17mod40, and Gauss pairing shows ten is a quadratic nonresidue.
Since p-1 is a power of two, its order must be the full group order.
Consequently every such Fermat prime divides R_k in exactly one index
class modulo p-1. For17,257,65537 the first indices are8,226,29252.
The primes3 and5 never divide R and require separate treatment.

This is an actual connection between the two families: for each prime
output of the tournament path-count construction other than3 and5,
ten generates the unit group modulo that prime and supplies a divisibility clock in R. The
first-hit phase remains necessary; group order alone does not give it.

**FINITE-EXACT:** every R_k through k=200 introduces a new prime,
checked by two exact gcd-stripping methods. This does not prove the
all-index statement. Nor does it put R_8 in the same exception class as
M_6=63: the latter introduces no new prime, whereas R_8 and F_5 do.

## 4. The rational three-cycle has a concrete twelve-element symmetry

For `x^2-29/16`, the cycle `-7/4 -> 5/4 -> -1/4 -> -7/4` is inherited.
Use the projective coordinate `x=X/(4Y)` and matrices

```text
B=(1/4)*[[1,-13],[1,3]],       R=[[-1,-2],[0,1]].
B^3=-I, B^6=I, R^2=I, RBR=B^-1.                       (S3)
```

**PROVED.** These matrices preserve the six integral points on
`X^2+2XY+13Y^2=48` and generate the dihedral group of order12, isomorphic
to `C_2 x S_3 = S_2 x S_3`. In coordinates U=X+Y,V=2Y, the norm is
`U^2+3V^2`; B acts as multiplication by `(1+i*sqrt(3))/2` and R as
negative complex conjugation. This is an exact hexagon with a specified
norm, group action, and projection.

The central sign is forgotten by projectivization, so its six vectors
give only three scalar points. Reflections reverse the cycle; only the
cyclic subgroup preserves its forward arrows. The same abstract lift
exists for every three distinct rational projective points. What is
special about c=-29/16 is the fixed affine reflection `x -> -x-1/2`:
invariance of a centered quadratic three-cycle under it forces this c.
The parameter's rational scalar map has no six-cycle, by the elementary
denominator and escape argument in the dynamics note.

There is a second exact construction:

```text
pi(z)=z+z^-1,          pi(z^2)=pi(z)^2-2.               (S4)
```

For z of odd order N>1, the trace period is the least m with
`2^m=+1 or -1 modN`. A minus return folds a squaring cycle of length2m
to a trace cycle of lengthm; a plus return pairs two length-m cycles.
This proves the complete period-six conductor list `13,21,63,65` for
`x^2-2`:54 real points in nine cycles, with multipliers ±64. These
algebraic cycles are not rational cycles at c=-29/16.

Here63 is structurally useful: `ord_3(2)=2`, `ord_9(2)=6`, and
`ord_7(2)=3`. Existing prime factors can support a larger period through
valuation depth and least common multiples. Conductor9 folds to trace
period3; conductor63 retains trace period6. Prime support alone loses
exactly the information needed to distinguish them.

**UNRESOLVED source identification:** the exact Anthropic `S_2 x S_3`
reference was not located. The dynamics note distinguishes permutation
groups, the manifold S^2 x S^3, and the nearby recent S^6 manuscripts.
No claim from those manuscripts is used as a proved dependency here.

## 5. The missing root repairs the arithmetic product

Let `T(t)=t(t+1)/2`, and `N=(A+1)(B+1)`, with A,B>=1. Then

```text
T(N-1)=T(A)+T(B)+T(AB)+AB(A+B+1),
T(N-2)=T(A)+T(B)+T(AB)+A(B^2-1)+B(A^2-1).              (S5)
```

The user's second polynomial is exactly the root-deleted count: the
difference is N-1, the star at the missing vertex. Decompose the pointed
Cartesian product into its root, A-axis, B-axis, and AB-interior to see
the mechanism. Edge counts alone do not specify a tournament orientation.

One coherent model for the proposed two-sided arithmetic is

```text
a boxplus b=a+b+1,       a boxtimes b=ab+a+b,
h(a)=a+1.
```

The shift h transports these operations to ordinary addition and
multiplication. On {-1,0,1,...}, they form a semiring with additive zero
-1 and multiplicative identity0. On natural labels the operations remain
closed but there is no additive identity. Their diagonals are2a+1 and
a^2+2a; their iterated power is `(a+1)^r-1`. This gives a precise model
of the requested arithmetic loop, without claiming to uniquely interpret
the undefined f and g notation.

**VERIFIED Lean scope:** the standalone package proves20 natural-number
statements about shift conjugacy, associativity, commutativity,
distributivity, diagonals, powers, and absence of an additive identity.
It has a clean public-root build, no custom axioms or proof placeholders,
and audited dependencies at most `propext` and `Quot.sound`. The graph,
prime, and dynamical theorems above are written proofs with exact controls,
not claims of Lean formalization.

## 6. A plane-tiling obstruction and its constructive repair

The strict summand graph misses exactly doubling arcs inside all ascending
pairs. The strict multiplicand graph misses exactly squaring arcs inside
proper divisibility pairs. Their ambient universes differ. Both underlying
graphs contain K_m for every m, so ordinary unrestricted minor membership
cannot distinguish them; neither full graph is planar.

As lattice incidence sets the multiplicative graph has only O(N logN)
points in an N by N box. **PROVED:** finitely many fixed invertible affine
copies of it, its reversal, and its removed square forest cannot cover
the lattice plane. This is not an obstruction to infinitely many
row-dependent copies or arbitrary nonlinear coordinate changes.

The constructive repair is to retain parents or remainders:
`(a,b)->(a,a+b)` is an integer shear, whereas `(a,b)->(a,ab)` reaches
only divisibility positions. For each fixed a>=1, restoring `z=ab+r`,
with b>=0 and `0<=r<a`, represents every nonnegative z uniquely.
Thus the user's micro remainder and macro
quotient are meaningful coordinates. A literal integer equation q=mq
has only q=0 for m!=1; modulo M it instead has `gcd(m-1,M)` solutions.
Naming the modulus changes the actual mathematical question.

## 7. Connection contracts and what remains useful for Collatz

| Source -> target and map | Preserved predicate | Lost information / needed sidecar | Decisive control |
|---|---|---|---|
| Pointed product -> shifted arithmetic via cardinality minus1 | Product identity | Distinguished root; orientation is extra | Root star has N-1 edges |
| E_n -> integer count H | Exact Hamiltonian-path count | Nearly all graph structure; no prime predicate is transported | F_5 is composite with unchanged construction |
| Fermat prime -> decimal residue orbit via powers of10 | Divisibility and period | First-hit offset | First hits8,226,29252 |
| Vector hexagon -> rational cycle via X/(4Y) | Cyclic successor | Central sign; reflection reverses time | B^3=-I but projective B^3=1 |
| Roots of unity -> Chebyshev cycle via z+z^-1 | Semiconjugacy | Inversion and return sign | N=9 folds, N=7 does not |
| Multiplication parents -> arc coordinates | Exact product | Remainder, divisibility ambient set | Pair(2,3) is absent |
| Finite residue description -> Collatz trajectory | Guarded finite word, if carry retained | Unbounded quotient and ordinary height | Arbitrarily long all-growth cylinders |

The [previous Collatz guard synthesis](collatz_guards_20260921_synthesis.md)
already proves why source density or a bounded residue automaton cannot
simply be read as descent on each orbit. The present work contributes
concrete carriers for some other analogies; it does not remove that
global quantifier. The strongest reusable move is to lift a proposed
finite picture until its action is exact, identify what the projection
forgets, and retain that coordinate in any proposed descent argument.

There is also a direct Collatz connection to the append-digit law. For
any positive odd x_0, set `x_j=4^j*x_0+(4^j-1)/3`. Then
`3x_j+1=4^j*(3x_0+1)`, so all x_j have the same accelerated Collatz image,
and `gcd(x_i,x_(i+t))=gcd(x_i,(4^t-1)/3)` for i>=0,t>=1. In particular
consecutive siblings are coprime. This recovers the inherited inverse-fibre
mechanism and adds its explicit repunit-return comparison; see the
[operations proof, section8](arithmetic_seams_20260921_operations.md).
The family indexes alternative predecessors of one node, not successive
forward iterates. Retaining that distinction prevents a false descent
or orbit-prime-novelty conclusion.

Three precise remaining directions replace an unspecified universal
exception principle:

1. **OPEN here:** prove or refute that every R_k, k>=1, has a primitive
   prime. A general append-digit theorem cannot be assumed: digit-one
   repunits already behave differently. The order/phase test supplies
   an exact starting criterion; the present evidence stops at200.
2. **OPEN research object:** classify which restrictions on chord
   interactions make an exact path-count recurrence depend on finitely
   many boundary states. The two solved families provide positive
   controls, and equal scores/triangle counts provide hostile controls.
3. **OPEN Collatz bridge:** construct a labelled directed carrier whose
   projection both preserves legal Collatz steps and retains a
   well-founded descent quantity. Neither ordinary graph minors, a
   shared group order, nor a prime-count pattern supplies such a map.

## 8. Reproduction and independent review

Run from the repository root:

```text
python 04-computation/experiments/arithmetic_seams_20260921_verify.py
```

The [manifest](arithmetic_seams_20260921_manifest.json) pins the inputs and
records normal/optimized exact replays plus the Lean clean build and
axiom audit. Finite universes are explicit in each lane:33,866 fixed-path
tournament configurations;18,720 general append-digit gcd controls;
125 rational marked-cycle parameters; and the operation-control boxes.
None is an enumeration of all ten-vertex tournaments or all positive
Collatz trajectories. Checks raise errors explicitly and remain active
under optimized Python.

The prime theorem, tournament bijections, quadratic lifts and operation
proofs received cross-lane read-only audits. They prompted an endpoint
repair: the squaring forest on positive vertices includes isolated1 as
well as the nonsquare roots at least2. The maintained mistakes ledger
records the mechanisms of the rejected implications, not a claim that
all of the user's intended interpretations have been uniquely determined.
