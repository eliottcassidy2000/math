# Recovering hidden orientations, growing prime sieves, and the Paley branch cover

2026-10-26. **PROVED elementary deductions, with the classical
Busch separator theorem CITED; FINITE-EXACT computational controls.**
No historical-priority claim. The synthesis distinguishes faithful encoding,
exact finite evaluation, and unresolved existence at unbounded height.

The main advance is an all-order reconstruction theorem. A tournament is
determined by the directed arcs used in its Hamilton paths, even though
some arcs never occur in any path. The missing orientations are forced by
detectable one-vertex cuts. This closes the named legal-guard injectivity
question left open by the recent strong-atom packet.

The arithmetic continuation has a different repair: growing prime tests,
with their prime-equality exceptions and a certified height bound, give an
exact recursive twin count. The F8 continuation supplies a third repair:
one determinant bit turns a conic point into an actual oriented seven-point
tournament. These are precise, different ways to recover a chosen target.

## 1. Hamilton-path support determines the tournament

For a labelled tournament T, let L(T) be its full directed Hamilton-path
language, and let U(T) be the spanning digraph containing the arcs that
occur in at least one such path. Then

    T -> L(T) -> U(T)

is lossless on its image: there is an explicit equivariant inverse from
U(T) to T. The syntactic loss of individual path words in taking their
union does not lose the underlying tournament.

The reconstruction has three steps:

1. Recover the ordered strongly connected components from U(T).
2. Orient every missing pair in different components forward in that order.
3. For a missing pair x,y inside one strong component C, find a vertex z
   such that U(T)[C]−z has a path in exactly one direction between x and y.
   Orient the pair in that direction.

Every required z exists and every qualifying z gives the same orientation.
The key new lemma is

    reach(U(T)−z) = reach(T−z),  when T is strong.

Its proof retains Hamilton-cycle edges after deletion by using two families
of paths: prepend z to paths starting in the first component, and append z
to paths ending in the last. If deletion leaves a strong tournament, breaking
different edges of one Hamilton cycle supplies the corresponding proof.
The strong hypothesis is essential: the transitive triple loses its long
arc from U, and deleting the middle vertex destroys its only support path.

The imported structural fact is Busch's 2005 thesis Theorem4.12: an unused
arc in a strong tournament has its endpoints separated in the ordered
strong components after a suitable vertex deletion. This necessary condition
is not treated as an iff test. The primary statement, proof pages, attribution
to Busch/Jacobson/Reid, and full reconstruction proof are in
[P1–P4](paths/RESULT.md), with the [independent proof audit](audit_paths/PROOF_AUDIT.md).

The first strong example with an unused arc has five vertices. Its adjacency
rows are `(20,1,2,7,14)`; U changes the last row to12, deleting `4->1`.
After deleting vertex0, its support has the order

    4 -> 3 -> 2 -> 1,

with additional forward arcs. This cut forces the missing orientation4->1.
Ordinary reachability in U alone cannot do this, because U is strong.

Equivariance has further exact consequences:

    Aut(U(T)) = Aut(T),       L(U(T)) = L(T),
    H(U(T)) = H(T),           H(U(T))/|Aut(U(T))| = tc(T).

The first identity follows by applying the inverse to a relabelling; the
second follows because U contains exactly the arcs used by T's paths and
is a subgraph of T. Thus isomorphism classes of tournaments correspond
bijectively to the image of U, preserving both raw and normalized path counts.
This is a faithful classification carrier, not an enumeration of that image
or an unrestricted A000568 recurrence.

In the fixed-path presentation cube X_n, all `2^binom(n−1,2)` named
one-step legal-guard vectors are therefore distinct at every order. The
coarsest quotient retaining these guards has exactly that many states.
This closes the all-order named-guard question in
[strong-atom C5](../strong_atom_continuation_vitali_20261026/README.md).
The fresh, unnamed successor-class observations of HYP-9416 are different
and their stated all-order questions remain open.

## 2. A growing sieve with a precise stopping horizon

Keep the previous twin-ray coordinate

    T(x)=(x(x+2),2(x+1),(x+1)^2+1),  odd x>=3.

For B>=2 let M_B be the product of primes at most B and define

    R_B(H)=#{3<=x<=H : gcd(x(x+2),M_B)=1}.

A strict wheel rejects an entry equal to one of its own sieving primes.
The corrected local test allows `ell|x` only when x=ell, and likewise for
x+2. Its accepted set is exactly the disjoint union of strict survivors
and actual twin pairs whose lower prime is at most B. This is an equality
proved from small proper factors, not an arbitrary exception list.

Let q be the next prime after B. If `H+2<q²`, every surviving pair in
the height-H box is prime. Writing pi_2(H) for twin lower primes at most H,

    pi_2(H) = R_B(H) + pi_2(min(H,B)).

In particular, for H>=3,

    pi_2(H) = R_floor(sqrt(H+2))(H)
                  + pi_2(floor(sqrt(H+2))).

The recursive argument decreases. The strict count has exact residue-floor
and Möbius inclusion-exclusion formulas, so this is a finite counting
certificate. It specializes classical sieve mechanics; no improvement of
analytic sieve bounds is claimed. [W1–W3 and the exact formulas](sieve/README.md).

The square boundary is sharp as a universal rule, but need not be attained
at every prime stage. At x=q²−2 the composite upper entry q² survives
exactly when the lower entry is prime. A stronger square-only pattern fails:
at cutoff61 the first false pair is `(4757,4759)=(67*71,4759)`, with4759
prime. The earlier possible boundary4487 is already divisible by7.

The limit statements must stay distinct:

- Every finite strict wheel has infinitely many positive survivors.
- No positive x>=3 survives all strict wheels, because a prime factor of x
  eventually enters the wheel, even when x itself is prime.
- The exception-restored tests stabilize pointwise to the twin predicate,
  uniformly on each bounded-height interval at the stated square threshold.

An explicit compatible wheel path is x congruent to−1 at every stage.
It represents the integer−1 and the degenerate signed triple `(−1,0,1)`,
outside the positive primitive-triangle domain. Its positive representatives
M_B−1 grow without a fixed height bound. Thus compatible residue paths do
not establish infinitely many positive twin pairs.

## 3. The F8 bit carries actual tournaments

The inherited paired F8 designs gave a two-sheet observation over
`P1(F7)`. Here is its literal arc formula. Choose projective representatives
v_t=(t,1), v_infinity=(1,0), and set

    K(p,u,v)=det(v_p,v_u) det(v_u,v_v) det(v_v,v_p).

At state (p,e), orient the seven points other than p by

    u -> v  iff  Legendre_7(K(p,u,v))=(-1)^e.

Representative rescaling changes K by a square. A projectivity g changes
its square class by det(g), so the exact transport is

    (p,e) -> (gp,e + determinant_square_bit(g)).

The sixteen objects form the homogeneous cover

    PGL2(7)/H21 -> PGL2(7)/B42,

where H21 is the square-determinant pole stabilizer and B42 the full one.
The integer root requires the explicit conjugation described in the
[field proof](field/README.md). U,D preserve the sheet; A flips it.
Each pole has two converse Paley tournaments, each with

    H=189,  |Aut|=21,  tc=9.

This connects branch parity to actual orientations and the order-21 arc
action. Moreover U(T)=T for each of these tournaments: the automorphism
group is transitive on arcs, and one arc of one Hamilton path can therefore
be carried to any other arc. The sixteen-state family sits inside the
reconstruction image above without any missing arcs.

The inherited raw-H exclusions7,21 transport to the whole image of U
because H is preserved. Their proof remains
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The separate through-eleven tc exclusions remain
[THM-4790 / presentation-fibres-7-and-21-absent-through-eleven-vertices](../../../01-canon/theorems/THM-4790-presentation-fibres-7-and-21-absent-through-eleven-vertices.md).
The Paley stabilizer's order21 does not prove either gap theorem.

## 4. Arithmetic predicates do not all descend through this cover

The branch observer can be refined by the Euclid vector modulo a prime
wheel, with reductions commuting with the branch matrices. At modulus15,
all3072 joint states are reached. This is the ambient transition system;
the subset passing a prime filter need not be closed under a branch move.

Even this refinement plus the true twin predicate does not decide an abc
hit in the **squared Pythagorean** packet `a²+b²=c²`:

| Euclid pair | Twin radii | Same finite state | Radical R | c² |
|---|---|---|---:|---:|
|(12,1)|(11,13)|t5, bit0, vector(12,1) modulo15|124410|21025|
|(432,1)|(431,433)|t5, bit0, vector(12,1) modulo15|8358844170|34828890625|

The first has c²<R; the second has c²>R. This does not contradict the
previous result for the different **centre** identity `1+x(x+2)=(x+1)²`,
whose twin-selected triples always have radical larger than height.
There are also two explicit branch words with equal full radical but
opposite determinant bits. The exact factors and independent checks are
in the field proof and [audit](audit_field/AUDIT.md).

The F1 language here consists of explicit pointed monoids and their
integral monoid algebras. On the free branch-word monoid, determinant is
evaluation of the commutative word monomial at `(1,−1,1)`; augmentation
evaluates at `(1,1,1)` and loses this character. Commutativization itself
retains the character but loses branch order and conic endpoint.
On prime-divisor polynomials, prime-weight evaluation gives aliquot sums,
whereas augmentation counts proper divisors. These are different maps
with explicit non-descent witnesses, not interchangeable F1 realizations.

## 5. What unifies the successful moves

The [mechanism note](METHOD.md) turns the common pattern into a precise
test: a target descends through an observation exactly when it is constant
on every observation fibre. To continue an operation, retain its guards
and the future distinctions that operation can create. The new results
identify three different repairs:

| Apparent missing information | Actual repair | Remaining obligation |
|---|---|---|
|Unused tournament orientations|Infer them from structural cuts; no extra state|Recognize or enumerate the image for a proposed taxonomy|
|Which oriented Paley object a branch carries|One determinant bit over the conic point|Arithmetic targets need their own retained data|
|Primality after finitely many congruence tests|Grow the test horizon and retain prime-equality exceptions|Positive survivors at unbounded height, not just finite wheels|

This explains a repeatable mechanism behind the cited advances without
declaring the different open problems equivalent. Exact maps organize
their obligations; they do not discharge image realization, termination,
uniformity or positivity automatically.

## Reproduction and concurrent integration

The [concept board](CONNECTION_BOARD.md) names inheritance, hostiles and
operation boundaries. The [audit record](AUDIT.md) separates proof review,
finite controls, source pins and status-only promotions. The incoming S60
audit merge restored several unqualified headlines despite retaining the
conditional proof; [exact repairs](REGRESSION_CORRECTIONS.json) restore the
scope and remove duplicate superseded rows while preserving its audit provenance.

    python3 -I -B 05-knowledge/results/growing_sieve_path_20261026/verify_packet.py

Twin-prime infinitude, abc, fresh re-marking minimality, unrestricted
tournament-spectrum completeness, all-order tc exclusions and black-core
connectivity remain outside the proved conclusions here.
