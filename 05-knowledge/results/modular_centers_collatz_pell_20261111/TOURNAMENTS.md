# From modular tables to genuine tournaments and extension obstructions

**11 November 2026. PROVED P9–P12 below; FINITE-EXACT controls in
verification.json. No literature-priority claim; Chan-Pak remains OPEN.**
The new integration strengthens the prior selected-set symmetry filter
to a coefficient-wise obstruction. No new global theorem ID is reserved.

## P9. A complete description of the Jacobi comparison graph

For odd m=3 mod4 take vertices Z/m and the intrinsic observable
chi_m(b-a), the Jacobi symbol. Orient a->b when it is +1; the value -1
gives the reverse arc, and zero is a missing edge. The sign gauge is
the choice of +1; switching it reverses every arc. Since
chi_m(-1)=-1, each unit difference has exactly one orientation.

Every vertex has outdegree and indegree phi(m)/2 and has
m-1-phi(m) nonneighbors. Indeed negation pairs the units into opposite
signs. Translations preserve all edges; multiplication by a unit u
preserves orientation for chi_m(u)=1 and reverses it for chi_m(u)=-1.
These claims retain the actual arithmetic vertex labels.

The graph is a tournament iff m is prime. For composite m the maximal
size of a subset inducing a tournament is exactly its least prime factor
p: such a subset must be injective modulo p, while {0,...,p-1} attains
the bound because every nonzero difference has absolute value <p.

CRT gives the full missing-edge rule: two vertices are adjacent iff
their residues differ at every prime divisor of m; the Jacobi sign is
the product of the local Legendre signs raised to the corresponding
prime exponents. Passing only to the radical loses multiplicities of
vertex fibres; retaining m/rad(m) recovers the blow-up size. This is an
exact classification of this comparison family, not of all tournaments.

For m=1 mod4 the same observable is symmetric, so it does not orient a
tournament. In particular 5,13 and33 have a different role from 3,7,11,
19 and1039. At25 all unit differences have Jacobi value +1. For the
actual prime-field tournaments at7 and1039 the special doubling
symmetry is inherited from
[PRIMES_AND_PAPERS.md](../odd_c_seed_symmetry_20261110/PRIMES_AND_PAPERS.md):
the orders of 2 are3 and519. A number used as a modulus or vertex count
has not thereby been identified with a seed slope of the same value.

At27 the Jacobi partial tournament is a directed triangle of three
independent nine-vertex fibres. It has243 arcs and108 missing unordered
pairs. This is not the Paley tournament on the field F27.

## P10. Resolving the missing levels produces a recursive family

Fix a prime p=3 mod4. For each e>=1 and nonzero difference d modulo p^e,
write d=p^v u with u a unit; its reduction modulo p is independent of
the chosen representative. Define a->b iff the leading unit of b-a is
a quadratic residue modulo p. Negation changes its sign, so this is a
genuine tournament R_(p,e) on all p^e vertices.

Writing vertices in base p, the first unequal digit, read from least
significant upward, determines the orientation. Therefore

    R_(p,e) = P_p[R_(p,e-1)],  R_(p,0)=K1,             (1)

where P_p is the prime-field Paley tournament and brackets mean
lexicographic substitution. This proves the recursion, rather than
assuming that integer-ring multiplication is finite-field arithmetic.
When e is odd this completes the Jacobi orientation; when e is even
the original Jacobi observable was symmetric, and (1) is a different
specified observable. The valuation and unit sign are the retained data.

Every vertex has (p^e-1)/2 outgoing arcs. The complete valuation profile
of directed arcs is

    A_(p,e,v) = p^e (p-1)p^(e-v-1)/2,  0<=v<e.        (2)

For p=3,e=3 the levels have243,81,27 arcs. Formula (2) totals
27*26/2=351. This supplies an exact triangular sequence at every depth.
Local translations inside each of the p^(e-1) bottom p-vertex blocks
are automorphisms supported on that block.

The ring-derived R_(3,3) and field P27 are not isomorphic: direct
computation gives common-outneighbor counts {4,10,12} in the former
and the single value6 in the latter. Their vertex counts alone do not
identify their multiplication, decomposition or symmetry structure.

The substitution operation and its need for incidence data were already
used in [Paley palettes and carry](../paley_multicoloring_and_carry_20261021/session23b-palettes/PALEY_MULTICOLORING_AND_CARRY.md).
Here the new role is to expose small automorphisms supported away from
the few endpoints touched by an extension.

## P11. Local endpoint symmetries force global divisibility of extensions

Let T be a tournament and k>=0. Suppose there is an integer d>1 such
that **for every** S subset V(T) of size at most2k, Aut(T) contains a
subgroup G_S of order d fixing S pointwise. The subgroup is allowed to
depend on S. Then

    d divides H(U) for EVERY tournament U extending T by k vertices. (3)

For k=0, this is the free action of an order-d automorphism group on
Hamiltonian paths. For k>0 use the inherited exact extension compiler
[EXTENSIONS.md](../odd_c_seed_symmetry_20261110/EXTENSIONS.md): a qualifying
permutation of old and named new vertices prescribes a partial assignment
of the orientation bits involving new vertices. Its monomial is the
indicator of that assignment. Group all qualifying permutations with
the same exact assignment c, and let N_c be the size of that group.

At most2k old vertices are incident to a new vertex along such a path.
These names are fixed by the assignment c; call their set S_c. The group
G_(S_c) preserves all old arcs and fixes every old endpoint appearing
in c, so it acts on this set of qualifying permutations. The action is
free: an automorphism fixing a permutation fixes every old vertex at
its specified position. Consequently d divides N_c for every c.

Each actual extension evaluates a sum of these coefficients times zero
or one. This proves (3). Equivalently every coefficient of the grouped
orientation-cylinder polynomial is divisible by d, and expanding the
indicators also gives an ordinary integer polynomial divisible by d.

This argument does **not** require one subgroup to preserve all of the
new attachments, or any nontrivial automorphism of U itself. It retains
the names of the old endpoints before taking the sum. This is the
additional coordinate that the single selected-set stabilizer misses.

### Consequence A: all nonprime finite-field Paley bases are excluded

Let q=p^f=3 mod4, f>1, and T=P_q over F_q. For every two distinct field
points a,b, the permutation

    g_(a,b)(x)=a+(b-a)*((x-a)/(b-a))^p                 (4)

is conjugate to Frobenius and has order f. It fixes a,b. It preserves
Paley orientation: its affine coefficient (b-a)^(1-p) is a square,
and Frobenius preserves squares. This is the inherited pair-fixer
construction, now applied **before summing path-cover coefficients**.
Single points and the empty set are covered by choosing additional
distinct points. Thus P11 applies with k=1,d=f:

    f divides H(U) for every one-vertex extension of P_q.            (5)

Also f divides H(P_q). Since q=3 mod4 forces f odd and f>=3, neither
3H(P_q)+2 nor3H(P_q)+4 is divisible by f. Hence these Paley bases cannot
realize either one-vertex Chan-Pak induction branch. This excludes a
specific infinite source family; it does not exclude other sources of
the same H or extensions with two or more new vertices.

For the seed ray, use the established exact formula

    H(T,O,t)=H(T)-W_T(O)+2^t U_T(O), t>=1,
    W_T(O)=sum_{a<b in O} W_ab, U_T(O)=sum_{a in O} U_a.

Here W_ab counts ordered two-path covers whose joining endpoints are
a,b in either order, and U_a counts the ordered covers whose second
path starts at a, allowing the first path to be empty. The same group
(4) acts freely on W_ab objects, and an order-f point fixer acts freely
on U_a objects. Hence

    f | H(T), f | every W_ab, f | every U_a.

It follows that **no selected set O whatsoever** in a nonprime Paley
tournament gives offset +1 or-1. Every count on the specified seed
ray is divisible by f. This strictly strengthens the earlier exclusion
of selected pairs and their complements. It does not say that every
actual extension is symmetric or that symmetry breaking is impossible.

### Consequence B: recursive prime blocks survive several added vertices

For R_(p,e) with e>=2 there are B=p^(e-1) bottom blocks. Any set of at
most2k vertices with2k<B misses one whole block. Translation on that
block supplies an order-p subgroup fixing the selected vertices.
Therefore

    p | H(U) for every k-vertex extension whenever 2k<p^(e-1).      (6)

For R_(3,3) on27 vertices this covers k=0,1,2,3,4. The conclusion
includes arbitrary orientations both among the new vertices and to the
old ones. For any two old endpoints there is an untouched bottom block,
so every W_ab and U_a is divisible by p. Thus no R_(p,e), e>=2, supports
a unit-offset seed at any transitive-tail length t either.

### An asymmetric hostile to a global-automorphism explanation

For R_(3,2) on9 vertices, exact counts are

    H=3159, every U_a=5085, W_ab in {1164,1242}.

All512 one-vertex attachments have H divisible by3. Attach new vertex9
with outgoing neighbors {0,1,2,3}; every other old vertex points to9.
The resulting10-vertex tournament has H=15033 and a **trivial** full
automorphism group. Divisibility survives even though no nontrivial
symmetry survives in the child. Endpoint-dependent groups in the
coefficient proof explain this example.

The test compares subset-DP counts with literal enumeration of all
9! permutations and every permissible two-run cut, exhausts all512
attachment masks, and directly recounts nine actual ten-vertex children.
Backtracking checks the full automorphism group of the displayed child.

## P12. An exact infinite sequence of symmetry-breaking subsets

For R_(3,e), retain the explicit iterated wreath subgroup W_e generated
by independent rotations of the three child blocks at every node of
its rooted ternary decomposition tree. This is a subgroup of its
automorphism group; no equality with the full group is needed below.
Put n=3^e. Then

    |W_e| = 3^((n-1)/2).

Exactly 2|W_e| subsets of the n vertices have trivial stabilizer under
W_e. They form two orbits, one of cardinality (n-1)/2 per subset, and
the other of cardinality (n+1)/2. In particular any unit-offset seed
would have to use one of these two selected cardinalities; P11 then
excludes all of them once e>=2.

**Proof by recursion.** At e=0 a one-point set has two subset types,
empty and full, each with trivial stabilizer. A subset at the next level
must first have trivial stabilizer in each child block. Each child
therefore has one of the two types. If all three types agree, a cyclic
permutation of the blocks, composed with appropriate child-group maps,
preserves the subset, producing a nontrivial stabilizer. If they are
not all equal, no nonzero cyclic rotation preserves their type word,
and trivial child stabilizers then force the full stabilizer to be
trivial. The six nonconstant binary triples form two rotation orbits.

If child sizes are (m-1)/2 and(m+1)/2, these two patterns have total
sizes (3m-1)/2 and(3m+1)/2. Also |W_e|=3|W_(e-1)|^3, proving the order
and the two free-orbit counts by induction.

The new tracked size/depth sequences are consequently explicit:

| e | n=3^e | number of W_e-free subsets | possible sizes |
|---:|---:|---:|---|
| 0 | 1 | 2 | 0,1 |
| 1 | 3 | 6 | 1,2 |
| 2 | 9 | 162 | 4,5 |
| 3 | 27 | 3188646 | 13,14 |

At n=9 the two free subset orbits have offsets -3903 and-8637, so
rigidity gives no unit offset. The counts and offsets through e=2 are
exhaustive; the general sequence follows from the proof. This makes an
infinite symmetry family manageable while explicitly retaining the
stronger arithmetic obstruction.

## What this changes for enumeration

For actual extension searches, group qualifying permutations by their
orientation demands and test the content of their coefficients before
enumerating all Boolean assignments. P11 supplies an all-size certificate
of a common divisor for the families above. An excluded target needs
no attachment enumeration. A divisor of a mean alone would not suffice.

The valuation profiles, local-fixer depths and wreath-free subset counts
are exact aspect sequences for specified infinite tournament families.
They are not formulas for A000568 or for all metagraph components.
The thin Pell prime branch in PROOFS.md supplies candidates and exact
A135244 predicates, not a faster general quadratic-prime sieve.
