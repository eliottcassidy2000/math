# Degree-57 correction: the involution branch is closed

2026-11-24. This supersedes the proposed involution search in
[the count-address write-up](README.md), not the exact tournament
counts or the physical +2 construction. [Ishida's June 2026 preprint,
Theorem 1.1](https://arxiv.org/html/2606.29183v1) already rules out
involutions, using the previously known fixed-star theorem. Below is a
short derivation of the trace test, a direct rejection of **both**
provisional branches in the prior write-up, and its odd-order use.

## Elementary projector trace test

For a prime p, take R=Z_(p) and a vertex permutation T of order p.
Let E be an idempotent over R commuting with T, and L=E R^V.
Set D=1−T, U=1+T+...+T^(p−1), L_+=ker(D|L), and L_0=ker(U|L).
On the permutation module ker U=im D: a free p-orbit has the
sum-zero differences and a fixed site contributes no kernel over R.
Applying E gives the same equality on L.

Let a=rank L_+ and (p−1)b=rank L_0. Over Q, the nontrivial
summands have polynomial Phi_p, so tr(T|L)=a−b and
det(D|L_0)=p^b. The maps D and U identify L/(L_++L_0)
with L_0/DL_0 and UL/pL_+, respectively. The latter thus has
order p^b. Since pL_+⊆UL, the quotient L_+/UL is an F_p-space
of dimension a−b. On the whole permutation module this
quotient is F_p on the fixed sites: the fixed vectors on free
orbits are norms. An idempotent splits this quotient and
its induced matrix is exactly E restricted to the fixed sites
and reduced mod p. Thus **tr(TE) equals the rank of that minor**.
This proof does not require classifying permutation lattices.

## Closure at 2 and the cost of a symmetric lift

Let A be adjacency, J the all-one matrix, f the number of fixed sites
and a1 the number of vertices adjacent to their image. The −8
projector is E=(455I+J−65A)/975, with 2-unit denominator, and
its trace with T is (7f+50−a1)/15. The previously derived
star branch has (f,a1)=(56,112), giving 22. Its principal
minor mod 2 is I+J+A_fixed: the center row and column vanish
and the 55-leaf block has rank 54. The provisional fixed
Hoffman–Singleton branch would give (50,700) and trace −20,
impossible for a rank. Neither branch survives, in agreement
with Ishida and the earlier fixed-star theorem he cites.

There is a precise encoding consequence. If a blue tiling decodes
a graph by merging reflection-invariant equivalence classes of
sites (including the earlier independent-twin recipe), reflection
induces a graph automorphism of order dividing 2. In degree 57 it
must be the identity. On an even spine the reflection is free,
so all 3250 decoded vertices must have even nonempty fibers:
at least **6500 sites**, not a 3300-ish fixed-site patch.
That bound is attained *conditionally* for the twin decoder:
duplicate every graph vertex, swap within each pair, and order
one representative from each pair along a Hamilton path of
the complement, then their twins in reverse. Dirac applies
since the complement minimum degree is 3192. The middle twin
edge is not in the graph mask, and reversing the graph edges
produces a blue tiling with exact recovery by merging twins.
This is not a bound for arbitrary nonlinear decoders, and does
not assert a graph exists.

## Four odd-order filters

The 7-projector E=(400I+50A−J)/750 is integral at the four
primes below and tr(TE)=(8f+a1−65)/15. Use the fixed graphs
from [Mačaj–Širáň, Lemma 4](https://www.sciencedirect.com/science/article/pii/S0024379509003735).
For a 7-star f≡2 mod 7 the center rows give zero and the
leaf block is I−J mod 7, rank f−2. For the pentagon mod 11,
the constant mode is killed and the nonconstant A-eigenvalues
are 7 and 3 of multiplicity two each; the latter is killed.
For a singleton mod 19, the numerator is 399=0 mod 19.

| Order | Fixed graph | Minor rank | a1 | Internally cyclic moved orbits |
|---:|---|---:|---:|---:|
| 7 | star of size f | f−2 | 7f+35 | 3(f+5) |
| 11 | pentagon | 2 | 55 | 25 |
| 13 | empty | 0 | 65 | 30 |
| 19 | singleton | 0 | 57 | 27 |

The rank applies to every nonidentity power. A free p-orbit
induces a circulant; two distinct unoriented nonzero step
sizes would give a square, so its internal graph is either
independent or a single p-cycle. Summing a1 over the p−1
powers counts 2p for every cyclic orbit. In particular for
order 11 all **25** cyclic orbits are in the 2970-site
residual graph: each of the five exterior pentagon columns
is independent (its vertices have a common fixed parent),
and the five fixed sites account for the rest. This gives a
separate count for the residual matching search.

## Exact Hamilton-path height filters

For a fixed pentagon the exterior return commutes with the
symmetry on each column. Every free p-cycle is a sheet orbit.
On one cycle of m quotient sheet orbits the composed offset is
in F_p. Zero offset lifts to p cycles of length m; otherwise
it lifts to one cycle of length pm. This also holds for the
rooted fixed-spoke triple. Apply the existing
[prime-packet height](../moore_mobius_prime_packets_20261124/)
with per-cycle factor 2^(k*m−2)+1 (k=5 for the pentagon,
k=3 for the triple) and deduplicate exact partitions and integer
products. These counts are candidate **integer heights** of
the relevant tournament encodings, not existence results:

| Symmetry and marked structure | Free orbit blocks | Fixed sheets | Distinct integer heights |
|---|---:|---:|---:|
| order 11, fixed pentagon | 5 | 0 | **36** |
| order 5, fixed graph a pentagon | 11 | 0 | 685 |
| order 5, fixed graph Hoffman–Singleton, mark a fixed pentagon | 10 | 5 | 444 |
| order 7, fixed-spoke triple | 8 | 0 | 74 |

The five fixed sheets of the third row have identity return
by the [Hoffman–Singleton exterior-cycle bound](../moore_mobius_prime_packets_20261124/RESIDUAL.md).
The last row excludes length-one cycles by triangle-freeness,
and requires at least three fixed spokes.
For order 11 the first three heights are
9^55, 9^44(2^53+1), and 9^33·257^11; the singleton count
c_1 is divisible by 11. The independently derived residual
pentagon count is therefore 36,383,149−c_1, with c_1 in
{0,11,22,33,44,55}. This is a direct odd-order successor
to the closed 55-matching involution branch.

Replay: run `node trace-replay.js` alongside this note.
