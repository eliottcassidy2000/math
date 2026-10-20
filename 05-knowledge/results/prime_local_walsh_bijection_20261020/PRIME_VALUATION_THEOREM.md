# A prime-local automorphism-preserving Walsh bijection

**Status: PROVED, conditional only on the proved elementary Walsh-basis mechanism reproduced below; FINITE-EXACT controls for every class through n=6.** This is a bijection on isomorphism classes. It preserves the order of a Sylow p-subgroup for one selected odd prime p, not its isomorphism type, conjugacy class, or action. It does not preserve all odd primes simultaneously and is not a vertex-relabeling-equivariant map.

## 1. Inheritance and scope

The input is THM-4721, `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`, Theorem W. Its full note is `05-knowledge/results/even_graphs_20261019_sylow_refinement_components_walsh_duality.md`. The existing explicit SCC–Lyndon codec (`05-knowledge/results/tournament_even_bijection_20261011c.md`) pairs irreducible atoms by equal rank. This note uses nonzero signed-subgraph coefficients instead and adds a prime-valuation invariant.

Here **even** always means Royle sign-even: every automorphism reverses an even number of edges relative to any reference orientation. It does not mean every degree is even. The classical equality of the two class counts is [Royle, Praeger, Glasby, Freedman, Devillers, *Tournaments and Even Graphs are Equinumerous*, Theorem 1.1](https://arxiv.org/pdf/2204.01947), published in *J. Algebraic Combin.* 57 (2023), 515–524, DOI [10.1007/s10801-022-01197-0](https://doi.org/10.1007/s10801-022-01197-0). Their paper asks for a natural bijection between unlabelled sets; the extra requirement of a labelled equivariant lift is a separate stronger notion. No claim of resolving their intended naturalness criterion is made here.

The new application is the full distribution of one odd-prime valuation, rather than only the product of automorphism orders or occurrence of an ambient Sylow subgroup. Its reusable mechanism is a weighted matrix pairing over Z[1/2].

## 2. The general weighted-pairing theorem

Let M be an N-by-N integer matrix invertible over R=Z[1/2]. Let a_1,...,a_N and b_1,...,b_N be positive integers and c a positive unit of R (thus a power of 2). Suppose

    M diag(b_j) M^t = c diag(a_i).                       (1)

For every odd prime p there is a permutation pi such that

    p does not divide M[i,pi(i)],
    v_p(a_i) = v_p(b_pi(i)) for every i.                 (2)

Consequently the multisets of individual p-valuations on the two sides are equal, and for every k>=0 the counts with p^k dividing the weight agree.

**Proof.** Invert (1) on the right to obtain

    (M^-1)[j,i] = b_j M[i,j] / (c a_i).                  (3)

All entries of M^-1 lie in R and hence have nonnegative p-valuation. Therefore

    v_p(M[i,j]) >= max(0,v_p(a_i)-v_p(b_j))              (4)

for nonzero entries. Since det(M) is a unit modulo p, the Leibniz determinant expansion has a nonzero term modulo p. Choose its permutation pi. Its entries are p-units, so (4) gives v_p(a_i)<=v_p(b_pi(i)) for each i. Taking determinants in (1) gives

    det(M)^2 product_j b_j = c^N product_i a_i.

The p-valuations of det(M) and c are zero. Thus the sums of the valuations of the a_i and b_j are equal. Summing the nonnegative differences above forces every difference to vanish, proving (2). □

This proves more than equality of products, but much less than equality of the joint prime-valued vectors. A nonzero entry modulo p need not itself have equal endpoint valuations. Equality is forced on a **complete matching**, not on every edge.

### The valuation filtration

Sort rows and columns separately by increasing p-valuation. The now-equal valuation multisets give aligned block sizes. Equation (4) makes M modulo p block upper triangular: an entry is zero if its row valuation is greater than its column valuation. Its determinant is nonzero, so every square diagonal block is invertible modulo p. One can find the bijection independently inside these diagonal blocks. Every perfect matching in the full p-nonzero support automatically remains inside them. Upper off-diagonal nonzero edges cannot belong to any complete matching.

Equivalently (1) is a congruence of diagonal forms over Z_(p); their p-primary Smith data agree. The direct determinant argument above establishes the needed conclusion without importing a classification theorem for forms.

## 3. Why the tournament Walsh matrix meets every hypothesis

Put m=binomial(n,2). Tournaments on [n] form the Boolean cube of choices of edge orientations. For an edge set X, let

    phi_X(T) = product of T(e) over e in X,

where T(e)=+1 for the increasing-label orientation and -1 for its reverse. The character matrix satisfies H H^t=2^m I, so its inverse 2^-m H^t lies in R. Thus the phi_X form an R-basis of all functions on labelled tournaments.

Vertex permutations act on this basis by signed permutation of the edge sets. On the line of phi_X the stabilizer Aut(X) acts by its edge-orientation sign. If X is odd, an invariant coefficient on its signed orbit must equal its own negative and hence be zero in R. If X is even, choose one representative X and transport coefficient 1 consistently across its orbit; this gives an invariant signed orbit sum Phi_X. These sums, one for each even graph class, form an R-basis of invariant functions. Orbit-indicator functions of tournament classes form a second R-basis. Hence their evaluation matrix

    M[T,X]=Phi_X(T)

is integral and invertible over R. This argument gives the equality of the two indexing cardinalities as well. Choice of the graph representative changes a whole column's sign, with no effect on support or the construction below.

Character orthogonality, grouped by vertex-isomorphism orbits, gives

    M diag(|Aut X|) M^t = 2^m diag(|Aut T|).             (5)

For completeness, the other orthogonality relation is

    M^t diag(|orb T|) M = 2^m diag(|orb X|),

and |orb Y|=n!/|Aut Y| converts one relation into the other. Thus (1) applies with a_T=|Aut T|, b_X=|Aut X| and c=2^m.

**Tournament consequence.** For every n>=0 and every odd prime p, there is an explicit bijection beta_p from tournament isomorphism classes to sign-even graph isomorphism classes such that

    p does not divide Phi_beta_p(T)(T),
    v_p(|Aut T|)=v_p(|Aut beta_p(T)|).

In particular, for all k>=0,

    #{[T]: p^k divides |Aut T|}
      = #{[X] sign-even: p^k divides |Aut X|}.           (6)

When k is the p-valuation of n!, (6) recovers the count of types whose automorphism group contains a Sylow p-subgroup of S_n. For smaller k it is an aggregate over subgroup shapes and actions; it does not establish HYP-9340(a) for a fixed non-Sylow subgroup H.

## 4. Fully specified construction, not equal-rank pairing

Encode each graph by its least edge mask under all vertex permutations; encode each tournament by its least reversal mask relative to i->j, i<j. Pairs are (0,1),(0,2),... in lexicographic order, low bit first. Order both finite class lists by these canonical masks. Form M modulo the chosen prime p, and its inverse Q.

At each stage use the first remaining row r. Choose the first remaining column j for which both M[r,j] and Q[j,r] are nonzero. Pair r to j and delete that row and column. The remaining minor is invertible: its determinant is det(M) times Q[j,r], up to sign. Existence of j follows from

    sum_j M[r,j] Q[j,r] = 1.

Repeat until empty. These are coefficients of the **original minor**, not a Schur complement used as a new matching support. This distinction prevents fill-in entries from becoming spurious permitted pairs. The inverse of the deleted minor is updated by

    Q_new[k,l] = Q[k,l] - Q[k,r] Q[j,l] / Q[j,r],

for k!=j, l!=r. The indices of Q are columns-by-rows. The standard block inverse identity verifies the formula, or multiplication by the retained minor does so directly.

Every chosen original coefficient is a p-unit; the theorem then guarantees equality of endpoint valuations. The finite table is a permutation and its inverse is read backwards. For an arbitrary labelled input, canonicalize its type, look up the image type, and return its canonical graph representative. Output labels do not covary with input labels.

This uses a canonical order to settle choices but does **not** pair equal list ranks. It is global and requires the full atlas. Given M, modular Gaussian elimination and the downdates use O(N^3) field operations and O(N^2) memory; these bounds do not include constructing the exponentially large class atlas or evaluating all coefficients. No polynomial-time-in-n claim is made. Column sign-gauge changes leave both nonvanishing tests unchanged, so the selected pairs are gauge-independent.

## 5. Sharp hostiles: one prime at a time

### Three vertices: equality is not edgewise

In row order (transitive, directed triangle) and column order (empty, P3),

    M = [[1,1],[1,-3]],  a=(1,3), b=(6,2), c=8.

Modulo 3 the upper-left entry is nonzero but joins valuations 0 and 1; it cannot appear in a full matching modulo 3. The unique such matching is transitive->P3 and triangle->empty. Over Q both matchings are supported, and the canonical cofactor recursion instead chooses transitive->empty and triangle->P3.

### A two-by-two general obstruction

The independent reviewer supplied the compact weighted-pairing example

    M = [[3,1],[-5,1]], a=(3,5), b=(1,15), c=8.

Indeed det(M)=8 and M diag(1,15) M^t=8 diag(3,5). The matching modulo 3 is off-diagonal and the one modulo 5 is diagonal. No single matching preserves both valuations. This isolates the arithmetic obstruction independently of graph enumeration.

### Five vertices: the actual first joint failure

The exact five-vertex census gives tournament automorphism orders

    1 (7 types), 3 (4 types), 5 (1 type),

and even-graph automorphism orders

    2 (5 types), 4 (1), 6 (1), 8 (2), 12 (1), 24 (1), 120 (1).

The empty graph has automorphism group S5, whose odd part is 15. No five-vertex tournament has automorphism order divisible by both 3 and 5. Thus no class bijection preserves both valuations simultaneously, even without the Walsh-support condition. This is already the first failure: the joint odd-part multisets agree for n<=4.

An elementary group argument independently excludes a tournament automorphism group of order divisible by 15 at n=5. Any such group has odd order dividing 120, hence order 15. A group of order 15 is cyclic (the Sylow counts force normal subgroups of orders 3 and 5), but S5 has no element of order 15. This proves the obstruction without trusting the complete census.

Correspondingly the matrix M modulo 15 is invertible, yet its support restricted to entries coprime to 15 has no perfect matching. Invertibility over a product of fields does not force a single determinant term to be a unit in every factor.

## 6. Exact reproduction and evidence

Files: `walsh_support.py`, `walsh_support.json`, `walsh_support.out`.

    python3 -B walsh_support.py --max-n 6 --out walsh_support.json

The script is independent of the inherited scripts and nauty. It enumerates all masks and all vertex permutations, builds signed graph orbits, rejects sign-odd stabilizers, and builds all tournament orbits. It verifies both weighted orthogonality identities exactly, determinant powers of two, rational and modular cofactor recursions, every prime-valuation inequality, and matching-preserved valuations for p=3,5,7. An independent augmenting-path routine checks support matchings. All 0<=n<=6 universes are complete, with N=1,1,1,2,4,12,56.

| n | det exponent | p=3 valuation histogram | p=5 valuation histogram |
|---|---:|---|---|
| 3 | 2 | 0:1, 1:1 | 0:2 |
| 4 | 8 | 0:2, 1:2 | 0:4 |
| 5 | 49 | 0:8, 1:4 | 0:11, 1:1 |
| 6 | 378 | 0:43, 1:12, 2:1 | 0:54, 1:2 |

Each histogram is identical on both sides. The determinant sign depends on the declared catalogue order and column sign gauges; the exponents are invariant. The script stores the complete matrices and explicit inverse-pair tables.

## 7. Exact remaining boundary

PROVED: all-order one-prime valuation distribution, a p-unit support-respecting constructive class bijection, and the abstract weighted-pairing generalization. REFUTED: simultaneous preservation of all odd-prime valuations, already at n=5. OPEN: a local graph operation implementing these matches without a complete atlas, preservation of a specified non-Sylow subgroup action, or a new meaningful naturalness criterion at the unlabelled level. The theorem neither supplies an equivariant lift nor contradicts THM-4531's scoped obstruction to such lifts.
