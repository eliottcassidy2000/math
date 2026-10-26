# The complementary quotient: score parity and an exact pair-flip chain

**Status:** PROVED elementary specialization of Theorems A–B in `EVEN_GRAPH_OPERATION_ENVELOPE.md`; FINITE-EXACT controls specified below. No historical priority claim. This quotient is **not** tournament switching. It divides by degree-even edge reversals rather than cuts.

## 1. Literal source map, not a counting analogy

Put V=F_2^{E(K_n)}, B=Cut(K_n), and Z=B^perp, the space of degree-even edge sets. For a tournament T write

    p(T)=(d_T^+(v) mod2)_(v in [n]),
    m=binom(n,2),
    P_n={p in F_2^n : sum_v p_v=m mod2}.

Reversing an edge uv toggles exactly p_u and p_v. Therefore reversing the edge set F changes p by its degree-parity boundary ∂F. The kernel of ∂ is Z, and for n≥1 its image is all even-weight vectors: a path/star pairs any even-sized boundary set. Consequently

    Tour/Z ≅ P_n

as affine S_n-sets. This is an actual equivariant correspondence. The n=0 case consists of one object on both sides. S_n acts on P_n by coordinate permutation, so the unlabelled quotient types are exactly

    r=|p|,  0≤r≤n,  r≡m (mod2).

No tournament with a distinguished ordering is retained or needed for this observable. This is different from the inherited THM-1470 score-parity law *within a cut-switching class*: there the cut action changes the parity vector; here equality of parity vectors defines the cycle-space quotient.

## 2. The exact weighted operation

On these types retain one contribution for every named arc reversal. The matrix A has entries

    A[r,r−2]=binom(r,2),
    A[r,r+2]=binom(n−r,2),
    A[r,r]=r(n−r),

with unavailable indices omitted. Choosing two odd vertices lowers their count by two, choosing two even vertices raises it by two, and a mixed pair leaves the count unchanged. Row sums are m. Loops and all multiplicities are part of this operator. It is reversible for weights proportional to binom(n,r); for n≥1 the probability normalization is binom(n,r)/2^(n−1) on the allowed parity sector.

The same source operation thus has three different faithful quotient operators: ordinary tournament isomorphism, cut-switching, and score parity. Nothing was inferred merely from matching state counts.

## 3. Full spectrum, including the exceptional middle mode

**Theorem S.** Let L=mI−A. Its complete eigenvalue multiset is

    {2a(n−a) : 0≤a≤floor(n/2),
                  excluding a=n/2 when a is odd}.

Every listed eigenvalue is simple. For n≥3,

    gap(L)=2(n−1),
    gap_algebraic(A/m)=4/n,
    gap((I+A/m)/2)=2/n,
    relaxation_time((I+A/m)/2)=n/2.

For n=0,1,2 there is one state and no positive gap is claimed. A concrete Fiedler observable is `n−2r`, where r is the number of odd outdegrees.

### Harmonic proof through the sign-even graph theorem

The dual of V/Z is B. Every nonzero cut graph is K_(a,n−a), where the two part sizes are considered without order. If a<n−a, every graph automorphism preserves the two degree-distinguished parts; orienting all cross-edges from one part to the other proves sign-evenness. In the balanced case n=2a, a part-swap reverses all a² edges. Thus the balanced cut is sign-even iff a is even. Permutations internal to the parts reverse no oriented cross-edge, so this condition is sufficient as well as necessary.

By Theorem A, these and only these cut types index invariant Fourier lines. They contain a(n−a) edges. The edge-flip Laplacian therefore has eigenvalue 2a(n−a) on each line. The function a(n−a) is strictly increasing for integer 0≤a≤floor(n/2), so all surviving eigenvalues are simple. For n≥3, a=1 survives and gives the least positive value 2(n−1). This proves all spectrum and gap assertions.

In particular the type count is

    floor(n/2)+1−[n≡2 (mod4)],

which equals the number of allowed weights r. The middle exclusion is forced by stabilizer phase, not an ad hoc correction to a count formula.

### Independent explicit eigenvectors

For each listed a define

    K_a(r)=[z^a](1−z)^r(1+z)^(n−r)
          =sum_j (−1)^j binom(r,j) binom(n−r,a−j).

This is the sum of the characters `(-1)^(sum_(v in A) p_v)` over all a-element vertex subsets A. Each pair flip has eigenvalue `m−2a(n−a)` on any one of these characters, so the same holds for K_a. On P_n, complementary character subsets satisfy

    chi_(A^c)=(-1)^m chi_A.

At balanced a=n/2, this forces cancellation exactly when m is odd, equivalently a is odd. Thus the excluded K_(n/2) is identically zero **on the allowed parity sector**. Other K_a can have individual zero entries (for example K_1(n/2)=0), but are not identically zero there. Distinct surviving cuts give the independent invariant lines established above. The first nonconstant polynomial is K_1(r)=n−2r.

## 4. A genuine three-regime consequence of one theorem

| Source quotient of the SAME arc-flip operation | Surviving dual graphs | First nonconstant mode | Unnormalized gap |
|---|---|---|---:|
| Tour/isomorphism | sign-even graphs | P3 plus isolates; centered cyclicity | 4, n≥3 (inherited THM-588) |
| (Tour/B)/isomorphism | sign-even Eulerian graphs | C4 plus isolates; centered tr(K^4) | 8, n≥4 |
| (Tour/Z)/isomorphism | sign-even cuts | K_(1,n−1); score-parity imbalance | 2(n−1), n≥3 |

The graph predicates and operations have not been swapped. These are **harmonic filters** of one labelled arc-flip operator. In particular B and Z are generally not nested. For odd n they are complementary submodules V=B⊕Z; neither quotient is asserted to refine the other, and no monotonicity comparison between their gaps is inferred. Both are quotients of the original tournament carrier, which is the valid common source.

At n=5 the last two gaps happen to equal eight, but their eigenfunctions are not the same: C4 and K_(1,4) have disjoint Fourier support. Indeed for odd n the two quotient invariant spaces intersect only in the constants, because B∩Z=0 on the dual side. Equal spectral gaps alone do not give an equivalence of observables or operations.

## 5. Exact finite evidence and boundaries

    python3 -B verify_score_parity.py > score-normal.json
    python3 -B -O verify_score_parity.py > score-optimized.json
    cmp score-normal.json score-optimized.json

The independent standard-library verifier uses no producer imports, assertions or floating point. It records **116,047 exact gates**. It enumerates every labelled tournament at n=0,...,6, checks the boundary map and all fibres, and checks the direct pair-flip rates for every allowed labelled parity vector. Separately it constructs the tridiagonal operator and every Krawtchouk eigenvector for n=0,...,64. Each eigenmatrix has an explicit full-rank certificate modulo the prime 1,000,000,007; a nonzero modular determinant/rank proves rational independence. Exact integer equations verify diagonalization and reversibility. The normal and optimized outputs byte-match.

Examples:

| n | allowed parity weights r | exact L spectrum |
|---:|---|---|
| 3 | 1,3 | 0,4 |
| 4 | 0,2,4 | 0,6,8 |
| 5 | 0,2,4 | 0,8,12 |
| 6 | 1,3,5 | 0,10,16 |
| 7 | 1,3,5,7 | 0,12,20,24 |
| 8 | 0,2,4,6,8 | 0,14,24,30,32 |

The all-order proof is the character/cut argument, not extrapolation from this table. This result does not reconstruct scores, Hamiltonian paths, cyclicity, the switching class, or the tournament type from score parity. Its source quotient is deliberately coarse, and the theorem states exactly the operation that remains valid after that loss.
