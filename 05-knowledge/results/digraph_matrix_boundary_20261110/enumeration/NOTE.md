# Pair states, cancellation fibres and Hamilton path marked classes

**PROVED elementary refinements of inherited Burnside and fixed-path
mechanisms; FINITE-EXACT independent checks through order four.** The
contribution is a degree-resolved interface between raw binary matrices,
their cancellation semantics and their Hamilton-path marked classes. It
does not classify the Hamilton-count spectrum or assert literature priority.

The scalar palette enumeration is already in
[THM-4720, k-valued tournaments](../../../../01-canon/theorems/THM-4720-k-valued-tournaments-by-burnside-the-parity-dichotomy-of-the-tie-state-and-three-hypertournaments.md),
[THM-4845, Gaussian cycle-index reformulation](../../../../01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md),
and [THM-5097, integral cycle-response recurrence](../../../../01-canon/theorems/THM-5097-factorial-gauge-integral-cycle-response-recurrence.md)
with its [palette proof](../../session239_sequence_engines_20261110/burnside/PALETTE_PROOF.md).
The matrix dictionary and fixed-path tiling gauge are inherited from
[THM-5041](../../../../01-canon/theorems/THM-5041-matrix-tiling-pair-dictionary.md)
and its cited THM-830. The following arguments retain orbit lengths and
the raw symmetric-pair decoration rather than only the number of colours.

## 1. Three degrees refine the five familiar class counts

Let D be a loopless digraph on n vertices. For each unordered pair there
are two asymmetric states, one absent state 00 and one reciprocal state11.
Write m(D), b(D), s(D) for their absent, reciprocal and asymmetric pair
counts. Thus m+b+s=N, where N=binomial(n,2). Define

    F_n(x,y,z) = sum over isomorphism types [D] of x^m(D) y^b(D) z^s(D).

Relabelling reverses the two asymmetric states when it reverses the pair,
while fixing both00 and11. This is the four-state reversal palette(q,t)=(4,2).
It is not the k=4 palette in THM-4720, which has two exchanged colour pairs
and no fixed state. The scalar palette formula is inherited; its weighted
version follows by retaining the number of pairs in every orbit.

For a partition lambda=(L_1,...,L_r) of n let z_lambda be the permutation
centralizer size, prod_l l^a_l a_l!. Define two multisets of orbit lengths:

- P_lambda contains floor((L_i-1)/2) copies of L_i for every cycle, and
  gcd(L_i,L_j) copies of lcm(L_i,L_j) for every i<j.
- R_lambda contains one copy of L_i/2 for each even L_i.

The first multiset contains pair orbits whose full return preserves the
endpoints; the second contains antipodal pair orbits whose return swaps
the endpoints. Then

    F_n(x,y,z) = sum_(lambda |- n) 1/z_lambda
       prod_(ell in P_lambda) (x^ell+y^ell+2z^ell)
       prod_(ell in R_lambda) (x^ell+y^ell).                 (1)

**Proof.** On a preserving orbit any of the four pair states can be chosen
at one position and is then propagated. Both asymmetric choices contribute
z^ell. On a reversing orbit the chosen state must be fixed by reversal, so
only00 and11 survive. Every fixed matrix has exactly one such assignment.
Weighted Burnside gives(1). The within-cycle and cross-cycle orbit lists
are those of THM-4720 with their lengths retained. QED.

The precise specializations are:

| Family up to isomorphism | Evaluation | Sequence | Terms at n=1,2,3,4,5 |
|---|---|---|---|
| Tournaments | F_n(0,0,1) | [A000568](https://oeis.org/A000568) |1,1,2,4,12|
| Oriented graphs | F_n(1,0,1) | [A001174](https://oeis.org/A001174) |1,2,7,42,582|
| Semicomplete digraphs | F_n(0,1,1) | A001174 |1,2,7,42,582|
| General loopless digraphs | F_n(1,1,1) | [A000273](https://oeis.org/A000273) |1,3,16,218,9608|
| Symmetric loopless matrices | F_n(1,1,0) | [A000088](https://oeis.org/A000088) |1,2,4,11,34|

The tournament boundary kills every reversing orbit; equivalently every
even vertex cycle has zero contribution. This is the structural reason
that Davis's A000568 sum runs only over odd cycle lengths. In the larger
families the symmetric pair states admit these otherwise impossible
symmetries. The sequence identifications are CITED to the linked OEIS
entries; the displayed terms are independently computed here.

Replacing 00 by11 and11 by00, while retaining each asymmetric direction,
commutes with relabelling. Hence F_n is symmetric in x,y. This explains the
oriented/semicomplete equinumeracy, but does not preserve ordinary H. For
example the empty order-two graph has H=0 and its reciprocal lift has H=2.

To retain loops, multiply each summand of(1) by
prod_i(1+t^L_i): a permutation-fixed diagonal is constant on each vertex
cycle. This gives a fourth degree counting looped vertices. At x=y=z=1 its
coefficients are [A384105](https://oeis.org/A384105), and the row sum is
[A000595](https://oeis.org/A000595), the number of binary relations up to
simultaneous row/column relabelling. At n=4 the loop-degree row is
218,752,1104,752,218, totalling3044. Loops do not alter Hamiltonian paths,
but these coefficients show that their isomorphism classes do not multiply
by2^n. Arbitrary independent row and column permutations are not the graph
isomorphism relation used here.

## 2. Cancellation changes both the object and the quotient

Define the effective oriented graph

    C(A)_ij=A_ij(1-A_ji) for i!=j,   C(A)_ii=0.           (2)

Thus 10 and01 retain their single arc, while both00 and11 become absent.
This map is idempotent and commutes with simultaneous relabelling.
Its fixed loopless matrices are exactly the oriented graphs. Restricted
to raw semicomplete matrices, it is a bijection onto oriented graphs:
the inverse fills every effective missing pair with11. Tournaments are
the common no-defect boundary, with no effective missing pair at all.

For a fixed labelled oriented graph O with m missing pairs, C has exactly
2^m loopless raw preimages, and2^(m+n) raw preimages when loops are allowed.
To count these preimages up to isomorphism, retain Aut(O)'s action on its
missing pairs and vertices. Let M(O) be its missing-pair set. Then the
polynomial counting raw preimage types by reciprocal pairs and loops is

    Phi_O(y,t) = 1/|Aut(O)| sum_(g in Aut(O))
       prod_(Q orbit of g on M(O)) (1+y^|Q|)
       prod_(V orbit of g on vertices) (1+t^|V|).          (3)

For loopless preimages put t=0. **Proof.** A raw preimage is exactly a subset
of M(O) marked11 instead of00 and a subset of vertices marked with loops.
An isomorphism between two such raw graphs must preserve their common
cancellation O, hence belongs to Aut(O). Conversely every automorphism of
O mapping the two decorations is a raw isomorphism. Burnside on these two
sets gives(3). QED.

This proves the finite decomposition

    A000273(n) = sum_[O oriented of order n] Phi_O(1,0),
    A000595(n) = sum_[O oriented of order n] Phi_O(1,1).   (4)

Replacing the fibre by a constant, by2^m after quotienting, or by2^m/|Aut(O)|
would generally give wrong answers. The group action is the missing sidecar.

The smallest positive-path control used here is the directed C4, with
out-neighbour bitmasks(2,4,8,1). It has H=4 and two missing diagonal pairs.
Its four labelled loopless raw preimages yield only three isomorphism types:
none, one, or both missing pairs changed to11. Thus Phi_O(y,0)=1+y+y^2.
Their ordinary H values are4,4,8, while all three cancelled totals remain4.
With loops, its64 labelled preimages give22 raw types. In particular even
strong connectivity does not make the first reciprocal addition increase H.

Under cancellation semantics the statement “H=7 or21 needs an effective
missing pair” follows from tournament exclusion. In a raw matrix that
missing pair can be00 **or11**. A claim requiring raw00 alone has changed
meaning and does not follow. Diagonal bits remain irrelevant to the count.

## 3. Marking one ordered Hamilton path removes all automorphisms

Let H(D) count ordered spanning paths under ordinary arc semantics and
H_c(D)=H(C(D)) under cancellation semantics. For n>=1, put k=N-n+1.
Then the following identities are exact polynomials:

    sum_[D] H(D)/|Aut(D)| x^m y^b z^s
        = (y+z)^(n-1) (x+y+2z)^k,                       (5)

    sum_[D] H_c(D)/|Aut(D)| x^m y^b z^s
        = z^(n-1) (x+y+2z)^k.                          (6)

**Proof.** An automorphism fixing an ordered Hamilton path fixes every vertex,
so its action on the path set is free. Consequently H/|Aut| is an integer,
counting isomorphism classes of pairs(D,P) with one selected ordered path.
Relabel P uniquely as1,2,...,n. For each of its n-1 consecutive pairs,
ordinary semantics permits either the forward single arc or11, contributing
z+y. Cancellation semantics permits only the forward single arc, contributing
z. Every other pair independently contributes x+y+2z. Conversely these
choices construct exactly all graph/path pairs in this fixed gauge. QED.

This is the inherited Hamilton-frame/tiling argument with the two defect
types retained, not a new canonical labelling algorithm. If diagonal bits
and raw loop-preserving isomorphisms are included, multiply each right hand
side by(1+t)^n. The selected ordered path fixes vertex positions, unlike
the unmarked class count in(1), so all loop bits now are independent.

Some useful consequences, always for uniform **labelled** graphs:

| Family | Expected ordinary H |
|---|---|
| Tournaments | n!/2^(n-1) |
| General loopless digraphs | n!/2^(n-1) |
| Oriented graphs | n!/3^(n-1) |
| Semicomplete digraphs | n! (2/3)^(n-1) |

For uniformly sampled raw binary matrices the expected cancelled count is
n!/4^(n-1). Conditioning away raw00 gives uniform oriented semantic images
and expected cancelled count n!/3^(n-1). These distinctions are consequences
of the elementary path-marked identities, not assumptions of independence
between different Hamilton paths.

The totals of isomorphism classes F_n and the automorphism-weighted path
mass in(5) answer different questions. Neither determines the support of H,
let alone whether a given odd integer occurs for tournaments. Equal scalar
class counts for oriented and semicomplete graphs coexist with different
path masses and spectra. Equation(1) generates related counting sequences;
equations(3),(5),(6) specify the extra data needed to interpret them.

## 4. Exact verification and reproduction

Run from the repository root:

    python3 -B 05-knowledge/results/digraph_matrix_boundary_20261110/enumeration/check.py > normal.json
    python3 -B -O 05-knowledge/results/digraph_matrix_boundary_20261110/enumeration/check.py > optimized.json

Saved `normal.json` and `optimized.json` are byte-identical. The checker:

- enumerates all1,4,64,4096 loopless labelled matrices at orders1 through4;
- canonicalizes under every vertex permutation and checks orbit-stabilizer;
- compares the complete three-degree orbit polynomial with(1), and both
  path-marked polynomials with(5),(6), using literal paths and subset DP;
- checks all139 partition pair-orbit profiles through order10 by directly
  following both endpoints around the permutation;
- checks104 oriented-source fibre polynomials, with and without loops,
  by explicit decoration and full canonicalization;
- recovers A000595 and its diagonal-degree refinement from these fibres;
- retains the C4 collision, the empty/singleton boundary in the sequence
  engine, and all source hashes.

This finite universe verifies the implementation independently of the
displayed algebra. The all-order statements are supplied by the proofs.
