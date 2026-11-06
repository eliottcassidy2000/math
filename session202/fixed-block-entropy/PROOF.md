# Arbitrary root incidences do not rescue one fixed strong block

Status: **PROVED, independently audited in the stated scope.** The all-order argument below is self-contained. `verify.py` supplies finite literal controls, not an all-target computation. No new eventual-coverage threshold has been proved.

## 1. Question and inheritance

The audited Session 201 twelve-template occupancy lemma is sufficient for every allowed odd target, but remains open. I tested a different physical reconstruction: keep one marked root, concatenate arbitrary numbers of copies of a strong block, and choose every root incidence independently. This is much freer than repeating a fixed incidence word or asking for an affine append clock.

The result is a quantitative obstruction for this entire reconstruction class. The source of the obstruction is the exact number of ordered two-path covers of a block, not a guessed affine recurrence. It is complementary to the fixed-apex regular-language entropy theorem in `session197/entropy-bridge/PROOF.md`, and uses neither its external strong-order bound nor a restriction on incidence words. The final multiplicative-closure argument is the same elementary Dirichlet-product mechanism as its Lemma 3, reproved below.

The closest positive tools are the actual root-amalgamation formulas in `session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md`. They also identify exactly how arbitrary O/Z singleton spacers escape the present obstruction. The hostile near miss is assuming that arbitrarily chosen root incidences already give enough arithmetic freedom. Live concepts are physical amalgamation, two-bin growth, scalar support entropy, strong-component intervals, and prime-atom occupancy.

## 2. Actual graph family and two-bin constant

Fix a strong tournament S on s>=3 vertices. Let h=H(S), and put

    beta=P_S(2)=sum_(A subset V(S)) H(S[A]) H(S[V(S)\A]),

with H(empty)=1. Thus beta counts ordered pairs of possibly empty directed paths partitioning S.

For every m>=0 take m disjoint labelled copies S_1,...,S_m, orient all edges S_i -> S_j for i<j, and adjoin one root r. Every one of the sm root incidences is freely chosen. Denote this family by G_S. The root may be chosen differently for every target. No fixed incidence blocks, automaton, periodicity, or seed is imposed.

Every induced tournament has a Hamiltonian path by the usual insertion argument, so every term in P_S(2) is at least one. The terms A=empty and A=V(S) both equal h. Strongness implies h>=s>=3: every specified vertex can be the endpoint of a Hamilton path. To see this without a cycle theorem, take a longest path ending at that vertex. An outside vertex beating a path vertex can be inserted immediately before the first path vertex it beats. If no outside vertex can be inserted, every path vertex beats every outside vertex, contradicting strongness. Reversing arcs proves the analogous starting statement.

Consequently

    beta >= 2^s+2(h-1) > 2^s.                         (1)

Define the exact exponent

    rho = log(2^s)/log(beta),       0<rho<1.            (2)

## 3. Strong outputs have two-bin exponential growth

**Lemma 1.** If a graph in G_S on m>=2 copies is strong, then

    H(T) >= beta^(m-2).                               (3)

**Proof.** The root must have an out-neighbor in S_1. Otherwise nothing outside S_1 points into S_1, contradicting strongness. Likewise it must have an in-neighbor in S_m. Choose a Hamilton path P_1 of S_1 starting at one such out-neighbor, and a Hamilton path P_m of S_m ending at one such in-neighbor.

For each interior block S_i, 2<=i<=m-1, choose an arbitrary ordered two-bin path cover (A_i,B_i), including empty bins. Form the vertex word

    A_2 A_3 ... A_(m-1) P_m r P_1 B_2 B_3 ... B_(m-1),

omitting empty pieces. Every edge between consecutive different blocks goes forward in block index. The two edges at the root are legal by construction. Thus this is a Hamilton path. Its restriction to each interior block recovers the ordered cover, so the beta^(m-2) choices give distinct paths. QED.

This argument retains every vertex and every prescribed root incidence. No formal profile is assumed realizable.

## 4. Weak outputs, arbitrary padding, and scalar support

Every S_i remains inside one strong component of T. The component containing r either is {r}, or contains precisely a consecutive interval S_a,...,S_b. Indeed, if it contains S_a and S_b, every intermediate block has paths to and from that component using the forward cross edges. No other component can contain two blocks, since a path back across block indices must use r. Hence all other components are whole copies of S.

Let ell be the number of blocks in the root component C. Then

    H(T)=h^(m-ell) H(C),                              (4)

where H(C)=1 when ell=0. This is an actual SCC factorization. Hamilton paths traverse tournament SCCs contiguously in their total order; their choices multiply.

Let B_S(X) count distinct H(C)<=X for strong root-containing graphs made from at least one copy. There are at most 2^(s ell) incidence presentations at length ell. Lemma 1 implies ell<=2+log_beta X when ell>=2; the same bound includes ell=1. Summing the geometric series gives the explicit estimate

    B_S(X) <= [2^(3s)/(2^s-1)] X^rho,       X>=1.     (5)

We count presentations only as an upper bound; there is no injectivity claim from presentations to values or isomorphism classes.

**Theorem 2.** The number of distinct scalar counts <=X generated by G_S is O_S(X^rho). In particular this family misses infinitely many odd integers.

**Proof.** By (4), every value is h^j or h^j b for some j>=0 and strong-root count b. The pure powers contribute O(log X). By (5), the other values contribute at most

    [2^(3s)/(2^s-1)] X^rho sum_(j>=0) h^(-j rho)
    = [2^(3s)/((2^s-1)(1-h^(-rho)))] X^rho.

The logarithmic term is O(X^rho) since rho>0. QED.

Adding arbitrary universal source/sink vertices outside a graph changes neither its count nor this conclusion. Arbitrarily many neutral vertices do not create a hidden order cutoff problem.

## 5. Arbitrary ordered joins still do not cover a tail

Let M_S be the scalar image after arbitrary finite ordered joins of members of G_S. These are real graph constructions and multiply H. For every epsilon>0,

    #{n in M_S : n<=X}=O_(S,epsilon)(X^(rho+epsilon)). (6)

Here is the elementary argument. Let A be the set of distinct generator counts at least two. For every t>rho, dyadic grouping and Theorem 2 show sum_(a in A) a^(-t)<infinity. Therefore

    Z(t)=product_(a in A) (1-a^(-t))^(-1)<infinity.

Expanding finite products and taking increasing limits overcounts all distinct products, because different exponent vectors can coincide. Thus sum_(n in M_S) n^(-t)<=Z(t), and the number of n<=X is at most Z(t)X^t. Take t=rho+epsilon. This uses unordered exponent vectors, not an incorrect count of ordered multiplication words.

Choosing epsilon<1-rho proves density zero. In fact sum_(n in M_S)1/n converges. Euler's divergence of the prime reciprocal sum then implies that M_S omits infinitely many odd primes. The density-zero conclusion itself needs no prime-distribution result.

## 6. Concrete physical cases and the scope boundary

1. For S=C3, rows (2,4,1), direct counting gives h=3 and beta=12. Hence rho=3 log2/log12, approximately 0.8368288369. All eight root patterns are allowed at every block, with no restriction on their chronology.
2. For the strong four-vertex tournament S with rows (8,1,3,6), h=5 and beta=32. Thus rho=4/5 exactly. Its entire arbitrary-incidence shared-root family has O(X^(4/5)) counts. Its arbitrary ordered-join closure has O_epsilon(X^(4/5+epsilon)) counts.

These classes are not equivalent to a fixed number of apexes above a transitive base: deleting the root leaves arbitrarily many nontransitive strong blocks.

**Necessary scope restriction.** This proof does not allow arbitrary binary singleton blocks between the copies. A singleton has h=1 and beta=2, while it has two root-incidence choices. Its entropy/growth ratio is exactly critical, so inequality (1) loses its strict gap. The audited prefix/core/suffix constructions permit arbitrarily many such singleton decisions, and are not obstructed here. Nor does the proof cover moving the root, editing old nonroot arcs, or an unrestricted bank of different strong block types.

The result therefore does not refute the twelve-template lemma, C1 eventual coverage, or the target threshold 23. It prevents a specific apparent shortcut: replacing the difficult binary source coordinates by arbitrary incidence choices on repeated copies of one stronger fixed block.

## 7. A precise surviving sufficient arithmetic problem

For comparison, the inherited arbitrary-core root-amalgamation compiler offers a genuine route outside the twelve central forms. For an actual non-singleton rooted core C, put

    h_C=H(C), H=H(C-r), beta=P_(C-r)(2),
    lambda=L_C-H, alpha=R_C-H,
    sigma=H-start_r(C), eta=H-end_r(C).

Here R_C,L_C are the labelled two-bin responses whose designated nonempty path starts or ends at r. Set F(q)=g(q+1) for q>=0, so

    F(0)=1, F(2q+1)=F(q),
    F(2q)=F(q)+q+(q mod2) for q>=1.

The audited formula and literal compiler give

    N_C(q,t)=h_C+H(F(q)+F(t)-2)
              +lambda q+sigma(q mod2)+beta qt
              +alpha t+eta(t mod2),       q,t>=0.     (7)

Construct the binary prefix from the rooted singleton, root-amalgamate C, and append the binary suffix, using the exact conventions in `session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md`, Sections 4-5. The graph order is |C|+bit_length(q)+bit_length(t). The two indices and their F-values are not independent coordinates.

**OPEN sufficient prime-atom occupancy assertion.** Choose a specified finite bank of actual cores. If its values (7), together with the no-core binary family, represent every odd prime except 7 and every 7p for odd primes p>=5, then every positive odd integer except 7 and 21 is a tournament path count.

This sufficiency uses the inherited multiplicative prime-atom reduction, not a new inference from the entropy theorem. For completeness, explicit C1 source words give the additional seeds

    49:  OJJOOJ,
    63:  OZZOOJ,
    343: OOJOZJJOJ.

The last J closes the source. Start at 0->1, root 0; O adds a sink, Z and J add a vertex beating only the old root, and J marks the newborn.

To check the reduction, write an odd integer as 7^a 3^b r with r coprime to 21. If a=0, multiply its allowed prime counts. If a>=2, factor 49 and 343 to remove the entire 7-power, and use the other primes. If a=1 and r>1, take one prime p|r and factor off the represented 7p. If a=1 and r=1, the cases b=0,1 are 7,21, and b>=2 uses 63 times 3^(b-2). Ordered joins implement every multiplication physically.

A finite range of prime atoms is not a proof of this assertion. This session has not proved it for any finite bank. Its purpose here is to state the remaining sufficient arithmetic obligation after the fixed-block shortcut has failed, without converting modular or approximate coverage into exact occupancy.

## 8. Reproduction and research outcome

Run in the packet directory:

    python3 verify.py > VALIDATION-normal.json
    python3 -O verify.py > VALIDATION-optimized.json
    cmp VALIDATION-normal.json VALIDATION-optimized.json

The verifier checks all root incidences on C3 blocks for m=0,1,2,3 and on the four-vertex strong block for m=0,1,2, plus a declared deterministic 38-presentation sample on three four-vertex blocks. It uses literal adjacency matrices, subset Hamilton-path DP, and independent graph reachability for SCCs. It checks the two-bin constants, the strong lower bound, contiguous root components and the exact weak factorization on 896 presentations. The seed compiler is checked separately on the three displayed words. Checks raise exceptions and remain active under optimization.

The exploratory files `probe.py` and `probe_parity.py` tested constant differences between small O/Z/J central matrices. They found no nonzero constant-offset pair through length twelve in their declared word census. This observation is not used in any theorem and does not imply an all-word obstruction.

**Outcome:** a fixed strong block with unrestricted root incidences cannot prove eventual coverage, even after arbitrary products. The positive two-index occupancy question and the possibility of threshold 23 remain open.
