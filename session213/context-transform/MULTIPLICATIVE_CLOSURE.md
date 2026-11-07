# SCC products and a finite seed bank do not repair the sparse cycle code

Status: PROVED elementary consequence; independent proof audit PASS. The counting lemma is a standard Euler-product majorant argument; no literature-priority claim. The application combines `PROOF.md` with the inherited dominance-join product law. It does not constrain tournaments outside the stated construction closure.

## 1. A sublinear atomic spectrum stays sublinear under all finite products

Let A be a set of positive integers, put A(X)=#{a in A:a<=X}, and assume

    A(X)=O(X^alpha),      0<alpha<1.

Write M(A) for the set of all finite products of elements of A, including the empty product 1. Duplicates and factors equal to 1 are ignored when listing generators.

**Lemma.** For every real beta>alpha,

    #{n in M(A):n<=X}=O_beta(X^beta).                  (1)

In particular M(A) has natural density zero, by choosing alpha<beta<1. Adjoining any finite set of positive integer generators does not change this conclusion.

**Proof.** Dyadic shells give

    sum_(a in A,a>1) a^(−beta)
      <= sum_(j>=0) #{a in A:2^j<a<=2^(j+1)} 2^(−j beta)
      = O(sum_(j>=0) 2^(j(alpha−beta))) < infinity.

As a>=2, the elementary inequality

    −log(1−a^(−beta)) <= a^(−beta)/(1−2^(−beta))

shows convergence of the positive Euler-product majorant

    E_beta=product_(a in A,a>1) (1−a^(−beta))^(−1) < infinity.   (2)

For each finite subset F of the generators, expanding its factors as geometric series gives

    product_(a in F)(1−a^(−beta))^(−1)
      =sum_(e in N_0^F) (product_(a in F) a^e_a)^(−beta).

All terms are nonnegative. Increasing F through all generators therefore gives the sum over all finitely supported exponent vectors. Distinct vectors may have the same product; this only increases the sum. Every member of M(A) has at least one such vector. Consequently

    sum_(n in M(A)) n^(−beta) <= E_beta.

For n<=X one has n^(−beta)>=X^(−beta), so (1) follows with constant E_beta. A finite extra generator set adds finitely many convergent factors to (2). QED.

The lemma does not assert the endpoint bound O(X^alpha), which need not follow from this argument. Nor does it assume that the generators are prime or multiplicatively independent. Allowing collisions is exactly why the Euler product is a majorant rather than an equality for distinct integers.

## 2. Application to the binary cycle family

Let F be the complete four-variant Session212 comb family, with all lengths and leading-zero padding. Fix any finite menu of deterministic marked recipes Phi_i, each retaining its input comb as an induced subtournament. Include the identity recipe if the original combs are also to be allowed. Let S be an arbitrary fixed finite collection of further seed tournaments.

Let C be the family obtained by taking any finite number of tournaments from S or from any Phi_i(F), arranging them in any dominance order, and orienting every arc from an earlier block to a later block. Arbitrary source and sink padding is allowed.

**Corollary.** With alpha=log(2)/log(3), for every beta>alpha,

    #{H(T):T in C, H(T)<=X}=O_beta(X^beta).            (3)

Hence this entire construction closure has density-zero path support and cannot realize an odd tail.

**Proof.** `PROOF.md`, Theorems 2 and Corollary 3, bound the union of the atomic path-count spectra by O(X^alpha). A finite seed collection adds finitely many values. Dominance-join multiplicativity gives

    H(T1 -> T2 -> ... -> Tr)=product_i H(Ti).

Indeed a Hamiltonian path cannot leave a later block and return to an earlier one, so it is exactly the concatenation of one Hamiltonian path per block. Conversely every such concatenation is valid. Thus the output values lie in M(A) for the atomic count set A, and the lemma applies. A singleton block has H=1, so source/sink padding changes no value. QED.

**Precise surviving boundary.** This excludes arbitrary *products of the specified atoms*, even with arbitrarily many factors. It does not exclude interleaving the blocks, changing cross-arcs, making a new connected interaction between independent comb codes, applying an operation to a product that is not already one of the fixed recipes on a single comb, or an infinite independent seed bank. Those operations can create new atomic values and are not covered by the majorant.

## 3. Dependencies and verification

- `PROOF.md`: the all-padding O(X^(log_3 2)) atomic count bound and finite-menu induced-extension corollary.
- `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`, section 1, citing `THM-1862-order-join-reduction-principle`: the inherited dominance-join multiplicativity. Its source identity is in `SOURCES.json`; the short concatenation proof above makes the application self-contained.

The all-order conclusion is the positive-series proof, not an extrapolation from finite product enumeration. No new computation is needed beyond the audited atomic bound. The important audit checks are the direction of the majorant, repeated factorisations, exclusion of generator 1 from the Euler product, and keeping beta strictly above alpha.
