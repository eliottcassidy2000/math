# Shared-port cycle indices and a full tournament-count recurrence

**PROVED** W1-W4 below, using standard Burnside/Pólya mechanisms. **FINITE-EXACT**
checks are in [check_weighted.py](check_weighted.py) and [its JSON](check_weighted.json).
IDs are local to this file, not global canon theorem IDs.

Concurrent provenance: [Session94's weighted cycle carriers](../session94_recursive_ports_crt_clocks_20261029/ports/WEIGHTED_CYCLE_CARRIERS.md)
arrived during this session with the same cycle-index mechanism, an extension
to all odd cyclic alphabets, and an unordered-subset polynomial hierarchy.
We credit that independent, broader development; the calculations here are
an alternate-path check and feed the canonical attachment carrier below.

## W1. The full recursive inventory

Let W_k be the depth-k iterated C3 tournament (3^k vertices), A_k=Aut(W_k),
and C_k the centralizer of digit negation. The inherited decomposition is

    A_k = A_(k-1) wr C3,    C_k = C_(k-1) x A_(k-1).

In C_k the second factor acts diagonally on two mirrored subtrees, after
folding one by digit negation. Write Z_G(p)=|G|^-1 sum_g product_j p_j^m_j(g).
Then Z_A,0=Z_C,0=p_1 and, for every k>=1,

    Z_A,k(p) = [Z_A,k-1(p)^3 + 2 Z_A,k-1(p_3,p_6,...)]/3,
    Z_C,k(p) = Z_C,k-1(p) Z_A,k-1(p_1^2,p_2^2,...).       (W1)

Proof: a root rotation is identity with probability1/3, giving three
independent child permutations. Either nonidentity rotation traverses the
three children; its return product is uniform in A_(k-1), and each length-j
return cycle becomes length3j. In the centralizer, the mirrored copies have
identical cycle lengths, each appearing twice. Thus the substitution is
`p_j -> p_j^2`, **not** `p_j -> p_(2j)`. The latter would join two separate
cycles and changes which attachments are fixed. The base cases include k=0.

For q colour choices at every vertex this specializes to

    a_0(q)=c_0(q)=q,
    a_k(q)=(a_(k-1)(q)^3+2a_(k-1)(q))/3,
    c_k(q)=c_(k-1)(q) a_(k-1)(q^2).                     (W2)

For binary attachments with z marking the number of selected vertices,

    a_0(z)=1+z,
    a_k(z)=(a_(k-1)(z)^3+2a_(k-1)(z^3))/3.

To retain the two port weights separately, put

    b_0(z,w)=(1+z)(1+w),
    b_k(z,w)=(b_(k-1)(z,w)^3+2b_(k-1)(z^3,w^3))/3,
    c_0(z)=1+z,   c_k(z)=c_(k-1)(z)b_(k-1)(z,z).        (W3)

For separate weights on each first-nonzero-digit stratum, the fully marked
answer is `(1+z_root) product_(j=0..k-1) b_j(z_(j,+),z_(j,-))`.
The two arguments are folded and acted on by the **same** A_j. Multiplying
their separate orbit counts would discard their joint alignment.

| k | vertices | binary A_k attachment orbits | binary C_k attachment orbits |
|---|---:|---:|---:|
|0|1|2|2|
|1|3|4|8|
|2|9|24|192|
|3|27|4624|887808|

These are attachment-pattern orbits under the specified groups, not arbitrary
tournament isomorphism classes. For actual named ports, retain the inherited
folded LCA depths **and cyclic directions**, or the relevant stabilizer action.
Cycle indices need not distinguish permutation groups or port orbits.

## W2. Weighted variance retains the shared source

Let X_j count fixed vertices of a uniform A_j element. Independently for each
j, the centralizer fixed-point statistic with weights constant on its single
port orbits is

    Y_w = w_root + sum_(j=0..k-1) (w_(j,+)+w_(j,-)) X_j.

The inherited critical ternary law gives E X_j=1, Var X_j=2j. Hence, for all
real weights (including negative weights),

    E Y_w^2 = (w_root+sum_j(w_(j,+)+w_(j,-)))^2
              +sum_j 2j(w_(j,+)+w_(j,-))^2.             (W4)

This is also the weight sum over ordered two-port orbits, by weighted
Burnside; fixed weights must be invariant under C_k. Its (2k+1)-square Gram
matrix is `11^T + sum_(j>=1) 2j v_j v_j^T`, where v_j selects the two signs
at level j. Its rank is k for k>=1 (rank1 for k=0): the vectors 1,v_1,...,
v_(k-1) are independent because the root coordinate singles out1. With all
weights1, this gives `(2k+1)^2+4k(k-1)=8k^2+1`.

This small rank concerns this named invariant statistic only. It is not a
bound on general return-Gram rank, nor a reduction of all joint-port data.
If local root identity is sampled with probability p instead of1/3, the
fixed-population mean offspring becomes3p. Criticality is supplied by the
uniform automorphism measure. It is not implied by the digit3 in Collatz.

## W3. A000568 has a closed multivariate recurrence

For all n>=0 define the formal weighted Burnside sum

    F_n(x_1,x_3,...) = (1/n!) sum_(sigma in S_n, all cycles odd)
                         2^E(sigma) product_j x_j^m_j(sigma),
    E(m) = [sum_(j,l) gcd(j,l)m_j m_l - sum_j m_j]/2.    (W5)

The exponent counts unordered-pair orbits. Even cycles fix no tournaments:
their halfway power swaps opposite vertices and would reverse their arc.
An odd cycle inventory fixes exactly2^E tournaments. The all-order recursion is

    F_0=1,
    n F_n(x) = sum_(odd j<=n) 2^((j-1)/2) x_j
                  F_(n-j)((2^gcd(j,l) x_l)_(l odd)).    (W6)

Proof: choose the cycle containing a designated vertex. Its length is j;
the usual labelled-cycle factor after division by n! is1/n. It contributes
(j-1)/2 internal pair orbits. Every remaining l-cycle has gcd(j,l) orbits
of pairs crossing to it. These give exactly the substituted weights. No
unseen child-isomorphism census is an input. Finally A000568(n)=F_n(1,1,...).

The scalar sequence alone does not close under this proof: future cycles
test the entire gcd-weighted cycle inventory. Equation(W6) retains exactly
the coordinates this recursion uses; no minimality among all possible
counting algorithms is claimed. It is an alternative organization of the
classical Burnside sum, not a newly discovered formula for A000568.

Exact regeneration n=0..20 matches the independent June implementation
[a000568_exact_burnside.py](../../../04-computation/a000568_exact_burnside.py).
The memoized recurrence uses1929 states for that declared run. Relevant rows:

| n | all tournaments | strongly connected tournaments |
|---|---:|---:|
|11|903753248|884464590|
|12|154108311168|152310149735|
|13|48542114686912|48234782263293|

The second column is A000568. The third follows the unique ordered sequence
of strong components: `s_n=t_n-sum_(i=1..n-1)s_i t_(n-i)`, equivalently
`T(x)=1/(1-S(x))`. This is a structural composition recurrence. It does not
resolve the different task of reconstructing every class from compressed
metagraph states, or establish a runtime improvement over existing code.

For any odd-order permutation group, evaluating its cycle index by
`product p_j^m_j -> 2^E(m)` counts all labelled tournament decorations modulo
that group. Thus(W1) and(W6) meet at the **same pair-orbit evaluator**.
The cross-gcd terms prohibit substituting independent per-cycle scalars.

## W4. The hostile pair remains separated

We re-import the actual [81/27 constructor](../relation_tower_frobenius_ports_20261028/check_root_switch.py),
check equal67 labelled cyclic triples and462915 Hamilton paths, and enumerate
all4096 binary subsets under each full automorphism group. Weighted subset
orbit coefficients, in weights0..12, are

    81: 1,2,6,14,21,28,32,28,21,14,6,2,1   (total176)
    27: 1,6,18,38,63,84,92,84,63,38,18,6,1 (total512).

Independent cycle-index evaluation gives the same polynomials. Already the
one-selected-vertex coefficient2 versus6 distinguishes them. The invariant
keeps information lost by cyclic triples, but this one successful control
does not prove completeness of weighted spectra for arbitrary tournaments.
