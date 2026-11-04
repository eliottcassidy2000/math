# Exact low-degree type dependence for every coarse modulo-eight response

## 1. Domain and exact run coefficients

Let Q be ANY tournament on m nonempty child blocks A_i. Attach actual ports x->y, uniformly on each whole child, with quotient masks a,c. The statement concerns the complete Hamilton count of this physical extension, not a truncation in its port bits.

Divide each child's complete ordered-path-cover profile by its odd Hamilton count modulo8. The five-type theorem gives

    r_i=(0,1, 2-2s_i+4b_i, 2c_i+4b_i,0,...),             (1)

with allowed actual states

    S:(s,b,c)=(1,0,0), P:(0,0,0), E:(0,1,0),
    O:(0,1,1), C:(0,0,1).

Here S is a singleton, P a pair, E an even-order child of order at least four, and O/C are odd-order children of order at least three with even/odd Hamilton-cycle parity. The odd H factors are retained separately.

For a positive visit vector k=(k_1,...,k_m), let K(k) count directed run words of the augmented quotient in which child i occurs k_i times and each of the two ports occurs once. Adjacent occurrences of the same child are forbidden, because they would be a single maximal run. Every other adjacent pair must be a directed arc of the actual augmented quotient. Prefixes and suffixes can start or end at any symbol, and both port orders are included when their actual intervening arcs allow them.

Filling a run word contributes exactly product_i q_(k_i)(A_i), one ordered child cover per occurrence list. Thus

    F_(Q[A_i])(a,c) = product_i H(A_i) · Phi_Q(a,c;r) (mod8),
    Phi_Q = sum_k K(k) product_i r_(i,k_i).               (2)

No physical child is relabelled or replaced in the comparison. Small model graphs introduced later are only a way to compute the universal coefficients.

## 2. A two-return involution, including a second repeated child

Suppose child i occurs exactly three times. A word has the unique form

    prefix i L1 i L2 i suffix,

where L1,L2 contain no i and are nonempty. Exchanging L1,L2 preserves every adjacent directed arc and every symbol multiplicity. Hence it is an involution on these words.

A fixed word would have L1=L2. Each return segment contains at least two DISTINCT non-i symbols. A single-symbol segment would require a directed two-cycle i->j->i, impossible in a tournament. A segment using only repeated copies of one symbol is also impossible, since adjacent equal child labels are forbidden and the ports occur only once.

If every other child occurs once, fixed points are therefore impossible. The same is true if exactly ONE other child j is repeated, whether twice or three times: two identical return segments would require at least two other repeated symbols, while only j is available. The two ports cannot supply a repeated symbol.

Consequently the integer coefficients satisfy

    K(i3) is even,
    K(i3,j2) is even,
    K(i3,j3) is even,                                    (3)

where all unlisted children occur once. This stronger second-child clause is essential; it is not inferred from the one-repeated-child case alone.

## 3. Complete modulo-eight formula

Write K0=K(1,...,1), K_i=K(i2), L_i=K(i3), and K_ij=K(i2,j2) for i<j. A coordinate of index at least four vanishes modulo8. Every repeated child contributes an even factor, so any word with three or more repeated children vanishes modulo8 as well.

The remaining patterns are: no repeated child, one child visited twice or three times, or two children each visited twice or three times. Equations(3) make every two-child pattern involving three visits vanish modulo8, because its profile product is divisible by four. The four-times-b_i part of the one-child three-visit term also vanishes, since L_i is even. Therefore the exact formula is

    Phi_Q = K0
      +sum_i(2+4b_i-2s_i) K_i
      +2sum_i c_i L_i
      +4sum_(i<j)(1-s_i)(1-s_j) K_ij            (mod8).  (4)

This proves that all five-type dependence is affine in b and c, with only singleton-singleton quadratic terms. There are no c_i c_j, c_i s_j, b_i c_j, or b_i s_j terms. Formula(4) holds separately for each full response, not just for a quotient-converse contrast.

When every child is nontrivial, s_i=0. Let F_P be the full response on Q[TT2,...,TT2]. Then

    Phi_Q = F_P+4sum_i b_i K_i+2sum_i c_i L_i (mod8),    (5)

and each L_i is even. Thus the complete coarse bank is AFFINE in the2m nontrivial type bits. In particular1+2m coefficient evaluations determine all4^m normalized profile banks. This is an exact structural reduction, not a proof that any contrast is nonzero.

## 4. Extraction from actual small model graphs

For the fixed quotient Q and fixed block masks, let F_0 denote the actual response when every child is a singleton. Let F_i^(2) and F_i^(3) replace only child i by TT2 or TT3. Let F_ij^(2,2) replace children i,j by TT2. The elementary actual profiles TT1=(1), TT2=(1,2), TT3=(1,6,6) give exact integer identities

    K0=F_0,
    K_i=(F_i^(2)-F_0)/2,
    L_i=(F_i^(3)-3F_i^(2)+2F_0)/6,
    K_ij=(F_ij^(2,2)-F_i^(2)-F_j^(2)+F_0)/4.           (6)

The numerator for L_i is divisible by twelve, by(3). All these are complete physical Hamilton counts, including every port position and every physical response contribution.

Thus1+2m+binom(m,2) actual SMALL-MODEL counts per block query determine the entire five-type function(4). Each augmented model has at most m+4 vertices, regardless of the original child sizes. The nontrivial-only formula can instead be determined from an all-pair baseline and2m type derivatives.

These are offline evaluations on known model tournaments, not a claim about the number of queries allowed to an oracle for an unknown original tournament. The scalar product of the original H(A_i) must be restored to recover its response residue. In a quotient-converse comparison with corresponding profiles retained, that scalar is common and odd, so zero versus nonzero contrasts is unaffected.

## 5. Scope and uses

The quotient order, orientation and child orders are unrestricted. The ports must be uniform on each whole child. A fixed ambient can be incorporated into the quotient with its exterior vertices as singleton blocks, but this formula alone is NOT an ambient non-annihilation theorem. It does not assert full named reconstruction or discard the possibility of a zero complete coarse contrast.

The complete named-response and marked-slot frontiers are separate: nonuniform individual-name attachments inside a child generally require marked cover data and are not captured by these three profile coordinates.

The verifier checks the exact extraction divisibilities, the return involutions, full actual canonical C3 banks, and additional arbitrary quotients and actual larger-child controls. Its C3 example verifies that the six nontrivial type bits have rank six in the complete normalized quotient-converse residue bank, consistent with the independently audited unique blind type(E,E,E). This finite example is not an all-quotient rank claim.
