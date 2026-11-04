# Independent small-physical replay recipe

Let E_j mean evaluation with one quotient child replaced by actual TT_j, holding all other profiles fixed. The coefficient functionals selecting that child's run multiplicities are

    L1=E1,
    L2=(E2-E1)/2,
    L3=(E3-3E2+2E1)/6,
    L4=(E4-6E3+11E2-6E1)/24,
    L5=(E5-10E4+35E3-50E2+24E1)/120.

These are exact integer identities obtained from q_k(TT_j)=k!S(j,k). Tensor products act independently on different children. They select precisely the corresponding run-word coefficient; no rounding of a division is permitted.

To reproduce the TT4 baseline modulo16, it is enough to select all multiplicity vectors with total sum(k_i-1)<=3, each k_i>=1. There are56 vectors on five quotient names. Every other product of TT4 run weights is divisible by16. For each selected vector, apply the tensor functional above and multiply by product_i q_(k_i)(TT4), then sum modulo16. Every physical evaluation required has at most five+three old names and two ports, hence at most ten names. This supplies an implementation structurally different from the producer's run-word dynamic program.

The first-order q2 corrections use L2 at a single child with all others L1. The q3 and q5 corrections can independently be checked using L3 and L5 with all others L1, requiring at most11 physical names. Their parity is zero for q3/q5 by the proved return-segment involution.

This recipe applies to complete responses at each actual block mask. It does not require sampling, a formal-profile realizability assumption, or a large22-name reconstruction.
