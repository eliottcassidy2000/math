# Twenty-four remaining center coefficients

**PROVED reduction conditional only on the explicitly inherited low-order data.** This is a support constraint, not a determination of the unknown coefficients.

Combine the independently reviewed degree47 theorem and its score-run restrictions with the two-shore insertion exclusion in this packet. For each possible center, choose any marked Hamilton cycle and write n=qm with exact cut index q and odd m.

- q=1 permits odd n<=47.
- 2<=q<=8 requires m in{1,3,5,7}.
- q=2,m>=5 is excluded by Delta>=2m-2.
- q=6 or8,m>=5 is excluded by the score-run theorem.
- q=5 or7,m>=5 requires regularity; regular n>=49 is excluded.
- Order2 cannot be strong.

The possible EVEN orders remaining from these necessary constraints are exactly

    {4,6,8,12,18,20,24,28}.

“Exactly” describes this parameter-filter output, not existence of centers at each order. The set of potential odd orders continues through47; no odd order in11..47 is removed merely by these constraints, since the q=1 lane is still permitted.

For the small coefficients we explicitly inherit the complete low-order data reported in `deficiency_three_20261020/centers_check.out` (and the preceding center packets): R6 has coefficients1 at orders5,6,7 and zero at the other orders through9. The degree47 packet freshly authenticated Paley7 but did not rerun this whole low-order census. The present support reduction likewise does NOT claim a fresh verification of those inherited coefficients, particularly order9.

Consequently

    R6(z) = z^5+z^6+z^7 + sum_{n in U} r_n z^n,

where every r_n is an unknown nonnegative integer and

    U = {11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47}
        union {12,18,20,24,28}.

There are24 unresolved coefficients: nineteen odd and five even. Some or all may vanish; the reduction gives no positive existence claim, uniqueness claim, or upper bound for an individual r_n. It remains consistent with R6=z^5+z^6+z^7, but does not prove that equality.

The small `support_check.py` independently enumerates the parameter filters, checks the even set and cardinality, and records its own source hash and actual clocks. It is arithmetic bookkeeping, not tournament enumeration.
