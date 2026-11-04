# S187 integration: exact cards, short supports, and stronger gate clocks

**PROVED inherited/scoped; FINITE-EXACT cross-checks.** Reviewed incoming
commits8103f2a359 through1b2caee0ce, plus compensation commits7228c8e534
and7018acb714. Universal prime12 exclusion and spectrum completeness remain
**OPEN**. This is an integration review, not a fresh independent review of
all incoming packets or a replay of their large finite searches.

## The exact current and the complete modular cubic agree

S187's [prime12 support theorem](../session187_named_response_20261104/prime12/PROOF.md),
with its [independent audit](../session187_named_response_20261104/prime12_audit/REPORT.md),
retains every named cubic of delta=(F_T-F_(T^op))/2 modulo2. At even order n,
only two disjoint nonempty paths with odd union size at most n-3 contribute.
At twelve the union sizes are3,5,7,9. These are support restrictions, not
permission to discard coordinates or replace integer equality by parity.

Combining this inherited theorem with [D1](PROOF.md) gives, for any
odd-order tournament Q of order m>=3 and any u,

    g_u(Q)/2 = sum S_u(Q[A]) (mod2),

where the sum is over A containing u with even size2<=|A|<=m-3.
To prove it, add any new name v to Q, then apply the cutoff to
[a_u a_v c_v]delta. Disjointness forces the middle path through v to be
the singleton {v}. Thus the first path has even size at most n-4=m-3.
D1 identifies the same coefficient with -g_u(Q)/2; the sign disappears
modulo2. The current is even by its exact response interpretation and
the oddness of every tournament Hamilton count. No attachment to v occurs
in the resulting formula.

For the eleven-vertex deletion cards of a twelve-core this requires only
even supports of sizes2,4,6,8. This is a coefficient-level compatibility
with the shared-label filter already used in this packet's SAT scout,
not a new stronger SAT condition or a solver conclusion.

`incoming.py` independently accumulates bit-packed path endpoints over
disjoint supports and compares **all7920 modular cubic coordinates** of
our four saved integer response contrasts and S187's prime control P.
It also compares528 exact deletion currents with the short-even-support
formula. No polynomial is reinterpolated from only those shared terms.

The incoming two radius-four balls contain1537384 source/edit
presentations,1407 endpoint-balanced presentations,1089 prime among those,
and zero combined endpoint/cover-skew survivors. This is credited bounded
evidence, not a census of all prime12 cores. One centre is our retained
H300453 core. The new prime P has H84599 and is already distinguished by
a singleton query; it does not prove that cubic coordinates are necessary
for its first detection.

## The seed interval has a fourth named-response exponential

The [S187 gate proof](../session187_named_response_20261104/gate/PROOF.md)
uses the inherited singleton score R_v=E_v+sum_w C_vw and our earlier
positive adjacent-twin contraction. For the kth module vertex,

    R_(v_k)(T_j(r))=18*3^(r-k)*4^(k-1)
                    +B_j*2^(r-k)*3^(k-1)
                    +C_j*2^(k-1)+D_j*1[k=r],

where (B_j,C_j,D_j) are(44,-5,-4),(42,-1,-3),(42,3,-3).
The actual query is H_j(r)+2R_(v_k). This retains a genuine4^(k-1)
channel absent from the scalar three-point interval's3^r/2^r law.
All27 rank queries and six exterior converse separators in our earlier
complete named tables at r=2,7 match these incoming formulas exactly.

The incoming theorem reconstructs the entire specified6r! family with
r+2 queries. Its [color-profile proof](../session187_named_response_20261104/collision/PROOF.md)
also separates any two transitive orders in every fixed quotient, even
when the compared tournaments have no common homogeneous pair. This
extends the closed transitive branch, without resolving prime cores.

For reversal only of5->6, the full response retains the quartic

    [a_0 a_3 c_0 c_6](F_T-F_T')=2^(r-2)+2-1[r=2]>0.

This all-order incoming result strengthens the previous finite hostile.
We retain the whole-converse cubic and the internal-reversal quartic in
their separate valid scopes. The original full-response tables are not
overwritten by the clock formulas.

## Exact SCC localization and the preparation-cost observer

The [S187 SCC theorem](../session187_named_response_20261104/converse/PROOF.md)
localizes an ambient hidden module converse exactly to the quotient SCC
containing its root. Every localized response is multiplied by the positive
product of H over the other SCCs. If the root is a singleton SCC this
implies exact intrinsic equality, so D1 forces exact deletion-current zeros.
For a nontrivial strong active quotient, it does not cancel that quotient;
our universal ambient mod4 scope is unchanged.

The [incoming preparation-cost identity](../compensation_returns_20261103/QUERY_COST.md)
uses the same physical singleton observer:

    H(T+z)-H(T)=(F_(T^op)(e_w,0)-H(T))/2-end_w(T).

The two stated children1409 and1427 have charges3632 and3536, producing
5041 and4963. This measures the missing preparation cost; it does not pay
it or give a count-preserving reset. The [cap-width theorem](../compensation_returns_20261103/CAP_WIDTH.md)
also shows endpoint mass1 alongside feedback number ceil(j/2), so a bounded
endpoint statistic cannot justify a fixed-width state model at arbitrary
depth. These incoming structural obstructions preserve the two infinite
prime/7p realization obligations in [BRIDGE.md](BRIDGE.md).
