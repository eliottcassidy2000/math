# Independent alternating-hub positivity and root-family review

Status: PASS. No mathematical correction required. This review covers Sections 4–6 of cyclic_child_disjoint/PROOF.md; the rank-one kernel and its marked-profile restoration in Sections 1–3 receive a separate independent review. Final producer proof 2c625e1a7c7cff166813d7f7243220b160fe04c3f32418e75bb07f267d9fb1be and verifier 1d727aac578e9aef9377136b7f6cd1d8899ed259fedb379e4d19a9a479cd7133 are pinned.

## Exact color-algebra and positivity review

The transitive chain run matrix is correct: a run word uses each color at most once, in increasing order. Marked first/last physical endpoint profiles replace ordinary q profiles at the two extreme children. Their prescribed-slot convention introduces no factor equal to the number of run occurrences. Full profiles may be arbitrary; only the two marked profiles must be balanced and have positive first coordinates.

For even u, the identity S_u*g=G_(u,last)*s gives rho_last*s=(beta-1)g. For odd v, E_v*g=e*G_(0,v) gives e*ell_0=(alpha-1)g. Consequently 2alpha beta g minus those two products is exactly (alpha+beta)g. The other terms yield precisely the stated expression (6), with all integer signs retained.

Project the final summand onto monomials containing every chain color. At k>=4, the last middle color is outgoing and the first middle color is incoming. Only E_last, S_0 and G_(0,last) can cover all intermediate colors in that summand. The two positive copies of b_0(2)b_last(2) cancel its negative two copies. The remaining full-support polynomial is exactly

    (z_0 z_last^2+z_0^2 z_last) product_(0<t<last) z_t.

At k=2, add z_0 z_last. This is stronger than merely checking a numerical sign: after discarding incomplete-support monomials, the final polynomial is coefficientwise nonnegative. The first two summands are already coefficientwise nonnegative. The all-once monomial in (alpha+beta)g has coefficient at least two, proving the claimed strict lower bound for every admissible profile.

No profile coordinate above three can occur here: the three exterior path slots each meet a chain child at most once. The verifier's h<=3 truncation therefore retains the complete disjoint quartic even for five-name extreme children.

## Family and physical-name review

The W quotient is prime, including the C3 boundary, and C3 children are prime. Hence its proper nontrivial original modules are exactly its children. If the hub child is selected along with the root, its orientation is restored; reversals of singly marked balanced extreme children and unmarked children are invisible to the selected quartic. If the hub is not selected, queries uniform on all other children agree exactly with whole converse. Its three physical names remain singleton blocks, so the exposed-anchor theorem applies.

For every orientation state, either cyclic order of the hub is admissible. Whole quotient converse is returned to the same W convention by reversing chain indices; this transports complete physical child sets and the resulting query masks, rather than identifying names. Every triangle remains endpoint-balanced. Thus the argument works against every source state, giving the stated orientation-cube distinction. Families with no root are the separately established disjoint-module case.

With m quotient children there are m+2 compressed blocks in the hub-unselected branch. The inherited bound gives at most m+4 anchor candidates, so adjoining sixteen quartic queries gives the safe m+20 detection bound. This is not asserted as one universal reconstruction battery. The all-order result is confined to the displayed core or inherited ordered SCC padding; arbitrary strong outer lifting is correctly excluded.

## Independent executable evidence

The audit imports no producer code. It expands the complete color polynomials at k=2,4,6,8 and verifies the original kernel reduction coefficient by coefficient, including the exact full-support last polynomial and strict all-once term. It explicitly enumerates marked child covers, reconstructs all six W literal examples and all forty retained root-containing families, checks every original module and XOR target, and rebuilds the compressed anchor masks.

A separate predecessor-DP C++ engine recomputes every retained physical query in these sections. It receives independently reconstructed graph rows, with all actual query values and cubic/quartic effects intact. Both independent modes are byte-identical and the producer's final paired receipts agree. Finite checks validate the all-order algebra and branch bookkeeping; they do not replace it. The kernel audit in Sections 1–3 is deliberately not subsumed by this report.
