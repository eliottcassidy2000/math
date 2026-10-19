# Directed defect screen under maximum majority span

Date: 2026-10-19. Status: proved sufficient exclusion screen; not a complete c=4 theorem.

Use the setup and directed-defect exchange lemma of one_defect_exchange_reductions_20261019.md. Assume the poset has no balanced pair and x,y is the selected maximum-index-span incomparable pair in the strict >2/3 majority order M, oriented x before y. This majority order is total and transitive: a directed triangle would have three comparison probabilities greater than 2/3, although the three cyclic comparison indicators sum to at most 2 in every linear extension. It extends the poset.

Every common neighbor z∈I(x)∩I(y) lies strictly between x and y in M. Otherwise z before x would make the incomparable pair z,y longer than x,y, and z after y would make x,z longer. Hence

    P(x<z)>2/3 and P(z<y)>2/3.

For each common z, let s(z)=7−d_Q(z). Its exterior incomparability counts q_A,q_B satisfy q_A+q_B≤s(z). The exact ambient directed defects are

    δ_P(x,z)=δ_Q(x,z)+q_B,
    δ_P(z,y)=δ_Q(z,y)+q_A.

Consequently either of the core-checkable inequalities

    δ_Q(x,z)+s(z)≤1,
    δ_Q(z,y)+s(z)≤1

contradicts the required majority orientation, by the one-defect exchange lemma. This is a sufficient core exclusion even if neither root/common pair is independently known to be balanced. No realizability assumption about arbitrary boundary weights is involved.

## Resolving the two-predecessor slack-one exceptional side

Suppose a common z has exactly two core comparable vertices, both predecessors b,w with w common. Then δ_Q(z,y)=0 and s(z)=1. Therefore the second directed screen excludes this core for every actual exterior. In the previously unresolved situation where the sole exterior incomparable element t lies in B, q_A=0 and δ_P(z,y)=0, giving the stronger P(z<y)≤1/2. In particular there is no need to migrate the witness to t in this case.

The other conditional slack-one case b<z<w with a missing lower-exterior relation is not settled by this argument in general. The note makes no claim that all such cores pass this screen.
