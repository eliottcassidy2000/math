# A recovered incidence recurrence for native partner path counts

**PROVED elementary inclusion-exclusion below; FINITE-EXACT through five
vertices.** This recovers a classical rook/path contraction mechanism for
the present complement problem; no literature-priority claim is made.

The route is
[THM-4225, bad-owner upper-zeta successor-rook hierarchy](../../../01-canon/theorems/THM-4225-bad-owner-upper-zeta-successor-rook-hierarchy.md).
It counts prescribed cyclic successors by injective, cycle-free rook
placements. Its factorial contraction proof survives for linear Hamilton
words, with one more available block position.

Let B be any target tournament and D=B^op its forbidden directed adjacencies.
For subsets of row vertices U and column vertices V define R_D(U,V;z) to
count partial injective maps from U to V whose chosen arcs lie in D, weighted
by z per chosen arc. For any i in U,

    R_D(U,V;z)=R_D(U\{i},V;z)
       + z sum_(j in V, i->j in D) R_D(U\{i},V\{j};z),       (R1)
    R_D(empty,V;z)=1.

This is an exact row/column recurrence; it keeps joint incidence, not just
row degrees. A rook placement may nevertheless contain closed directed
cycles, which cannot occur inside a linear Hamilton word. Let C range over
collections of pairwise vertex-disjoint directed cycles in D, including the
empty collection, and write V(C) for their vertex union. Then the cycle-free
rook polynomial is

    A_D(z)=sum_C (-1)^|C| z^|V(C)| R_(D-V(C))(z).            (R2)

Proof: in any injective successor placement the cyclic components are
disjoint. Choosing a subcollection of them contributes a sum of signs
(1-1)^c. It vanishes unless there are no cyclic components. Distinct
directed cycles on the same vertex set retain their multiplicities in R2.

Write a_e=[z^e]A_D. Inclusion-exclusion over forbidden consecutive arcs gives

    H(B)=sum_(e=0)^n (-1)^e a_e (n-e)!.                     (R3)

A cycle-free placement is a disjoint union of ordered paths and isolated
vertices. Contracting its e arcs leaves n-e blocks, in (n-e)! linear orders.
All other chosen-arc sets force an impossible successor relation. This
proves R3. For the actual metagraph take B=C_P(T), not T^op.

The minimal hostile is D=C3:

    R_D=1+3z+3z^2+z^3,
    A_D=1+3z+3z^2.

Using R_D directly in R3 gives 2; excluding its closed cycle gives the actual
3 Hamilton paths. rooks.py tests R1--R3 on all 74 fixed-path presentations
of sizes 3,4,5, compares R2 to literal arc-subset enumeration, and compares
R3 to an independent Hamilton-path dynamic program on the native partner.

For general strong-atom composition, R1 suggests a practical remaining
target: a child response should expose the used incoming/outgoing ports and
which path endpoints are connected. The cycle rejection in R2 is a global
realization guard. Scalar rook coefficients are sufficient for the final
unattached H count but are not a proved substitution-closed attachment
carrier. The cheapest next test is to compose two child rook responses
through a revisiting quotient walk and compare with literal global R2.
The inherited six-vertex two-endpoint obstruction must remain active.
