# Exact block rescaling, and why it cannot preserve a unit gap

For the rooted state (h,g,R,S), let I=h−2S+2g and Z=R+S−g.
Appending a whole nonempty tournament B as a near-sink block
means all old nonroot vertices beat it, and it beats the old root.
Put h_B=H(B) and let q_B count ordered two-slot path covers of B
with an empty slot allowed. Then its literal output and new state are

    H_out = h_B I + q_B Z
    (I,Z,g) maps to (h_B I, q_B Z, h_B g).

Proof: let U,V sum H(B[A])H(B−A) over nonempty A of odd, even
size respectively. Then U−V=h_B: in the alternating sum over all
A, expand both H factors by disjoint odd-cycle packs. Each odd
component, including any singleton, can be assigned independently
to either side, cancelling by sign. Also h_B+U+V=q_B.
Root paths can traverse an old nonsingleton root path and then
a final B run; grouping by parity gives

    R' = h_B R + V R + U(S−g),
    S' = h_B S + U R + V(S−g),
    g' = h_B g.

The h count includes the ordinal product h_B h plus
2UR+2V(S−g) from odd cycles through both old and new vertices.
Such a mixed cycle must contain the root, and no two such cycles
can occur in one packing. Substitute to obtain the formulas
and note S'−R'=g'. Thus this operation resynchronizes root
parities even if they did not start synchronized.

For a fully ordinary sink block (all old vertices beat B)
the same argument gives h'=h_B h, g'=h_B g,
R'=S'=q_B(R+S)/2, and

    Z' = q_B Z + (q_B−h_B) g.

In either common block operation a difference of exactly two
between synchronized readers is multiplied by h_B. To preserve
that gap with a common block, one must have h_B=1, hence B
transitive and q_B a power of two. At fixed h_B there are only
finitely many normalized odd parts of q_B: the nonsingleton
strong SCCs have bounded orders and number by ordinal
factorization and Busch's order bound; singleton SCCs multiply
the two-slot count by two. So replacing each binary letter by
one large uniform block does not by itself give a variable
odd radix while retaining +2. A useful next move would have to
use different blocks on different readers, or compensate a
scaled gap with a separate physical port.

This barrier is for the specified uniform append architecture,
not arbitrary interactions or arbitrary tournament growth.
