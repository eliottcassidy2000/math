# A dual-query separator robust under increasing child profiles

**PROVED:** a block-uniform quotient-converse separator can be made strict in the direction that survives every coordinatewise increase of the target child cover profiles. This yields a fixed two-query battery for all root-containing original-module XOR families whenever the quotient has a nonzero Hamilton endpoint-parity bit. One of the two queries, chosen from the baseline source and quotient converse, works simultaneously for every such family, also in every fixed ambient tournament.

## 1. Uniform-port monotonicity

Fix disjoint nonempty physical child sets A_i and a quotient tournament Q. Write T=Q[A_i] and U_0=Q^op[A_i], retaining all source child interiors in U_0. The target U=Q^op[B_i] uses the same child sets and orders, and suppose

    q_k(B_i) >= q_k(A_i) for every child i and every k>=1.       (1)

Here q counts ordered covers by nonempty directed paths. Attach actual ports x->y, with a_v=[v->x] and c_v=[v->y], uniformly on each child. Each child remains a module of the actual augmented graph. Its Hamilton count is a polynomial in the complete child q profiles with nonnegative integer coefficients, obtained by fixing the word of maximal child runs and the two singleton port positions. Equivalently replace one child at a time using the nonnegative module-run substitution formula. Hence

    F_U(a,c) >= F_(U_0)(a,c).                                  (2)

Therefore any query with F_T(a,c)<F_(U_0)(a,c) also separates T from every U satisfying (1), with the same strict direction and at least the baseline gap. No equality of the target profiles is required.

The same argument holds in a fixed rooted ambient tournament R, replacing its root by T,U_0,U. Both ports must receive the root's exterior incidences. They remain uniform on every child, so (2) continues to hold for the complete physical augmented ambient graphs.

## 2. The core converse-complement involution

For a block-uniform core query z=(a,c), let

    Jz=(1_M-c, 1_M-a),

where M is the core name set and subtraction means Boolean complement inside M. Converse every arc of the augmented core and interchange the names of its two ports. This restores x->y and transforms z into Jz. Each child interior is conversed, but it remains a module and has exactly its original complete q profile. Restoring child interiors changes no ordered-cover count. Thus, for every h>=1,

    q_h(T+ at z) = q_h(U_0+ at Jz),
    q_h(U_0+ at z) = q_h(T+ at Jz).                            (3)

In particular the entire augmented-profile contrast at Jz is the negative of the contrast at z. This is a statement about all q_h, not merely Hamilton counts or coefficients through degree three.

In a fixed ambient, lift both queries by giving BOTH ports the original root incidences to the exterior. Complement only the core bits. The exterior bits are not complemented. The two augmented cores are modules with the same collapsed ambient R. Applying the same nonnegative run-substitution coefficients to (3) proves

    Delta_R(Jz) = -Delta_R(z),                                (4)

where Delta_R is the baseline source-minus-U_0 physical response. The ambient tournament itself is never conversed.

If either baseline contrast is nonzero, exactly one is negative. Choose that member once. By Section 1 it separates the source from every profile-increasing target (1). The choice may depend on the source interiors and ambient, but never on which permitted target changes occur. Without making that choice, the fixed two-query set {z,Jz} is a separating battery for all such targets simultaneously.

## 3. An explicit all-order sufficient hypothesis

Let S_i(Q),E_i(Q) count quotient Hamilton paths starting and ending at i. Suppose

    S_i(Q)+E_i(Q) is odd.                                     (5)

Use the two core queries

    z:  a=1_(A_i), c=0;
    Jz: a=1_M,     c=1_(M minus A_i).                         (6)

The independently audited endpoint-profile transfer theorem in `../parity_profile_transfer/PROOF.md` proves that the first baseline contrast is 2 modulo 4, for arbitrary nonempty children, including singletons. Its fixed-ambient lift retains this residue. Hence both contrasts in (4) are nonzero, and Sections 1–2 apply.

**Theorem.** Under (5), the fixed pair (6), with the specified exterior lift, separates Q[A_i] from every Q^op[B_i] satisfying (1), in every fixed ambient tournament. A single member chosen to have negative baseline contrast works for all these targets at once. The final contrast need not retain the baseline modular residue when child profiles increase; the guarantee is a strict integer sign.

This is a two-actual-query upper bound for the given source/partition/anchor and family, not a universal reconstruction bound or a claim that a predetermined one of the two masks always has the correct sign. Selecting the single member from known baseline counts is an additional computation, not a free oracle query to an unknown tournament.

## 4. Root-containing original-module families

Assume Q is prime and satisfies the endpoint hypothesis. A directed triangle has zero endpoint bits and is not an instance of this sufficient criterion. Every proper original module of Q[A_i] lies inside one child. Select any family of original modules whose whole root is reversed oddly, and apply all edge reversals simultaneously by XOR. The final quotient is Q^op. Within a child, the final tournament is the converse of the result of its selected proper-module reversals. Converse does not change q. The hereditary original-module monotonicity theorem in `../module_decomposition/PROOF.md` gives exactly (1), including crossing or nested selected modules inside a child.

Consequently (6) separates EVERY such root-containing family, including the families whose original Hamilton count increases. There is no scalar-query fallback in this strengthened result. The older scalar-query-plus-endpoint battery remains valid; this is a separate stronger corollary, not a modification of its frozen proof.

The same conclusion holds when all selected operations are contained inside a fixed ambient core. It does not permit arbitrary changes to the exterior or a different target quotient. No exact response collision is asserted.

More generally, Sections 1–2 apply whenever one already has a negative baseline block-uniform separator. The equal-profile cyclic SOS query and its own bounded-order ambient extension are complementary applications developed separately; this packet does not duplicate or enlarge their hypotheses.

## 5. Controls and complete named data

The verifier enumerates every original-module family inside several small child systems, including overlapping intervals of a transitive three-vertex child that create a directed triangle. It checks complete child-profile monotonicity, both dual baseline signs, both actual queries on every target, and the common fixed physical exterior lift. A separate literal C++ Hamilton DP retains actual values rather than truncating the response.

For one seven-old-name source and a scalar-strict crossing-family target, it retains the complete named response bank: every bidegree-at-most-(2,2) coordinate, including zeros, its actual source/baseline/target values and all three exact Boolean coefficients. In particular all cubic and quartic data are preserved. Random full-mask queries independently check reconstruction from that complete polynomial. The universal proof is Sections 1–4; the finite controls are not an extrapolation to arbitrary size.
