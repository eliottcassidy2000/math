# Exact restart costs for both branches of the stored doubling gate

2026-11-28. This follows the [independent +2 / dyadic switch](../chanpak_unit_gap_dyadic_gate_20261128). After storing the high switch, its two branches have h=303 (or 305 after the old unit flip), e=g=9, and u/g=55 or 63. Keeping the old root, we can rearm either branch by attaching a marked nontransitive suffix, but not by the same suffix core. The **minimum possible suffix Hamilton count is 9 on the lower branch and 45 on the upper branch**. The lower result has a short structural proof, and the upper minimum uses a complete finite catalogue of the relevant strong components.

## 1. The exact contact subset sum

Use a nonempty tournament B, with old nonroots beating all of B and arbitrary contacts between the retained root a and B. Write b=H(B); let C be B plus a. The marked-block law from the previous note gives
\[
g'=gb,\qquad e'=g e_C,\qquad u'=(u-g)Q(B)+g u_C,
\quad Q(B)=\sum_{S\subseteq B} H(B[S])H(B-S).
\]
The directed contact a→x is a binary choice that contributes a nonnegative weight to \(u_C-e_C\):
\[
w_x = \#\{\text{full paths in }B\text{ ending at }x\}
+\sum_{\substack{S\subseteq B\\x\in S}}
 H(B-S)\,\#\{\text{paths in }B[S]\text{ starting at }x\}.
\]
The singleton root contributes b to \(u_C\); other root-starting parts are classified by their first B site. The count \(e_C\) sums the full B paths whose last site beats a, so the two b terms cancel in the difference. Hence \(u_C-e_C=\sum_{a\to x}w_x\). Summing over x gives \(Q(B)\): the first part contributes b and the second contributes \(Q(B)-b\). For \(U=u/g\), the exact new guard is
\[
\frac{u'-e'}{g'}=\frac{(U-1)Q(B)+\sum_{a\to x}w_x}{b}.
\]
A power of two on the right rearms the gate, and each solution of the subset sum gives literal arc contacts. The old flat differential \((2,0,0,0)\) becomes \((2b,0,0,0)\): the old arc stays independent but may cease to have unit amplitude.

## 2. Two minimum-cost witnesses

Masks encode forward bits for pairs i<j in lexicographic order, least significant bit first. Contact bits are the sites of B beaten by a. Start with the h=303 member of each stored branch.

| Branch U | order, forward mask | contact mask | b, Q(B) | new (h,e,u,g) | new guard |
|---:|---|---:|---|---|---:|
| 55 | 5, 10 | 29 | 9, 84 | (13113, 27, 41499, 81) | 512 |
| 63 | 8, 238981866 | 127 | 45, 1464 | (241569, 135, 829575, 405) | 2048 |

The corresponding flat companions have exactly the same side registers and h larger by 18 and 90, respectively. Each witness may now invoke the original storage gate again within its branch.

For minimality, the weaker necessary interval is already enough:
\[
[(U-1)Q/b,\ UQ/b]\quad\text{must contain a power of two}.
\]
Both b and Q multiply in an ordinal join, by factorization of paths in every induced subgraph. Singleton components only multiply Q by 2 and do not affect whether such an interval meets the powers of two. For \(b<9\) the nonsingleton core is empty, a directed triangle, or the strong four-vertex tournament; their respective \(Q/b\) are \(1,4,32/5\). This follows from the strong-order bound below, the unique small types and oddness; none of the intervals for U=55 hits a power of two. The lower witness attains 9.

For the upper branch it suffices to exclude \(b\le43\), since b is odd. Use the [strong tournament bound \(H\ge 1+2^{n-2}\)](../nemenyi_stress_sparse_path_range_20261127/README.md) on strong components. Each component of a B with b≤43 then has at most seven vertices. The replay generates the full isomorphism catalogue through order seven by adjoining a vertex to each representative with every possible contact orientation. Its exact canonical code permutes within equal-outdegree blocks, which is sufficient because any isomorphism preserves the outdegree partition. The successive class counts are 1,1,2,4,12,56,456. Count b and Q directly on the strong representatives and close their profiles under multiplication for b≤43, dropping all singleton factors. No resulting \([62Q/b,63Q/b]\) contains a power of two. The order-eight witness above has b=45 and attains the minimum. The finite enumeration is over bounded strong components, but the exclusion holds for suffixes of **every order**, including unlimited singleton padding and any contacts.

## 3. No common restart core, even with branch-dependent padding

Even if contacts differ and the target powers of two differ, a common B would require one power in \([54,55](Q/b)\) and another in \([62,63](Q/b)\). The ratio of any number in the second interval to one in the first is strictly between 1 and 2. This cannot be a ratio of powers of two. Moreover, if the two suffixes have the same nonsingleton strong components, their Q/b differ only by a power of two. Dividing each proposed guard by its own dyadic singleton factor gives the same obstruction. Thus branching contacts or arbitrary transitive padding between a common sequence of nontrivial components cannot renew both branches. Several successive fixed-root suffix blocks reduce to this same calculation by taking their ordinal join.

This is a limitation of the fixed-root precedence interface, not of arbitrary tournaments or of root-changing gadgets. It also makes the cost distinction operational: a single-branch restart is feasible, but it cannot preserve the old +2, and renewing both stored branches requires genuinely different cores or a different interface.

## Replay

Run node replay.js. It verifies the actual suffix registers and contact weights independently by subset DP, produces the complete strong-core profile catalogue needed for the upper minimum, and checks the interval exclusions after strong-component factorization.
