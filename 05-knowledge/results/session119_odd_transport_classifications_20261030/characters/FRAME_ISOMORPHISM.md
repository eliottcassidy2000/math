# All prime-power valuation frames form one tournament isomorphism class

**PROVED all-parameter extension; no novelty or full-automorphism claim.**
The p=3 all-positive case is inherited explicitly from
[THM-4754, the 3-adic Sierpinski tournament](../../../../01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md).
[THM-4755](../../../../01-canon/theorems/THM-4755-blow-up-gains-skew-spectrum-lexicographic-law-ternary-tower.md)
already gives the general lexicographic spectral law. The present statement
identifies every valuation frame from the character classification with one
lexicographic family for every prime p≡3 mod4. Ring residues Z/p^a and the field
F_(p^a) are different objects; this is not a claim about field Paley tournaments.

## 1. Digitwise isomorphism for every frame

Fix p≡3 mod4 and a>=1. For frame ε=(ε_0,...,ε_(a-1)) in {±1}^a orient Z/p^a by

    x→y iff ε_v χ_p((y-x)/p^v)=+1,  v=v_p(y-x)<a.

Write x=Σ_{j=0}^{a-1} x_j p^j with digits x_j in F_p, least significant first.
For x≠y, let v be the first differing digit. The unit part of y-x modulo p is
exactly y_v-x_v; any borrow or carry affects higher digits only. Consequently
arcs are determined by the first differing low digit, with the Paley orientation
or its reverse prescribed by ε_v. This is an iterated lexicographic product of
those p-vertex tournaments, with low digits labelling the outer blocks.

Given two frames ε and ε', choose t_j in F_p^× with

    χ_p(t_j)=ε_j/ε'_j.

Such t_j exist because both square classes are nonempty. Define the bijection

    F(x)=Σ_j [t_j x_j mod p] p^j.

Its first differing digit on a pair is still v, and its difference there is
 t_v(y_v-x_v). Therefore ε'_v χ_p(t_v(y_v-x_v))=ε_v χ_p(y_v-x_v).
F preserves every arc. This proves:

**Theorem.** All 2^a labelled valuation-frame symbols produce a single tournament
isomorphism class, namely the a-fold lexicographic power P_p^[a] of the p-vertex
Paley tournament. The labelling count and the isomorphism count are different.

The digitwise map is generally not a single multiplication modulo p^a. It does
not assert equivalence under the narrower affine multiplier action. For example,
a global unit multiplier changes every frame sign by the same χ_p(u), whereas
independent digit maps can change each layer separately.

## 2. An explicit iterated subgroup and a divisibility consequence

On P_p, every map t→ut+b with χ_p(u)=+1 preserves arcs. These maps form a group
of order g=p(p-1)/2. At every internal node of the rooted p-ary digit tree, apply
such a permutation to its next digit, leaving all later digits unchanged; the
choice may depend on the preceding low-digit prefix. The first differing digit
argument again proves every resulting tree permutation is a tournament automorphism.

There are

    E=1+p+...+p^(a-1)=(p^a-1)/(p-1)

internal nodes. The independent local choices give a faithful iterated wreath
subgroup of order g^E. One can verify the count recursively: at the root choose
g local permutations and independently choose an automorphism of each of the p
subtrees, giving G_a=g G_(a-1)^p with G_1=g. Choices are unique because the
permutation on first digits and then on each rooted subtree can be recovered
from the resulting vertex permutation. Conjugate by the digitwise frame map to
obtain the same subgroup for every frame.

Let H(T) count ordered directed Hamiltonian paths. Any automorphism fixing one
such ordered path fixes every vertex in its list, hence is the identity. The
displayed subgroup therefore acts freely on the Hamiltonian paths, proving

    [p(p-1)/2]^((p^a-1)/(p-1)) divides H(T).

This is a divisibility statement, not an exact path-count formula. It does not
identify the full automorphism group. For p=3 the known ternary tower has the
same displayed group order; that special case and earlier computed path counts
are credited to THM-4754/4755 rather than rediscovered.

## 3. Literal checks

`check_frames.py` checks every map between every frame pair for the listed prime
powers, including 3^4,7^2,11^2, against every ordered pair of vertices. It separately
checks the tree-node translation/square-multiplier generators on every prefix.
For orders at most11 it computes Hamiltonian paths by an independent subset DP
and verifies the displayed group-order divisibility. Its receipt distinguishes
these finite controls from the proof and from an unproved full-Aut assertion.
Normal and optimized output agree with `frames_receipt.json`.
