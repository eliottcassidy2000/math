# The alternating-group lift remembers a parity that the merge forgets

This is a direct A5 action at order five, and a related two-sheet result at every order. It does not put A5 inside a tournament automorphism group (those groups have odd order), and it does not assert Bott periodicity.

Label vertices 1,...,n and let chi(T) be the product of all arc signs above the diagonal: + for i→j when i<j and − otherwise. Under a relabeling sigma and global converse, respectively,

    chi(sigma T)=sgn(sigma)chi(T),
    chi(T^op)=(-1)^binom(n,2) chi(T).

To see the first formula, reorder the unordered pairs; the only sign changes come from inversions of sigma. The second reverses all pairs. Let G⁺ be the kernel of (sigma,epsilon)↦sgn(sigma)(−1)^[epsilon·binom(n,2)] on S_n×C2. Every converse-merged class splits into exactly two G⁺ orbits, one at each chi sign. Indeed if a transformation relates two tournaments of the same sign, the formulas force it to be in G⁺; any odd relabeling switches signs. For n≡0,1 mod4, G⁺ is literally A_n×C2. At n=5 the two sheets are the two A5 labelling orbits with converse merged inside them.

For edges, keep labelled path-pointed tournaments (T,P) **only modulo A_n**, not modulo converse: merging converse at the dart level would erase the actual reflected black twins. Quotienting by S_n recovers precisely the original oriented tiling half-edges. Since the action on a directed path-pointed tournament is free, each original dart splits into two lifts, identified by chi(T). Their node projections go to the two G⁺ orbits above. The path-reversal involution on darts reverses n−1 arcs and traverses the reversed path, hence multiplies chi by (−1)^(n−1). It commutes with relabeling and preserves blue/black. Thus for odd n this lift is the disjoint union of two copies of the original metagraph, and for even n it is the canonical bipartite double cover. These statements retain loops and multiplicities as statements about dart involutions. For a color adjacency matrix A (diagonal counted twice), the two lifts have matrices respectively

    [ A 0 ]     [ 0 A ]
    [ 0 A ]  or [ A 0 ].

This also gives a sharp contrast with the nonpath-complement convention. Its sign factor is (−1)^[binom(n,2)−(n−1)], rather than (−1)^(n−1). The two conventions are indistinguishable on the merged graph (THM-4729), but the sign effect upstairs differs exactly at n≡2,3 mod4. The result is a genuine period-four effect of the orientation character, not an eight-period statement about the spectrum or the graph. At n≡2,3, whenever the base graph has an odd closed walk, the two covers cannot be isomorphic even as uncolored graphs: one has such a walk and the other is bipartite. For example, the known order-six loops already force the distinction. No all-order existence of loops or connectivity is used.

The corresponding closed-walk determinant functions are exactly

    det(I−z A_lift)=det(I−zA)^2      (same-sheet move),
    det(I−z A_lift)=det(I−z² A²)    (sheet-switching move).

More generally put any formal weights on the blue and black line colors and replace A by the weighted sum of their adjacency matrices. The block identities remain exact over that polynomial ring. In the switching cover all odd closed-walk coefficients vanish; in the same-sheet cover each coefficient doubles. Thus the orientation-character lift gives a precise periodic rule for derived generating functions even without a scalar recurrence for A itself.
