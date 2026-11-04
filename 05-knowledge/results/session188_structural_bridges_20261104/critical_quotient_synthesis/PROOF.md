# Every critical quotient with a pair child is separated

## Verified classification input

A tournament of order at least five is critically indecomposable if it has no proper nontrivial module but deleting any vertex creates one. Schmerl and Trotter classify these tournaments into three families in Corollary 5.8(3), printed page 203; their definitions are on page 198. The primary author-hosted PDF was inspected directly, including both rendered pages: [Critically indecomposable partially ordered sets, graphs, tournaments and other binary relational structures](https://trotter.math.gatech.edu/papers/83.pdf), Discrete Mathematics 113 (1993), 191–205, [DOI](https://doi.org/10.1016/0012-365X(93)90516-V).

Their family 1 is our cyclic quotient. Their family 2 maps to our staircase quotient by a_i -> B_(r-i) and b_i -> A_(r-i). Their family 3 maps to our alternating-hub quotient by a_i -> i-1 and b -> hub. Thus their classification covers exactly the three explicitly defined quotient families used here, at every order 2r+1>=5. The checker verifies these label conventions independently at finite orders; the formulas establish them generally.

## Consequence of the three response theorems

Let Q be any critically indecomposable tournament with at least five vertices. For every quotient name v, let D_v be an arbitrary tournament with at least two physical vertices. Suppose at least one D_z is exactly a directed pair p->q. Then

    F_(Q[D_v]) != F_(Q^op[D_v])

as complete named two-port responses, with all child interiors and all physical names unchanged.

Indeed choose the pair child at z. The classification identifies Q with one of the cyclic, staircase or alternating-hub quotients. The corresponding theorem supplies disjoint unions X,Y of exterior child blocks such that X->Y, M=X union Y is a module after deleting D_z, and |M|>=3. The exact grouped quartic is

    [a_q a_X c_p c_Y](F_(Q[D_v])-F_(Q^op[D_v]))
      = C_YX(Q[D_v]-D_z)-C_XY(Q[D_v]-D_z) > 0.

The sixteen associated actual queries therefore contain a separator. Relabeling the classification's quotient names transports entire physical child sets; it does not identify child names or port bits. The only internal nonuniform port attachment is at the designated pair. The exact finite difference retains every cubic and quartic contribution.

Dependencies are ../ordered_module_skew/PROOF.md (cyclic family and inequality), ../alternating_prime_quotient/PROOF.md and ../staircase_prime_quotient/PROOF.md, with the exact coefficient in ../marked_pair_quotient/PROOF.md.

## Boundaries

This is a classification-based corollary for critical prime quotients, not for every prime quotient. Arbitrary other child orders and orientations are allowed subject to nontriviality. It does not cover the case with every child of order at least three and no pair child, does not claim arbitrary outer-ambient lifting, and does not itself prove injectivity for all module-reversal families.

The directed triangle is handled separately by the elementary ordered-cut case: with a pair root, exterior total order at least three suffices. A zero quartic when both exterior children are singletons is not a full-response collision.

Publication excludes the downloaded source PDF and rendered pages; cite the verified primary links instead.
