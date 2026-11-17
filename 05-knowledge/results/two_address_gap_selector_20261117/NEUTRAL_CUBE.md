# Moving the center inside an unchanged gap state

Fix the nonport graph N and the signed mask difference f=I_x−I_a from [the five-coordinate gap state](SELECTION.md). Let p=f^+, m=f^−, and let g be the common-contact indicator on vertices where f is zero. Then I_x=p+g, I_a=m+g. If K=M−M^T, skew symmetry gives

    I_x^T K I_a = p^T K m + f^T K g.

The only gap-state coordinate affected by g is D, and its dependence on g is exactly the affine term f^T K g. All four boundary coordinates use f alone. Thus at fixed f the **exact same continuation state** is a literal 0/1 subset-sum slice in the common contacts. In particular a common-contact vertex with coefficient f^T K e_v=0 may be flipped at both ports without changing any future gap read. This is a concrete second route after the cost of copying suffix digits: move inside the kernel rather than spend another address cancelling drift.

This does not make the center linear. For the literal nonport rows (0,53,57,83,225,201,135,15) and masks (34,68), the common contacts at v=0 and v=3 have zero coefficients (bit v is 2^v). The four configurations are:

| Common contacts flipped | Masks | Literal path counts | Center |
|---|---|---|---:|
| Neither | 34,68 | 2325,2327 | 2326 |
| 0 | 35,69 | 12263,12265 | 12264 |
| 3 | 42,76 | 2805,2807 | 2806 |
| Both | 43,77 | 11623,11625 | 11624 |

Every pair has the entire gap state (2,0,0,0,0). Vertex 3 is genuinely internal to the nonport graph (outdegree four). Flipping it raises the center by 480 with the other contact unchanged but lowers it by 640 after flipping vertex 0. The mixed difference is −1120. Even switches independently neutral for the *difference* cannot be presumed independently additive for the *sum*.

This refines the selection problem to a tangible search: make the exact neutral slice large, then control its center interaction terms so that actual intervals overlap. A proof of gap neutrality alone cannot justify a greedy digit selector. For this example, reconstruct each pair with the profile/host/count functions in two_address.py; both path counts were independently computed by literal subset DP.
