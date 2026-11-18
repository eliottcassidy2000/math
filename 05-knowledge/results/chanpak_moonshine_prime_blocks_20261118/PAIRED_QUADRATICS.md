# Direct and nested paired quadratics

The [paired rooted grammar](PAIRED_RECURSION.md) scales a state under repeated attachments. A complementary way to obtain p and 7p together is to match the **visit tensor** of two different quotients. This controls two independently chosen input blocks or, in one branch, just their normalized profile.

With ports 0,1 insert independent blocks L,R. Each below can have just one or two visits per port; as the two marks alone are transitive the absolute max is bounded by the number of ordinary vertices plus one. The replay counts all visits up to that bound and certifies the shown tensor. In visit order (1,1),(1,2),(2,1),(2,2), the lower and upper tensors are:

| Lower masks | Upper masks | Lower tensor | Upper tensor / 7 |
|---|---|---|---|
| (34,32,11,19,39,12) | (6,4,24,51,195,87,143,47) | (25,12,12,6) | same |
| (62,16,26,50,96,70,15) | (6,4,24,51,195,343,143,303,95) | (67,32,32,16) | same |
| (66,64,3,23,39,79,28) | (6,4,336,295,459,343,395,303,3) | (53,30,30,15) | same |
| (190,16,2,22,100,142,175,30) | (6,4,152,115,195,471,391,11,159) | (61,41,23,16) | (61,32,32,16) |

Use H for one visit and B−2H for two. The first three pairs read p=alpha B_L B_R+beta H_L H_R and 7p, for (alpha,beta)=(6,1),(16,3),(15,−7). The fourth reads the analogous (16,−3) values exactly when H_L/B_L=H_R/B_R, since the lower minus upper/7 is 9(H_L B_R−B_L H_R). In particular an identical block in both ports always works. Indefinite and definite integer forms thus arise on the diagonal. No rotation symmetry is inferred from tensor scaling in these four templates.

The matched independent readings have the same numerical form as an affine reader on the ordinal join L→R. Their extra information is the actual two-port geometry and entrywise visit matching, not an assertion of a larger numerical spectrum by themselves. One can group any number of physical inputs into two ordinal groups and obtain p=alpha·product(B_i)+beta·product(H_i) and 7p. In the (16,3) row, independent positive transitive widths yield every 2^m+3 for m≥6, not just even exponents. In the (6,1) row a strong core with (h,B)=(19,256) in both ports gives the certified prime 393,577 and 2,755,039.

## A quadratic stage retaining a genuine seven-rotation

To make the upper symmetry literal, compose a quadratic *seed compiler* with the [paired carousel reader](DIHEDRAL_HEIGHT.md). Given arbitrary Q, put X=C3[Q,Q,singleton]. In the odd-cycle packing formula, either the singleton is uncovered or its cycle visits one directed path from each marked block, with total even support. Keeping even support by signed evaluations q+1,q−1 gives
```
f_X(q) = q [f_Q(q)^2 + (f_Q(q+1)−f_Q(q))^2
                      + (f_Q(q−1)−f_Q(q))^2].
```
Indeed the nonempty support reads subtract the unvisited term, and the factor 2q per odd cycle cancels the half from selecting even parity. At q=1,2 this gives
```
H_X = 2h²+(B−h)²,
B_X = 2[B²+(C−B)²+(B−h)²],   C=f_Q(3).
```
Append one ordinal singleton to X and apply the paired 1036B−209H reader. The result is p and 7p with p=2072B_X−209H_X and an actual C7 action on the upper. The pad matters: without it, X has odd order and the unpadded affine read is divisible by 3 by slot reflection. A literal Q=(126,252,56,112,96,64,132,61) with (h,B,C)=(71,1976,21591) gives prime p=1,624,855,313,581 and 7p=11,373,987,195,067.

These are ways of composing shared input, compare whether one can reach more absolute targets with different input factorizations or complete one of the rooted reachability grammars. They do not show that an arbitrary prescribed prime reaches one of the physical profiles. Replay: `python3 paired_quadratics.py`.
