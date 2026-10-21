# Independent bounded review: tetrahedral Fano duality

Verdict: PASS, with no mathematical correction required. This is an elementary, exact duality bridge between the two specified tetrahedral constructions, not a novelty claim about Fano planes. Reviewed on 2026-10-21. No repository or Git operations were performed.

## Independent proof

Write the four affine vertex labels as x_0=1, x_1=2, x_2=4, x_3=0 in V=F2^3. Put omega=sum_i x_i=7 and v_ij=x_i+x_j. The six edge labels are precisely V minus {0,omega}. Edge-plane lines are nonzero kernels L_a={v != 0 : a(v)=0}, for nonzero a in V*. These are exactly the four face triples and the three opposite-edge-pair triples containing omega.

For a in V* define C_a by its edge indicator a(v_ij+omega). At any vertex i, the sum of the three incident vectors v_ij+omega is zero: the sum of their v_ij terms is omega, and the three added omega terms sum to omega. Every C_a therefore has even degree at each vertex. The map a -> C_a is linear and injective, since the vectors v_ij+omega are the six nonzero vectors other than omega and span V. Its image is all of the dimension-three binary cycle space of K4.

Let s=a(omega). If s=1, C_a consists of the edge labels in L_a; omega is not on L_a, so its three edges form a face triangle. If s=0, L_a contains omega and two opposite edges, and C_a consists of the other four edges. Thus the coordinate map gives exactly the stated face rule and omitted-opposite-pair rule, not an identification of quadrilateral supports with three-point lines.

For distinct nonzero a,b, the three cycle points C_a,C_b,C_(a+b) form a cycle-plane line. Their edge-plane kernels meet in the unique nonzero vector annihilated by the two independent covectors a,b. Conversely, three distinct kernels concurrent at v have normals equal to the three nonzero elements of the two-dimensional annihilator of v, and those normals sum to zero. This proves both directions of the incidence reversal. It also identifies all seven cycle lines with all seven edge-point pencils. The purely combinatorial description is the same: three quadrilaterals give the omega pencil, and each of the six face pairs gives the pencil at its common edge.

For a vertex permutation p, its action on edge differences extends to an invertible linear map T_p on V, and T_p fixes omega. One direct verification uses the affine map carrying the four affinely independent x_i to x_(p(i)); its linear part is T_p. Consequently the induced covector action is a -> a composed with T_p^{-1}, not generally a -> T_p a. The identity T_p(v_e+omega)=v_(p(e))+omega shows that C_a and L_a transform under precisely this contragredient action. The cycle-point-to-edge-line map and the cycle-line-to-edge-point map are therefore S4-equivariant.

The natural edge-point action has one fixed point omega and one six-element edge orbit. The cycle-point action has a four-element face orbit and a three-element quadrilateral orbit, with no fixed point. Thus there is no S4-equivariant bijection of the two point sets. This does not contradict the duality, whose point-to-line direction is essential.

## Controls and audit limitations

The candidate's checker passes 276 explicit checks in both normal and optimized Python. Its reported orbit sizes are literal output values rather than computed orbits; this is not a flaw in the written proof, but the independent checker computes the orbits directly.

`check_independent.py` does not import or execute candidate code. It derives the cycle space by testing even degree on all 64 edge subsets, derives edge-plane lines as coordinate kernels, verifies all 35 triples and all 49 incidence flags, and checks all 24 tetrahedral actions, including the contragredient coordinate action and equivariance of the pencil map. It also enumerates all 64 tournament orientations, finding all 24 strong orientations and verifying the stated three-support noncollinearity claim and its dual nonconcurrence. No circulation or tournament enumeration claims outside that one sentence are re-audited here.

Three hostile controls pass: omitting the omega translation gives four noncycles; each uncomplemented quadrilateral support fails to be an edge-plane line; using the ordinary vector action on covectors instead of the contragredient action fails in 108 of 168 permutation/covector cases. The independent checker passes 1,364 explicit checks, and its normal/-O receipts are byte-identical.

The proof carefully separates incidence from extra orientation/circulation information. It claims no canonical transfer of Fano-line cyclic orientation to directed-cycle feasibility, and makes no metric, sphere-feature, unit-distance, clock, or trijection identification. No further research expansion is needed for this bridge.

## Reviewed source binding

- `session23b-duality/PROOF.md`: SHA-256 `af297f77570dc82643e1453eaccc3f7b6976001e4dc016996a07a01d96d43729`.
- `session23b-duality/check.py`: SHA-256 `bc531f4099b72dce0671c692158b06d46870a8175d9a7ee68cba8e25af12f695`.

The independent checker rejects drift from these reviewed bytes. `manifest.json` binds the review and controls. A final wording-only promotion from awaiting review to proved with independent review was verified by replacing the exact new status sentence with the old one and recovering the prior SHA-256 `74d88ad82898d6740eeb41d1ebdbd0ade08a37c5bd19c5a170736d29e0a52ef4`. No mathematical text or candidate-checker bytes changed. The binding above is to the promoted final version; any subsequent source edit requires another explicit rebind.
