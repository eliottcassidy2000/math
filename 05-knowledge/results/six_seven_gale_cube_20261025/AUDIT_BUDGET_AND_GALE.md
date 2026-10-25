# Independent audit: canonical budget and Gale return

Reviewer: arithmetic_shapes. Date: 2026-10-25.
**Verdict: PASS in the stated scope.** No mathematical blocker found.

Reviewed `BUDGET_AND_GALE.md` and `verify_budget_gale.py` in the current
`six_seven_geometry_20261025` packet. Independently reran the checker:
305 gates pass, including all 720 permutations at u=v=1, with eight
matroid-compatible permutations and exactly two projective matches.

Reviewed SHA256 values: proof
`168ee9b57f09015551ce084e8dce5012057341868f1edf9b779644ee96d4b83f`;
checker `975f86c34a7958ed340c1681bcf5101c086d899a0d84d391d29aa0415464a3d8`.

1. The smooth general-position hypothesis of BG1 is sufficient: at most d
   branch hyperplanes pass through a point, and their local equations are
   independent parameters. Rooting those parameters and taking etale roots
   of units preserves regularity. Distinct hyperplane valuations prove full
   generic Kummer degree, including composite root order r in characteristic
   zero. The monic defining equations are a regular sequence. Adjunction,
   degree, and the Hilbert series therefore have the stated values.

2. Intermediate H^i(O) vanishing is the standard complete-intersection
   Koszul calculation for 0<i<d. Integrality supplies H^0(O)=k, and Serre
   duality at rc=n supplies H^d(O)=k. No special singular arrangement is
   silently put under this smooth theorem. For the actual six/seven-root
   surfaces the inherited A1/crepancy theorem is separately needed.

3. The symbolic Gale return is verified column by column:
   M C has columns -2 A_5, -u A_2, u A_1, -v A_4, v A_3, 2 A_0.
   Its determinant is 2uv. Thus u,v nonzero are exactly the displayed
   invertibility conditions; balanced dimensions by themselves are not
   used as a substitute for projective equivalence.

4. The arithmetic twist is retained correctly. Substituting lambda=M^t nu
   gives C_i dot lambda=D_i L_(p(i))(nu). Their product is
   -4u^2v^2 product_i L_i(nu). Dividing the covering coordinate by 2uv
   leaves w^2=-product_i L_i. Over an algebraic closure this twist splits;
   over Q(u,v) the displayed map retains the minus-one squareclass. No
   assertion about all conceivable alternative isomorphisms is made.

5. The degree-eight complete-intersection polarization, degree-two
   discriminant cover, and degree-16 character quotient are correctly
   distinguished. A repeated K3 label does not establish isomorphism,
   equality of integral cohomology lattices, or derived equivalence.
   The discriminant determinant is elementary and does not require a
   smooth-net theorem at the singular special arrangement.

One wording clarification was suggested to the author: the connection
ledger's root-equation preservation should not suggest that every diagonal
column rescaling lifts over the original field without adjoining square
roots. The detailed twist paragraph already prevents the false double-plane
inference. This is a precision improvement, not a failed theorem.

The Madonna–Nikulin arXiv abstract was independently opened and confirms
the classical net/discriminant correspondence and that isomorphism questions
require extra lattice data:
https://arxiv.org/abs/math/0206158
The Huybrechts PDF timed out during this independent source check; this audit
does not certify its exact example-number citation. It is not needed for
the elementary matrix identity or the proved character quotient.

Integration amendment: the root reviewer made the suggested ledger clarification: individual column scalings require their square roots to lift to root coordinates. The equations, proof, checker and receipt are unchanged.
