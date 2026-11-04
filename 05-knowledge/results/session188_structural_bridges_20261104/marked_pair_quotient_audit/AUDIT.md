# Independent marked-pair quotient-converse audit

**PASS.** Reviewed the complete proof and checker; independently recovered grouped coefficients, the converse signs, the ordered-cut formula and the retained literal witnesses. The result applies to the stated marked-pair setting and does not assert arbitrary strong outer-module transfer.

## Alias and converse bookkeeping

The four retained Boolean variables are a_q,a_A,c_p,c_B, with p,q outside the disjoint exterior unions A,B. The original response has bidegree at most(2,2). Therefore any physical monomial that becomes this grouped degree-four monomial must contain exactly a_q, exactly one a_i with i in A, exactly c_p, and exactly one c_j with j in B. There are no extra terms from two aliased A names, from two aliased B names, from lower degree, or from omitted higher degree. This justifies summing the named quartic coefficients once each.

After reversing only p->q to form W, the whole converse W-op has the target quotient and the desired p->q interior. Every other child is conversed. On the prescribed grouped queries each such child is seen uniformly by both ports, so restoring it preserves the actual Hamilton count exactly. Thus F_U=F_(W-op) on the entire four-variable family. Whole-converse quartic invariance can then be applied before grouping; substitution cannot raise degree.

The twin derivative's sign is correct. In T a contributing word contains the marked segment j,y,p,q,x,i, producing C_ji with two complement-factor minus signs. In W it contains i,x,q,p,y,j, producing C_ij with two minus signs. Both quartic signs are positive, so their contrast is C_ji-C_ij. The pair's exterior incidences are not tested by these marked segments. Boundary singleton cover paths are included and both cover slots are nonempty.

## Ordered-cut positivity and scope

For R=A->B, the first cover count is H(A)H(B). In the opposite endpoint order each of the two ordered slots may cross the cut or remain on its required side. The four cases give H(A)H(B),q2(A)H(B),H(A)q2(B),q2(A)q2(B), without an extra ordering factor. Their difference is positive precisely when at least one side has at least two vertices. The two-singleton zero case is retained and is not called a response collision.

An ordered SCC cut of Q-z lifts to full unions of child blocks. Its exterior order is therefore genuinely A->B even when the child interiors are nontransitive. The example Q=(14,28,8,16,5) has the stated physical orientation, is strongly connected and prime, and deletes to TT4 at vertex4. With all children TT2, the grouped coefficient is224. The regular-five example only lies outside this particular ordered-cut corollary; it is not asserted to survive every other same-session separator.

## Independent computation

`independent.py` uses predecessor Hamilton-path DP and, separately, direct Hamilton-word cuts for ordered two-cover matrices. The latter is independent of the producer's subset endpoint convolution: concatenate an ordered cover, retain every required cut, and count the actual ordered endpoints without testing the connecting arc.

The independent census covers every labelled quotient on3 or4 vertices, every choice of the marked pair child, and every ordered choice of disjoint nonempty exterior groups, leaving any remaining child unmarked. All other children in this exhaustive small census are singleton blocks. This gives:

- 3,120 grouped sixteen-query batteries;
- 4,656 underlying named-quartic contributions, checked against complete retained physical response vectors;
- exact agreement of the grouped coefficient, pair-flip derivative and literal oriented cover skew.

Every one of the producer's18 additional batteries is independently reconstructed from its quotient and physical children, including nontrivial child restoration and the ten-name prime witness. The equality between U and W-op is checked on every actual query, and the grouped T-W coefficient is checked separately. Across the small census and retained larger batteries there are101,568 literal integer count comparisons, with repeated queries cached. The18 larger grouped coefficients include0,positive and negative controls, and the224 witness.

All121 ordered-cut pairs of arbitrary labelled tournaments of sizes1 through3 are checked by direct word cuts, including the sharp singleton boundary. Normal and optimized independent receipts are byte-identical, as are the frozen producer receipts. Full cubic/quartic response information is used; no low-degree equality is substituted for a full-function comparison.

## Frozen pins

- `marked_pair_quotient/PROOF.md`: baba1aeb86fddf83f000a63f4da0ccff1418a3f6add077c7495c35bfe85f6e09
- `marked_pair_quotient/verify.py`: dd20526244816c4477099d034be062a8a3b42d89512c61ba226c65d5e8157cea
- Producer receipts: 41a5992ca0c3678bd6a701f633fee313cbe15391a62bce544abf16dbb2a89882
- Incoming `response/PROOF.md`: dcd7a33664f63039b937582f241910e98ccb3fa341f3ba3999f626f5c61b51ec
