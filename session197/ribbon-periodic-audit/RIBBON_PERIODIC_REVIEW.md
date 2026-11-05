# Independent review of the fixed-block ribbon obstruction

Verdict: **PASS**, bound to the frozen `session197/ribbon-periodic` snapshot:

- `PROOF.md`: `69d23585eea361ab1b7d104cd80384eeb3ece659bd147b483efa268d0c6577f7`
- `verify.py`: `8f9a536a86f563ab79d54031e587539af2fef1d114bd6ddb144c4cf437def3ef`
- `representation.json`: `d418468ec6869b2eab56a8406b48b19e78ac64993da686fbe767e5e264cd8b5a`
- `eigenvalue2-counterexample.json`: `6baabd5b2a4dd5e907e854f34c5753be2d2293eaff08dc2b08fe96cf7cab92dc`
- `VALIDATION.json`: `c8cbef1f039b8efe2f23ca7bee86a3ae11549357359f354c8049d7400ba4d3d0`

## All-order proof review

The forward-sink reset argument is sound for nonempty old tournaments. Two distinct old path components cannot be joined through a new sink and would leave two forgotten starts; a single old Hamiltonian path has continuation options independent of its interior. The verifier's complete nonempty frontier-basis check independently confirms the required rank-one response.

For the periodic triangle case, a translated directed triangle has diameter at most 2w. When the period length L>w, overlap can occur only between consecutive translates and only by identifying the last vertex of one triangle with the first of the next. In that case every adjacent pair in the displayed union is one of the forced backward edges, and every nonadjacent pair has natural index separation at least L>w. There are therefore no unspecified short cross-arcs: the union really is R_(2m+1). The tribonacci inequality used to obtain 3^m is correct, including its n=1 boundary. In the disjoint case, the separately proved deterministic path-merging injection applies regardless of the cross-arcs.

In the triangle-free suffix case, the finite set U is correctly taken through the last of the first w suffix vertices in the suffix's transitive order. No suffix vertex after natural position 2w can precede a member of that set. Every residual vertex has natural position greater than w, so it is dominated by the entire old prefix as well as by U. Thus the ordered-join factorization and eventual constancy hold even for a hostile nontransitive prefix.

The finite representation certificate is used correctly: entrywise M_a B=B R_a and Bf=F suffice for every future output, without claiming minimality or reachability of every signed reduced coordinate. Nonzero modular determinants prove that 2 is not an eigenvalue. Cyclic rotation reduction is valid for square matrix products, including singular ones. Cayley-Hamilton then excludes an eventual nonconstant affine-doubling sequence exactly as stated.

For W=001476, the exact vector relation y_3-77y_2+150y_1=0 propagates by right multiplication with the block matrix. The initial values determine the two-exponential closed form. Thus the factor-two mode is physically present, but its simultaneous factor-75 mode prevents a scalar binary append rule. The counterexample does not contradict the longer-block graph-growth theorem.

## Reproduction and additional independent graph checks

I ran the full frozen verifier with assertions enabled. It reproduced all 952 intertwining row checks, 220 short-block rotation classes, 118 reset basis rows, 504 overlap geometries, 92 hostile-prefix transitive-tail cases, the boundary-deletion counterexample and the exact mixed-exponential certificate. The output JSON retained its frozen digest.

I additionally constructed the actual width-three graphs for W and W^2 and counted them by the separate Held-Karp graph routine: H(W)=9 and H(W^2)=543, agreeing with the formula. The higher-k conclusion uses the reviewed exact vector recurrence, not extrapolation from these two checks.

The verifier deliberately uses Python assertions. It must be run without `-O`, as the packet states. This review does not claim an optimized-mode validation.

No correction or unresolved blocker was found. The result excludes fixed repeatable binary-affine append blocks, not state-dependent controllers or the full span-three spectrum.
