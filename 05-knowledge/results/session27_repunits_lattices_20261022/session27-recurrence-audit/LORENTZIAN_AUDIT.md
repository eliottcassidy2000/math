# Independent audit: 25-vertex Lorentzian score lift

Date: 2026-10-22 UTC. Verdict: PASS. No mandatory mathematical correction. No priority claim. No repository/Git writes. This is a separate addendum to AUDIT.md and does not broaden its prior review scope.

## Audited inputs and reproducibility

Reviewed LORENTZIAN_SCORE_BRIDGE.md, check_lorentzian_score.py, and lorentzian_score.json in the root's session27-root directory. The requested filename check_lorentzian_score.json does not exist; the actual preserved result is lorentzian_score.json. Re-running the script reproduced the actual JSON byte for byte. Input and output hashes are in LORENTZIAN_SHA256.txt.

The independent lorentzian_independent.py imports no reviewed code. It directly checks all 40320 permutations of each eight-vertex tournament, builds both padded 25-vertex tournaments, counts all their triangles directly, checks their identical scores and norm, and exercises 26 explicit root/quotient representatives. It does not enumerate 25-vertex Hamiltonian paths. Output: lorentzian_independent.json.

## Score identity and lattice membership

Every transitive triple has exactly one source, so its count is sum_i binom(s_i,2). Therefore c3=binom(25,3)−sum_i binom(s_i,2). Since sum_i s_i=300 and 2 binom(25,3)=4600, sum_i s_i²=4900−2c3. Subtracting 70² gives Q(v)=−2c3 exactly.

The coordinate model is correctly stated: all integral or all half-integral coordinates, with integral pairing against (1/2,...,1/2;1/2). For every tournament lift that pairing is (300−70)/2=115. Thus the vector belongs to the actual even Lorentzian lattice. Time coordinate 70 is nonzero, so nullness is a genuine nonzero null-vector statement. Absence of a directed triangle is equivalent to transitivity for tournaments. Hence every nontransitive tournament lift is timelike and every transitive lift is null. No inference about arbitrary graphs is warranted or asserted.

For a transitive tournament the scores are precisely 0 through 24 in some order. Spatial coordinate permutations are lattice isometries. The cannonball sum is 4900, so all such lifts lie in the orbit of the displayed w. The vector w is primitive: a lattice vector r0=e0−e1 pairs with it as −1.

## Classical source and quotient identification

Directly inspected [Borcherds–Conway–Queen–Sloane, A Monster Lie Algebra?, 1984](https://neilsloane.com/doc/Me99.pdf), printed page 76. It explicitly gives the coordinate model, w=(0,1,...,24;70), and identifies roots of norm 2 pairing −1 with w, with the Leech lattice under the difference metric. Its reference [6], printed page 79, attributes this to Conway–Sloane, Lorentzian forms for the Leech lattice (1982). The addendum's attribution is accurate; the 1982 paper itself was not separately inspected here. No Monster-group conclusion follows for tournaments.

The quotient argument is correct, including its sign. Write q=r0+y+a w with y in w-perp. Since (r0,w)=−1, q²=2+y²+2(r0,y)−2a. Thus q²=2 exactly when a=(y²+2(r0,y))/2. This is integral by evenness and integral pairing. Replacing y by y+k w replaces a by a−k and leaves q unchanged. Conversely every relevant root q yields y=q−r0 in w-perp. These maps are inverse between the roots and the quotient w-perp/Zw. Differences have precisely the quotient's induced norm. Choosing r0 supplies the origin of the resulting affine metric lattice. The construction consequently yields Leech, not an arbitrary rank-24 Niemeier lattice.

If v²<0, v is not in v-perp. Therefore Zv is not a subgroup of v-perp, and the displayed quotient by Zv is not defined. One can study the timelike orthogonal complement by other constructions, but it is not the same null quotient and would not have the asserted rank-24 interpretation. The addendum correctly refuses that extension.

## Padded collision

Both supplied eight-vertex bitmask lists encode valid tournaments and have identical labelled scores (4,4,4,3,3,3,3,4). Independent full-permutation counts give 657 and 643.

Attach vertices 8,...,24 with arcs from each new vertex into all eight old vertices, and with i→j for 8≤i<j≤24. The old scores stay unchanged. New scores are 24,23,...,8. Direct enumeration of the 2300 triples gives c3=20 for each padded tournament, and Q=−40. No mixed triple is cyclic, and the new chain is acyclic.

Every Hamiltonian path must start with that entire chain in its unique forward order: no old vertex points to a new vertex, and no new vertex can precede an earlier chain vertex. Its remaining suffix is an arbitrary old Hamiltonian path. This is an explicit bijection, proving preservation of H without any 25-vertex Hamiltonian enumeration. Both 25-vertex lifts are exactly the same labelled Lorentzian vector, not just equal in norm or orbit. Therefore neither that vector nor its norm can determine H.

Wording improvement applied and audited: the manuscript now says “a transitive chain of17 vertices dominating the original tournament.” Only the first new vertex is a universal source in the final tournament. Reversing precisely this replacement reproduces the previously audited input SHA256 b2af3aef6c8ca0963c950258af1984066f61c49383064f39e7bebd9b18dcd0b6, proving this was the sole input change. No computational rerun was necessary; the PASS verdict and prior controls remain valid.

## Provenance and limits

The manuscript explicitly attributes the concurrent proposal to incoming MSG-1913 at commit a08f25ecbf8abd97959305f0836c1c25890edf65 and does not use reserved theorem IDs as proof dependencies. This audit verifies that attribution is present and the algebra stands independently; it does not independently authenticate the incoming message or remote commit history. The root supplied that provenance. No claim of first discovery is certified.

The concluding reference to a finite-context reconstruction theorem is contextual and was not re-audited here. The audited score/lattice and collision arguments do not depend on it.
