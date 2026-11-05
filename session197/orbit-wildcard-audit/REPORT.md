# Independent audit: the 14-vertex two-orbit 7p witness

**PASS. No correction required.** This report audits the explicit graph/certificate part of `session197-orbit-wildcard/PROOF.md`. The complementary parity-context theorem audit is in `session197-lowcycle/WILDCARD_REVIEW.md`; its reviewer reports PASS. No original wildcard file was edited.

## Frozen snapshot

- PROOF.md: `69497d610f4ffdb3530ed260cec3618af3d973aadaadbbddd23aaff343424f4a`
- verify.py: `9f354c23dfe8616f73f58cc15f56b6054a5f197bdb51fc944c6550b04d6ceba9`
- receipt.json: `50f7e00c1375e818bdc01f419d0c13b2340d59caee197cd2e2aaeae3e5c31aef`
- witness.json: `b4358dec35824a42a5756bd73a2cd180a3b7fd242870628fb1e54ae5430d36d7`

The independent checker verifies these hashes on every run. Normal and optimized runs pass and have identical output.

## Independent computations

`check_graph.py` imports none of the author's verification code.

1. Reconstructed both Paley blocks and the cross-difference set {1,3} from the definition, then matched all fourteen outgoing rows and checked the tournament axioms.
2. Computed H by a different counting identity from either author DP:

       H(T)=sum_(S subset V) (-1)^(n-|S|) 1^T A_S^(n-1) 1.

   Inclusion-exclusion retains exactly length n-1 directed walks visiting all n vertices. Such a walk visits each vertex once and is precisely a Hamilton path. The implementation first matched literal permutation counting on all 1,099 labelled tournaments through order five. Applying all 16,383 nonempty-subset terms to the fourteen-vertex graph gives exactly

       H=15,144,269=7*2,163,467.

3. Independently trial-divided 2,163,467 through floor(sqrt(p))=1,470. It is prime. Every recursive Lucas certificate node was also trial-checked, its complete factorization of p-1 independently recomputed, and its modular full-order witness checked.
4. Enumerated every nontrivial proper vertex subset as a possible module, using equality of all members' exterior out-neighborhoods rather than the author's outside-vertex criterion. None exists. There are16,368 such candidate subsets in the universe of16,384 subsets.
5. Checked strongness by independent bitset transitive closure.
6. Enumerated the full automorphism group by joint degree-constrained backtracking on all fourteen vertices, not by first assuming the two Paley automorphism banks. All partial assignments were checked against previously assigned arcs. The search visited106 nodes and returned exactly the seven simultaneous translations. Their orbits are exactly A={0,...,6} and B={7,...,13}.
7. Verified outdegrees5 on A and8 on B, H(A)=H(B)=189 by the independent walk formula, and mixed incidences from every opposite-orbit vertex to each full orbit.
8. Verified the resulting path-orbit mass2,163,467 and labelled mass14!/7=12,454,041,600.

## Mathematical review

The structural proof of module-indecomposability is sound. In a regular seven-tournament, a module induces a regular tournament, so its possible nontrivial proper sizes are3 or5. Size5 contradicts the exterior degree3 condition. For size3, an exterior vertex beating the module must lose to all other exterior vertices, and one losing to the module must beat all of them. Each class has size at most one, so four exterior vertices cannot be accommodated.

A module of the full graph intersects either Paley block in a module of that induced block. Containing a full block but omitting an opposite vertex contradicts the mixed cross incidences. Therefore any remaining nontrivial proper module would be a cross-pair; its two total outdegrees would differ by exactly one, whereas the actual difference is three. This excludes all cases.

The short automorphism argument is independently valid: simultaneous translation has order7. Free action on Hamilton words makes the full group order divide7p. It also divides14!, while the certified prime p exceeds14. Hence the group has order exactly7. The direct backtracking check confirms this without relying on that divisibility argument.

The finite parity-context gcd proof was also read for consistency. Its monic recurrence applies to positive lengths1,3,... and2,4,..., not formal length zero. The grouped-square criterion correctly weights coefficients by b for odd input and b^2 for positive even input; bases zero modulo the prime vanish. The paired reviewer owns the detailed context/profile controls.

## Scope

The explicit graph is a strong module-indecomposable symmetric7p tournament whose two moving orbits both induce Paley tournaments. It refutes an unrestricted module-reduction extrapolation, not the earlier threshold-bounded seven-orbit result. No minimum order, minimum H for two-orbit examples, infinite prime family, or full-spectrum closure follows or is claimed.

Reproduce with `python3 check_graph.py` and `python3 -O check_graph.py`, comparing output. The report's PASS applies to the exact frozen hashes above.
