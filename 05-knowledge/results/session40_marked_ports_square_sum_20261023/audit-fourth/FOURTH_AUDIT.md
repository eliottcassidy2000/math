# Fourth audit: late incoming boundary correction

Incoming commit d45865c3323af6d1b3ce44525d6d4789b22778d7 was read directly, including its owned-boundary frontier, aligned-port proof, and MISTAKES correction. This addendum does not alter the first three frozen audits.

## Independently verified strengthening

In an undirected Hamiltonian graph of order at least three, every edge is optional among unrestricted Hamiltonian paths. For a cycle edge, delete that edge to omit it and a different edge to include it. For a chord uv, cutting u–succ(u) and v–succ(v) splits the cycle into two paths; inserting uv joins them into one spanning path containing the chord. Cutting the original cycle omits the chord. The construction does not preserve endpoints or simultaneous required edges.

The root's explicit witness implementation was read and rerun normally and under optimized Python. Its canonical cycles of orders3..64 cover every possible cycle edge and added chord, validating each resulting vertex list and edge membership. Its four-vertex graph with chord0–2 has a Hamiltonian path from0 to1, but none using that chord with those endpoints. The fixed-endpoint firewall is therefore valid. These computational controls support the elementary all-size proof rather than replacing it.

The primary OEIS page https://oeis.org/A090461 was opened during this audit. Robert Gerbicz's comment dated21January2018 explicitly records Hamiltonian cycles at every n≥32. Thus the all-size square-sum edge-role consequence is properly labelled elementary proof plus CITED input. This audit did not reproduce the cited gluing algorithm. Each fixed distinct positive square-sum pair a,b is includable and omittable for every n≥max(32,a,b).

## Aligned ownership and boundary-state scopes

The signed path identity A+B−(J−I)=F−F^T and the joint equitable-refinement peeling proof are correct for common actual vertex alignment. Endpoints appear in round one and successive singleton incidences expose inward neighbors, giving at most floor(n/2) rounds for n≥2. The nonnegative doubly stochastic commutant proof fixes endpoint rows and columns using unique extrema, then propagates along the signed path. Neither theorem finds the missing target or alignment from compressed source-only data.

The agreement-matrix equivalences require the separate in/out-degree guards. Under those guards it is one directed endpoint path plus detached cycles, so nilpotence, determinant1 for I−D, and the full start-to-end walk are equivalent ways to reject detached cycles.

The connected affine-F3 three-boundary section count is structurally correct: one empty relation,27 singleton triples,36 full affine three-point graphs. A first boundary owner fixes its entire connected completion. Empty section relations can conceal different holonomy centralizers, so the inherited420 empty-section systems do not contradict differing residual symmetries. The order31 counts and full composition census remain inherited finite results, not independently replayed here.

## Final integration review

LATE_INHERITANCE.md, SYNTHESIS.md and README.md were read as actual files. They correctly credit the incoming work, close unrestricted coarse single-edge roles, and preserve the unresolved prescribed-endpoint, joint-edge, and compressed source-only wiring questions. The path-clock formula is understood for n≥2 and s≥1; order-one paths are a separate trivial case. The incoming claim of a user's4–12 correction is not presented here as a newly received message; the earlier13–12 parsing remains a valid edge statement.

Verdict: PASS for the integration, with cited and inherited boundaries retained. Exact reviewed source hashes and replay commands are recorded in MANIFEST.json.
