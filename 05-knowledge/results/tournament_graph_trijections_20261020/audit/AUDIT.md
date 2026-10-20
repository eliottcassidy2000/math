# Independent audit: tournament/graph trijections and cut-Cayley aggregation

2026-10-20. **PASS. No mathematical correction or scope narrowing is required for the audited packet.**

The audited proof note is the 369-line `TRIJECTIONS_AND_TREE_GAUGES.md` with SHA-256
`4a19202298e373dff64cb9ca4b887642c415b560102ff40d34d01735ce0660e7`.
The verifier hash is
`f489e1d8f7074d9854ede3b1b8c624be71135c44ed865052e6c546a7fc4e7c31`.
The producer packet was not modified. This is a proof-and-finite-computation audit,
not a formal proof-assistant certification or an independent literature-priority claim.

## 1. Verdict and publication boundary

The following pass as stated:

- The arbitrary spanning-constraint-graph theorem, with cycle rank
  `beta=m-n+c`, nonempty fibre size `2^(c-1)`, image size
  `2^(N-m-c+1)`, and full one-to-one coverage exactly for spanning trees.
- The order-four path/triangle free-host chart, its inverse, the obstructed
  star free-host chart, and the transported single-edge-flip operation.
- Directed-agreement path incidence and the inherited score-fibre Gram
  interpretation, with no path-witness multiplicity accidentally discarded.
- The all-order weighted cut-partition aggregate, its connected-component and
  bipartiteness classification for `n>=3`, uniform positive weights, and
  generator injectivity from `n=8`.
- The explicit non-equitable/nonlifting limitation. This is essential, not
  merely cautious wording: an independent order-three counterexample below
  shows that connectedness of the aggregate does not imply connectedness of
  the original labelled graph.

Do not advertise these results as original/black/regular-return metagraph
connectivity, a natural unlabelled tournament/all-graph bijection, or a new
proof of inherited spanning-tree/Gram mechanisms. The packet does not do so.
No new theorem identifier is needed for publishing this results packet.

## 2. Adversarial proof audit of the new Cayley theorem

### Incidences, factors of two, and loops

Fix a labelled switching class and an ordered Hamiltonian word. Its prescribed
oriented tree has exactly one realizing tournament in that class. Reversing the
complement of the underlying path adds `[E minus P]` in the cut quotient.
There are exactly two ordered orientations for each undirected Hamiltonian
path. Thus the coefficient `2` in equation (24) and the row sum `n!` are correct.
The agreement graph recovers the underlying path and its direction, so different
words cannot represent the same ordered tournament pair.

For `n>=3`, a zero generator would mean that the complement of a spanning path
is a cut. The path would then be a disjoint union of the two within-shore cliques,
contradicting connectedness unless one shore is empty. That remaining case
would require the path to equal `K_n`, which is impossible for `n>=3`.
Consequently there are no loops, including no loops introduced by aggregation.
The qualifier `n>=3` is indispensable and is already explicit: at `n=2`, the
path complement does nothing and the one quotient vertex has diagonal weight 2.
Order one is outside the packet's declared domain.

Relabelling acts linearly on displacement vectors modulo the cut space and
transitively on undirected Hamiltonian paths. Hence every nonempty generator
fibre has the same size. This action is distinct from the affine relabelling
of tournament state coordinates, and the note distinguishes them correctly.

### Difference span and connected components

For any three distinct vertices `a,b,c`, choose a Hamiltonian word starting
`a,b` and ending `c`. Its cyclic rotation from `a,b,...,c` to `b,...,c,a`
changes precisely the two incident edges `ab,ac`. Such adjacent-edge pairs
span the even-edge-weight subspace because the line graph of `K_n` is connected.
All path differences have even weight, so the equality of spans is proved in
both directions. This includes the boundary `n=3`.

At even `n`, a singleton cut has odd cardinality, so the image of the even
subspace is the whole cut quotient. At odd `n`, every cut has even cardinality,
so parity descends and its kernel has codimension one. The generator weight
`binom(n-1,2)` is odd exactly when `n=3 mod 4`. These facts prove full span
except for `n=1 mod 4`, where the span is exactly the parity hyperplane.
Therefore there are exactly two, not merely at least two, exceptional components.

### Bipartiteness

In a Cayley graph on an elementary abelian 2-group, bipartiteness of a connected
component is equivalent to a linear functional taking the value one on every
generator. Such a functional must vanish on generator differences. For even
`n`, and on the exceptional `n=1 mod 4` component subgroup, those differences
already span the subgroup, ruling out the functional. At `n=3 mod 4`, parity
is the required functional. This proves every direction of the classification;
there is no inference from connectivity alone and no reliance on finite ranks.

### Injectivity and sharp threshold

Equality of two generators means that the symmetric difference of their paths
is a cut. Every vertex of this difference has degree at most four. At `n>=9`,
a nonzero cut has a shore of size at least five and hence has a vertex of degree
at least five. At `n=8`, only `K_(4,4)` survives that degree bound; however,
an endpoint of either path has difference degree at most three. Hence that
cut is also impossible. The empty cut forces equality of the two path edge sets.
The generator degree `n!/2` and edge weight 2 follow, without hidden multiplicity.

The order-seven collision in the independent receipt proves sharpness directly:

- `P = (1,6,2,4,7,3,5)`;
- `Q = (2,3,1,4,5,6,7)`.

Their symmetric difference is the complete cut between `{3,4,6}` and
`{1,2,5,7}`. They are distinct undirected paths with the same generator.

## 3. Constraint, chart, and incidence checks

The component normalization count divides `2^c` potential assignments by exactly
one global constant ambiguity, not by one ambiguity per component after global
cut formation. This is why the multiplicity is `2^(c-1)`. Every slice is nonempty
before quotienting, and the stated image exponent is nonnegative. Requiring
surjectivity gives `beta=0`; requiring injectivity gives `c=1`; together these
are precisely the spanning-tree condition, including at order two.

At order four, the free host and its prescribed complement are never interchanged
silently. Path and triangle free hosts have tree complements. A star free host
has a triangle-plus-isolate complement, with `(beta,c-1)=(1,1)`, hence four covered
classes with two representatives each. Changing the prescribed triangle pattern
changes the affine parity sector but cannot change the rank or multiplicity.

The map `(a,b,c) -> (a,b+c,a+b)` and inverse
`(u,v,w) -> (u,u+w,u+v+w)` are mutually inverse over F2. Sending `(1,1,1)` through
the linear map gives `(1,0,0)`, verifying that path-free complement transports
to one triangle-edge flip, not triangle complement. An independent subset-DP
Hamiltonian count verifies the displayed same-cut-class hostile `H=5 -> H=1`.

The observer graph chart is invertible only with the observer retained; its
relabel action is affine, while its displacement action is linear. The cyclic
triple polynomial is correct and must not be mistaken for ordinary graph
triangle count. The empty-graph fixed-point obstruction to an equivariant
all-labelled-graph/all-labelled-tournament bijection is valid for every `n>=2`.

Agreement edges determine the directed Hamiltonian word uniquely. This remains
true at `n=2`, although the original operation then has loops. The endpoint score
identity gives the same ordered endpoints for two score-return paths, and the
net reversal set is their path symmetric difference. The inherited Gram identity
therefore retains raw `H` on the diagonal. Its order-five icosahedral interpretation
uses the verified existing source, not the unrelated twenty three-edge graphs.

## 4. Independent non-equitable and nonlifting hostile

At `n=3`, the original labelled path-complement graph is exactly two disjoint
copies of `K_(1,3)`. Each labelled switching class contains three degree-one
transitive tournaments and one degree-three cyclic tournament. Thus the cut
partition is not equitable.

Its aggregate is nevertheless `K_2` with weight 6 on its only adjacency. It is
connected, while the original labelled graph is disconnected. This gives a
concrete minimal illustration of the packet's no-connectivity-transfer warning.
It is supplementary audit evidence, not a requested change to the frozen note.

## 5. Computational and provenance evidence

The original verifier was replayed in ordinary and optimized Python modes.
Both produced the exact existing receipt, byte-for-byte:

- checks: **416,131**;
- receipt SHA-256:
  `59a9cdf7d35f6e79571be2fb33b686f0a8065938477e00256454cb9116f66e03`;
- ordinary replay approximately 4.24 seconds; optimized replay approximately
  3.79 seconds in this workspace.

Independent controls import no producer code. They use elimination against an
explicit cut basis rather than the producer's direct triangle encoder, direct
BFS through order six, independent constraint-equation traversal through order
four, and subset dynamic programming for the literal Hamiltonian-count hostile.
They pass **12,612** live checks in both normal and optimized modes, with
byte-identical receipts. They enumerate generators only through order eight,
not all tournaments there. No large new experiment was run.

For `n=3,4,5,6,7,8`, respectively:

- generator degrees: `1,6,15,180,1260,20160`;
- positive weights: `6,4,8,4,4,2`;
- quotient ranks: `1,3,5,10,15,21`;
- generator-difference ranks: `0,3,5,10,14,21`.

All six entries of the producer `SHA256SUMS` pass. All ten inherited-source hashes
match the recorded checkout `f3e658d7cc105d5e367e4dd565812a81e127a6fc`.
All eight Markdown source links resolve when published at the destination in
`PUBLICATION_MANIFEST.json`. That destination was absent when checked.

The repository remained read-only and clean during this audit. The targeted
message-peek command could not run because `.machine-id` is absent; no registration
or repository change was made to work around this administrative prerequisite.
This does not affect the proof, exact replay, or source-hash findings.

## 6. Reproduction and optional editorial additions

From a directory containing the two staged folders:

    python3 session17-trijections/verify_trijections.py > replay.json
    cmp replay.json session17-trijections/verification.json
    python3 -O session17-trijections/verify_trijections.py > replay.optimized.json
    cmp replay.json replay.optimized.json
    python3 session17-trijection-audit/independent_controls.py > independent.json
    cmp independent.json session17-trijection-audit/independent_controls.json
    python3 -O session17-trijection-audit/independent_controls.py > independent.optimized.json
    cmp independent.json independent.optimized.json

No edits are required. If the producer chooses a later revision, a one-sentence
`n=2` loop note and the explicit `n=3` disconnected-original/connected-aggregate
hostile would be helpful pedagogically, but neither repairs a current defect.
Any revision should receive fresh hashes and another exact manifest check.
