# Fixed-deficiency kernels, endpoint-changing detours, and exact boundary games

Research continuation, October 20, 2026. Initial integrated base:
`7f8df95a96ea`; incoming work through `9852737522` was integrated at close-out.
The isolated worktree preserves the unrelated shared checkout.
This packet distinguishes proofs, cited inputs, finite computations, and
remaining conjectures. No claim of literature novelty is made.

## Strongest result: every fixed deficiency has a finite kernel formula

[General fixed deficiency](general_fixed_deficiency.md) proves a structural
extension of the D3 mechanism. For fixed D>=0, put m=floor(D/2). In the stable
range k>=m+1, the number N_D(k) of capped weighted pair-free chain cores with
feedback at most k and exactly 2k+1-D heavy coordinates satisfies

    N_D(k) = (2k-2m-1)!! P_D(k),
    degree(P_D) = floor(5D/2),
    leading coefficient = 2^floor(3D/2)/D!.

The polynomial is computable from kernels with at most 4D endpoints and
2D defect edges. The owner now compares q, the **full** vertex-weight word,
and the edge list on full vertex positions. This preserves ownership under
arbitrary matching insertion without dropping ports or permitting clone
mergers. Large light vertices have fixed positions equal to their indegrees,
which proves independent specialization of their weights.

This is a **PROVED** fixed-D theorem with independent proof audits. It is
not an unrestricted scalar recurrence for A000568. At the zero-heavy
boundary D=2k+1, arbitrary pair-free cores can still occur. Those boundaries
and the growth of the interface with D are the remaining counting problem.

The new owner was compared on all 36,942 small D3 records: 392 individual
decisions change while the class counts agree. Its controls also cover
7,105 matching-insertion/joint tests and 313 independent large-light
substitutions, including unequal values. These finite tests supplement
the general proof rather than supply its unbounded quantifier.

[The completed D3 arithmetic](coefficients.md) gives, for every k>=2,

    N3(k) = (2k-3)!!/6 *
       (16k^7+64k^6+8k^5-140k^4-22k^3+135k^2+61k-110).

The separate value is N3(1)=5. C4 is a **FINITE-EXACT** owner computation;
C5/C6 use proved activation shortcuts. A complete separate C5 owner replay
agrees on all 1,988,280 primitive records. The formula uses the inherited
D3 owner; individual kernel coefficients are not identified with those of
the newly chosen general owner.

## Strong mixed blue detours: a literal three-move family

The incoming [nonstrong theorem](../spectrum_basins_20261020e/outer.md) and
[pure-blue common hubs](../spectrum_basins_20261020e/pure_hubs.md) reduce the
black-connectivity problem to strong mixed-to-mixed blue lines, relative
to the stated **CITED** combined-connectivity input.

[Endpoint-changing detours](detours.md) proves a three-black bypass for the
aligned blue word `(a,V,z)` in `a -> V -> z -> a` whenever |V|>=4. It preserves
the complete middle, every orientation guard, and the literal labelled
target. Open blue middles of order three also work. The closed C3 middle
reaches the pure-blue regular pentagon and is a necessary exclusion.
All 4,693 applicable blue middle presentations through order eight are
replayed exactly. The proof is uniform in middle size.

The inherited fixed-endpoint hub fails on genuine strong five-vertex
examples. The new words release those endpoints. A separate even-order
hostile corrects a tempting symmetry shortcut: consecutive whole-cycle
reflections compose to a **two-step** rotation, not a one-step rotation.
General strong blue detours and all-order black-core connectivity are OPEN.

## The clarified complement game is completely classified

Player 2 finds a path among the uncoloured edges. Let H mean Hamiltonian
path and W mean spanning walk, allowing repetitions. W is connectivity.
For n>=4 the sharp minimum number of red edges is:

| Red must lack / blue must lack | Minimum red edges |
|---|---:|
| H / H | n-1 |
| H / W | n-1 |
| W / H | b(n) |
| W / W | impossible |

Here `b(2m)=min(4m-5,binom(m+1,2))` and
`b(2m+1)=min(4m-3,binom(m+2,2)-1)`. At n<=3 no variant admits a win.
[Game proofs](game.md) use elementary complement/closure arguments and one
**CITED** sharp Erdős bound. All labelled graphs through order six are
checked; literal path permutations independently agree through five.

The incoming [independent game analysis](../complement_path_games_20261020/AUDIT_AND_PROVENANCE.md)
agrees on these budgets and adds exact rooted path-cover join criteria.
Its stronger hostile has **both colours connected and nontraceable** at
every n>=7; our independent replay confirms the construction, and the
two censuses exclude it through six. This is a concrete further connection
to the repository's full path-cover substitution interface.

The undirected star and triangle-plus-isolate also expose a port loss:
the line graph of K1,3 is Hamiltonian, but the star has no Eulerian trail.
[Walks and knots](walks_knots.md) supplies the exact directed-line/dart repair
and a separate torus-family residue graph with integer winding labels.
Its cycles carry `(p/d,q/d)`, d=gcd(p,q); forgetting carries loses q beyond
its residue modulo p. This recovers the precise theme of THM-2646 without
identifying arbitrary tournament path counts with knot parameters.

## Six-colour boundary forcing beyond the triangular lattice

[Exact palette forcing](forcing.md) constructs a 16-point unit-distance
gate in Q(sqrt3,sqrt11). Nine prescribed ports force a central vertex into
their three-colour palette, with 24 extensions for each output colour.
Two extra private ports give an 18-point graph where **every proper**
six-colour boundary assignment extends. A third gives a 19-point graph
with an explicit nonextending proper boundary assignment.

All three graphs are exactly four-chromatic. The result is a boundary
relation, not a six-chromatic graph. The extension threshold follows from
an elementary induced-diamond list-colouring proof; all geometric contacts
are checked exactly. An independent SymPy audit reconstructs all 171 pair
distances and confirms the palette counts and four-colouring. The remaining
task is an actual geometric palette synchronizer, not abstract port merging.

## Inheritance and research board

The anchor is the fixed-deficiency recurrence and its exact owner. The
niche is strong endpoint-changing transport. The wildcard is the geometry
of list constraints and the arithmetic hidden by finite graph quotients.
The board has six concepts: exact outside ports; canonical owners; bulk
insertion; endpoint charges; boundary palettes; integer carries.

The closest proved mechanism is the D3 matching-completion theorem; the
canonical hostile is equal local gaps with different outside partners.
The corrected near misses are deficiency-changing clone contraction,
fixed endpoints without a black alternative, and a cyclic rotation with
two parity orbits. The least-used relevant sidecars recovered here are
full vertex positions, degree-list colouring, and THM-2646's integer braid
height. The packets state each map's source, target, preserved predicate,
lost information, and cheapest hostile control.

The concrete next targets are:

1. Sum or classify the **zero-heavy boundary**, where the all-fixed-D law
   leaves arbitrary core complexity. Counting relocation-connected fibres
   still requires exact fibre sizes and transport; connectedness alone
   supplies no isomorphism-class recurrence.
2. Extend the endpoint detour beyond the closed score pair `(n-2,1)` to
   `(n-3,2)`, retaining the extra mirrored boundary incidences.
3. Realize the palette gate's input relation geometrically while retaining
   the actual induced boundary graph and a proper-colouring hostile.
4. For the Collatz comparison, keep exact integer affine carries on the
   excursion forest. A finite carrier with unbounded labels can encode the
   arithmetic, but no reset-paying decreasing rank or universal certified
   entry theorem follows from the present constructions.

## Reproduction

Each proof links its own checker and gives its exact universe and controls.
The main checkers use the Python standard library; `independent_audit.py`
uses SymPy 1.14.0 as a separate algebra/geometry path. No Lean proof is
claimed. Normal and optimized runs agree for the principal receipts.

    python 05-knowledge/results/kernels_games_detours_20261020f/verify.py
    python 05-knowledge/results/kernels_games_detours_20261020f/verify.py --replay

The manifest pins local sources, receipts, and imported owner dependencies
using SHA-256 after CRLF-to-LF normalization. The complete C5 owner audit
is separately reproduced with `coefficients_check.py audit-five`.
