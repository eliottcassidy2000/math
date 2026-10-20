# Source audit: the supplied six-colour claim

2026-10-20. **Status: primary-source scope inspection, NOT a Lean build or complete proof audit.** No demonstrated refutation of the present manuscript is asserted.

## Pinned primary source

The GitHub connector verified `openai/math` main as commit
`1581d4cbdd69a0faf34b2089614fe2a814316a6a`. `source_manifest.json` records
read file paths, Git blob hashes, line counts, and selected declaration positions.
All repository access was read-only. No remote code was executed or installed.

The supplied PDF is
https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf

Its actual TeX source is under `build/sections/`, not a top-level `overview.tex`.
The introduction, transfer, and rigidity source were inspected, alongside the
formal definitions, final theorem, forward-transfer module, comparator metadata,
and root import. The complete geometric dependency graph and Lean kernel build
were not independently checked.

## What the primary source actually claims

The introduction's Theorem `thm:main` concerns every colouring of R², not merely
measurable colour classes. Its intermediate `thm:transfer` asserts, for every
positive finite palette, equivalence between unrestricted proper colourability
and weak measurable colourability. The latter excludes same-colour unit pairs
almost everywhere for the product of planar Lebesgue and angular Haar measures.
This is not pointwise properness of the measurable representative.

The paper then claims no weak measurable five-colouring. Its finite graph is
existential, obtained by compactness after the unrestricted theorem. No explicit
six-chromatic coordinate set or five-colour UNSAT certificate is supplied by the
inspected source. Prior repo packets already distinguish that consequence from
an exhibited finite witness.

Primary locations:
- Introduction, lines 62 onward and 104 onward:
  https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/build/sections/introduction.tex
- Scope document:
  https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/docs/158.md

## The definition-baking allegation does not match this revision

Both the challenge comparator and the actual solution define proper colouring
by the ordinary unit-edge condition: whenever two complex points have distance
one, their colours differ. There is no extra clause assuming five-colour
impossibility or the analytic transfer.

- Actual definition at `Coloring.lean:7--8`:
  https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/OAI/Geometry/PlaneColoring/Coloring.lean#L7-L8
- Challenge comparator:
  https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/ComparatorChallenges/EuclideanFiveColor.lean

The challenge ends with `sorry`; it is a comparison target, not the proof module.
That placeholder must not be described as the solution's proof. Conversely, the
fact that an actual proof module exists is not evidence that this session built it.

The JSON explicitly lists `ProperColoring` among definitions to compare, points
to solution module `OAI.Geometry.PlaneColoring.Five`, and permits only `propext`,
`Quot.sound`, and `Classical.choice`. `enable_nanoda` is false. These are declared
validation settings, not a successful validation result from this session.

https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/ComparatorChallenges/EuclideanFiveColor.json

The actual final theorem at `Five.lean:124--128` takes a proper five-colouring,
applies `proper_to_borel_weak`, and then applies
`no_measurable_weak_five_coloring`. The latter constructs graph-valued
approximants and invokes a planar cover/mesh contradiction. The main import
file `lean/OAI.lean:479` imports this solution module.

https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/OAI/Geometry/PlaneColoring/Five.lean#L105-L128

No `lean` or `lake` executable was found on the current PATH. No build, full
transitive axiom audit, comparator execution, or independent kernel check was
performed. Therefore this packet records **SOURCE CLAIM / NOT BUILD-VERIFIED**,
not PROVED, and does not turn an apparently outdated criticism into a refutation.
Whether any historical revision had a different definition was not investigated.

## The analytic load-bearing interface

The forward transfer uses the algebraic plane E and its algebraic unit rotations
K. Their semidirect product is amenable. Averaging gives an invariant law on
proper labelings of E. This step alone does not make the translation
representation strongly continuous for the Euclidean topology.

The source next identifies the continuous-vector subspace as a probability
factor, so orthogonal projection is conditional expectation and preserves the
nonnegative label coordinates and their sum. To preserve zero unit correlations
it invokes its central Haar-rigidity statement: a K-invariant probability on the
discrete character dual of E that assigns zero mass to the ordinary continuous
characters is Haar. The residual coefficients then vanish off zero. Strong
continuity, jointly measurable field representatives, Tonelli, and a positive
coordinate selection produce a weak measurable colouring.

The relevant locations are `transfer.tex:145--258`, especially the explicit
rigidity invocation around line 185, and the formal declaration
`BorelTransfer.lean:436`:

https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/build/sections/transfer.tex

https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/OAI/Geometry/PlaneColoring/BorelTransfer.lean#L436-L440

The inspected rigidity source develops a compact-factor tower, conditional
multiple averages, a radial difference inequality, a general-position
number-field basis selection, and a bounded-residue argument. This session does
not certify that entire argument. Its presence does show that the source is not
simply replacing 'amenable' by 'measurable' without an intervening theorem.

For the converse, the source proves density-point classes are unit-independent,
translates each finite configuration into their conull union, and then uses graph
compactness. This is a separately recognizable density/finite-configuration
mechanism; it should not be confused with the difficult forward direction.

## Exact unjustified transfers to avoid in our research

These are invalid shortcuts, not alleged newly located errors in the source:

1. **Measurable exclusion ⇒ unrestricted exclusion** without a valid forward
   transfer. Compactness in the opposite direction does not repair this.
2. **Amenability ⇒ ordinary measurable sample.** The usual-topology continuous
   factor and preserved unit correlations are additional obligations.
3. **Invariant radial PSD kernel ⇒ an actual colouring law.** The earlier q17
   finite cut-polytope hostile refutes this implication.
4. **Finite residue/toroidal graph not five-colourable ⇒ Euclidean plane not
   five-colourable.** Reduction/quotient maps have the wrong direction for this
   inference; G2 gives an exact K7 quotient of a three-colourable lattice.
5. **A quotient closed walk ⇒ a planar closed walk.** G3's nonzero winding and
   G4's acyclic carousel lift give explicit witnesses.
6. **Every finite tested patch succeeds ⇒ all patches succeed.** Compactness
   requires every finite subgraph, not a bounded census. Conversely a genuine
   finite obstruction is decisive only with exact realizable unit edges.
7. **An existential finite obstruction ⇒ a practical finite certificate.** No
   size, coordinate-degree, or running-time bound follows just from compactness.

## Remaining verification that would change status

A genuine promotion of the formal result requires a pinned clean build of the
solution and comparator, full transitive axiom audit, agreement of compared
statements and definitions, root-import reach, and explicit review of trusted
code. A finite real-algebraic coordinate set plus independently checked
five-colour refutation would establish the lower bound without relying on this
analytic chain. Neither deliverable is supplied by this packet.
