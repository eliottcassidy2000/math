# Source and target status, checked 2026-10-25

**CITED external theorem, not independently formalized in this session.**
OpenAI, *The Euclidean plane is not five-colorable*, 23 September 2026,
Theorem 1.1, states the unrestricted ZFC bound 6 <= chi(R^2) <= 7.
Theorem 1.3 gives a transfer to weak measurable colorings; this is not
an assertion that arbitrary proper colorings become measurable proper ones.
The source was downloaded and its introduction/theorem statements inspected.
Its exact source commit is adc7f1241b42e322a6451854ab7e4b4c146bf78a:
https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf

The linked Lean root, `EuclideanFiveColor.no_proper_five_coloring`, was read:
https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/lean/OAI/Geometry/PlaneColoring/Five.lean
No Lean build, transitive axiom audit, or complete independent proof review
was performed here. Google's Formal Conjectures catalogue independently
records the six lower bound and points to that external formal proof;
its own displayed `sorry` is a catalogue stub, not the cited proof.
https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/ErdosProblems/508.lean

The present session contributes no improved plane bound. Its mathematical
bridges and obstructions below do not depend on adopting the new lower bound.
Where a dichotomy between six and seven is mentioned, the external theorem
is the premise. A seven-color quotient or a residue-coloring obstruction is
not a seven-chromatic Euclidean unit-distance witness.

## Classical triangle/Fano source

Philippe Caldero and Jerome Germoni, *Equilateral Triangles and the Fano
Plane*, American Mathematical Monthly 123(8) (2016), 789-801, supplies the
classical finite-field equilateral-triangle/Fano/Heawood context.
https://doi.org/10.4169/amer.math.monthly.123.8.789
Primary PDF read:
https://kam.mff.cuni.cz/~spring/2017/papers/equilateral_triangles_in_Fano.pdf

The torus and Fano ingredients in this session are classical and independently
reproved here in explicit coordinates. The application selecting the
four-element Campedelli-code orbit from the hourglass arrangement is a
new-to-this-route comparison, not a literature-wide novelty claim.

Do not republish the downloaded source PDFs or full extracted text as part
of the research packet; retain URLs, hashes and our own deductions only.
