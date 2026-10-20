# Connected golden basins, outer blue detours, and the raw path spectrum

20 October 2026, continuation e. This is the maintained synthesis.
**Global black-core connectivity and raw-spectrum completeness remain OPEN.**
The statements below distinguish actual path-complement moves, regular
score-preserving moves, and generating functions that count different objects.

## 1. Feedback inventory no longer separates golden terminal basins

[The four-black creation theorem](inventory.md) proves, for any nonempty
source context A,

    A → 1 → 1 → 1 → C3  ↝  A → C3 → C3.

The four moves restore all context arcs and change the endpoint feedback
inventory by one. Their blackness follows from explicit opposite-position
score inequalities, including a nonendpoint pair that eliminates the
apparent exceptional context size three.

Together with the inherited fixed-inventory normalization, this proves:
**all active SCC words over `{1,C3,R4}` at any fixed order lie in one
ambient black component.** They all reach `1^(n-3)→C3`. This is stronger
than fixed-feedback connectivity and does not assert that intermediate
tournaments stay in the SCC language.

[Every strong five-vertex end atom](five_atoms.md) can also be normalized:
`A→S5 ↝ A→1→1→C3` in four black moves, for every nonempty A with a source.
Six finite templates cover the six strong types. Their uniform proof
uses known source arcs and untouched cross arcs; the small base has one
explicit `A=1→R4` hostile, repaired by changing the context path. The
[independent C5 audit](five_atom_audit.md) checks the regular atom separately.
An arbitrary nonempty context on *both* sides is not covered.

## 2. Every nonstrong mixed source has blue-edge detours

[The endpoint-rigidity proof](outer.md) establishes a black detour of at
most five moves around **every blue line incident to a nonstrong mixed
tournament**, with no bound on its outside SCCs.

An endpoint-flexible outer component uses the existing common-hub move.
For an endpoint-rigid component, either a consecutive cyclic triangle
supplies a three-move ledger, or two forbidden prefix rotations force a
near-source vertex. Removing it exposes a rigid SCC word and a smaller
terminal strong component. Recursion selects a terminal three- or
five-move ledger; it does not accumulate five moves at each level.
The all-order rigid-strong alphabet used here is the incoming source-backed
classification `{1,C3}`; it is not a conjectural pure-blue alphabet.

Consequently a genuine mixed-to-mixed blue bridge, if one exists, must
have **two strong endpoints**. This does not prove that the black core
is connected. The older black descent still ends at arbitrary nonstrong
terminals, many of which have larger SCCs than our normalized language.
Here a bridge is tested in the full metagraph, or after contracting its
black components. The blue-only subgraph can have bridges that acquire
black detours in the full graph; it is a different question.

[Combining this with the recovered pure-blue hubs](pure_hubs.md) gives an
exact all-order reduction, for n>=4: **the active black graph is connected
if and only if every strong-to-strong mixed blue line has a black detour**.
This uses Reid's CITED combined connectivity. Merely excluding blue
bridges is weaker: several black components could be joined by a blue cycle.
The same integration note gives an elementary replacement for the
rigid-strong classification input, using the incoming circulant reroutings.

## 3. The sequence has an exact matching interpretation

[The arithmetic bijection](arithmetic.md) maps a golden SCC word of
weight n to a matching of the graph on `1,...,n` with edges `(i,i+2)`.
That graph is two disjoint paths. The map preserves vertex count,
directed-triangle count, and converse as reflection. Hence

    a_(2k)=F_(k+1)^2,       a_(2k+1)=F_(k+1)F_(k+2),
    sum a_n z^n = 1/(1-z-z^3-z^4).

Thus the user's sequence is exactly

    1,2,6,15,40,104,273,... = F_j F_(j+1),
    273 = 13·21 = 1+1+4+9+25+64+169.

The sum of squares comes from a recursive matching partition, not just
an algebraic coincidence. Marking the second edge of the shorter P6 at
order13 splits the matchings into 63 containing it and 210 avoiding it.
Adding its reflected fourth edge gives a reflection-invariant union of
105 matchings: 21 contain both and 84 exactly one. This is an actual
`21+84` partition, but no equivariant identification with the earlier
105-junction atlas is claimed. The restricted converse quotient has
54 types; 105 is its ordinary, unmerged count. Two reflection-fixed words
are pure blue, leaving 103 active ordinary words and 52 active merged types.

The marked primitive pair `(13,21)` gives a neighboring Farey direction,
the torus-surface framing 273 of `T(13,21)`, and the primitive Pythagorean
triple `(136,273,305)`. The preceding pair `(8,13)` gives word count 104
and triple `(105,208,233)`. These are maps of family cardinalities and
marked arithmetic data, not graph isomorphisms or a Petersen-minor theorem.
The limiting Euclid leg ratio is 1:2. Finite piecewise bijections preserve
cardinality, so no finite Banach–Tarski duplication is involved.
The two-step Fibonacci update also has a precise upper-half-plane action:
`z -> (2z+1)/(z+1)`, with determinant one and boundary fixed points
`phi,-1/phi`. It preserves Farey adjacency, but increases the tournament
order by four rather than moving within one metagraph.

The [local interaction refinement](arithmetic_interaction.md) retains the
raw Hamiltonian count too. If chi_i marks edge `(i,i+2)`, put
`E=sum chi_i` and `C=sum chi_i chi_(i+1)`. Then C counts R4 atoms,
`c3=E`, minimum feedback-arc number is `E-C`, and
`H=3^(E-2C)5^C`. A four-state transfer records this adjacent-row interaction;
weighting the two parity paths independently would give the wrong R4 count.

## 4. The forbidden raw counts are a known open conjecture

[The spectrum note](spectrum.md) identifies the exact question as
[Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf):
whether all positive odd integers except 7 and 21 occur as raw H(T).
Both omissions are proved. Completeness is still open; it does not follow
from the metagraph counts or from finite realizations up to 80,557.

Using Busch's strong minimum, deleting singleton SCCs compresses every
H=N witness to at most `floor(3 log_3 N)` vertices. Thus each individual
membership question has a finite logarithmic-order core test, and the
raw counting function is concise. This does not supply every missing
positive construction uniformly.

For each fixed N the ordinary-class series F_N(z) is rational, with
denominators powers of `(1-z)`. The note determines the strong raw-count
polynomials through 23 and gives, for example,

    F_9(z)=2z^5/(1-z)^2 + z^6/(1-z)^3,
    [z^n]F_(3^a5^b)(z) ~ n^(a+b)/(a!b!).

These are raw-H strata, distinct from D_beta, blue maxima M_k, and the
earlier orbit-deficit rational functions. Historical false common-
intersection and conflict-graph sentences were repaired with their
lineage retained in the canon and MISTAKES ledger.

## 5. What regular edge changes can and cannot explain

[The regular-tournament note](regular.md) proves the classical fixed-score
connectivity under directed triangle reversals. At order seven their
isomorphism-class quotient connects

    H=175  —  H=171  —  H=189.

All three have outdegree 3 at every vertex and 14 directed triangles.
Thus H and isomorphism type can change while these quantities remain fixed.

A single arc reversal destroys regularity. So does every actual R_P or
K_P Hamiltonian-path move at nontrivial regular order: its endpoints have
opposite score charges. Regular vertices form an independent set in the
actual metagraph. On fixed labels at odd order, R moves also preserve
orientation parity whereas a triangle reversal changes it; converse or
odd relabelling can erase that parity. No conclusion about merged black
components follows from the labelled obstruction alone.

## 6. Reconciliation with incoming work

The [all-order rigid-strong classification](../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md)
is used in the endpoint-rigid recursion above. It classifies one-path-orbit
strong tournaments; it does not identify rigidity with being pure blue.
The [independent pure-blue closure](../rotation_owners_spectrum_20261020e/purity_closure.md)
now gives exactly the strong atoms `{1,C3,C5}` and the ordinary class series
`(1+z+z^3+z^5)/(1-z^2-z^6)`. Its cutoff and certificate composition was
independently replayed during integration.
The [incoming Reid translation](../paired_black_bridges_order_nine_20261020/BRIDGES_THROUGH_NINE_AND_LOCAL_CYCLE_LAWS.md),
Section 6, gives combined blue-and-black connectivity and a black path
from every nonselfconverse node to some mixed self-converse node.
Neither input identifies all mixed anchors with each other.

The separately audited deficit-two classification gives `D(z)=z^4`.
The [higher-center refinement](deficit_centers.md) uses the incoming cubic
order bound to make every fixed-defect center series D_d a polynomial.
It also supplies two small strong self-converse witnesses proving
`D_6 >= z^5+z^6`: the quartic result does not continue as a monomial ladder.
For the paired-black quotient Q, the incoming complete census excludes
nonleaf bridges through order nine; the all-order nonleaf contribution
remains unknown. This is distinct from the mixed blue-edge question.
The [new maxima packet](../blue_maxima_sqrt_e_20261020/README.md) proves
`liminf M_k/k! >= sqrt(e)` and records exact `M_6=1185, M_7=8295`.
It does not prove a limiting equality or a general formula for M_k.
The [complete wiring certificate](../merged_metagraph_realizability_20261020/REALIZABILITY.txt)
retains all halfedge assignments and endpoint isomorphisms; its exact
realizability criterion is not a replacement for a black connectivity proof.

## 7. Reproduction and remaining targets

All new sources and receipts are colocated. Run a checker from the repo
root, e.g. `python3 05-knowledge/results/spectrum_basins_20261020e/inventory.py`.
Normal and `-O` modes are checked; independent proof reviews and replay
records distinguish source-dependent theorems from finite controls.

| checker | declared finite universe |
|---|---|
|`inventory.py`|1,101 source/path contexts through order7; abstract macrograph through order21|
|`five_atoms_verify.py`|six end atoms, 77 small contexts per atom, repaired hostile, parametric indices|
|`five_atom_audit.py`|independent C5 replay: 1,100 labelled contexts, 7,888 context Hamiltonian paths|
|`outer_verify.py`|31,731 strong fixed-path presentations, each in three middle contexts|
|`outer_audit.py`|separate boundary-ledger and restricted-rigidity checks|
|`arithmetic_probe.py`|inverse matching bijection and reflection through order18; exact recurrences through200|
|`regular_probe.py`|all regular labelled tournaments through7, all36,960 order7 triangle moves|
|`spectrum_probe.py`|all1,100 fixed-path presentations through6; raw-count polynomials and one-apex controls|
|`pure_hubs.py`|117 inherited padded hub routes, three C5 boundary routes, 126 circulant witnesses|
|`arithmetic_interaction.py`|all7,920 words through18; local weights, symbolic and integer transfers|
|`deficit_centers_check.py`|two specified positive defect-six witnesses; no large-defect census|

The next structural targets are now narrower: normalize interior SCCs
beyond R4; connect terminals with arbitrary larger atoms to the golden
component; and analyze blue lines whose two endpoints are strong. No
mixed-to-mixed blue bridge has been found, and neither full black-core
connectivity nor the `{7,21}` completeness conjecture is claimed solved.
