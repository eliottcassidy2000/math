# Session137 independent proof-frontier audit

**Verdict: ACCEPT the reviewed claims in their stated scopes.** The global
determinant closure has its own portable first-checkpoint report,
[`GLOBAL_DETERMINANT_AUDIT.md`](GLOBAL_DETERMINANT_AUDIT.md). This file records
the five concurrent lanes and the independent finite replays. No global odd
value-spectrum realization or LRC(14) closure is asserted.

Read-only baseline: `bb7cb2d1f61a14632389814469317131d897b640`. The startup
packet and guardrails were read. The source worktree remained clean; this
reviewer did not fetch, commit, push, or edit it. Message-processor inspection
could not run because this isolated checkout had no `.machine-id`; no identity
file was fabricated. Coordination used direct parent/producer messages.

Reviewed source files and exact hashes are in `SOURCE_PINS.json`. A hash pins
the version reviewed, not any later unreviewed producer edit. All independent
mathematical kernels below avoid producer imports. The certificate-verifier
audit intentionally imports that producer implementation as the program under
test, while computing expected gap counts independently. The explicit
adjacency witness snapshot `width_witnesses.json` is producer-derived data,
carries its own source-output hash, and is independently evaluated.

## 1. Substitution-width and residue contexts

Reviewed `width/PROOF.md`.

* The free action on internal return excursions proves the universal divisor
  `d_c=c!(c-1)!`; nonempty disjoint labelled return segments make the action
  free. The upper-triangular matching detector has unique diagonal matching,
  leading coefficient exactly d_c, and the optimal exterior budget 2c-2.
* Triangular induction proves the exact contextual congruence, without division
  by modular nonunits. The minimal exterior detecting size is conditional on
  an actual pair with that first weighted difference.
* The independent-copy substitution kernel is correctly typed as an oriented
  graph, not a tournament. Excursion divisibility and detector expansion still
  hold there. Divisibility by each input d_c separately is enough for modular
  representative independence; a product-divisibility claim is not needed.
* Reversal reciprocity, parity of set partitions, and the stated even-order
  OCF congruence justify the all-tournament upper state counts 5,3,14,11 at
  moduli 3,4,8,12. Actual witnesses establish the lower bounds.
* The seven-state mod12 closure is an exact all-depth statement for its
  specified J/C grammar, not for every width-three construction. The infinite
  branching-resource separation and its width-two boundary are sound.
* The late inactive-repeat lemma is valid: each return through an inactive
  repeated color must contain a globally unique singleton color, so return
  permutations remain free after other colors are collapsed. Only lcm/max
  divisibility is claimed, not a product of the factorials.

Independent replays: 1,099 labelled quotient contexts through order five;
coefficient gcds 1,2,12; sharp detectors through c=4. Separate unordered
set-partition profiles verify every supplied state witness, all 532 tournament
classes through order seven, all four upper residue restrictions, and the
seven-state grammar closure. Both programs match normal and optimized runs.

Files: `audit_width.py`, `width.normal.json`, `width.optimized.json`,
`audit_width_states.py`, `width_witnesses.json`, `width_states.normal.json`,
`width_states.optimized.json`.

The later `WIDTH_UNIVERSALITY_AUDIT.md` accepts the full two-slot-triangle
grammar's all-modulus scalar local universality, with strong width-three
witnesses and arbitrarily large fibres. Independent controls include 121
literal factor pairs, 3,562 simultaneous landings, and 132 complete positive
polynomial-function recipes. This does not close exact integer realization.

## 2. Cross-parent coverage and physical-move barriers

Reviewed `coverage/RESULT.md`, `PHASE2.md`, and `SYMMETRY_EDIT.md`.

* The strong-parent overlap graph has the claimed radius n-4 about the
  transitive blowup of the unique strong four-vertex class. Replacing deleted
  outer vertices by clones leaves a strong common deletion core at each step.
  Its scalar-hull gap corollary does not assert solidity.
* Both inter-ear seam derivations have the correct endpoint orientation and
  signs. The gcd extraction from singleton/pair/triple cuts needs m>=4 exactly
  as stated. Positive nested cuts do not imply feasible unit steps.
* The arbitrary-parent ear criterion uses the source and sink SCCs, not
  nonconstant cuts alone. This restores the explicit H25/H27 seam but does
  not repair the later all-order full-star obstruction.
* Common deletion cards exactly characterize arbitrary one-vertex replacement
  on isomorphism classes. The explicit H147 orbit is closed under zero/unit
  star moves, because every equal-response successor is in the same orbit.
* The r-vertex edit divisibility argument is a valid cellwise free action,
  not an assertion that the altered tournament retains those automorphisms.
  Every cell fixes R and its at most 2r old path-neighbors. At least k-3r
  triangle blocks remain rotatable. The explicit warning about possible
  zero-response escapes is essential and retained.

Independent finite controls:

* all 33,867 labelled tournaments through order six for Hamilton counts and
  strongness; direct permutation/OCF cross-checks on all 1,099 through five;
* 16,932 inter-ear gap identities through base order four, including singleton
  bases; 4,323 strict nested positive-capacity floors;
* all 22,320 labelled strong order-six children: exactly five strong-parent
  unit components of sizes 4320,4320,7440,3600,2640, with the stated H images;
* exactly 1,440 H25-to-H27 two-arc bridges; every bridge has a common pivot,
  same-direction cut change, and nonstrong deletion parent;
* all 353 strong order-seven classes and all deletion cards: exactly eleven
  full-star unit components; all 448 literal star presentations of the H147
  witness, with eleven equal-H successors and no H145/H149 successor;
* the H147 witness is rigid, all seven deletion cards are strong, their H
  values are 43,29,29,41,43,43,43, and each response gcd is two;
* arbitrary arc-flip distances one/two/three reach H145 in exactly 0/0/4 ways.

The large substituted-family zero-response censuses in the producer's
`SYMMETRY_EDIT.md` were not independently rerun here; their all-order
divisibility theorem was audited mathematically. Do not relabel that finite
producer computation as this reviewer's complete independent census.

The later `COVERAGE_CONTINUATION_AUDIT.md` accepts the family-wide canonical
card completion theorem, independently checks 16,215 cards and 18,828 literal
completion cuts, reproduces one full order-eighteen row without orbit
compression, and independently verifies the formal-capacity firewall. These
selected additions do not claim to replicate the large complete censuses.

Files: `audit_tournaments.py`, `normal.json`, `optimized.json`,
`audit_bridge.py`, `bridge.*.json`, `audit_fullstar.py`, `fullstar.*.json`.

## 3. Exact integer constructors and source geometry

Reviewed `constructors/THEORY.md`.

* Four actual sources generate Z4 by integer differences. The two transitive
  source rays approach the missing coordinate axes, proving the full closed
  cone. Axis multiples 3^k and vanishing-weight spikes prove the stated
  closed convex hull, not realization of every lattice point inside it.
* Universal linear actual-source maps therefore have nonnegative integer
  entries and the stated height parity. This closes exactly the inherited
  finite-bank positivity caveat for source-universal fresh-block operations.
* The fixed-height L=g=0 face classification is sound. A nonempty cut is a
  singleton; the selected vertex cannot lie in the last SCC. Ordered
  two-path covers multiply over ordinal SCCs. Removing unselected suffix
  singletons removes factors of two; at most two singleton SCCs need remain.
  The elementary bound of h+2 parent vertices gives a finite complete base
  universe. The reverse sink-padding construction is literal and exact.
* The original-port collapse holds through arbitrary allowed re-marking:
  original A stays in the parent and dominates every new vertex. If any new
  vertex beats the original port, the nonnegative q coefficient is positive
  and terminal F exceeds q. Otherwise the history seals to an ordinal product.
* The overflow tag is an exact capped terminal-output congruence, not an
  assumption that current q is monotone. A sealed fresh re-mark depends only
  on original F; its empty local-parent unit is algebraic, with a literal
  actual global realization. The stated inverse-prefix and affine-digit
  corollaries follow within this operation class only.

Independent finite controls: all 1,098 marked sources through parent order
four, literal path/ordered-cover counts, and 210 local recipes through local
order four, reconstructed on the four unimodular witnesses. Their 38 distinct
matrices are entrywise nonnegative. Separate complete strong-core enumeration
and literal face reductions confirm the minimal dyadic banks
`B1={1}`, `B3={3}`, `B5={5,6,7,9,16}`. The formal state `(5,0,0,1)` is a
source-height hostile even though its F=11 is globally attainable.

Files: `audit_constructor_cone.py`, `constructor_cone.*.json`,
`audit_budget_and_faces.py`, `budget_faces.*.json`.

**Later accepted continuation:** `constructors/BOUNDED_INTERIOR.md` extends
the normal form to every fixed (h,L), allowing base g>0. Its independent
proof review and 4,054-source, all-class-through-six state replay are in
`BOUNDED_INTERIOR_AUDIT.md` and `bounded_interior.*.json`. It proves an exact
finite-base affine-dyadic inverse and O(log X) output growth for every finite
height/interior bank, not only the zero-interior face.

## 4. Natural cycle-packing carriers

Reviewed `cycle_carriers/THEORY.md`; Moon's theorem is a **cited** input,
not independently reproved by this reviewer.

* The parity-sensitive integer inversion of Moon's floor gives the claimed
  source-order bound. The inherited SCC correspondence makes recognition
  componentwise and exact; invisible singleton SCCs must first be discarded.
* Packing k forces source order at least 3k. The consecutive-edge tournament
  has exactly consecutive directed intervals, proving the sharp bound and
  its balanced multi-component extension. The interval generating function
  and coefficient formula follow by ordered gaps and even interval extensions.
* Universal cycles are exactly those with transitive complement. Natural
  antineighborhoods are literal deletion carriers. The link double-count is
  correct, including coefficient one, not two, on N14.
* The complete connected budget table plus ordinal products justifies global
  minimum cycle costs when the minimum lies within the table's budget. It
  does not identify cycle budget with tournament order or substitution width.

Independent finite controls: all classes through order eight, with counts
1,1,2,4,12,56,456,6880 and 6,404 strong classes through eight. Direct
subset/permutation cycle enumeration reproduces all 29 connected carriers
of at most 15 cycle-vertices, and the all-source-order conclusion uses the
proved source bound rather than extrapolation. Independent min-plus closure
gives minimum disjoint-pair counts 2,2,1 at budgets 3,10,14. Independent
multiplicative budget closure has exactly 300 H values through cost 15 and
reproduces the selected global minimum costs, including q(31)=9 and q(35)=13.

Files: `audit_cycles.py`, `cycles.*.json`, `audit_cycles15.py`,
`cycles15.*.json`, `audit_budget_and_faces.py`, `budget_faces.*.json`.

## 5. Harmonic first-missing-speed owners

Reviewed `harmonic_owners/README.md`.

* Circle-gap confinement correctly separates strict gaps from the complete
  primitive equality grid. A nonprimitive q-center violates an existing
  smaller speed. No zero-wrap branch is lost in the owner audit.
* The first missing index d controls the new closed floor; the unit-owner
  set alone does not. Primitive boundary-grid overlap forces owner e=d and
  an integer numerator of absolute value one, giving precisely the two
  reflected overlaps in the coprime case.
* All endpoint owner ties, the harmonic maximum dichotomy, every ceiling,
  and the primitive-grid correction are complete in their stated chamber.
* The fastest-tail owner proof handles the exceptional tail speed q at m=1
  by a residue argument. The mass correction and period formulas retain the
  isolated boundary maxima rather than discard null sets.
* The generic denominator bound for a maximum comes from opposite active
  slopes (including a single function's half-integer cusp); it supports the
  exact positive/negative rational certificate, not an arbitrary sampling
  optimizer.

**Correction caught and repaired:** the draft unpromoted-owner inequality
omitted a factor p. Correct is `(v-p)epsilon<=p(c-1)/q`. The hostile
q=5,d=p=2,c=2,v=9,r=1/7 makes the incorrect sides 2/5 and 1/5. The next
displayed algebra already proved the corrected inequality, so the theorem
was unchanged. The final reviewed manuscript retains the repair/lineage.

Independent interval intersection: every promotion set at q=3..10 and five
rational thresholds, totaling 5,100 cases and 1,020 full boundary component
classifications. A second continuation audit verifies 447 maximum/boundary
decisions, 464 complete ceiling sets, 292 exact terminal-mass slopes, and
428 additional primitive-grid points. The tied owner {2,9} and isolated
{1/5,4/5} hostiles are explicit. All normal/optimized outputs match.

Files: `audit_harmonic.py`, `harmonic.*.json`,
`audit_harmonic_tails.py`, `harmonic_tails.*.json`.

## 6. Replay and manifest

The later counting-complexity consequence of the global determinant proof
is separately accepted in `COUNTING_GAPS_AUDIT.md`. Its implementation-under-
test audit uses independently computed expected gaps and tests all 23,834
single-bit mutations of 184 valid certificates on 78 small/control inputs.
It does not assert arbitrary nonnegative GapP-to-#P closure.

The separate phase-two `WEIGHTED_STABILITY_AUDIT.md` accepts the quantitative
weighted determinant-gap bounds, including their exact pair-deletion form,
and gives an independent 857-matrix rational replay. The first-checkpoint
global determinant packet is unchanged.

The averaging refinement in `WEIGHTED_POLYNOMIAL_AUDIT.md` proves a
polynomial-in-order coefficient asymptotic to 24/n^2, rather than only an
exponentially small variance coefficient. The n=3 identity and n=4 sharp
example are exact; general-order optimality is not claimed. Its separate
856-matrix rational replay checks both root inequalities and coefficient
recurrences, retaining the false n=2 strong-gap extension as a hostile.

Run a Python checker from this directory, redirect output, and compare it
with the corresponding `*.normal.json`; repeat with `python3 -B -O`.
`audit_tournaments.py` uses `normal.json` without a filename prefix.
The scalar discriminant witness checker additionally uses SymPy (version
1.14.0 in this audit); the universal global proof checker does not.

The global C++ reproduction commands and exact modular reconstruction bound
are in the separate determinant audit. All paired runs listed above are
byte-identical. Checks use explicit exceptions; no Python `assert` is
silently removed by optimization. `PUBLIC_FILES.txt` and `SHA256SUMS` identify
the complete public packet, excluding binaries and caches.

All claims remain typed: algebraic proofs, cited classical inputs, complete
finite universes, and finite producer-only auxiliary checks are not merged.
The value-spectrum, all-order scalar interval coverage, and LRC(14) remain
open beyond the explicitly proved results.
