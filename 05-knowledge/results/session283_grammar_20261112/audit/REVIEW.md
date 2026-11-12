# Session283 independent review

**Status: PASS WITH MINOR SCOPE/PRESENTATION REPAIRS (applied).**

The new mathematical claims reviewed here are sound within their stated
scopes. This review does not re-audit every inherited theorem, the historical
prime/cycle-index census, or the whole Chan–Pak program.

## 1. Paths: one positive physical code

Independently derived the triangle-cap coefficients by enumerating cyclic
quotient run words; these agree with the inherited formula. The union-bound
surjection estimate, factorial digit bounds, superincreasing inequality,
greedy decoder and reverse-lex ordering argument are valid at every stated
order. All blocks are nonempty, so the physical context is strong. The
polynomial vertex/bit bounds do not assert efficient recovery of an arbitrary
input's path-cover profile.

Independent controls use matrix Hamilton-path DP and permutation/cut profiles,
not the producer's subset path-cover recursion. They check 150 literal small-cap
counts and seven distinct-profile code round trips. The two H9 witnesses have
identical complete profiles but the displayed different full extension menus;
64 attachment counts were recomputed independently. Hence the intended loss
boundary is real, not merely a labelled-cut artefact.

No claim that a decoded formal digit vector is tournament-realizable, or that
one scalar code is a complete isomorphism/growth invariant, is accepted.

## 2. Weighted-core grammar

Reviewed the forest-rank identity as an inherited dependency. The sharp
four-rank-change proof is valid: only homogeneous edges incident to the
reversed arc's endpoints can change, and the mutual endpoint pair does not.
The split-template presentation retains every positive chain piece, which
justifies its exact-object rewrite; dropping lengths preserves rank only.
The internal-chain polynomial and endpoint-type compiler count every arc
once. Nonselected block contraction cannot create a false endpoint
homogeneity because any affected third-block separator belongs to the other
selected block, whose endpoint type and nonempty sides remain represented.

The all-non-singleton bulk polynomial is valid: a distinct-core pair has at
least two separating vertices, so one flip cannot create a new homogeneous
pair across old chains. Cross-chain flips therefore destroy exactly the old
homogeneous edges at their endpoints. The entire row is shape-independent in
this sector, not a multi-step Markov quotient.

Independent controls recompute components of the homogeneous graph directly,
not merely n-e0. They cover all 1,099 labelled tournaments through order five,
10,650 arc reversals, and 27 weighted triangle families for the bulk polynomial.
The two non-lumpability examples and the one-step versus arbitrary-step
qualification are necessary and correctly retained.

## 3. Separator heat flow

The generator is correct: for a fixed vertex pair, each exterior separator
bit is toggled by exactly two possible reversed arcs. Expectation closes even
though the bits for different vertex pairs are correlated. Independent full
trajectory distributions verify 288 spectral moments at orders three and four.
The Krawtchouk eigenvalues, Fourier inversion, stationary expected rank and
Poissonization all have the stated normalization. This is uniform labelled
arc dynamics, not uniform sampling of ordinary isomorphism classes.

The bounded-core scaling proof uses e0=n-b and e1=O(n), and controls all
higher coefficients by Mq²/(1-q). It proves convergence of expected normalized
rank. The initial wording suggesting a typical trajectory was narrowed to
expectation; no concentration result is claimed. Derivative notation was
also separated from the random tournament T'. These are minor scope and
presentation repairs, not failed identities.

## 4. Exact A000568 correction packets

The complementary elimination identity is correct because residual and
backbone cycles have disjoint permitted lengths. All positive correction
defects reside within the residual core; backbone completions retain the
full gcd field, so their moment cell and elementary polynomial are unchanged.
The inherited projective gauge is compatible with those cells. The new
factorial-scaled integer recursion is valid: multiplying by N! converts each
short-cycle choice into the integral factor N!/[s^m m!(N-ms)!]. Final raw
integer, raw rational and projective rational controls agree.

An independent direct odd-partition implementation reconstructed all cells
for 100 cases (orders0..24 and four backbones), plus orders64 and70. This is
stronger than comparing interval endpoints. A separately written affine
moment consumer recovers the 10,303-bit order150 integer from the B35 packet
and agrees with the separately produced modular completion. The B35 packet
was not given that integer as input. A freshly generated order200/B49 packet
also gives the same 18,655-bit integer under this independent affine consumer.
Order250/B59 recovery is covered by the producer normal/optimized replay and
its comparison with the inherited counter, not an additional independent
order250 recount by this reviewer.

The threshold proof is sound: a constructive unique-minimum defect witness
occurs by4a+28; since a>=29, no five-part core can precede it. Repeated parts
cannot leave five positive shifted defects. The remaining four distinct
lengths lie in exactly53 offset patterns. Their gcd comparisons are periodic
modulo6930. Independent enumeration of all3465 odd classes reproduces the
complete offset histogram. The fifth finite-difference perturbation annihilates
moments0..4 and changes the zero atom; strict positivity of the core's support
makes this a genuine rational ambiguity. It is not necessarily integer
rounding ambiguity at that same order.

## 5. Positive contraction decks and their exact limitation

The TT3/C3 profile difference is twice the singleton profile, so multilinearity
proves the displayed Boolean contraction-deck expansion exactly. The mod-2^r
truncation, mod8 criterion and positive-interaction range obstruction follow
without approximation. The range argument excludes a complete odd interval
from the cube's own minimum to maximum only; it does not exclude internal
subintervals, unions of cubes or compensated edits. This is correctly stated.

Independent literal adjacency-matrix checks cover 74 quotient families and
237 switch outputs, including the twelve-vertex mixed-residue example. They
also check the displayed counterexample to the discarded multiplicative
mod8 shortcut. No producer path-counting implementation is imported in these
controls. The full path-cover identity follows from the inherited multilinear
kernel; that historical general kernel is not re-audited in full here.

## 6. Limits and verdict boundary

The session does not prove a finite full core alphabet, faster-than-existing
full exact A000568 counting, or Chan–Pak. It improves a rich certificate
producer and classifies several exact interfaces and their evolution.
A000568 counts already known through the inherited cache are recovered, not
claimed as first numerical discoveries. Runtime comparisons must retain their
workload and machine qualifications, including when the old full counter is
faster than this richer producer.

## Independent reproduction

Run from this directory:

```
python3 independent.py > independent.json
python3 -O independent.py > independent.optimized.json
python3 enumeration_independent.py > enumeration_independent.json
python3 -O enumeration_independent.py > enumeration_independent.optimized.json
python3 switch_independent.py > switch_independent.json
python3 -O switch_independent.py > switch_independent.optimized.json
```

All three normal/optimized receipt pairs agree exactly. `SHA256.json` records
the reviewed producer source and proof snapshot and the independent controls.
It excludes this report and the digest file itself. Only the listed files are
covered by that snapshot; absence of a historical dependency from it is not a
retraction or a claim that the entire inherited repository was reviewed.

Novelty review is scoped: generic substitution, Burnside enumeration,
homogeneous-core decomposition, projective gauges and earlier ground-moment
inequalities are inherited. The packet producer, sharp periodic obstruction
classification, weighted bulk jump law and particular physical code/deck
consequences are checked here against their stated dependencies, not certified
as historically unprecedented throughout mathematics.
