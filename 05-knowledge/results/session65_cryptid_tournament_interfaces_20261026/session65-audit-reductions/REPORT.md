# Independent audit — Session 65 typed reduction atlas

**Verdict: PASS, scoped. No unresolved mathematical or execution blocker found.**

Audited on 2026-10-26 by the independent `audit_typed_reduction_atlas` lane.
The reviewed producer is `session65-reductions`; the accepted freeze is
identified below. This is a proof review and finite verification receipt,
not a proof of any headline conjecture, an external formalization build,
or a new completeness theorem for a natural mathematical subfamily.

## 1. Accepted scope and the central conclusion

The report supports a **typed directed atlas of reductions and interfaces**,
not the assertion that all major open problems are mutually equivalent.
It distinguishes fixed statements from uniform languages, and logical
implication from computable or polynomial many-one reduction, simulation,
conjugacy, graph isomorphism, and operation-relative observation equivalence.
The main negative boundaries are correct and useful:

* Two families sharing a compiler into NONHALT do not thereby reduce to
  each other. A common target supplies no reverse arrow.
* An arithmetical formula supplies an upper-bound presentation; it does not
  by itself prove completeness, minimal quantifier level, or a separation
  result for a fixed sentence.
* NP-complete languages are decidable. Finite checking neither establishes
  nor rules out NP-hardness. It does exclude HALT-hardness under total
  computable many-one reduction to a decidable target.
* A common tape/code carrier, an arithmetic address tree, and a truth-preserving
  predicate compiler are three separate constructions.

The explicit NP-decidability warning was requested during this audit and is
present in the frozen report. This was preventive clarification, not a
refutation of the producer's earlier one-way statement.

## 2. Proof review

### Compiler propositions and fixed/uniform typing

Proposition 1's reduction chooses a known yes/no target instance after
running the source decider. It is valid for every decidable source and
nontrivial target; it carries no efficiency claim. Treating one sentence
as a singleton-domain decision problem gives a constant extensional
language, not a uniform undecidability theorem.

Proposition 2 searches finite counterexamples to a uniformly decidable
predicate: false instances map to HALT, true instances to NONHALT. The
assumed total algorithm is explicit, rather than a claim that arbitrary
program totality can be tested.

Proposition 3 sequentially searches witnesses and beeps after each success.
If every search succeeds there are infinitely many beeps; the least failed
search instead causes finite beeping without halting. Applying this compiler
to the complement of a Sigma-2 normal form proves FINBEEP completeness. The
separate NONHALT-to-FINBEEP construction proves FINBEEP cannot uniformly
reduce to HALT: simulate silently, then beep forever exactly if the original
machine halts. This lower-bound argument does not exclude a different Pi-1
formulation of one fixed Collatz sentence. The distinction agrees with
[Aaronson, section 5.10](https://www.scottaaronson.com/papers/bb.pdf).

The full positive Collatz and weak cycle-exclusion rows have the correct
quantifiers and implication direction. Strong Goldbach's prime-pair search
is bounded. Rule 30's non-eventual-periodicity and frequency statements are
kept distinct. The displayed ABC formulation uses reciprocal positive integer
epsilons, which suffice for the usual every-positive-epsilon formulation.
The real/set-coded Mahler, Kakeya and Hopf problems are honestly left outside
this elementary atlas.

### LRC finite checking

The cited theorem really assumes the preceding runner case, primitive
positive integral speeds, and a strict subset-gcd sum bound. Singleton
terms bound every speed in a remaining counterexample. The report discharges
the lower-case hypothesis sequentially and claims a decision procedure for
**through r runners**, not an unsupported isolated-r consequence. The
primitive normalization and replacement of negative relative velocities by
absolute values are legitimate. The individual-tuple interval method is
exact and not polynomial-time in binary speed length as presented.
[Malikiosis–Santos–Schymura, Theorem A and section 1.1](https://arxiv.org/html/2411.06903).

This audit independently checked the theorem statement and hypotheses, not
the whole geometric proof. No feasible LRC(14) enumeration was attempted.

### Bounded JC/DC, and the status boundary

For fixed rank and input degree, actual automorphisms have the cited
inverse-degree bounds: polynomial rank n uses d^(n-1); Weyl rank n uses
d^(2n-1). Thus the existence of inverse coefficients can be expressed in
finitely many coefficient equations. Universal quantification over input
coefficients gives an effective first-order ACF0 sentence; quantifier
elimination decides it without assuming the conjecture. PBW normalization
makes Weyl composition and relation checks finite polynomial computations.
[Bavula, Theorems 1.1–1.2](https://arxiv.org/pdf/math/0512215),
[Adjamagbo–van den Essen, Theorem 2](https://arxiv.org/html/math/0608009#S3).

The bounded degree argument also removes the apparent uncountable-coefficient
coding gap: ACF0 truth is independent of choosing C or its algebraic
subfield. Unbounded degree is still a universal quantifier over these
decidable sentences. No global degree cutoff is supplied.

The historical theorem chain is exactly JC_(2n) => DC_n => JC_n, hence
all-dimensional equivalence. It is not fixed-rank iff, a Poisson-shadow
identification, a complexity reduction, or a dynamical conjugacy.
[Belov-Kanel–Kontsevich, introduction and Theorem 1](https://arxiv.org/html/math/0512171).

The report does not replace the repository's newer higher-dimensional JC
record with stale background knowledge, and does not independently promote
[Gao's 2026 preprint](https://arxiv.org/abs/2608.00222) to an audited theorem.
This audit inspected its metadata/claimed classical definition, not its
counterexamples or proofs. The repo's THM-1300 certificates and live JC(2),
DC(1), DC(2) status were not re-audited.

### RH, Diophantine solvability, and machine-source labels

Lagarias's reverse argument gives a strict inequality violation if RH is
false. Certified rational upper intervals therefore suffice to semidecide
that failure; no decision of equality of arbitrary computable reals is
needed. The producer uses precisely this safe direction.
[Lagarias, Theorem 1.1 and reverse proof in section 3](https://arxiv.org/pdf/math/0008177).

The general Diophantine completeness benchmark is legitimate and does not
transfer to a fixed shared-root variety. The source's equation (3) represents
r.e. subsets by a polynomial with nonnegative integer witnesses; unrestricted
integer witnesses can be obtained by replacing each nonnegative witness by
a sum of four integer squares. Membership is semidecidable by enumeration.
The report's cited primary author source was independently opened here.
[Matiyasevich–Robinson, section 3, equation (3)](https://arxiv.org/pdf/0802.1052).

Cryptids verification labels are accurately attributed as source-reported,
not as fresh builds. Full-Collatz beeping, weak-Collatz ordinary halting,
and the verification caveats are not conflated.
[Cryptids index](https://wiki.bbchallenge.org/wiki/Cryptids).

### Graphs, symmetry, incidence codec

Proposition 4 is the arc-count obstruction N(N-1)/2 <= Nk. With k=1,
a deterministic configuration graph on at least four vertices cannot
literally be a tournament. A halt vertex sharpens the necessary bound to
N<=2k. Neither condition is asserted sufficient.

Proposition 5 correctly classifies a single rooted deterministic orbit as
a terminating path, a tail entering a cycle, or a ray. Nonrepeating infinite
runs all have the same unlabelled ray, while time-indexed infinite runs
are rays even when configurations repeat. This makes bare orbit isomorphism
a nearly vacuous classification of the interesting nonhalting runs.

Proposition 6 correctly excludes a **faithful symmetry-preserving functor**
from all stated control graphs to finite tournaments: an order-two
source automorphism cannot inject into an odd-order tournament automorphism
group. Colored tournaments have no exception. The report explicitly does
not claim this refutes an arbitrary graph-isomorphism many-one reduction,
or the ordered/reachable deterministic subcategory.

Proposition 7 is a valid n+2m ordered incidence codec. The three cross-class
pair rules are disjoint; within-class order fills every remaining pair.
Unique incidences decode loops and repeated endpoints, and a port payload
retains read/write/move labels. Preserving state/start/output decorations is
a straightforward additional carried label in the proof, not something
exhaustively exercised by the producer's simpler edge-table test input.
Dense filler arcs are not legal computation steps. No native K move,
Hamiltonian predicate, or whole computation-graph equivalence is claimed.

### Native continuations and black core

The recovered THM-3121 substitution formula and Session 33 context-budget
converse retain path-cover profiles, not merely present Hamiltonian counts.
The factor-return budget and triangular nonzero diagonal argument are sound.
The sparse-port and deterministic-versus-signed-linear distinctions are
preserved. The general continuation quotient is an output-preserving right
congruence; the report correctly requires two-sided contexts for a two-sided
composition claim and gives no infinite-TM minimization algorithm.

New section 4.0 fixes a finite native target: tournament classes, path-free
arc complement K_P, path reversal as the anti-automorphism color test,
and non-pure-blue eligible vertices. It matches the pinned THM-643 and
HYP-9346 conventions. Converse merging is explicitly a separate quotient
target. A quotient path is not assumed to lift, and pure-blue nodes are
not accidentally inserted into the black core.

Session 64's valuation-word correspondence is kept on its punctured 2-adic
domain. Constant word 2 sends Collatz 1 to twisted 3/5, so positive-integer
realizability already fails. The native Antihydra guard/debt is retained
as an additional obstruction. The three-carrier summary keeps this apart
from the sibling finite tape code and uniform predicate compiler.

## 3. Exact finite checks

The producer script was fully read before execution. Its normal and `-O`
runs reproduce both frozen JSON files byte-for-byte. Every gate uses an
explicit runtime exception, not an assertion disabled by optimization.

The independent `check.py` in this directory additionally verifies:

| Test | Exact universe/result |
|---|---|
| Cross-codec tests | 25 one-state and 6,561 two-state partial deterministic tables, varying target, write bit and move independently in every read slot |
| Independent implementations | Tagged-vertex/set-arc encoder and incidence decoder, cross-tested in both directions with producer functions |
| Mutation controls | All 45 single-arc corruptions of the canonical 10-vertex sample rejected; changed payload remains a valid code and decodes to changed source data |
| Continuation hostile | Independent permutation-and-cut profiles (15,78,198,240,120) and (15,90,210,240,120); direct seven-vertex counts 123 and 135 |
| LRC interval controls | 774 tuple/threshold cases, lengths 1–3, speeds 1–6; independent exact common-denominator integer-grid method; 446 nonempty |
| Integrity | Eleven input snapshots match both SHA-256 pins and exact git blobs at the pinned commit; 17 producer manifest rows verified |
| Optimization safety | Auditor and producer contain zero AST Assert nodes; independent normal and optimized JSON are byte-identical |

Reproduce from the session workspace, keeping the producer and auditor
folders as siblings and the `session65-repo` checkout available:

```bash
python3 -B session65-audit-reductions/check.py > session65-audit-reductions/normal.json
python3 -B -O session65-audit-reductions/check.py > session65-audit-reductions/optimized.json
cmp session65-audit-reductions/normal.json session65-audit-reductions/optimized.json
python3 -B session65-reductions/check.py > session65-audit-reductions/producer-normal.json
python3 -B -O session65-reductions/check.py > session65-audit-reductions/producer-optimized.json
cmp session65-audit-reductions/producer-normal.json session65-reductions/normal.json
cmp session65-audit-reductions/producer-optimized.json session65-reductions/optimized.json
```

These are finite controls, not proofs by testing of the unbounded statements.
The independent script deliberately rejects a changed producer freeze.

## 4. Freeze, provenance, and execution boundary

Source commit: `6a8ad5d3bf01a232da4de41bcbff68390a9337ee`.
The current checkout was `fcc112a69d7080db4cd619d6e7927be1992bc51e`;
source pins were not silently moved forward.

| Accepted producer file | SHA-256 |
|---|---|
| REPORT.md | `434cba60cf4760eb9ddb0d9cf25500e54b1db321376eeefc6b4bdd5d946fc898` |
| check.py | `a5bf1c7b398393416d660008b9f7804f77de2ca39dc9592c3beaa80e68626396` |
| normal.json and optimized.json | `153b7925cf69f895b634e94388e3f8f5a8708cad49f06b2f7c501e65c6cf16e3` |
| SOURCE_PINS.json | `1c203205287a884f5eb38091b9e0503b8bdb7ce21823c678950412087993e192` |

External primary documents were read through the web tool; no downloaded
external code was executed, and no external paper bodies were copied into
this audit. Links and concise paraphrases are retained instead. Only the
inspected local producer script and newly authored audit code were executed.
No repository file was edited, staged, committed or pushed by this lane.
The parent owns integration and publication. No permission denial occurred.

**Remaining frontier:** establish useful reductions in both directions
between explicitly chosen natural input families, with initialization,
guards, witness transport and resource bounds. The present audit validates
the proposed map's boundaries; it does not fill those missing arrows.
