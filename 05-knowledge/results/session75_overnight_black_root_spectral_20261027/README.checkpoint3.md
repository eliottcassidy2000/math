# Session75 checkpoint 3 — audited boundaries and finite resource tests

**Current checkpoint: independently audited results in the scopes below.
The overnight session continues; this is not session closure.**

The [C1 spectral packet](README.md) and [C2 core/ROOT packet](README.checkpoint2.md)
are preserved byte-for-byte and mode-for-mode. This addition was prepared
against incoming main `86d49e6bf47` plus the semantic content of local C2
`edebe48cd565`; publication commit identities are supplied by the publication
owner, not invented here.

## 1. Black exit and core targets are different

**PROVED and independently audited:** explicit strong, black-active examples
at every order `n >= 10` have `lambda=1`, minimum semidegree `delta=2`, and
exact black distance `2` to the non-strong locus. The even rooted-converse
construction has one balanced minimum shore; the odd extension has two
minimum shores and the same nine blue decreasing words. Separate literal
black routes prove the upper bound. A different, separately proved route
attains one singleton and a genuinely strong deleted core in exactly two moves.

- [Family, finite minimum-order gate, and conventions](session75-black-boundary/RESULT.md)
- [Even rooted-converse proof](session75-black-boundary/family/ROOTED_CONVERSE_COUNTEREXAMPLE.md)
  and [odd middle-vertex extension](session75-black-boundary/family/ODD_MIDDLE_EXTENSION.md)
- [Separate core survivor and unique-unbalanced-minimum-shore theorem](session75-black-boundary/survivors/SURVIVORS.md)
- [Independent adversarial audit](session75-black-boundary-audit/AUDIT.md)
  and [acceptance](session75-black-boundary-audit/ACCEPTANCE.json)

The minimum order ten is **FINITE-EXACT**, authenticated over every ordinary
class through order nine, not a converse-merged sample. Pure-blue `C3` and
`R5` are explicitly separate. Native complement and path-reversal gauges
are related by converse; named source/sink endpoints are not silently equated.
General singleton-plus-strong-core distance `delta` and all-order black-core
connectivity remain **OPEN**. The unique-unbalanced-shore theorem is not a
claim for balanced or nonunique minima.

## 2. Equal H and every named cut do not determine black exit distance

The two fixed order-ten matrices have the same `H=657`, the same named score
vector, and all 1,022 identical named cut capacities. Their automorphism
orders are `1` and `9`; black exit distances are `1` and `2`. Both have core
distance two. Their full odd-cycle-collection polynomials differ; even the
sets of labelled cyclic triples differ despite equal triangle count fourteen.

- [Exact matrices, counts, paths, cycle polynomials and incoming isomorphism](session75-black-twins/RESULT.md)
- [Independent permutation-cover audit](session75-black-twins-audit/AUDIT.md)
  and [final acceptance](session75-black-twins-audit/ACCEPTANCE.json)
- [Root's separate counts, all named cuts and literal core routes](session75-root-twins/normal.json)

**Concurrent credit and temporal scope:** the E10/rooted-converse object agrees
by an explicit orientation-preserving isomorphism with the incoming THM-4811
recipe at `87d61f67673`. No first-discovery claim is made. The
[source-bound odd11 correction](session75-black-twins-audit/C3_CORRECTION_NOTE.md)
refutes the original HYP-9421 even-only classification and original even-family
exclusivity. It does not itself prove or refute the universal `lambda+1` bound.
Later incoming main separately claims the corrected all-order classification
and bound; those stronger claims are preserved, not independently certified
by this frozen C3 family audit.

## 3. Joint depth/clock cutoffs: exact finite coverage, not R3

At precisely nineteen registered scales, all **245,826** depth×integer-clock
breakpoint rectangles satisfy `D+ <= S_N/2`. This is not every integer scale
between 4096 and 262144. The denominator retains all **1,000,908** original
source labels across the separate scale-specific families. The largest
recorded ratio is `4132/44381` at N18275.

The independent auditor rebuilds every source, checks the entire rectangle
universe, recomputes **42,425** rectangles numerically at N4096/N18275/N65536,
and directly checks **43** further cutoffs. It does **not** independently
recompute all 245,826 numerical rectangles. The producer normal/optimized
full replay is separate evidence, not another independent implementation.

- [Definitions, exact pair-birth identity, finite result and exclusions](session75-joint-cutoffs/README.md)
- [Registered coverage](session75-joint-cutoffs/COVERAGE.md)
  and [finite summary](session75-joint-cutoffs/SUMMARY.json)
- [Independent audit with precise coverage](session75-joint-cutoff-audit/AUDIT.md)
  and [acceptance](session75-joint-cutoff-audit/ACCEPTANCE.json)

**PROVED clipping repair:** in the same fixed terminal-labelled tree, with
removed labels retained as zero leaves, put
`delta_k=Delta_k(full)-Delta_k(retained)`. Then
`D+(retained)-D+(full) <= sum_k (-delta_k)_+`.
At authentic N4096, removed norm22 raises debt by56; the complete signed
boundary-pair invoice is82. This refutes payment by removed self-norm alone,
not the half-source target. Sharpness on a generic scalar family is NON-NATIVE.
Both depth-only R3 and its stronger joint-cutoff version remain **OPEN**.
No Collatz, convergence, universal debt, or all-source theorem is inferred.

## 4. Independent all-order tc7/21 exclusion audit

The [analytic/group audit](session75-tc-audit/AUDIT.md) accepts the exclusion
with [explicit proof repairs](session75-tc-audit/PROOF_REPAIRS.md): positive
factors only, an exact integer tail argument rather than false floor-margin
monotonicity, infinite-order recurrence, index-five elimination, and actual
`G`-orbit versus full-automorphism wording.

Its independent acceptance path replaces the original large order12 element-
type scan with six configurations and **45 jobs / 46,968 strong tournaments**.
The [window replay](session75-tc-window-audit/AUDIT.md) separately reproduces
**83 jobs / 185,642 strong tournaments**, with no target, plus 526 Cayley inputs
and analytic/group controls. The window checker retains source algorithms
while independently replacing the automorphism gate and auditing every job;
it is not wholly separate code. The original multi-million order12 scan was
**not rerun by Session75**. Historical n<=11, OCF and classical dependencies
remain inherited. Later THM-4812(g) is outside this C3 acceptance and unchanged;
newer incoming audit evidence remains credited separately.

## 5. Frozen bytes, source bindings and replay

All original native-manifest entries are retained. Large data files alone
use deterministic gzip transport; scripts and proofs remain readable.
See [curation](CURATION.checkpoint3.md), [transport index](TRANSPORT.checkpoint3.json),
[exact source origins](SOURCE_ORIGINS.checkpoint3.json), and
[bindings](PROVENANCE.checkpoint3.json). Forty-four snapshot origins were
compared with targeted `git show` at their declared commits during assembly.
The portable verifier checks those exact bytes/blob identities offline; it
does not claim a fresh GitHub-origin lookup. With the historical Git objects
available, `python3 -B c3_check_git_origins.py --repo /path/to/math`
repeats just these 44 exact source comparisons.

From this directory, using Python3, Bash and a C++17 compiler:

```sh
python3 -B verify_packet.py --manifest FILE_MANIFEST.checkpoint3.json
python3 -B c3_transport.py
python3 -B -O c3_transport.py
python3 -B c3_provenance.py
python3 -B -O c3_provenance.py
python3 -B test_c3_corruption.py
python3 -B -O test_c3_corruption.py
python3 -B test_c3_provenance_corruption.py
python3 -B -O test_c3_provenance_corruption.py
T=$(mktemp -d)
python3 -B c3_transport.py --materialize "$T/native"
python3 -B c3_replay.py --materialized "$T/native" --work "$T/normal" --mode normal
python3 -B -O c3_replay.py --materialized "$T/native" --work "$T/optimized" --mode optimized
python3 -B c3_replay.py --materialized "$T/native" --work "$T/ubsan" --mode ubsan \
  --component boundary --component boundary-audit --component tc --component tc-window
python3 -B check_cpp_rejections.py --normal "$T/normal" --optimized "$T/optimized" --ubsan "$T/ubsan"
python3 -B -O check_cpp_rejections.py --normal "$T/normal" --optimized "$T/optimized" --ubsan "$T/ubsan"
```

Compilation and generated ledgers are outside the immutable packet. The replay
adapter only reroutes an original absolute sibling lookup and four exact
hash-pinned source reads; it does not alter mathematical algorithms. Native
reports, including their historical “audit pending”, “reserved” or “OPEN”
phrases, are unchanged snapshots. The current acceptance and incoming-source
scope above supersede those time-local labels, not their mathematical content.

[C3 verification summary](VERIFICATION.checkpoint3.md) records the actual runs,
coverage and exclusions. Original C1/C2 manifests remain unchanged historical
freezes; the latest whole-packet manifest includes all older files plus C3.
