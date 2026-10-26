# Independent audit — Session66 presentation-count scope appendix

**Verdict: PASS.** The sealed `session66-tc-scope/REPORT.md` correctly
distinguishes raw Hamiltonian-word counts from automorphism-orbit counts.
Its elementary free-action and graph-kernel proofs are valid. No correction
is requested. This is a bounded audit, not a novelty claim, new spectrum
theorem, source-repository edit, or publication.

Audited producer report SHA-256:
`4b0e7acc63ee313ded19ada74e0ed3025ebcfbb49e2ce8d94b5baa238ae56c22`.
Exact producer-file and source checks are recorded in `SOURCE_REVIEW.json`.

## 1. Proof and quantifier checks

**Hamiltonian words.** For a finite digraph, an automorphism sends a
Hamiltonian vertex word to another such word. Fixing a word fixes each of
its entries and therefore every vertex. The action is free, including the
empty-word convention, so `tc=H/|Aut(T)|` is its orbit count. This is the
inherited elementary LEM-003, not a new tournament-specific theorem. The
directed-cycle caveat is essential: C3 has one cycle modulo rotation, fixed
by all three automorphisms.

**Joint 7n codec.** The frame tree's unique degree-three center and unequal
arms of lengths 1, 2, 3 fix every internal vertex. A joint automorphism is
therefore determined by its component permutation; preserving the second
layer is exactly preserving tournament arcs. Thus
`Aut(G1,G2) ≅ Aut(T)`. The component-word decoder intertwines these actions,
so it preserves the actual word-orbit invariant. Both the common carrier
identification and fixed order of the two layers matter.

**Gamma kernel.** Degree-one leaves are distinguished from degree-at-least-
three cores. The cores' numbers of leaf neighbors, respectively 1, 2, 3,
distinguish centers, tails and heads, while each center identifies its two
ports. The induced action on centers is exactly `Aut(T)`. Its kernel fixes
every core and independently permutes the two tail leaves and three head
leaves at each vertex, giving `(S2 × S3)^n`, of order `12^n`. Choosing leaf
slots gives a splitting. Hence the correct numerical transport is

```
tc(T) = I(Omega(T),2) / |image(Aut(Gamma(T)) on centers)|
      = 12^n I(Omega(T),2) / |Aut(Gamma(T))|.
```

The kernel fixes every decoded word, so the full Gamma group does **not**
act freely on those words when `n>0`. At `n=1`, `H=tc=1` and
`|Aut Gamma|=12`; the naive quotient is `1/12`. At C3 it is `1/1728` rather
than one. Neither quotient concerns Hamiltonian paths of the gadget graph.
The inherited OCF is used only as a cardinality identity; no equivariant
independent-set/path bijection or division by `|Aut Omega|` is asserted.

## 2. What remains open or unverified

The source lineage is correctly separated:

- LEM-003, `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`,
  supplies free action.
- THM-4741,
  `01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md`,
  distinguishes the proved ordinary-`H` gaps from finite `tc` evidence.
- HYP-9352(2),
  `05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md`,
  leaves the all-order non-rigid `tc` exclusions **OPEN**.

Excluding `H=7,21` does not exclude `H=7|Aut(T)|,21|Aut(T)|` when the
automorphism group is nontrivial. A faithful recoding cannot close this
gap by itself. The appendix explicitly retains this boundary.

The six-family order-ten JSON table is faithfully described as **REPORTED,
NOT INDEPENDENTLY REPLAYED OR CERTIFIED**. Neither its coverage reduction
nor the through-nine census was re-established here. Candidate counters
are not disjoint class counts; the source loop checks both targets
separately. No order-ten source script was executed in this audit. No
all-order `tc` exclusion or new canon status is approved by this audit.

## 3. Independent finite verification

`check.py` was independently authored without importing producer code.
It has no input files, Git, network, or non-standard-library dependency.
Its generic graph routine enumerates actual vertex automorphisms by
color refinement and backtracking, **without false-twin compression or
codec decoding**. It is separately checked against brute permutation
search on all 76 simple graphs of orders 0 through 4. The seven-vertex
frame's asymmetry is also checked by all `7!` permutations.

The normal and `python3 -O` outputs are byte-identical:

| Bounded test | Exact result |
|---|---:|
| All labelled tournaments of orders 0–4 | 76 |
| Explicit Hamiltonian-word orbits | 172 |
| Permutations checked for joint-lift equivalence | 1,590 |
| Full joint automorphism sets compared with unique tournament lifts | 76 |
| Literal odd-cycle/independent-set evaluation versus word count | 76 |
| Direct full Gamma automorphism counts | empty 1; K1 12; TT2 144; C3 5,184 |
| Gamma automorphism multiplicity over each center action | 1; 12; 144; 1,728 |

All gates use explicit exceptions and remain active under `-O`. The
negative control changes the expected singleton quotient to `1/13`;
both modes reject it with the intended explicit error, as recorded in
`NEGATIVE_CONTROL.json`. The
all-order conclusions in Section 1 rest on proofs, not extrapolation from
these small tests. The producer's own normal and optimized checkers were
also rerun; both are byte-identical to its sealed result files.

## 4. Integrity and portable replay

All 15 producer allowlisted files were hashed, and every producer seal
entry passed. All eight copied sources matched their pinned SHA-256.
The six repository sources also matched their calculated Git blob IDs
and `git show` bytes at
`67e6ca529b9ca594a7ddc042abe9800a6fdb0276`. The two graph-packet copies
matched the sealed originals. These are audit-time local-object integrity
checks, not a new independent remote-authenticity assertion and not a
certificate of the source computations' mathematical correctness.

The assigned read-only source boundary was preserved; no startup/fetch,
checkout edit, commit, push, publication, or canon update was performed.
Only this audit directory was written. Replay from any directory with
Python 3.9 or newer:

```
python3 /path/to/check.py > /tmp/tc-audit-normal.json
python3 -O /path/to/check.py > /tmp/tc-audit-optimized.json
cmp /tmp/tc-audit-normal.json /tmp/tc-audit-optimized.json
```

`ALLOWLIST.txt` is the exact delivery set. `SHA256SUMS` seals all listed
files other than itself.
