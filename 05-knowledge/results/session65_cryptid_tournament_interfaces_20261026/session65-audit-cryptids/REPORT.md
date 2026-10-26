# Independent audit: three total Turing machines and ordered-join tape encoding

**Verdict: PASS for the stated guarded macros, native clocks, and tournament tape observer.**
No asymptotic morphology theorem, Cryptid halting resolution, autonomous tournament dynamics,
K-move simulation, or equivalence with a full Turing-machine configuration graph is certified.

Audit date: 2026-10-26. Audited repository base:
`6a8ad5d3bf01a232da4de41bcbff68390a9337ee`.
Input packet: sibling `session65-cryptids`, including its post-SCC-decoding data and
the narrowly authorized Coq-inclusion wording amendment recorded below.
This auditor made no repository changes, commits, or pushes. Incoming research reported by
other agents was not silently substituted for these pinned dependencies.

## 1. Identity, scope, and evidence

The exact audited programs are:

| Name | Literal transition string | Macro section |
|---|---|---|
| Spaghetti sandwich | `1RB1LC_1LA1RD_1RA0LC_1LB0RD` | A, B |
| Bouncer + bell | `1RB0LA_0RC1LA_1RD0RA_0LB1RB` | A |
| Typical irregular bell | `1RB0RC_1LC1RA_1RA1LD_0LC0LA` | A, C |

Their eight entries all write a binary symbol, move left/right, and enter A/B/C/D.
Thus their transition maps are total on **every** two-sided binary tape and head/state,
not just finite-support tapes or the blank orbit. Induction on native time gives an
instruction at every finite time. This elementary nonhalting fact does not prove
recurrence, nonperiodicity, envelope growth, or the advertised visual morphology.
The word “Cryptid” in the input directory is not a theorem that these total machines
have a hard unresolved blank-tape halting problem.

The [zoology page](https://wiki.bbchallenge.org/wiki/Non-halting_Turing_machine)
identifies precisely these examples, calls the third a typical irregular bell, and
uses an informal description for spaghetti behavior. Its retrieved-from permalink
identifies oldid 8776. The [Cryptids page](https://wiki.bbchallenge.org/wiki/Cryptids)
uses a different, halting-hardness-oriented definition. These are source-supported
identifications, not independent morphology classifications.

**Independent implementation:** `independent_check.py` imports no producer code.
It parses the literal tables with a separate token validator, uses a dictionary of
both zero and one cells, and advances native instructions until a control/read
**stopping event**. It does not execute for the macro's predicted time and then
assume the endpoint is the first intended return. Its macro specification reads
an immutable original support and constructs a separate final patch. Graphs are
adjacency matrices; SCCs use two-pass DFS rather than the producer's reachability
bitsets. All assertions use explicit exceptions and survive `python -O`.

## 2. All-parameter native proof, independently rederived

Let a tape be a binary function on Z, with head h and state q. Cells not listed
below remain unchanged. Finite support guarantees each maximal run terminates;
for arbitrary infinite tapes the formulas apply whenever the indicated run is
finite. A cell-write counts as an edit only when its symbol changes. Times count
literal machine instructions, including unchanged writes.

### Sandwich

Put d(A)=-1 and d(B)=+1.

* On 0, the literal instruction writes 1, enters the other section state, and
  moves to h-d(q). Time 1, one edit, population change +1.
* On 1, let k>=0 be maximal with ones at h+d,...,h+kd and zero at h+(k+1)d.
  The initial instruction keeps the current one and enters C (q=A) or D (q=B).
  The exact instruction word is `A1 C1^k C0` or `B1 D1^k D0`.
  Each repeated instruction erases one distinct cell and moves one place in d;
  the final instruction writes the terminal zero to one and moves back.
  Therefore the final head is h+kd, state q, time k+2, edits k+1, and population
  change 1-k. The k=0 case still takes two instructions and edits the adjacent cell.

This proves the formula for every k: one more repeated C1/D1 deletes the next
one, adds exactly one native instruction and one displacement, and leaves the
same terminal calculation. The A/B section returns after finite positive time.

### Bouncer + bell

At A on 1, let k>=1 count the maximal leftward run including the current cell.
The word `A1^k` erases those k cells and stops at the first zero in state A,
head h-k, time k, edits k, population change -k.

At A on 0, write one and enter B at h+1. In B, the following **disjoint,
exhaustive** prefix guards describe the native instruction words:

| Guard | Instructions | Tape effect | Exit |
|---|---|---|---|
| `001` | B0 C0 D1 | `001 -> 011` | B, three cells right |
| `1` | B1 | no edit | A, one cell left |
| `01` | B0 C1 | erase the second cell | A, two cells right |
| `000` | B0 C0 D0 B1 | write the middle cell to one | A, same starting cell |

After k>=0 maximal `001` blocks, j=h+1+3k is the start of exactly one terminal
prefix `1`, `01`, or `000`. The initial cell becomes one and each traversed block
becomes `011`. The full terminal data are respectively:

| Terminal guard | Final (head,state) | Time | Edits | Population change |
|---|---|---:|---:|---:|
| `1` | (j-1,A) | 3k+2 | k+1 | k+1 |
| `01` | (j+2,A) | 3k+3 | k+2 | k |
| `000` | (j,A) | 3k+5 | k+2 | k+2 |

Induction simply appends one three-instruction `001` loop on disjoint cells.
Maximality excludes `001` at termination. Finite support precludes infinitely
many such blocks. Weakening terminal `000` to `00` is invalid: it admits the
continuation guard `001` and predicts a false endpoint/clock.

### Irregular bell

At A, the exact instruction words and data are:

| Original local symbols | Instructions | Edit | Exit (head,state) | Time |
|---|---|---|---|---:|
| current 1 | A1 | erase h | (h+1,C) | 1 |
| current 0, right 1 | A0 B1 | write h to one | (h+2,A) | 2 |
| current 0, right 0 | A0 B0 | write h,h+1 to one | (h,C) | 2 |

At C, the two-step word C1 D0 traverses a `01` pair to the left, beginning with
the head on its right-hand one. Both writes preserve the existing symbols.
After k>=0 maximal pairs let j=h-2k. If j is zero, C0 writes it to one and
exits (j+1,A), with time 2k+1 and exactly one edit, population +1.
Otherwise maximality forces cells (j-1,j) to be `11`; C1 D1 erases j-1 and
exits (j-2,A), with time 2k+2 and exactly one edit, population -1.
Each additional pair contributes two native instructions, no edits, and a
left displacement of two. This proves both branches for every k, not just
for sampled k. Finite support forces termination.

All three macros preserve the full tape, signed head, control, and positive
native clock on their section. They need not be shortest return maps; the
bouncer A1 branch intentionally combines a maximal sequence of A self-loops.
None reduces the retained variable-length word to a fixed collection of counts.

## 3. Ordered-join encoding and exact theorem hypotheses

G0 is the transitive triple with arcs 0->1,1->2,0->2. G1 reverses only 0->2,
so it is C3. Define E(w)=G(w1) join ... join G(wm), with earlier blocks beating
later blocks. For equal-length words u,v,

`d_iso(E(u),E(v)) = Hamming(u,v)`.

The pinned dependency is
`01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md`.
Its metric is the minimum arc disagreement over vertex bijections, **isomorphism
only**, without complement folding. Its hypotheses are equal corresponding
factor orders and an order-join with a uniform forward cross-cut. It requires
neither strong factors nor a uniquely marked/maximal cut. Every paired factor
here has order three, so successive prefix orders also match. G0 and G1 are
nonisomorphic and differ by one arc, hence have distance one; induction gives
the formula, including the empty-word boundary with distance zero.

The dependency proof was inspected: any block-mixing bijection has a source
pair a in the first block, x in the second whose images cross. Swapping those
two images repairs their own cross-pair and never increases the sum of the
other two incident costs, in any of the four source/target block cases.
Thus mixing strictly loses to a block-respecting map; the latter costs the
sum of independent factor distances. This argument explicitly covers merged
transitive zero blocks. It is not an inference from SCC sizes alone.

The separately pinned
`01-canon/theorems/THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md`
correctly marks the limitation: cyclic quotient substitutions do not enjoy
this law universally. No cyclic-substitution generalization is made here.

### Decoder and the zero-block issue

A path can never go backward across an order-join cut. Each one-block is
therefore a C3 strong component; each zero-block contributes three singleton
strong components. The condensation is canonically ordered. Between consecutive
C3 components, and at either end, a run of 3r singletons decodes to r zeros.
This proves injectivity of the word encoding up to isomorphism, including all
zero words and adjacent zeros; it does **not** claim that original zero triples
remain maximal modules. Invalid singleton-run lengths or strong blocks other
than C3 are outside the image and must be rejected.

### Invariants and genuine information loss

Every directed Hamiltonian path visits join blocks in their order. G0 has one
such path and G1 has three. No triangle spans join blocks. Consequently
`H(E(w))=3^(number of ones)` and `c3(E(w))=number of ones` for all finite words.
We count directed vertex sequences, not isomorphism-orbits of paths; H(empty)=1.

For words 01 and 10 the scalars coincide, while distance is two. Bouncer+bell
at head 0 in state A maps 01 to 11 in two instructions and leaves head 0;
it maps 10 to 00 in one instruction and leaves head -1. Population, H, and c3
therefore do not determine the next macro. The symbols' ordered positions matter.

### Literal arc operation, actual clone changes, and clock separation

Flipping bit i reverses exactly arc (3i,3i+2), never a cross-block arc. An
unchanged native write reverses no arc. For 00, the tournament is T6 and its
twin pairs (equal orientations to every third vertex) are
`(0,1),(1,2),(2,3),(3,4),(4,5)`. After writing the first bit to one, the word
is 10 and only `(3,4),(4,5)` remain. In particular the literal cross-boundary
pair (2,3) loses its twin relation: vertex 0 distinguishes them. Erasing that
one restores the relation. These are actual split/merge effects, not a visual
analogy; the first write is realized by A0 on the two-zero window.

For every k>=0, place C on the rightmost symbol of `0(01)^k`. The k pair loops
make no edits; C0 then changes only the initial zero. Thus tournament distance
is one while native time is 2k+1. For k=0 add a zero padding cell to contain
the outgoing head. This all-k proof rules out any universal upper bound on
macro time that depends only on edit distance. It does not rule out bounds
that also depend on window length or other retained data.

A finite trajectory can be represented by E(window word) plus its origin,
head index, control, and program, provided the window includes every read,
write, and head position. A wrapper must detect any attempted exit. This
recovers the configuration **within that window**, not arbitrary remote tape
cells. The bare tournament forgets head/control/origin/program; native dynamics
is consequently a controlled marked local rewrite, not an autonomous graph map.
No fixed finite window can be presumed to contain an arbitrary infinite orbit.
The rewrite is not K_P, and none of the Hamiltonian-word/successive-source
blackness/literal-K-target obligations for the endpoint-four theorem are supplied.

## 4. Exact finite controls and independent receipts

All normal and optimized checks passed; the two independent outputs are byte
identical. The separate producer reruns also match its sealed receipts exactly.

* **156,934 macro cases / 599,010 native instructions:** the exact producer
  universe has 112,640 short cases plus 1,005 long cases = 113,645. The short
  universe is every binary word of length 0..10, every head -1..n, and every
  stated section state. Independent branch counts match the producer receipt.
* Added 25,580 translated signed-head cases (origins -17,23; length 0..7;
  heads two places beyond both ends), retaining remote one-cell sentinels at
  +/-1,000,000; 5,000 deterministic sparse signed tapes; 49 long templates
  including k=1024; 12,658 blank-orbit macro calls; and two position hostiles.
  The additional direct hostiles are outside this macro-case count.
* Each blank orbit was replayed to its first macro endpoint at/after 10,000
  native steps. Exact outcomes match the producer: 10,005 sandwich, 10,001
  bouncer+bell, 10,056 irregular steps. These are not infinite-orbit evidence.
* **21 equal-length word pairs through length two:** all 11,545 vertex
  permutations, including every order-six bijection, give Hamming distance.
  All Hamiltonian paths and all triples were enumerated for the seven words.
* **511 words through length eight:** decoded after both reversed and seeded
  shuffled labels, yielding 1,022 independent decoder checks. Empty words and
  merged zero runs are included. Three nonimage tournaments are rejected.
* **3,586 one-bit changes** in these words alter only the distinguished internal
  arc. **3,096 actual native transitions**, across all four states of all three
  programs, length 1..5 windows, all words and all head cells, match the marked
  arc rule. Of these, 2,352 keep the head inside and 744 explicitly exit.
* Hostiles separately retain signed coordinates, reject forgetting head or
  control, reject replacing guard 000 by 00, and check exact differing outputs
  of 01 and 10. Remote sentinels test preservation of untouched tape cells.

Finite tests are not substituted for the all-parameter proofs in Sections 2–3.

## 5. Source and provenance audit

The sealed producer hashes and both theorem hashes match; see `normal.json`.
`SOURCE_AUDIT.json` records exact source URLs and inspection outcomes.

The [Bouncer page](https://wiki.bbchallenge.org/wiki/Bouncer), whose rendered
permalink is oldid 8027, gives an inductive analysis for a different two-state
machine. The [Bell page](https://wiki.bbchallenge.org/wiki/Bell), oldid 6789,
gives ordinary/cubic examples. Neither is silently substituted for the requested
four-state machine. The primary [Yongyi Chen discussion](https://groups.google.com/g/busy-beaver-discuss/c/lcr_6buFz_8/m/e2BiAob8AgAJ)
contains the dated finite recurrence report and tentative sandwich equivalence.
Those historical computations were not independently rerun and prove no
nonperiodicity assertion here.

The primary [pinned Bouncers.v](https://github.com/meithecatte/busycoq/blob/1c23203761602c199a7ac5ea6f7947fb80851d93/verify/Bouncers.v)
was retrieved through the GitHub connector and inspected as source text only.
Its 20,685 bytes equal the fetched master version. A read-only local snapshot is
`Bouncers.source.txt`, SHA256
`2a2455657fa4eabcbfeff11517307580e5fb656a1f6664baa12203d10dbf0d17`.
Its pinned upstream copyright/permission notice is retained in `Bouncers.LICENSE.txt`.
The source contains positive-repeat semantics, `shiftrule_left_exec`,
`steps_some`, and the theorem `verify_bouncer_correct`. The source theorem
requires successful certificate verification and concludes blank-tape nonhalting.
No certificate for the three audited machines is asserted.

**Adopted wording refinement:** the source's alignment/subsumption lemmas preserve
matching membership through forward implication. Subsumption explicitly states
`t ~ xs -> u ~ xs`; it is not equality of the matching languages. At the parent
agent's explicit request, the sole producer prose amendment replaced “preserve
matching languages” with “preserve membership by one-way inclusion, not language
equality.” Producer code, source manifest, and both computational receipts are
unchanged. `AMENDMENT.json` records old/new report pins; the producer checksum
manifest and this checker's input pin were refreshed, and both audit modes rerun.
No macro, metric, or decoder theorem changes.

Coq was **not built**; imports/axioms and root reach were not audited. Therefore
formal verification remains CITED, not newly VERIFIED. No downloaded software
was executed. Only the authored Python checkers were run.

## 6. Reproduction and remaining boundary

From this audit directory:

```sh
python3 -B independent_check.py > normal.json
python3 -B -O independent_check.py > optimized.json
cmp normal.json optimized.json
python3 -B ../session65-cryptids/check.py > producer-rerun-normal.json
python3 -B -O ../session65-cryptids/check.py > producer-rerun-optimized.json
cmp ../session65-cryptids/normal.json producer-rerun-normal.json
cmp ../session65-cryptids/optimized.json producer-rerun-optimized.json
sha256sum -c SHA256SUMS
```

A genuine further dynamical result would need a closed blank-orbit word grammar,
with exact boundary guards, control/head ports, and native clock. The present
static carrier and complete variable-word macros do not supply such a smaller
closed recurrence. Recent Hydra/A6 work on other, partially defined machines
cannot bypass those obligations or change these three machines' totality.

**Blockers:** none for this audit's stated deliverables. Direct web requests for
some oldid pages and GitHub source returned cache/access errors, but the rendered
wiki pages supplied the same revision IDs and the GitHub source connector
successfully returned the exact pin. No action was denied, no unconfirmed
mutation was retried, and no approval is pending. The unbuilt Coq stack is an
explicit scope limitation, not a passing build claim.
