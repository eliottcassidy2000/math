# Independent audit: the remaining order-ten endpoint-three case

2026-10-21. **PASS / FINITE-EXACT.** The frozen producer packet proves the
stated unfiltered order-ten theorem. No mathematical repair is required.
This audit is independent of the producer's two Python implementations and
checks the literal examples printed in the README, not just their JSON copies.

## Approved mathematical scope

For every tournament with a specified directed Hamiltonian word P on ten
vertices, if P's actual positional reversal is an anti-automorphism,
the word closes from its last vertex to its first, and the endpoint
outdegrees are (6,3), there exist three successive directed black K moves
whose final tournament is the literal labelled K_P(T). K retains the word's
path arcs and reverses every other arc. Black means actual failure of that
word's positional anti-reversal in its actual source, not failure of a
sufficient score test.

Combining this finite theorem with the separately proved and pinned
endpoint-three structural theorem for n>=11 gives the endpoint-score
(n-4,3) conclusion for n>=10. This audit does not re-prove that entire
all-order predecessor. It establishes the missing n=10 boundary, without
source/target mixedness, strong-partner, or isomorphism filters. Neither global
black-core connectivity nor a classification of other endpoint layers follows.

## 1. Complete parameter universe

Normalizing the marked word to (0,...,9) is a bijective change of names for
this fixed marked presentation. Since 9->0 and 0 has three inneighbors,
N^-(0)={9} union I with |I|=2. Vertex 1 cannot be in I because 0->1 is a
path arc. Thus I ranges over all 21 pairs in {2,...,8}. Anti-reversal then
forces N^+(9)={0} union {9-i:i in I}; overlaps of the two port sets impose
no new restriction and are not removed.

On {1,...,8}, the seven adjacent pairs are fixed by directedness. The 21
nonadjacent pairs have three reflection-fixed members and nine two-element
orbits, hence twelve independent bits. The printed orbit representatives
are exactly the canonical representatives, in the same order as the
certificate variables. The independent checker starts from 45 unconstrained
upper-triangle arc bits, imposes path, reflection, closing-arc and port
constraints, and obtains consistent rank 33 for each I. It then substitutes
the proposed twelve-bit parametrization in every raw constraint and checks
that the twelve representative coordinates remain independently free.
This proves surjectivity and injectivity of the stated parametrization,
not just correctness of some constructed examples.

The three cell coordinates are T(1,8), T(2,8), T(2,7), corresponding to
three different middle-orbit bits. Their guard equations have rank three.
Each of the 21 port choices therefore has eight nonempty, disjoint cells
of dimension nine. The checked universe contains exactly 168 cells and
21*2^12=86,016 normalized marked presentations. No complete assignment
census was run, and this is not a claim to count unlabelled tournaments or
all Hamiltonian markings on a fixed external label set.

The ordered-graph compression also has the correct orientation: a backward
nonadjacent tournament pair j->i becomes {i,j-1}. Thus the graph bits are
complements of the forward source coordinates, K_P complements the graph,
and the specific order reflection acts by s->8-s. No arbitrary graph
relabeling was substituted for changing a marked tournament path.

## 2. Certificate soundness and completeness

The new checker represents affine forms by sets of variable names, using
symmetric difference for GF(2) addition. It imports neither producer checker
nor inherited research kernels. It verifies every supplied contradiction and
unit certificate directly and also runs its own elimination, with a different
pivot convention, to establish each system's rank/consistency independently.

For every certificate cell, the checker verifies:

- A unique valid port pair and cell, with no omitted or duplicate cell.
- Exactly three permutation words, retaining every literal vertex label.
- All 27 joins identically directed after the cell guards are substituted.
- Each of the 45 anti-reversal equations at each word's actual intermediate
  source, rather than at the initial source, final target, or a score proxy.
- Every listed XOR proof vector against those exact 48 equations (three
  guards plus 45 anti equations), including index range and no duplicate indices.
- The labelled undirected path-edge identity
  E(Q1) XOR E(Q2) XOR E(Q3)=E(P), and separately the final symbolic arc identity.

Exactly 501 of the 504 systems are inconsistent. A listed XOR yielding
0=1 is a direct certificate that the corresponding word is black for every
source point in that cell. The other three systems have rank twelve; their
36 unit equations pin all twelve coordinates, and the claimed assignment
satisfies all original equations. Thus these are exactly three singletons,
not merely three examples among potentially larger failure spaces.

The printed residual table was parsed and compared with these exact outcomes:

- I={3,6}, cell (0,1,0), closed_epsilon_one, step 3, mask 2202.
- I={5,6}, cell (0,0,0), closed_epsilon_zero, step 3, mask 2446.
- I={7,8}, cell (0,0,0), closed_endpair_yx, step 2, mask 3525.

No other source point can have a blue step in its assigned template.
The counts are 4,536 symbolic joins, 168 literal XOR targets, 501
contradiction certificates, and 36 unit certificates. The certificate's
branch labels are metadata for the proof; correctness rests on the literal
words. As an additional provenance check, the producer's ledger function
has an identical Python AST to the pinned predecessor's ledger function.

## 3. Literal repairs and actual-source witnesses

All three README headings, incoming/outgoing port sets, masks, ten adjacency
rows, and nine displayed words were parsed from the prose. They agree with
the JSON data, the parameter point reconstructed from each mask, and the
actual reflected ports. Every row list is a loop-free complete tournament.

All nine displayed words are directed at the source reached by the preceding
moves. Each has a verified anti-reversal failure at that same source. All
three repaired routes end at the literal stored target rows and at the
independently computed K_P(T), without isomorphism testing. The receipt
records each witness's two arc values explicitly. A valid failure may be
0 versus 1 or 1 versus 0; the producer's first receipt specifically supplies
an existing arc whose anti-image is absent, while its second verifier can
supply either orientation. Both conventions are valid and all 21 stored
witness tuples (9+9+3, including the hostile) were rechecked at their actual
intermediate sources.

The first repair really has endpoint score sums (9,9,9), yet all three
moves are black. This confirms the printed warning that equal endpoint
scores do not establish blue colour.

The separate printed hostile rows
(190,44,760,176,994,704,907,770,559,11) form a closed blue source with
endpoint degrees (6,3). Its one-step cyclic cut is directed and genuinely
blue. The third displayed replacement triple nevertheless repairs it with
three black moves and the exact labelled target. The method's sharp cyclic
cut boundary is therefore not mistaken for a metagraph obstruction.

## 4. Reproduction and hostile corruption controls

Both independent audit runs, normal and python -O, pass and produce
byte-identical receipts. Every check uses an explicit exception, not a
Python assert. Five deliberate corruptions are rejected: omitted cell,
duplicate cell, empty contradiction vector, changed singleton assignment,
and an invalid word.

The producer's check.py was also run in both modes with output paths outside
the frozen source directory. Both regenerated certificates and both receipts
are byte-identical to the frozen originals. The producer's second verifier
was separately rerun in both modes; both receipts match their frozen copies.
The audit never overwrote the producer packet.

From the published packet directory, the portable audit commands are:

    python3 audit/independent_audit.py --packet . --output /tmp/order10.audit.normal.json
    python3 -O audit/independent_audit.py --packet . --output /tmp/order10.audit.optimized.json
    cmp /tmp/order10.audit.normal.json /tmp/order10.audit.optimized.json

To reproduce the stored dependency-verifying audit receipts exactly, add
--repo /path/to/the/exact/inherited/repository/snapshot to both commands.
That snapshot is commit c4a7875652a5ea516e93a31c3c6410d2204773b2. A later
checkout may legitimately fail the optional source-pin check after shared
navigation or correction files change; this does not change the frozen
order-ten certificate. Without --repo, the mathematical checks still run,
but the receipt honestly records that inherited source hashes were not
rechecked in that invocation.

## 5. Dependency and publication boundary

All ten files listed in the producer MANIFEST.sha256 match their hashes.
All eleven producer source pins match the inspected inherited snapshot.
The supplemental dependency manifest also pins AGENTS.md, the research
protocol, and the predecessor's open-end-pair proof probe_new_ledgers.md.
The last file is linked from the pinned predecessor README but not listed
in the producer's source-pins file; adding its hash here improves inheritance
traceability without changing any producer byte or creating a finite-proof
dependency on the all-order argument.

Relevant corrected near misses were checked in the correction ledger:
reflection-fixed pairs can be independent; empty-block XOR can retain
extra labelled edges; singleton reversals cannot be exported to longer
blocks; and K-complement is not the full tournament converse. The present
certificate separately checks precisely those lost coordinates.

The publication manifest is add-only and contains the ten frozen source
files, their original manifest, and this audit's checker, dependency pins,
normal/optimized receipts, report and reproduction summary. Exploratory
search scripts and scratch JSON outside the frozen source manifest are not
needed by the proof and are not proposed for publication. No repository or
remote write was made by the auditor. Publication remains reserved to the
sole publisher following exact root approval.
