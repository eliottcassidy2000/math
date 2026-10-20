# Scope and validation receipt

**PROVED:** the elementary all-order arguments in the four notes.
**FINITE-EXACT:** only the universes listed below. The computations are
controls and literal-witness audits; they are not the basis for extrapolated
all-order claims. No Lean formalization is claimed.

From the repository root, run:

    python3 05-knowledge/results/trijection_reset_20261020/validate.py
    python3 agents/check_docs.py
    git diff --check

The first command runs every checker normally and with optimization, fails
on an error or a receipt mismatch, and saves canonical JSON receipts and
[a source/receipt hash manifest](validation_manifest.json). Every mathematical
guard uses an explicit exception, so optimization does not erase it.

| Checker | Exact universe / filters | Independent path or hostile |
|---|---|---|
| `independent_check.py` | All path-normalized tournaments n=1..4; all20 labelled three-edge graphs on4 vertices; no extra filters | Arc permutations for H and Aut, direct reachability versus interval coverage, direct cycles versus score formula; same unlabelled one-edge graph lifts to three types |
| `arithmetic_check.py` | k=2..40 symbolic exponent boxes; three marked prime assignments at k=2..8 (21 numeric cases); all10 tetrahedral corner involutions | Trial-division divisor sets and numeric divisibility versus coordinates; explicit inverse and removed-endpoint collision; p², p^q*r and p^4qr fail balance |
| `square_counts.py` | Ordinal powers Q^r, r=0..4; every path and every cell of the25-by25 square; no extra filters | Direct full-tournament path DFS versus block products; actual anti-automorphism; same-endpoint blue/black words; diagonal deletion versus pointed deletion |
| `reset_composition_check.py` | All1365 words of length0..5 in valuations1..4, all7737 cuts, first three positive cylinder sources per word (4095 replays); odd sources3..1001 with first-descent endpoint>1 (492 pairs) | Affine recurrence versus suffix sums and native odd-map replay;1113 expanding hostile,13 mixed-cone hostile; descent search cap10000 is an explicit error, never a negative result |

Proof review is separate from the tests. The graph referee independently
checked the proposed coordinate map and supplied the interval-cover proof.
The arithmetic agent derived the divisor chart from recovered sources;
the parent checked the kernel/inverse, action, and parity boundary. The
arithmetic agent separately audited the parent's square construction from
a hand-listed arc set and all24 permutations. The graph referee separately
audited the Collatz agent's gauge laws, exact-cylinder proof, source pullback,
boundary parity and source-payment inequalities. These are independent-agent
reviews within this session, not external refereeing.

One draft error was corrected before publication: path words from an
isomorphic four-tournament were initially used with the wrong arc labelling.
The direct legal-arc checker exposed the impossible30 step. The actual paths
and reflection-fixed words are now listed in `square_counts.md`; the
mechanism is recorded in `01-canon/MISTAKES.md`. The count5 and all ordinal
product statements survived, and the independent review used the repaired
literal representative.

Existing direct predecessors were checked rather than rediscovered:
LEM-003 and THM-4729's elementary tiling layer; the defect-boundary
tetrahedral graph calculation; the divisor-layer10-point and complement
obstruction; the balanced-family square and actual4,2,1 deletion map;
the positive-charge first-exit and exact-cylinder notes. Provisional
THM-4730 and THM-4729 L4 are not dependencies. No literature priority
claim is made. The parent chat's closed-blue endpoint-three proof was not
duplicated.

Before publication, incoming commits70cb917742 through43e4f5c6a9 were
read as mathematical input. Their forward debt cone overlaps R1 and adds
terminal payment and unbounded-reset examples; these are attributed and
integrated as R3. The zero-slack-coefficient c=0 proof boundary was checked
separately. Incoming endpoint-two work remains separate from the parent
chat's endpoint-three target. The new palette note's576-element extension
product is compared with, not identified with, our Hamiltonian-path product.

The final clean rebase also incorporated8dcc679ae3 and1e07ac8243. These
independently repaired the same c=0 strictness issue and sharpened the
separate unit-gap mass layer. The current charge source already contains
the dyadic-denominator proof, so no duplicate source repair was needed.
The parent checked R3's terminal payment proof against that incoming
version. The only merge conflict was reconciled in `MISTAKES.md`, retaining
all incoming corrections and our own labelled-path correction below the
maintained opening text. Normal/optimized receipts remain unchanged.
