# Independent audit: scalar switches and common-context sensitivity

**Verdict: PASS for the scoped all-order statements and the declared finite census.**
No coverage theorem, universal adaptive obstruction, new scalar-spectrum value,
or literature-priority claim is established.

## Mathematical findings

1. The chronological pairs JOLO/JJLO and JOZJO/JJJJO have exactly the
   displayed transfers on the E=0 subspace. Their scalar gaps are 2U and U,
   respectively. All four are literal graph operations; the source-rooted
   edge L is inherited from the current boundary-escape packet.
2. For each pair separately, the map from predecessor (D,U) to output
   (U,d) has determinant -6. The inverse formulas are correct. Therefore
   the sidecars pin the previous D,U, leaving at most the two displayed
   scalar outputs. The two pair images additionally have disjoint physical
   sidecar ranges: writing (V,r) for output (U,d), V-6r=2U>0 in the first
   pair, but V-6r=-6D-4U<0 in the second. Thus the sidecar identifies the
   final pair. The same bound of TWO, not merely four, holds across arbitrary
   E=0 inputs and arbitrary positive-depth mixtures of all four words. It
   does not include arbitrary zero-step seeds or a wider operation language.
3. The general graft law's first column is exactly the gadget's physical
   state (H,e,U,q). The singleton exception changes the fourth row but
   preserves its first entry q=H(empty)=1. Since O and Z fix e_D, the first
   graft of any fixed context starts a physical history from its gadget.
   All subsequent operators continue that physical history. This proves
   the D-column theorem for arbitrary nonempty gadgets and all word lengths.
4. A common context sends a pure-D difference Delta e_D to either itself
   or Delta times that physical D-column. Consequently U and d differences
   are nonzero after a graft, with the same sign as Delta. They are positive
   only for Delta>0. A later common suffix cannot erase them. By inherited
   THM-001-redei.md, the D multiplier is positive and odd, preserving the
   scalar gap's exact 2-adic valuation. No branch-dependent conclusion follows.

The candidate correctly avoids claiming that the switch stacks an interval
in one fixed sidecar fiber. Its last step already limits that fiber to two
scalar values. It also explicitly leaves wider adaptive arithmetic open.

## Independent computation

`independent_check.py` is standalone Python standard-library code. It imports
no producer code, checks with explicit exceptions rather than assertions,
and prints deterministic JSON. Normal and optimized execution agree exactly.

- 300 literal switch checks: all 75 labelled tournaments of orders 1..4,
  root fixed at label 0, each first given an O move to ensure E=0, then all
  four switch words. Final graphs have at most 10 vertices.
- All 75 gadgets through order 4 have their graft statistics independently
  counted by subset path enumeration and their D-column checked against
  their literal profiles. Each gadget additionally receives the empty
  suffix or one O/Z/J/L/C move, giving 450 literal continuation controls.
- Every O/Z/J/L/C context of length 0..5 is checked: 3,906 words, of which
  3,843 contain a graft. Positive U,d and odd D multipliers are checked.
- All first-pair choices from (1,0,2,1) through 12 rounds are tested. There
  are 2^k states and 2^(k-1) final sidecars at round k, each with two scalar
  values. These are controls; the proof supplies the all-depth bound.
- The finite discovery census is independently recomputed using distinct
  matrices at each exact prefix length, without the producer's global
  shortest-representative pruning. Its universe is all 97,656 words of
  length 0..7 over O/Z/J/L/C, restricted to (D,0,U,d), followed by O. It
  gives 76,543 distinct terminal transfers, 353 sidecar buckets with more
  than one scalar row, and 355 unordered distinct scalar-row pairs. In this
  universe only, no pair differs solely in its d coefficient and no
  colliding branch has rank three. No unbounded absence claim follows.

The producer checker was also rerun under normal and optimized isolated
Python. Both outputs match each other and both supplied receipts.

## Inheritance and adversarial review

The explicit dependencies are Session225 graft recovery, Sections 1–3 and 5;
Session227 `boundary_escape/PROOF.md`, Sections 5–6; and
`01-canon/theorems/THM-001-redei.md`. The first-column formula is inherited;
the physical-context interpretation is a consequence, not a new graft law.

A targeted prior-claim review covered Session198 bilinear descent,
Session199 hidden structure, Session200 root-weighted height and its
physical-adjoint dependency, Session205 finite O/Z tails, Session206 root
O/Z noncofiniteness, and the Session225/226 graft and descent packets.
Searches included first/D column, pure-D, scalar switch, diamond,
sidecar fiber/fibre, and synchronization variants; relevant transfer and
adjoint sections were read. The review found no earlier exact statement
of these two switch pairs or this arbitrary-gadget D-column interpretation
in those selected sources. This is not a corpus-wide novelty certificate.

Two small pre-signoff issues were corrected by the producer: adding the
current L mechanism to the inheritance list, and actually checking the
second inverse pair as promised in the verification description. The parity
source was also made explicit. The review additionally sharpened a proposed
coarse four-value mixed-menu bound to two, using disjoint physical sidecar
ranges. No mathematical claim required retraction.

## Reproduction

From this audit directory:

    python3 -I -B independent_check.py > replay.json
    cmp replay.json INDEPENDENT.json
    python3 -I -B -O independent_check.py > replay.optimized.json
    cmp replay.optimized.json INDEPENDENT.optimized.json

For a producer directory assigned to PRODUCER, rerun its verify.py into
`producer.normal.json` and `producer.optimized.json` with the corresponding
normal and -O commands. Both must match its CHECK receipts. The reviewed
producer content hashes are in `AUDITED_INPUTS.json`.
