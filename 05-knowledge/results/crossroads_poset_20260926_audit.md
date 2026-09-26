# Audit: families beyond 27, arithmetic posets, and fixed-integer carriers

**Status: PROVED elementary additions + FINITE-EXACT controls + independent
proof audits. Submitted global poset results remain claims under review.
Collatz remains OPEN.**

Run from the repository root:

    python -X utf8 04-computation/experiments/crossroads_poset_20260926_audit.py
    python -X utf8 04-computation/experiments/crossroads_poset_20260926_audit.py --verify-manifest

The first command replays five lanes normally and with `-O`, and compares
both with retained stdout. The second verifies the content pins in
[the manifest](crossroads_poset_20260926_manifest.json), using LF-normalized
bytes. `--record` intentionally refreshes outputs; `--manifest` intentionally
refreshes pins after replay. Load-bearing predicates use explicit errors,
not assertions removed by optimization. Scripts use the Python standard
library and exact integer/rational comparisons. Decimals are display only.

## Finite universes and independent paths

| Lane | Universe and checks | Positive / hostile controls |
|---|---|---|
| bridge |119 fixed-count word universes through length14;3471 legal swap checks;200 ternary lifts;actual family through6833-bit input;exact counts through X=2^320;all odd sources<=100000;closure covers through2^256 |independent minimal-element recursion vs binary words;direct integer trajectories vs affine formula;27 exception;27/31 height selection |
| integer |12800 address checks;8174 capped-cylinder checks;exact weighted mass through depth9;finite tail atH64;all starts2..1000000 for first descent with10000-step failure cap;hover/drop m6..96 |negative fixed-point shadow;sparse computable infinite address;positive3n-1 and5n+1 cycles;Z2=Z3;27 plateau;actual distinct large hover/drop sources |
| barrier |5230 naturally labelled posets on2..6 elements;252265 extensions;4376 displacement-triple checks;1089374 variance triangles;91196 ideal cuts;96426 antichain windows;289 rigid triple gluings;1230 rational case-D triples;63000 Fibonacci boundary controls |direct extension enumeration vs matching formulas;fixed central pair vs moving boundary maximum;normalized moments vs maximum |
| width |5231 naturally labelled posets on1..6 elements;252266 extensions;601302 ordered XYZ triples;76193 pair-window bounds;96427 each ideal and antichain checks;8190 parity words through12;151 integer-law checks |full extension law vs induced-subset restriction;four-order negative covariance vs selected positive covariance;nonconvex chamber support;exact many-period TV |
| Moran integration |50 exact finite stopping cases, including17 empty-hit events;2000 greedy threshold prefixes |undivided stopping identity vs undefined empty-event conditional;finite hitting vs unattained supremum3;prefix scale vs integer-height scale |

The poset census is every transitive relation compatible with a fixed
natural topological ordering, not every labelled poset. Every isomorphism
type at those sizes has such a labelling. The six-event height obstruction
is first only in the explicitly probed universe, not under every conceivable
encoding. The first-descent record list is finite-exact through one million;
no asymptotic density or recurrence for those records is inferred.

## Proof audits and scope

Root derived the family and two-chain bridge. Automata independently derived
the bridge and audited the written family, carry and closure statements.
Geometry and flow independently checked family counting, injection, phase
recursion, closure and the O(J) cover. The exception r=0,n=27 is explicit.
An initial overbroad generic-tail stopping-time reading was narrowed: a=1
can already start at1, and generic r=0 exceptions depend on v3(a+1).

Geometry audited the fixed-integer address, cap, endpoint refutation and
record repair; automata independently accepted those and the atomic/Abel
reduction. Root re-derived the hover/carry and suffix bounds and checked
the finite record interpretation. Geometry and automata reciprocally audited
the elementary Fibonacci boundary and height-conditioned law additions.
The wording about order statistics was repaired: the underlying uniform
draws are independent, the order statistics themselves are dependent.

Promoted results are [THM-4502](../../01-canon/theorems/THM-4502-collatz-41-tail-growth-family.md)
and [THM-4503](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md).
They have elementary in-repo proofs, no dependency on accepting the new
submitted papers, and no formal Lean claim. The other elementary statements
and exact reductions are retained in the lane notes with complete proofs.

## Source audit

- KSBFT_v7.pdf:31 pages, SHA256
  `48bfecea2d4b47698861bd3eeb544f350617e5828a018e1511c67efdf1ddb759`.
- Kahn_Saks_Proof.pdf:25 pages, SHA256
  `83bc04b19ca1790913f16183322764a5537ea1aa1101d1e4b8a0cbbd43681897`.

Both local PDFs were text-extracted and material pages visually checked;
the lane notes enumerate these pages. The full internal proof structures
were read. No fatal internal error was found, but cited external inputs
were not all independently re-proved. The KSBFT global epsilon arithmetic
and Kahn–Saks quantitative rate are not independently certified. The
[author announcement](https://igorpak.wordpress.com/2026/09/25/how-to-make-small-improvements-and-break-barriers/)
corroborates source provenance; it does not certify truth or byte identity
between a public version and the local attachments.

## Corrections and residual obligations

- HYP-9161's literal arbitrary-segment average is REFUTED by a parameterized
  actual positive distinct-orbit construction. The O(k) endpoint repair and
  an entire-injective-orbit version remain OPEN. THM-4476/4499 are unchanged.
- The record-peak proof is repaired by suffix carry at the large endpoint.
  The original denominator lacked the necessary lower bound. Final/initial
  missing windows contribute O(k); the asymptotic result survives.
- XYZ covariance reversal refutes only the height-law transfer. The
  submitted paper's hypotheses exclude that selected nonconvex law.
- Atomic-mass decay is equivalent to Collatz first descent; no decay bound
  is supplied. Computability of each rational Z_k is not convergence of Z_k.
- Recursive convergent examples, finite census records, entropy dimensions,
  and the Fibonacci/Zeckendorf233 encoding are distinct statements. No
  arithmetic map to the earlier tournament233 is asserted.
- The incoming THM-4504 integration retains its original large census as
  attributed finite evidence; it was not re-run here. New exact controls
  repair empty-hit conditioning, the open-set definition, and height-scale
  exponent conversion, and restore model/empirical labels in its summary.
  The common ceiling-critical word is excluded even from rational starts:
  scaling an odd denominator gives a fixed-odd-b integer orbit of linear
  growth, contradicting THM-4476. Root and flow audited this corollary,
  detailed in the board; it is not a claim about arbitrary infinite words.

Current synthesis and next targets: [board](crossroads_poset_20260926_board.md).
