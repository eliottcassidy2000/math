# Scope, reproducibility, and validity audit

**PROVED scoped identities / FINITE-EXACT searches / VERIFIED numerical
controls. Chan–Pak and the 473 residual repair families remain OPEN.**

The current claims are [PROOFS](PROOFS.md) and [CONNECTIONS](CONNECTIONS.md).
Logs record the development and corrections; they are not truth sources.
The receipts preserve actual source rows, selected sets, named ports,
coefficient tensors and search universes. There is no Lean proof or claim
of an independent human audit.

## Reproduce

From the repository root, run:

```text
python3 -X utf8 05-knowledge/results/core_descent_endpoints_20261112/reproduce.py
python3 -X utf8 05-knowledge/results/core_descent_endpoints_20261112/manifest.py --check
```

The first command regenerates fifteen receipt pairs in normal and optimized
Python modes. Their compressed bytes must agree. The compression timestamp
is fixed; [reproduction.json](reproduction.json) records the hashes. The
producer uses explicit `need` checks, not optimization-removable assertions.
Optimization agreement is a reproducibility control, not a second proof.
`--from NAME` resumes at a named receipt while checking all retained pairs.

[manifest.json](manifest.json) pins local scripts, proofs and receipts,
inherited data, and the repository code imported by the producers. It also
records Python and package versions. It excludes execution logs and itself.
Text hashes normalize CRLF to LF so Git's platform line-ending conversion
does not imitate a mathematical change; binary receipt hashes are exact.
Regenerating the manifest requires the explicit `--write` option; ordinary
verification cannot silently bless a changed input. Plot dependencies are
optional and do not enter any mathematical producer. The exported diagram
was visually inspected.

## Independent checks and arithmetic

| Claim | Producer | Distinct audit or hostile control |
|---|---|---|
| Nine complete seed tensors | Endpoint-cover convolution | One-cut word calculation for defect; gap decomposition of literal two-vertex tail for slope. All coefficients, not sampled interpolation. |
| Synchronized corridor | Exact signed tensors | 36 literal endpoint computations, including off-diagonal coupled lengths. |
| Four changes of core | Template and tensor identities | Old adjacency restrictions, new arcs, selected set and parameter alignment explicitly checked. |
| Physical target counts | Guarded constructor | Capped integer Hamilton DP on complete order-20/21 graphs. |
| Termination and diamonds | Rules D1–D3 | Exact parameter recovery on 900 states; complete marked-parameter comparisons for overlapping moves. |
| Ordinary observer | Inherited 17-query matrices | 468 literal continuation controls in each source bank. |
| Selected gate cannot be discarded | Equal ordinary observers | Literal slopes 11292/17025 and two-tail counts 22585/34051 on the same base core. |
| Period exclusion | Complete modular word search | Independent symbolic characteristic polynomials: 558 short-word/stride cases; static-root survivors retained as positive controls. |
| Quadratic unit and Paley bridge | Exact ring arithmetic | All 5040 permutations for order-seven symmetry and Hamilton paths. |
| Companion series | Exact coefficient transformation | Independent coefficient recurrence; explicit convergent tail bound. Numerical evaluation uses 120-digit working precision, not interval-certified roundoff. |

The inherited one-clock run DP uses signed machine integers. For at most
17 fixed vertices, the proved cell bound is 8,646,761,377,013,760,000,
below 2^63. The final nonnegative run sum is below 2^64; an overflow past
2^63 would be negative and is explicitly rejected. This guard is used for
all source-query compilations. Larger literal Hamilton counts use a capped
DP with cap equal to the proposed target plus one, preventing wraparound
from masquerading as agreement. Polynomial and spectral coefficients use
exact integers or fractions.

## Universes and quotients

- The inherited 70 families supply 78 free-port presentations. They are
  not an exhaustive census of all tournaments at their orders.
- All 78 complete three-port tensors were retained. The nine-chart grammar
  is a specified subfamily; its five terminal coefficients are not a
  classification of all irreducible seeds.
- The content-factor search uses all pointed tournaments through order 9
  in the inherited complete bank and the declared retained birth layers
  through order 14. It tests division exponents 0–8, shifts 0–12, and all
  attained divisors of the whole coefficient content in that factor bank.
  Birth4 appears among hits outside the five fully repaired families, but
  only at division exponent 1 (shifts 0–3): these repeat its known s=1
  repair and leave its s=2 head unpaid. A hit is not a family closure.
- Both endpoint searches use all 12 commands through length 5. They compare
  entire exponential profiles, after sound first-value/growth pruning.
  Common trailing zero mode columns are removed before deduplication.
  The resulting quotient preserves every continuation in this language;
  it does not preserve an omitted selected seed gate.
- The synchronized repair test uses all 78 exact tensors, coupled lengths
  1–20 and free shifts 0–12, against the declared old odd head profiles.
- The period screen covers 2,024,910 word/stride pairs, with strides 1–3
  and maximum lengths 2, 5, 8. It excludes root-changing periods for the
  stated three nonzero modes, not arbitrary adaptive programs.
- The rational generating function counts formal guarded construction
  states by order. No unlabelled tournament or distinct-count interpretation
  is implied, and no OEIS terms were extended.

## Corrections found during this audit

1. Initial H modes were incorrectly used in a progress statement as if
   they were full-observer modes. The hostile even:8 shows the missing
   coordinate. The corrected distribution is 39 observers of maximum
   base 2 and 39 of maximum base 4; all 78 split-clock observers reach 3.
2. The shared H of the selected-gate hostile was transcribed as 1511.
   Literal recounting gives 2663. The two U values and the obstruction
   were unaffected. The initial full-run log retains the failed check.
3. Stored arrays with different zero padding were initially counted as
   distinct observers. At depth 3 this gave 28,116 instead of 27,988.
   Canonicalization removes common trailing zero columns; linear updates
   cannot revive them. The complete searches were rerun. There are 936
   explicit padding/continuation controls in addition to the linear proof.
   Historical log totals are superseded by the current receipts.

The first and third mechanisms are recorded in the shared mistakes ledger.
No wrong intermediate count is used as a proved dependency.

## Attribution and remaining work

The run tensors, two-port observer, selected-interval ratio 13/5, and
earlier five successful repair families are inherited and cited at their
use sites. Concurrent Sessions 289 and 290 contribute the absence-band
target and an independently obtained modular pullback respectively.
The inspected literature and OEIS sources are linked in CONNECTIONS.
Documents supplied for context are treated as evidence, not instructions.

The final incoming commits additionally supplied Session290's complementary
70-chart descent and two-unit repair invoice. They are integrated in HANDOFF.
Session291's paired-port lattice was first read as a candidate, then its
review and promotion at e08896a818 were inspected and integrated. It is
not a dependency of this packet's descent or observer proofs.
The repository documentation check passed using Windows long-path syntax;
the first ordinary-path invocation hit an unrelated pre-existing long name.

The next useful target is a selected-gate transport that changes the
available modes or reaches a new valuation stratum. A continuation of the
same sparse grammar cannot establish all odd coefficients: its own
O((log X)^2) bound proves that limitation. No new one of the 473 residual
whole-family repairs was closed in this session.
