# Scope and reproducibility

Same-agent proof review with independent computational paths. The complete
positive-odd spectrum and a closed physical ternary source family are OPEN.

## Explicit universes

* `integer_census.py`: direct odd partitions at every order0..30;
  integral and rational projective algorithms at60,100,150,200; original
  implementation at60,100. Known float-rounding failures at14/15 are exact
  controls. Every division checks its remainder; Python -O removes no checks.
* `bucket_census.py`: both integral gauges through30 and four exact
  engines at100,150,200. Same projective state counts; all power-of-two
  divisions checked. The older full-prefix baseline remains faster at200
  in the first recorded run; no overall speed record is claimed.
* `exterior_search.py`: every196417 rooted unique-start parent class at
  orders3..15, with literal DP on the two4277 positives and all eligible
  matching-one exterior clones. The all-order consequence additionally
  uses the proved parent reduction and cited Busch bound.
* `offset_fibres.py`: the same complete grammar, five declared offsets,
  seven positive bases and34 clone queries. The proof separately handles
  clone marks and the excluded unique in-port, including moving marks.
* `exterior_two.py`:1387 saved two-entry core presentations through9 and
  all legal tails through parent15. This is a bounded sublanguage; no
  complete two-ready exclusion is claimed.
* `audit_feedback.py`: all2730 four-letter words of lengths0..5 and both
  feedback-edge choices, compared with literal Hamilton subset DP.
* `feedback_search.py`: all prefix-clone/two-feedback tails through13,
  using monotone insertion pruning. Zero hits, so the100-hit cap never
  truncated the search. The cited strong-order bound justifies the
  all-size conclusion only for the architecture in FG1.
* `next_fibres.py`: complete unique-start parents through17,1346268
  classes, for12833/12835. Uses the incoming graded C++ compiler with only
  its target list changed, recorded in `next_bank.generated.cpp`.
* `family_certificate.py`: five positive laws at m=1,2,3, all-size extension
  by SF2, a full seed eighteen-polynomial observer, literal two-birth12833,
  and five Q4/charge phases at root and growth port.
* `source_pins.py`: the actual18979 source, both incoming native traces,
  retained/deleted old names, root frames, next counts, fresh phases,
  actual full automorphism orders81/27 and eight named-cut responses.
  The full source observer is inherited unchanged and hash-pinned; its
  root projections are independently checked. Full source rows are kept.
* `splice_mobius.py`: all endpoint guards for the two fixed-anchor square
  catalogues, each checked by independent subset matching; all proper
  intersections in parity controls r=2..5. SP2 proves the general-r claim.
* `first_ready.py` and `audit_first_ready.py`: complete parent banks
  through17 for11735/11737, seven literal positives, the actual repaired
  m=2 source, a22-native return with minimum old-name loss2, and next35207.
  The independent audit uses separate matrix and Hamilton routines;
  all completed repairs are strong, with no native-cost minimum claimed.

## Quantifier and type boundaries

One-ready closure is narrower than two-ready closure and much narrower
than arbitrary realization of a value. A fixed uniform clone module is
narrower than a split module or a growing exterior. The new SF3 result
allows arbitrarily complicated fixed exteriors and moving ready marks,
but does not drop these module hypotheses.

The actual18979 receipts are incoming work, not newly discovered returns.
Old deletions are explicit. Native counts do not assert intermediate
helper-expanded strongness. Golden/field readouts are not fabricated
interval-growth ancestry; the saved original histories remain separate.
The failed55->111 one-birth certificate is inherited and source-pinned,
not enlarged into an unrestricted impossibility.

A000568 counts isomorphism classes. It does not recover their path-count
fibres, full automorphism actions or metagraph edges. Exact counts at200
are an implementation extension of the comparison suite, not a literature
record. Runtime ratios are observed once on this machine and stored
separately from deterministic receipts.
The final Session105 recovery shows that both integer gauges and their
packed bucket realization are inherited too. IC1/IC2 are explicitly
credited reconstructions; the numerical controls remain valid. Final
integration changed reference-document pins only, not computational inputs.

SP1 counts paths for a fixed endpoint skeleton, and its full square
application is complete only for the specified anchors and minimum
connector count. Pairwise intersection data fail in the explicit parity
hostile; no implication from that failure to an arbitrary mathematical
conjecture is asserted.

## Reproduction

From the repository root:

    python -X utf8 -B 05-knowledge/results/second_feedback_18979_20261103/reproduce.py
    python -X utf8 -B -O 05-knowledge/results/second_feedback_18979_20261103/reproduce.py

Both modes compile their own C++ tools, regenerate all thirteen deterministic
receipts, and compare LF-normalized SHA256 pins. The dependency set
includes the inherited failed55 source certificate and previous81/27
construction. Timings regenerate but are intentionally not hash-checked.
The original measurement quoted in the summary is retained in
`benchmark_observed.json`.
`--pin` is a maintainer action after review, not part of verification.

The stored normal and optimized verification files cover the first eleven
programs. After adding the constructive11735 branch, only its two new
programs needed rerunning, in both modes; their separate receipts are
`verification.extra.normal.json` and `verification.extra.optimized.json`.
The default command now replays all thirteen. For that last targeted batch:

    python -X utf8 -B 05-knowledge/results/second_feedback_18979_20261103/reproduce.py --only=first_ready,audit_first_ready
