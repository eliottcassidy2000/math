# Handoff and validity readback

2026-10-26. This is the producer's readback, not an independent referee report.

## Deliverable

`REPORT.md` gives the typed atlas, seven elementary propositions, the
principal open-problem quantifier ledger, the native tournament inheritance,
and the three-carrier distinction requested by the parent session.

## Exact replay

```
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
```

Both modes passed, with byte-identical JSON. The corpus contains:

* 4,181 two-symbol partial deterministic target tables on 1–3 states;
* 49,758 total encoded vertices, with maximum 15 in one code;
* 45/45 rejected one-arc mutations of the designated canonical codec;
* 1,099 labelled tournaments through order five, automorphism sizes in
  `{1,3,5}`, and orbit counts `1,1,2,4,12`;
* a same-source/read/target machine skeleton with a two-step configuration
  repeat versus a three-instruction halt after changing a write symbol;
* exact rational constant-valuation domain hostiles and finite LRC interval
  controls, including a genuinely empty intersection at threshold 1/2 for
  speeds (1,2).

All tests use runtime exceptions rather than Python assertions. They validate
only these finite constructions and examples, not any famous conjecture.

## Mathematical readback boundaries

* LRC finite decisibility is for **through r** runner cases, with Theorem A's
  lower-case hypothesis discharged sequentially. It does not claim a feasible
  full enumeration or prove LRC(14).
* Decidable does not mean non-NP-hard. NP-hardness requires its own reductions;
  no decidable target can be HALT-hard under computable many-one reductions.
* JC/DC inverse-degree bounds concern actual automorphisms and permit bounded
  algebraic decision; they do not assume the conjectures being tested.
* The stable JC/DC arrows are not fixed-rank iff, polynomial reductions, or
  dynamical conjugacies. The higher-dimensional 2026 record is not re-audited.
* FINBEEP and HALT differ. A natural Pi2 presentation for fixed Collatz does
  not prove that sentence has no hidden Pi1 reformulation.
* The tournament automorphism obstruction concerns faithful symmetry
  transport, not arbitrary graph-isomorphism reductions; colored tournament
  automorphisms still have odd order.
* The incidence tournament is an ordered **codec**. Dense filler arcs are
  not legal machine moves; program labels, tape, head, guards and output
  semantics remain explicit.
* Class-level and converse-merged black-core targets are separately defined.
  Pure-blue nodes are excluded, and no quotient path is automatically lifted.
* Session 64's conjugacy is on punctured odd 2-adics; positivity/integrality
  and the native Antihydra barrier do not survive automatically.

## Literature and provenance

Primary theorem statements were opened for MSS Theorem A, Bavula inverse
bounds, Adjamagbo–van den Essen Theorem 2, Belov-Kanel–Kontsevich Theorem 1,
Lagarias Theorem 1.1/reverse proof, and Aaronson section 5.10. The current
Cryptids page supplies source-reported machine verification labels only.
The original MRDP record/author metadata were found, but two direct full-page
fetches timed out; the exact Diophantine representation was then verified in
Matiyasevich–Robinson's primary paper section 3, equation (3), at
`https://arxiv.org/pdf/0802.1052`. No fresh formal build is claimed.

Eleven exact repository source snapshots are hash-pinned in SOURCE_PINS.json.
The original worktree was read only; source statuses and historical warnings
were not changed. The parent owns integration, independent audit and Git
publication. No unresolved tool-permission denial or mathematical execution
blocker remains in this lane.
