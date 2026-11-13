# Scope, proof audit and exact replay

2026-11-13, Codex. **Self-audited proof notes, independently implemented
finite controls.** No separate agent or external mathematical review was
performed. Chan--Pak and ordinary Collatz remain OPEN.

## Reproduction

From the repository root:

```text
python3 -X utf8 05-knowledge/results/prime_ternary_lift_20261113/replay.py
```

`REPLAY.json` pins23 local executable dependencies and the inherited
ordinary-class input bank. Eight scripts run in ordinary and optimized
Python, with explicit exceptions rather than removable assertions.
Eight matching receipt pairs are byte-identical. Do not silently reseal
changed dependencies when reproducing an existing claim; `--seal` is an
authoring operation. Third-party environment requirements are Python3,
NumPy and Numba; no network access or external solver is needed.
The repository documentation checker, all local packet links, all23 pinned
inputs, and the final eight receipt-pair comparisons also passed.

| receipt | explicit universe and check |
|---|---|
| surplus | all ordinary classes2..8, every vertex occurrence;307 delta1/2 occurrences,47 of indegree2 |
| frames | all strong frames3..8, every two-run port;2,477 responses,2,074 normalized rays,264 transitions,9,908 literal expansions |
| moves | three displayed frames; independent permutation-slots coefficients; actual arc edits and Hamilton counts for m1..7 |
| visits | all strong frames3..8, every port, every possible visit count up to4;22,817 responses,312 transitions,91,268 literal expansions |
| triage | attained parents3..73 with a missing direct digit and a prime/7p child; search only the stated normalized-ray catalog, retaining lower exponent thresholds |
| renewal | all85 root occurrences through5 compared against independent skeleton enumeration; six selected higher-mode lifts, m1..5 and cyclic-module hostile controls;8,190 digit words; known renewable +2 clock |
| surplus_audit |85 direct permutation extractions of F and85 literal triangle-context children;600 ordinal-block factorizations; complete small-head missing-digit lists |
| clone_parity |all20 ordinary frames through5, all subsets of doubled ports, each selected coordinate changed2->3;1,117 comparisons; cyclic three-port controls at widths1..4 |

The finite integer DP never exceeds15 vertices here; even15! fits signed
64-bit arithmetic. Symbolic power formulas and cycle-word arithmetic use
Python integers. Infinite-family results use finite coefficient identities
and the run-expansion bijection, not numerical extrapolation.

## Proof obligations checked

- **Types:** Hamilton paths are simple spanning vertex words, not matrix
  walks. K_j means ordered nonempty path covers; P(2) in the padding
  lemma allows empty slots. The two conventions are not interchanged.
- **Universe:** ordinary isomorphism classes and labelled root occurrences
  are distinguished. A response-vector quotient forgets other ports and
  endpoint statistics. No all-frame absence is inferred from order eight.
- **Construction:** exact source/target adjacency masks and maps are
  retained. Every chosen move checks tournament validity, strongness,
  Hamilton count, old module retention and the claimed old-arc budget.
- **All-parameter scope:** complete skeleton coefficients give the universal
  arbitrary-module identity. Substitution of C3 verifies that the digit
  becomes dH(M), preventing a false recursive promotion.
- **Finite-phase theorem:** the cycle argument includes the constant and
  zero-growth cases; exponential independence, factorial divisibility,
  oddness and the extremal additive term force all digits to be2. The
  known amplitude60 loop is a positive equality control.
- **Parity theorem:** the one-port proof works with any finite exterior;
  applying it coordinatewise avoids assuming a finite simultaneous
  expansion when module ports contain cycles. The corollary requires
  eventually fixed singleton status. Alternating singleton status is
  explicitly outside it.
- **Small-surplus all-order upgrade:** root position in the strong-component
  order is proved before invoking the cited Busch floor. It applies to
  the delta1/2 existence question at h<=73, not all positive surplus minima.
- **Collatz direction:** the plus branch has odd r and a literal constructor;
  the first halving has no constructor in the stated fixed-frame classes.
  Both-sign arithmetic descent is not ordinary Collatz, and a forward
  scalar descent is not an inverse graph-construction theorem.
- **Prime lanes:**49 is excluded from the infinite7p lane and remains a
  finite seed. No primality infinitude is inferred for a displayed family.

## Corrections and hostile probes retained

The initial exploratory clone expression omitted internal transitive
arcs. It was discarded and logged in `01-canon/MISTAKES.md` before any
theorem claim. The literal validity guard later caught a test fixture
that encoded a directed path rather than TT4; the fixture is now
(14,12,8,0). A small-count validation initially limited its bank to order6,
which omitted39; it now stores actual order7 witnesses for39 and79.
Neither failure was treated as a global absence result.

The direct triangle-context failure at15->47 is a proved obstruction
to that operation, using the all-order padding reduction. Failure to
match a larger frame in a finite catalog is not. The tentative bound
delta>=2^indegree-2 remains OPEN despite no counterexample through8.

## Inheritance and incoming work

The packet builds on the Session193 prime-atom reduction; the universal
run calculus and excursion symmetry; the identified-root endpoint
compiler; the factor-five theorem; established dyadic lines; the
Session304 physical fork and signed arithmetic; and the inherited
ordinary-class bank. Full file links are provided alongside each use.

Incoming commits on the63 repair and fibre-absorption programs were read
after the first checkpoint. They reinforce the distinction between
count supply and physical transport. The retained-star three-arc repair,
paired-module primitive guard and finite activation ideal are credited
in `COLLATZ_AND_TARGETS.md`; no automatic entry into our frame rules is
inferred from them.

After the second checkpoint, Session305's synthesis and Session306's
exact-observer note were read. Their physical-versus-arithmetic and
finite-state boundaries are recorded in the final section of
`COLLATZ_AND_TARGETS.md`. They are not dependencies of the proofs or
computations in this packet.

No new theorem ID, global tournament sequence term, universal seed,
all-prime selector, complete spectrum classification, or Collatz
convergence result is claimed. The remaining construction obligation
is explicit in the README and connection ledger.
