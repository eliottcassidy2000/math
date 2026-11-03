# Audit and reproduction

Status: **same-agent proof audit with independently implemented finite checks**;
no independent researcher review, no Lean formalization, no all-size minimality
claim. All validations use exceptions rather than optimization-sensitive asserts.

## Reproduce

From the repository root:

```text
python -B 05-knowledge/results/square_birth_windows_20261103/verify.py
python -B 05-knowledge/results/square_birth_windows_20261103/incoming_check.py
python -B 05-knowledge/results/square_birth_windows_20261103/seal.py --check
```

The first command runs both complete calculations normally and with `-O`,
requires identical generated artifacts and stdout, and updates `checks.json`.
This includes actual controller tables, not just the printed counts. The seal
uses LF-normalized text hashes; `checks.json` also records platform-native
artifact hashes from the replay. The literal 56-path is an input certificate,
validated for all vertices and every edge, not an asserted search result.

## Universes and independent paths

- Birth universes: ALL orientations of Q25,Q26,Q27, respectively 2^32,
  2^34,2^36. Active sources are complete unions of path-prescription cubes.
  The remaining cube is handled by the proved residual reduction, including
  old inactive sources which may become active later.
- Future relevance is also checked by projecting each future path pattern
  first and enumerating only its old-edge cube. This independently reproduces
  both support sets without completing free future coordinates.
- Path-bank exhaustiveness: forced endpoint 18; exhaustive traversal and
  both directions. The growth checker independently constructs graph edges
  from square roots and paths with set-based traversal; the birth compiler
  uses pair testing and bitmask traversal.
- Minimality: forward signature refinement versus separately implemented
  inverse-image partition splitting. Equality is checked as a partition of
  the full explicit domain, not merely equality of state counts. The symbolic
  residual has a proved uniform behavior and exact cardinality.
- Closure: every explicit member agrees on all outputs, old command targets,
  supplied-bit birth targets and converse targets. Exact source masses sum
  to the whole cube, and every target mass is recovered by summing all source
  masses over all supplied new-bit words. No orientation averaging resets a
  source between commands.
- Hostile: a one-edge pair equal in the complete fixed-26 observer but
  different after the identical 27 birth. The differing path legalities
  are replayed by literal arc membership, without the mask-matching routine.
- Positive inherited controls: fixed counts 33/71/637, active source counts
  4864/11232/62952, and the full 25 -> 26 count 97. No three-run matrix or
  81/27 transporter result is recomputed or extrapolated.
- Growth: direct pairs versus square-root edge construction at every n<=300;
  literal parents versus the inherited floor formula; incoming cumulative
  formula versus direct edges. Complete small Hamilton bank through 28.
- Orientation: topological elimination on ALL active orientations at 25,26,27
  gives 20/24/70 acyclic ones, and actual complement of the forward path
  orientation fails acyclicity as the general proposition predicts.
- Arithmetic controls: all 25 prime-value prefix graphs, Ogg bipartition,
  the user's induced components, exact 55/56 edge counts and the 56-path.

## Proof audit

The future-support reduction works for every finite nested guard family, not
just square sums. It relies on the selected path surviving K_P and on births
retaining ALL old bits. Losing either property would invalidate the residual
state argument. A residual of zero size is omitted in the general statement;
the present computations check it is nonzero.

Refinement proves necessity through distinguishing permitted words and
sufficiency through stable transitions. Converse is an additional allowed
operation because it descends and conjugates every named command, not because
it has been silently identified with a path complement. No graph relabeling
or vertex isomorphism quotient is taken.

The majority threshold has an all-size monotonicity proof, with finite exact
checks only for the first 56 entries and a nonvacuous 56-path. The square-root
asymptotic is credited to the incoming stronger formula. Equal cardinality
with 15 Ogg primes is not a map: the explicit 6/9 bipartition is a hostile
to a proposed Hamilton-path-preserving identification.

## Remaining boundaries

- 761 and 903 are minima for the named-path, complement, converse and supplied
  birth-bit language at the specified horizons. They do not establish an
  unnamed, count-only, arbitrary-re-marking or infinite-horizon minimum.
- Class cardinalities transport uniform source counts exactly. A prescribed
  correlated source distribution requires its own pushforward; no claim of
  independent joint factors is made.
- Guarded Hamilton paths differ from unrestricted tournament Hamilton paths.
  No new A000568 value, solution of the 7/21 gap conjecture, or moonshine
  equivalence follows from these calculations.
- The first25-prime graph is a finite object; its isolated primes do not stay
  isolated when larger prime vertices are admitted. No all-prime obstruction
  follows from that finite disconnection.
- No mathematical error requiring a new MISTAKES entry was found. Incoming
  overlap is credited in RECOVERY.md; frozen earlier packets were preserved.

Repository `agents/check_docs.py` passed before publication, including bounded
startup, status/frontmatter, correction sentinels and router smoke. Explicit
path staging and a clean rebase preserve concurrent packets.

Post-rebase portability repair: incoming Session148 `alignment.csv` had
26 CRLF records and was byte-for-byte identical to HEAD, but the repository's
text-cleaning rule reported it dirty. An exact-path `-text` exception in root
`.gitattributes` preserves that frozen CSV's audited bytes. Its companion
`paths.json` initially failed the incoming frozen-byte checker because Windows
checkout expanded LF to CRLF. After verifying equality up to line endings,
that file was restored to the exact HEAD bytes and given an exact-path LF
attribute. No incoming record, path witness or mathematical claim was changed
or staged; the two attributes preserve the existing manifest conventions.
The incoming wrapper still compares platform-native stdout bytes to its LF
hash and therefore rejects Windows CRLF stdout. `incoming_check.py` checks
every artifact against that unchanged manifest, then replays its verifier
normally and with `-O`, normalizing ONLY stdout newlines before comparison.
Both match the published output hash. This portability limitation and its
exact normalized checks are retained in `incoming_check.json`; the frozen
incoming wrapper itself was not edited.
