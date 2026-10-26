# Source receipt, inheritance, and computation scope

Read on 26 October 2026. The mathematical note is an elementary exact
specialization, with no historical-priority or analytic-improvement claim.

## Classical attribution

Adam J. Harper, *Lecture Notes 1 for Cambridge Part III Course on Elementary
Methods in Analytic Number Theory*, Lent 2015, dated 5 February 2015,
[author-hosted PDF](https://warwick.ac.uk/fac/sci/maths/people/staff/harper/elementarymethodsnotes1.pdf).
Read Section 1, Definition 1.1, Proposition 1.2 (basic Eratosthenes–Legendre
version), and its proof on printed pages 1–3. The source gives the Möbius
indicator for avoiding a fixed set of prime factors and its exact floor sum;
it also explains the square-root cutoff and why a rapidly growing remainder
prevents a useful bound from following automatically. These are the classical
mechanisms used here. All two-form identities, equality exceptions, and
integer-height claims needed by this note are proved explicitly in README.
No Selberg-sieve estimate, parity-barrier theorem, or distributional claim is
imported as an unstated dependency.

## Repository inheritance

- `twin_pythagorean_124_20261026/tree/RESULT.md`: the exact integer T(x) map,
  D-ray address, primality restriction, and the line-graph interpretation.
  Those results are inherited, not reintroduced as new findings.
- `twin_pythagorean_124_20261026/residue_obstruction/README.md`: every fixed
  residue observer has both-composite false positives. The new height gate
  explains the additional datum needed by a growing observer.
- `prime_pairs_local_symmetry_20261019/prime-pairs-structural-note.md`, §§2,
  5: fixed twin-center wheels, explicit warning about self-prime divisibility,
  and the distinction between roughness and primality. The present exact
  self-hit partition makes that warning operational for a single pair.
- `prime_pairs_local_symmetry_20261019/power-depth-rough-pairs.md`: already
  supplies a stronger analytic lane for four rough linear forms and their
  almost-prime count. Its carrier and target differ from this two-form finite
  count, and no implication upgrades its roughness conclusion to primality.
- `THM-2447-twin-gap-local-prime-harmonics-and-detrended-spectrum.md`:
  local roots and their exact multiplicities are distinct from measured
  global prime spectra.

The research protocol, active guardrails, CORE-PAPERS and its analytic-number-
theory sidecar were read before choosing this lane. A targeted mistakes search
recovered the warning against calling a weaker-sieve survivor globally
admissible (MISTAKE-427), and the requirement for a controlled sieve tail rather
than finite CRT compatibility. The prior four-form note was recovered before
deriving the exact two-form count, to avoid presenting CRT wheel counts as a
new analytic sieve result.

## Exact finite universes

The companion `check.py` imports no producer and uses only the standard
library, with explicit runtime gates rather than Python assertions.

- Every least prime factor through 10002 is compared with direct trial
  division. Prime pairs are counted by their lower member through 10000.
- Strict gcd counts are computed for every height through 10000 and every
  prime cutoff through 101.
- Complete wheel lift sets are compared with literal enumeration through
  modulus 30030 (cutoff 13). Allowed-residue floors and full Möbius-root sums
  are compared on the declared small interval and selected larger heights.
- The proper-factor exception test is evaluated independently for all
  x=3..200 and B=2..101, and its exact disjoint patch is checked.
- The square-root recurrence is checked at every height 3..10000; the
  next-unsieved-prime form is also checked throughout each available box.
- Every first-false candidate is searched through 5000 for prime cutoffs
  through 61. A found witness and rejection of every smaller candidate
  certify minimality in that stage. The nonsquare failure at 61 is independently
  factored, and its companion is tested by trial division.

Reproduce from this directory:

```bash
python3 -B check.py > replay.json
python3 -B -O check.py > replay.optimized.json
```

Both must be byte-identical to the stored receipt. The Möbius formula is a
finite exact identity; its exponentially growing number of root contributions
is not suppressed in a claimed running-time or asymptotic bound.

Independent audit [PASS](../audit_sieve/REPORT.md) is complete. The all-parameter proof and the finite experiment
remain separately typed. Exact source and artifact hashes are in
`SOURCE_MANIFEST.json`.
