# Validity and reproduction audit

Status: **VERIFIED checkpoint**, with scoped proofs and finite exact computations.
This is a same-agent audit using independent mathematical and computational
paths. No external referee, separate-agent audit, or literature-priority claim
is implied. The inherited 65-symbol bank is fixed, not reselected after testing.

## Proof obligations and their boundaries

| Claim | Mechanism | Boundary checked |
|---|---|---|
| THEORY T2: common 3H-10 and 3H-20 families | An exact row-error bridge, a conserved charge, and a terminal triangle | These replace an initial source block; they do not retain every old edge |
| THEORY T3: no append-only return from odd-first primitive errors | First-column parity and divisibility of row content | Does not exclude re-marking, deletion, or source replacement |
| THEORY T4 / ARITHMETIC A1-A4: loop laws and minimal arithmetic carriers | Exact integer matrices, unique macro parsing, valuation cone, separating contexts | The six-macro language is specified; its gain is not H |
| ARITHMETIC A5: three linear Hamilton coordinates | A nonzero actual-source Hankel minor and the inherited upper bound | The modular quotient is on ambient integer states, not every graph fibre |
| PORT_GUARD P1-P3: minimal extra legality coordinate | A bijection between paths and two-piece covers, plus explicit separating sources | Observes only H and start-at-port feasibility/count under these joins |
| BIRTHS B2: degree at most min(2r,n+r-1) | Each path uses at most two incident arcs at each new vertex | Exhaustive evaluation is still exponential in the new-edge count |

The all-order claims do not follow from extrapolating the tests. Their proofs
are the displayed identities, bijections, and invariant arguments. The finite
tests independently check definitions, orientation conventions, and named
witnesses. No result proves full prime/7p coverage, the unrestricted {7,21}
spectrum conjecture, a complete metagraph carrier, or A000568 enumeration.

## Exact universes and controls

`reproduction.json` records a successful complete replay against the pinned
input hashes (including incoming integration sources). Every Python probe is run normally and with `-O`; explicit `need`
checks remain active, and the LF-normalized JSON outputs agree.

- **channels.py:** all 65 atoms by permutations; three inherited seeds by
  subset DP; all three new families for k=0,1,2,3 by literal adjacency
  endpoint-forest DP, reaching 38 vertices. This path count does not read
  the three-coordinate matrix output. The two equal-gain/equal-clock but
  different-H loop orders are checked by the same literal method.
- **Error census:** all words of at most five bank letters, merging identical
  signed rows at each depth. Nonnegative matrices justify pruning one-sided
  rows. Mixed-row counts are 20,113,539,2462,9570. No all-length classification
  of the possible gain primes is inferred.
- **New recurrence clocks:** exact pair iteration modulo 7,13,169,2197,91,1183;
  the new seed satisfies the inherited all-power lifting proof's hypotheses.
  Trial division independently certifies the two named primes. The composite
  quotient at k=2 passes both small-prime index guards and is a negative
  control against confusing these guards with primality.
- **arithmetic_probe.py:** exact macro-word multiplicities through added
  vertex cost 48, compared with rational-series coefficients for both sign
  evaluations. All seven Mobius states and all eight joint states have
  explicit distinguishing responses. Nine literal graph counts verify the
  Hamilton Hankel minor.
- **port_guard.py:** 1,270 joins: every pointed tournament through order four
  with four selected blocks, plus all 65 blocks on an actual zero-start and
  positive-start source. A permutation two-cover computation supplies kappa;
  independent subset DP supplies the joined graph's start count. The sharp
  rank-four witness consists of actual sources, not formal rows alone.
- **context_entry.py:** all absent states in the physical cone with odd
  h<=15 and 0<=x,y<=2h, each followed by all 65 letters: 942 inputs and 61,230
  tests, with 49 entries into the observer image. The exact inverse tests
  all word lengths for each state. A literal four/five/four-vertex control
  and the complete H<=3 word census distinguish image entry from graph entry.
- **two_birth.cpp / birth_verify.py:** all 196,608 two-birth extensions of
  the six old H55 presentations; all 33,554,432 three-birth extensions of
  their two ordinary classes. Explicit isomorphisms transport the result
  to all six presentations. Every monomial/zeta output equals its separate
  Walsh/Hadamard output. The common path-constraint input is independently
  tested by subset DP on 585 literal extensions, comparing the exact value
  at the same orientation in both two- and three-birth tests. Python also reconstructs
  all two-birth outputs from sparse coefficients. The external H55 source
  supplies 4 and 18 positive witnesses, respectively.
- **81/27 hostile:** the actual adjacency rows reproduce 67 identical cyclic
  triples and H=462915, but the specified sink/port extension gives 2390013
  versus 1754919. Automorphism orders are recomputed, not taken from labels.

- **integration_probe.py:** exact quotient intertwining, determinant and
  characteristic polynomial for D=L2 b L2; all eight nonzero normalized F9
  phases and their source-relative Frobenius map; literal graph counts for
  all three seed families at k=0,1,2,3 through 66 vertices. The uniform level
  step is proved by invertibility of R modulo 3. The universal parity
  extension and the exact index-100 return lattice are proved separately;
  the three original seeds fail the integral reset gate. The universal parity
  extension credits the incoming all-block calibration. Other incoming
  censuses are cited with their scope, not claimed as fresh independent checks.

## Reproduce and inspect

From this directory run:

```text
python -X utf8 -B reproduce.py
python -X utf8 -B seal.py
```

Python needs SymPy; the extension census needs a C++17 compiler. The reproducer
locates the compiler and deletes its own temporary executable. The largest
three-birth transform uses about 400 MB of arrays. Runtime is measured in the
receipt, not asserted as a platform-independent bound. The file manifest uses
LF normalization to survive Git's Windows newline conversion.

The inherited inverse is a terminating decision procedure for a specified
state, not a proof that a prescribed scalar target has a successful state.
All four next ternary states still fail after the new bridge. The exact
zero-start invariant is a permanent obstruction to one tempting source
substitution. These failures are retained alongside the positive families.
