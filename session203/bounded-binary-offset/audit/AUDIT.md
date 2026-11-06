# Independent audit of bounded binary offsets

Verdict: **PASS in the stated scope**. No correction required. This audits
`session203-occupancy-tail/PROOF.md` at SHA-256
`ccc0f66fbdfc2ed19b5767b8e12ffbbc02dfc7e81d1c3c7ac0307bf21447bf8e`.
The immutable input copy is `PROOF-frozen.md`.

## Main quantifiers and mathematical audit

The theorem counts parent scalar values having at least one actual source
representative and at least one successful append word. It does not fix a
preferred representative. Consequently its union over source profiles really
does include arbitrary target-preserving numerical-source rebasing. It does
not cover rebuilding with a different parent value, deleting source letters,
or an arbitrary tournament operation.

All identities in Section 1 were reconstructed in positive coordinates and
checked on the four coordinate basis vectors. In particular `JZ=OJ` is a
numerical profile identity, not a graph isomorphism assertion. The lower
readout factor 3/2 holds coefficientwise; after two letters the factor is at
least 9/4. The strict hypothesis n>4C is exactly what excludes two or more
appended letters for a target 2n+r with |r|<=C.

The bounded-D proof strips only O, which preserves D. The bounded-c proof
strips only Z, which preserves c. The remaining J case uses the bounded-D
lemma and rewrites all trailing Z letters across J. In the ell argument a
Z output is itself a bounded core, while a J output reduces to bounded c.
For h<0, stripping Z leaves a seed or J, never O. For h>0 it leaves O;
the E=0, preceding-Z, and preceding-J cases are exhaustive. The latter uses
the vanishing f coordinate at a J output. These steps prove containment in
the six suffix families without assuming their converse.

Each normal representative is still an actual word. The universal F growth
therefore bounds its total length, not merely its original source recipe.
There are three families with two free exponents and three with one. Counting
all choices, including redundant representations and representatives whose
F exceeds N, is a valid upper bound. This proves the claimed exact
`6 b_C (L_N+1)^2`, hence O_C(log^2 N).

Section 4.1 is also valid. Rewriting JZ to OJ preserves length and profile,
and strictly decreases the total number of Z letters; termination needs no
confluence. The number of JZ-avoiding words of length l is Fib_(2l+2), and
the cumulative count is Fib_(2L+3)-1. Since F<=2U, this yields the stated
polynomial core bound with an absolute implied constant. For N a multiple
of four, (N/2,N] has exactly N/4 odd integers. The condition C<N/8 puts all
of them above 4C. Applying the parent bound gives the claimed necessary
growing-offset lower bound. This is an interval of parents each admitting
some move, not an assertion about all possible children of a sparse parent
set; the proof keeps that distinction correctly.

## Independent computation

`check.py` imports no author verifier. It uses a separate positive-coordinate
matrix implementation and enumerates the complete orbit through U<=1024:
12,100 profiles. For eight C values it constructs the union of all six
normal-form families from every bounded core, rather than calling the
author's normalization algorithm. Every qualifying profile with F<=1024
is in that union. Scalar-first unions agree with the direct append condition
after retaining every profile in each source fiber.

Additional checks:

- all 9,841 recipes of lengths zero through eight satisfy the length bound
  and physical sidecar parity/positivity;
- the complete parent-25 fiber has two profiles and minimum absolute
  one-letter offset five;
- an independent two-state bigram automaton verifies Fibonacci counts and
  cumulative counts through length sixty;
- all 256 bounded-core census checks with 1<=C<=256 satisfy Section 4.1's
  exact Fibonacci bound.

Run `python3 check.py > NORMAL.json` and `python3 -O check.py > OPTIMIZED.json`.
The two deterministic outputs agree byte for byte. All check gates use
exceptions, not removable assertions. These are finite implementation
controls; the all-order conclusion rests on the audited proof.

No eventual Hamilton-path spectrum threshold, occupancy assertion, or
all-order absence of a target is proved by this result.
