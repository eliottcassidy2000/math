# The 223 recovery: private certificates and shared Boolean ownership

Status: **PROVED, root-audited, exact computation reproduced**. This note proves a
finite and asymptotic compatibility statement, not Collatz, an LRC result, or
a lower bound comparing global and private edit prices.

## Inheritance and corrected target

The closest proved mechanism is the requirement to retain joint guard state:
[THM-2204](../../01-canon/theorems/THM-2204-scalar-depth-223-thirteen-lift-capacity-law.md)
preserves the total capacity of the thirteen lifts of depth profile `(2,2,3)`
but does not identify the location of a capacity deficit. Its least-used
sidecar is the full labelled guard/hole correlation vector.
[THM-2233](../../01-canon/theorems/THM-2233-guard-danger-hidden-state-bellman-profile-exclusion.md)
retains a guard-parent bit in the Bellman state, while
[THM-2232](../../01-canon/theorems/THM-2232-same-core-signed-eigen-markov-dual-exclusion.md)
and [THM-2239](../../01-canon/theorems/THM-2239-unrestricted-multicore-signed-dual-profile-exclusion.md)
retain signed residuals that unsigned moments discard. These are conditional
LRC profile exclusions, not LRC(14).

The other literal `223` is a different object:
[THM-3174](../../01-canon/theorems/THM-3174-projected-k3-z223-terminal-descent-and-cap222.md)
closes the projected `z1=223` atlas layer and obtains a cap of 222. Neither
this integer nor the depth profile has an identified arithmetic map to the
Collatz starting value 223.

The corrected near miss is private-to-global gluing.
[THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md)
and the [pairpeak note](procgen_pairpeak_20260926_pairing_peak_price.md)
refute the old polynomial comparison to the unweighted bad-word population
for arbitrary edits and for private pairing prices. Private pairing price
has a stretched-exponential improvement. One globally consistent pairing is
a different problem. The canonical hostile example developed below is the
pair of starting values `7,11`: their unique cheapest private certificates
require opposite values of the same saved bit.

Live concept board:

1. Guard state becomes ownership of the actual pair bit shared by sources.
2. A signed residual becomes the opposite load changes on the two members
   of a flipped pair; cancellation needs their joint weights.
3. Private trajectories become Boolean clauses, with conflict graphs before
   any attempt to construct a global map.
4. The 223 source becomes an explicit short conflict chain, explained by a
   2-adic valuation rather than by numerical coincidence.
5. Peak-discounted prices remain the correct private benchmark; compatibility
   is an additional coordinate and does not reverse their established gain.

The relevant predecessor is section 2.4 of the
[pairing-transitions note](procgen_brackets_20260924_pairings_transitions.md):
it already derives the global two-step landing clauses, the multiplicative
path graph, and a `v3` vertex cover for a weaker necessary condition. The
new statement here identifies the *unique privately cheapest certificates*
and gives their exact finite compatibility optimum. It is not a claim to
have discovered the multiplicative graph or the two-step clauses anew.

## Exact private-certificate compatibility theorem

Partition the positive integers into pairs `P_i={2i-1,2i}`. A bit
`epsilon_i` selects the pairing map

```
epsilon_i=0:  2i-1 -> 3i-1,   2i -> i;
epsilon_i=1:  2i-1 -> i-1,    2i -> 3i.
```

The zero-bit map is the shortcut Collatz map. For a source `n`, a private
two-step certificate specifies the bits encountered before the first value
strictly below `n`, which must occur in at most two steps. Its private cost
is `sum(n/i)` over its flipped pair indices. Specifying additional unused
positive-cost flips cannot improve this certificate.

**Proposition.** For every even `i>=2`, the source `n=2i-1` has the unique
cheapest private two-step certificate

```
epsilon_i=0, epsilon_(3i/2)=1,
2i-1 -> 3i-1 -> 3i/2-1,
cost = (2i-1)/(3i/2).
```

For sources whose even pair index is at most `X`, the maximum number of
these unique optimum certificates that a single bit assignment can realize
simultaneously is exactly

```
floor(X/2) - R(X),
R(X) = sum_(j>=0) [floor(X/(6*3^(2j))) - floor(X/(18*3^(2j)))].
```

Consequently `R(X)=X/8+O(log(X+2))`: exactly one quarter asymptotically of
these privately optimal certificates must be abandoned. Here abandonment
means changing the certificate, not failing to descend and not paying an
additional globally distinct flip.

**Proof of the private optimum.** Flipping the source pair produces
`i-1<n` immediately and costs `n/i`. Leaving it unflipped produces `3i-1`,
an odd integer at least `n`. Its pair index is `j=3i/2`. Flipping `j`
produces `j-1<n`; leaving `j` unflipped does not descend. The second option
costs `n/j<n/i`, and these exhaust the possible paths through two steps.

**Proof of exact compatibility.** Two optimum certificates conflict if and
only if the zero-bit index of one is the one-bit index of the other. Thus
the conflict graph has vertices the even integers `i<=X` and edges
`{4r,6r}` for `6r<=X`. These are paths: on a component the integer
`i=2^a 3^b u`, with `gcd(u,6)=1` and `a>=1`, moves from `(a,b)` to
`(a-1,b+1)`.

The set `R={i even: i<=X, v3(i) odd}` is a vertex cover, since the two ends
of every edge have opposite valuation parity. The edges

```
M = {{4r,6r}: 6r<=X, v3(r) even}
```

form a matching. Two such edges sharing an endpoint would force the
3-adic valuations of their parameters to differ by one. Each point of `R`
is the right endpoint of exactly one edge in `M`, and conversely. Hence
`|M|=|R|`, proving that this cover is minimum without any stationarity,
density, LP-duality, or scaling assumption. Counting even integers of odd
3-adic valuation gives the displayed formula. Summing its geometric main
terms gives `X/8`; only `O(log X)` floor terms are nonzero.

An independent family is obtained by keeping the vertices with even `v3`.
It is realized by the one assignment `epsilon_i=v3(i) mod 2` for all `i`.
This proves attainability of the compatibility optimum. It does **not**
assert global two-step descent: this assignment sends `6 -> 9 -> 14`.
The previously good sources therefore supply additional global clauses.

**Corollary (an exact compatibility penalty for a specified functional).**
For a common bit assignment that rescues every `n=3 mod 4` within two steps,
let its executed cost on such an `n` be the sum of `n/i` over flips actually
encountered up to first descent. Average this cost over all source integers,
putting zero on the other three residue classes. The infimum of the lower
asymptotic mean, over these assignments, is exactly `3/8`. By comparison,
the unrestricted individual private optimum has mean `1/3`.

Indeed, each such bad source has exactly the two options in the proof. An
abandoned optimum must flip its own pair and pay the extra
`(2i-1)/(3i)=2/3-1/(3i)`. At least `R(X)` optima are abandoned for even
`i<=X`. The harmonic correction is `O(log X)`. The total individual optimum
is `(2/3)X+O(log X)`, and the extra cost is at least `X/12+O(log X)`.
Dividing by the source cutoff `2X` gives `1/3+1/24=3/8`. The `v3` parity
assignment rescues every bad source: at odd `v3(i)` it flips immediately,
and otherwise it follows the cheapest certificate. Its abandoned set has
size `R(X)`, so it attains the bound. No density assumption is needed for
the lower bound. The `9/8` penalty is for this executed, source-weighted
functional. Shared flips are counted once per source that uses them, and
previously good sources are excluded; hence it is **not** a `9/8` lower
bound for global flip density divided by private price.

## The actual 223 experiment

For source 223 the consecutive conflicting private optima are

| Source | Its pair | Cheapest flipped pair | Private path |
|---:|---:|---:|---|
| 223 | 112 | 168 | 223 -> 335 -> 167 |
| 335 | 168 | 252 | 335 -> 503 -> 251 |
| 503 | 252 | 378 | 503 -> 755 -> 377 |
| 755 | 378 | 567 | 755 -> 1133 -> 566 |

Every two adjacent rows prescribe opposite bits at the shared pair.
The `v3` cover retains the first and third optima. The chain has four
vertices because `v2((223+1)/2)=4`; after four upward source transitions
the next source is `1133=1 mod 4`. No property of 223 being prime is used.

For all bad sources `n<=223`, the even pair indices satisfy `i<=112`:
there are 56 optimum certificates and the exact minimum abandoned count
is 14. For the different cutoff `i<=223`, there are 111 certificates and
the count is 28. These two universes must not be interchanged.

The smallest conflicting pair is `7,11`: source 7 prefers bits
`epsilon_4=0,epsilon_6=1`, whereas source 11 prefers
`epsilon_6=0,epsilon_9=1`. Yet the single flip set `{6}` rescues both:
`7 -> 11 -> 5` and `11 -> 5`. Source 11 loses its private optimum, but
the repair adds no distinct global flip. Thus the proposition gives no
constant-factor lower bound for global price divided by private price.

## What the transfer preserves, and what it loses

Source: labelled guard-capacity or guard-state Bellman data.
Target: labelled private Collatz rescue certificates.
Map: replace each local capacity choice by its constraints on the actual
shared pair indices. Preserved predicate: simultaneous consistency of
specified local choices. Lost by summing private costs: the identity and
sign of the Boolean requirement, overlap of charged flips, and constraints
created for previously descending sources. Needed sidecar: a signed
source-by-pair incidence relation, plus the selected trajectory itself.
The cheap decisive test is two private minima that prescribe opposite bits.

For one flipped pair the two changes in image are exactly `-2i,+2i`, and
the pair image sum remains `4i-1`. This signed neutrality helps only if the
actual source weights put the required joint mass on both endpoints.
It does not imply cancellation along a positive orbit. That is the same
logical loss that makes scalar capacity totals insufficient to locate a
guard hole; it is not an arithmetic identification of the two systems.

## Reproduction, controls, and next boundary

Run `python3 04-computation/experiments/crossroads223_20260926_flow.py`, or
the same command with `-O`. The retained output is
[crossroads223_20260926_flow.out](crossroads223_20260926_flow.out).
Normal and optimized outputs are identical; explicit checks remain enabled.

The experiment enumerates all private branch assignments for 1,000 sources
`3<=n<=3999`, checks 1,000 incompatible pairs, and compares cover/matching
counts with independent exhaustive search on every connected path for all
cutoffs `2<=X<=300`. The `X=10^6` value is the exact floor-sum formula,
not an exhaustive million-source search. Positive controls realize every
selected even-valuation policy through source 9,999 and verify immediate
rescue for the other bad sources. Exact rational sums verify the trace-cost
gap identity at pair cutoffs 112, 223, and 1,000; printed means are decimal
renderings of those exact rational sums. The hostile control
`6 -> 9 -> 14` ensures that compatibility is not mistaken for global descent.

**OPEN, narrower question.** At longer horizons, retain the actual signed
pair incidence of near-cheapest private certificates and their alternatives.
Can changing a conflicting certificate preserve the peak discount while
controlling the newly affected partners? Even at horizon two a count of
abandoned optima does not measure that repair cost. This is the exact
missing implication; a scalar price comparison or a private success count
does not supply it.
