# Sevenfold selection: exact guards and the first live one-vertex lift

**Exact:** a count-5 source has a literal one-vertex extension of
count 35, and order seven is minimum for that extension. The two
sevenfold holes at 7 and 21 are inherited. Universal sevenfold
coverage and Chan–Pak are not proved.

## Two distinct operations

UNIT_VERTEX_RETURN.md says the only one-vertex gains of two are
1→3 and 3→5. On the sevenfold axis the same three inputs ask for
7,21,35; their first two values are excluded and the third occurs.
This is alignment of the first available states, not transport of
a +2 operation by a sevenfold multiplier. In particular it is
conjectural that no other positive odd input gives a sevenfold hole.
On the tripling ray q=3^k there are exactly two holes, since the known
63 seed gives 63·3^j by ordinal sums with directed triangles.

There are two distinct ways a sevenfold selector could fail:
1. **Source/interface failure:** a given count-p tournament exposes
   no suitable attachment, even if 7p has a witness elsewhere.
2. **Exact target failure:** a proposed approximate or modular
   repair has not actually landed at 7p.

The e-means-plus study separates two failures of its heuristic:
a seed-loss check protects against deterioration relative to the
seed, and a pseudocount floor controls its forgetting counterexample.
The analogy here is architectural, not a transfer of that theorem:
maintain alternatives keyed by their full marked responses rather
than averaging them into one source (a diversity safeguard); accept
only after the full exact path-count test (an independent landing
safeguard). If the test fails, switch source or repair class, not
just the numerical tolerance. This neither proves the search
terminates with a witness nor excludes another sevenfold hole.
Source: https://github.com/rohanarun/emeans-plus .
The accompanying Zenodo link did not resolve in this pass.

## An exact physical guard

For a tournament U with added vertex v, the inherited odd-cycle
deletion identity gives

    H(U) = H(U−v) + 2 sum_(odd directed cycles C through v) H(U−V(C)).

The sum counts cycles up to cyclic rotation. Thus if the retained
source has count p, landing at 7p requires the exact weighted
sum 3p. Keeping only the number of marked triangles, or the
sum modulo a base, is not sufficient as a general certificate.

Literal seven-vertex example, outgoing row masks on 0,…,6:

    U   = [120,113,115,54,32,64,24]
    U−6 = [56,49,51,54,32,0].

The parent is the strong four-vertex count-5 core followed by
two sinks in transitive order. An all-labels scan verifies that
no output of order <=6 admits this 5→35 deletion pair.

| cycle length through vertex 6 | number of cycles | sum of complement counts |
|---|---:|---:|
| 3 | 4 | 8 |
| 5 | 6 | 6 |
| 7 | 1 | 1 |

The gain is exactly 2(8+6+1)=30, and subset DP independently
verifies both total counts. The two longer cycle lengths are
essential to the accounting on this literal witness, not
asserted necessary in every possible lift.

The marked 5→35 witness immediately gives one-vertex bridges
5h→35h for every realized h: append the same H=h tournament as
an ordinal summand to both sides, so deleting v still recovers
the source, and both path counts multiply by h. In particular
5·3^k→35·3^k holds physically for all k>=0. This continuation
does not place arbitrary primes: the scalar source now
contains the factor 5.

## Selector experiment suggested by the guard

For prime inputs, construct a bank of physical count-p sources
with *different* marked weighted-cycle responses; seek a proof
one response lands at 3p, or fall back to a different construction
if no one-vertex extension works. The test identifies the exact
resource to preserve, not a numerical proxy for it. The successes
at input 5 and the tripling ray do not establish that 1 and 3
are the only missing sevenfold inputs.

Run python3 sevenfold_guard_replay.py for the literal gains,
cycle contributions, and classification of the minimal output order
for an unchanged-old-arcs 5→35 pair.
