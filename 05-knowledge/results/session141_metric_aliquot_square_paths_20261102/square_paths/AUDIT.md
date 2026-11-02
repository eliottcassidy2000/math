# Verification and scope

## Proof review targets

The new mathematical content is the common-block normal form (including
both directions of its singleton-conflict criterion), the `d=1` cycle
boundary, and source-preserving directed block blow-up. "New" here means
this continuation beyond the explicitly inspected pinned sources; no
literature-priority claim is made. The order-25 formula, general indicator
expansion and one swap witness are inherited and credited.

The first cubic response at 26 is a finite exact result. Its complete
24-path bank, all twelve compatible pairs and both triples are retained
in `normal.json`; canonical block decompositions certify their union arc
counts. Order 25 retains its full 20-path bank and four disjoint pairs.
No all-tournament enumeration or probabilistic sampling is involved.

## Algorithms and controls

1. All-start DFS exhaustively visits every square-compatible permutation;
   it uses no endpoint/leaf pruning that could silently discard a case.
2. Compatible subsets are enumerated as cliques; actual arc unions yield
   the binomial moments and exact inversion to the response histogram.
3. Independently from the clique recursion, each path's free guard edges
   are oriented explicitly and the active sources counted. The unvisited
   guard orientations have response zero. Both histograms agree for all
   `n=1,...,30`; neither algorithm enumerates all `2^binom(n,2)` tournaments.
4. 7,432 pairs of distinct permutations through order five check both
   directions of the canonical block criterion, including incompatible
   examples. A singleton swap is a negative control; a two-block cyclic
   rotation is a positive `d=1` control. Three literal block blow-ups check
   the source-preserving bijection. A directed triangle/four-cycle control
   rejects dropping the source hypothesis.
5. Every histogram has total `2^m` and first moment
   `number_of_directed_paths * 2^(m-n+1)`. Normal and optimized outputs
   match byte-for-byte. `PORTABILITY.json` records a fresh run from an
   isolated directory containing ONLY `verify.py`, plus source-snapshot
   digest checks.

The two producer counting algorithms share the path enumeration. The
parent session's independent reviewer separately reconstructs the path
bank and n25/n26 subset census; that review is a separate packet, not
represented here as a second proof assistant or formal verification.

## Intended boundaries

- Full tournaments remain labelled on the fixed integer vertex set.
- The square-sum guard is retained; an even guarded response says nothing
  contrary to parity of unrestricted tournament Hamilton paths.
- The n25 matching hypothesis fails at n26. No inherited complement-state
  formula that assumes a matching is used at n26.
- The abstract blow-up theorem preserves directed paths but does not
  promise square-sum labels on an arbitrary expansion.
- The finite maximum responses decrease from five at n28 to three at n29;
  there is no monotonicity or eventual-generalization claim.
- Counts equal to unrelated divisor counts, or to `704-694`, are not a
  typed correspondence and play no role in the proofs.
