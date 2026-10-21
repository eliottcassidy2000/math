# Degree-local odd-factorial transfers for tournament Burnside enumeration

Status: PROVED arithmetic identity; FINITE-EXACT controls and same-host benchmark receipt. Independent review is reported separately. This packet is new executable work, not a restoration of the lost Session21 packed-prime source or timings.

## Inheritance and live board

Closest available mechanism: incoming `tournament_enumeration_20261021a/check.py`, frozen verbatim as `baseline.py`, uses descending odd cycle lengths, packed future-gcd responses and a common odd(N!) scalar. The hostile is dropping an unaffordable coordinate before proving it can never return; the existing descending-order/remaining-degree mask is preserved exactly. The corrected near miss is the prior lost/reconstructed prime-order hybrid: its old evidence is not claimed here. The least-used sidecar is the degree-dependent denominator lattice.

The live concepts are (1) cycle-index histories, (2) future-response quotient, (3) odd factorial denominator lattice, (4) reusable integral transfers, and (5) pair-free-core substitution. The anchor is faster exact A000568; the niche is signed pair-free-core coefficients; the wildcard is degree-local rather than ambient normalization. The useful connection is between the denominator lattice and integer permutation-class counts. It changes arithmetic work, not the response equivalence relation or transition graph.

## The exact identity

Let O(d) be the odd part of d!, and let a state have used degree d and rational weight w. Store W=O(d)w rather than O(N)w. At descending odd length ell, select k actual cycles, and in core mode an excess e from substituting x^ell/(1+x^ell). Put q=d+ell(k+e). The original rational transition is

  w' = w · 2^[k(ell−1)/2 + ell*k(k−1)/2 + k*c] / (ell^k*k!)

multiplied by (−1)^e binomial(k+e−1,e), where c is the retained response at ell. For k=0 only e=0 is permitted. Write k!=2^v O(k). Then

  W' = (W*T) << [k(ell−1)/2 + ell*k(k−1)/2 − v + k*c],
  T = O(q)/(O(d)*ell^k*O(k)) · (−1)^e binomial(k+e−1,e).

The rational-looking transfer T is always an integer. Indeed

  q!/(d!*ell^k*k!)
  = binomial(q,ell*k) · [(ell*k)!/(ell^k*k!)] · (q−ell*k)!/d!

is a product of integers: the middle factor counts permutations consisting of k cycles of length ell, and q−ell*k=d+ell*e≥d. Taking odd parts of numerator and denominator proves the divisibility defining T. The nonnegative shift is the inherited factorial-two cancellation; v2(k!)≤k(k−1)/2, including k=0,1. The core sign does not affect integrality.

Thus transfers depend only on (stage,degree,k,e), not response. Precompute and verify each integer division once, then use multiply/shift/add on every response-state contribution. The previous engine performs a checked large-integer division at every (state,k), before any core-excess expansion. The new one uses no divisions on response-state contributions. Precomputation includes divisions and is included in the measured time.

Induction from W=1 at d=0 proves exact equality to the rational series throughout; merging is valid because all histories merged at degree d have the same O(d). At termination divide degree d by O(d), not O(N). For unrestricted raw histories O(d) remains the least common odd scalar because the all-fixed-point history realizes denominator O(d). No claim of this minimality is made for no-fixed-point restrictions or after signed cancellation.

## Quotient/loss ledger

Source: partial odd-cycle Burnside histories with total degree and future gcd sums. Target: (degree,packed response) and locally scaled aggregate weight. Map: retain exactly the incoming packed responses and sum O(degree)-scaled weights. Preserved: every future transition, all four coefficient modes, stage and degree. Lost: individual cycle partitions, individual tournament objects, automorphism identities, metagraph ports and adjacencies. Sidecar: bound, mode, stage, packing width, remaining-degree masks, exact normalization O(d). This is a counting acceleration, not a classification or smaller response-state space.

## Reproduction and controls

Run `python check.py --bench` from this directory (or via its path). Run `python -O check.py --receipt optimized-receipt.json` for assertion-independent controls. The checker uses explicit exceptions, not assert. It compares all four core/no-fixed combinations for every bound 0..26 with the frozen incoming kernel, including identical stage state/transition counts. Independently unmerged rational odd partitions and direct core substitution control all modes through 20. A separate permutation action on unordered labelled edges computes Burnside counts through 7 without the odd-cycle formula. A known A000568 prefix through 12 and invalid-bound controls are included.

The receipt includes interpreter/platform, SHA-256 source hashes, actual samples, coefficient hashes, mathematical transition counts and precomputed transfers. Benchmarks alternate execution order over three repetitions, compare exact complete coefficient vectors after each pair, and include setup/precomputation/output Fraction conversion. They are timings on this shared host, not claims about other implementations, languages, orders, machines or universal fastest performance. Tiny bounds may lose to overhead. Memory is not benchmarked here.

## Remaining boundary

This removes state-dependent rational division and shrinks arithmetic operands; it does not remove future-response distinctions. A structural reduction of the response graph remains open. Prime-order/channel-sharing improvements could be combined later but are not implemented or credited in this packet. For pair-free cores, the number of signed excess contributions can still dominate runtime.

## Recorded finite timing results

- Bound 80, ordinary: median speedup 1.880x; 36,384 contributions; 6,753 transfer entries.
- Bound 80, core: median speedup 1.416x; 645,752 contributions; 93,858 transfer entries.
- Bound 100, ordinary: median speedup 2.016x; 90,219 contributions; 11,109 transfer entries.
- Bound 100, core: median speedup 1.394x; 1,697,975 contributions; 186,770 transfer entries.
- Bound 150, ordinary: median speedup 2.041x; 601,461 contributions; 27,274 transfer entries.
- Bound 200, ordinary: median speedup 2.269x; 2,994,814 contributions; 51,374 transfer entries.

These are the three-sample alternating-order measurements in receipt.json, including concurrent-host variability.

Timing attribution caveat: the actual incoming API also collects response-coordinate-update totals, tuple-cell counts and packed-response-bit counts. The new API retains row states/transitions but omits those extra telemetry passes. The end-to-end speedups above therefore combine the arithmetic change with telemetry reduction; they are not an isolated measurement of the arithmetic identity. A separately hashed telemetry-matched baseline ablation is included with the independent review when complete.

## Polynomial transfer-family bound

Even if every degree is visited at every stage, the number of ordinary transfer entries is at most

  sum_(odd ell≤N) sum_(0≤d≤N) [1+floor((N−d)/ell)] = O(N² log(N+2)).

For core mode put m=floor((N−d)/ell). There is one k=0 entry and exactly m(m+1)/2 possible pairs k≥1,e≥0,k+e≤m, giving

  sum_(odd ell≤N) sum_(0≤d≤N) [1+m(m+1)/2] = O(N³).

Only encountered degrees are materialized. These are bounds on the number of exact transfer divisions, not bit complexity or the full algorithm: response-state contributions remain unchanged and integer operands grow. They explain why factoring the degree sidecar out of response states is a real arithmetic reuse mechanism rather than a new state quotient.

The independent review's telemetry-stripped incoming ablation also retains a speed advantage in its own alternating three-sample runs: ordinary bounds80/100/150/200 give1.833x/1.971x/2.158x/1.855x; core bound100 gives1.326x. See the separate `session22-enumeration-review` packet's `telemetry.patch`, hashed `baseline_light.py`, `benchmark.json` and `benchmark_extra.json`. These samples reduce the identified telemetry confound but do not attribute every implementation difference to the identity. Ratios are rounded; receipts retain full samples.
