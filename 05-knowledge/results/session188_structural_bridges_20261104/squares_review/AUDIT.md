# Independent review: square-sum edge and parity inventory

**Verdict: PASS. No mathematical issue found.**

I independently checked the proof and replayed the producer checker in normal and optimized Python. Its outputs are identical and match the supplied receipt.

## Proof review

- At the birth of vertex n, a neighbor i must satisfy 1<=i<n and n+i=k^2. Thus n<k^2<=2n-1, giving exactly isqrt(2n-1)-isqrt(n). Both the lower cutoff and exclusion of the loop n+n are correct, including n=1.
- Same-parity edges have even-square sum 4k^2 and the unique center-radius parametrization {2k^2-t,2k^2+t}. The available radius is min(2k^2-1,n-2k^2). For odd n both bounds are odd; positivity gives 2k^2<=n-1. For even n the unsaturated radius is even, while saturation occurs at 4k^2<=n+1 and gives an odd radius. These yield the two stated floor-square-root formulas.
- For a Hamilton path, summing degrees over the two parity classes gives 2(p_oo-p_ee)=2(V_o-V_e)-(e_o-e_e). The endpoint count is two incidences for n>=2, so the displayed allowed intervals are correct.
- For a cycle every vertex has degree two and there is no endpoint term. Subtracting used from available edges gives the stated lower bounds on omitted odd-odd edges. These are necessary inventory facts, not sufficient conditions for Hamiltonicity.

A more detailed retained identity is

    omitted_oo - omitted_ee = (O_n-E_n) - (n mod 2) + (e_o-e_e)/2

for paths. This specializes to the packet's universal lower bounds after maximizing the used parity difference. No correction is needed.

## Independent checks

`independent.py` does not import producer functions. It adds edges by their later endpoint, directly generating square values, and verifies both births and cumulative parity inventory for every n through 10,000.

It separately enumerates all parity arrangements with the required odd/even vertex counts through n=18: 99,293 path patterns and 99,291 cycle patterns. These check the endpoint and degree identities without relying on the square-edge formula.

It also verifies actual square-sum Hamilton paths on 15, 16, and 17 vertices and checks both the endpoint-specific omitted-edge identity and the weaker universal omission bound. These witnesses are finite controls, not an existence theorem at other orders.

Normal and optimized independent outputs are identical. All checks use explicit exception-raising conditions.

Artifacts: `normal.json`, `optimized.json`, `independent.py`, `independent.json`, and `independent.optimized.json`.

## Added order-25 application: separate audit

The subsequently added `path25.py` and order-25 paragraph also PASS. The producer's normal and optimized outputs agree byte for byte with `path25.json`.

A separate memoized endpoint-rooted computation in `path25_independent.py` builds neighbors from square values rather than pairwise square tests. It verifies that 18 is the unique leaf, so every Hamilton path has endpoint 18. It counts exactly ten paths starting there; reversal gives all twenty directed paths, with no duplicate orientation or missing start case.

The independent histogram and every ordered endpoint count agree exactly with the supplied receipt. The available same-parity inventory is O_25=9, E_25=6. Every observed path has p_oo-p_ee equal to 1 or 2, and the omitted odd-odd count is between 1 and 4, as claimed. Normal and optimized independent receipts agree. The memoized computation uses 2495 distinct states.

Additional artifacts: `path25.normal.json`, `path25.optimized.json`, `path25_independent.py`, `path25.independent.json`, and `path25.independent.optimized.json`.
