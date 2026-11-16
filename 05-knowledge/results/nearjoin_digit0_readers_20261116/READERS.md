# Can the chosen reader break the unwanted divisor without overshooting the absolute target?

Session klein-2026-11-16-nearjoin. Script `scripts/probes.py`; receipt `probes.out`. The question is the two-sided
design rule of `metagraph_flow_observers_20261116` §5 ("first test whether the chosen core can leave the unwanted
residue through the permitted reader, then test whether activation already overshoots the absolute budget"), put to the
constructions of this programme.

## 1. Homogeneous readers are module insertions (PROVED, THM-5260 (A))

A homogeneous reader attaches exterior vertices each of which beats all of the core `T` or loses to all of it: the
core is a module, the reader is the gadget `Q` with the core substituted for `τ`. Every directed cycle of the result
meets the core in one run, so the count is a nonnegative combination of the core's run sums — the ordered path-cover
numbers `Q_j(T)` of the flow-observers packet — and a reader whose `τ` lies on no cycle of `Q` keeps every divisor of
`H(T)` (the count is `c·H(T)`). So "breaking the unwanted divisor" is only possible through readers whose module vertex
lies on a cycle, and only if the divisor does not divide the relevant run sums (the packet's readout gcd).

## 2. The cheapest breaking readers and their cost (FINITE-EXACT; `probes.out`)

The cheapest readers of that kind are a single cycle through the module: `l` exterior vertices `a_1 → … → a_l` in a
transitive chain, the core beating `a_1, …, a_(l−1)` and `a_l` beating the core. Every mixed cycle is a run `P` of the
core followed by an increasing segment of the chain that ends at `a_l`, so

    R_l(T) = H(T) + 2 Σ_{m=2}^{l} C(l−1, m−1) · Π_{|P| ≡ m+1 (mod 2)}(T),   Π_par(T) = Σ_{runs P of parity par} H(T − P),

checked against direct subset-DP counts on all strong cores of orders 3–7 and `l ≤ 5` (1,231 cases, no failure).
For the `+2` question the absolute target is the five-window `{3H, …, 3H + 4}` around the forced midpoint `3H`:

| `l` | least `R_l − 3H` over strong cores of orders 3–7 | median | probes landing on `3H + 2` or `3H + 4` |
|---|---|---|---|
| 2 | 6 | 1,038 | none |
| 3 | 24 | 3,276 | none |
| 4 | 60 | 7,752 | none |
| 5 | 132 | 5,356 | none |

No single-cycle homogeneous reader lands on `3H + 2` or `3H + 4` for any strong core of orders 3–7; the least overshoot
over the midpoint is 6, attained by the two-vertex probe (`H + 2Π_odd`) on the smallest cores. The readout gcd of
`{H, R_2, …, R_5}` equals `H` for 9 of the 396 cores (no divisor of `H` is ever broken by these probes) and exceeds 1
for 166 of them.

## 3. The answer

- **Can the chosen reader break the unwanted divisor?** A homogeneous reader can, exactly when the divisor does not
  divide the readout gcd of the core; the flow-observers packet's H49 witnesses show both outcomes.
- **Without overshooting the absolute target?** Not in the `+2` window: on every strong core tested, breaking the
  divisor 3 of the midpoint costs at least `+6`, and typically a thousandfold more, because the cheapest breaking reader
  adds `2Π_odd(T) ≥ 2 Σ_v H(T − v)`. The exact `+2` is reached only by an *inhomogeneous* reader — the height-two tower
  at a digit-1 vertex, which reads one vertex of the core — and only when THM-5260 (C)'s criterion `2^k B_t = h + 1`
  holds. The port, not the exterior, carries the digit; this is the same conclusion as THM-5260 (A) from the other side.
- The near-join reader at a digit-0 vertex (`B = h`) never breaks the divisor `h`: its outputs are `h(2^k + 1)`.
