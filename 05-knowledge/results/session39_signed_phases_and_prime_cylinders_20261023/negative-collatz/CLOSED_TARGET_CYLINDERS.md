# Closing any finite signed valuation prefix to a chosen target

Session 39 extension, 23 October 2026. **PROVED elementary construction; FINITE-EXACT replay controls.** This extension does not alter the frozen signed-prime packet. It makes no novelty claim for inverse-orbit density and does not prove Collatz convergence.

## 1. Exact closing-step theorem

Use Uε(n)=(3n+ε)/2^v2(3n+ε), ε∈{−1,+1}, on positive odd integers. Given an arbitrary finite word w=(a0,...,a(r−1)) of positive valuations, let

A=sum ai, B=sum_(j=0)^(r−1) 3^(r−1−j)2^(sum_(i<j) ai).

For the empty word set r=A=B=0. Let y be any positive odd integer with 3∤y. There is precisely one exponent class b modulo L=2·3^r satisfying

2^(A+b)y ≡ ε(2^A+3B) (mod 3^(r+1)).                       (1)

For every sufficiently large positive b in that class,

n_b=[2^(A+b)y−ε(2^A+3B)]/3^(r+1)                         (2)

is a positive odd integer whose first r+1 exact valuations are (a0,...,a(r−1),b), and Uε^(r+1)(n_b)=y. In particular, every finite valuation cylinder contains infinitely many starts that reach the chosen y in one extra odd step after the prescribed prefix.

### Proof

The residue ε(2^A+3B)y^(−1) is a unit modulo 3^(r+1). The number 2 generates that unit group and has order 2·3^r, hence (1) has exactly one exponent phase. For completeness, the order assertion follows from v3(4^(3^r)−1)=r+1: the base case is 4−1=3, and binomial expansion of (1+3^(r+1)u)^3 increases the valuation by exactly one when 3∤u. The order of 4 modulo 3^(r+1) is 3^r, while the order of 2 is even because 2≡−1 modulo 3. Thus it is 2·3^r, the size of the unit group.

Equation (1) makes (2) integral. The numerator is odd: if r=0 its subtracted constant is ε; if r≥1, B is odd and 2^A+3B is odd. The denominator is odd, so n_b is odd. As b increases through its phase, its numerator tends strictly to infinity, proving eventual positivity and distinctness.

The carry of the extended word (w,b) is 3B+2^A. Formula (2) says exactly that its formal affine endpoint is y. Since n_b is integral and y is odd, the finite-cylinder sufficiency proof applies. Explicitly, for each prefix of the extended word its numerator 3^j n_b+εB_j is divisible by 2^A_j: reduce the full endpoint equation modulo that power, using the carry decomposition into prefix plus 2^A_j times suffix. All intermediate values are thus integers. Working backward from the odd endpoint in 3n_j+ε=2^a_j n_(j+1), each previous n_j is odd. Each valuation is exact, and forward positivity follows from n_b>0. This authenticates the actual orbit, not merely rational branch composition. ∎

The target restriction is necessary for a nonempty forward path: every odd-map image is nonzero modulo 3. No universal target assertion for y divisible by 3 is made. The theorem includes the empty prefix: then 2^b y≡ε modulo 3, giving the familiar one-step inverse family.

## 2. Three basins inside every finite observation

For ε=−1, take y=1,5,17. These lie on the three different directly verified positive 3x−1 cycles. The theorem supplies infinitely many starts in each basin inside every finite valuation cylinder. This result holds without assuming that these are the only basins or that every input eventually enters a cycle.

An explicit common prefix is (1,2), corresponding to the class 5 modulo 16:

- 38837 has exact word (1,2,17) and lands at 1.
- 6069 has exact word (1,2,12) and lands at 5.
- 21 has exact word (1,2,2) and lands at 17.

The last-exponent phase period is 18 in all three examples. It is the exponent phase, not an assertion that the start values form an ordinary arithmetic progression. The start families are exponential in b.

**Finite-observer consequence.** No procedure using only the valuation stream, and always stopping after finitely many entries, can correctly name the eventual one of these three cycles for all positive odd inputs, even under the promise that the input reaches one of these cycles. On any purported terminal observation, the same finite prefix occurs for a start in each of the other basins. The procedure receives indistinguishable information and cannot give three different correct answers. This is an observational limitation, not algorithmic undecidability with access to the input value: exact orbit values immediately supply additional information, and a checked finite path to a cycle is a valid certificate.

This strengthens the distinction between a local modular observer and a terminal certificate. It does not identify different basins, merge their cycles, or show that basin membership has positive natural density in each cylinder.

## 3. Inheritance and primary-source comparison

Inverse-orbit ubiquity in congruence classes is established literature, not a newly discovered phenomenon here. Roger E. Zarnowski, “The Congruence Structure of the 3x+1 Map,” *Fibonacci Quarterly* 46/47 (2008/2009), 115–125, Theorem 3, shows that classes modulo 2^a3^b generate all targets not divisible by 3 infinitely often under the shortcut map. Primary paper: https://www.fq.math.ca/Papers1/46_47-2/Zarnowski_1-09.pdf . Our self-contained formulation retains a prescribed accelerated valuation prefix, the sign, and the exact one-extra-odd-step exponent phase. No historical priority is asserted for this formulation either.

Closest inspected repository context at snapshot e5d18e55:

- `05-knowledge/results/collatz_basin_observer_20261009.md` distinguishes formal affine state, evaluated equality, terminal cycle and clock phase; proves finite-observer limitations in a different pair-state setting; and carefully separates Haar-null target basins from ordinary integer questions. The present three-basin cylinder result concerns actual positive 3x−1 starts and a valuation-only observer, not its base-p pair chain.
- `05-knowledge/results/collatz_terras_inverse_clock_bridge_20261007.md` treats Haar count-clock/inverse-clock recurrence, explicitly without concluding universal ordinary convergence. No probabilistic assertion from it is required for the present elementary construction.
- `05-knowledge/results/collatz_minus17_return_20261003.md` authenticates repeated valuation cylinders and the seven-step negative 17 word. Its cycle certificate is inherited.

The searches recover relevant predecessors but do not establish an exhaustive priority census.

## 4. Non-implications and remaining obligations

Do NOT combine “each cylinder contains infinitely many primes” from the signed-prime packet with “each cylinder meets each basin” to conclude that each basin-cylinder intersection contains infinitely many primes. Those are separate existence statements. The constructed target families (2) are exponential sequences; their primality is not supplied by Dirichlet's theorem for arithmetic progressions.

The exponent clock L=2·3^r grows with prefix length. This construction does not give a uniformly bounded closing cost, a useful uniform density lower bound, first-exit payment for arbitrary given starts, or the annulus 153N estimate from other work. It chooses new starts; it does not prove that an already specified start has a closing step. Global convergence and complete cycle classification remain open.

## 5. Executable controls

`verify_closed_targets.py` has no repository imports and uses exceptions rather than removable assertions. Normal and optimized Python have identical receipts:

- 1,210 unique exponent-phase checks: all 121 words of lengths 0 through 4 over {1,2,3}, both signs, and targets 1,5,7,17,25.
- 3,630 exact trajectory replays: three successive admissible exponents for every case.
- 2,000 checks that actual images of the first 1,000 odd starts, both signs, avoid multiples of 3.
- The three distinct basin examples above are retained in JSON.

Reproduce with `python3 verify_closed_targets.py` and `python3 -O verify_closed_targets.py`, comparing output with `closed-targets-receipt.json` and `closed-targets-receipt.optimized.json`. These finite checks validate the implementation; the all-word result is the proof in Section 1.
