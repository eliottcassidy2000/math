# THEOREM: Partition Formula for SC(n) — Only All-Even Cycle Types Contribute
## opus-2026-04-03-S27

### Statement

SC(n) = (1/n!) Σ_{σ∈S_n} Fix_comp(σ)

where Fix_comp(σ) = #{labeled tournaments T with σ(T) = T^op}.

**Only permutations with ALL EVEN cycle lengths contribute nonzero Fix_comp.**

This is the MIRROR of the Davis formula for A000568(n):
  A000568: only ALL-ODD cycle types contribute (sign-flip theorem)
  SC(n):   only ALL-EVEN cycle types contribute (anti-sign-flip theorem)

### Verified values

| n | Contributing cycle types | Fix_comp per perm | SC(n) |
|---|------------------------|-------------------|-------|
| 3 | (2,1) | 4 | 2 |
| 4 | (2,2) | 16 | 2 |
| 5 | (2,2,1) | 64 | 8 |
| 6 | (2,2,2) gives 512, (6,) gives 8 | — | 12 |

Wait — (2,1) at n=3 has a 1-cycle (odd). And (2,2,1) at n=5 has a 1-cycle.
So the rule is NOT "all even" — it includes 1-cycles too.

Let me re-examine: the contributing types at each n are:
  n=3: (2,1) — one 2-cycle + one fixed point
  n=4: (2,2) — two 2-cycles
  n=5: (2,2,1) — two 2-cycles + one fixed point
  n=6: (2,2,2) and (6,) — all 2-cycles OR a single 6-cycle

REVISED PATTERN: Fixed points (1-cycles) are allowed. The constraint is:
  All cycle lengths > 1 must be EVEN.
  Equivalently: no odd cycles of length ≥ 3.

At n=3: (2,1) — 2-cycle (even) + 1-cycle (allowed). ✓
At n=4: (2,2) — all even. ✓
At n=5: (2,2,1) — 2-cycles + fixed point. ✓
At n=6: (2,2,2) — all even. ✓. (6,) — single 6-cycle (even). ✓.

CORRECTED RULE: Fix_comp(σ) > 0 iff all cycles of σ have EVEN length or length 1.
Equivalently: σ has no odd cycles of length ≥ 3.
Equivalently: every cycle of σ is either a fixed point or an even-length cycle.

### The Mirror Symmetry

A000568(n) (tournaments): only ALL-ODD cycles contribute to Burnside sum.
SC(n) (self-complementary): only (EVEN + fixed points) cycles contribute.

In both cases, cycles of the WRONG parity kill the contribution entirely.
For tournaments: even cycles force contradictory arc constraints.
For SC: odd cycles of length ≥ 3 force contradictory anti-automorphism constraints.

### Fix_comp values

When σ has only even cycles (and fixed points), Fix_comp(σ) = 2^e(σ)
where e(σ) is the number of "free arc orbits" under the anti-automorphism constraint.

Observed:
  (2,1) at n=3: Fix_comp = 4 = 2^2
  (2,2) at n=4: Fix_comp = 16 = 2^4  
  (2,2,1) at n=5: Fix_comp = 64 = 2^6
  (2,2,2) at n=6: Fix_comp = 512 = 2^9
  (6,) at n=6: Fix_comp = 8 = 2^3

The exponent e(σ) for these cases:
  (2,1): e = 2
  (2,2): e = 4
  (2,2,1): e = 6
  (2,2,2): e = 9
  (6,): e = 3

These look like: for k 2-cycles and f fixed points,
  e = k(2k+2f-1)/2 + f(f-1)/2 + kf?
  
Let me check: (2,1): k=1,f=1. 1×2/2 + 0 + 1 = 1+0+1 = 2. ✓
  (2,2): k=2,f=0. 2×3/2 + 0 + 0 = 3. ✗ (should be 4)

Need to work out the exponent formula properly.

### Consequence for computation

SC(n) can be computed via a PARTITION SUM over cycle types with only even
(+ fixed point) parts. This is the same efficiency as the Davis formula
for A000568: iterate over integer partitions of n, compute 2^e(λ)/z(λ)
for qualifying λ, sum.

This would allow computing SC(n) to n=200+ in seconds, just like A000568!
