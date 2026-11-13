# Session 295: Compare the constructions by what they let you do next

**Chan–Pak remains open.** Thirteen construction families are compared
by actual inputs, retained state, legal continuation and the first
boundary where the comparison breaks. The new results make some
of those boundaries exact.

- [Construction atlas and selective-port filtration](atlas.md).
  With at most q nonuniform exterior vertices, the k-run sector is
  divisible by max(k−1−q,0)!·max(k−q,0)!. Glue each exceptional
  excursion to its adjacent core runs before permuting anything.
  At modulus 27 and q=1 the cutoff k≥7 is sharp. A single selective
  port changes the previously 27-obstructed core to H=1147, prime to 3.
- [Exact observer envelope for uniform contexts](ARITHMETIC.md).
  Weighted evaluations of the odd-cycle polynomial generate precisely
  the same integer ideal as weighted cover coordinates. Under ordinal
  sums the polynomials multiply. The p-adic valuation is a finite,
  effectively computable minimum of affine costs for a fixed list of
  cores. Repeated copies eventually have a correction given by the
  factorial tariff at the first minimizing evaluation. This makes the
  earlier 8/9/27 calculations instances of a general mechanism.
- [Continuation equivalence and distant precedents](bridges.md).
  Uniform substitution has an exact finite continuation algebra
  modulo any fixed modulus; modulo 3 it has five states, all realized.
  Swapping the order of a directed triangle and a singleton is invisible
  to *every* uniform context, but the same selective attachment
  distinguishes the two nonisomorphic cores by counts 9 and 5.
  Ore-grade correction, Rule 30, lonely-runner observability, the
  planar Jacobian and square-sum work supply distinct controls on
  what an observer is allowed to forget.

**Design implication.** Separate absolute target, physical port guards,
and continuation response. They are not interchangeable. A uniform
observer is the right comparison for uniform substitution, but a
selective repair must carry more information. The selective theorem
still removes all sufficiently high run sectors at fixed arithmetic
precision; it does not guarantee that the remaining low-run sectors
can be tuned to a prescribed count or renew the interface. That is
the next constructive obligation.

## Reproduce and review

Run `python check_arithmetic.py`, `python verify_selective.py` and `python verify_bridges.py`
from this directory. Compare with [arithmetic.out](arithmetic.out)
[selective.out](selective.out), and [bridges.out](bridges.out). Normal and optimized runs agree.
The arithmetic transform and finite Pareto recursion were independently
reviewed. The selective filtration and sharp model were independently
derived using compound blocks; the sharp sector is also checked by
subset DP. The controls test identities, not an extrapolation to
coverage. [Initial plan](PLAN.md).
