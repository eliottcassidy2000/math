# Irrational Collatz sources selected by positive Mahler integers: contraction, logarithmic memory, and a complexity dichotomy

2026-10-14, mac-mini-2026-10-14-wilson.

- **PROVED:** [THM-4654](../../01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md).
- **FINITE-EXACT:** controls on 4,980,199 pairs and on an equality family.
- **OPEN:** Mahler's Z-number problem and Collatz.

This answers the target "stronger restrictions on the irrational Collatz sources selected by positive Mahler integers". That target is the next obligation of [mahler_carry_interface_20261014i.md](mahler_carry_interface_20261014i.md) §6.

## 1. What was known

The interface note built the code-preserving isometry `Ψ` with `ΨM = TΨ`. It proved that a positive Mahler integer `A` with an all-depth no-11 code selects an **irrational** Collatz source `z_A = Ψ(A)`. It gave an explicit separation from rationals of height `Q`: `v₂(z_A − u/v) < 14Q + 2 + bitlength(A + 1)`.

Two explicit series make the duality concrete. Let the ones of the carry code sit at `n₁ < n₂ < …`. Then
- `A = −(1/3) Σ_m (2/3)^(n_m)` (2-adically);
- the Z-number is `ξ = A + (1/3) Σ_m (2/3)^(n_m)` (really).

So a Z-number is the difference between the real and the 2-adic values of **one** series. Meanwhile the Collatz source is `z_A = −Σ_m 2^(n_m) 3^(−m)`.

## 2. What is new (THM-4654)

1. **The safe language is a one-dimensional automaton,** `m ↦ min((3/2)m − e, 1)`.
   - This gives exact counts: `N(ℓ) ~ 1.551 (3/2)^ℓ`.
   - It gives the block law: no `11`; a gap 2 forces a gap of at least 6; `(2, 6)` forces at least 3. This is Parry's criterion against the greedy expansion of 1, `d = 1010000010010010100000000010…`.
   - A safe word of length `ℓ` has at most `⌊ℓ/3⌋ + 1` ones.
2. **Real safety contracts the Collatz maps in the real metric.** Along any safe word the slope is at most `3(3^(1/3)/2)^ℓ` and the offset is below `11/10`.
   - So the orbit of every rational shadow enters the window `[−1/2, 8/5)`.
   - Two shadows with the same denominator `v` that read the same `≈ 2.12 log₂ v` code bits then **coincide**: logarithmic memory.
3. **Complexity–approximation inequality.**
   `v₂(z_A − u/v) < log₂3 · (n₁ + min(2.1v + 1, ℓ(v) + N_e(ℓ(v)))) + log₂(A + 1)`.
   - Unconditionally this improves the interface bound: linear in `v`, logarithmic in `|u|`.
   - If the carry code has factor complexity `O(ℓ)`, then `z_A` has a **finite 2-adic irrationality exponent**, `μ <= 3.36(K + 2)`.
   - The periodic-run ingredient (the interface's eq. (8)) is sharp to 0.507 bits on the family `A_t = (2^(18t+8) − 9)/19`, whose codes shadow the Collatz cycle 1/5 for exactly `18t + 8` steps.

## 3. Reading: one duality per scale

- The Mahler condition is **real** (tails below 1). The positive-integer condition is **2-adic** (the series equals `−3A`). The Collatz source is **2-adic**.
- Real safety enters Collatz only through contraction. That is enough to turn the 2-adic approximation problem into a finite-state one with logarithmic memory.
- The new restriction (THM-4654 (d)) is a **dichotomy**:
  - either a hypothetical Z-number's carry code has super-polynomial factor complexity;
  - or the irrational Collatz source it selects is 2-adically badly approximable to polynomial order.

  Low-complexity Z-numbers are the natural "structured" candidates (Sturmian, morphic). They are therefore exactly the ones whose Collatz sources would have finite irrationality exponents. A transcendence-type input for such `z_A` would exclude them.
- The same theme appears in the OpenAI Chowla paper's deletion of minimal prohibited words before a trace estimate; see [chowla_ordinary_exponent_ledger_20261014.md](chowla_ordinary_exponent_ledger_20261014.md) §3. In both settings, a suffix condition defines forbidden factors, and deleting them buys structure on the surviving paths.

## 4. Reproduction

```text
python3 04-computation/experiments/mahler_collatz_sources_20261014/safe_language.py
python3 04-computation/experiments/mahler_collatz_sources_20261014/collatz_sources_check.py 200000
```

The first script uses exact `Fraction` arithmetic throughout; the second took under 4 minutes on the mac-mini. Outputs are saved alongside the scripts.

## 5. Next obligations

1. **Factor complexity of Mahler codes of positive integers.** The periodic-run inequality gives only aperiodicity (`N_e(ℓ) >= ℓ + 1`). A super-linear lower bound for every positive `A` with long safe prefixes would remove the linear-complexity branch of the dichotomy.
2. **2-adic transcendence input.** Is there a Mahler-method (functional-equation) proof that `z_A` is not a quadratic irrational when the code is morphic? The series `−Σ 2^(n_m)3^(−m)` with morphic `n_m` is of Mahler-function type.
3. **The marked sheet.** Combine with the triangular quotient of `recursion_carry_bridges_20261014i.md`. There the complementary sheet `1 − e` turns safety (`Y < 1`) into `Y > 1`. Does it select the negatives `−1 − z_A`?
