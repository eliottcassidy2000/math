# Irrational Collatz sources selected by positive Mahler integers: bounded entry, 2-adic memory, and a complexity-approximation inequality

2026-10-14, mac-mini-2026-10-14-wilson. Audited by audit D (corrections recorded as MISTAKE-612).

- **PROVED:** [THM-4654](../../01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md).
- **FINITE-EXACT:**
  - controls on 4,980,199 pairs (session);
  - 2.1 million worst-case cases (audit D);
  - an equality family.
- **OPEN:** Mahler's Z-number problem and Collatz.

This addresses the target "stronger restrictions on the irrational Collatz sources selected by positive Mahler integers", the next obligation of [mahler_carry_interface_20261014i.md](mahler_carry_interface_20261014i.md) §6. The honest yield is a sharper, checkable separation of `z_A` from rationals, plus a clean account of where it comes from. It is not a restriction that by itself bears on the existence of Z-numbers.

## 1. What was known

- The interface note built the code-preserving isometry `Ψ` with `ΨM = TΨ`. It proved that a positive Mahler integer `A` with an all-depth no-11 code selects an **irrational** Collatz source `z_A = Ψ(A)`. It gave an explicit separation from rationals of height `Q`: `v₂(z_A − u/v) < 14Q + 2 + bitlength(A + 1)`.
- Flatto (Contemp. Math. 135, 1992) frames Mahler Z-numbers through the `β = 3/2` shift.
- Dubickas (Glasgow Math. J. 51, 2009) proves that the Mahler code of every positive integer has factor complexity `liminf N(ℓ)/ℓ >= log 2/log(3/2) ≈ 1.71`.
- Andrieu–Eliahou–Vivion (arXiv:2510.11723) identify Mahler's binary words with minimal words in base 3/2.

Two explicit series make the real/2-adic duality concrete. Let the ones of the carry code sit at `p₁ < p₂ < …`. Then
- `A = −(1/3) Σ_m (2/3)^(p_m)` (2-adically);
- the Z-number is `ξ = A + (1/3) Σ_m (2/3)^(p_m)` (really).

So a Z-number is the difference between the real and the 2-adic values of **one** series. The Collatz source is `z_A = −Σ_m 2^(p_m) 3^(−m)`.

## 2. What THM-4654 adds

1. **The safe language as a one-dimensional slack automaton,** `m ↦ min((3/2)m − e, 1)`. This is THM-4072's follower graph in slack coordinates.
   - Exact counts: `N(ℓ) ~ 1.55104519·(3/2)^ℓ`, with constant `3/Σ_j j d_j (2/3)^j` (THM-3848's renewal law).
   - The block law: no `11`; a gap 2 forces a gap of at least 6; `(2, 6)` forces at least 3. This is Parry's criterion.
   - The exact density bound: at most `⌊ℓ/3⌋ + 1` ones.
2. **Real safety contracts the Collatz maps.** Along any safe word the slope is at most `3(3^(1/3)/2)^ℓ` and the offset is below `11/10`. So every rational shadow enters the bounded window `[−1/2, 8/5)`.
3. **The memory is 2-adic, not real.** Inside the window, a shadow's state is determined by its next `⌊log₂(2.1v)⌋ + 1` code bits, because Collatz parity vectors determine residues (Terras, Everett).
   - The pre-audit text attributed this memory to real contraction, with a longer window of about `2.12 log₂ v` bits. That lemma was true but vacuous (audit D, M2).
4. **Complexity–approximation inequality.** With `N` the number of distinct future `ℓ`-factors (`ℓ = ⌊log₂(2.1v)⌋ + 1`):

   `v₂(z_A − u/v) < min( log₂3·(n_in + ⌊2.1v⌋ + 1) + log₂(A+1),  max(n_in + N + ℓ, log₂3·(n_in + N) + log₂(A+1)) )`.

   - Unconditionally the bound is linear in `v` and logarithmic in `|u|`, improving the interface's `14Q`.
   - Under linear factor complexity (automatic or primitive morphic codes; Sturmian codes are excluded by Dubickas), `z_A` has a finite 2-adic irrationality exponent, `μ <= log₂3·(K + 1) + 3.361`.
   - The periodic-run ingredient (the interface's eq. (8)) has an optimal additive constant.

## 3. Reading

- **Three metrics, three roles.**
  - The Mahler condition is **real** (tails below 1). Through contraction it gives boundedness.
  - The Collatz source is **2-adic**. The parity bijection gives finite memory.
  - The positive-integer premise (the series equals `−3A`) forbids periodic shadows.
- **Dichotomy, stated correctly.** For a hypothetical Z-number, either `limsup N_e(ℓ)/ℓ = ∞`, or `z_A` has a finite 2-adic irrationality exponent.
  - This is not a restriction by itself: a finite exponent is the generic, expected behaviour, proved for real automatic numbers by Adamczewski–Cassaigne.
  - Excluding low-complexity Z-number codes would need a Liouville-type (infinite-exponent) input for `z_A`, or an argument on the Mahler side.
  - The pre-audit version said that low-complexity codes would be "exactly" the ones with finite exponents, and that a transcendence input would exclude them. That ran backwards (audit D, M1).
- **ANALOGY only.** The OpenAI Chowla paper deletes sites where minimal prohibited words start; see [chowla_ordinary_exponent_ledger_20261014.md](chowla_ordinary_exponent_ledger_20261014.md) §3. Its words are defined by a suffix **divisibility** condition, and deleting them costs density. THM-4654 deletes nothing: its constraint is the `β`-shift **inequality**, which Z-number codes satisfy by hypothesis.

## 4. Reproduction

```text
python3 04-computation/experiments/mahler_collatz_sources_20261014/safe_language.py
python3 04-computation/experiments/mahler_collatz_sources_20261014/collatz_sources_check.py 200000
```

The first script uses exact `Fraction` arithmetic throughout; the second takes a few minutes on the mac-mini. Audit D's independent scripts are in `04-computation/experiments/audit_D_20261014/`.

## 5. Next obligations

1. **Super-linear factor complexity.** Dubickas gives `liminf N(ℓ)/ℓ >= 1.71` for Mahler codes of positive integers. A super-linear lower bound, at least along safe prefixes, would remove the linear-complexity branch.
2. **A Mahler-side argument.** For automatic or morphic `p_m`, the 2-adic series `−(1/3)Σ(2/3)^(p_m)` is of Mahler-function type. Can it be a positive integer while the real series stays safe?
3. **The marked sheet.** Combine with the triangular quotient of `recursion_carry_bridges_20261014i.md`, where the complementary sheet `1 − e` turns safety (`Y < 1`) into `Y > 1`. Does it select `−1 − z_A`?
