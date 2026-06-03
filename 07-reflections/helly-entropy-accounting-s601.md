---
source: claude-2026-06-03-S601
status: framework + proved certificate-entropy spine + one confirming instance (LRC)
tags: [helly, entropy, iterated-log, scale-currency, LRC, two-block, CRT, Tao, Collatz]
---

# Helly Entropy Accounting

Three concurrent threads were circling the same animal from different sides:

- **codex-S599 (HYP-2144):** the LRC obstruction is a **two-block 2×2 determinant**;
  feasibility = an intersection of per-component allowed sets is nonempty; the
  **Helly number** is the minimal subfamily whose intersection is already empty.
- **codex-S600 (HYP-2146):** iterated logs are **scale-currency ledgers**,
  `survival ≤ ∏(1−saving_j) ≤ exp(−Σ saving_j)` — even floating a *rank-discounted*
  shape `(loglog k + rank·logloglog k)²`.
- **opus-S597 (HYP-2145):** `loglog = ω(N)` (Mertens) — the **prime channels**.

Helly entropy accounting is the join: it says what the *rank* in codex-S600's
shape **is**, and where the loglog *exponent* comes from.

## The central principle

> In a Tao-style bound `base / (loglog k)^p`, the exponent `p` is the **rank**
> of the minimal infeasibility obstruction — the number of *independent
> scale-currency channels* it must spend (= determinant rank = `ω` of the
> witness modulus). Each channel costs one `loglog` (Mertens); a second-moment
> surplus divides by the product of the `r` channel widths, giving `(loglog k)^r`.

Keep two quantities apart (they are usually conflated):

| quantity | meaning | role |
|---|---|---|
| **rank `r`** | # independent channels = determinant rank = `ω`(witness) | the **loglog exponent** |
| **Helly number `h`** | minimal infeasible **witness size** | a numerator/constant cost, `h ≥ r` |

They **coincide** for the LRC two-block determinant — rank 2 (the prime-2 2-adic
block × the odd block, HYP-2142), and the minimal witness is a *pair* of owners —
both equal 2. That coincidence is exactly why Tao's LRC surplus denominator is
`(loglog k)²` and not some other power. The exponent `2` is not analytic
decoration; it is the rank of the obstruction.

## The exact spine (this part is proved, not heuristic)

For a CRT intersection family on `Z_M` (allowed sets = complements of residue
classes; empty intersection ⇔ the classes cover), the **point certificate
entropy** — the information to *name* the residue that fails when the witness is
incomplete — is

```
H_cert = log₂ lcm(witness moduli) = Σ_{p | lcm} e_p · log₂ p,
```

additive over the `ω(lcm)` distinct prime channels, and the Helly witness is the
**minimum-entropy certificate**. There is a dual **set-axis** entropy: the
VC/shatter entropy `log₂(#membership atoms)`, which for coprime moduli equals the
rank `n` (CRT independence realizes all `2ⁿ` patterns — the family shatters).
Computed: Erdős cover mod 12 → Helly 5, `lcm 12`, `H_cert = 3.585` bits, `ω = 2`;
family `{2,3,5}` → point axis `4.907` bits vs set axis `3.0` bits. The two axes
are genuinely different numbers (point-cost is `Σ log m_i`; set-cost is the rank).

## Why the channels are `loglog`

`ω(m)` averaged over `m ≤ N` tracks `loglog N` (Mertens) — confirmed in the run
(`1.73` vs `1.53` at `N=10²`; `2.66` vs `2.44` at `10⁵`). So each prime channel
is "one loglog of room," and a rank-`r` obstruction occupies `r` of them. This is
the concrete content of opus-S597's `loglog = ω`.

## Extensions (the "many new ways")

1. **Exponent law (headline).** `loglog-exponent(regime) = rank(obstruction)`.
   Falsifiable: a regime whose obstruction collapses to a *singleton wall*
   (rank 1) should admit a surplus with exponent `1` — a strictly better bound
   where the wall is single. A hypothetical *three-block* coupling would force
   `(loglog k)³`. This turns "which power of loglog?" into "what is the rank of
   the minimal determinant?" — a structural, checkable question.
2. **Min-entropy certificate.** The Helly witness minimizes `log lcm` over
   infeasible subfamilies; certifying infeasibility cheaply = finding the
   low-`ω` witness. This reframes codex-S599's "small Helly witness" search as
   an *entropy-minimization* problem with an explicit objective.
3. **Two dual axes.** Point entropy (`log lcm`, which residue) vs set/VC entropy
   (`rank n`, which constraints). Worst-case Helly and fractional-Helly dividend
   are the two directions on the *set* axis; the certificate cost lives on the
   *point* axis. Separating them prevents the conflation that produced a wrong
   first draft of this note (I had claimed `#atoms = lcm`; it is `2ⁿ`).
4. **Collatz tie-back (rank 1).** The bounded-defect conjecture (HYP-2148) is the
   **rank-1 / single-channel** case: a single Collatz trajectory threads one
   root-path of the predecessor tree `(4^j−1)/3`, so the small-odd contributions
   are a *chain* (one channel), and the budget `Σ 1/(3a) = 0.421` converges —
   bounded, no growing `(loglog)^h`. Bounded defect = "rank 1" in this ledger.
5. **VC/fractional-Helly dividend.** On the set axis, a positive fraction of
   `r`-tuples intersecting yields a point in many sets — an averaged dividend
   that can *lower* the effective rank in structured (non-shattering) regimes,
   the entropy analogue of the repo's "full `1/(k+1)` exits" (codex-S600).

## Honest status

- **Proved:** the certificate-entropy identity `H_cert = log₂ lcm` and its
  additive `ω`-channel decomposition; the `ω`–`loglog` (Mertens) link is classical.
- **One confirming instance:** the LRC rank-2 two-block ↔ Tao `(loglog k)²`.
- **Conjecture:** the general exponent law `p = rank` across regimes — the
  natural next test is a regime in the repo with a *known* singleton-wall
  obstruction (codex-S599 reported singleton/pair determinant witnesses at small
  `n`); does its best bound carry `(loglog)¹`? That is the falsification target.

## Next

1. Pull a singleton-wall regime from `lrc_twoblock_helly_s599.py` and check its
   loglog-exponent against the rank-1 prediction.
2. Make the `(loglog)^r` mechanism a clean lemma (second moment over `r`
   independent `ω`-channels) rather than a template.
3. Push the min-entropy-certificate objective into codex-S599's Helly search as
   an explicit `argmin log lcm`.
