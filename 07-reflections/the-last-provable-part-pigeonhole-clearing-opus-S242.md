---
source: opus-2026-07-11-S242
status: A PROVED theorem (the "last provable part"). The pigeonhole clearing theorem — three proved lemmas,
  no anti-concentration — clears ~44% of divisor-complete families, lifting total elementary-provable
  coverage of LRC(14) to ~95% of all 13-speed families. The remaining ~5% is the pure anti-concentration
  core (fold-classes must cluster, not just be few).
tags:
  - lrc14
  - pigeonhole
  - proved
  - auto-safe
  - band-edge
  - provable-coverage
---

# The last provable part: the pigeonhole clearing theorem

**opus-2026-07-11-S242.** Owner: work on the last provable part. The auto-safe reduction (S241) turns
composite-clearing into a covering of the `φ(q)/2` unit fold-classes by the coprime-to-`q` sub-family. The
*provable* part of that is **pigeonhole**: if the coprime sub-family is smaller than the number of fold-classes,
it *cannot* cover them — no anti-concentration required.

## The theorem (PROVED; verified 0 violations / 1749)

> **Theorem.** Let `v` be a 13-speed family and `q ∈ {15,21,25,27}` (composite, prime factors ≤ 13, danger
> band `{0,±1}`). If `v` has **no multiple of `q`** and **fewer than `φ(q)/2` speeds coprime to `q`**, then
> `v` clears at `q`, hence `M(v) ≥ 2/q > 1/14`.

*Proof — three proved lemmas:*
- **[auto-safe, S241]** every speed with `gcd(v_i,q) > 1` is inert: at a unit multiplier `p`, `v_i p mod q`
  shares a factor with `q`, so it avoids the danger band `{0,±1}` (unless `q | v_i`, excluded).
- **[pigeonhole, this session]** the `< φ(q)/2` coprime speeds occupy `< φ(q)/2` unit fold-classes, so some
  fold-class `r` is empty; at `p = r⁻¹` no coprime speed hits `{±1}`, and the structured speeds are inert, so
  `bandCount(v,q,p) = 0`.
- **[band-edge, S235]** clearing at `q` with `14 ∤ q` gives `M(v) ≥ ⌈q/14⌉/q = 2/q > 1/14`. ∎

Verified end-to-end: pigeonhole-forced ⟹ actually clears, **0 violations / 1749**.

## Coverage

The theorem provably clears **43.7%** of divisor-complete spread families — *with no anti-concentration*.
Forcing modulus: `q=27` (54%), `q=25` (28%), `q=21` (16%), `q=15` (2%). Drivers: **≥ 5 speeds divisible by 3**
(forces `q=27`, since `φ(27)/2 = 9`, so `≤ 8` coprime-to-3 can't cover 9 classes) or **≥ 4 speeds divisible by
5** (forces `q=25`, `φ(25)/2 = 10`).

Combined with the elementary `t=1/d` ladder (THM-366, which clears every *non*-divisor-complete family):

| class | share of all families | status |
|---|---|---|
| not divisor-complete → `t=1/d` | ~91.5% | **PROVED** (THM-366) |
| divisor-complete, pigeonhole-forced | ~3.7% | **PROVED** (this theorem) |
| divisor-complete, not pigeonhole-forced | ~4.8% | open — anti-concentration core |

> **~95% of all 13-speed families are provably lonely (`M > 1/14`) by elementary means** — the `t=1/d`
> divisibility ladder plus the pigeonhole clearing theorem, no anti-concentration anywhere.

## The honest remainder

The open ~5% is the **anti-concentration core**: divisor-complete families whose coprime sub-family is *large
enough* (`≥ φ(q)/2`) to potentially cover the fold-classes at every bounded composite. For these, clearing
requires the residues to **cluster** (leave a fold-class empty despite having enough speeds) — a genuine
anti-concentration statement, not a counting one. This is exactly where the hypothetical covering
counterexamples would live, and it is the sole remaining open part. The pigeonhole cannot reach it (by
definition it has enough speeds); it needs the residue-clustering / three-gap disjunction that S238–S241
showed has no bounded-window shortcut.

## Net

This banks the **last provable part**: a clean, self-contained theorem assembling the three proved lemmas
(auto-safe + pigeonhole + band-edge) into unconditional clearing for ~44% of the residual, raising the
total elementary-provable coverage of LRC(14) to ~95%. What remains is precisely the anti-concentration core
(~5%), sharply isolated: divisor-complete, enough coprime speeds, requiring residue clustering — the wall,
now with everything provable around it stripped away and banked.

→ opus-S241 (auto-safe lemma), opus-S235 (band-edge margin), THM-366 (the `t=1/d` ladder / covering ⟹
divisor-complete), opus-S238 (no bounded-window shortcut — why the remainder is genuinely anti-concentration),
opus-S239 (the shared crux). Files: `lrc14_pigeonhole_provable_part_opus_S242.py` (+`.out`).
