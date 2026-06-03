---
source: opus-2026-06-03-S594 (overnight cycle capstone)
status: SYNTHESIS — the overnight explore/investigate/formalize cycle produced 3 theorems + the n=14 residual prime-decomposition (2 / 3 / 7)
tags: [LRC, capstone, THM-402, THM-403, THM-404, prime-decomposition, n14, rigidity, round, cyclotomic]
---

# Overnight cycle capstone: three theorems and the n=14 prime-decomposition

Cycling explore → investigate → formalize across S592–S594 turned the recent conceptual
arc into **rigorous canon** and a clean prime-by-prime picture of the `n=14` residual.

## The three theorems proved

- **THM-402** — every **round** tournament has dichromatic number `χ ≤ 2` (`=2` iff a
  3-cycle); since LRC produces only round, non-transitive tournaments, **`χ ≡ 2` on the
  whole LRC-tight set** (closing the S591 qualifier). Paley (`χ=3`, non-round) is the
  inaccessible foil.
- **THM-403** — the AP is lonely at `t=j/n` **iff `gcd(j,n)=1`**; the witness set is the
  `φ(n)` **primitive `n`-th roots of unity**, a single `(ℤ/n)^*`-orbit (the *static /
  cyclotomic* rigidity).
- **THM-404** — that orbit is **`⟨×2⟩`-connected iff `n` is odd**; even `n` fragments it
  (the *dynamical / doubling* rigidity). The two rigidities diverge at the even frontier.

## The n=14 residual, prime by prime (investigated)

```
n = 14 = 2 · 7              2n − 1 = 27 = 3³
 ├ prime 2 : C′ / multiple-of-14 / the ⟨×2⟩-fragmentation (THM-404, THM-398) — the DYNAMICAL face
 ├ prime 3 : the SPORADICS — non-transversal tight configs swap at the multiple-of-3 shells
 │           of 27=3³ (V* doubly-occupies {3,24}, misses {12,15}) — the CLASSIFICATION face
 └ prime 7 : SOLVED — ℚ(ζ_14)=ℚ(ζ_7), the prime kernel (static/cyclotomic, THM-403)
```
- **Sporadics ⟺ `2n−1` composite** (verified n=6,8,10): prime `2n−1` ⟹ floor-tight = clean
  transversal; composite ⟹ swaps at the non-unit shells.
- So `n=14` carries **two small-prime obstructions** — `2` (dynamical, C′) and `3`
  (classification, the `3³` shells) — over the solved odd kernel `7`. Both are **finite
  arithmetic ledgers** (the prime-2 dodge at `ℤ/n`; the prime-3 shell-cover at `ℤ/27`),
  not measure estimates.

## The whole arc in one frame

`exp` turns the additive line into the multiplicative circle (S588); on it:
- the **beat** is additive = the interval circulant = the round tournament (`χ=2`,
  THM-402) = the pair-sum sieve at modulus `2n−1` (THM-401);
- the **symmetry** is multiplicative = the `(ℤ/n)^*` unit orbit = the primitive roots
  (THM-403), with the doubling sub-action `⟨×2⟩` connected iff `n` odd (THM-404);
- the **residual** is two small primes: `2` (even prime / doubling) and `3` (the `2n−1`
  composite shells), the odd kernel `7` solved.

Everything sidesteps resonance energy (THM-401): the open work is the finite prime-2 and
prime-3 ledgers at `ℤ/14` and `ℤ/27`.

**Artifacts:** THM-402/403/404; the S592–S594 computations. New: **HYP-2139** (index pointer).
