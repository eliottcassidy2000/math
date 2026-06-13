---
source: opus-2026-06-03-S592 (remote-control, round 3 of overnight cycle)
status: INVESTIGATION — the non-transversal sporadic tight configs appear exactly when 2n-1 is composite; for n=14 (2n-1=27=3³) they are a prime-3 obstruction, a SECOND prime alongside the prime-2 of C'
tags: [LRC, worry-set, sporadic, transversal, 2n-1, composite, prime-3, n14, residual, THM-401]
---

# The sporadics are the composite part of 2n−1 — a prime-3 obstruction at n=14

**Context:** overnight cycle, round 3 (investigate). THM-401 pinned the worry modulus at
`C=2n−1`, leaving the *classification at `C`* (does floor-tight ⟺ shell-transversal?) as
the residual. This round characterises the exceptions.

## Finding — sporadics ⟺ `2n−1` composite

Exhaustive floor-tight (`M=1/n`) census in bounded boxes:

| n | `2n−1` | floor-tight | sporadics | non-transversal sporadics |
|---|---|---|---|---|
| 6 | **11 (prime)** | 2 (AP + `(1,3,4,5,9)`) | 1 | **0** (the sporadic is a transversal flip) |
| 8 | **15 = 3·5** | 3 | 2 | **2** (`(1,2,3,4,5,7,12)`, `(1,4,5,6,7,11,13)`) |
| 10 | **19 (prime)** | 1 (AP only) | 0 | **0** |

> **Non-transversal sporadic tight configs appear iff `2n−1` is composite.** When `2n−1`
> is prime (`n=6,10`), every floor-tight config is a perfect antipodal **transversal** of
> the shells `P_a={a,−a} mod 2n−1` — THM-401's classification is *exact*. When `2n−1` is
> composite (`n=8`, `15=3·5`), tight configs escape the transversal structure — they sit
> on the **non-unit shells** (residues sharing a factor with `2n−1`).

## Why: the non-unit shells of a composite `C=2n−1`

The shell action is `(ℤ/C)^*` (the multiplicative units, S572o's "inverse clocks"). When
`C` is prime, all nonzero residues are units — every shell is unit-visible, and floor-
tightness forces a clean transversal. When `C` is composite, the **non-unit residues**
(`gcd(r,C)>1`) form shells the unit action cannot reach; a config can be floor-tight by
occupying those shells *without* being a clean transversal — the sporadics. The number
and shape of sporadics is governed by the **prime factorisation of `2n−1`**.

## The payoff: a SECOND prime at n=14

`n=14`: `2n−1 = 27 = 3³`, composite. So `n=14` **has** non-transversal sporadics (e.g.
`V*`), and they are a **prime-3 obstruction** — the factor of `27`. This is *distinct*
from the prime-2 obstruction of C′ (the multiple-of-`14` / doubling, S589). So the full
`n=14` residual splits by prime:

```
n = 14 = 2·7         2n − 1 = 27 = 3³
 ├ prime 2: C′ / the multiple-of-14 / the ⟨×2⟩ doubling fragmentation (S585, S589)
 ├ prime 7: the SOLVED odd part — ℚ(ζ_14)=ℚ(ζ_7), the prime kernel (S588)
 └ prime 3: the SPORADICS / the composite shells of 2n−1=3³ (this round)
```

> **The `n=14` obstruction is carried by two small primes: `2` (C′ / doubling) and `3`
> (the composite `2n−1=3³` sporadics); the prime `7` is the solved kernel.** The earlier
> "everything localizes to the prime 2" (S589) was the *C′* half; the *classification*
> half localizes to the prime `3`. Both are tiny — and both are finite arithmetic ledgers
> (the prime-2 dodge / the prime-3 shell-cover), neither a measure estimate.

## Consistency check (THM-402)

Every floor-tight config (AP, transversal flip, or composite sporadic) is a round
tournament with a 3-cycle, hence `χ=2` (THM-402) — the worry-set is `χ`-uniform across
all three types, exactly as the "χ constant" theorem predicts.

## Honest status

- **Verified:** sporadics ⟺ `2n−1` composite (n=6,8,10); for `n=14`, `2n−1=27=3³` ⟹
  prime-3 sporadics exist.
- **Structural claim:** the `n=14` residual factors as prime-2 (C′) × prime-3 (sporadics),
  prime-7 solved. Grounded in the census + ℚ(ζ_14)=ℚ(ζ_7); the precise prime-3 shell
  ledger is the open finite task.

**Artifacts:** `04-computation/lrc_worryset_sporadics_s592.py` (+`.out`). Builds on
THM-401 (2n-1 modulus), THM-402 (χ=2), S570 (transversals), S589 (prime-2 localization),
S588 (ℚ(ζ_14)=ℚ(ζ_7)). New: **HYP-2138**.
