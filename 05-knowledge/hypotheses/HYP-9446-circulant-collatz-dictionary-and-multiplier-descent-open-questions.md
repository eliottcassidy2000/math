---
id: HYP-9446
title: "After THM-4848/4849. (1) REFUTED (pre-registered): locked tournament cycles are rotational. It holds for q <= 19 and fails at q = 23, where 14 of the 18 locked cycles are not rotational. (2) OPEN, Ailon-Rudnick type: infinitely many primes q have an interval cycle locked at no prime, i.e. gcd(2^((q+1)/2)-1, 3^((q+1)/2)-1) = 1; this holds for 299 of the 782 primes below 6000. Locking at p happens iff lcm(ord_p 2, ord_p 3) | (q+1)/2, e.g. at 5 iff q = 7 (mod 8), which is the owner's 35. (3) OPEN, pre-registered and passed so far: the Paley determinant beta_k(P_q) in Q(sqrt(-q)) has no split odd prime in its rational content, i.e. it is one-sided at every split odd prime. Passed at P7 (11, 5861), P11 (5) and P19 (5); 2 is two-sided. (4) OPEN: at most two tournament base points per cycle for every odd q (holds for q <= 27)."
status: >
  OPEN (mac-mini-2026-10-29-cycle101, 2026-10-29), apart from (1), which is REFUTED, and the proved parts marked
  PROVED. PRE-AUDIT checkpoint: independent audit pending. Evidence is FINITE-EXACT unless marked otherwise.
  Scripts: 04-computation/experiments/cycle101_tournament_collatz_20261029/ (locked_tournaments.py, chirality.py,
  composite_basepoints.py, locking_iff.py).
session: mac-mini-2026-10-29-cycle101
source: 05-knowledge/results/cycle101_tournament_collatz_dictionary_20261029.md
related:
  - 01-canon/theorems/THM-4848-multiplier-descent-for-run-observers-of-circulant-tournaments-paley-quadratic-norms-and-the-multiplicity-law.md
  - 01-canon/theorems/THM-4849-circulant-tournaments-are-rational-collatz-cycles-shift-decimation-reversal-two-base-points-and-the-locking-theorem-behind-35-and-175.md
---

# HYP-9446 — the circulant–Collatz dictionary and multiplier descent: open questions

Terms are those of THM-4849.
- A *tournament cycle* is a shape-`((q−1)/2, q)` rational cycle with a tournament base point.
- A cycle is *locked at `p` under `v`* if every marking cancels `p` on appending `v`. By THM-4849 (E) this happens
  iff every letter is `≡ v (mod ord_p(2))`, given `p | D(wv)` and `p ∤ 6(2^v − 3)`.

## 1. Locked tournament cycles are rotational — REFUTED (pre-registered)

**The prediction**, written into `locked_tournaments.py` before running it: every locked tournament cycle is a
cycle of `R_q`.

**The test.** Every tournament cycle for prime `q ≤ 23`, every prime `p` with `ord_p(2) ≥ 3` dividing a difference of
letters, and every appended letter class.

**The outcome.**

| `q` | tournament cycles | locked | not rotational |
|---:|---:|---:|---:|
| 7 | 5 | 1: the interval cycle, at 5 | 0 |
| 11 | 27 | 2: the interval and `6·I`, both at 7 | 0 |
| 13 | 58 | 0 | 0 |
| 17 | 248 | 0 | 0 |
| 19 | 503 | 1: interval, at 11 | 0 |
| 23 | 2037 | 18 | **14** |

So the prediction FAILS at `q = 23`. For example `(1,1,1,1,1,1,1,4,1,1,10)`, a non-rotational tournament word, has
all letters `≡ 1 (mod 3)` and is locked at 7. Locking is a congruence condition on the gaps; the tournament class
is invisible to it.

**What survives** (PROVED, THM-4849 (E)). The interval cycle of `R_q` is locked under `(1)` exactly at the primes of
`gcd(2^m − 1, 3^m − 1)`, where `m = (q+1)/2`.

## 2. Infinitely many unlocked interval cycles (OPEN; a prime-restricted Ailon–Rudnick question)

**PROVED.** For a prime `p ≥ 5`, the interval cycle of `R_q` is locked at `p` under `(1)` iff
`lcm(ord_p 2, ord_p 3) | (q+1)/2`. By Dirichlet, each such `p` locks infinitely many `R_q`:
- `p = 5` iff `q ≡ 7 (mod 8)`: `q = 7` (the owner's 35), 23, 31, 47, 71, 79, …;
- `p = 7` iff `q ≡ 11 (mod 12)`;
- `p = 13` iff `q ≡ 23 (mod 24)`;
- `p = 11` iff `q ≡ 19 (mod 20)`.

**Conjecture.** For infinitely many primes `q`, `gcd(2^((q+1)/2) − 1, 3^((q+1)/2) − 1) = 1`, so the interval cycle is
locked at no prime.
- *Evidence.* This holds for 299 of the 782 primes `q < 6000`: the first are 3, 5, 13, 17, 29, 37, 41, 53, 61, 67, 73,
  89, 97, 101, 113. The pattern "`m` odd implies `gcd = 1`" holds below 109 and fails there (`m = 55`, `gcd = 23`).
- *Literature.* Over all `n` (not restricted to `n = (q+1)/2` with `q` prime), the statement
  `gcd(2^n − 1, 3^n − 1) = 1` infinitely often is the Ailon–Rudnick conjecture for `(a, b) = (2, 3)`. It is open:
  Ailon–Rudnick, Acta Arith. 113 (2004) 31–38, prove the function-field analogue. Bugeaud–Corvaja–Zannier (Math. Z.
  243 (2003)) bound the gcd by `exp(εn)`.

## 3. Paley one-sidedness (OPEN; pre-registered, passed so far)

**Setting.** For Paley `P_q` (`q ≡ 3 mod 4`) write `β_k = det B_1(M_k(P_q)) = a + bη` in `Z[η]`, where
`η = (−1 + √−q)/2` (THM-4848 (iv)). Then `min(v_𝔭(β), v_𝔭̄(β)) = v_p(gcd(a, b))` at a split prime `p`.

**Conjecture**, written into `chirality.py` before computing it. `gcd(a, b)` has no odd prime factor that splits in
`Q(√−q)`. Equivalently, at every split odd prime the nontrivial part of the run-observer cokernel lives over a
single prime of `Q(√−q)`. This is CLOCK_MINIMAL_PAIRING C5's "one cubic, zero of its dual".

**Evidence** (`chirality.out`).

| case | `gcd(a, b)` | split odd primes with exponent ≥ 2 in `N(β)` |
|---|---|---|
| P7, k = 2 | `2^6·3` (3 inert) | none |
| P7, k = 3 | `2^9·5^2` (5 inert) | 11: `(0, 2)`; 5861: `(2, 0)` |
| P11, k = 2 | `2^6·11` (11 ramified) | 5: `(2, 0)` |
| P19, k = 2 | `2^4` | 5: `(0, 2)` |

- *Strength.* All four non-automatic instances are one-sided. Primes with exponent 1 are automatically one-sided.
  Heuristically, one-sidedness at an exponent-2 prime has probability about 2/3, so four passes are weak evidence:
  about `(2/3)^4 ≈ 0.2` by chance.
- *The exception at 2.* 2 splits in `Q(√−7)` and is two-sided there: `(8, 6)` for `k = 2`, `(9, 13)` for `k = 3`,
  because `2^6` and `2^9` divide `gcd(a, b)`.
- *Next test.* P23 (`k = 2`, a 506 × 506 observer, 22 × 22 blocks over `Q(ζ_23)`), and P7, P11 with `k = 3, 4`.

## 4. At most two base points for every odd `q` (OPEN)

**PROVED (THM-4849 (B')).** Two base points `S`, `S + t` force `T(S) = C[R_(q/g)]`, `g = gcd(t, q)`.

**Conjecture.** A shape-`((q−1)/2, q)` cycle has at most two tournament base points for every odd `q`, not only for
primes.

**Evidence** (`composite_basepoints.out`).

| `q` | base points per cycle (count: cycles) | doubly-based cycles that are not `m·{1..k}` |
|---:|---|---:|
| 9 | `{1: 6, 2: 5}` | 2 |
| 15 | `{1: 104, 2: 12}` | 8 |
| 21 | `{1: 984, 2: 20}` | 14 |
| 25 | `{1: 4060, 2: 18}` | 8 |
| 27 | `{1: 8130, 2: 31}` | 22 |

Every pair `(S, S + t)` has the lexicographic-product structure of (B'). Three base points `S, S + t, S + u` would
make `T(S)` such a product with respect to both `⟨t⟩` and `⟨u⟩`. Whether that is impossible is the open step.

## 5. Not conjectured (typed; recorded so that it is not re-derived)

- *Reduced denominators of the Paley pair.* `5 | D_q` iff `q ≡ 3 (mod 8)`. At `q = 11` only the QR cycle is
  5-integral (reduced denominator 361), at `q = 19` only the NQR cycle, and at `q = 43` neither. This is consistent
  with a `1/5` chance per cycle, so GENERIC.
- *Consecutive primes.* The five (3,7) minima 19, 23, 29, 31, 37 are consecutive primes: NUMEROLOGY.
- *Cycle traces.* Decimation does not preserve them (R7's three cycles: 267, 177, 155 over 101). Reversal does
  (five-mirrors Proposition 6).
