# The windows around 10007, 10009, 30011, 30013, 30021, 30027, and the theorems where an approximation selects the integer

Session klein-2026-11-16-force. Receipts `numbers.out`, `tribonacci.out`.

## 1. The owner's numbers (PROVED arithmetic)

| number | factorisation | reading |
|---|---|---|
| 10003 | 7 · 1429 | a hidden generator `7p` (THM-5258 (B)) |
| 10007, 10009 | twin primes | |
| 30009 | 3 · 7 · 1429 | the midpoint `3·10003`, digit 0 |
| **30011, 30013** | **twin primes** | `3·10003 + 2` and `3·10003 + 4`: the digit-1 and digit-2 children of the hidden generator 10003 |
| 30021 | 3 · 10007 | the midpoint of 10007 (forced by `C_3 → T`) |
| 30023 | 7 · 4289 | `3·10007 + 2`, a hidden generator |
| 30025 | 5² · 1201 | `3·10007 + 4` |
| 30027 | 3 · 10009 | the midpoint of 10009 |
| 30029 | prime | `3·10009 + 2` |
| 30031 | 59 · 509 | `3·10009 + 4` |

So the three five-windows `{3q − 4, …, 3q + 4}` of the consecutive odd numbers `q = 10003, 10005, 10007, 10009`
overlap (`30023 = 3·10007 + 2 = 3·10009 − 4`), and the two twin-prime pairs sit in them as children: the twin primes
30011, 30013 are the `+2, +4` children of the hidden generator 10003, while the twin primes 10007, 10009 are parents
whose midpoints 30021, 30027 are forced and whose `+2` children 30023 (hidden generator) and 30029 (prime) are the
cases §3 of FORCING.md is about. In the address `t = (H − 1)/2`: `10003 ↦ 5001 = 7·714 + 3`, i.e. the image of the
prime address `714 = (1429 − 1)/2`.

## 2. Where an approximation selects the integer (PROVED)

- **One row is Tribonacci, exactly rounded.** The counts of a single Collatz row built as a one-parent chain (wave 46;
  compositions with no two adjacent 1s) satisfy `C(m) = C(m−1) + C(m−2) + C(m−3)` and
  `C(m) = round(c·τ^m)` for every `m ≥ 1`, with `τ = 1.839286755…` the tribonacci constant and
  `c = 0.435616389351097…` (`tribonacci.out`, checked to `m = 40`; pre-registered P-E **confirmed**). The real function
  `c·τ^m` is within `1/2` of the count, so it selects it.
- **The selection principle.** If a construction yields a count `H` with `H ≡ N (mod M)` and `|H − N| < M`, then
  `H = N`. Hadamard's one-apex placement (`one_apex_all_residues_20261116`) gives, for fixed `M`, `H ≡ N (mod M)` with
  `|H − N| = O_M(√N)`; the principle says that Chan–Pak above `N_0` would follow from a placement with error `≤ C√N`
  and residue control modulo some `M > 2C√N`, with `C` not depending on `M` — which is exactly what that packet
  states is not provided (its constant grows with `M`). This is the honest form of "close enough to select the
  integer": the approximation must be uniform in the modulus.
- **The windows select the digit.** For `3q + 2σ` with `σ ∈ {0, 1, 2}`, the residue modulo 9 determines `σ` once
  `q (mod 3)` is known (`3q + 2σ ≡ 3(q mod 3) + 2σ`), so within a five-window the count's residue mod 9 reads the
  digit — the mod-9 clock of TWO_APEX_HOST (2 is a primitive root mod 9) is the selector.
