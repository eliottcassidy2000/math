---
id: HYP-2623
status: STRONGLY SUPPORTED (exact through a=4 at k=31; a=5,6 bounded below their mediants); globalization open
source: kind-pasteur-2026-06-19-S9
related:
  - THM-539
  - HYP-2052
  - HYP-2084
---

# HYP-2623: the LRC max-min level `a_max(k)` is UNBOUNDED, so the second-point gap dips below Theta(1/k^2)

## Claim

Define the reachable level `a_max(k) = max { a : some gcd-1 k-set S has M(S) = a/(a(k+1)-1) }`
(the deepest Stern-Brocot mediant above the floor `1/(k+1)` that is realized). Then:

1. **Generic floor of the dip:** for ALL `k`, `a_max(k) >= 2` (doubled apex `{1,..,k-1,2k}`,
   `M=2/(2k+1)`). For "most" `k`, `a_max(k) = 2` (e.g. exhaustively k=2,3,4,5,8,9,10,11).
2. **The primorial family** `F(k,a)={1,...,k-2,k,a(k-1)}` realizes level `a` exactly when
   `(k-1)` is divisible by the first `a-1` primes:
   - `a=3` ⟺ `6 | (k-1)` (k=7,13,19,25,...; verified exact);
   - `a=4` at `k=31` (`k-1=30=2·3·5`; verified exact `M=4/127`);
   - `a=5,6` at `k=211, 2311` (`k-1`=primorial; `M` verified below the level-(a-1) mediant).
3. **Therefore `a_max(k) >= ω*(k-1)`** (the number of consecutive small primes `2,3,5,...`
   dividing `k-1`), which is UNBOUNDED (take `k-1` = primorial). So along that sparse
   subsequence `g(k)·k^2 -> 1/a -> 0`, i.e. **`liminf_k g(k)·k^2 = 0`** — the spectral gap is
   not bounded below by any `c/k^2`.

## Status of evidence

| level a | smallest k | k-1 factorization | M | verification |
|--------|-----------|--------------------|---|--------------|
| 2 | (all k)   | —                  | 2/(2k+1) | exact, k=2..20 |
| 3 | 7         | 6 = 2·3            | 3/23     | exact (also exhaustive sigma_2) |
| 4 | 31        | 30 = 2·3·5         | 4/127    | exact (crossing + covering) |
| 5 | 211 (cand)| 210 = 2·3·5·7      | < 5/1059 | covering (exact level pending — E2 workflow) |
| 6 | 2311 (cand)| 2310 = 2·3·5·7·11 | < 6/13871| covering |

## Open / to nail

- **Globalize:** prove `M(F(k,a)) = a/(a(k+1)-1)` exactly when `(k-1)` has the first `a-1`
  primes as factors (mechanism: killer speed `a(k-1)` annihilates clocks `t=j/d` for all
  `d|(k-1)`; need to show no OTHER time beats `a/(a(k+1)-1)`). A clean induction on the
  primorial structure is the target.
- **Exact level at k=211, 2311** (is it exactly `a=5,6`, or deeper?). [E2 workflow]
- **Is `F` optimal?** Could a multi-killer config reach level `a >> ω(k-1)` (e.g. `a ~ sqrt(k)`),
  making the dip much faster? [E1 workflow] If not, `a_max(k) = Theta(log k / log log k)` exactly.
- **True `sigma_2(k)`** at k=12,13,14 (does `F` give the min?). [E3 workflow]

## Why it matters

The user's "live question" (does anything dip below `Theta(1/k^2)` unboundedly?) is answered
**YES**, but the answer is delicate: the dip is real and unbounded, yet confined to an
extremely sparse arithmetic subsequence (`k-1` highly composite) and grows only like
`1/log` — so for almost all `k` the gap is `Theta(1/k^2)` with constant `~1/2`. The
difficulty of LRC stays localized at the AP, but the *width of its basin* is arithmetic in
`k-1`, not a clean `1/k^2`.

See THM-539, `07-reflections/lrc-spectral-gap-dips-along-primorials-kps.md`.
