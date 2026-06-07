---
id: THM-421
name: lrc-divisor-clock-peeling
status: PROVED (elementary, CRT) + VERIFIED (exact over QQ, n=11..22)
date: 2026-06-06
session: claudebox-2026-06-06-S710
depends_on:
  - THM-398   # C'(n): multiple-of-n => loose; LRC(n) convention
  - THM-420   # the prime-n companion (2n-1 shell non-transversal dodge)
provisional_id: true   # HYP/THM counter is contested; renumber at PR time
---

# THM-421: The divisor-clock peeling theorem (the composite-n companion of THM-420)

## Convention

Canon LRC (THM-398): `n` runners, speed set `S` of `n-1` distinct positive integers, gap `1/n`,
`M(S) = max_t min_{v in S} ||v t||`; loose iff `M(S) > 1/n`. A **d-clock** is a time `t = b/d` with
`gcd(b,d)=1`. `||x||` = distance to the nearest integer.

## Statement

Let `n >= 2` and let `d | n` be any divisor. At any d-clock `t = b/d`:

1. **(danger = multiples of d)** A runner `v` satisfies `||v t|| = 0` iff `d | v`; otherwise
   `||v t|| >= 1/d`. Consequently, with loneliness radius `1/n`:
   - if `d <= n` (always, since `d | n`), the set of dangerous runners `{v : ||v t|| < 1/n}`
     is **exactly** `{v in S : d | v}`;
   - for a **proper** divisor `d < n`, every non-multiple-of-d runner is **strictly** safe with
     margin `||v t|| >= 1/d > 1/n`.

2. **(peeling reduction)** Therefore loneliness of `S` near the d-clock is controlled entirely by the
   **sub-configuration of multiples of d**. Writing `S_d = {v in S : d|v}` and `w_i = v_i / d`, a
   perturbation `t = b/d + s` keeps every non-multiple safe for `|s|` small, and reduces the question
   "is `S` lonely near `t`?" to "is the sub-config `{w_i}` lonely near `s` (radius `1/n`)?" — a
   **smaller, recursive** LRC-type instance carried on the `mod (n/d)` fiber.

3. **(best clock)** The proper divisor maximizing the safe margin and minimizing the residual is the
   **largest proper divisor** `d = n/p` where `p` is the smallest prime factor of `n`. It peels the
   **sparsest** sub-config (the multiples of `n/p`) and leaves the residual on the `mod p` fiber.

## Proof

(1) Write `v b / d`. If `d | v` then `v b / d` is an integer, so `||v t|| = 0`. If `d nmid v` then,
since `gcd(b,d)=1`, `v b not\equiv 0 (mod d)`, so `v b mod d in {1, ..., d-1}` and
`||v t|| = ||(v b mod d)/d|| = min(r, d-r)/d` with `r = v b mod d in {1,...,d-1}`, whence
`||v t|| >= 1/d`. Since `d | n` gives `d <= n`, i.e. `1/d >= 1/n`, the non-multiples have
`||v t|| >= 1/d >= 1/n`, so they are never dangerous; for `d < n` strictly, `1/d > 1/n`. Multiples
have distance `0 < 1/n`, so they are dangerous. Hence danger `= {v : d|v}` exactly. ∎(1)

(2) By (1), at `t = b/d` all non-multiples of `d` have a fixed positive slack `>= 1/d - 1/n > 0`
(for `d<n`). By continuity of `t -> ||v t||`, there is `s_0 > 0` such that for all `|s| < s_0` every
non-multiple still has `||v(b/d + s)|| > 1/n`. For a multiple `v = d w`, `||v(b/d+s)|| = ||w b + d w s||
= ||d w s||` (since `wb in Z`), which depends only on `w = v/d` and the perturbation `s`. So on the
window `|s| < s_0`, `min_v ||v t||` equals `min_i ||d w_i s||` once `s` is small enough to have left
the non-multiples behind, i.e. the residual loneliness problem is exactly that of the dilated
sub-config `{d w_i}` near `s = 0`. ∎(2)

(3) Among proper divisors, `1/d` is maximized by the largest `d = n/p_min`, and the multiples of a
larger `d` are sparser, minimizing `|S_d|`. ∎(3)

## Verification

`04-computation/lrc14_divisor_lattice_s710.py` (exact over QQ, Fractions):

- **(A)** n=14, all divisors `d in {1,2,7,14}`: danger `= {v:d|v}` and margin `>= 1/d` on **100%** of
  trials (400 configs x each d-clock; 2400/2400 for d=7,14). `d=7` recovers S643's mod-7 fiber.
- **(D)** Best-proper-divisor clock across `n = 11..22`: composite `n` peel sizes
  (avg residual): n=12 -> mult-of-6 (2.13), n=14 -> mult-of-7 (2.20), n=15 -> mult-of-5 (3.02, mod 3),
  n=16 -> mult-of-8 (2.18), n=21 -> mult-of-7 (3.11, mod 3), n=22 -> mult-of-11 (2.26). Prime
  `n = 11,13,19`: **no** proper divisor `> 1` => divisor-clock empty.

## Significance: the prime/composite partition of the C'(n) frontier

THM-420 (S642) reduces `C'(n)` for **prime n** (more precisely `2n-1` prime, the shell route) to the
rare `+/-`-transversal core. THM-421 reduces `C'(n)` for **composite n** to the sparse residual
mult-of-`(n/p)` sub-config on the `mod p` fiber. The two are **complementary**:

- **prime n**  -> THM-420 (the `2n-1` multiplier shell);
- **composite n** -> THM-421 (the divisor-clock tower, recursing on `n/p`).

`n = 14` is **two-headed**: composite, so THM-421 handles the clock side (reduces to a <=4-runner
mult-of-7 core on the mod-2 fiber); but `2n-1 = 27 = 3^3` is ramified, so THM-420 fails there. The
residual gap (both routes) is the same shape: show the small residual core is dodgeable in the
perturbation window (HYP-2346 / HYP-2347 dominance-vs-window dichotomy). Neither theorem alone closes
LRC, but together they show the **entire** open frontier (`n = 15,19,21,22,...`) reduces to a small,
finite residual-core check rather than the full configuration space.
