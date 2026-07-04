---
id: HYP-4054
title: THE CAPACITY ARGUMENT -- the creative reason the log-census kernel closes. Blocking a modulus q (making a covering family have NO lonely witness at any a/q) requires the 13 speeds to sit in a "no-witness" residue configuration of density f_q < 1 (bounded away from 1, verified f_q in [0.0008, 0.14]). By CRT the blocking constraints at distinct primes are independent, so blocking primes p_1..p_r requires the 13-tuple of speeds to hit a joint config of density prod f_{p_i}; the smallest such tuple has size ~ e^{(1/13) sum log(1/f_{p_i})}, so a blocking tuple with max|v| <= M exists only if sum_{blocked p} log(1/f_p) <= 13 log M. Since log(1/f_q) >= c > 0 for every q, AT MOST O(log M) primes can be blocked, hence the first WITNESS prime q* = O(log M loglog M): loneliness at a bounded denominator follows because the speeds' MAGNITUDE M limits the INFORMATION they can carry, and each blocked modulus costs a fixed quantum of it.
status: CREATIVE ARGUMENT (mechanism verified numerically; the geometry-of-numbers step (smallest blocking tuple >= e^{sum log(1/f)/13}) and the f_q>=c-lower-bound are the two rigor gaps). Explains HYP-4040's q~3.6 ln M. Not a full proof.
source: mac-mini-2026-07-03-S28
related:
  - HYP-4040   # no uniform band / q ~ 3.6 ln M for compressed lcm-blockers -- THIS is the reason WHY
  - THM-515    # klein: L_C measure = singular series (the circle-method side; f_q is its finite-q shadow)
  - HYP-3982   # kps star-safe measure route (13/7 > 1 union bound = the same wall)
  - MISTAKE-098 # the weak-adversary trap this argument's f_q correctly captures
results:
  - 04-computation/resonance_kernel_macmini_20260703.py
  - 04-computation/capacity_fp_macmini_20260703.py
  - 05-knowledge/results/resonance_kernel_macmini_20260703.out
---

# HYP-4054 -- the capacity argument (the creative reason)

## The kernel, restated
Every covering `gcd=1` family of 13 speeds `<= M` is lonely at some `a/q`. The empirical law (HYP-4040) is
`q* ~ 3.6 ln M` -- unbounded but logarithmic. The open kernel: *why does a witness appear by `q* = O(log M)`,
despite the union bound failing (`13 * (2/14) = 13/7 > 1`)?* Here is the reason.

## The character sum and `f_q`
The number of lonely numerators at denominator `q` is
```
    N(q) = #{ a in (0,q) : min_i ||v_i a / q|| >= 1/14 }
         = q * (6/7)^13  +  q * E(q),   E(q) = sum_{m != 0, q | sum m_i v_i} prod_i shat(m_i),
```
with `shat(m) = -sin(pi m/7)/(pi m)` (`m != 0`), `shat(0) = 6/7` (the safe-arc Fourier coefficients). The main
term `q (6/7)^13 > 0` is the count you'd get if the 13 danger arcs were **independent**.

Define `f_q = ` (fraction of 13-tuples of nonzero residues mod `q` that are **no-witness**, `N(q)=0`). This is
the finite-`q` shadow of the singular series (THM-515). **Verified numerically:** `f_q in [0.0008, 0.14]` for
`q = 17..113` -- always `< 1`, and `>= ~10^-3`, so `log(1/f_q) in [2, 7]`, bounded below by `c ~ 2`. As
`q -> infinity`, `f_q ->` the probability that 13 random arcs of measure `1/7` cover the circle, a fixed
constant in `(0,1)` -- so `log(1/f_q) >= c > 0` for **all** `q`.

## Why `f_q < 1`: no-witness = commensurability
`N(q) = 0` requires `E(q) = -(6/7)^13`, i.e. the error must kill the main term. Because `shat(m)` decays like
`1/|m|` and its large-`|m|` mass is bounded, the kill needs **small-`m` resonances** `q | (m_i v_i + m_j v_j)`
with `|m| <= K ~ 7`. Verified: for 57/57 generic compressed families, every no-witness modulus divides such a
small combination -- runners `i, j` are **commensurate** mod `q` (`v_i/v_j ≡ -m_j/m_i`), their danger arcs
LOCK, and 13 pairwise locks let the arcs TILE `Z/q`. Tiling is a *special* (density `f_q`) configuration; a
*generic* family leaves `(6/7)^13 = 0.135` of the circle safe -- a witness.

## The capacity bound (the creative core)
To block prime `p` (force no-witness), the speed vector `(v_1,...,v_13) mod p` must land in the no-witness set
`C_p`, of density `f_p`. The blocking conditions at **distinct primes are CRT-independent**, so to block
`p_1, ..., p_r` the vector mod `prod p_i` must land in `prod C_{p_i}`, of density `prod f_{p_i}`. A lattice of
density `d` in `[0, prod p)^13` has its smallest point at scale `~ (prod p) * d^{-1/13}`; a **blocking vector
with `max|v_i| <= M` exists only if**
```
    (1/13) * sum_{blocked p} log(1/f_p)  <=  log M,     i.e.   sum log(1/f_p)  <=  13 log M.
```
Since `log(1/f_p) >= c > 0`, **at most `r <= 13 log M / c = O(log M)` primes can be blocked**. The first
`O(log M) + 1` primes therefore contain an unblocked one -- a **witness prime** `q* <= p_{O(log M)} =
O(log M * loglog M)`. QED (heuristic).

The physical statement: *the speeds carry `~13 log M` bits (their magnitude); each blocked modulus spends a
fixed quantum `>= c`; so only `O(log M)` moduli can be blocked, and loneliness reappears by `q = O(log M)`.*
This is exactly HYP-4040's `q ~ 3.6 ln M` -- the `3.6` is `13/c` with `c ~ ln(1/max f_q) ~ 13/3.6 ~ 3.6`... the
constant is the mean blocking cost.

## Lemma (i) is EASY (`f_q` bounded away from 1, unconditionally)
Over a uniformly random speed vector mod `q`, the 13 danger events `{v_i a in danger}` at a fixed `a` are
INDEPENDENT (distinct `v_i` uniform), so the expected safe (witness) fraction is
```
    E[ N(q)/q ]  =  prod_i (1 - |danger|/q)  =  (6/7)^13  =  0.135   (as |danger|/q -> 1/7).
```
Since `N(q)/q in [0,1]`, Markov gives `P(N(q) > 0) >= E[N(q)/q] = 0.135`, hence
```
    f_q  =  P(no witness)  <=  1 - 0.135  =  0.865 < 1,   so  log(1/f_q) >= 0.145 =: c.
```
No circle method needed for `c > 0`. (Empirically `f_q` is far smaller, `<= 0.14`, so the true `c ~ 2`; but
`0.145` already suffices for the capacity bound `r <= 13 log M / c`.) **So lemma (i) is done -- the ONLY real
gap is (A) below / lemma (ii).**

## The single genuine gap: (A) no-witness => small resonance
Restated as the capacity's teeth: a prime `q` can be blocked (no-witness) ONLY IF `q | (bounded-order small
combination sum m_i v_i)` (`|m_i| <= K`, `<= t` nonzero terms, `t` bounded). VERIFIED necessary for 57/57
generic families (2-term already accounts for all). Granting it, the elementary count finishes:

**(B) [elementary]** the `<= (13 choose t)(2K+1)^t` small combinations are integers `<= t K M`, so their
DISTINCT prime divisors number `<= (# combos) * log(tKM) = O(log M)`. Hence at most `O(log M)` primes are
no-witness, and the first witness prime `q* <= p_{O(log M)} = O(log M loglog M)`. 

So the WHOLE compressed crux rests on the one circle-method lemma (A): *if `q` divides no small speed-
combination, a lonely witness exists at `a/q`.* That is the clean, isolated open kernel.

## What is rigorous, what is not
- **Verified:** `f_q < 1` and `>= c` (numeric); no-witness => small resonance (57/57 general); the whole
  picture reproduces `q* ~ 3.6 ln M`.
- **Rigor gaps (2):** (i) `f_q >= c` for ALL `q` uniformly -- needs the arc-covering-probability lower bound
  (`13` arcs of measure `1/7` fail to cover with probability `>= c`; plausibly `(6/7)^13`-ish but must be an
  a.s.-type bound, not just in-mean); (ii) the geometry-of-numbers step -- "smallest blocking vector
  `>= e^{sum log(1/f)/13}`" -- needs a genuine successive-minima / large-sieve bound, because a special
  blocking class could a priori have an atypically small representative.

## Why this matters
It converts the wide-residual crux from "prove the log-census bound" into two SEPARATE, classical-flavoured
lemmas: an **arc-covering lower bound** (probability 13 arcs miss a point `>= c`) and a **geometry-of-numbers
capacity bound** (blocking `r` independent congruence conditions of density `f_p` forces `max|v| >=
e^{sum log(1/f_p)/13}`). Both are the kind of statement the sieve / large-sieve literature handles. The
`13/7 > 1` union-bound wall (kps star-safe, klein measure floor, opus tower) is the SAME `f_q`: it is exactly
the fraction that survives, and the capacity argument says the adversary can only suppress it `O(log M)` times.
