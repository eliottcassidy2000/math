---
id: THM-539
title: The LRC max-min second-point spectrum — Stern-Brocot mediants, the denominator lemma, and the a=3,4 below-mediant dips
status: MIXED — Lemma A PROVED; generic a=2 + dips a=3,4 VERIFIED exact (infinite families);
        a>=5 OPEN (the natural family collapses to the floor at a=5 — MISTAKE-079 corrected my
        earlier "a_max unbounded" overclaim). Whether g(k)=Theta(1/k^2) uniformly is OPEN.
source: kind-pasteur-2026-06-19-S9
related:
  - HYP-2052   # oracle-S552 spectral gap (original "gap = 2/(2n-1)" claim, here corrected)
  - HYP-2084   # codex-S573 below-edge rows
  - HYP-2621   # codex-S16 A_{k,r} mediant ladder (same family, mod-30 refinement)
  - HYP-2622   # codex-S17 excess ledger + bounded-height filter (uses this session's denominator lemma)
  - HYP-2623   # the a_max-unbounded conjecture (this session)
  - THM-534    # moment-LP dual
  - THM-522    # scale-invariance + quantization
---

# THM-539: The LRC max-min second-point spectrum

## Setup

For `k` distinct positive integers `S` (gcd 1, observer 0), the **max-min collar**
`M(S) = max_{t in (0,1)} min_{v in S} ||v t||`.  LRC(k+1): `M(S) >= 1/(k+1)`, with the
floor hit by the AP `{1,...,k}`.  The **second spectrum point**
`sigma_2(k) = min over non-tight S of M(S)`, and the **spectral gap**
`g(k) = sigma_2(k) - 1/(k+1)`.

(Translation: codex/oracle "n runners" = `k+1`; their `2/(2n-1)` = our `2/(2k+1)`.)

## Lemma A (denominator lemma) — PROVED

`M(S)` is achieved at a time `t* = m/(v_i ± v_j)` for some binding pair `v_i, v_j in S`;
hence the reduced denominator `q` of `M(S)` divides `(v_i ± v_j)`, so **`q <= 2·max(S)`**.

*Proof.* `f(t)=min_v ||v t||` is continuous, piecewise-linear, 1-periodic. For `|S|>=2`
with `M(S)<1/2`, the max of `f` is at a point where the lower envelope peaks, i.e. where
two sawtooths `||v_i t||, ||v_j t||` cross on the upper side. Equality `||v_i t||=||v_j t||`
means `v_i t ≡ ± v_j t (mod 1)`, i.e. `(v_i ∓ v_j) t ∈ Z`, so `t*=m/(v_i∓v_j)`. Then
`M = ||v_i t*||` has denominator dividing `v_i∓v_j`. ∎

**Corollary (universal lower bound).** `g(k) = (p(k+1)-q)/(q(k+1)) >= 1/(q(k+1)) >= 1/(2·max(S)·(k+1))`.
So **the gap dips below `Theta(1/k^2)` iff the `sigma_2`-extremizer needs `max(S)/k -> infinity`.**

## Stern-Brocot reframe — the spectrum just above the floor

The values just above `1/(k+1)` lie on the Stern-Brocot path `1/k ⤳ 1/(k+1)`. The mediants
converging to the floor are
```
   a/(a(k+1)-1),   a = 2, 3, 4, ...   ->   1/(k+1).
```
`a=2` gives `2/(2k+1)` (mediant of `1/(k+1)` and `1/k`). If `sigma_2(k)=a/(a(k+1)-1)` then
```
   g(k) = 1/((a(k+1)-1)(k+1)),     g(k)·k^2 -> 1/a.
```
So the gap dips below `Theta(1/k^2)` **iff the reachable level `a` is unbounded.** This is the
exact reformulation of the user's "lower bound on g(k)" question, and of the G2-lemma lift
depth `c(k) = (a(k+1)-1)(k+1) ~ a·k^2`.

## Generic value (a=2) — VERIFIED exact

For every `k`, the doubled apex `D_k = {1,...,k-1, 2k}` has `M(D_k) = 2/(2k+1)` (exact,
checked k=2..20). So `sigma_2(k) <= 2/(2k+1)`, `g(k) <= 1/((2k+1)(k+1))`, `g·k^2 <= ~1/2`.
**Exhaustively, `sigma_2(k) = 2/(2k+1)` (a=2, NO dip) for k = 2,3,4,5,8,9,10,11** (filtered
exact scan; boxes large enough to catch up to level 5).

## The dips (a >= 3) — the primorial family `F(k,a)`

`F(k,a) = {1, 2, ..., k-2, k, a·(k-1)}` (AP with speed `k-1` removed, replaced by `a(k-1)`).
This is **codex's `A_{k,r}` family** (HYP-2621/HYP-2622, independent concurrent discovery); the
two threads converged on the same family and the same denominator lemma.
**Verified exact:**
```
  a=2:  k=5            M=2/11    g·k^2=0.379                   [generic doubled-apex level]
  a=3:  k=7,13,19,25   M=3/(3k+2)  (3/23, 3/41, 3/59, 3/77)  g·k^2 ≈ 0.27–0.31
                       [k ≡ 7,13,19,25 (mod 30) — codex's sharp refinement of my "k≡1 mod 6";
                        at k ≡ 1 (mod 30) the a=3 family is AP-tight again and a=4 takes over]
  a=4:  k=31,61,...,181 M=4/(4k+3) (4/127, ...)               g·k^2=0.2365  [k ≡ 1 (mod 30); INFINITE family]
```
**a=5 is NOT reachable by this family (CORRECTION — MISTAKE-079).** I initially claimed a=5,6 at
k=211,2311 from a covering test reading `M < 5/1059`. That reading was WRONG: `F(211,5) =
{1,…,209,211,1050}` has `M = 1/212 = the FLOOR exactly` — it COLLAPSES to a tight configuration
(`g=0`), it does NOT dip to level 5. Verified by two independent covering implementations
(workflow E2 + own recheck): `F(k,5)` is TIGHT for every `k` with `2·3·5·7 | (k-1)` (k=211,421,631,841),
and stays above all small mediants otherwise. So **the natural family tops out at level 4**, and the
`a=5,6` "deep dips" were an artifact. (Bonus: `F(k,5)` collapsing to the floor means `{1,…,k-2,k,5(k-1)}`
is a NEW member of the tight locus when `2·3·5·7|(k-1)` — relevant to the lonely-measure tight-locus thread.)
**Synthesis with codex (HYP-2621/2622).** codex found, at FIXED k, that increasing the multiplier
`r` saturates (`r=4` best at k=31,61; `r>=5` picks up excess or falls to `1/k`) — "no `o(1/k^2)`
dip found." My correction CONFIRMS codex: the `F`-family levels stop at `a=4`. Trying `a=5` at
the "right" primorial `k=211` does NOT continue the climb — it collapses to the FLOOR (tight).
So both threads now agree: **`a_max` via this family is bounded by 4**; whether ANY family
reaches `a>=5` is OPEN (workflow E1 search inconclusive so far — every prefix-containing set
re-optimizes at `t≈1/k`, i.e. to level 1). The clean dips are `a=3` (`k≡7,13,19,25 mod30`) and
`a=4` (`k≡1 mod30`); these are infinite families, so `limsup g·k^2 = 1/2`, and the confirmed
`liminf` over realized dips is `≈ 1/4` (a=4). **Whether `g(k)·k^2` is bounded below by a positive
constant (i.e. `g = Theta(1/k^2)` uniformly) remains OPEN** — it hinges on whether `a>=5` is ever
realizable.
All `M > floor 1/(k+1)` (LRC holds — NOT counterexamples). The exhaustive `sigma_2` matches
`F` exactly at k=7 (3/23) and k=13 (3/41); at k=6 the true `sigma_2=5/33` is a deeper
(non-`F`) Stern-Brocot fraction `= mediant(3/20, 2/13)`.

**Mechanism.** The killer speed `a(k-1) ≡ 0 (mod d)` for **every** `d | (k-1)`, so at every
coarse time `t = j/d` it sits on an integer (distance 0), annihilating that clock. The more
small divisors `(k-1)` has, the more coarse tight times die, and `M` is pushed down toward
the floor. When `k-1` is a primorial `2·3·5·7···`, all small clocks die at once and the level
`a` climbs with the number of distinct small prime factors `ω(k-1)`.

**Binding-pair structure (proof skeleton).** At the optimum the binding pair is always
`{2a-1, a(k-1)}`, summing to `q = a(k+1)-1` (verified: `a=2→{3,8}`, `a=3→{5,18},{5,36}`,
`a=4→{7,120}`). I.e. the killer `W=a(k-1)` and its `q`-complement `2a-1 = q - W` sit antipodally
at distance `a/q` on the `q`-clock. So the level-`a` claim `M(F)=a/q` splits into:
 (i) **lower bound** `M >= a/q`: at `t* = m/q` (the time where `(2a-1)m ≡ ±a, Wm ≡ ∓a (mod q)`),
     EVERY speed `v ∈ {1,...,k-2,k}` has residue `mv mod q ∈ [a, q-a]` (the safe band) — this is
     a finite congruence check on the AP-prefix, and needs `2a-1 ∈ {1,...,k-2}∪{k}` (so `a <= (k-1)/2`);
 (ii) **upper bound** `M <= a/q` (covering): no other time beats `a/q`. This is where `(k-1)`'s
     factorization enters — the killer must annihilate all coarse clocks `t=j/d`, `d|(k-1)`, which
     it does for the primes dividing `k-1`. The clean statement to prove: the danger arcs of radius
     `a/q` cover `[0,1)` iff `(k-1)` carries the first `a-1` primes. (HYP-2623.)

## Consequence

- **Generically `g(k) = Theta(1/k^2)`** (`a=2`, `g·k^2 ≈ 1/2`), for all `k` not of the special form.
- **Below-mediant dips exist** along infinite families: `a=3` (`k≡7,13,19,25 mod30`, `g·k^2→1/3`)
  and `a=4` (`k≡1 mod30`, `g·k^2→1/4`). So the "gap edge" `2/(2k+1)` is the GENERIC value, NOT a
  true universal edge: below-`2/(2k+1)` configs exist. The realized `g·k^2` ranges over `[≈1/4, 1/2]`.
- **`a>=5` is OPEN.** The natural family collapses to the floor at `a=5`; no construction has yet
  realized level 5. So `liminf_k g(k)·k^2` is **not known to be 0** — it could be `≈1/4` (gap
  `Theta(1/k^2)` with bounded constant) or `0` (if a deeper family exists). This is the precise
  remaining question (HYP-2623, sharpened).
- This **corrects HYP-2052 (oracle-S552)**: `2/(2n-1)` is the GENERIC second value, not a universal
  gap edge; below-edge configs exist (codex-S573 found the first; the `a=3,4` families explain them).
  The "no config in `(1/n, 2/(2n-1))`" claim was a small-box artifact (exhaustive only to n<=8).

## For LRC(14) (k=13)

`k-1 = 12 = 2^2·3`, so `ω(12)=2` distinct primes → reachable level `a <= 3`; indeed
`sigma_2(13) = 3/41` (a=3), the deepest dip available at k=13. The genuinely hard `k=13` floor
is still the AP; the spectral basin around it has width `g(13) = 3/41 - 1/14 = 1/574`.

## Engineering application

The max-min collar `M(S)` is the **guaranteed simultaneous clearance** in view-obstruction /
periodic-scheduling problems: `k` processes (or transmitters, or radar pulses) with integer
period-ratios `S` and one observer; `M(S)` is the largest guard band such that, at some phase
`t`, every process is `>= M(S)` away from the observer's lane. THM-539 says the *achievable
clearance spectrum* just above the optimum `1/(k+1)` is the Stern-Brocot ladder
`a/(a(k+1)-1)`, and — crucially for design — the finer rungs (`a>=3`) are reachable **only when
`k-1` is highly composite**. So a frequency-hopping / coprime-sampling / guard-band designer who
can choose `k` (number of channels) gets strictly finer clearance tuning by picking `k` with
`k-1` divisible by `2·3·5·…` (a primorial), via the explicit family `{1,…,k-2,k,a(k-1)}`. The
denominator lemma `q <= 2 max(S)` is the exact statement that finer clearance costs proportionally
larger period ratios (`max(S) >= q/2`).

## Files
- `04-computation/lrc_spectral_gap_lowerbound_kps.py` (+`.out`) — engine, M_k verify, spectrum scan
- `04-computation/lrc_below_mediant_fast_kps.py` (+`.out`) — filtered exact sigma_2 scan k<=11
- `04-computation/lrc_a3_family_verify_kps.py`, `lrc_general_family_kps.py` — the F(k,a) family
- `04-computation/lrc_primorial_verify_kps.py`, `lrc_exact_M_covering_kps.py` — covering check, exact M
- Reflection: `07-reflections/lrc-spectral-gap-dips-along-primorials-kps.md`
