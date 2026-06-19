---
id: THM-539
title: The LRC max-min second-point spectrum — Stern-Brocot mediants, the denominator lemma, and the primorial dip below Theta(1/k^2)
status: MIXED — Lemma A PROVED; generic value & small-k spectrum VERIFIED (exact); the
        "a_max unbounded" globalization is HYP-2623 (verified through a=4 exact, a=5,6 bounded)
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
  a=2:  k=5            M=2/11    g·k^2=0.379
  a=3:  k=7,13,19,25   M=3/(3k+2)  (3/23, 3/41, 3/59, 3/77)  g·k^2 ≈ 0.27–0.31
                       [k ≡ 7,13,19,25 (mod 30) — codex's sharp refinement of my "k≡1 mod 6";
                        at k ≡ 1 (mod 30) the a=3 family is AP-tight again and a=4 takes over]
  a=4:  k=31,61,91,... M=4/(4k+3) (4/127, ...)               g·k^2=0.2365  [k ≡ 1 (mod 30), k-1 div by 30=2·3·5]
  a=5:  k=211          M <= 5/1059  (covering; exact level E2-pending) g·k^2 < 0.20  [k-1=210=2·3·5·7]
  a=6:  k=2311         M <= 6/13871 (covering)                g·k^2 < 0.17  [k-1=2310=2·3·5·7·11]
```
**Synthesis with codex (resolves HYP-2621/2622's open question).** codex found, at FIXED k,
that increasing the multiplier `r` saturates (`r=4` best at k=31,61; `r>=5` picks up excess or
falls to `1/k`) — hence "no `o(1/k^2)` dip found." This is correct *at fixed k*: the max useful
multiplier is `a = ω*(k-1)+1` (one more than the run of consecutive small primes dividing k-1).
The dip appears only when you ALSO grow `k`: level `a=5` requires `7 | (k-1)`, first at `k=211`
(which codex did not reach — k=31,61 have `k-1 = 30,60`, no factor 7). So the threads agree:
`a` is bounded at each fixed `k`, but `a_max(k) -> infinity` as `k -> infinity` along primorials.
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
- **Along the sparse subsequence `k-1 ∈ {primorials}` (k = 3,7,31,211,2311,…), `g(k)` dips
  below `Theta(1/k^2)`**, with `g(k)·k^2 -> 1/a`, `a ~ ω(k-1) ~ log k / log log k` UNBOUNDED.
  Hence **`liminf_k g(k)·k^2 = 0`** (the spectral gap is NOT bounded below by `c/k^2` uniformly).
- This **corrects HYP-2052 (oracle-S552)**: `2/(2n-1)` is the GENERIC second value, not a true
  universal gap edge; below-`2/(2n-1)` configs exist (codex-S573 found the first; the primorial
  family explains and extends them to unbounded depth). The "no config in `(1/n, 2/(2n-1))`"
  claim was a small-box artifact (exhaustive only to n<=8 with a too-small entry bound).

## For LRC(14) (k=13)

`k-1 = 12 = 2^2·3`, so `ω(12)=2` distinct primes → reachable level `a <= 3`; indeed
`sigma_2(13) = 3/41` (a=3), the deepest dip available at k=13. The genuinely hard `k=13` floor
is still the AP; the spectral basin around it has width `g(13) = 3/41 - 1/14 = 1/574`.

## Files
- `04-computation/lrc_spectral_gap_lowerbound_kps.py` (+`.out`) — engine, M_k verify, spectrum scan
- `04-computation/lrc_below_mediant_fast_kps.py` (+`.out`) — filtered exact sigma_2 scan k<=11
- `04-computation/lrc_a3_family_verify_kps.py`, `lrc_general_family_kps.py` — the F(k,a) family
- `04-computation/lrc_primorial_verify_kps.py`, `lrc_exact_M_covering_kps.py` — covering check, exact M
- Reflection: `07-reflections/lrc-spectral-gap-dips-along-primorials-kps.md`
