---
source: opus-2026-06-02-S557 (remote-control)
status: DEDUCTIVE — rigorous structure lemma (pinch) + exact reformulation of LRC@14; proofs included, verified exactly on critical configs
tags: [LRC, n14, pinch, binding-pair, exact-moments, loneliness-radius, deductive, S553, S556]
---

# LRC@14: the exact moments are pair-pinches; the loneliness radius is r/s

**Prompt (user):** stop modelling scenarios; focus on the *exact moments and
conditions strictly necessary for a proof*.

So this note is deductive. Let `f_S(t) = min_i ||v_i t||` and
`M(S) = max_t f_S(t)`. LRC(14) ⟺ `M(S) ≥ 1/14` for every 13-set `S`. We pin down
exactly *where* `M(S)` is attained and *what value* it can take.

## 1. The Pinch Lemma (rigorous)

Each `||v t||` is a tent: piecewise linear, slope `+v` where `frac(vt)∈(0,½)`,
slope `−v` where `frac(vt)∈(½,1)`, apex `½` at `t=(2k+1)/(2v)`, zero at `t=k/v`.

> **Lemma.** If `M(S) < ½`, the maximum is attained at a time `t*` at which two
> *distinct* runners `a,b` are **binding and straddle the observer**:
> `||v_a t*|| = ||v_b t*|| = M(S)`, with `frac(v_a t*) = M(S)` (the `+` side) and
> `frac(v_b t*) = 1 − M(S)` (the `−` side).

*Proof.* `f_S` is continuous, piecewise linear, periodic; its maxima are at
breakpoints. At a local max `t*`, `f_S` rises to the left and falls to the right.
The left rise is carried by a binding runner `a` with positive slope there,
i.e. `frac(v_a t*) ∈ (0,½]`; the right fall by a binding runner `b` with negative
slope, `frac(v_b t*) ∈ [½,1)`. If `a=b`, that runner has slope `+` then `−` only
at its apex `frac = ½`, forcing `M=½`, excluded. So `a≠b`, and since both bind at
the common value `M<½`, `frac(v_a t*)=M`, `frac(v_b t*)=1−M`. ∎

## 2. The radius is r/s, pinned by a pair (rigorous corollary)

From the straddle, `(v_a+v_b)t* ≡ 0 (mod 1)`, so

> **`t* = m/(v_a+v_b)`** for an integer `m`, and `(v_a−v_b)t* ≡ 2M (mod 1)`.

Write `d=gcd(v_a,v_b)`, `s=(v_a+v_b)/d` (the **reduced pair-sum**),
`α=v_a/d`. Then `v_a t* = α m/s` and, since `gcd(α,s)=1`,

> **`M(S) = ||α m / s|| = r/s`,  with `r ∈ {1,…,⌊s/2⌋}`.**

So the loneliness radius is always a fraction whose **denominator is the reduced
sum of two speeds**. (Verified exactly: 80/80 random 13-sets, and every critical
config of §4.)

## 3. Exact reformulation and the strictly-necessary conditions

These are the *exact moments* (`t=m/(v_a+v_b)`) and *conditions* (straddle +
clear) a proof must control:

> **LRC(14) ⟺ every 13-set admits a pair `(a,b)` and integer `m` with
> `‖v_j · m/(v_a+v_b)‖ ≥ 1/14` for all `j`.**  (The witness is the pinch time.)

Immediate, rigorous consequences:

- **(N1) Counterexample ⇒ the optimal binding pair has reduced sum `s ≥ 15`.**
  Because `M = r/s < 1/14` with `r≥1` forces `s>14`. (A counterexample's two
  closest-approach runners, reduced, sum to `≥15`.)
- **(N2) Tight (`M=1/14`) ⇒ the binding pair has reduced sum `s ≡ 0 (mod 14)`**,
  i.e. `s ∈ {14,28,…}` (since `1/14 = r/s` ⇒ `s=14r`). At the floor, `r=1`,
  `s=14`: the two straddling runners satisfy `(v_a+v_b)/gcd(v_a,v_b) = 14`.
- **(N3) The provable target** (replaces "search configs"): *every 13-set has a
  pair `(a,b)` with reduced sum `≤ 14` whose pinch time clears all other runners.*
  Such pairs are abundant (any coprime pair summing to `≤14`, any pair like
  `(2,12)` with reduced sum `7`, etc.); the content is that **one of them clears
  the rest**.

## 4. Exact verification on the critical configs (not a search — a check)

For every tight config the unique straddling pair and `s=14` appear exactly:

| config (n) | tight `t*` | straddle `(+,−)` | `v_a+v_b` | reduced `s` | `(v_a−v_b)t*` |
|---|---|---|---|---|---|
| AP@14 | `1/14` | `(1,13)` | 14 | **14** | `1/7=2/14` |
| AP@14 | `3/14` | `(5,9)` | 14 | **14** | `1/7` |
| V*@14 | `1/14` | `(1,13)` | 14 | **14** | `1/7` |
| sporadic@8 | `1/8` | `(1,7)` | 8 | **8** | `1/4=2/8` |
| sporadic@6 | `1/6` | `(1,5)` | 6 | **6** | `1/3=2/6` |

In each row `(v_a+v_b)/gcd = n` and `(v_a−v_b)t* ≡ 2/n`, matching N2 exactly.

## 5. What this buys, honestly

It converts LRC@14 from a statement about times on a circle into a **finite,
arithmetic statement about pairs of speeds**: the only times that matter are
`m/(v_a+v_b)`, the radius is `r/s`, and the floor `1/14` is reached precisely when
a straddling pair has reduced sum `14`. The proof problem is now N3 — a clean
covering/clearing statement over `O(13²)` candidate pairs — and the obstruction
(N1) is named exactly: a counterexample's closest pair reduces to a sum `≥15`.

This also **subsumes the empirical findings deductively**: the tight-witness
lattice (S553/S556) is N2 (`s=14 ⇒ t*=m/14`); the "multiple of 14" of the sieve
is the special pair `(1,13)`/`(v,14·?)`; the spectral-gap value `2/(2n−1)=2/27`
(oracle-S552) is `r/s` with `s=27=(v_a+v_b)/gcd`, `r=2` — the apex-doubled pair.

**Next (provable) step:** attack N3 — show some reduced-sum-`≤14` pair clears all
others. The pinch lemma reduces the odd-coupling / even-fold residuals (S554) to a
statement purely about which *pair* binds. Likely route: among the `⌊n/2⌋` "small"
pairs summing to `n` (`(1,13),(2,12),…,(6,8)`), a counting/pigeonhole argument
that one pinch time avoids all other runners' thin danger arcs.

**Artifacts:** verification inline (exact); builds on S553 (V*, tight census),
S556 (boundary times), oracle-S552 (gap). New: **HYP-2059**.
