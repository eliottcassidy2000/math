---
source: opus-2026-06-03-S599i/j/k/l (remote-control)
status: STRUCTURAL PROGRESS on C2b (the whole remaining content of LRC(14)) — five angles. VERIFIED: C2b ⟺ (multiple of n ⟹ p_0>0); the sharp floor M_min=2/(2n−1); the discrete-(2n−1)-tick mechanism (jS avoids {0,±1} ⟹ M≥2/(2n−1)); a two-regime dichotomy (critical-tick vs coarse); minimizers carry a shell-partner pair containing the multiple. NOT a proof.
tags: [LRC, C2b, C-prime, n14, discrete-witness, 2n-1, shell-collision, positive-measure, p0, farey, two-regime, multi-angle]
---

# C2b multi-angle assault: localising the last obstruction of LRC(14)

**Prompt (user):** work hard and creative on C2b this long session; try many angles.

**C2b** = the residual of C': configs `S={v_1,…,v_{n−1}}` with some `v_j ≡ 0 (mod n)`, no
dominance, prove **loose** (`M(S) > 1/n`, sharply `M ≥ 2/(2n−1)`). This is the *entire*
remaining content of LRC(14). Five angles; each result labelled.

## Angle 2 — the clean reformulation (verified)

> **C'(n) ⟺ TIGHT (`M=1/n`) ⟹ no speed `≡0 (mod n)`.** Via THM-406 (`tight ⟺ p_0=0`) this is
> **multiple of `n` ⟹ `p_0 > 0`** (the lonely set `{min‖v_i t‖≥1/n}` has *positive measure*).

**Verified** (`…s599j.py`): over the window, **zero** tight configs contain a multiple of `n`
(`n=3..7`; tight counts `1,1,2,2,1`). So C2b is exactly *"a multiple of `n` forces a
positive-measure laminar channel."* This is the cleanest statement of the target: not "a lonely
point exists" but "a lonely **interval** exists", separating C2b from the measure-zero worry-set.

## Angle 1 — the sharp floor and the discrete `2n−1` tick (verified)

The minimum of `M` over C2b configs is **exactly `2/(2n−1)`** (`2/5,2/7,…,2/15` for `n=3..8`),
the Farey neighbour above `1/n` (HYP-2163). Mechanism:

> **Sufficient witness:** if some dilate `jS` avoids `{0,±1} (mod 2n−1)` (i.e. `v_i j mod (2n−1) ∈
> {2,…,2n−3}` for all `i`), then `t=j/(2n−1)` gives `M(S) ≥ 2/(2n−1)`.

The forbidden `j` for runner `v_i` is `{±v_i^{-1}}` (plus zero-divisor `j` if `gcd(v_i,2n−1)>1`),
so `≤ 2(n−1)` forbidden values among `2n−2` nonzero `j` — a **tight** budget. A good `j` exists
iff the forbidden set does **not** cover, which needs *collisions*.

## Angle 1b — shell collisions free the tick (verified, partial)

Two runners share a forbidden `j` iff `v_i^{-1}=±v_k^{-1}`, i.e. `v_i ≡ −v_k (mod 2n−1)` —
**shell partners** (`v_i+v_k ≡ 0`, THM-401). Each such collision frees one `j`.

**Verified** (`…s599k.py`): "discrete tick works" correlates with "has a shell pair" but is **not**
equivalent — `works∧no-shell` and `fails∧shell` both occur (zero-divisor speeds and triple
coincidences also free `j`). So shell-collision is a *mechanism, not the criterion*. The
discrete tick **fails** on many configs — but precisely the **coarse** ones (next angle).

## Angle 3 — the two-regime dichotomy (verified)

Every C2b config falls in one of two regimes:
- **Critical:** the would-be-tight configs (near the floor). The `2n−1` discrete tick fires
  (`M ≥ 2/(2n−1)`); these are the minimizers.
- **Coarse:** a small-denominator witness fires. **Verified:** the discrete-tick *failures* all
  have `M ∈ {1/2,1/3,1/4,…}` — **strictly above `2/(2n−1)`** (e.g. `n=7`: failures start at
  `1/6 > 2/13`). They are loose by *more*, via a coarse clock `t=a/q`, `q` small.

> **C2b dichotomy:** a multiples-of-`n` config either (i) admits the `2n−1` tick (critical,
> `M=2/(2n−1)`), or (ii) admits a coarse clock `t=a/q` with `q` small (`M=1/q > 2/(2n−1)`).
> In both, `M ≥ 2/(2n−1)`. **No C2b config was found below the floor.**

## Angle 4 — minimiser structure (verified, refined)

The `M=2/(2n−1)` minimizers contain a **shell-partner pair mod `2n−1` in which the multiple of
`n` participates**. The *canonical* minimizer is `{n,n−1}` (since `n+(n−1)=2n−1≡0`; `n≡2^{-1}`,
`n−1≡−2^{-1}`, so their forbidden `j`'s `±n^{-1}=±2` and `±(n−1)^{-1}=∓2` **coincide**). Literal
`{n,n−1}`-containment holds for `n=3,4,5,7`; **`n=6` shows `(1,3,4,5,18)`** with `18≡7`, `4+7=11`
— the multiple `18` is the shell-partner of `4`, not the literal `n`. So the invariant is the
**shell-pair-carrying-the-multiple**, up to the scaling/`⟨2,−1⟩` symmetry (THM-407).

## The synthesis — where C2b now stands

> **C2b is reduced to a finite avoidance problem at modulus `2n−1`:** *for every config with a
> multiple of `n`, either some dilate `jS` avoids `{0,±1} (mod 2n−1)`, or `S` lies in a coset with
> a small-`q` clock.* Equivalently (Angle 2): **a multiple of `n` forces a positive-measure lonely
> interval.** The bound it yields is the sharp `M ≥ 2/(2n−1)`.

For **n=14** (`2n−1=27=3³`): the critical regime is the 3-shell residual (THM-407, gcd `1,3,9`);
the discrete tick must fire on each, with the multiple of 14 (`≡7·2^{-1}=7·14… mod 27`) supplying
the shell collision. The coarse regime is the all-odd / coset configs (`t=1/2`, etc.).

## What would close it (the open gap)

The dichotomy is **verified, not proven**. To finish C2b one must prove **either** branch always
holds:
1. *(tick)* the forbidden inverse-set `{±v_i^{-1}}∪Z` never covers `ℤ/(2n−1)\{0}` when `S` has a
   multiple of `n` **and** no small-`q` coset — i.e. the multiple always forces a collision; **or**
2. *(measure)* directly, `p_0>0` whenever a multiple of `n` is present — an explicit lonely
   interval around a non-clock point.
The Garsia–Milne involution (T2, HYP-2160) whose fixed points are the 3 n=14 shells is the route
to (1); the laminar-channel construction (S599g) is the route to (2).

## Honest status

- **Verified:** Angle 2 (tight ⟹ no multiple, `n=3..7`); the `2/(2n−1)` floor; the discrete-tick
  sufficiency; the two-regime dichotomy (no config below floor); the minimizer shell-structure.
- **Established (rigorous):** the sufficient-witness lemma (`jS` avoids `{0,±1}` ⟹ `M≥2/(2n−1)`);
  the `p_0>0` reformulation (via THM-406).
- **Open:** proving the dichotomy is exhaustive (no C2b config evades both regimes) — i.e. C2b
  itself. This is structural progress (a finite `2n−1`-avoidance problem + a positive-measure
  reformulation), **not** a proof of C2b/LRC(14).

**Artifacts:** `04-computation/lrc_c2b_angle1_discrete_witness_s599i.py`,
`…angle2_tight_nomultiple_s599j.py`, `…angle1b_shell_collisions_s599k.py`,
`…minimizer_struct_s599l.py` (+`.out`s). Builds on THM-398/401/406/407, HYP-2160/2163. New:
**HYP-2168**.
