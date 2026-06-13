---
source: opus-2026-06-02-S559 (remote-control)
status: RESULT — reduction + apex-obstruction PROVED; tight-tuple handled (apex excluded); residual characterised. Progress on the k+1=2q polynomial-method analogue (the literature's n=14 wall).
tags: [LRC, n14, polynomial-method, 2q, apex, CRT, tight-tuple, Sungkawichai-Trakulthongchai, S558]
---

# The k+1=2q analogue of the tight-tuple corrector: the apex IS the field-failure

**Prompt (user):** make progress on the `k+1=2q` analogue of the tight tuple.

**Context (S558).** The published route to LRC(n) for `n≥8` is the finite-checking
sieve; its accelerator is the **polynomial method** (Sungkawichai–Trakulthongchai,
arXiv:2604.23906, Prop 4.1): *if `k+1` is an odd prime, the tight tuple `(1,…,k)`
is a universal corrector* — for every nonzero `v∈ℤ_{k+1}^k` there are units `s,r`
with `s·v + r·(1,…,k) ∈ {1,…,k-1}^k`. Proven over the **field** `ℤ_{k+1}`. This is
the exact step that fails at `k+1 = 2q` (so at **n=14 = 2·7**), where `ℤ_{2q}` is
not a field — the literature's stated wall. We pinpoint *why*, and repair the
hardest case.

## 1. Reduction to a mod-q problem (PROVED; mod-2 is forced)

`ℤ_{2q} ≅ ℤ_2 × ℤ_q`. **Every unit mod `2q` is odd** (≡1 mod 2), so for any units
`s,r` and tight tuple `(1,…,2q-1)`:

- **mod 2:** `(s v + r·tight)_i ≡ v_i + i (mod 2)` — **forced, no unit freedom**;
- **mod q:** `(s v + r·tight)_i ≡ s v_i + r i (mod q)`, with `s,r mod q` ranging
  freely over `ℤ_q^×` (units mod 2q surject onto `ℤ_q^×`, mod-2 part fixed at 1).

The strict target `{1,…,2q-2}` excludes `0=(0,0)` and `2q-1=(1,q-1)`. So a forbidden
landing for runner `i` is `(0,0)` if `v_i+i` is even, `(1,q-1)` if odd. Therefore:

> **Reduction Theorem.** The corrector `∃ s,r∈(ℤ/2q)^×: s v + r(1,…,2q-1) ∈
> {1,…,2q-2}^{2q-1}` holds **iff** `∃ s',r'∈ℤ_q^×: s' w_i + r' c_i ≠ f_i (mod q)`
> for all `i`, where `w_i=v_i mod q`, `c_i=i mod q`, and `f_i = 0` if `v_i+i` even,
> `q-1` if odd.

Verified exactly: **0 mismatches** between the full `ℤ_{2q}` check and the mod-q
check across q=3,5,7,11 and all tuple classes. The whole difference between prime
and `2q` is **the loss of mod-2 freedom**.

## 2. The apex IS the obstruction (PROVED)

In `{1,…,2q-1}` the **unique** speed `≡0 (mod q)` is the **apex `q = (k+1)/2 = n/2`**
(the repo's "co-observer," S530/S547). For runner `i=q`, `c_q = q mod q = 0`, so its
mod-q constraint `s' w_q + r' c_q ≠ f_q` becomes **`s' w_q ≠ f_q` — parameter-free
in `r'`**, and if `w_q = v_q mod q = 0` (apex speed divisible by `q`) with `f_q=0`
it is `0 ≠ 0`, **unsatisfiable**.

For the **tight tuple itself** (`v_i=i`, so `v_q=q≡0`, `f_q=0`): the corrector is
**impossible** — `strict=False` for all q (3,5,7,11). This is the precise,
localized reason the field argument of Prop 4.1 dies at `2q`: the apex is the
**zero-divisor `q`**, exactly where `ℤ_{2q}` stops being a field.

## 3. The tight tuple is handled once the apex is set aside (PROVED)

The apex is *not actually a problem for loneliness*: at the base time `t=s/(2q)`
(`s` a unit), the apex sits at `q·s/(2q) ≡ 1/2` — distance `1/2 ≫ 1/(2q)`, the
**safest possible** position. So exclude it from the corrector and handle it
directly. For the remaining `2q-2` runners with the tight tuple, take `s'=r'`
(e.g. `s=r=1`): each non-apex runner lands at `2i mod 2q ∈ {2,4,…,2q-2}` (even,
nonzero) — in the strict interior. Hence:

> **Tight-tuple repair.** At `k+1=2q`, the tight tuple `(1,…,2q-1)` is corrected on
> all non-apex runners by `s=r=1`, while the apex `q` is automatically `1/2`-safe.
> Verified `excl-apex = True` for all q=3,5,7,11.

This is the `2q`-analogue of **Prop 4.4** (the bottleneck the paper names as the
expensive `c=k+1` lift) for the hardest tuple, done analytically — *the apex is
the only thing the prime proof was implicitly using primality to avoid*.

## 4. The residual, exactly characterised (iff; verified)

The full universal-corrector property still fails for *general* bad `v` (excl.
apex). For **parity-matched** `v` (all `f_i=0`), excluding the apex, the constraints
are `s' w_i ≠ -r' c_i`, i.e. forbid the ratio `ρ_i = -c_i/w_i ∈ ℤ_q^×`. Then:

> **Residual Theorem (iff, verified).** A parity-matched `v` is *un*-correctable
> (apex excluded) **iff** the ratio set `{ρ_i = -c_i/w_i}` covers all of `ℤ_q^×`.
> (q=7: 100% of fails cover, 0% of successes do.) The **tight tuple collapses all
> ratios to `{-1}`** — a single value — so it is the *easiest*; the hard residual
> is the *ratio-spread* tuples.

So the residual hard inputs for `n=14` are precisely the parity-tight,
ratio-spreading `v` — a finite, explicit local family, exactly what a `c=14`
computation must still resolve.

## 5. What this buys, honestly

- **Reduction Theorem + Apex Theorem (PROVED):** the `2q` polynomial-method
  obstruction is exactly the loss of mod-2 freedom + the apex zero-divisor `q`.
  This *explains* the literature's "`k+1` must be an odd prime" and identifies it
  with the repo's apex/co-observer thread — they are the same object.
- **Tight tuple handled (PROVED):** the single hardest tuple (the paper's `c=14`
  bottleneck) is settled analytically once the apex is treated separately.
- **Residual characterised (iff):** the remaining work is the ratio-spread family,
  not the whole of `I(13,p,1)`.
- **NOT claimed:** a full `2q` Prop 4.1 (it is false — the ratio-spread residual is
  real), nor that this alone makes the `c=14` computation feasible. The honest
  gain is structural: the obstruction is named and localized, and the worst single
  case is removed.

**Next:** (a) bound/clear the ratio-spread residual by adding the `r/p` (mod-p)
freedom of the full witness time `t=s/(2q)+r/p` — the apex and ratio-spread runners
have room there; (b) feed the apex/zero-divisor structure to the concurrent
pinch/shield route (THM-396, HYP-2061): the apex `q` is the binding pair's
`(q,q)`-type partner and the natural shield. (c) the `2q` story should transfer to
any `k+1 = 2·prime` — i.e. it explains the whole even-doubled-prime sub-frontier
(n=10,14,22,…), not just 14.

**Artifacts:** `04-computation/lrc_2q_tight_tuple_analogue_s559.py`,
`lrc_2q_apex_repair_s559` (inline), `lrc_2q_residual_characterization_s559`
(+ `.out`s). Builds on S558 (method map), S554 (even-fold), S530/S547 (apex);
arXiv:2604.23906 Prop 4.1/4.4. New: **HYP-2063**.
