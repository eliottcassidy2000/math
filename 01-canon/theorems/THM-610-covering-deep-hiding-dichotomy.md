---
id: THM-610
title: The covering deep-hiding dichotomy — a covering family hides only at denominators q*≥n+1 (elementary), and a family tight at 1/n hides at a denominator divisible by n; hence a TIGHT COVERING family must hide at q*≥2n (a multiple of n, ≥28 for n=14). Elementary complement to THM-566 (no uniform bounded witness), sharpening the covering-route/non-covering-route split.
status: PROVED (Lemmas 1,2 elementary; verified 0/213 covering families violate q*≥15, AP q*=14 & even block q*=28 for the tight cases). The uniform-looseness upgrade (covering ⟹ M≥1/n+c) it points to is CONJECTURAL (HYP-2566).
source: mac-mini-2026-07-03-S30
depends_on:
  - THM-523   # the q-witness reduction (non-covering ⟹ shallow witness at q≤n); this is the dual
related:
  - THM-566   # codex: covering sets have no UNIFORM bounded-denominator witness (deeper cousin; this is the elementary fixed lower bound)
  - HYP-4060  # kps-S36: PRIMITIVE covering-min=14/183, unique extremizer=deep well, gcd refinement, tight-locus=AP mechanism. Lemma 2 here RIGORIZES its "tight config lands on 14th-roots / principal branch" step.
  - HYP-4061  # opus-S58: CV floor gatekeeper on covering families (deep well passes R'=0.818; drop-7 bottoms the floor R'=0.315). The measure-floor axis; THM-610 is the denominator-structure axis.
  - HYP-2566  # residual covering-set hard core uniformly loose (the conjecture this sharpens the structure of)
  - HYP-2561  # inf L = 1/1260; tight locus = {AP, Goddyn–Wong}
  - HYP-4043  # mac-mini LRCDilation (WLOG gcd=1) — load-bearing here: it is what makes the imprimitive even block collapse to the non-covering AP
  - CASE-convergent-not-covering-min  # the covering-min trajectory (2/13,2/15,4/33,3/31,…), a separate open thread
results:
  - 04-computation/deep_hiding_lemma_macmini_20260703.py
  - 05-knowledge/results/deep_hiding_lemma_macmini_20260703.out
  - 04-computation/covering_min_anneal_macmini_20260703.py
  - 05-knowledge/results/covering_min_anneal_macmini_20260703.out
external: Lonely Runner Conjecture; classical covering/small-divisor reductions (Bohman–Holzman–Kleitman, Tao arXiv:1701.02048).
---

# THM-610 — the covering deep-hiding dichotomy

**Setup.** `n=14`, `S` a set of `n-1=13` distinct positive integers, `M(S)=max_{t}min_{v∈S}||vt||`
the LRC view (`||x||`= distance to nearest integer). `S` is **covering** if for every
`q∈{2,…,n}` some `v∈S` has `q|v`. The **hiding denominator** `q*(S)` is the (lowest-terms)
denominator of a maximizer `t*=a/q*` of `min_v||vt||`. LRC(n) = `M(S)≥1/n` for all primitive `S`.

## Lemma 1 (covering ⟹ deep hiding). *A covering `S` has `q*(S) ≥ n+1`.*

**Proof (elementary).** Take any `t=a/q` with `2≤q≤n` (any `a`). Since `S` covers `q`, some
`v∈S` has `q|v`, so `v·(a/q)=(v/q)·a∈ℤ`, hence `||v·a/q||=0` and `min_{v∈S}||v·a/q||=0`.
Thus the objective is `0` at every rational of lowest-terms denominator `≤n` (denominator `1`
is `t∈ℤ`, also `0`). As `M(S)>0`, no such point is a maximizer; every maximizer has
`q*≥n+1`. ∎

This is the exact **dual** of the q-witness lemma (THM-523): *non-covering ⟹ a shallow
hiding spot at some `q≤n` (whence `M≥1/q≥1/n`)*, while *covering ⟹ all hiding is deep
(`q*≥n+1`)*. It is the elementary, fixed-lower-bound shadow of THM-566 (codex: covering sets
admit no **uniform** bounded-denominator witness); Lemma 1 gives the clean unconditional floor
`q*≥n+1` with a one-line proof.

## Lemma 2 (tight ⟹ `n | q*`). *If `M(S)=1/n` attained at `t*=a/q*` (lowest terms) then `n | q*`.*

**Proof.** With `gcd(a,q*)=1`, `||v·a/q*|| = ||v·a||_{q*}/q*` where `||·||_{q*}` is the integer
circular distance on `ℤ/q*` (`||m||_{q*}=min(m\bmod q*, q*-(m\bmod q*))∈ℤ_{≥0}`). Tightness gives
`min_{v}||v·a||_{q*} = q*·M = q*/n`. The left side is a nonnegative **integer**, so `q*/n∈ℤ`,
i.e. `n | q*`. ∎

## Corollary (tight covering ⟹ `q*≥2n`, a multiple of `n`).
Lemma 1 gives `q*≥n+1`; Lemma 2 gives `n|q*`; the least multiple of `n` exceeding `n` is `2n`.
**For `n=14`: a covering family tight at `1/14` must hide at `q*∈{28,42,56,…}`.**

## Verification (`deep_hiding_lemma_macmini_20260703.py`)
- **Lemma 1:** over **213** random primitive covering 13-sets, `q*≥15` in **all** (min `q*` observed
  `=15=n+1`, the bound is attained); min-dist `=0` confirmed at every `t=a/q`, `q≤14`. **0 violations.**
- **Lemma 2 / Corollary:** the two tight `M=1/14` families sit exactly where predicted —
  the **AP `{1,…,13}`** (primitive, **non-covering**) at `q*=14`; the **even block `2·{1,…,13}`**
  (imprimitive, covering) at `q*=28=2n`. The **deep well `{1,…,12,182}`** is *not* tight
  (`M=14/183`, `q*=183`, `14∤183`) — consistent.

## Why this matters (the route it sharpens)
1. **Clean covering/non-covering split.** LRC(14) after WLOG-`gcd=1` (LRCDilation): *non-covering*
   families get `M≥1/n` from a shallow witness (THM-523, tight at the AP); *covering* families
   hide only at `q*≥15` (Lemma 1). The two buckets are separated by hiding depth.
2. **The tight locus avoids covering.** The `M=1/14` tight primitive families (AP, Goddyn–Wong;
   HYP-2561) are non-covering; the only covering tight family is the **imprimitive** even block,
   which under WLOG-`gcd=1` (HYP-4043) collapses to the non-covering AP. So *primitive covering ⟹
   `M>1/14`* (strict). The **uniform** upgrade `M≥1/14+c` is HYP-2566.
   - **Lemma 2 rigorizes kps HYP-4060's mechanism.** kps derives the tight locus by taking the
     *principal branch* (`v_i∝k_i`, runners on the nonzero 14th-roots `k/14`) of the exact-landing
     condition — a heuristic branch choice. Lemma 2 makes it unconditional: **every** tight family,
     on any branch, has `14|q*`, so its runners lie on the `q*`-th roots with min circular distance
     `q*/14` — i.e. on a `(q*/14)`-dilated copy of the 14th-roots. At `q*=14` this is exactly kps's
     principal config; `q*=28,42,…` are the higher branches, and there the even block (`q*=28`) is
     the imprimitive representative. So "tight ⟹ 14th-root config" is a theorem, not a branch choice.
3. **Quantified looseness (evidence, `covering_min_anneal`):** annealing over primitive covering
   families n=7..14 gives the covering-min ratio `M/(1/n) ∈ [1.06, 1.11]` at every n — bounded
   away from `1`, **not** shrinking. (Reproduces the known winners `2/13,2/15`; finds a new
   sub-construction family `3/31<11/111` at n=11; is unreliable at large n — it misses the known
   `4/33` at n=9 — so it does not settle the exact n=14 covering-min, only bounds the margin.)

## What is NOT claimed
- Not a proof of LRC(14). Lemma 1,2 are structural; the looseness bound `c>0` (HYP-2566) and the
  tight-locus finiteness (HYP-2561) remain the open core.
- The exact covering-min trajectory (`2/13,2/15,4/33,3/31,…`) is a separate open problem
  (CASE-convergent-not-covering-min); Lemma 1 only says every entry hides at `q*≥n+1`.
