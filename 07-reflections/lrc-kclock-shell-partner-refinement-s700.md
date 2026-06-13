---
source: opus-2026-06-07-S700 (attack the LRC — proof or major refinement/reframe)
status: MAJOR REFINEMENT (not a full proof of LRC(14)). Two NEW explicit-witness lemmas reduce LRC(n)
  to a small, structured residual: (A) the k-clock witness (∃k≤n with no v_i≡0 mod k ⟹ M≥1/k≥1/n —
  generalizes THM-369 to ALL k, larger bound for smaller k); (B) the shell-partner lemma (NEW, rigorous
  for speeds coprime to C=2n−1: a pair v_i+v_j≡0 mod C shares its ±-forbidden set, dropping |F\{0}| from
  ≤2(n−1) to ≤2n−4<2n−2 ⟹ ≥2 good ticks ⟹ M≥2/(2n−1)>1/n). Together A∪B cover all but the residual
  R(n) = {all clocks fail AND no shell-partner}, which is SMALL (|R|=18,50,62,233 at n=5..8) and VERIFIED
  loose — witnessed at PAIR-SUM resolutions t=m/(v_i+v_j). The reframe: a witness hierarchy by resolution
  (k-clock 1/k ⊃ shell-partner 2/(2n−1) ⊃ pair-sum). Honest caveat: B's coprime hypothesis excludes the
  n=14 prime-3 strata, so it does not close LRC(14).
tags: [lonely-runner, LRC, k-clock-witness, shell-partner, 2n-1, discrete-witness, forbidden-set,
  residual, pair-sum, witness-hierarchy, THM-369, THM-401, refinement, honest-status]
---

# Attacking the LRC: the k-clock witness, the shell-partner lemma, and the pair-sum residual

**Prompt (user):** "attack the LRC and try to come away with a proof or a major refinement and reframe."

No proof of LRC(14) — that remains open. What this session produced is a genuine **major refinement**:
two clean explicit-witness theorems that cover all but a small, structured residual, plus a reframe of
LRC as a *witness hierarchy graded by resolution*. Recorded as **THM-420** + **HYP-2283**.

## 1. The k-clock witness (Lemma A — generalizes THM-369)

THM-369 used `t=1/n`: if no runner is `≡0 mod n`, then `‖v_i/n‖≥1/n`, done. The trivial but
underused generalization: **any** `k∈{2,…,n}` works, and *smaller `k` gives a larger margin*.

> **If ∃ `k∈{2,…,n}` with no `v_i≡0 (mod k)`, then `M(S) ≥ 1/k ≥ 1/n`** (witness `t=1/k`).

So the only configs not handled by *some* clock are the **divisibility-complete** ones: every
`k∈{2,…,n}` divides at least one speed. This is a strong constraint — at `n=14` it forces the speed
set to hit a multiple of each of `2,3,5,7,11,13` (and `4,8,9`), pinning much of the structure.

## 2. The shell-partner lemma (Lemma B — NEW, rigorous)

The repo's pair-sum modulus is `C=2n−1` (THM-401); a **shell-partner** is a pair `v_i+v_j≡0 (mod C)`.
The discrete witness `t=m/C` is good (gives `M≥2/C`) iff `v_k m ∉{0,±1} (mod C)` for all `k`. For
speeds coprime to `C`, the forbidden set is `F={0}∪{±v_k^{-1}}`. The key cancellation:

> A shell-partner `v_j≡−v_i` ⟹ `v_j^{-1}≡−v_i^{-1}` ⟹ `{±v_j^{-1}}={±v_i^{-1}}` — the partner pair
> contributes **2 forbidden values, not 4**. So `|F\{0}|≤2(n−1)−2=2n−4`, hence `|F|≤2n−3<2n−1`, so
> **≥2 good ticks survive** ⟹ `M≥2/(2n−1)>1/n`.

This is the same `2/(2n−1)` floor (S599h) reached not by min-over-multiples but by an **overlap count
in the forbidden set** — a structural reason a shell-partner is *automatically* loose with margin.
Verified: all shell-partner configs at `n=5,6,7` (121, 572, 2470 of them) have a good tick and
`M≥2/(2n−1)`.

## 3. The refined reduction and the pair-sum residual

Combining: LRC(n) ⟺ every config in the **residual**
```
   R(n) = { every k∈{2,…,n} divides some v_i  ∧  no shell-partner (v_i+v_j≢0 mod 2n−1) }
```
is loose. **Verified `R(n)` is small and entirely loose** (`…s700c.py`):
```
   n :   5    6    7    8        (of 205, 786, 2996, 11432 gcd-1 configs)
  |R| :  18   50   62  233
  minM : 3/13 1/5 3/19 3/20  — all > 1/n, in fact ≥ 2/(2n−1)
```
And the residual's loose witnesses sit at **pair-sum resolutions** `t=m/(v_i+v_j)`: e.g.
`M(1,3,4,10)=3/13` at `n=5`, and `3+10=13` is a pair-sum. So `R(n)` — the divisibility-complete,
shell-free core — is exactly the regime where the *additive/signed* face (S699: pair-sums, the
`2n−1` shells, the 2-coloring) does the work the multiplicative clocks cannot.

## 4. The reframe: a witness hierarchy graded by resolution

> **LRC's witnesses form a hierarchy by denominator (resolution):**
> - **`k`-clock** `t=1/k`, margin `1/k` — Lemma A, kills the divisibility-incomplete configs;
> - **shell-partner / discrete** `t=m/(2n−1)`, margin `2/(2n−1)` — Lemma B, kills the shell-partnered;
> - **pair-sum** `t=m/(v_i+v_j)` — kills the residual `R(n)` (verified n≤8).
>
> Each finer level handles what the coarser one cannot. LRC(n) is precisely: *the
> divisibility-complete, shell-free configs are loose, and a pair-sum witness does it.*

This is the reframe — LRC as a descent through resolutions, with the hard core being the additive
(pair-sum) face, not the multiplicative (clock) one. It explains *why* `n=14` is hard: `27=3³` makes
the prime-3 strata divisibility-complete with gcd>1 speeds (so coprime Lemma B doesn't apply), pushing
the whole difficulty into the pair-sum residual.

## 5. Honest status

- **PROVED:** Lemma A (all `n`); Lemma B (speeds coprime to `2n−1`, all `n`).
- **VERIFIED:** every config outside `R(n)` is loose by A/B; `R(n)` is small and loose (n≤8 window),
  via pair-sum witnesses.
- **NOT proved:** that `R(n)` is always loose (no uniform pair-sum witness theorem yet) — this is the
  irreducible core. **Lemma B's coprime hypothesis excludes the `n=14` prime-3 strata**, so the result
  does NOT close LRC(14). The honest gain is the reduction + reframe, not a proof.
- **Next probe:** a *pair-sum witness lemma* for `R(n)` (does every divisibility-complete shell-free
  config have a good `t=m/(v_i+v_j)`?) — and an extension of Lemma B to the gcd>1 (prime-3) strata via
  the `C/d` sublattice count, which is what `n=14` actually needs.

**Artifacts:** `04-computation/lrc_shell_partner_lemma_s700.py`, `lrc_residual_no_shell_s700b.py`,
`lrc_kclock_shellpartner_residual_s700c.py` (+`.out`s). Theorem: **THM-420**. Builds on THM-369
(clock witness), THM-398 (`C'` reduction), THM-401 (`2n−1`/shells), S599h/i (the `2/(2n−1)` floor +
discrete witness), S699 (signed LRC / pair-sums). New: **HYP-2283**.
