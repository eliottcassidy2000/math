---
source: opus-2026-06-02-S569 (remote-control)
status: REDUCTION (stride) — the paired/anchored dichotomy splits LRC(n) so the HARD worry-set is the EASY route; one sharp open lemma remains
tags: [LRC, private-obstruction, paired, anchored, shield, dichotomy, proof-split, THM-396, THM-397, n14, even-n]
---

# Every private obstruction is paired or anchored — and it splits the proof

**Prompt (user):** see if every private obstruction is paired or anchored; that
will help the proof.

It does. The right notion of "private obstruction" is a runner that **blocks a
small pair's pinch**, and the dichotomy (THM-396/397) is exactly **paired ∨
anchored**. Used correctly it flips the proof so the *hard* worry-set is the
*easy* route. Verified across even n=4..14.

## 1. The dichotomy, correctly stated (paired ∨ anchored) — verified

Fix a pair `(a,b)`, `D=a+b`, reduced sum `s=D/gcd(a,b) ≤ n`. Its **pair-safe
pinch times** are `t=m/D` (`m≢0 mod s`). A runner `c` **blocks** a pinch if
`‖cm/D‖<1/n`. Then (THM-396/397):

> a small pair `(a,b)` is **BLOCKED** (every pair-safe pinch killed) **only via**
> a **PAIRED** blocker — a *sum-multiple shield* `D | c` — **or** an **ANCHORED**
> blocker — an *endpoint* `‖c/D‖ < 1/n` (`c mod D` in the `±⌊(D-1)/n⌋` window).
> For `D ≤ n` the anchor window is empty, so the blocker must be a shield.

Verified: over 3428 blocked small pairs, **0** had neither a shield nor an anchor.
So **a small pair is UNBLOCKED ⟺ no runner is a shield and none is an anchor**, and
then its pinch is a witness.

## 2. The split this produces (the stride)

Classify every speed set by safe-set measure:

- **positive measure ⟹ `M>1/n` ⟹ lonely.** (trivial — the IGNORE regime.)
- **measure zero (the WORRY set, where a counterexample would live).**

The dichotomy lets us decide the worry set with the **small-pair** test, and the
data is decisive (exhaustive small boxes, `lrc_private_obstruction_split_s569.py`):

| n | WORRY (measure-0) | with an UNBLOCKED small pair | block-all configs | of those positive-measure |
|---|---|---|---|---|
| 4 | 1 | **1/1** | 65 | 65/65 |
| 6 | 2 | **2/2** | 50 | 50/50 |
| 8 | 3 | **3/3** | 19 | 19/19 |
| 10 | 1 | **1/1** | 4 | 4/4 |
| 12 | 1 | **1/1** | 1 | 1/1 |
| 14 | 1 | **1/1** | 2 | 2/2 |

**Two clean facts, every even n, no exceptions:**
1. **Every measure-zero (worry) config has an UNBLOCKED small pair** — so it has a
   witness pinch `t=m/D` with all runners `≥1/n`. The *hard* set is handled by the
   *cheap* route.
2. **Every config that blocks all its small pairs is positive-measure** (lonely).
   The block-all "residual" is entirely the *easy* spread set.

## 3. The reduction

> **LRC(n) ⟸ every measure-zero speed set has an unblocked small-reduced-sum pair.**

Proof of the implication: a positive-measure set is lonely; a measure-zero set, by
the hypothesis, has an unblocked pinch `t=m/D` where every runner is `≥1/n` — a
witness — so `M=1/n≥1/n`. Hence no set has `M<1/n`: **no counterexample.** ∎

And by the dichotomy, "unblocked small pair" = a pair `(a,b)`, `s=(a+b)/gcd≤n`,
with **no runner divisible by `a+b`** and (if `a+b>n`) **no runner in its endpoint
window**. So the reduction is purely arithmetic.

## 4. Why this is the right split (the worry-set is paired, the residual anchored-free)

The worry-set is the regular rotational encirclement (S566) / perfect transversal
(S568) / polygon face (HYP-2091). Its tight extremiser is the AP-like config whose
**straddling pair `(a, n-a)` sums to `n`** (S557 N2). That pair has `D=n`, so its
anchor window is empty, and the config has **no multiple of `n`** (C′, S564) — so
the pair is **unblocked**. *The paired/anchored dichotomy is exactly why the
worry-set's own straddle pair cannot be blocked:* to block a sum-`n` pair you would
need a multiple of `n` (paired shield), which the tight config lacks; anchors don't
exist at `D=n`. This is the hint made precise: **the worry-set's private
obstruction is its straddle pair, which is neither paired-away (no shield) nor
anchorable (D=n) — hence open.**

## 5. The one remaining lemma (sharp, open)

> **LEMMA (open, verified n=4..14): every measure-zero speed set has an unblocked
> small-reduced-sum pair.** Equivalently: *if every small pair is shielded-or-
> anchored, the safe set has positive measure.*

This is now a finite-flavored, arithmetic statement about the shield/anchor cover —
far sharper than the original conjecture, and the swarm's THM-396/397 + the pinch
lemma (S557) supply the machinery. Two sub-cases:
- **no small pair at all** (all pairwise reduced-sums `>n`, the spread/pairwise-
  coprime configs): prove positive measure (spread ⟹ low resonance ⟹ `M` large,
  S563).
- **small pairs all shielded/anchored:** a covering-design constraint — `n-1`
  runners providing a shield (`D|c`) or anchor for every small pair-sum `D`. Show
  this forces the speeds "spread" enough for positive measure, or is outright
  impossible for a measure-zero config.

## 6. Honest status

The dichotomy (paired ∨ anchored) is THM-396/397 (proved) and re-verified (0/3428).
The **reduction** (LRC ⟸ measure-0 has an unblocked small pair) is rigorous given
the dichotomy. The **split** — worry-set always route-1, residual always positive-
measure — is verified exhaustively for even n=4..14 with no exceptions. The
**remaining lemma** is open but sharp and arithmetic. This is a genuine
simplification: it concentrates LRC(n) into "shielding/anchoring every small pair
forces positive measure," handing the hard worry-set to the cheap unblocked-pinch
route via exactly the paired/anchored structure the hint named.

**Artifacts:** `04-computation/lrc_private_obstruction_split_s569.py`,
`05-knowledge/results/lrc_private_obstruction_{dichotomy,correct}_s569.out`,
`lrc_proof_split_even_n_s569.out`. Builds on THM-396/397 (shields/anchors), S557
(pinch), S564 (C′/no-multiple-of-n), S568 (worry-set), HYP-2091 (parity). New:
**HYP-2095**.
