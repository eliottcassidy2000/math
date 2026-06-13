---
source: opus-2026-06-03-S577 (remote-control)
status: SYNTHESIS — the additive-circuit lemmas (A: circuit-free randomness, B: 3-term fold) and the multiple/Φ framework share ONE residual; the fold IS the δ-clock; safe measure is depressed only by 3-term relations
tags: [LRC, additive-circuit, 3-term, fold, sieve, circuit-free, Lemma-A, Lemma-B, C-prime, Phi, delta-clock, n14, synthesis]
---

# Going back and forth: the fold is the δ-clock, and only 3-term relations depress the gap

**Prompt (user):** focus mostly on Lemma B; bounce to Lemma A when stuck; use the
particularities of one as noise for the other. **Lemma A** (randomness): circuit-free
⟹ `G ≥ 1/(k+1)` via measure/equidistribution. **Lemma B** (structure): a 3-term
relation `v_c=v_a+v_b` is a literal fold `t v_c = t v_a + t v_b` (fold+sieve); that it
delivers `1/(k+1)` is unproven (2025 prime-divisibility machinery plugs in here).

Four rounds of bouncing (k = #nonzero speeds, δ = 1/(k+1)) converged the additive
picture onto the multiple-of-(k+1) / `Φ` framework I had built. **The fold is the
δ-clock; the open delivery of δ is exactly C′.**

## Round 1 — fold = shield; circuit-free is far from the floor
- **The fold is a shield.** `v_c=v_a+v_b` ⟹ at the `(a,b)`-pinch `t=m/v_c`, `c` sits
  at integer `m` (on the observer). A 3-term relation *blocks* the `(a,b)` pinch (= a
  THM-396 shield). Verified.
- **Lemma A premise (de-risked):** circuit-free (no 3-term) `min M ≈ 0.21–0.25` —
  bounded below by a constant `≈1/5`, far above δ, margin **growing** with k
  (+0.05→+0.13, k=4..10). 4-term-rich circuit-free configs share that min: **4-term
  structure alone never makes a hard config.**

## Round 2 — hardness ⟺ 3-term richness
Binning `M−δ` by 3-term count (k=6): `#3t=0: +0.088`, `1: +0.045`, `2: +0.034`,
`≥3: +0.011`. **Hardness is monotone in 3-term richness.** And among near-tight
configs (`M−δ<0.02`): mean #3-term `= 4.7, 7.7, 12.5` (k=6,8,10), **circuit-free among
them = 0%**. The contrapositive of Lemma A holds empirically: *a hard config always
has a 3-term relation.*

## Round 3 — the convergence: the fold is the δ-clock
Near-tight configs split exactly by **multiple of (k+1)**:

| k (n=k+1) | near-tight | no multiple of n | δ-clock `j/n` witness |
|---|---|---|---|
| 6 (7) | 28 | 18 | **18** |
| 8 (9) | 10 | 4 | **4** |
| 12 (13) | 1 | 0 | **0** |

`#(no multiple) == #(δ-clock witness)` **exactly** — i.e. **no-multiple-of-(k+1) ⟺
the δ-clock `t=j/(k+1)` is a witness** (THM-369), reached from the additive side. The
3-term relations (folds) strip away every pair-pinch *except* the straddle pair
`(a, n−a)` summing to `n=k+1`; that pinch *is* the δ-clock. So **Lemma B's "fold+sieve
delivers δ" = the straddle pinch at `t=j/(k+1)`** — a witness iff no runner is a
multiple of `k+1`. The remaining `had-multiple` configs are loose off-clock (C′/Φ>0).

## Round 4 — only 3-term relations depress the measure
Safe measure `μ` at level δ:

| k | circuit-free min μ | 4-term-rich-but-cf min μ | 3-term-rich (≥3) min μ |
|---|---|---|---|
| 6 | 0.0806 | **0.0806** | 0.0037 |
| 8 | 0.1015 | **0.1015** | 0.0185 |

**The safe measure is depressed *only* by 3-term relations.** Circuit-free (even
4-term-rich) keeps `μ ≳ 0.08`; 3-term relations push `μ → 0` (the worry-set boundary).
So the relevant "energy" for hardness is the **3-term count**, not the additive energy
(4-term) — exactly the user's de-risking, now quantified on the measure.

## The 2×2 synthesis (two slicings, one residual)

|                | no multiple of (k+1) | multiple of (k+1) |
|----------------|----------------------|-------------------|
| **circuit-free** | δ-clock; `M ≥ ~0.2` (Lemma A) | `M ≥ ~0.2` (Lemma A) — loose |
| **3-term-rich** | δ-clock; **can be TIGHT** (worry-set) | **C′ / Φ>0** (the residual) |

- **Worry-set (M=δ) lives only in `[3-term-rich ∩ no-multiple]`.** Tightness *requires*
  folds (circuit-free ⟹ `M≥0.2`); the maximally-folded config (the AP, every
  `i = 1·(i−1)+1`) is the extremal tight one.
- **The open residual is the same under both slicings: `[multiple of (k+1)]` = C′ =
  `ker Φ` ∩ {n|v}.** The additive route (Lemma A+B) and the multiple route (THM-398)
  cross, and both bottom out at C′.

## What the cross-pollination delivered (the "noise")

- **additive → multiple.** The fold explains *why* the worry-set falls on the δ-clock:
  3-term relations remove every other pair-pinch, leaving only the straddle pair = the
  δ-clock. And tightness *needs* folding — circuit-free can't reach the floor.
- **multiple → additive.** Lemma B's unproven "fold delivers δ" **is** C′
  (multiple-of-(k+1) ⟹ loose), for which my ECCP machinery gives the *exact gap*
  `G(v)=Φ(C)` (HYP-2112) and proved peeling criteria (B′/C/E/F, S581) covering 100% of
  sampled n=14.
- **measure coupling.** `μ(safe) ≈ const − (3-term contribution)`: the two lemmas meet
  in one inequality — circuit-free keeps `μ` up; 3-term relations drive it to 0.
- **prime machinery.** For `k+1=p` prime, `p∤v_i j ⟺ p∤v_i` (any `i`), so no-multiple
  ⟹ `t=1/p` works trivially; the residual `p|v_i` needs a speed `≥ p` (large, often
  dominant ⟹ dodgeable, THM-398 Lemma B). For `k+1=14=2·7` composite the δ-clock has
  *extra* `j` witnesses but the `2·7` apex residual (HYP-2063) is the genuine C′ core.

## Honest status / skeleton

```
LRC(k)  ⟸  [ no multiple of (k+1)  ⟹  δ-clock witness ]      (THM-369, PROVEN; = the fold/straddle pinch)
         ∧  [ multiple of (k+1)     ⟹  loose, Φ(C)>0  ]       (C′, OPEN residual; both slicings end here)
   with  [ circuit-free ⟹ μ(safe) ≥ ~0.08 ⟹ M > δ ]          (Lemma A, VERIFIED, proof open: discrepancy from 3-term count)
```
- **Lemma A** (circuit-free ⟹ `M ≥ const > δ`): **verified** (min μ ≈ 0.08–0.10,
  depressed only by 3-term); proof open (a discrepancy bound keyed to the 3-term count).
- **Lemma B** (3-term fold ⟹ δ): **reduces to** δ-clock (proven, no-multiple) + C′
  (open, multiple) — *not* a separate open problem.
- **The residual is C′ = `ker Φ` ∩ {multiple of (k+1)}**, attacked by the ECCP
  translator / `Φ` / peeling lemmas. The additive picture confirms it from the other
  side and explains the fold→δ-clock mechanism.

**Artifacts:** `04-computation/lrc_additive_circuit_AB_round{1,2,2b,3,4}_s577.py`
(+`.out`). Builds on THM-398/HYP-2102 (C′), HYP-2112 (Φ), S569 (unblocked straddle),
THM-369 (δ-clock), HYP-2063 (2q apex). New: **HYP-2114**.
