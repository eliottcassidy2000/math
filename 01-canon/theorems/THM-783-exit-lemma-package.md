---
id: THM-783
title: The corrected r=8 exit package — anchored φ-extension, visitor balance, conditional factor-two de-phasing, and a conditional metric-extent theorem; raw wall count and unconditional de-phasing are withdrawn
status: PROVED at the stated hypotheses (anchored simple-wall φ-extension; period-sum; single-visitor break; cluster balance; factor-two drift under one-step paired indices; no-companion extent bound) + VERIFIED (bounded-bank maximum 6 only; laws 40/40, 45-period, 6/6 batteries) + CORRECTED (unconditional de-phasing and every absolute raw-wall K0 are withdrawn) + OPEN (varying-index balanced co-landings and a universal metric/core-sensitive exit)
source: opus-2026-07-14-S303 package, scope-corrected and strengthened by codex-2026-07-14-S10
renumber_note: claimed as THM-782 before the author observed codex-S9's earlier-pushed
  THM-782 phase-cell theorem; renumbered to THM-783 under the first-pusher protocol.
depends_on:
  - THM-779   # the token-walk criterion this analyzes
  - THM-773   # the token algebra
related: [THM-767, THM-771, THM-784, THM-786, THM-788, HYP-6840, HYP-6845, MISTAKE-147, MISTAKE-148]
verification: 04-computation/lrc14_exit_lemma_decision_opus_S303.py
  (+ 05-knowledge/results/lrc14_exit_lemma_decision_opus_S303.out)
---

# THM-783 — the corrected r=8 exit package

> **Superseding correction to the historical S304 banner:** the ratio-boxed
> census diagnosis was right, but the claimed de-phase/serving law and general
> extent completion were not. MISTAKE-148 withdraws them. The anchored
> recurrence, period-sum/single-visitor laws, cluster balance, and the
> no-companion conditional extent theorem are the surviving proved package.

> **Correction (codex-S10, THM-784/MISTAKE-147).** Sections (1)--(3) retain
> the corrected local content below; the unconditional claim in (4) is
> withdrawn but a factor-two lemma survives under explicit one-step pairing;
> (5) is corrected
> below without changing its metric conclusion. The finite census in (6) never
> established a uniform bound, and the absolute wall-count conjecture in the
> original (7) is false. A fixed slow rainbow chamber supports arbitrarily many
> fastest-owner walls. Raw wall count and metric extent are not equivalent.

**Frame.** THM-779's setting: lens 7, r = 8 owners W (7 ∤ w, distinct), wall of
owner o at x = (m + ½)/w_o, blocking run = a maximal streak of walls all
satisfying the wall condition (the seven non-walling tokens exactly F_7). Write
φ_i = w_{o_i}^{-1}·m_i (mod 7) for the i-th wall of the run.

## (1) The anchored φ-extension lemma (PROVED)

> Let wall `i` be a simple covered wall, and let the next global event coordinate
> be a simple wall `i+1`. Then wall `i+1` is covered **iff
> `φ_{i+1} ≡ φ_i + w_{o_i}^{-1} (mod 7)`.** Same-owner steps satisfy the
> congruence identically; every owner switch is one determined mod-7 equation.

*Proof.* After wall i the non-o_i tokens are F_7 and o_i's new token
−w_{o_i}^{-1}(m_i+1) duplicates its unique holder; no token changes before wall
i+1; wall i+1's condition (its non-walling tokens = F_7) holds iff o_{i+1}'s
token −w_{o_{i+1}}^{-1}·m_{i+1} equals that duplicated value. Rearranged, this
is the recurrence; for o_{i+1} = o_i it is m ↦ m+1, an identity. ∎

Both hypotheses matter. Without a covered anchor the congruence can hold while
both walls fail (for example `W={1,2,3,4,5,6,8,9}`, from the owner-9 wall at
`1/6` to the owner-8 wall at `3/16`). Without a simple next coordinate, a
simultaneous wall is uncovered even if one selected owner's congruence holds
(`W={1,2,3,4,5,8,10,30}`, from the covered owner-30 wall at `19/60` to the
owner-10/30 double wall at `7/20`). Calling both candidate walls “walls of the
run” in advance would make the original iff circular.

## (2) The period-sum law (PROVED; battery 40/40)

> If owner b walls twice inside a run (consecutive b-walls, m ↦ m+1),
> telescoping (1) gives **Σ_{a≠b} n_a·w_a^{-1} ≡ 0 (mod 7)**, where n_a is a's
> wall count strictly inside b's period — a Beatty count taking two adjacent
> values. **Flip-break:** if exactly one owner's count flips between
> consecutive periods, the law fails (its w^{-1} ≢ 0) and the run ends.

## (3) THE SINGLE-VISITOR BREAK (PROVED — the unconditional tooth; 45-period battery, 0 occurrences in-run)

Let f be the FASTEST owner. Every other owner's mesh exceeds 1/w_f, so each
complete f-period contains **at most one wall of each other owner**: visitor
counts n_a ∈ {0,1}. The period-sum law then reads Σ_{visitors} w_a^{-1} ≡ 0:

> **An f-period wholly inside a blocking run can never have exactly one
> visitor** (w_a^{-1} ≡ 0 mod 7 is impossible). Unconditional, all tuples.
> **Cluster balance:** any visitor set S must satisfy Σ_{a∈S} w_a^{-1} ≡ 0
> (mod 7); for pairs, **w_c + w_{c'} ≡ 0 (mod 7)** (battery 6/6).

## (4) CORRECTION — the unconditional de-phase lemma is withdrawn

The original argument assumed that paired `c`- and `c'`-wall indices both
advance by one and keep a fixed relative order. Neither follows from
co-visitation. For `f=8`, `c=2`, `c'=5`, every owner-2 wall co-visits an
f-period with an owner-5 wall forever (pairs near `.25/.30`, then `.75/.70`,
and their integer translates), although the displayed bound would be `17/12`.
The owner-5 index advances alternately by two and three and the relative order
flips.

Under the **extra** hypotheses of one-step paired indices and fixed order, the
elementary drift estimate is valid. Without fixed order the available window
has width `2/w_f`, giving at best a factor-two version. No blocking-specific
principle forcing those hypotheses was proved here. Therefore the unconditional
bound and the asserted mandatory `>=3`-cluster handover are withdrawn; a sound
varying-index/Beatty formulation is an open part of the cascade.

The factor-two statement under the weaker explicit **one-step paired-index**
hypothesis is nonetheless rigorous and does not require fixed order. Suppose
two owners `c!=c'` make `L` paired co-visits, each paired index advances by one,
and both walls of each pair lie in one `f`-period. Their signed separation
drifts each time by

```text
|1/w_c-1/w_(c')|=Delta/(w_c w_(c')),
Delta=|w_c-w_(c')|.
```

Every signed separation lies in `(-1/w_f,1/w_f)`. Hence the `L-1` equal drifts
fit in a window of width `2/w_f`, proving

```text
L < 2 w_c w_(c')/(w_f Delta)+1.                         (D)
```

The factor two is necessary: `(w_f,w_c,w_(c'))=(11,6,8)` has four successive
one-step paired co-visits with signed separations

```text
1/16, 1/48, -1/48, -1/16.
```

The old bound is `35/11<4`, while (D) gives `59/11>4`. The pair is balanced
because `6+8=0 (mod 7)`. This lemma controls a fixed one-step strand only; it
does not repair the varying-index counterexample `(8,2,5)` or prove a mandatory
handover.

THM-786 obtains the same factor-two scale without a one-step assumption in the
special but important case where one fixed companion serves consecutive
second-fastest-owner walls: comparing the first and last companion walls gives
the bound directly. Alternation among several companions remains open.

## (5) The corrected conditional extent theorem (PROVED)

Here a g-wall has a **balanced companion cluster** when the other visitors in
its complete f-period have inverse-residue sum `-w_g^{-1}`; equivalently, the
full visitor set containing g has sum zero.

> **If no g-wall in a run co-lands with a balanced companion cluster** (g = the
> second-fastest owner), then all g-walls sit within 1/w_f of the run's ends,
> and the run's extent is
> **< 1/w_g + 2/w_f.** Consequently, under the same no-companion hypothesis for
> its putative covered streak, any closed core-safe component longer than that
> contains a wall where blocking fails — a full 1/14-witness moment.

*Proof.* A g-wall at least `1/w_f` from both ends is bracketed by two consecutive
f-walls that are still inside the run. The intervening complete f-period has a
visitor set containing g. By the period-sum law this set has inverse-sum zero;
since `w_g^{-1}` is nonzero, it contains a balanced companion to g, contrary to
the hypothesis. Thus the closed central interval contains no g-wall. It is
strictly shorter than g's wall mesh `1/w_g`, giving the displayed bound. ∎

The original sentence “every complete in-run f-period is visitor-free” was too
strong and is withdrawn. Exact counterexample:

`W={2,8,17,18,19,20,24,29}`, with `f=29`, `g=24`.

The consecutive covered walls at `11/58, 7/36, 7/34, 13/58` have owners
`29,18,17,29`. The complete f-period has the balanced visitors `{18,17}`
because `18^{-1}+17^{-1}=5+2=0 (mod 7)`, but it contains no g-wall. This does
not affect the proof above, which needs only complete periods containing g.

**Why the companion hypothesis cannot simply be dropped (the Davenport
observation):** the seven non-f inverses are seven nonzero residues of Z_7, and
the Davenport constant D(Z_7) = 7 means a zero-sum (balanced) subset ALWAYS
exists algebraically. The obstruction to long runs is therefore not algebraic
scarcity of balanced clusters but GEOMETRY: balanced owners must co-land in one
f-window at every interior g-wall. A correct varying-index analysis of repeated
co-landings is still missing. THM-784 shows that this geometry governs only the
visitor-rich mode: a visitor-free fast refinement can persist inside a slow
rainbow chamber without invoking the cascade at all.

On a strand whose paired indices really do advance one-by-one, (D) bounds the
life of a fixed pair. General strands can skip indices and reverse order, so
their handovers still require a varying-index/Beatty analysis. After persistent
stalks are factored out, the remaining owner-switch problem is exactly how long
this geometric-arithmetic conspiracy can persist.

## (6) The decision census (VERIFIED on its stated bounded banks only)

Anchored-run maxima over sampled windows of length `1/5`, heights to `10^4`:

| family | n | median | max |
|---|---|---|---|
| generic random (w ≤ 10^4) | 80 | 3 | 5 |
| synchronized {2}+odds, Σw^{-1} ≡ 0 (the S303 loophole candidate) | 80 | 2 | 4 |
| near-ratio packets n·base ± 1 | 80 | 2 | 4 |
| all w ≡ r (mod 7) | 60 | 2 | 6 |
| annealed peak (600 steps from the best) | — | — | **6** |

The synchronized packet — the one tested family whose period-sums cancel by
design — performs worse than the generic samples. The maximum **within these
banks** is 6 (THM-779's earlier height-500 census found 5). THM-784 shows why
this cannot be called a working universal constant: the banks omitted static
seven-owner rainbow chambers refined by one arbitrarily fast owner.

THM-779 supplies a stronger arithmetic adversary. A fixed seven-owner
permutation stalk on a divisor-complete core-safe interval, together with
`A_m=182m+1`, has `2m` consecutive covered walls. Every complete fast-owner
period there is visitor-free, so all the proved visitor laws coexist with the
diverging count. The finite generator contained neither this family nor the
independent THM-784 slow-rainbow refinement.

## (7) Status of the exit lemma after this package

- Proved at the displayed hypotheses: anchored simple-wall extension (1); no
  single-visitor periods and cluster balance (3); the factor-two drift bound on
  one-step paired strands (4); and the extent bound under the no-companion
  hypothesis (5).
- The r = 8 pierce is now UNCONDITIONAL for every run configuration in which
  the second-fastest owner's walls fail to attract balanced co-landing clusters
  — while repeated companion co-landings still need a correct varying-index
  analysis.
- **REFUTED by THM-784:** “every blocking run has at most 6 walls.” The claimed
  equivalence with `extent < C/w_g` was also false: changing only the fastest
  mesh makes wall count diverge while the blocked interval stays fixed.
- REMAINING (sharp, geometric): split off visitor-free fast refinements of a
  slow seven-owner rainbow chamber. For the visitor-rich complement, bound
  consecutive balanced co-landings using a sound Beatty/varying-index cascade;
  the factor-two lemma applies only to one-step strands and no mandatory-
  handover theorem is available. Globally,
  prove a metric/core-sensitive statement that forces a free time in every
  relevant core-safe component. A normalized extent bound or a theorem that
  static slow rainbows cannot contain the required core component would suffice.
- Everything here is Lean-friendly (mod-7 arithmetic + Beatty counting).
