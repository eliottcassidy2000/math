---
id: THM-782
title: The exit-lemma package — the φ-recurrence, the period-sum law, the SINGLE-VISITOR BREAK (unconditional), cluster balance, the de-phase bound, and the conditional extent theorem; the loophole packet eliminated; Davenport shows balance is algebraically unavoidable so the surviving gap is geometric
status: PROVED (φ-recurrence; period-sum; single-visitor break; cluster balance; de-phase lemma; the no-companion extent bound) + VERIFIED (run cap 6 at heights to 10^4 incl. all targeted packet families; laws 40/40, 45-period, 6/6 batteries) + OPEN (the absolute run constant; the balanced co-landing cascade, named)
source: opus-2026-07-14-S303 (owner directive: prove the unconditional exit lemma)
depends_on:
  - THM-779   # the token-walk criterion this analyzes
  - THM-773   # the token algebra
related: [THM-767, THM-771, HYP-6840, HYP-6845]
verification: 04-computation/lrc14_exit_lemma_decision_opus_S303.py
  (+ 05-knowledge/results/lrc14_exit_lemma_decision_opus_S303.out)
---

# THM-782 — the exit-lemma package

**Frame.** THM-779's setting: lens 7, r = 8 owners W (7 ∤ w, distinct), wall of
owner o at x = (m + ½)/w_o, blocking run = a maximal streak of walls all
satisfying the wall condition (the seven non-walling tokens exactly F_7). Write
φ_i = w_{o_i}^{-1}·m_i (mod 7) for the i-th wall of the run.

## (1) The φ-recurrence (PROVED)

> Blocking propagates from wall i to wall i+1 **iff φ_{i+1} ≡ φ_i + w_{o_i}^{-1}
> (mod 7).** Same-owner steps satisfy it identically; every owner switch is one
> determined mod-7 equation.

*Proof.* After wall i the non-o_i tokens are F_7 and o_i's new token
−w_{o_i}^{-1}(m_i+1) duplicates its unique holder; no token changes before wall
i+1; wall i+1's condition (its non-walling tokens = F_7) holds iff o_{i+1}'s
token −w_{o_{i+1}}^{-1}·m_{i+1} equals that duplicated value. Rearranged, this
is the recurrence; for o_{i+1} = o_i it is m ↦ m+1, an identity. ∎

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

## (4) The de-phase lemma (PROVED)

Two owners c ≠ c′ co-visiting f-periods drift apart at rate
|1/w_c − 1/w_{c′}| = Δ/(w_c w_{c′}) per c-period (Δ = |w_c − w_{c′}| ≥ 1, distinct
integers), while co-landing requires staying within one f-window 1/w_f:

> consecutive co-visits of a FIXED pair ≤ **w_c·w_{c′}/(w_f·Δ) + 1.**

Balanced clusters are temporary; sustaining a run past a fixed pair's de-phase
needs a handover through a ≥3-cluster with Σ w^{-1} ≡ 0 co-landing at the
handover moment — each a further determined coincidence, cascading with the
two-owner congruences (n ≡ w_f·w_b^{-1} count-locks; the (f,c,c′)
tri-consistency forces the f-count n ≡ 4 or 0 mod 7).

## (5) The conditional extent theorem (PROVED) and why it cannot be unconditional as stated

> **If no g-wall in a run co-lands with a balanced companion cluster** (g = the
> second-fastest owner), then every complete in-run f-period is visitor-free,
> all g-walls sit within 1/w_f of the run's ends, and the run's extent is
> **< 1/w_g + 2/w_f.** Consequently any closed core-safe component longer than
> that contains a wall where blocking fails — a full 1/14-witness moment.

*Proof.* A g-wall ≥ 1/w_f from both ends lies in a complete in-run f-period,
which then has a visitor; by (3) it needs ≥ 2 balanced visitors — excluded by
hypothesis. So the interior (minus two 1/w_f buffers) contains no g-wall, hence
is shorter than g's mesh 1/w_g. ∎

**Why the companion hypothesis cannot simply be dropped (the Davenport
observation):** the seven non-f inverses are seven nonzero residues of Z_7, and
the Davenport constant D(Z_7) = 7 means a zero-sum (balanced) subset ALWAYS
exists algebraically. The obstruction to long runs is therefore not algebraic
scarcity of balanced clusters but GEOMETRY: balanced owners must co-land in one
f-window at every g-wall, fixed pairs de-phase by (4), and handovers demand
cascading coincidences. The absolute run bound is exactly the statement that
this geometric-arithmetic conspiracy cannot persist.

## (6) The decision census (VERIFIED; the loophole eliminated)

Anchored-run maxima over quarter-period windows, heights to 10^4:

| family | n | median | max |
|---|---|---|---|
| generic random (w ≤ 10^4) | 80 | 3 | 5 |
| synchronized {2}+odds, Σw^{-1} ≡ 0 (the S303 loophole candidate) | 80 | 2 | 4 |
| near-ratio packets n·base ± 1 | 80 | 2 | 4 |
| all w ≡ r (mod 7) | 60 | 2 | 6 |
| annealed peak (600 steps from the best) | — | — | **6** |

The synchronized packet — the one family whose period-sums cancel by design —
performs WORSE than generic: the within-period φ-equations it cannot tune kill
it, exactly as (4)–(5) predict. **Working constant updated: K0 = 6** (THM-779's
census said 5 at heights ≤ 500).

## (7) Status of the exit lemma after this package

- UNCONDITIONAL and proved: no single-visitor periods (3); cluster balance (3);
  the de-phase bound (4); the extent bound under the no-companion hypothesis (5).
- The r = 8 pierce is now UNCONDITIONAL for every run configuration in which
  the second-fastest owner's walls fail to attract balanced co-landing clusters
  — and quantitatively constrained (de-phase + cascade) when they do.
- REMAINING (sharp, geometric): bound consecutive balanced co-landings
  absolutely — the cascade (n ≡ 4 or 0 tri-consistency, handover triple-balance,
  simultaneity exclusion) is the named machinery; the census says the true
  constant is 6 at all tested scales. Conjecture: **every blocking run has at
  most 6 walls; equivalently extent < C/w_g with small explicit C.**
- Everything here is Lean-friendly (mod-7 arithmetic + Beatty counting).
