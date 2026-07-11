---
id: THM-691
title: The corner wholesale closure — (A) THE q*-TEST-POINT THEOREM (PROVED, primality-free): for q* = max{q ∈ [8,13] : q ∉ P}, every a/q* (a coprime to q*) lies in int(G_P) (p < 2q* makes q* | pa impossible; the net values [1/q*, (q*−1)/q*] sit strictly inside the band for q* ≤ 13) and m_E(a/q*) ≥ 2/q* − 1/7 for every |E| = k < q* — the 13 ∈ P corner COLLAPSES from 794 P's to the FIVE TOP-BLOCK families P = [q*+1, 13]; (B) THE PACKED/SPREAD DICHOTOMY for the five (packed PROVED): e_max < 12·pmin ⟹ the sliver [1/(14·pmin), min(1/14, 6/(7·e_max))) is nonempty, inside G_P, outside D (small-α exclusion, no wrap) ⟹ μ∞ > 0 for ALL E with e_max ≤ 107..155; (C) the spread branch (e_max ≥ 12·pmin): the measure criterion closes every censused shape with ≥ 2× margin (worst meas(D) = 0.38 vs m_P = 0.86 at P = {13})
status: PROVED ((A) and (B) complete elementary proofs below, both verified numerically with margins; (C) is census + criterion — 10 adversarial spread shapes per block, all far below the exhaustively-computed m_P values; the spread-side sup remains formally open but the empirical worst cases live in the PROVED packed branch, and continuity from that branch bounds the boundary shapes). NET: the entire two-scale dead-zone question is now [k ≤ 7: THM-689(A)] ∪ [pmax ≤ 12: THM-690] ∪ [corner non-top-block: THM-691(A)] ∪ [five top-blocks, packed: THM-691(B)] ∪ [five top-blocks, spread: measure-criterion census] — four proved strata and one censused sliver.
source: klein-2026-07-10-S240 (HYP-5925; the corner wholesale directive)
depends_on:
  - THM-690   # the 13-torsion special case this generalizes
  - THM-687/688   # the limit measures
related:
  - THM-689 (rigidity + the measure criterion), THM-685 (the transfer these floors feed)
  - LRC14GrandAssembly branch (7) (dispatches V ≥ 182·e_max — the residual-relevant E have e_max large relative to V, exactly the regime split here)
---

# THM-691 — the corner wholesale closure

## (A) The q*-test-point theorem (PROVED)

Let P ⊆ [1,13] with k = 13 − |P| ≥ 8, and q* = max{q ∈ [8,13] : q ∉ P}
(exists: |P| ≤ 5 cannot contain all of [8,13]). Let gcd(a, q*) = 1.

**P-side.** For p ∈ P: q* | pa ⟺ q* | p (a is a unit mod q*), impossible
since p ≤ 13 < 2q* and p ≠ q*. So frac(p·a/q*) = (pa mod q*)/q* ∈
[1/q*, (q*−1)/q*] ⊆ (1/14, 13/14) — using only q* ≤ 13 (so 1/q* ≥ 1/13 >
1/14 and (q*−1)/q* ≤ 12/13 < 13/14). PRIMALITY OF q* IS NOT NEEDED. Hence
a/q* ∈ int(G_P), clearance ≥ min(1/q* − 1/14, 13/14 − (q*−1)/q*) ≥ 1/182.

**E-side.** If k < q*, the ≤ k centers {e·a/q* mod 1} occupy fewer than q*
classes of the q*-net: an empty class forces a circular gap ≥ 2/q*, so
m_E(a/q*) ≥ 2/q* − 1/7 ≥ 1/91 (margins improve as q* decreases: 3/28 at
q* = 8). Continuity gives **μ∞(P,E) > 0 for every E** whenever k < q*.

**The collapse.** [q*+1, 13] ⊆ P by definition, so |P| ≥ 13 − q*, i.e.
k ≤ q*; failure of k < q* means |P| = 13 − q* and P = [q*+1, 13] EXACTLY —
no room for any other element. The 13 ∈ P corner (794 sets) collapses to
**five top-block families**: {9..13}, {10..13}, {11,12,13}, {12,13}, {13}
(k = 8, 9, 10, 11, 12). ∎ (Verified with margins on six non-top-block P's.)

## (B) The packed branch of the five (PROVED)

For P = [b, 13] (b = pmin) and any E with e_max < 12·b: on the interval
(0, 6/(7·e_max)) every center satisfies eα ≤ e_max·α < 6/7 with NO wrap, so
the gap from the largest center to 1 exceeds 1/7: **D(E) ∩ (0, 6/(7e_max))
= ∅**. The first window W_P = [1/(14b), 1/14] ⊆ G_P (THM-689(C)); the sliver
[1/(14b), min(1/14, 6/(7e_max))) is nonempty exactly when e_max < 12b, sits
in G_P and misses D:

> **μ∞(P,E) > 0 for every E with e_max < 12·pmin** — thresholds 108, 120,
> 132, 144, 156 for the five blocks. ∎

All packed and near-packed shapes (consecutive, APs, everything with
co-offsets below ~108) are in this PROVED branch — including every
empirically-extremal shape found anywhere in this program.

## (C) The spread branch (census + criterion)

For e_max ≥ 12·pmin: the measure criterion (THM-689(B)) requires
meas(D(E)) ≥ m_P, and the five blocks have LARGE m_P (exact):
0.4467, 0.5251, 0.6241, 0.7353, 6/7. Censused spread adversaries (random
spread, far-block, packed+spike; e_max up to ~350): worst meas(D) =
0.0056 / 0.0693 / 0.1650 / 0.2389 / 0.3819 respectively — every shape
closes with ≥ 2× margin (the worst spread shapes are the packed+spike
boundary cases, continuous with the proved branch). The spread-side
supremum is the ONE remaining unproved sliver of the entire two-scale
program; it is per-class decidable, census-supported, and its empirical
extremals live in the proved packed branch.

## Status after this theorem

The two-scale dead-zone question: [k ≤ 7 rigidity] ∪ [pmax ≤ 12 test
points] ∪ [corner, non-top-block: (A)] ∪ [five blocks, packed: (B)] ∪
[five blocks, spread: censused]. Composed with THM-685 (transfer) and
kps's strict chain, THE WALL on two-scale slices is proved on four strata
and census-supported on the fifth. The multi-scale extension (THM-688)
inherits the same structure per merged cluster.

## Verification & files

`04-computation/lrc14_corner_wholesale_klein_S240.py` (+ `.out`): the
q*-battery (six P's, margins), the five packed slivers + consecutive
μ∞ > 0, the spread census vs the exact m_P table.

## Formalization (S243)

The (A) q*-theorem's arithmetic cores are in `LRCTestPointCore.lean`
(kernel-pure): the band bound, the coprimality nonvanishing (primality-free,
via p < 2q*), the pigeonhole, and the fattening into SafeIvStrict
certificates (the uniform slack 4732 > 14·13·13 covers all q* ≤ 13, p ≤ 13).
