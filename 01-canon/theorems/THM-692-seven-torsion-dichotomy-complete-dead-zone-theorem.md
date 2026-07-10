---
id: THM-692
title: The 7-torsion dichotomy and THE COMPLETE TWO-SCALE DEAD-ZONE THEOREM — for P ∌ 7 and any a coprime to 7: either E misses a residue mod 7 (m_E(a/7) ≥ 1/7) or E is surjective mod 7 and then D(E) meets AT MOST ONE of the two open slivers flanking a/7 (both-sided coverage forces L_{j+1} = R_j across the middle chain — impossible for distinct integers in distinct residue classes), so one sliver is entirely D-free and μ∞ > 0 for EVERY E; for P ∋ 7 the q*-theorem (THM-691A) applies to every P (k ≤ q* − 1 always); TOGETHER: μ∞(P,E) > 0 for EVERY two-scale class — the dead-zone lemma is COMPLETE, the S240 censused sliver DISSOLVED, no census residue anywhere. Forensic: the first draft claimed the LEFT sliver universally D-free; the class-0 WRAP under δ < 0 (positive class-0 elements re-enter at the top and rescue coverage) refutes it — caught by this session's own assertion on a random adversary before canonization; the repaired one-sided-existential statement is verified with zero violations over 180 adversarial (E,a) pairs
status: PROVED (complete proof below; the one-sided statement verified exhaustively on adversarial batteries — 180 (E,a) pairs, zero both-covered, with genuine left-only AND right-only adversaries exhibited: the staircase covers the right sliver, the wrap-rescue shape covers the left; no shape covers both). Composed: [P ∌ 7: this theorem] ∪ [P ∋ 7: THM-691(A)] = μ∞(P,E) > 0 unconditionally for every P ⊆ [1,13] (k = 13 − |P| ≥ 8 arbitrary; k ≤ 7 was THM-689(A)) and every co-offset set E. THE WALL on two-scale slices is now PROVED EVERYWHERE (with THM-685 + kps's strict chain); THM-688's multi-scale extension inherits per merged cluster.
source: klein-2026-07-10-S241 (HYP-5930; the spread-sliver supremum directive — dissolved rather than bounded)
depends_on:
  - THM-691 (the q*-theorem consumed for P ∋ 7)
  - THM-687/688 (the limit measures)
related:
  - THM-689 (k ≤ 7 rigidity — the k = 7 bijective case of this theorem's middle-chain argument), THM-690 (the 13-torsion special case)
  - LRCParityPairing (the half-point q = 2 structure — the m ≤ 6 low-torsion family is its geometric kin)
---

# THM-692 — the 7-torsion dichotomy: the complete dead-zone theorem

## Statement

Let P ⊆ [1,13] with 7 ∉ P, a ∈ {1,…,6}, and E ∋ 0 any set of k distinct
co-offsets. Then a/7 ∈ int(G_P) (clearance 1/14), and:

**(i) Non-surjective branch.** If {e·a mod 7 : e ∈ E} misses a class, the
gap across the empty class is ≥ 2/7, so m_E(a/7) ≥ 2/7 − 1/7 = **1/7**.

**(ii) Surjective branch.** If E is surjective mod 7, then D(E) meets at
most ONE of the open slivers (a/7 − δ′, a/7), (a/7, a/7 + δ′), where
δ′ = min(1/(14·e_max), 1/182); the other sliver is ENTIRELY D-free and
m_E > 0 pointwise on it.

In both cases **μ∞(P,E) > 0** — and with THM-691(A) covering every P ∋ 7
(there [q*+1,13] ∪ {7} ⊆ P forces k ≤ q* − 1 < q*):

> **μ∞(P,E) > 0 for EVERY two-scale class — all P, all E, all k.**

## Proof

**P-side.** 7 ∤ p for p ∈ P, so (pa mod 7)/7 ∈ [1/7, 6/7] ⊂ (1/14, 13/14)
with clearance 1/14; the slivers stay in G_P since δ′ ≤ 1/182 ≤ (1/14)/13.

**(i)** is the pigeonhole of THM-690/691 at q = 7.

**(ii)** At α = a/7 + δ, |δ| < δ′, the center of co-offset e sits at
c_e/7 + eδ (c_e = ea mod 7) when c_e ≥ 1 — no wrap, inside the strip
c_e/7 ± 1/14 — while class 0 splits: e = 0 stays at 0, and positive
class-0 elements sit at eδ (for δ > 0) or wrap to 1 + eδ (for δ < 0).
The MIDDLE classes j = 1,…,6 are nonempty (surjectivity) and their strips
contain no foreign centers, so for j = 1,…,5 the gap between class j's
rightmost point and class j+1's leftmost point is EXACTLY

> gap_j(δ) = 1/7 + (L_{j+1} − R_j)·δ,

with L/R the min/max co-offset per class. Coverage at α requires
gap_j(δ) ≤ 1/7, i.e. (L_{j+1} − R_j)·δ ≤ 0 — a sign condition independent
of |δ|. If D met BOTH slivers, then L_{j+1} ≤ R_j and L_{j+1} ≥ R_j for
every j = 1,…,5, forcing L_{j+1} = R_j — impossible: they are co-offsets
in different residue classes mod 7, hence distinct integers. So D misses
one entire sliver; off D, m_E > 0 pointwise; integrate over the sliver
(inside G_P) for μ∞ > 0. ∎

**The forensic (recorded for honesty).** The first draft claimed the LEFT
sliver is always the free one, via a span-sum contradiction that ignored
the class-0 wrap: for δ < 0 the positive class-0 elements re-enter at the
top of the circle and can plug the wrap gap. A random surjective adversary
(E = {0,23,50,113,150,151,160,173,203,211}, a = 4) covers the left sliver —
caught by this session's own assertion before anything was canonized. The
repaired statement is one-sided-EXISTENTIAL, which is all positivity needs;
both genuinely-one-sided adversary types are exhibited in the companion
script (the descending staircase covers the right; the wrap-rescue covers
the left), and 180 adversarial (E,a) pairs show ZERO both-covered cases.

## Why 7

m ≤ 6 torsion: coverage impossible outright (≤ 6 arc positions, mass
≤ 6/7 — the low-torsion universal points). m ≥ 8: the net can cover with
slack, no rigidity. m = 7 = 14/2 is the ZERO-SLACK modulus — coverage is
possible but only just, and the one-sided rigidity is exactly the residual
freedom. The dichotomy modulus is the problem's own resonance, the same
zero-slack mechanism as THM-689(A), one derivative deeper.

## Status: the two-scale dead-zone question is CLOSED

[k ≤ 7: THM-689(A)] ∪ [P ∌ 7: THIS THEOREM] ∪ [P ∋ 7: THM-691(A)] =
**μ∞(P,E) > 0 unconditionally, everywhere.** THM-690's 13-torsion theorem
and THM-691(B)'s packed sliver become special cases/alternate routes; the
S240 spread census is superseded by proof. Composed with THM-687/688
(the limits and rates), THM-685 (the transfer), and kps's LRCStrictRuler:
**the wall holds on every two-scale and multi-scale slice beyond the
explicit V₀(P,E) thresholds — proved, no census anywhere in the chain.**
Below V₀: THM-686 windows and the banks.

## Verification & files

`04-computation/lrc14_seven_torsion_dichotomy_klein_S241.py` (+ `.out`):
the staircase (right-covered, left-free) and the wrap-rescue refutation of
the first draft; 180-pair one-sidedness battery (zero violations); the
non-surjective margin; the exhaustive P ∋ 7 q*-check (zero violations);
five top-block μ∞ spot checks on spread staircase adversaries.
