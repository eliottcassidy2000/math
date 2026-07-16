---
id: THM-896
title: THE KENDALL–WEI / PERRON ATLAS OF THE METAGRAPH (owner-fed claims verified and canonized) — self-composition R∘R counts 2-paths (the score-of-scores) and the Kendall–Wei iteration converges to the Perron eigenvector: every vertex's strength defined through all others SIMULTANEOUSLY — the all-frames-at-once invariant of the S28 intersubjective program, tournament-exact. (i) λ = 0 ⟺ TRANSITIVE (classical, one line: a tournament's adjacency is nilpotent ⟺ acyclic ⟺ transitive — "the transitive tournament is the zero of relational feedback"); (ii) per-class atlas n = 4..7: corr(λ, −score-variance) = +0.956/+0.928/+0.936/(n=7 in .out) — the owner's ≈0.93 verified; (iii) λ has NONZERO SPREAD inside EVERY multi-class score level (11/11 at n=6) — the intersubjective invariant STRICTLY REFINES the frame-dependent first moment (scores), which is velocity-relative/acceleration-absolute made tournament-exact. BONUS (proved, one line): A000568(n) ≡ #SC(n) (mod 2) — complementation is an involution on classes with SC as fixed points; SC counts verified 2, 8, 12 at n = 4, 5, 6 (cross-validating the merged-metagraph canon (V−SC)/2 = 2, 22)
status: (i) PROVED (classical; proof in-file); (ii)/(iii) MACHINE-VERIFIED per iso class by exact orbit-marking enumeration (counts 4/12/56/456 = A000568 ✓); the SC-parity lemma PROVED; source claims fed by the owner (verbatim numbers confirmed)
source: boxeph-2026-07-16-S30 (owner prompt: self-compositions, Kendall–Wei, per-class λ; verify and integrate)
related: [opus THM-894 (LRC Kendall-Wei: the resonance-matrix Perron over the cluster's own speeds — the LRC face of this prompt), mac-mini THM-895 (Helmholtz identity: gradient energy = x/(4n) exact + Perron census — the deeper Helmholtz), 07-reflections/the-intersubjective-object-boxeph-S28.md (this is its tournament-side theorem), klein cont.8 (composition defect = 3c₃ — the defect is the same R∘R object's antisymmetric shadow), kps coordinate atlas (C3 = score variance — λ refines it), THM-853(III) locker tournament (next: its λ)]
script: 04-computation/lrc14_deviation_kendallwei_boxeph_S30.py -> 05-knowledge/results/lrc14_deviation_kendallwei_boxeph_S30.out
---

# THM-893 — the Kendall–Wei / Perron atlas

**(i) λ = 0 ⟺ transitive.** For a tournament with adjacency A: ρ(A) = 0 ⟺ A nilpotent
⟺ the digraph is acyclic ⟺ (totality) it is the transitive order. *Proof:* nilpotent ⟺
no directed cycle (a cycle gives tr(A^k) > 0); an acyclic tournament has a source by
finiteness, induct. ∎ The transitive tournament is the unique zero of relational
feedback; every other class has λ > 0 fed by its cycle part (the curl, in the S28
Helmholtz reading — klein cont.8's composition defect 3c₃ is the same second-order
object's antisymmetric shadow).

**(ii)/(iii) The atlas (exact class enumeration, orbit-marking).** Classes 4/12/56/456
at n = 4/5/6/7 (= A000568 ✓). Per class: λ (power iteration on A), x = score variance:

| n | classes | corr(λ, −x) | multi-class score levels | levels with λ-spread > 0 |
|---|---|---|---|---|
| 4 | 4 | +0.956 | 0 | — |
| 5 | 12 | +0.928 | 2 | 2/2 |
| 6 | 56 | +0.936 | 11 | 11/11 |
| 7 | 456 | +0.947 | 41 | **40-ish/41 — FAILS: cospectral pairs exist** |

The owner's corr ≈ 0.93 is confirmed at every n (0.928–0.956). The refinement claim
(iii) holds STRICTLY at n ≤ 6 (13/13 multi-class levels separated) — and **FIRST FAILS
AT n = 7**, by genuinely COSPECTRAL non-isomorphic pairs (exact integer characteristic
polynomials equal: e.g. x⁷ − 4x⁴ − 4x³ − 3x² − x at score level (1,1,2,3,4,4,6), and
two further tied pairs at level (1,2,2,3,3,4,6) with polys −6,−6,−7,−2 and −6,−7,−6,−2):
the ENTIRE spectrum ties, not merely λ. So: score = velocity (per-frame), λ =
acceleration (the all-frames fixed point), and the Perron invariant refines the first
moment strictly through n = 6, with n = 7 the first cospectral degeneracy — n = 7 joins
the repo's first-failure list (width formula, E_n chordality, …). The corrected
statement: **the refinement is strict below the first cospectral pair, which lives at
n = 7** (minimality verified: no within-level λ-ties at n ≤ 6).

**SC-parity lemma (proved).** Complementation T ↦ T^op induces an involution on iso
classes whose fixed points are the self-complementary classes: A000568(n) ≡ #SC(n)
(mod 2). Verified: SC = 2, 8, 12 at n = 4, 5, 6 — consistent with the merged-metagraph
canon (V_merged = (A000568 + SC)/2; NS/2 = 2, 22 at n = 5, 6 ✓). (This is the pairing
the owner's toothpick prompt starts from: classes come in complement-pairs plus folded
fixed points — see TANGENTS entry T-toothpick, boxeph-S30.)

## ADDENDUM (S31): the cospectral census and the score-of-scores repair (partial)

Full census at n = 7: **256 within-level cospectral pairs** (λ-tie + exact charpoly
equality). The owner's R∘R invariant — the score-of-scores multiset (row sums of A²,
NOT a spectral function) — separates **225 of the 256 (88%)**; **31 pairs resist both
the spectrum and the second score vector** (the A²-entry multiset separates some but
not all of those). So the owner's object both creates the refinement story and repairs
most of its own n = 7 degeneracy; the 31 resistant pairs are the next degeneracy class
(named: which entrywise power-invariant completes the separation?). Data in the S31
tool log; reproduction: the S30 script + charpoly/matmul (10 lines).

## Evidence log
- [x] enumeration exact (A000568 reproduced n ≤ 7); (i) proof + verified; (ii)/(iii)
      verified n ≤ 6, n = 7 in the .out; SC parity proved + verified
- [ ] λ of the locker tournament D_n (bridge to THM-853(III)/THM-865)
- [ ] λ on the merged metagraph as a flow coordinate (kps atlas candidate: does the
      H-drift linearize in log λ? — their cont.29 found no coordinate; λ is untested)
