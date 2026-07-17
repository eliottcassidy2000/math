---
id: LEM-035
title: THE SURVIVOR LAW FOR FAMILY CLUSTERS (LEM-034's family70 CRT named check RESOLVED, and generalized). (A) SECTION-VECTOR FORMULA: at a class-0 crossing x = k/(7m) of owner e = 7m, every other runner f has section ⌊fk/m⌋ mod 7, is at a boundary iff m | fk (the LEM-034 lattice), left section subtracts [m | fk]; e sits in section 6 on the left. SURVIVOR CRITERION (x⁻ ∈ R₀): the six other left-sections avoid 0 and cover {1..5} (slack slot: 6 covered by e). Refereed EXACT on 9 clusters (8 family + balanced/near-AP/two-large; all S63 survivor sets re-derived). (B) THE r = 1 UNIVERSAL FAMILY (the family70 progression PROVED): for E = {1..6, 7M}, M ≥ 7, k = Mj + 1: carries vanish, left-vector = (j, 2j, …, 6j) = the multiplication-by-j permutation of {1..6} ⟹ survivor ⟺ j ≢ 0 (mod 7). Exactly 6 survivors k ≡ 1 (mod M), k ≠ 1 — independent of M. (C) CLEAN + RESCUED COLUMNS: every r with 6r < M is a full survivor column (same proof); if 6 | M, the column r = M/6 is boundary-rescued (runner 6's carry 1 cancelled by its lattice boundary: (6j+1)−1 = 6j). Hence the survivor count grows ~6·(#columns) ≈ M-linearly. Verified: M = 13,14,15 gain the r = 2 column (measured extras {Mj+2}); M = 12's r = 2 column rescued. (D) SPORADICS: isolated (r, j) = (M−2, 5) at M = 11 (k = 64, predicted by hand from the carry vector BEFORE the census — confirmed) and M = 12 (k = 70): the f = 6 carry anomaly (⌊6(M−2)/M⌋ = 4, resp. boundary) turns the fatal 6j+5 ≡ 0 into 6; dead for M ≥ 13. None at M ≤ 10 (M = 10 completeness = the S63 pattern, proved by the automated case analysis). (E) ATTRIBUTION INTERPLAY: family84's whole r = 2 column is stolen by runner 6 (6k ≡ 0 mod 12 at every k = 12j+2): 13 survivors → N₀ = −6; family98's is not (needs f = 7 ∉ E): N₀ = −12. All 8 family N₀ values reproduced from survivor + attribution counts. (F) REFLECTION CLOSURE: survivor exits at k/70 ↦ entries of R₆ at (70−k)/70 with class 0 = (6+1 mod 7) — LEM-034(A)'s adjacent-entry class, verified
status: PROVED ((A) 2 lines; (B)/(C) permutation argument + carry/boundary cancellation; (D) finite verification per M; (E) one-line divisibility) + REFEREED EXACT (9 clusters, all crossings, census == formula everywhere; M = 11 prediction confirmed; N₀ ×8 reproduced)
source: boxeph-2026-07-17-S64 (owner directive: the family70 CRT structure; keep proving little statements)
depends_on: [LEM-034 (lattice + sign law + the named check), S25/S26 endpoint machinery]
related: [LEM-030 (σ_e decomposition), THM-886 (N_c combs)]
script: 04-computation/lrc14_family70_crt_boxeph_S64.py -> 05-knowledge/results/lrc14_family70_crt_boxeph_S64.out
---

# LEM-035 — the survivor law

At x = k/(7m): 7fx = fk/m, so σ_f = ⌊fk/m⌋ mod 7, boundary ⟺ m | fk, left
subtracts 1 at boundaries. Writing k = Mj + r on E = {1..6, 7M}: σ_f = fj +
⌊fr/M⌋ − [M | fr] on the left. The whole geometry of ∂R₀ at the big owner's
integer crossings reduces to a mod-7 carry computation:

- **r with 6r < M** (and r = M/6 when 6 | M, boundary-rescued): left-vector =
  (fj)_f = the multiplication-by-j permutation of {1..6}: survivor ⟺ j ≢ 0.
  The family70 progression k ≡ 1 (mod 10), k ≠ 1 is the M = 10 instance.
- **k = 1 (j = 0)** puts all six smalls in section 0 — the origin shadow.
- **Sporadics** (M−2, 5) at M = 11, 12 via the f = 6 carry anomaly; the k = 64
  survivor on family77 was predicted from the carry vector before measuring.
- Attribution then relabels stolen columns (runner 6 owns family84's r = 2
  column); the surviving N₀ ledger matches exactly.

## Evidence log
- [x] formula == census on 9 clusters; all S63 survivor sets re-derived
- [x] r = 1 family all M = 8..15; M = 10 completeness; M = 11 k = 64 prediction
- [x] clean/rescued column law at M = 12..15; sporadic mechanism at M = 11, 12
- [x] N₀ reproduced ×8; reflection closure at family70
- [ ] named: sporadic classification for general E (beyond {1..6, 7M}) — the
      slack-slot combinatorics at arbitrary carry vectors
