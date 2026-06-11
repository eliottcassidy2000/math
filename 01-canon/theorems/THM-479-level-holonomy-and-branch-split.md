# THM-479 — Level permutations via pair-orbit holonomy; the branch-split Burnside for A049313

**Status:** PROVED (claudebox-2026-06-11-S2) — an independent re-derivation, by a new
elementary argument, of Babai–Cameron Lemma 7.1 + Theorem 7.2 (EJC 7 (2000) #R38; proved
1981-82, manuscript lost, published 2000), found here BEFORE the literature check and
verified exhaustively (all cycle types n ≤ 10, zero discrepancies; closed form matches
OEIS A049313 terms 1..16 exactly). The branch split (parts 3-4) appears to be new.
**Provenance:** claudebox-2026-06-11-S2, dispatch: RM(r,m)⊥ = RM(m−r−1,m) duality as a lens.
**Companions:** THM-474 (tilings = labeled switching classes), THM-477 (blue code),
THM-430 (antipodal involution), OPEN-Q-060. **Literature:** Babai–Cameron EJC 2000
(Lemma 7.1 "level", Thm 5.1, Thm 7.2, Lemma 3.1); Mallows–Sloane SIAM JAM 28 (1975);
Cameron, Math. Z. 157 (1977); OEIS A049313 (Cameron 1999), A002854.

## Setup

Fix the reference tournament T₀ (i→j for i<j) and identify tournaments on [n] with flip
vectors x ∈ 𝔽₂^E, E = E(K_n). S_n acts AFFINELY: x ↦ P_π x + t_π, where P_π permutes pair
coordinates and the cocycle t_π(e) = 1 iff transporting T₀'s orientation along π reverses it
at e. Switching translates by the cut space C (dim n−1). The graph case is the same picture
with t ≡ 0: tournaments are the odd affine torsor over the even/linear graph world (LEM-004).
π fixes some switching class of tournaments ⟺ ∃x, c: (P_π+1)x = t_π + c ⟺ t_π ∈ Im(P_π+1)+C.

## Statement

1. **(Level law; = BC Lemma 7.1.)** π fixes a switching class of tournaments iff π is
   **level**: all cycle lengths have the same 2-adic valuation. (Equivalently: all cycles
   odd — i.e. π of odd order — or all cycles even with equal v₂.)
2. **(Burnside; = BC Thm 7.2.)** With o₂(μ) = Σᵢ⌊ℓᵢ/2⌋ + Σ_{i<j} gcd(ℓᵢ,ℓⱼ) (pair-orbits)
   and k(μ) = #cycles, the number of switching classes of tournaments up to isomorphism is
   A049313(n) = Σ_{μ all-odd} 2^{o₂−k+1}/z_μ + Σ_{μ all-even, equal v₂} 2^{o₂−k}/z_μ.
3. **(Branch split — new.)** Write A049313(n) = N_odd(n) + N_lev(n) for the two sums in (2).
   N_lev vanishes for odd n (no all-even partitions); for even n ≥ 4 both branches are
   separately INTEGERS (verified n ≤ 16; only n = 2 splits fractionally as 1/2 + 1/2):
   N_odd = 1, 1/2*, 1, 1, 2, 4, 12, 64, 792, 19296, 886288, 75346496, … and
   N_lev(2t) = 1/2*, 1, 2, 15, 280, 23464, 6819360, 6981258440 (t = 1..8). Neither branch
   sequence, nor the count of level permutations itself, is in OEIS (checked 2026-06-11) —
   flagged for contribution. Branch integrality for n ≥ 3 is observed, not yet explained.
4. **(Where the obstruction lives.)** The proof shows the cocycle's holonomy is supported
   EXACTLY on the diameter pair-orbits {i, π^{ℓ/2}(i)} of even cycles — the antipodal /
   half-turn pairs (THM-430's σ-self-partners; THM-447's twin pairs; the same pairs whose
   tiling avatar underlies the blue-line canon). Mixed cycle types die because an even
   cycle demands odd cut-mark-parity while any odd cycle in its company forces it even.

## Proof

**Orbit reduction.** (P+1) acts blockwise on the π-orbits of pairs; on a cyclic orbit of
length ℓ_O its kernel is the constants and its image the even-weight vectors. Hence
t ∈ Im(P+1) + C ⟺ ∃ cut c with Σ_{e∈O}(t+c)_e = 0 for every pair-orbit O.

**Holonomy of t.** Σ_{e∈O} t_π(e) = 1 iff the orbit returns the pair REVERSED. Within a
cycle of length ℓ, the pair {i, i+d} (0<d<ℓ) returns reversed iff d = ℓ/2 (the diameter
orbit, length ℓ/2; even ℓ only); all other within-cycle orbits (full length ℓ) and all
cross-cycle orbits (length lcm) return the pair unswapped. So Σ_O t = 1 exactly on the one
diameter orbit of each even cycle.

**Cut sums.** For c = cut(W) with cycle parities sᵢ = |W ∩ cycleᵢ| mod 2:
Σ_O c = 0 on within-cycle non-diameter orbits; = sᵢ on the diameter orbit of even cycle i;
= (ℓⱼ/g)sᵢ + (ℓᵢ/g)sⱼ on cross orbits (g = gcd(ℓᵢ,ℓⱼ)).

**The system.** Solvability ⟺ ∃ s ∈ 𝔽₂^k with sᵢ = 1 for every even cycle and
(ℓⱼ/g)sᵢ + (ℓᵢ/g)sⱼ = 0 for all i<j. Odd–odd pairs: both quotients odd ⟹ sᵢ = sⱼ.
Even(i)–odd(j) pairs: g odd ⟹ ℓⱼ/g odd, ℓᵢ/g even ⟹ sᵢ = 0, contradicting sᵢ = 1 ⟹
mixed types unsolvable. All odd: s ≡ const works. All even: s ≡ 1 forced; cross condition
becomes ℓᵢ/g + ℓⱼ/g ≡ 0, which holds iff v₂(ℓᵢ) = v₂(ℓⱼ) (equal v₂ ⟹ both quotients odd;
unequal ⟹ exactly one even). This is the level law. ∎(1)

**Counting.** For solvable π, #fixed classes = 2^{o₂ + dim(C∩Im) − (n−1)}, and cut(W) ∈
Im(P+1) ⟺ s(W) satisfies the homogeneous system: sᵢ = 0 on even cycles, all odd-cycle sᵢ
equal — giving dim(C∩Im) = n − k if some cycle is odd, n − k − 1 if all are even, whence
the exponents o₂−k+1 and o₂−k. Burnside over S_n gives (2). ∎

## Verification (script switching_classes_level_burnside_cbx2.py)

- Direct 𝔽₂ linear algebra (bitmask Gaussian elimination, no closed form) over ALL cycle
  types n ≤ 10: solvability ≡ level law, 0 violations; Burnside totals 1,1,1,2,2,6,12,79,
  792,19576 (n=1..10).
- Closed form (2) reproduces OEIS A049313 terms n = 1..16 exactly, including
  1938818712501985736 at n = 16 (independent validation against Jovovic's 2000 terms).
- Targeted all-even types at n = 12: (6,6), (6,2,2,2), (4,4,4) fixable; (8,4) not —
  exactly the equal-v₂ pattern.
- Graph control (t ≡ 0, same code path): reproduces A002854 = 1,1,2,3,7,16,54,243,2038,
  33120 (Mallows–Sloane).

## Remarks

1. **The RM dictionary (the dispatch).** The graph world is the 𝔽₂-LINEAR story: cut space
   ⊥ cycle space inside 𝔽₂^E is the K_n-avatar of RM(r,m)⊥ = RM(m−r−1,m) (repetition↔SPC
   ends; in repo coordinates: base-path arcs = cut/score data, wiggly arcs = cycle data —
   already canon). Mallows–Sloane's theorem (#two-graphs = #Euler graphs) is the duality's
   Burnside shadow. The tournament world is the AFFINE TORSOR over it: the cocycle t_π is
   the odd twist, its cohomology class obstructs fixing, and the obstruction is 2-adic —
   "level" is the seam. Graphs : Euler graphs :: tournaments : (OPEN-Q-060 — the branch
   split in (3) says the answer must itself split by the 2-adic branch).
2. **Self-duality and n = 4.** dim C = n−1 equals dim Z = C(n−1,2) only at n = 4 (K₄ is
   self-dual as a planar graph): the k = n/2 self-dual point of the dispatch, sitting at
   the first even-branch order with nontrivial level structure.
3. **Contrast with graphs (BC Lemma 3.1).** Every single automorphism of a graph switching
   class fixes a member graph; a tournament-class automorphism fixes a member iff it has
   odd order. The even-level branch is symmetry existing ONLY at class level — perspective
   without a fixed representative, the σ-exile motif (S643) at the switching layer.
