---
id: THM-852
title: THE BLACK SELF-LINE LAW, STRUCTURED — 2·selfK = SC(n) (opus S311's law) carries a KLEIN-FOUR and INVOLUTION skeleton. Setup: κ = complement tiling (T(κt) = P ∪ Ā = "op modulo the base path": it differs from T(t)^op = P^op ∪ Ā exactly by re-reversing the base), g = the grid reflection (realizes op at class level: cls(gt) = cls(t)^op); quasi-fixed X = {t : cls(κt) = cls(t)}, SC-ness of cls(t) ⟺ cls(gκt) = cls(t). PROVED HERE: (i) the Klein four-group {1, g, κ, gκ} acts on X, κ freely; (ii) THE INVOLUTION LEMMA: any witness coset σ·Aut(T) with σ² ∈ Aut(T) and |Aut| ODD (always, tournaments) contains an element of 2-power order; when Aut is trivial the unique witness is an INVOLUTION — at odd n ALL quasi-fixed classes have trivial Aut (verified n=5,7: Burnside-affine W = |X| = 8, 88), so every black self-line carries a unique involutive κ-twisted symmetry, the exact analogue of the SC classes' involutive anti-automorphism; (iii) the Klein-4 orbits on X_non-gridsym ALL have full size 4 (n=5: 2 orbits; n=6: 3; the law reads SC(n) = 4·#orbits); (iv) THE n=8 CHECK REFUTES THE LAW: non-gridsym quasi-fixed = 404 ≠ 176 = SC(8) (total X = 412, gridsym-qf = 8) — the 'all n ≥ 5' law was a small-n coincidence killed at the first even n beyond its data; the S311 total-law SC + 2selfB also fails (412 ≠ 184). Pipeline certified: the identical two-pass code reproduces 8/12/88 at n=5,6,7 and gridsym-qf = 0 at odd n
status: STRUCTURE PROVED ((i),(ii) odd-order-coset; (iii) verified n=5,6) + **THE LAW REFUTED AT n=8**: 404 ≠ 176 (run complete, pipeline certified against 8/12/88 at n=5,6,7). The law 2selfK = SC holds at n = 5,6,7 ONLY as far as is known — there is NO all-n bijection to prove; the 5,6,7 agreement is either coincidence or the shadow of a modified law (404 = 4·101 orbits-worth; the Klein-orbit granularity SURVIVES: 404 ≡ 0 mod 4). The open question is now: what does 2selfK count in general? (n=8 data: selfK = 202, selfB = 4, X = 412.)
source: kind-pasteur-2026-07-15-S128 (cont.13, overnight; owner: prove the 2selfK = SC bijection, check n=8)
depends_on: []
related:
  - opus S306/S310/S311 (HYP-6860/6885/6890): the law's discovery arc (odd-n restriction was an artifact; carriers not all SC ⟹ weighted bijection needed)
  - THM-791 (the even-n fixed-point-free involution mechanism — the SAME odd-|Aut| forcing produces involutions here)
  - THM-466/OCF (Aut odd via Rédei-adjacent parity), THM-793 (Mode-B tower: the κ/g actions descend it)
---

# THM-852 — structure of the black self-line law

**The dictionary.** T(κt) = P ∪ Ā and T(t)^op = P^op ∪ Ā: the complement tiling is "converse
modulo the base path". Hence: quasi-fixed (self-line) ⟺ cls(P∪A) = cls(P∪Ā); SC ⟺ cls(P∪A) =
cls(P^op∪Ā) = cls(gκt). The law 2selfK = SC equates the sizes of two twisted diagonals of the
same Klein-four action — off the grid-symmetric locus.

**Involution lemma.** If σ(T) = T' with T' ≅ T via σ and σ² ∈ Aut(T), and |Aut(T)| is odd (every
tournament, by Rédei-parity of automorphisms), then the coset σ·Aut(T) contains an element of
2-power order; if Aut(T) = 1, σ itself satisfies σ² = 1. (⟨σ⟩ has order 2^a·b, b odd dividing
into Aut-data; σ^b is a witness of 2-power order.) Consequently at odd n — where the referee finds
every quasi-fixed class has TRIVIAL Aut (W = |X| exactly at n=5, 7) — each black self-line carries
a UNIQUE involutive witness σ, mirroring the unique involutive anti-automorphism π of a trivial-Aut
SC class. The conjectured bijection should match (t, σ) ↔ (K, π) data; the Klein-orbit granularity
(SC = 4·#orbits) says each orbit {t, gt, κt, gκt} accounts for exactly four SC classes.

## Evidence log

- [x] Klein-4 action + free κ + orbit sizes (all 4): exact n=5,6
- [x] Burnside-affine W(n) = Σ_{t∈X}|Aut| = 8, 20, 88 (n=5,6,7): odd-n triviality of Aut on X
- [ ] n=8: #non-gridsym quasi-fixed =? 176 = SC(8)   [two-pass run in flight]
- [ ] the orbit ↦ 4-SC-classes bijection (open crux)
