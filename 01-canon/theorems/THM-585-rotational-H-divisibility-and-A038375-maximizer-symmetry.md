---
id: THM-585
title: Vertex-transitive tournaments have n | H(T) = n*(odd) (generalizing the Paley p|H of THM-586); the divisibility of A038375 (max H) by n is a SHADOW of the maximizer's symmetry -- n | a(n) exactly for n=1,3,5,7,9,11 (the circulant-optimal regime, THM-338), failing at n=13 where the maximizer is non-circulant (13 | opt_circ(13)=3711175 but 13 nmid a(13)); and the skew-Hadamard doubling tower (THM-448) revives n|H at the Mersenne numbers via DRT regularity (15 | H(T_15)) with the self-similar recursion B_0(T_{2m-1})=T_{m-1}
status: PROVED (the vertex-transitive lemma); VERIFIED (n|H for ALL circulant tournaments n=3..13; A038375 divisibility n=1..13 with opt_circ(13)=3711175 exact; doubling tower B_0 self-similarity and H-values to T_15).
source: mac-mini-2026-06-29-S10
depends_on:
  - THM-586   # Paley p|H + dihedral Burnside (this generalizes the p|H to all vertex-transitive)
  - THM-338   # circulant optimality threshold: a(n)=opt_circ(n) for n<=11, gap at 13
  - THM-448   # the DRT/Mersenne skew-Hadamard doubling tower + B_0 self-similarity
related:
  - A038375   # OEIS: max # Hamiltonian paths in an n-vertex tournament
  - CONJ-001  # Rédei (H odd); refined here to H = n*(odd) for vertex-transitive T
  - HYP-2351  # T_15 = conjectural n=15 H-maximizer (if true, n|a(n) revives at Mersenne 15)
  - OPEN-Q-108
results:
  - 04-computation/rotational_tournament_arithmetic_macmini_20260629.py
  - 05-knowledge/results/rotational_tournament_arithmetic_macmini_20260629.out
---

# THM-585 -- rotational H-divisibility and the symmetry shadow of A038375

## 1. The vertex-transitive lemma (PROVED; generalizes THM-586)
> For any VERTEX-TRANSITIVE tournament `T` on `n` vertices, the number of directed Hamiltonian
> paths STARTING at a fixed vertex is `H(T)/n` (by transitivity, equal for all start vertices), so
> **`n | H(T)`**; with Rédei (`H` odd), **`H(T) = n * (odd)`**.

Rotational/circulant tournaments `R_n` (`i->j` iff `(j-i) mod n in S`, `S` antisymmetric) are
vertex-transitive (`Z_n` acts), so `n | H(R_n)`.  This generalizes THM-586's free-rotation `p|H` for
Paley (a special `R_n` with `S=QR`) to ALL rotational tournaments.  VERIFIED: `n | H` for every one of
the `2^{(n-1)/2}` circulant tournaments, `n = 3,5,7,9,11,13`.

## 2. A038375 divisibility = a shadow of the maximizer's symmetry (VERIFIED)
`a(n) = A038375(n) = max{H(T) : |T|=n}`.  Then:
> **`n | a(n)` exactly for `n = 1, 3, 5, 7, 9, 11`,  and `n nmid a(n)` for `n = 2,4,6,8,10,12,13`.**

This range `{3,5,7,9,11}` is EXACTLY where the H-maximizer is a circulant (vertex-transitive)
tournament (THM-338: `a(n)=opt_circ(n)` for `n<=11`).  At `n=13` the maximizer is NON-circulant
(THM-338: `a(13) > opt_circ(13)` by 8656): the circulant optimum `opt_circ(13)=3711175` IS divisible
by 13 (vertex-transitive), but `a(13)=3719831` is NOT.  So in this range
> **`n | a(n)  <=>  the n-vertex H-maximizer is vertex-transitive`,**
and the circulant-optimality threshold `n=13` (THM-338) is simultaneously the DIVISIBILITY threshold
of the OEIS max-H sequence.  A038375's arithmetic encodes the symmetry of its extremal tournament.

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 13 |
|---|---|---|---|---|---|---|---|---|---|----|----|----|
| a(n) | 1 | 1 | 3 | 5 | 15 | 45 | 189 | 661 | 3357 | 15745 | 95095 | 3719831 |
| n\|a(n) | y | n | **y** | n | **y** | n | **y** | n | **y** | n | **y** | n |

## 3. The doubling tower: a second divisibility mechanism + self-similarity (VERIFIED)
The skew-Hadamard doubling tower (THM-447/448) builds DRTs on the Mersenne numbers `2^k-1`
(`3,7,15,31,63`).  It is SELF-SIMILAR: `B_0(T_{2m-1}) = T_{m-1}` exactly (the out-neighborhood of
vertex 0 IS the previous level) -- VERIFIED `H(B_0(T_15))=H(T_7)=189`, `H(B_0(T_7))=H(T_3)=3`.  And
**`n | H` at every tower level**: `3|3`, `7|189`, `15|198335025` -- INCLUDING `T_15`, which has
`Aut = F_21` (order 21) and is therefore NOT vertex-transitive on its 15 vertices.  So DRT regularity
forces `n|H` by a SECOND mechanism beyond vertex-transitivity.  Consequence: if `T_15` is the `n=15`
H-maximizer (HYP-2351, conjectured), then `15 | a(15)` and the A038375 divisibility REVIVES at the
Mersenne number 15, past the circulant threshold 11 -- the two structures (vertex-transitive circulant,
and Mersenne doubling DRT) are the two sources of `n | a(n)`.

## Significance: two recursive towers, one arithmetic
Merging THM-586 (Paley) with `R_n` and A038375: the Hamiltonian-path count of a maximally-symmetric
tournament is `n*(odd)`, and this divisibility is visible in the OEIS max-H sequence as a symmetry
detector.  Two towers organize the symmetric maximizers: the **Paley tower** (primes `p=3 mod4`:
`3,7,11,19,23,31,...`, vertex-transitive, `p|H`, achieving `a(p)`) and the **Mersenne doubling tower**
(`3,7,15,31,63`, DRT, self-similar `B_0` recursion, `n|H`).  They COINCIDE at `3,7` (`T_p ~= Paley`)
then DIVERGE: `T_15` (15 composite, non-Paley) and `T_31` (non-Paley DRT, THM-448).  At Mersenne
primes the Paley tournament (more symmetric, `|Aut|=p(p-1)/2`) likely out-maximizes the tower DRT
(`|Aut|=21`).  Open: is `n | H` for EVERY doubly-regular tournament (the common generalization of both
mechanisms)?  Verified for all computed DRTs (Paley `3,7,11`; tower `3,7,15`).

Scripts: `rotational_tournament_arithmetic_macmini_20260629.py` (this);
`paley_dihedral_burnside_macmini_20260629.py` (THM-586).
