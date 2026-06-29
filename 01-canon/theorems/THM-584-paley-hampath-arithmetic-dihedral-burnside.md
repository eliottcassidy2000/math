---
id: THM-584
title: The arithmetic of the Paley-tournament Hamiltonian-path count -- for p prime = 3 mod 4, the free Z_p rotation gives p | H(T_p) with v_p=1, so H(T_p) = p*(odd) (a Paley refinement of Rédei); the full dihedral D_{2p} gives the Burnside identity #orbits = (H(T_p) + p*f(T_p))/(2p) in Z, where f is the palindromic count (THM-583), forcing H(T_p)/p == f (mod 2); R(p)=H*2^{p-1}/p! -> e; and the conjecture H = |Aut|*3^{(p-3)/2} is FALSE (fails at p=11, cofactor 1729=7*13*19, not 81)
status: PROVED (Z_p free action and Burnside on D_{2p}; the equal-fix of all p reflections is conjugacy in D_p, p odd). VERIFIED exhaustively p=3,7,11 (H=3,189,95095; f=1,9,185; #orbits=1,18,4415) and H/f data computed to p=31.
source: mac-mini-2026-06-29-S9
depends_on:
  - THM-582   # f = #palindromic Ham paths is odd
  - THM-583   # f via the half-system (the reflection-fixed count)
related:
  - CONJ-001  # Claim A / Rédei (H odd); this refines it for Paley
  - why-the-paley-path-ratio-is-e-the-cherry-is-the-unique-cluster   # R(p)->e
  - H-impossibility-the-multiplicative-mechanism-s599   # H-spectrum multiplicativity (Moon)
  - OPEN-Q-108
results:
  - 04-computation/paley_tournament_number_theory_macmini_20260629.py
  - 04-computation/paley_dihedral_burnside_macmini_20260629.py
---

# THM-584 -- the arithmetic of the Paley Hamiltonian-path count

## Setup
Paley tournament `T_p`, `p` prime `= 3 mod 4`: vertices `Z_p`, `i->j` iff `j-i` is a nonzero QR.
`H(T_p)` = # directed Hamiltonian paths (Rédei: odd).  `f(T_p)` = # palindromic Hamiltonian paths
(THM-582/583, the reversal-fixed count; `f=1,9,185,573057,...` for `p=3,7,11,19,...`).  The symmetry
group is `D_{2p} = <r: v->v+1 (automorphism, order p), s: v->-v (involutory anti-automorphism)>`,
acting on directed Hamiltonian paths: rotations directly (`r^k.P=(v_1+k,...,v_n+k)`), reflections by
reversal-conjugation (`rho_g(P)=g(reverse(P))`).

## 1. `p | H(T_p)` and `H(T_p) = p*(odd)` (PROVED)
The rotation `Z_p` acts FREELY on Hamiltonian paths (no directed path is invariant under a nonzero
rotation: rotation moves the endpoints), so orbits have size `p` and **`p | H(T_p)`**, with
`v_p(H)=1` (verified: `v_p=1` at `p=3,7,11,19`).  Combined with Rédei (`H` odd):
> **`H(T_p) = p * (odd number)`** -- a strict refinement of Rédei for Paley tournaments.
VERIFIED: `H/p = 1, 27, 8645, 61720828785` (all odd) at `p=3,7,11,19`.

## 2. The dihedral Burnside identity (PROVED + verified)
All `p` reflections of `D_{2p}` are involutory anti-automorphisms, conjugate in `D_p` (`p` odd), so
each fixes the SAME number `f(T_p)` of palindromic paths.  Burnside (`Fix(id)=H`, `Fix(r^k)=0` for
`k!=0`, `Fix(reflection)=f`):
> **`#{Hamiltonian-path orbits under D_{2p}} = (H(T_p) + p * f(T_p)) / (2p)`  in `Z`.**
VERIFIED: `(3+3)/6=1`, `(189+63)/14=18`, `(95095+2035)/22=4415`.  Consequence: `2p | (H + p f)`, and
since `p|H` this gives `H/p + f == 0 (mod 2)`, i.e.
> **`H(T_p)/p == f(T_p) (mod 2)`** -- both ODD (`f` odd by THM-582).
This ties the path count, the palindromic count, and the prime `p` through the dihedral symmetry.

## 3. The asymptotic `R(p) -> e` (consistency)
`R(p) := H(T_p) * 2^{p-1} / p! = 2.0, 2.4, 2.4395, 2.5272, ...` (p=3,7,11,19) `-> e = 2.71828...`
(the cherry-cluster result, reflection `why-the-paley-path-ratio-is-e`).  These `H`-values are the
exact inputs; the Burnside structure (1) is the `p|H` factor that survives normalization.

## 4. The `3^{(p-3)/2}` conjecture is FALSE (correction)
The repo conjecture `H(T_p) = |Aut(T_p)| * 3^{(p-3)/2}` (`|Aut|=p(p-1)/2`) holds only for `p=3,7`
(`H=3,189`).  At `p=11` it predicts `55*81=4455` but `H(T_11)=95095=55*1729`, cofactor
`1729 = 7*13*19` (the taxicab number), NOT `81=3^4`.  At `p=19` the cofactor is `5*7*11*23*774463`.
So the clean `3^k` form is a small-`p` coincidence; the true cofactor `H/|Aut|` is a product of
assorted primes (and `R(p)->e`, not a `3`-power growth).  DO NOT USE the `3^{(p-3)/2}` formula.

## Significance
A self-contained tournament<->prime result on the canonical bridge (Paley `T_p` exists exactly for
`p=3 mod 4`).  It refines Rédei (`H=p*odd`), realizes the dihedral symmetry as a Burnside identity
linking `H`, the palindromic `f` (mac-mini's recent THM-582/583 odd index), and `p`, and corrects a
standing conjecture.  It sits in the number-theoretic web mapped this session: the `H`-values that
`->e` (Paley ratio), the H-spectrum multiplicativity (Moon `H=prod H(C_i)`), and the totient/`zeta(2)`
floor.  Open: a closed form / linear recurrence for `H(T_p)/p` and `f(T_p)` (`1,27,8645,...` and
`1,9,185,...`); their QR/Gauss-sum governance.

Scripts: `paley_tournament_number_theory_macmini_20260629.py`,
`paley_dihedral_burnside_macmini_20260629.py`.
