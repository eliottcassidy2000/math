---
id: HYP-3594
title: THE TRUTH ACROSS ALL FRAMES (the odd-cycle reading, complementing klein-S14's value reading) -- every frame the project called fundamental is a measurement of ONE object, and that object is a single ODD CYCLE: the apex p-cycle C_p (p=7 for LRC14). Its VALUE is the spectral floor 4cos^2(3pi/7)=2+lambda_min(C_7)=0.198 (klein-S14's atom); its IDENTITY is an odd-length cycle in the R-odd cycle space; its MULTIPLICITY is the genus (1 at N=14 = ONE odd cycle); and the parent theory is literally the ODD-CYCLE Collection formula -- the project began with odd cycles (Redei/OCF, H=I(Omega,2)) and the LRC14 obstruction collapsed back to ONE odd cycle. ESSENTIAL x BOUNDED = "the odd cycle EXISTS (non-separable/non-bipartite, anti-Littlewood) AND its spectrum is positive (4cos^2(3pi/7)>0)". PLUS the two open tasks resolved: (T1) the rib SC-side-even theorem PROVED generally (R is an involution-automorphism; non-fixed neighbors of a fixed vertex pair up); (T2) the NS-rib defect (0,0,0,6 by tournament-n) is NOT the genus jump (0,0,1,2,2 by apex-prime-p) -- different index variables, shared R/even-odd root only (honest non-connection)
status: META-SYNTHESIS (complements klein-S14 HYP-3593) + (T1) PROVED general + verified n<=6 + (T2) COMPUTED + honest non-connection. The odd-cycle identity is grounded (4cos^2(3pi/7)=2+lambda_min(C_7), HYP-3590; cyclicity=3-cycle, THM-588; OCF=odd cycles). Not a proof of LRC14; the synthesis names the one object.
source: mac-mini-2026-06-29-S33
related:
  - HYP-3593  # klein-S14: the truth = one cyclotomic atom 4cos^2(3pi/7); frames = monotone zoom-in; exact inf R'=114382/332563
  - HYP-3590  # mac-mini S31: 4cos^2(3pi/7)=2+lambda_min(C_7) (the atom IS the 7-cycle, even-graph dual)
  - THM-588   # the metagraph's unique invariant = the 3-cycle (minimal odd cycle / cyclicity)
  - HYP-3587  # klein: genus = dim cusp forms = the obstruction multiplicity
  - HYP-3595  # mac-mini S32: the rib SC-side-even theorem (T1) + equinumerosity
  - THM-584   # complement R = the involution-automorphism (the engine of T1 and the R-odd cycle space)
results:
  - 04-computation/equinumerosity_and_eulerian_metagraph_macmini_20260629.py
---

# HYP-3594 -- the truth is a single odd cycle (the full circle)

The owner asked: across all the frames we called most fundamental, what truth were we aiming at?
klein-S14 (HYP-3593) answered with the VALUE: one cyclotomic atom `4cos^2(3pi/7)=0.198` at the apex-7
cusp, surrounded by the bulk floor `3/pi^2`; the 16+ frames are a monotone zoom-in on it. This is the
COMPLEMENTARY answer -- the atom's IDENTITY and its lineage.

## The truth: the atom is a single ODD CYCLE
The 17 frames (covering / metagraph / triangle / 2-adic / CV-variance / R-eigenspace / Euler-product /
Gamma_0(N) / relations / ESSENTIAL x BOUNDED / transitive-Z_7 / finite-cyclotomic / X_0(14)-cusps /
genus-cusp-form / odd-boundary / blue=SC / density-of-states) all measure ONE object, and that object is
a single **odd cycle** -- the apex `p`-cycle `C_p` (`p=7` for LRC14):
- **VALUE** (klein-S14): `4cos^2(3pi/7) = 2 + lambda_min(A(C_7)) = 0.198` -- the spectral floor of the
  7-cycle (HYP-3590). The doublet's autocorrelation IS `2I+A(C_7)`; the binding number is a 7-cycle
  eigenvalue.
- **IDENTITY**: an odd-LENGTH cycle, living in the R-odd / cycle space / even-graph dual `E_n` (it is a
  2-regular even GRAPH but an odd-LENGTH cycle). Detected by the complement involution `R` (the `-1`
  eigenspace = the obstruction).
- **MULTIPLICITY**: the **genus** of `X_0(2p)` = dim of cusp forms = the number of independent odd-cycle
  obstructions. Verified `0,0,1,2,2` (p=3,5,7,11,13). **LRC14 has genus 1 = exactly ONE odd cycle.** This
  is why N=14 is the first hard, first genuinely global case: one irreducible odd-cycle mode the
  local/bulk/Eisenstein data cannot determine.
- **LINEAGE / THE FULL CIRCLE**: the parent theory is literally the **Odd-Cycle Collection** formula
  (`H(T) = I(Omega,2)`, `Omega` = the directed ODD cycles). The metagraph's unique quadratic invariant is
  the **3-cycle** (minimal odd cycle, THM-588). The LRC14 floor binds at **C_7** (apex odd cycle). The
  project began with odd cycles and the obstruction collapsed back to a single odd cycle.

## Why an odd cycle is THE obstruction (the ESSENTIAL x BOUNDED reading)
An odd cycle is the obstruction to bipartiteness / 2-colorability / separability. So the two halves of the
proof (HYP-3570) are two facts about the SAME odd cycle:
- **ESSENTIAL** = the danger relation `D(v,t)` does not factor `f(v)g(t)` (rank > 1, non-separable, the
  bilinear `vt` is anti-Littlewood) = **D contains an odd cycle** (D is non-bipartite). The odd cycle IS
  the essentiality.
- **BOUNDED** = that odd cycle's spectral floor is positive = `4cos^2(3pi/7) > 0` (it does not reach the
  `gap=0`-only-at-full-`Z_7` disproof boundary).
So LRC14 = "one odd cycle exists and its spectrum is bounded away from zero." Both clauses, one cycle.

## The three fingerprints (= klein-S14's, read as the cycle)
1. bulk/obstruction split = the cut space (no odd cycle, Eisenstein, `3/pi^2`) vs the cycle space (the odd
   cycle, cusp form).
2. obstruction is 1-dim at N=14 = ONE odd cycle (genus 1).
3. cyclotomic at apex 7 = the cycle has odd length 7 (`Q(cos 2pi/7)`); the apex prime IS the cycle length.

## The two open tasks (resolved)
**T1 -- the rib SC-side-even theorem, PROVED in general.** `R` (complement) is an involution-automorphism
of the metagraph (THM-584). For ANY involution-automorphism and any `R`-fixed vertex `v`, the non-fixed
neighbors of `v` come in `R`-orbits of size 2; so `deg(v) ≡ #(R-fixed neighbors) (mod 2)`. An SC class is
`R`-fixed; its NS neighbors are the non-fixed ones -> **every SC class has even SC-NS (rib) degree**, all
`n` (verified SC-side-odd count `= 0` at n=3,4,5,6). This is a corollary of `R` detecting the cycle space
(it is the same `R` whose `-1`-eigenspace is the odd-cycle obstruction).
**T2 -- the NS break is NOT the genus jump (honest non-connection).** The NS-rib defect (#NS classes of
odd rib degree) is `0,0,0,6` by tournament size `n` (always even, by the handshake from T1); it first
appears at `n=6`. The genus is `0,0,1,2,2` by apex prime `p` (N=2p). These are **different index
variables** (tournament vertex count vs LRC modulus) with no natural identification -- the n=4 metagraph
(=X_0(14)'s 4 cusps, klein-S10) is FIXED across the genus family. They share only the deep `R` / even-odd /
cut-cycle root. Documented as a non-connection, not forced.

## What it buys
Names the one object the project kept rediscovering: a single odd cycle. klein-S14 gives its value and the
zoom-in; this gives its identity (odd cycle), its multiplicity (genus), and the full-circle narrative (OCF
-> one odd cycle). ESSENTIAL x BOUNDED = the odd cycle exists and is spectrally positive. The LRC14 is the
first case with exactly one such global odd cycle, and the proof is bounding its floor `4cos^2(3pi/7)` off
the disproof boundary.
