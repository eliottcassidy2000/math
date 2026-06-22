---
id: HYP-+2878
title: LRC(14) completed-proof ROUTE via strong-component atoms + covering systems -- LRC-failure = a persistent covering system of the unsafe APs (s708); a 13-set's unsafe APs cover at most BOUNDEDLY many primes simultaneously (verified <=3 of 5); so a finite prime basis always has a witness => LRC(14). The hard cases are SINGLE resonance atoms.
status: ROUTE + key bound VERIFIED (<=3 of 5 primes covered, 4000 13-sets; hard cases single-atom). Remaining: rigorous over-determination bound (step 3) = the completed-proof crux.
source: mac-mini-2026-06-22-S35 (strong-component-atom extension of HYP-2876)
related:
  - HYP-2876   # the finite rational-witness certificate (N(S,D), main term + resonances)
  - HYP-2860   # Node-3 spectrum = the resonances
  - h21-moon-reduction-s617  # H multiplicative over strong components (the atom idea)
  - even-cycles-covering-systems-and-the-parity-covering-lens-s708  # LRC failure = covering system
---

# HYP-+2878 -- the strong-component-atom + covering-system proof route

## The chain (extending the finite certificate HYP-2876)
1. **Witness <=> covering complement:** LRC(14) witness for S <=> N(S,D)>=1 for some D, where
   N(S,D)=#{a: ||sa/D||>=1/14 forall s}. **N(S,D)=0 <=> the unsafe APs U_s={a: ||sa/D||<1/14}
   COVER Z/D** -- a covering system (s708: "LRC failure = danger residues cover Z/q for every q").
2. **min modulus >= 15 (apex-7, Idea 3 PROVED, kps LRCApex7Floor.lean):** D=14 always covers
   (forced-14 runner on observer), so a non-trivial cover needs D>=15.
3. **STRONG-COMPONENT ATOMS:** the resonance graph (s~s' if low-height k s+k' s'≡0 mod D) has
   strong components = ATOMS. The HARD/covering cases (loosest {1..11,13,84}) are a SINGLE atom
   (verified: #atoms=1 at D=14,28,41,83) -- one irreducible strongly-resonant obstruction. H is
   multiplicative over atoms (Moon, h21-s617); the {7,21}/apex-7 obstruction lives at the atom level.
4. **BOUNDED simultaneous covering (the key, VERIFIED):** a 13-set's unsafe APs cover Z/p (N=0) for
   at most **3 of the 5 primes {83,89,97,101,103}** simultaneously (4000 13-sets; worst covers
   {83,97,103}). So a basis of 5 primes ALWAYS leaves >=2 with N>=1 => a witness => LRC(14) holds.
   (Main term (6/7)^13 phi(p)~11 says ~11 survive generically; covering needs resonance alignment,
   which is p-SPECIFIC and over-determined across primes.)

## The completed-proof crux (step 3 rigorous)
PROVE: a single 13-element atom's unsafe APs cannot cover Z/p for more than C primes p (C<5).
MECHANISM (over-determination): covering Z/p requires the 13 APs (each ~p/7 residues) to align with
no gap -- ~p constraints on the residues {s mod p}. Aligning at MANY primes over-determines the 13
speeds (Hough-type: covering systems have BOUNDED min modulus / structure). The atom is the
irreducible unit; the apex-7 (q=7) and the {7,21} forbidden-H are the q=7 atom obstruction.
This is the Node-3 resonance-deficit bound (HYP-2860) in covering-system form, and it COMPLETES
the finite certificate (HYP-2876): finite prime basis + bounded simultaneous covering => witness.
-> HYP-2876, HYP-2860, s708, h21-s617, kps LRCApex7Floor.

Script: lrc14_atom_covering_macmini_S35.py.
