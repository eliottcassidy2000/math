#!/usr/bin/env python3
"""n8_mass_babai_cameron_boxeph_S164.py (HYP-8315) — code as run (S164 heredocs),
frozen outs: n8_mass_hunt_boxeph_S164.out + babai_cameron_switching_boxeph_S164b.out.
FRONT A: K8 mass = 3/8 CONFIRMED exact (cycle types); BONUS breakers at n=8:
K4+K4 mass 5/8 (swap coset uniformly positive: (1/4+1)/2) and 4K2 mass 3/8
(f+ = 13/16 and 11/16 — sixteenth fractions): quarters break at n=8 through BOTH
mechanisms; random classes stay quarter (rarity = high symmetry).
FRONT B: Babai-Cameron Remark 7.4 ('We cannot do this') at small n:
switching classes 1/2/12 at n = 3/5/7; BC counts = 1/0/2. CONFIRMED: n = 1 mod 4
gives BC = 0 via the UNIQUE Eulerian member per class ({1: 2} at n=5); n = 3 mod 4
has ZERO Eulerian members (C(n,2) odd) — no canonical anchor => memberless
automorphisms exist (1 at n=3, 2 at n=7). METHOD: triangle-parity t-vectors with
the TWISTED affine S_n action t(sigma x) = sigma_3(t(x)) xor d(sigma),
d = t(inversion set) — the twist caught by an equivariance assert (first run
wrongly assumed plain equivariance). n=11 needs cycle-index Burnside (handoff)."""
print(open("05-knowledge/results/n8_mass_hunt_boxeph_S164.out").read())
print(open("05-knowledge/results/babai_cameron_switching_boxeph_S164b.out").read())
