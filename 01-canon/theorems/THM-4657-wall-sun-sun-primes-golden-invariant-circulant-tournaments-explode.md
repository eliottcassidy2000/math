---
id: THM-4657
title: "CEDED to THM-4656 (duplicate). Independent confirmation of THM-4656 parts 1 and 4 for split primes: Wall-Sun-Sun primes are the Wieferich primes to the base phi, the circulant tournaments on Z_(p^2) invariant under the odd part of phi explode exactly there, and the golden Paley primes 11, 19, 31, ... (52 below 2000)"
status: >
  SUPERSEDED BY THM-4656 (mac-mini-2026-10-14-tourn reserved THM-4656 at 2b6d2c3abf, before this number was reserved at
  d00166b48d, and pushed the full text at 448c0d367e). The number is ceded; no priority is claimed. This file is kept only
  as an independent confirmation. CITED: Wall-Sun-Sun = Wieferich to the base phi (Bragman-Rowland 2022, Theorem 10;
  the name "Fibonacci-Wieferich prime", McIntosh-Roettger 2007). PROVED in THM-4656 (independently re-proved here):
  the split-prime orbit count. FINITE-EXACT (this session): the orbit count at all 40 (p, phi) pairs with split p <= 199;
  no Wall-Sun-Sun prime below 2*10^6. CITED: none below 2^64 (PrimeGrid, 2020-2022).
source: opus-2026-10-14-S36
depends_on: [THM-4656, THM-4644 (part 2), Wall 1960, Sun-Sun 1992, Bragman-Rowland 2022]
verification: 04-computation/experiments/fibonacci_pisano_wall_sun_sun_20261014.py (output 05-knowledge/results/fibonacci_pisano_wall_sun_sun_20261014.out); audit 04-computation/experiments/s36_audit_20261014/audit_pisano_wss_orbits.py
note: 05-knowledge/results/strong_balance_pisano_golden_eleven_20261014.md
related: [THM-4656, THM-4644, THM-4631]
---

# THM-4657 — ceded to THM-4656

**Use THM-4656.** mac-mini's THM-4656 (Wall–Sun–Sun primes as golden Wieferich primes) contains everything this number was reserved for, and more.
- Its part 1 is the split-prime statement below, for every metallic unit `μ_a`.
- Its part 2 is the inert case, which this session left unwritten.
- Its part 4 lists the same 52 golden Paley primes below 2000.

This session found the same statements independently, in parallel, before reading THM-4656. The independent agreement is the only reason the file is kept.

## What was confirmed independently

Let `p ≡ ±1 (mod 5)` be a split prime and `φ ∈ Z/p²` a Hensel lift of a root of `x² − x − 1`. Put `d = ord_p(φ)`, let `d'` be the odd part of `d`, and let `g = φ^(d/d')` be the odd part of `φ`.

1. **Wieferich form (CITED).** `p` is Wall–Sun–Sun ⟺ `φ^(p−1) ≡ 1 (mod p²)`, the same for `ψ = −φ^(−1)`.
   - This is known: Bragman–Rowland (arXiv:2202.00704, Theorem 10) and the name "Fibonacci–Wieferich prime" (McIntosh–Roettger 2007).
   - The proof is one line. `F_(p−1) = (φ^(p−1) − ψ^(p−1))/(φ − ψ)`, and `ψ^(p−1) = φ^(−(p−1))`, so `F_(p−1) ≡ ±2pt/√5 (mod p²)` when `φ^(p−1) = 1 + pt`.
2. **Orbit count (THM-4656 part 1).** The circulant tournaments `S` on `Z_(p²)` with `gS = S` number `2^E`, where
   - `E = (p − 1)/d'` if `p` is not Wall–Sun–Sun;
   - `E = (p + 1)(p − 1)/(2d')` if it is.

   *Proof.* `⟨g⟩` has odd order `d'p^ε`, with `ε = 0` exactly at Wall–Sun–Sun primes. So `−1 ∉ ⟨g⟩` and negation pairs the orbits freely. There are `p(p−1)/(d'p^ε)` orbits on units and `(p−1)/d'` on the non-zero non-units. A `g`-invariant tournament picks one orbit from each pair.
3. **Golden Paley primes (THM-4656 part 4).** Let `p ≡ 3 (mod 4)` split, and suppose the root that is a square mod `p` generates `QR_p`. Then the invariant circulants on `Z_p` are exactly `P_p` and its reverse.
   - Below 2000 there are 52 such primes: `11, 19, 31, 59, 71, 79, 131, …`.
   - At `p = 11` the square root is `4`, with `⟨4⟩ = {1, 3, 4, 5, 9}`. Here `⟨x + 1, 4x⟩ = Aut(P_11)` has order 55 and acts regularly on the 55 arcs.

**Wording corrections from the audit.**
- "Golden-invariant" means invariant under the odd part `g` of `φ`. When the order of `φ` is even, no tournament is invariant under `φ` itself, because `−1 ∈ ⟨φ⟩`.
- At golden Paley primes, `ψ = −1/φ` is a primitive root.
- The PrimeGrid bound is `2^64 ≈ 1.8·10^19` (2020–2022), not "past `4·10^16`" (a stale OEIS A244801 comment).

**Independent numerics.** The orbit count was enumerated at all 40 `(p, φ)` pairs with split `p <= 199`, in agreement with the formula. The S36 audit re-derived the count with its own code. It also tested the explosion branch on 39 artificial Wieferich-type units, since no Wall–Sun–Sun prime is known.
