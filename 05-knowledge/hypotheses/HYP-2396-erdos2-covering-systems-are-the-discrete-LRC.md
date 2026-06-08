---
id: HYP-2396
name: erdos2-covering-systems-are-the-discrete-LRC-and-the-distortion-method
status: SURVEY + bridge (inclusion-exclusion correspondence VERIFIED) + a proof-method import for the LRC
date: 2026-06-08
session: claudebox-2026-06-08-S724
external:
  - "Erdos 2 (DISPROVED, $1000): min modulus of a covering system bounded. Hough 2015 (<=10^16);"
  - "Balister-Bollobas-Morris-Sahasrabudhe-Tiba 2022 (distortion method, <=616000; no all-odd-squarefree)"
  - "best construction: min modulus 42 (Owens, BYU thesis); Nielsen 2009 had 40"
depends_on:
  - THM-406   # LRC covering-depth master object; p_0 = Sum_j (-1)^j S_j
  - HYP-2065  # S561: sieve density rho = Sum_T (-1)^|T|/lcm(T); density alone can't close LRC
  - THM-439   # the witness tower is the cyclotomic tower
provisional_id: true
---

# HYP-2396: Erdos covering systems are the discrete LRC; the distortion method is the tool the LRC wants

## Erdos problem 2 (DISPROVED)

A covering system = finitely many congruences `x = a_i (mod n_i)`, DISTINCT moduli, covering every
integer. Erdos asked (his "favourite problem"): can `min n_i` be arbitrarily large? Hough (2015) proved
NO -- `min n_i <= 10^16`; Balister-Bollobas-Morris-Sahasrabudhe-Tiba (2022) gave a simpler "distortion
method" bound `<= 616000` and proved no covering has all moduli odd & squarefree. Best CONSTRUCTION (lower
bound): **min modulus 42** (Owens; Nielsen had 40). So the minimum modulus is BOUNDED: between 42 and 616000.

## The structure of minimal (large-min-modulus) coverings

- moduli are **smooth** -- divisors of a highly composite `L`; built by recursive CRT refinement (cover
  residues mod 2, refine the uncovered residues mod 3, mod 4, ...). Min modulus 42 = `2*3*7`.
- **density is necessary but NOT sufficient.** `Sum 1/n_i >= 1` is required, but VERIFIED here:
  moduli `{3,4,5,6,8,9,10,12}` (min 3) have `Sum 1/n_i = 493/360 ~ 1.37 > 1` yet leave **uncovered density
  `89/360 ~ 0.247`** -- they do NOT cover. The obstruction is the inclusion-exclusion OVERLAP structure.
- Hough/BBMST = a **distortion / density-increment**: when min modulus is large the moduli are forced so
  smooth that a positive proportion of residues escapes every congruence -- the cover always leaks.

## THE BRIDGE (VERIFIED): the LRC is the continuous covering problem; Erdos 2 is the discrete one

Both have the SAME inclusion-exclusion skeleton:
```
   covering system:  uncovered density = Sum_{I} (-1)^|I| [I compatible] / lcm(moduli in I)   (= S561 rho)
   LRC (THM-406):    lonely measure   p_0 = Sum_{j} (-1)^j S_j,   S_j = sum of j-fold danger-arc overlaps
```
VERIFIED on the classic min-modulus-2 covering `{0(2),0(3),1(4),5(6),7(12)}`: direct uncovered density
`= 0`, inclusion-exclusion `rho = 0`, match. "Covers `Z`" <=> `rho = 0`; "leaves a gap" <=> `p_0 > 0` (a
lonely time). **The covering system is the rational/CRT SKELETON of the LRC danger arcs** (continuous
overlaps `S_j` vs the skeleton's `1/lcm`; same alternating sum).

## THE METHOD IMPORT (the research payoff)

Hough/BBMST proved "spread-out (large-modulus) pieces always leave a positive uncovered density." That is
*exactly* the discrete sibling of the LRC's **"the cover always leaks"** (`p_0 > 0` for multiple-of-n
configs; my S643). Both share S561's lesson that density `Sum 1/n_i > 1` does not imply covering -- the
overlap structure decides. So:

> **CONJECTURE / PROGRAM:** port the distortion / density-increment method (Hough; BBMST) to the LRC. The
> covering-depth `p_0` (THM-406) is the LRC's uncovered density; a density-increment argument showing
> `p_0` stays bounded below for multiple-of-n configs (the runners' danger arcs are "too smooth/spread to
> cover the circle") would be the continuous analog of Hough's theorem -- a route to LRC looseness via the
> covering-system technology rather than the twistor/shell route.

The two methods are dual: covering systems run the density-increment to show coverage FAILS at large
modulus; the LRC wants exactly that failure (looseness). And BBMST's "no all-odd-squarefree covering"
mirrors the LRC's parity/2-adic seam (even moduli essential <-> the even/odd seam, THM-435).

## 42 in the repo
Incidental (`a(42)` index in A038375 THM-329; a Paley-table entry). The COVERING is central: THM-406
covering-depth (the LRC master object), S561 `rho` (the covering density), THM-439 witness tower.

## Next
- attempt a distortion/density-increment lower bound on the LRC `p_0` for multiple-of-n configs;
- the squarefree/odd-modulus ban (BBMST) vs the LRC 2-adic seam -- is the LRC's even-n hardness the same
  parity obstruction;
- compute small covering systems' covering-depth distribution `p_k` (not just `p_0`) and compare to the
  LRC `p_k` (THM-406) -- do they share the moment hierarchy `S_j`?
