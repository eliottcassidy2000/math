---
id: HYP-3715
title: WORKING THE LRC->HEXAGONAL BRIDGE -- the covering-min's binding configuration IS a zeta_6-direction LINE in the hexagonal lattice: M(S)=gap=covering radius j/D (HYP-3551), covering-min 14/183=n/Phi_6(n) at {1..12,182} (densest core + lcm(n-1,n) killer); at the binding rotation a=zeta_6=14 the speeds map to the n-spaced AP {-14,14,28,...,168} = the runners equally spaced (spacing n) along the zeta_6 (hexagonal 60deg) direction in Z/Phi_6(n), covering radius n. SUBCONDITIONS: (1) DENSE covering sets -- the covering-min lives at the densest coverable core, and the densest core {1..n-2} alone has gap 1/(n-1) > 1/n, so the dense-core family satisfies M > 1/n (partial conjecture, HYP-2566/3551); (2) ANTIPODAL binding pair (1, n(n-1)=-1 mod D); (3) the lcm(n-1,n) MINIMAL KILLER (coprime->largest->equidistributes->minimal perturbation, gives D=Phi_6(n)); (4) DISCRETE KERSHNER -- the tightest covering uses maximally-spread speeds. The covering-min is the zeta_6-line; full optimality (no covering beats it) = the hexagonal-line optimality, connected to Kershner
status: VERIFIED (M=14/183 & 7/89 reproduced via covering radius; densest-core-wins scan; the zeta_6-AP binding config; discrete-Kershner spread minimizes covering radius). The dense-core partial bound (M>1/n on dense-core coverings) is established (HYP-2566/3551). The global covering-min (no exotic covering beats 14/183) + the cyclic-line<->2D-Kershner step remain OPEN.
source: klein-2026-06-29-S26
depends_on:
  - HYP-3551   # M = covering radius j/D; 14/183 = the covering-min; densest core + minimal killer
  - HYP-3706   # the hexagonal/Eisenstein wallpaper bridge (zeta_6 = mult-by-n = the 60deg rotation)
related:
  - HYP-2566   # inf M over covering sets > 1/14 (uniform looseness) -- the dense-core bound
  - HYP-3710   # the Coxeter-Catalan A1->A2 ladder (hexagonal = A2)
  - THM-523    # covering = kill all resonances b<=14
results:
  - 04-computation/hex_bridge_subconditions_klein.py
  - 05-knowledge/results/hex_bridge_subconditions_klein.out
---

# HYP-3715 — the covering-min is a zeta_6-line in the hexagonal lattice; the subconditions

## M = covering radius; the covering-min reproduced
`M(S) = the lonely-runner gap = max over rational t=a/D of min_s ||s a/D|| = j/D` at the binding modulus
`D` (HYP-3551). Verified: `M({1..12,182}) = 14/183` (binding rotation `a=14`, `D=183=Phi_6(14)`);
`M({1..11,13,84}) = 7/89`. Densest-core scan (`{1..13}\{skip}` + killer `lcm(skip,14)`): the minimum is
`skip=13` -> `14/183` (the densest coverable core `{1..12}` + minimal killer `182=lcm(13,14)`).

## The binding configuration IS a zeta_6-direction line (the concrete hexagonal realization)
At the binding rotation `a = 14 = zeta_6` (`14^6 = 1 mod 183`, the 6-fold hexagonal rotation), the speeds map
to
> `S * 14 mod 183 = {14, 28, 42, ..., 168} u {169}  =  {-14, 14, 28, ..., 168}`
the `n`-spaced arithmetic progression (spacing `14 = n`) along the `zeta_6` (hexagonal `60`-deg) direction in
`Z/Phi_6(n) = Z[zeta_6]/(n-zeta_6)`, missing `0`; the two points nearest `0` are `+-14`, so the covering
radius is `n = 14` and `M = n/Phi_6(n)`. **So the LRC covering-min is, concretely, the runners equally
spaced along a line in the `zeta_6` direction of the hexagonal lattice** -- a 1D sub-lattice (the
`zeta_6`-line), spacing `n`, modulus `Phi_6(n)`. The abstract "hexagonal quotient" (HYP-3706) becomes this
explicit line.

## The subconditions (the structural skeleton of the covering-min)
1. **DENSE covering sets** (the strongest lever). The covering-min lives at the *densest coverable core*
   plus the minimal killer (densest-core-wins, verified). And the densest core `{1..n-2}` ALONE has lonely
   gap `1/(n-1)` (`= 1/13 > 1/14` for `n=14`); the killer only perturbs it slightly. So the **dense-core
   family satisfies `M > 1/n`** -- the conjecture holds on dense coverings (HYP-2566/3551 uniform
   looseness). This is exactly where the covering-min sits AND where a partial bound is provable.
2. **ANTIPODAL binding** `(1, -1)`. The covering-min's binding pair is `(1, n(n-1))` with `n(n-1) = -1 mod
   Phi_6(n)`; modulus `D = Phi_6(n) = 1 + n(n-1)`. The simplest possible binding (speed `1` vs its negative).
3. **The `lcm(n-1, n)` MINIMAL KILLER**. The killer `n(n-1)` (`= lcm` since `gcd(n-1,n)=1`) is the largest
   minimal killer -> equidistributes -> perturbs the dense core least -> `D = Phi_6(n)` (Eisenstein norm).
4. **DISCRETE KERSHNER**. For `k` speeds in `Z/D` the tightest covering (min covering radius `~ D/2k`) uses
   MAXIMALLY-SPREAD speeds (verified: minimizers are spread sets). For the LRC covering-min the spread is the
   `zeta_6`-line in `Z/Phi_6(n)` -- the discrete analogue of the hexagonal (Kershner) optimum.

## Creative special cases to consider next
- the `zeta_6`-line vs the full 2D hexagonal lattice: the covering-min is a *line* (1D AP) in the hexagonal
  lattice, not the 2D Kershner covering -- the bridge is "the optimal LRC line is the `zeta_6` line";
- prime-power `n-1` (Singer/projective-plane exists) vs general `n`;
- the antipodal-binding family `(1, -1 mod D)` for general `D = 1 + n(n-1)`;
- the renormalization fixed points (the doublet attractor, HYP-3700's apex column) vs this covering column.

## Honest status
- RIGOROUS/VERIFIED: `M` = covering radius; `14/183 = n/Phi_6(n)` at the densest-core + `lcm` killer; the
  `zeta_6`-line binding configuration; discrete Kershner (spread minimizes covering radius); the dense-core
  partial bound `M > 1/n` (HYP-2566/3551).
- OPEN: that `14/183` is the GLOBAL covering-min (no exotic, non-dense-core covering beats it); and the step
  from the cyclic `zeta_6`-line covering radius to the 2D hexagonal (Kershner) optimality. The bridge is now
  a concrete geometric claim: **the optimal LRC covering is the `zeta_6`-line in the hexagonal lattice.**
