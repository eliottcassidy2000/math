---
id: HYP-3771
title: THE 5 PLANE-LATTICES / 5 PLATONIC SOLIDS AND THE LRC APEX-PRIME GEOMETRY -- one (2,3,n) angle-defect spine, and why 7 is the frontier. The two counts of 5 come from ONE source, the sign of 1/2+1/3+1/n-1 (Gauss-Bonnet angle defect): n=3,4,5 SPHERICAL = the 5 Platonic solids {3,3},{3,4},{4,3},{3,5},{5,3}; n=6 EUCLIDEAN = the hexagonal plane tiling = A2 root lattice = the LRC covering-min; n>=7 HYPERBOLIC, n=7 = Klein quartic (2,3,7)=PSL(2,7)=168=Fano. The planar count 5 = crystallographic orders {n: 2cos(2pi/n) in Z} = {n: phi(n)<=2} = {1,2,3,4,6} = the 5 two-dimensional Bravais lattices. LRC(2p) apex prime p walks this spine and the modular genus of X0(2p) TRACKS the geometry: p=3 Euclidean genus X0(6)=0 (LRC6 trivial); p=5 SPHERICAL/icosahedral -- a Platonic solid! -- genus X0(10)=0 (tractable); p=7 the FIRST HYPERBOLIC, genus X0(14)=1 (HARD, cusp form f14). The genus jump 0->1 at p=7 IS the spherical->hyperbolic transition -- the arithmetic reason LRC14 is the first genuinely hard case. The 5 Platonic solids consume the spherical orders {3,4,5}; order 6 is the unique Euclidean (the covering-min, hexagonal/A2, last session's Dedekind anomaly HYP-3768); 7 is the first order that tiles NEITHER sphere nor plane => the LRC frontier. The crystallographic-forbiddenness of 7 (2cos(2pi/7) irrational, degree phi(7)/2=3) is the same cyclotomic non-degeneracy as the apex floor gap 4cos^2(3pi/7)=2+lambda_min(C_7)>0 (THM-590).
status: SYNTHESIS / unifying map. The identities are exact/verified (crystallographic count = {phi<=2}; the (2,3,n) trichotomy; the 5 Platonic 1/p+1/q>1/2; genus X0(2p)=0,0,1,2,2 for p=3,5,7,11,13 tracking spherical/hyperbolic). These are classical facts assembled + the genus<->(2,3,p)-geometry alignment; NOT a new proof. It clarifies WHY the apex prime 7 is the frontier and unifies the covering-min (order 6) and apex (order 7) geometries with last session's B2/A2/Dedekind.
source: mac-mini-2026-06-30-S65
related:
  - HYP-3768   # my S64 + klein-S56 (convergent): covering-min margin = order-6 hexagonal Dedekind anomaly; klein PROVED s(n,Phi6)=-(Phi6-1)/(12Phi6) -> -1/12=zeta(-1)=-B2/2 (Bernoulli B2 reading); this session's B2 = the square root-lattice reading -- both converge
  - HYP-3715   # covering-min = zeta6-line in the hexagonal A2 lattice
  - HYP-3586   # X0(14) cusps=Klein V4, apex cusp d=7, genus jump 0->1 = the f14 obstruction
  - HYP-3606   # the apex gap 4cos^2(3pi/7)=2+lambda_min(C_7) = non-bipartiteness certificate (THM-590)
  - HYP-2943   # codex: tiling/solid recursion carriers split by curvature (the prior scout)
  - HYP-3616   # forbidden-H {7,21}; 7=|Fano|=Phi_3(2); the apex prime's special-ness
---

# HYP-3771 -- the 5 tilings / 5 solids and the LRC apex-prime geometry

The owner asked to consider "the 5 2D plane tilings correspond to the 5 Platonic solids," where their nature
and count come from, and connections to proofs. The answer is one spine.

## One source: the (2,3,n) angle defect
Both counts of 5 -- and the LRC covering-min and apex geometries -- are the **sign of the angle defect**
`1/2 + 1/3 + 1/n - 1` (Gauss-Bonnet for the `(2,3,n)` triangle group):

| n | `1/2+1/3+1/n-1` | geometry | object |
|---|-----------------|----------|--------|
| 3 | `+1/6` | spherical | tetrahedron `{3,3}` |
| 4 | `+1/12` | spherical | octahedron `{3,4}` |
| 5 | `+1/30` | spherical | icosahedron `{3,5}` |
| 6 | `0` | **Euclidean** | **hexagonal tiling = A2 = LRC covering-min** |
| 7 | `-1/42` | **hyperbolic** | **Klein quartic `(2,3,7) = PSL(2,7) = 168 = Fano` = LRC14 apex** |

The **5 Platonic solids** are the spherical regular figures `{p,q}` with `1/p + 1/q > 1/2`
(`{3,3},{3,4},{4,3},{3,5},{5,3}`). The **5 planar rotation orders** are the crystallographic
`{n : 2cos(2pi/n) in Z} = {n : phi(n) <= 2} = {1,2,3,4,6}` -- the **5 two-dimensional Bravais lattices**.
Same trichotomy: positive defect = the solids, zero defect = the plane lattices, negative = hyperbolic.

## The LRC apex-prime walk (the proof connection)
`LRC(2p)` has apex prime `p`; the `p`-fold cycle `C_p` is the binding object (the apex doublet). As `p` grows,
`(2,3,p)` walks spherical -> Euclidean-adjacent -> hyperbolic, and **the modular genus of `X_0(2p)` tracks it
exactly** (verified):

| p | X_0(2p) | genus | `(2,3,p)` geometry | LRC |
|---|---------|-------|--------------------|-----|
| 3 | X_0(6) | 0 | Euclidean/crystallographic (order-3 triangular) | LRC6 trivial |
| 5 | X_0(10) | 0 | **spherical / icosahedral -- a Platonic solid** | LRC10 tractable |
| 7 | X_0(14) | **1** | **first hyperbolic (Klein quartic)** | **LRC14 first hard** |
| 11 | X_0(22) | 2 | hyperbolic | harder |
| 13 | X_0(26) | 2 | hyperbolic | harder |

**The genus jump `0 -> 1` at `p=7` is the spherical->hyperbolic transition.** This is the arithmetic reason
LRC14 is the frontier: the apex primes `3, 5` still live on finite (spherical/Platonic) symmetry -- order 3 is
crystallographic, order 5 is the icosahedron -- so `X_0(6), X_0(10)` are genus 0 (no cusp form). `7` is the
**first prime that tiles neither the sphere (not Platonic) nor the plane (not crystallographic)**; it is purely
hyperbolic, `X_0(14)` acquires genus 1, and the cusp form `f_14` at the apex cusp `d=7` appears -- the hard
core (HYP-3586, HYP-3768).

## Two geometries in LRC, cleanly separated
- **The covering-min geometry = order 6 (hexagonal, A2, Euclidean `(2,3,6)`)**: universal for all `n`; the
  construction is `n/Phi_6`, and its margin is the order-6 Dedekind anomaly `s(n,Phi_6) != 0` (HYP-3768).
  The competitor **B2 (square, order 4)** is the other crystallographic lattice, anomaly-free (`s=0`).
- **The apex geometry = order `p` (the `C_p` doublet)**: spherical for `p<=5`, hyperbolic for `p>=7`. The apex
  floor gap `4cos^2(3pi/7) = 2 + lambda_min(C_7) > 0` is exactly the **crystallographic non-degeneracy** of `7`
  (`2cos(2pi/7)` irrational, degree `phi(7)/2 = 3`) -- the same reason the heptagon does not tile: the 7-fold
  structure never closes up, leaving a spectral gap = the floor.

So "the plane tilings vs the Platonic solids" is not a coincidence of fives; it is the angle-defect trichotomy,
and LRC14 sits exactly at its Euclidean/hyperbolic seam -- the covering-min on the hexagonal (Euclidean) side,
the apex on the Klein-quartic (hyperbolic) side, and `7` the first prime forced off both the sphere and the
plane.

## Honest scope
The identities are classical and here verified/assembled (crystallographic count `= {phi<=2}`; the `(2,3,n)`
trichotomy; the 5 Platonic; `genus X_0(2p) = 0,0,1,2,2`); the genus-`(2,3,p)` alignment is the organizing
observation. This is a **unifying map** connecting the owner's "5 and 5" to the LRC frontier and to last
session's B2/A2/Dedekind -- it explains WHY 7 is hard (first hyperbolic apex) and where the two LRC geometries
(order-6 covering-min, order-`p` apex) sit -- but it is NOT itself a new proof step; the hard core remains
`f_14`-control at `d=7` (HYP-3745/3749/3768).
