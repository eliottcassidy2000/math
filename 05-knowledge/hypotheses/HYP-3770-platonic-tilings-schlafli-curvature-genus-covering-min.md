---
id: HYP-3770
title: PLATONIC SOLIDS / PLANE TILINGS / BRAVAIS on the LRC -- the counts 5 (Platonic) and 3 (regular plane tilings) are the integer solutions of ONE Schlafli inequality (p-2)(q-2)<4 (=4), with combinatorial curvature kappa(p,q)=1/p+1/q-1/2; and this CURVATURE TRICHOTOMY (sphere chi=2 g=0 / plane chi=0 g=1 / hyperbolic chi<0 g>=2) IS the X_0(2p) GENUS trichotomy of the LRC family (genus 0,0,1,2,2 for p=3,5,7,11,13 = n=6,10,14,22,26). SO n=14 (LRC target) = genus 1 = the FLAT / EUCLIDEAN / PLANE-TILING case, whose lattice is the HEXAGONAL A_2 = the covering-min's home (Phi_6=Eisenstein, HYP-3715/3706): the Platonic solids are the genus-0 small-n analogues, hyperbolic the genus>=2 large-n. THREE proof-bridges: (1) covering-min = KERSHNER's thinnest hexagonal covering (the plane-covering theorem, HYP-3706); (2) DUALITY (p,q)<->(q,p) = the antipodal iota (THM-584/HYP-3767): self-dual {3,3},{4,4}=iota-FIXED, dual pairs octa/cube, icosa/dodeca, triangular/hexagonal = the SC-spine / complement structure; (3) GAUSS-BONNET: total combinatorial curvature = chi = the danger-cover NERVE Euler characteristic (codex M102) = the genus obstruction (S56) = the iota-odd degree (S55). CRYSTALLOGRAPHIC RESTRICTION (2D orders {1,2,3,4,6}) explains why the covering-min is 6-fold (Phi_6) not 5-fold; 5-fold (golden/icosahedral) is lattice-forbidden -> Fibonacci/quasicrystal thread
status: SYNTHESIS / exploration (rigorous classification + conceptual bridges, not a new proof). PROVED/standard: the Schlafli counts (5 Platonic, 3 plane) from (p-2)(q-2) vs 4; curvature=1/p+1/q-1/2; chi=2-2g; crystallographic restriction {1,2,3,4,6}; Kershner (hexagon = thinnest plane covering). VERIFIED (computed): the count check, the duality self-dual/pair split, the curvature=genus trichotomy table (n=14 <-> genus 1 <-> flat/hexagonal). CONJECTURAL/bridge: the curvature<->genus correspondence as a LRC organizing principle; Gauss-Bonnet chi = danger-nerve Euler char = the coverage-crux obstruction. HONEST count clarification: the 5 Platonic are the 5 regular SPHERICAL tilings; there are exactly 3 regular PLANE tilings and 5 planar BRAVAIS lattices -- the '5 plane tilings <-> 5 Platonic' folklore conflates spherical with plane; the rigorous content is the shared {p,q}/curvature classification.
source: klein-2026-06-30-S57
depends_on:
  - HYP-3587   # genus X_0(2p) = the obstruction dimension (the genus trichotomy)
  - HYP-3715   # covering-min = hexagonal A_2 / Phi_6 = Eisenstein (the flat-case lattice)
related:
  - HYP-3706   # Kershner thinnest hexagonal covering (proof-bridge 1)
  - HYP-3767   # antipodal iota / sign theory (duality = iota; Gauss-Bonnet chi = iota-odd degree)
  - HYP-3768   # covering-min = E2/Eisenstein bulk, residual = genus cusp form (the curvature=genus refinement)
  - THM-584    # complement = antipodal (duality)
  - HYP-3726   # the hexagonal margin number
---

# HYP-3770 — Platonic solids, plane tilings, and the LRC

## Where the counts come from (one Diophantine inequality)
A regular `{p,q}` (p-gons, q at each vertex) has **combinatorial curvature** `kappa = 1/p + 1/q - 1/2`,
and `sign(kappa) = sign(4 - (p-2)(q-2))`:
```
  (p-2)(q-2) < 4   kappa>0   SPHERE      chi=2   genus 0   -> the 5 PLATONIC solids
  (p-2)(q-2) = 4   kappa=0   PLANE       chi=0   genus 1   -> the 3 regular EUCLIDEAN tilings
  (p-2)(q-2) > 4   kappa<0   HYPERBOLIC  chi<0   genus>=2  -> infinitely many
```
- **5 Platonic** = `{3,3},{3,4},{4,3},{3,5},{5,3}` (tetra, octa, cube, icosa, dodeca) = the 5 regular
  SPHERICAL tilings.
- **3 plane** = `{3,6},{4,4},{6,3}` (triangular, square, HEXAGONAL).
The counts `5` and `3` are the integer points of a single inequality. **Honest clarification:** the
"5 plane tilings <-> 5 Platonic" folklore conflates *spherical* (5) with *plane* (there are exactly 3
regular plane tilings; the "5" in 2D is the 5 planar **Bravais lattices**: oblique, rectangular,
rhombic, square, hexagonal). The rigorous content is the shared `{p,q}` / curvature classification.

## Duality = the antipodal iota (sign theory, S55/S56)
`{p,q} <-> {q,p}` is an involution -- the *same* antipodal `iota` = complement `R` (THM-584, HYP-3767):
- **self-dual** `{3,3}` (tetrahedron), `{4,4}` (square tiling) = `iota`-FIXED;
- **dual pairs** octahedron/cube, icosahedron/dodecahedron, triangular/hexagonal.
This mirrors the tournament **SC-spine** (self-complementary = `iota`-fixed) + complement pairs, and the
covering-min's `iota` structure. Platonic duality IS the sign involution of the LRC/tournament program.

## Curvature = Genus (the S56 bridge): LRC-14 is the FLAT case
The tiling trichotomy `chi = 2/0/<0` (`g = 0/1/>=2`) **is** the genus trichotomy of `X_0(2p)`
(HYP-3587/3041): `genus = 0,0,1,2,2` for `p = 3,5,7,11,13` (`n = 6,10,14,22,26`). Hence
```
  n=6,10  genus 0  SPHERE      <-> Platonic (small-n analogues)
  n=14    genus 1  PLANE (FLAT)<-> HEXAGONAL A_2 = the covering-min (Phi_6=Eisenstein)   *** LRC target
  n=22,26 genus>=2 HYPERBOLIC  <-> hyperbolic {p,q} (large-n)
```
**LRC-14 is the Euclidean / plane-tiling case.** Its lattice is the hexagonal `A_2` (`Phi_6 = Eisenstein`,
HYP-3715/3706) -- which is exactly *why* the covering-min is `n/Phi_6` and hexagonal, and why `n=14` is
the borderline-tractable LRC (the flat genus-1 boundary, the "last doublet-tractable case" HYP-3587,
elliptic curve `14a` = the genus-1 cusp form of HYP-3768).

## Three proof-bridges
1. **Kershner.** The covering-min `= n/Phi_6` is the LRC instance of **Kershner's theorem** (the regular
   hexagon is the thinnest covering of the plane, Fejes Toth): the covering-min lower bound is a
   plane-covering-density statement on the hexagonal `A_2` lattice (HYP-3706).
2. **Duality = iota** (above): the proof's `iota`-equivariant structure (S55) is Platonic duality; the
   self-dual cases are the `iota`-fixed / self-complementary spine.
3. **Gauss-Bonnet = the coverage crux.** Total combinatorial curvature `= chi` (Gauss-Bonnet). The
   danger-cover **NERVE**'s Euler characteristic (codex atlas M102, "danger-cover nerve / Euler hole
   certificate") is the coverage obstruction; `chi = the genus (S56) = the iota-odd degree (S55)`. So the
   crux "the danger sheaf covers" is a **Gauss-Bonnet / Euler-characteristic** statement -- a lonely
   point is a curvature/genus defect. This unifies S55 (iota-odd degree), S56 (genus cusp form), and the
   Platonic curvature into ONE invariant `chi`.

## Crystallographic restriction: 6-fold, not 5-fold
2D lattices admit only rotation orders `{1,2,3,4,6}` (`2cos(2pi/k) in Z`) -- **5-fold is FORBIDDEN**. So
the covering-min is **6-fold** (hexagonal, `Phi_6`, the densest lattice), never 5-fold. The 5-fold
symmetry (icosahedron/dodecahedron, the golden ratio) is lattice-illegal and surfaces instead in the
project's **Fibonacci/Zeckendorf/quasicrystal** thread (the CF/Farey ladder of the covering-min rung).
The "5" of the Platonic solids and the "6" of the covering-min are the sphere-vs-plane faces of the
curvature classification.

## Net
The 5 Platonic solids (spherical `{p,q}`), the 3 regular plane tilings (Euclidean `{p,q}`), and the 5
Bravais lattices all come from one Schlafli/curvature classification; that curvature trichotomy is the
`X_0(2p)` genus trichotomy, placing **LRC-14 at the flat / hexagonal / genus-1 boundary** -- the reason
the covering-min is `n/Phi_6` (hexagonal `A_2`, Kershner) and the reason `n=14` is borderline-tractable.
Platonic **duality** is the antipodal `iota` (sign theory), and **Gauss-Bonnet** identifies the coverage
crux's obstruction with `chi` = the genus = the `iota`-odd degree, tying S55/S56 to the Platonic
curvature. Synthesis, not a new proof -- but it locates the whole LRC-14 program at curvature zero.
