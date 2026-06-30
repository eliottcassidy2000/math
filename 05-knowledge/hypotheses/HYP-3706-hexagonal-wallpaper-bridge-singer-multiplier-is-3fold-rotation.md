---
id: HYP-3706
title: THE DESIGN-COVERING BRIDGE via the HEXAGONAL WALLPAPER GROUP -- Q(sqrt-3)=Z[zeta_6]=the hexagonal lattice (wallpaper p6m); Phi_6(n)=norm=index, Z/Phi_6(n)=hexagonal-lattice quotient; the rotation tower in the quotient is zeta_6=mult-by-n (order 6, the ambient 60-deg), zeta_3=zeta_6^2=mult-by-n^2=mult-by-(n-1)=q=the SINGER MULTIPLIER (order 3, the 120-deg = the Frobenius of GF(q^3)/GF(q)), -1=n^3 (order 2); the Singer difference set carries p3 (3-FOLD) symmetry exactly (mult-by-q IS a multiplier; mult-by-n is NOT) -- a union of 3-fold orbits (n=4: {0}u{1,3,9}); the 3-fold zeta_3 is the Condorcet/3-cycle atom. Bridge: design-optimality (projective plane = optimal Steiner covering, p3-symmetric Singer set) <-> continuous covering-optimality (hexagonal lattice = THINNEST plane covering, Kershner 1939, p6m), transferred by the shared Eisenstein/hexagonal symmetry; covering-min n/Phi_6(n) = the hexagonal-symmetric covering density
status: STRUCTURE VERIFIED (n order 6 / q=n-1=n^2 order 3 / -1 order 2 mod Phi_6(n), n=3..14; mult-by-q a multiplier of the (v,n,1) Singer set, mult-by-n NOT, n=3,4,6). The wallpaper identification (p6m lattice, p3 difference set, multiplier=3-fold=Frobenius) is solid; the LRC-covering = hexagonal-covering claim (the actual bridge) + Kershner-transfer remain OPEN.
source: klein-2026-06-29-S24
depends_on:
  - HYP-3705   # n/Phi_6(n) covering-min = the projective-plane optimality (this bridges it)
related:
  - HYP-3700   # mac-mini: the apex column (the complementary Q(cos2pi/p) side)
  - HYP-3610   # the multi-axis poles (difference-set = spread pole = this construction)
  - HYP-3602   # the 3-cycle / Condorcet (= the 3-fold zeta_3 multiplier here)
results:
  - 04-computation/hexagonal_wallpaper_bridge_klein.py
  - 05-knowledge/results/hexagonal_wallpaper_bridge_klein.out
---

# HYP-3706 — the design-covering bridge is the hexagonal wallpaper group; the Singer multiplier is the 3-fold rotation

## The hexagonal lattice and its rotation tower (verified)
`Q(sqrt-3) = Q(zeta_6) = Q(zeta_3)` and `Z[zeta_6]` = the **hexagonal (triangular) lattice**, wallpaper
group **p6m**. `Phi_6(n) = |n - zeta_6|^2` = the Eisenstein norm = the index `[Z[zeta_6] : (n-zeta_6)]`, and
`Z[zeta_6]/(n-zeta_6) = Z/Phi_6(n)` (cyclic). In the quotient `zeta_6 = n`, so the hexagonal point group
`C6` acts on `Z/Phi_6(n)` as multiplication:
```
zeta_6  = mult-by-n        order 6   (60-deg rotation; the AMBIENT hexagonal symmetry)
zeta_3  = mult-by-n^2 = mult-by-(n-1) = q   order 3   (120-deg; = the SINGER MULTIPLIER = Frobenius)
-1      = mult-by-n^3      order 2   (180-deg)
```
Verified `n=3,4,6,8,14,...`: `ord(n)=6`, `ord(q)=ord(n-1)=3`, `ord(-1)=2`. The identity `q = n-1 = n^2 mod
Phi_6(n)` (since `n^2 - (n^2-n+1) = n-1`) is the key: **the classical Singer multiplier `q` IS the 3-fold
rotation `zeta_3`**, the Frobenius of `GF(q^3)/GF(q)` (Galois group `C3`).

## The Singer difference set has p3 (not p6) symmetry
For the `(Phi_6(n), n, 1)` Singer difference set `D`:
- **mult-by-`q` = the 120-deg rotation IS a multiplier** (`qD = D + s`): verified `n=3` (shift 1), `n=4`
  (shift 0), `n=6` (shift 3). A centered translate of `D` is a union of **3-fold orbits** (e.g. `n=4`:
  `D={0,1,3,9} = {0} u {1,3,9}`, the orbit of `1` under `x -> 3x` mod 13).
- **mult-by-`n` = the 60-deg rotation is NOT a multiplier** (no shift works, all cases).
So the difference set carries the **`p3`** wallpaper symmetry (3-fold `zeta_3` = Frobenius = Singer
multiplier), embedded in the ambient **`p6m`** hexagonal lattice. The full 6-fold is the lattice's, not the
design's. (The 3-fold `zeta_3` is exactly the project's Condorcet / odd-3-cycle atom, HYP-3602 -- the
minimal intransitivity is the difference set's rotational symmetry.)

## The bridge
- **Discrete side:** the projective plane `PG(2,n-1)` / the `(v,n,1)` Singer difference set is the OPTIMAL
  covering DESIGN (a Steiner system `S(2,n,v)` attains the covering number exactly), carrying `p3` symmetry.
- **Continuous side:** the hexagonal lattice is the **THINNEST covering of the plane** (Kershner 1939) -- the
  unique optimal LATTICE covering, with full `p6m` symmetry.
- **The transfer:** both live on `Z[zeta_6]` (the same hexagonal lattice / Eisenstein field). The covering-min
  `n/Phi_6(n)` is the density of the hexagonal-symmetric covering; the discrete design-optimality and the
  continuous covering-optimality are two faces of the hexagonal lattice's extremality, linked by the shared
  `p6m`/`p3` symmetry. The wallpaper group is the bridge: it is the symmetry under which BOTH the discrete
  Steiner optimality and the continuous Kershner optimality are the extremal configuration.

## Honest scope
- RIGOROUS: the hexagonal-lattice / Eisenstein identification; the rotation tower (`n` order 6, `q=n-1=n^2`
  order 3 = the Singer multiplier = Frobenius, `-1` order 2); the difference set's `p3` symmetry (`q` a
  multiplier, `n` not); Kershner (hexagonal = thinnest covering) and Steiner optimality are known theorems.
- OPEN (the actual bridge): that the LRC continuous covering-min IS achieved by the `p6m`-hexagonal covering
  (so its optimality = Kershner), i.e. the LRC floor equals the design/hexagonal covering number. This needs
  M's exact LRC definition (floor owners') and the LRC -> hexagonal-torus reduction. The wallpaper frame
  SHARPENS the open bridge to a clean geometric claim: "the optimal LRC covering is the hexagonal one."

## Net
The design-covering bridge runs through the hexagonal wallpaper group `p6m`: `Phi_6` is the Eisenstein/
hexagonal norm, the Singer multiplier is the 3-fold rotation (`zeta_3` = Frobenius = the Condorcet 3-cycle),
the difference set is `p3`-symmetric, and the optimality transfers between the discrete Steiner design and
the continuous Kershner covering because both are the hexagonal lattice's extremal configuration. Open: that
the LRC covering-min is that hexagonal covering.
