---
id: HYP-3704
title: THREE COORDINATED ROUTES to the zeta_6-line covering radius (proving the construction {1,..,n-2,n(n-1)} = the cyclic-Kershner covering-min, M=n/Phi_6(n), for n>=7) -- (1) TORUS LIFT: Z[omega]/(n-omega) ~ Z/Phi_6(n) with omega (the 60deg hexagonal rotation) = mult-by-n (verified n^2=n-1, n^3=-1, ord(n)=6 mod Phi_6); the construction's residues v*n mod Phi_6 are a hexagonal AP, covering radius n/Phi_6, reducing to KERSHNER 1939 (hexagonal = thinnest lattice covering, PROVEN) for lattice/AP coverings; the binding pair at t*=n/Phi_6 is {1, n(n-1)} (smallest & outlier); (2) THREE-DISTANCE: the AP residues have gaps EXACTLY {1,n,2n} (the missing origin k=0 = the -1-mod-Phi_6 'lcm-killer' gives gap 1; the origin-straddle gives 2n), so M = min|residue|/Phi_6 = n/Phi_6, computable in closed form for any {1}+AP, and minimizing the max-gap among APs picks the construction; (3) SPECTRAL/LP: the covering radius = a Delsarte/Cohn-Elkies LP (sup over lonely t = inf over Fourier-positive certificates vanishing on the danger frequencies), with R-hat = the AP's Dirichlet/Gauss-sum; the EISENSTEIN-SYMMETRIC (omega=xn-invariant) optimal certificate function is the open node (a Viazovska-style Fourier-extremal problem) and handles ALL coverings, not just lattices. NET: the three reduce the cyclic-Kershner to ONE Eisenstein-symmetric Fourier-positivity certificate
status: COMPUTED + verified (the lift iso n^2/n^3/ord; the three-distance gaps {1,n,2n} n=7,14; the binding pair {1,n(n-1)}; R-hat AP/Dirichlet structure). Routes 1,2 are concrete and reduce APs/lattices to Kershner (proven); route 3 (the spectral/LP for ALL coverings) is framed (Delsarte/Cohn-Elkies) with the Eisenstein-symmetric extremal function the open node. Not a complete proof; a coordinated attack reducing to one certificate.
source: mac-mini-2026-06-30-S45
related:
  - HYP-3703  # the covering-min IS a tiling-optimality; the 6/7 boundary, A2/Eisenstein, Kershner
  - HYP-3706  # klein-S24: omega=xn, p6m, Singer multiplier; design<->Kershner bridge
  - HYP-3702  # the covering-set taxonomy; the {1}+AP hard core (route 2's family)
  - HYP-3606  # the least-eigenvalue/non-bipartiteness certificate (the CUSP spectral; this is the OFF-CUSP/Eisenstein spectral)
results:
  - 04-computation/three_routes_to_zeta6_covering_radius_macmini_20260630.py
  - 05-knowledge/results/three_routes_to_zeta6_covering_radius_macmini_20260630.out
---

# HYP-3704 -- three routes to the zeta_6-line covering radius

The owner asked to work the torus lift, the three-distance gaps directly, and a spectral route to the
`zeta_6`-line covering radius -- the three attacks on the cyclic-Kershner hard core (the construction
`{1,..,n-2,n(n-1)}`, `M=n/Phi_6(n)`, is the covering-min for `n>=7`).

## Route 1 -- the TORUS LIFT (the geometry)
`Phi_6(n) = N(n-omega)` in the Eisenstein integers `Z[omega]` (`omega = e^{i pi/3}`, the primitive 6th
root), and `Z[omega]/(n-omega) ~ Z/Phi_6(n)` as rings, with `omega ≡ n`. VERIFIED: `n^2 ≡ n-1`, `n^3 ≡ -1`
(mod `Phi_6`), so `ord(n)=6` -- **multiplication by `n` IS the 60-degree hexagonal rotation `omega`**. The
LRC danger configuration, scaled by `Phi_6` and lifted, becomes a configuration on this hexagonal torus;
the construction's speed-residues `v*n mod Phi_6` form a hexagonal arithmetic progression, and the covering
radius is `n/Phi_6`. The binding pair at `t*=n/Phi_6` is `{1, n(n-1)}` (the smallest speed and the outlier --
the two framing the AP). For LATTICE / AP coverings this reduces to **KERSHNER 1939** (the hexagonal lattice
is the thinnest plane covering -- PROVEN), so the construction is optimal *within lattice/AP coverings*.

## Route 2 -- the THREE-DISTANCE gaps (the explicit M)
The residues `r_v = v*n mod Phi_6` are the AP `{kn : k=1,..,n-2}` plus `-n` (`= n(n-1)*n mod Phi_6`),
i.e. `{kn : k=-1,1,2,..,n-2}` -- the full AP `{kn}` with the **origin `k=0` MISSING**. By the three-distance
(Steinhaus) theorem its gaps take exactly the three values `{1, n, 2n}`: `n` between consecutive `kn`; `1`
at the `(n-2)n -> -n` step (the `-1 ≡ n^3 mod Phi_6` "lcm-killer" that closes the AP); `2n` straddling the
missing origin. The covering radius is half the origin-straddle distance to the nearest residue, `= n`, so
`M = n/Phi_6`. This is a CLOSED FORM, and for any `{1}+AP` covering (the hard-core family, HYP-3702) the
three-distance theorem gives `M` directly; minimizing the max-gap among AP coverings selects the
construction.

## Route 3 -- the SPECTRAL / LP route (the optimality machinery)
The covering radius `= sup{r : exists a lonely t}` is a linear program; by Delsarte/Cohn-Elkies duality it
equals `inf` over **Fourier-positive certificate functions** `F >= 0` whose Fourier support lies on the
danger frequencies (the speeds and their multiples) and which is positive at a lonely `t`. The construction's
residue Fourier transform `R-hat(j) = sum_{r in R} omega^{jr}` is the AP's **Dirichlet/Gauss sum** (a
geometric series, peaked), and the max gap `2n` is read off it (Beurling-Selberg). The optimal certificate
is **EISENSTEIN-SYMMETRIC** -- invariant under `omega = mult-by-n` -- a Viazovska-style Fourier-extremal
function on the hexagonal torus. Unlike routes 1-2, this route handles **ALL coverings**, not just
lattices/APs; constructing the `omega`-symmetric optimal `F` is the precise open node. (This is the OFF-CUSP
/ Eisenstein analogue of the CUSP least-eigenvalue certificate HYP-3606, which was `Z_7`-cyclotomic; here it
is `Z[omega]`-Eisenstein.)

## The synthesis -- one certificate
The three routes are a coordinated reduction:
- **lift** gives the GEOMETRY (the hexagonal torus, `omega=xn`) and, with Kershner, optimality among
  lattice/AP coverings;
- **three-distance** gives the EXPLICIT `M = n/Phi_6` and the closed form for AP coverings;
- **spectral/LP** gives the OPTIMALITY MACHINERY for all coverings, reducing to one `omega`-symmetric
  Fourier-positive certificate.
So proving the construction is the cyclic-Kershner covering-min reduces to constructing ONE
Eisenstein-symmetric Fourier-positivity certificate (Cohn-Elkies/Viazovska-style) on `Z[omega]/(n-omega)` --
the single, modern, well-posed open node. Kershner (1939, lattice) is the proven backstop; the certificate
extends it to the general (non-lattice) covering.

## Convergence with klein-S27 (HYP-3716): the 2D side is a THEOREM; the spectral route is a Jacobi operator
klein-S27 independently closed the **2D covering side**: the rhombic-lattice covering density
`delta(theta) = pi R^2 / sin(theta)` is minimized EXACTLY at `theta=60deg` (`= 2pi/sqrt27` = Kershner; the
rhombus = 2 equilateral triangles = hexagonal), beating `theta=90deg` (the square, `pi/2`); the Voronoi cell
of a lattice covering is a centrally-symmetric convex `<=6`-gon (hexagon or rectangle, a FINITE check,
hexagon wins -- Reinhardt-Rao); and the Socolar-Taylor decorated hexagon gives NO covering gain (Kershner is
for all coverings). So **'the optimal plane covering is hexagonal' is a THEOREM** -- exactly route 1's
backstop, now fully closed. The remaining gap is precisely the **LRC -> 2D reduction** (does the `zeta_6`-line
in `Z/Phi_6(n)` inherit the 2D hexagonal optimality) = my torus lift. And klein's SPECTRAL lead makes route 3
concrete: **the A2/triangular lattice 'tridiagonalized' = the JACOBI operator whose moments are the Catalan
family** (HYP-3710); its spectrum is the semicircle (free-probability), and **the `zeta_6`-line covering
radius = a spectral EDGE of this A2-Jacobi operator**. So route 3's "Eisenstein-symmetric certificate" has a
concrete spectral form: the edge of the Catalan-moment (semicircle) Jacobi operator. The two accounts
dovetail exactly: 2D side = theorem (klein), the lift + the Jacobi-edge spectral route = the open node.

## What it buys
Turns the covering-min hard core into three concrete, mutually-reinforcing attacks, and pins the open node
to a single object: an `omega = mult-by-n` (60-degree hexagonal) symmetric Fourier-positive certificate on
the Eisenstein torus. The torus lift and three-distance are done (APs/lattices, via Kershner); the spectral
LP certificate is the one remaining piece, and it is exactly the off-cusp/Eisenstein twin of the proven
apex-cusp non-bipartiteness certificate (HYP-3606).
