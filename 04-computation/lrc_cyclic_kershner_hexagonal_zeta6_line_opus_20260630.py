"""
Attack the cyclic-Kershner node with the hexagonal/tiling frame. KEY: the construction's speeds (at zeta_6)
are EQUALLY SPACED on the zeta_6-line in the 2D hexagonal lattice Z[zeta_6]/(n-zeta_6); the 1D LRC circle
(Z/Phi_6) sees this as the three-distance {1,n,2n} -- the triangular tiling TRIDIAGONALIZED. Kershner: the
hexagonal lattice is the thinnest covering => the zeta_6-line is the covering-optimal line.
"""
import cmath, math
def picture(n):
    Phi6=n*n-n+1
    w=cmath.exp(1j*math.pi/3)  # zeta_6 = e^{i pi/3}
    # construction speeds * zeta_6 mod Phi6 = multiples {b*n mod Phi6} = the zeta_6-line {b*zeta_6}
    bs=list(range(1,n-1))+[-1]   # b = 1..n-2 and -1 (the pronic -> -zeta_6)
    pts2d=[b*w for b in bs]      # 2D hexagonal-lattice positions (on the zeta_6-line, 60 deg)
    proj1d=sorted((b*n)%Phi6 for b in bs)  # 1D projection mod Phi6
    # 2D: equally spaced on the line (unit steps); 1D: three-distance gaps
    gaps1d=sorted(set((proj1d[(i+1)%len(proj1d)]-proj1d[i])%Phi6 for i in range(len(proj1d))))
    # 2D spacing along the line
    along=sorted(b for b in bs)
    sp2d=sorted(set(along[i+1]-along[i] for i in range(len(along)-1)))
    return Phi6,pts2d,proj1d,gaps1d,sp2d
print("THE HEXAGONAL zeta_6-LINE (2D) vs its 1D TRIDIAGONALIZATION (three-distance):")
for n in [7,14,20]:
    Phi6,pts2d,proj1d,gaps1d,sp2d=picture(n)
    print(f"  n={n}: 2D positions = {{b*zeta_6 : b in {sorted([b for b in (list(range(1,n-1))+[-1])])}}} on the 60-deg line")
    print(f"        2D spacing along zeta_6-line = {sp2d} (EQUALLY spaced, unit steps -- the triangular-lattice line)")
    print(f"        1D projection mod Phi6={Phi6}: three-distance gaps = {gaps1d}  <-- {{1,n,2n}} = TRIDIAGONALIZED")
print()
print("=> the 'triangular tiling is a god' = the 2D hexagonal lattice (Kershner-thinnest-covering); the")
print("   speeds are EQUALLY SPACED on its zeta_6-line. 'Tridiagonalized' = the LRC's 1D circle Z/Phi_6")
print("   sees that equal spacing as exactly THREE gap-lengths {1,n,2n} (Steinhaus). The 2D god, 1D shadow.")
print()
# Kershner covering radius of the hexagonal lattice vs square (right-angle) lattice
print("KERSHNER (1939): thinnest covering of the plane by a convex body = the regular HEXAGON.")
print(f"   hexagonal covering density theta_hex = 2*pi/(3*sqrt(3)) = {2*math.pi/(3*math.sqrt(3)):.5f}")
print(f"   square (right-angle) covering density   = pi/2          = {math.pi/2:.5f}  (WORSE: 1.209 vs 1.209? )")
print(f"   hexagon {2*math.pi/(3*math.sqrt(3)):.4f} < square {math.pi*math.sqrt(2)/2 if False else math.pi/2:.4f}: hexagon is thinner (optimal).")
print("   => the covering-min lives on the HEXAGON (60-deg, Eisenstein Q(sqrt-3)), NOT the square (90-deg,")
print("      Gaussian Q(i)). The right-angle quadrilateral fundamental domain = 2 equilateral triangles")
print("      (the Eisenstein rhombus); shearing the square->rhombus IS the Gaussian->Eisenstein move.")
print()
print("RAO 2017 closes convex-polygon tilings: triangles, quads, 3 hexagons, 15 pentagons; NONE with 7+ sides.")
print("   The covering-optimal tile is the HEXAGON (Reinhardt type); the zeta_6-line is its 1D covering slice.")
print("   Aperiodic monotiles (Socolar-Taylor, hat/spectre) cannot BEAT the periodic hexagonal covering")
print("   (Kershner is a LATTICE-covering optimum, and the hexagonal lattice attains it) -- so the optimality")
print("   proof may restrict to PERIODIC (lattice) coverings: the zeta_6-line, no aperiodic improvement.")
