import cmath, math
print("=== the n-1 <-> n <-> n+1 mapping: degree, roots, coefficients (FTA + Vieta), across the arc ===")
print("  a degree-n polynomial has: n+1 COEFFICIENTS (incl the constant) -> n ROOTS in C (FTA, with multiplicity).")
print("  monic: n free coefficients = the n elementary symmetric functions of the n roots (VIETA = a bijection).")
print("  the constant term (the +1th coefficient) = (-1)^n * (product of roots) = the BASE/vacuum/apex.")
print("  the DERIVATIVE: n roots -> n-1 critical points (Gauss-Lucas: in their convex hull) = the n -> n-1 map.\n")
rows=[
("FTA / Vieta", "n+1 coefficients <-> n roots; constant term = product of roots (the +1 = the base)"),
("Hamiltonian path", "n VERTICES, n-1 EDGES (a path) -- H counts them; the n <-> n-1 map (the spanning tree)"),
("LRC", "n runners, gap delta = 1/(n+1); the n <-> n+1 (the stationary runner / the origin = the +1)"),
("perspective key", "persp(n) = #structures(n+1) for small n -- the n <-> n+1 lift"),
("independence poly / H", "degree = #objects (n), n+1 coefficients alpha_0..alpha_n, alpha_0=1 (empty set = the base); roots = resonance/Lee-Yang/CM spectrum"),
("cyclotomic Phi_n", "degree phi(n), roots = primitive n-th roots of unity; Phi_3=x^2+x+1, 2 roots = w,w^2 (the cube root)"),
("derivative / Gauss-Lucas", "n roots -> n-1 critical points (in the convex hull) = the n -> n-1 contraction"),
("Euler characteristic", "V - E + F: the alternating n,n-1,... = the deletion-contraction / Tutte (S633)"),
]
for a,b in rows: print(f"  {a:22}: {b}")

print("\n=== the cube root: x^3 - 1 has 3 roots (FTA, degree 3 = 3 solutions) = the AG_n / resonance eigenvalues ===")
for k in range(3):
    w=cmath.exp(2j*cmath.pi*k/3)
    print(f"  cube root {k}: e^(2pi i*{k}/3) = {w.real:+.3f}{w.imag:+.3f}i")
roots=[cmath.exp(2j*cmath.pi*k/3) for k in range(3)]
print(f"  sum of cube roots = {sum(roots).real:+.3f} (=0, the x^2 coefficient: n roots sum to the trivial = Vieta)")
print(f"  product of cube roots = {(roots[0]*roots[1]*roots[2]).real:+.3f} (=1 = -(-1) = constant term: the +1/base)")
print("  x^3-1 = (x-1)(x^2+x+1) = (x-1)*Phi_3; the 3-cycle's char poly; eigenvalues 1,w,w^2 = the AG_n generators (S635).")

print("\n=== the deep principle: the +1 (constant term / vacuum / base / apex / stationary runner) vs the n nontrivial ===")
print("  the n+1th coefficient (constant) = the GROUND STATE: alpha_0=1 (empty set), the lonely origin, the apex, the identity.")
print("  the n roots/objects = the EXCITATIONS: the resonances, the runners, the cycles, the non-identity group elements.")
print("  FTA = 'the excitation count matches the degree'; Vieta = 'the symmetric functions of excitations = the coefficients';")
print("  the partition function (indep poly / H / chromatic / Tutte) IS this polynomial; coefficients<->roots = counts<->spectrum.")
