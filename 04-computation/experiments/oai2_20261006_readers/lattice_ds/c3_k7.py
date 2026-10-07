# FINITE-EXACT: the 7-point Eisenstein torus Z[w]/(2+w), w = e^{i pi/3} (w^2 = w - 1)
# and the 21-point torus Z[w]/(4+w): minimal-distance graphs, 120-degree orientation, faces, symmetries.
from itertools import permutations, combinations
def ring_map(a, b, N):
    # Z[w]/(a+b w) -> Z/N via w -> root of x^2 - x + 1 with a + b*x = 0 mod N (b=1 => x=-a)
    assert b == 1
    x = (-a) % N
    assert (x*x - x + 1) % N == 0
    return x
for (a, N) in ((2, 7), (4, 21), (3, 13)):
    x = ring_map(a, 1, N)
    units = [1, x, (x*x) % N, (x**3) % N, (x**4) % N, (x**5) % N]   # w^k, k=0..5
    print(f"N={N}: w -> {x};  images of the six units w^k:", units)
    nb = set(units)
    edges = {frozenset((i, (i+u) % N)) for i in range(N) for u in units}
    print(f"   minimal-distance graph: circulant C_{N}{sorted(min(u, N-u) for u in nb)}; edges = {len(edges)}; degree 6")
    even = sorted({units[0], units[2], units[4]})   # 1, w^2, w^4 : the 120-degree directions
    print(f"   120-degree directions 1,w^2,w^4 -> {even}")
    # up-triangles {0,1,w}, down-triangles {0,w,w^2}
    up = sorted({0, 1, x}); dn = sorted({0, x, (x*x) % N})
    for nm, tri in (("up", up), ("down", dn)):
        diffs = sorted((p - q) % N for p in tri for q in tri if p != q)
        print(f"   {nm}-triangle base {tri}: differences {diffs}  perfect difference set: {sorted(diffs) == list(range(1, N)) }")
    if N == 7:
        QR = sorted({(k*k) % 7 for k in range(1, 7)})
        print("   QR_7 =", QR, " orientation = Paley P_7:", even == QR)
        # automorphisms of P_7 among affine maps x -> m x + t
        aut = [(m, t) for m in range(1, 7) for t in range(7) if all(((m*q) % 7) in QR for q in QR)]
        print("   affine maps preserving P_7:", len(aut), "(F_21 order 21);  multipliers:", sorted({m for m, t in aut}))
        # rotation by w (60 degrees) maps orientation to its converse:
        print("   rotation by w sends {1,2,4} to", sorted((x*q) % 7 for q in QR), "(= non-residues -> converse tournament)")
        # Fano: the 7 up-triangles as lines
        lines = [sorted({(p + t) % 7 for p in up}) for t in range(7)]
        print("   up-triangle lines:", lines)
        pairs = [frozenset(p) for L in lines for p in combinations(L, 2)]
        print("   every pair of points on exactly one line:", len(set(pairs)) == 21 and len(pairs) == 21)
    if N == 21:
        print("   multiplicative order of w^2 image:", next(k for k in range(1, 30) if pow(units[2], k, 21) == 1),
              "; subgroup <4> =", sorted({pow(4, k, 21) for k in range(3)}))
        # how many triangles (faces) = 2N; is the orientation a tournament? out-degree 3 < 10
        print("   faces (triangles):", 2*N, "; Euler V-E+F =", N - len(edges) + 2*N)
