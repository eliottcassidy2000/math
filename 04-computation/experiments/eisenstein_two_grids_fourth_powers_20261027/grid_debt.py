"""Catalytic merge debt by lattice shape (THM-4824).  A lattice is given by its integer Gram form Q(a, b) on Z^2 and its
unit (minimal-vector) set; g(x) = min over binary merge trees with unit leaves and root value x of sum beta(y,z)_+,
beta = Q(y+z) - Q(y) - Q(z).  Relaxation over all pairs in the window Q <= W (values only grow by a-priori bounds; the
window is enlarged until stable).  Compared: the square grid Z[i]; the first grid Z[sqrt-3] alone (units +-1); the
Eisenstein grid Z[w] = Z[sqrt-3] glued with its (1/2,1/2)-shift (written in the basis 1, w: Q = a^2 - ab + b^2)."""
import sys
def run(name, Q, units, W, report):
    R = 2
    while any(Q(R, b) <= W for b in range(-4*R, 4*R + 1)) or any(Q(a, R) <= W for a in range(-4*R, 4*R + 1)): R += 1
    pts = [(a, b) for a in range(-4*R, 4*R + 1) for b in range(-4*R, 4*R + 1) if Q(a, b) <= W]
    INF = 10**9; g = {p: INF for p in pts}; g[(0, 0)] = 0
    for u in units: g[u] = 0
    for rnd in range(200):
        changed = False
        fin = [p for p in pts if g[p] < INF]
        for i, y in enumerate(fin):
            qy = Q(*y); gy = g[y]
            for z in fin[i:]:
                x = (y[0] + z[0], y[1] + z[1]); qx = Q(*x)
                if qx > W: continue
                beta = qx - qy - Q(*z); c = gy + g[z] + (beta if beta > 0 else 0)
                if c < g[x]: g[x] = c; changed = True
        if not changed: break
    unreached = sum(1 for p in pts if g[p] >= INF)
    print(f"{name}: window Q <= {W}, {len(pts)} points, converged after {rnd + 1} rounds, unreachable points {unreached}")
    for x in report:
        print(f"   g{x} = {g.get(x, 'outside')}   (Q = {Q(*x)})")
    return g
W = int(sys.argv[1]) if len(sys.argv) > 1 else 200
gi = run("square grid Z[i]", lambda a, b: a*a + b*b, [(1,0),(-1,0),(0,1),(0,-1)], W, [(2,0),(3,0),(5,0),(7,0),(11,0),(3,4),(13,0)])
print("   max g over the whole window:", max(v for v in gi.values() if v < 10**9))
gr = run("first grid Z[sqrt-3] alone (units +-1)", lambda a, b: a*a + 3*b*b, [(1,0),(-1,0)], W, [(2,0),(3,0),(5,0),(11,0),(0,1),(1,1)])
print("   on the real line g(n) = n(n-1):", all(gr[(n,0)] == n*(n-1) for n in range(1, 14)))
ge = run("Eisenstein Z[w] = first grid + (1/2,1/2)-shift", lambda a, b: a*a - a*b + b*b, [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)], W, [(2,0),(3,0),(5,0),(7,0),(11,0),(13,0)])
print("   g = floor(Q/3) except at the exceptional points; real n: ", [(n, ge[(n,0)], n*n//3) for n in range(1, 15)])
