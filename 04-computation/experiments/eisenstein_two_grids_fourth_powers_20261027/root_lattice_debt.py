"""Is the Eisenstein catalyst floor special?  Catalytic merge debt g on other glued ("checkerboard") lattices, with the
minimal vectors as unit charges and beta(y,z) = Q(y+z) - Q(y) - Q(z):
  * the centred rectangular family in the plane: points (u, v) in Z^2 with u = v (mod 2), Q = (u^2 + c v^2) [scaled];
    c = 1 is the square grid (rotated), c = 3 the Eisenstein grid, c = 2 the intermediate rhombic grid, c >= 4 collinear;
  * the root lattices A3 = D3 (FCC) and D4 (checkerboard in R^3, R^4: even coordinate sum, Q = x.x/2, roots +-e_i +-e_j).
Relaxation over all pairs inside the window Q <= W."""
import sys, itertools
def relax(points, Q, units, W):
    pts = [p for p in points if Q(p) <= W]
    INF = 10**9; g = {p: INF for p in pts}; zero = tuple(0 for _ in pts[0]); g[zero] = 0
    for u in units: g[u] = 0
    for rnd in range(100):
        ch = False; fin = [p for p in pts if g[p] < INF]
        for i, y in enumerate(fin):
            qy = Q(y); gy = g[y]
            for z in fin[i:]:
                x = tuple(a + b for a, b in zip(y, z)); qx = Q(x)
                if qx > W: continue
                be = qx - qy - Q(z); c = gy + g[z] + (be if be > 0 else 0)
                if c < g[x]: g[x] = c; ch = True
        if not ch: break
    return g, pts
if __name__ == '__main__':
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print("centred rectangular family, points u = v mod 2, Q = u^2 + c v^2")
    for c in (1, 2, 3, 4):
        pts = [(u, v) for u in range(-40, 41) for v in range(-40, 41) if (u - v) % 2 == 0]
        Q = lambda p, c=c: p[0]*p[0] + c*p[1]*p[1]
        m = min(Q(p) for p in pts if p != (0, 0)); units = [p for p in pts if Q(p) == m]
        g, P = relax(pts, Q, units, W * m)
        reach = [p for p in P if g[p] < 10**9]
        vals = sorted({(Q(p), g[p]) for p in reach if p != (0, 0)})[:14]
        ratio = max((g[p] / Q(p) for p in reach if Q(p) > 0), default=0)
        print(f"  c={c}: minimal norm {m}, {len(units)} unit charges, reachable {len(reach)}/{len(P)}, max g/Q = {ratio:.4f}; (Q, g) samples: {vals}")
    print("root lattices with Q = x.x/2 (roots have Q = 1)")
    for d in (3, 4):
        R = 5 if d == 3 else 4
        pts = [p for p in itertools.product(range(-R, R + 1), repeat=d) if sum(p) % 2 == 0]
        Q = lambda p: sum(a*a for a in p) // 2
        units = [p for p in pts if Q(p) == 1]
        g, P = relax(pts, Q, units, 8 if d == 3 else 5)
        reach = [p for p in P if g[p] < 10**9]
        print(f"  D{d}{' = A3' if d == 3 else ''}: {len(units)} roots, window Q <= {8 if d == 3 else 5}: {len(P)} points, reachable {len(reach)}, max g = {max(g[p] for p in reach)}")

