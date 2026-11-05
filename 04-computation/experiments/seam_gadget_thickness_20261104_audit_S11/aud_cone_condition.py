"""aud_cone_condition.py -- audit S11: the rigorous counting form of the HYP-9471 "necessary density condition".

Let C = {(H, R): t0 <= R/H <= t1} and let a finite family of plane maps A_i = [[a, b], [c, d]] (H' = aH + bR,
R' = cH + dR, a, b >= 1) self-cover the lattice points (H odd, R even) of C beyond some height: every such point is
A_i q for some lattice point q of C.  Counting points with H <= X (area / 4 + O(X)) gives, as X -> infinity,

    1 <= (1 / (t1 - t0)) * sum_i  integral over {t in [t0, t1] : A_i . t in [t0, t1]} of (a_i + b_i t)^(-2) dt,

where A.t = (c + d t) / (a + b t) is the slope action (monotone on t >= 0).  Without the restriction A_i.t in C the
right side is the slope-uniform average of D(t) over [t0, t1]; the restriction only lowers it.
Prints both sides for the order <= 8 alphabet (distinct nonsingular restrictions of S-row-zero matrices) on a grid of
cones, and the singular fraction of the gadget matrices by order (note section 2 says "most" are singular).
Usage: python3 aud_cone_condition.py DIR KMAX"""
import sys
import numpy as np
D = sys.argv[1]; K = int(sys.argv[2])
maps = set(); sing = []
for k in range(1, K + 1):
    ns = tot = 0
    for line in open(f'{D}/aud_mats{k}.txt'):
        m = list(map(int, line.split())); tot += 1
        det = m[0]*(m[4]*m[8]-m[5]*m[7]) - m[1]*(m[3]*m[8]-m[5]*m[6]) + m[2]*(m[3]*m[7]-m[4]*m[6])
        ns += det == 0
        if m[3] == m[4] == m[5] == 0 and m[0]*m[8] - m[2]*m[6] != 0:
            maps.add((m[0], m[2], m[6], m[8]))
    sing.append((k, ns, tot))
print("singular gadget matrices by exact order (singular / distinct):", ", ".join(f"{k}: {s}/{t} ({s/t:.1%})" for k, s, t in sing))
A = np.array(sorted(maps), dtype=float); a, b, c, d = A.T
assert (a >= 1).all() and (b >= 1).all()
print(f"plane maps (order <= {K}): {len(A)}")
def prim(lo, hi):   # integral of (a + b t)^-2 over [lo, hi] (vectorised, lo <= hi)
    return (1.0 / (a + b * lo) - 1.0 / (a + b * hi)) / b
def slope(t): return (c + d * t) / (a + b * t)
best = None
rows = []
for t0 in (2.0, 2.25, 2.5, 3.0, 4.0):
    for t1 in (t0 + 0.05, t0 + 0.25, t0 + 0.5, t0 + 1.0, t0 + 2.0, t0 + 4.0, 12.0, 50.0):
        if t1 <= t0: continue
        full = prim(t0, t1).sum() / (t1 - t0)
        # sub-interval of [t0, t1] mapped into [t0, t1]; slope map is monotone: solve on a fine grid robustly
        lo = np.full(len(A), np.nan); hi = np.full(len(A), np.nan)
        grid = np.linspace(t0, t1, 4001)
        inside = np.array([(slope(t) >= t0 - 1e-12) & (slope(t) <= t1 + 1e-12) for t in grid])   # shape (grid, maps)
        anyin = inside.any(axis=0)
        idx_first = np.argmax(inside, axis=0); idx_last = len(grid) - 1 - np.argmax(inside[::-1], axis=0)
        lo[anyin] = grid[idx_first[anyin]]; hi[anyin] = grid[idx_last[anyin]]
        restr = np.where(anyin, prim(np.where(anyin, lo, t0), np.where(anyin, hi, t0)), 0.0).sum() / (t1 - t0)
        rows.append((t0, t1, full, restr))
        if best is None or restr > best[3]: best = (t0, t1, full, restr)
print("cone [t0, t1] | slope-average of D_8 (all maps) | rigorous necessary-condition value (images kept in the cone)")
for t0, t1, full, restr in rows:
    print(f"   [{t0:.2f}, {t1:.2f}] | {full:.4f} | {restr:.4f}")
print(f"largest necessary-condition value on this grid: {best[3]:.4f} at [{best[0]}, {best[1]}]")
