"""Li-bigraded Hilbert series of J_infty(R_q): dims by (weight w, Li degree l), l = number of derivatives."""
import sys
import flint
from arc_hs import monos, mul, deriv, poly, RELS, WT, P
def bigraded(rels, W):
    gens = [poly(r) for r in rels]
    gw = [sum(WT[v] + j for (v, j) in next(iter(g))) for g in gens]
    ders = []
    for g, w in zip(gens, gw):
        lst = [g]
        while w + len(lst) - 1 < W:
            lst.append(deriv(lst[-1]))
        ders.append((w, lst))
    li = lambda m: sum(j for (v, j) in m)
    table = {}
    for w in range(0, W + 1):
        M = monos(w)
        for l in sorted({li(m) for m in M}):
            Ml = [m for m in M if li(m) == l]; idx = {m: i for i, m in enumerate(Ml)}
            rows = []
            for (gw0, lst) in ders:
                for n, dg in enumerate(lst):
                    rest = w - gw0 - n
                    if rest < 0: break
                    for m in monos(rest):
                        if li(m) + n != l: continue
                        row = {}
                        for mo, co in dg.items():
                            k = idx[mul(mo, m)]; row[k] = (row.get(k, 0) + co) % P
                        rows.append(row)
            rk = flint.nmod_mat(len(rows), len(Ml), [r.get(j, 0) for r in rows for j in range(len(Ml))], P).rank() if rows else 0
            table[(w, l)] = len(Ml) - rk
    return table
q, W = int(sys.argv[1]), int(sys.argv[2])
T = bigraded(RELS[q], W)
Lmax = max(l for (w, l) in T)
print(f"q={q}: rows = Li degree l, columns = weight w (0..{W})")
for l in range(0, Lmax + 1):
    print(f"l={l:2d}:", " ".join(f"{T.get((w, l), 0):3d}" for w in range(W + 1)))
