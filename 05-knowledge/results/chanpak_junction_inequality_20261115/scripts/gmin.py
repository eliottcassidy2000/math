"""g_m = least H(T[u]) over factors u of m vertices whose operator has rank >= 2 (rank mod p), for small m,
and the longest such factors with H(T[u]) <= G.  usage: python3 gmin.py G MAXDEPTH CAP"""
import pickle, sys, os, numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'badtree.py')).read().split("# level 1:")[0].replace("RHO = float(sys.argv[1]); CAP = int(sys.argv[2]); MAXD = int(sys.argv[3]) if len(sys.argv) > 3 else 200", "G = float(sys.argv[1]); MAXD = int(sys.argv[2]); CAP = int(sys.argv[3])"))
level = [((l,), col_basis(Mm[l]), v0f.copy()) for l in range(8)]
for depth in range(1, MAXD + 1):
    keep = [(w, L, x) for w, L, x in level if L.shape[1] >= 2 and hf @ x <= G + 1e-9]
    hs = sorted({int(round(hf @ x)) for w, L, x in keep})
    print("m = %2d: %7d factors of rank >= 2 with H <= %g ; least H among them: %s ; example %s" %
          (depth, len(keep), G, hs[0] if hs else '-', keep[0][0] if keep else ''), flush=True)
    if not keep or len(keep) > CAP: break
    nxt = []
    for w, L, x in keep:
        m = depth + 1
        for l in range(8):
            le = l & ((1 << min(3, m - 1)) - 1)
            nxt.append((w + (l,), col_basis(modmul(Mm[l], L)), Af[le] @ x))
    level = nxt
