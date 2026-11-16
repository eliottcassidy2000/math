"""Classify digit-0 vertices t (B(T,t) = H(T)) of strong tournaments of order <= 7:
out-degree of t, and the shape of T - t (its last strong component, and whether t's out-neighbours sit before it).
Near-join tip: out-degree 1 and T - t = A => (transitive tail of k >= 1 singletons) with t0 in A."""
import sys
sys.path.insert(0, '.')
from w52_digit0 import classes, strong, B, near_join_tip
from w45_hk import ham_split, sub
def comps(n, out, verts):
    """strong components of the subtournament on verts, in dominance order (first = top)"""
    vs = list(verts)
    def reach(v):
        r = {v}; st = [v]
        while st:
            u = st.pop()
            for w in vs:
                if w not in r and out[u] >> w & 1: r.add(w); st.append(w)
        return r
    R = {v: reach(v) for v in vs}
    cs = []; seen = set()
    for v in vs:
        if v in seen: continue
        c = {w for w in vs if v in R[w] and w in R[v]}; seen |= c; cs.append(c)
    cs.sort(key=lambda c: -len(R[next(iter(c))]))     # top component reaches the most
    return cs
rows = []
for n in range(3, 8):
    for out in classes(n):
        if not strong(n, out): continue
        H = ham_split(n, out)[0]
        for t in range(n):
            if B(n, out, t) != H: continue
            d = bin(out[t]).count('1'); tip = near_join_tip(n, out, t)
            cs = comps(n, out, [v for v in range(n) if v != t])
            rows.append((n, H, d, len(cs[-1]), tip))
from collections import Counter
c = Counter((d, last, tip) for (n, H, d, last, tip) in rows)
print("digit-0 vertices of strong tournaments, n <= 7, by (out-degree, size of last component of T - t, tip?):")
for k, v in sorted(c.items()): print("  ", k, v)
print("non-tip digit-0 counts:", sorted(set(H for (n, H, d, last, tip) in rows if not tip)))
print("all digit-0 counts:", sorted(set(H for (n, H, d, last, tip) in rows)))
