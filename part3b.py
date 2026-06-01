from fractions import Fraction as F

def all_points(speeds, t):
    pts = {F(0)}
    for v in speeds:
        pts.add((v * t) % 1)
    return sorted(pts)

def observer_gaps(speeds, t):
    p = all_points(speeds, t)
    if len(p) == 1:
        return F(1), F(0)  # degenerate: only 0
    right = p[1] - p[0]
    left = (p[0] + 1) - p[-1]
    return left, right

def lonely_gap(speeds, t):
    """min distance from 0 to any other runner = the standard LRC gap to nearest."""
    best = F(1)
    for v in speeds:
        x = (v*t) % 1
        if x == 0:
            return F(0)
        best = min(best, x, 1-x)
    return best

def walls(speeds):
    W = set()
    vs = list(speeds)
    for v in vs:
        for k in range(0, v+1):
            W.add(F(k, v))
        for k in range(0, 2*v+1):
            W.add(F(k, 2*v))
    for i in range(len(vs)):
        for j in range(i+1,len(vs)):
            for d in (abs(vs[i]-vs[j]), vs[i]+vs[j]):
                if d>0:
                    for k in range(0, d+1):
                        W.add(F(k,d))
    return sorted(w for w in W if 0 <= w <= 1)

def maximize(speeds, func, require_distinct=True):
    W = walls(speeds)
    cands = set()
    for i in range(len(W)-1):
        a,b = W[i],W[i+1]
        if b<=a: continue
        cands.add((a+b)/2)
    best=F(-1); bt=None
    for t in cands:
        if not (0<t<1): continue
        if require_distinct:
            p = all_points(speeds,t)
            if len(p) != len(speeds)+1:  # require all points distinct (non-collapse)
                continue
        val = func(speeds,t)
        if val>best:
            best=val; bt=t
    return best,bt

print("Standard Lonely-Runner gap: max_t min_i dist(0, v_i t)  (loneliness)")
print("(non-degenerate cell midpoints)\n")
for m in [4,5,6,7]:
    speeds=list(range(1,m+1))
    b,bt = maximize(speeds, lonely_gap, require_distinct=False)
    print(f"  AP m={m}: max loneliness = {b} = {float(b):.5f}  at t={bt}   1/(m+1)={F(1,m+1)}")

print("\nmax_t min(observer-adjacent gaps), EXCLUDING point-collapse (all m+1 distinct):")
for m in [4,5,6,7]:
    speeds=list(range(1,m+1))
    b,bt = maximize(speeds, lambda s,t: min(observer_gaps(s,t)), require_distinct=True)
    print(f"  AP m={m}: {b} = {float(b):.5f}  at t={bt} (den {bt.denominator})   2/(m+1)={F(2,m+1)}={float(F(2,m+1)):.5f}")
