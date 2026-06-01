from fractions import Fraction as F

def observer_gaps(speeds, t):
    pts = {F(0)}
    for v in speeds:
        pts.add((v * t) % 1)
    p = sorted(pts)
    if len(p) == 1:
        return F(0), F(0)
    right = p[1] - p[0]
    left = (p[0] + 1) - p[-1]
    return left, right

def walls(speeds, T=1):
    """Collect all wall times in (0,T) where combinatorial structure changes.
       Walls: k/v_i, k/(2 v_i)  (covered by k/v_i and crossings), and k/|v_i-v_j|.
       We use: t where v_i t in Z (t=k/v_i), and v_i t ≡ v_j t mod 1 (t=k/|v_i-v_j|),
       and v_i t ≡ -v_j t i.e. t=k/(v_i+v_j) (point coincides with reflection through 0 ordering)."""
    W = set()
    vs = list(speeds)
    for v in vs:
        for k in range(0, v*T+1):
            W.add(F(k, v))
    for i in range(len(vs)):
        for j in range(len(vs)):
            if i==j: continue
            d = abs(vs[i]-vs[j])
            if d>0:
                for k in range(0, d*T+1):
                    W.add(F(k,d))
            s = vs[i]+vs[j]
            for k in range(0, s*T+1):
                W.add(F(k,s))
    W = sorted(w for w in W if 0 <= w <= T)
    return W

def max_min_observer(speeds):
    """Exact max over t in (0,1) of min(g_left, g_right) via cell midpoints.
       min(g_left,g_right) is piecewise of form (linear) but the points are
       linear in t within a cell; the two neighbors of 0 are fixed within a cell,
       so each observer gap is linear in t. min of two linear functions -> max
       on a cell is at an endpoint OR where they cross. We sample dense candidate
       t: all walls, midpoints, and pairwise-equality crossing points.
       To be safe & exact we evaluate at midpoints of cells AND at wall points,
       AND at internal crossings where g_left=g_right within each cell."""
    W = walls(speeds)
    cands = set()
    for i in range(len(W)-1):
        a, b = W[i], W[i+1]
        if b<=a: continue
        mid = (a+b)/2
        cands.add(mid)
        cands.add(a); cands.add(b)
        # within this cell neighbors fixed; find crossing of the two linear gaps.
        # evaluate gaps at two interior points and solve linear g_left=g_right.
        t1 = a + (b-a)/3
        t2 = a + 2*(b-a)/3
        l1,r1 = observer_gaps(speeds, t1)
        l2,r2 = observer_gaps(speeds, t2)
        # f(t)=l-r linear: f(t1),f(t2). root:
        f1 = l1-r1; f2 = l2-r2
        if f1 != f2:
            # t* where f=0
            tstar = t1 + (t2-t1) * (F(0)-f1)/(f2-f1)
            if a < tstar < b:
                cands.add(tstar)
    best = F(-1)
    bestt = None
    for t in cands:
        if not (0 < t < 1):
            continue
        l,r = observer_gaps(speeds, t)
        val = min(l,r)
        if val > best:
            best = val; bestt = t
    return best, bestt

print("="*70)
print("PART 3: max_t min(g_left, g_right)  [observer-adjacent gaps]")
print("="*70)
print("\nAP speeds v=1..m  (points k t, k=0..m):")
for m in [4,5,6,7]:
    speeds = list(range(1,m+1))
    best, bt = max_min_observer(speeds)
    print(f"  m={m}: max min(obs gaps) = {best} = {float(best):.6f}   at t={bt} (den {bt.denominator}),  1/(m+1)={F(1,m+1)}={float(F(1,m+1)):.6f}")

print("\nRandom speed sets:")
import random
random.seed(7)
for _ in range(6):
    m = random.choice([4,5])
    speeds = sorted(random.sample(range(1,15), m))
    best, bt = max_min_observer(speeds)
    print(f"  speeds={speeds}: max min = {best} = {float(best):.6f}  at t={bt} (den {bt.denominator})")
