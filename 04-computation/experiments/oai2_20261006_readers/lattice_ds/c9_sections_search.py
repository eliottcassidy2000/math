# Search (simulated annealing) for dense N-point sections of the triangular lattice Z[w], unit = sqrt(D).
# Lower bounds only (a found count is FINITE-EXACT, recounted exactly); compare with u(N) (AMP24 table).
import random, math, sys
def nrm(a, b): return a*a + a*b + b*b
def shell(D):
    R = int(math.isqrt(4*D)) + 2
    return [(a, b) for a in range(-R, R+1) for b in range(-R, R+1) if nrm(a, b) == D]
uN = {1:0,2:1,3:3,4:5,5:7,6:9,7:12,8:14,9:18,10:20,11:23,12:27,13:30,14:33,15:37,16:41,17:43,18:46,19:50,20:54,21:57}
ubounds = {22:(60,61),23:(64,66),24:(68,72),25:(72,78),26:(76,84),27:(81,90),28:(85,96),29:(89,103),30:(93,110)}
def count_edges(S, sh):
    Sset = set(S); return sum(1 for (a, b) in S for (s, t) in sh if (a+s, b+t) in Sset)//2
def anneal(N, sh, iters=20000, seed=0):
    rnd = random.Random(seed)
    S = {(0, 0)}
    while len(S) < N:  # greedy growth
        front = {}
        for (a, b) in S:
            for (s, t) in sh:
                p = (a+s, b+t)
                if p not in S: front[p] = front.get(p, 0) + 1
        mx = max(front.values()); cand = [p for p, d in front.items() if d == mx]
        S.add(rnd.choice(cand))
    def deg(p): return sum(1 for (s, t) in sh if (p[0]+s, p[1]+t) in S)
    E = count_edges(list(S), sh); best = E; bestS = set(S); T = 1.0
    for it in range(iters):
        T = max(0.05, 1.5*(1 - it/iters))
        v = rnd.choice(list(S)); dv = deg(v)
        S.remove(v)
        front = {}
        for (a, b) in S:
            for (s, t) in sh:
                p = (a+s, b+t)
                if p not in S and p != v: front[p] = front.get(p, 0) + 1
        if not front: S.add(v); continue
        items = list(front.items())
        mx = max(d for _, d in items)
        cand = [p for p, d in items if d >= mx - (1 if rnd.random() < 0.3 else 0)]
        w = rnd.choice(cand); dw = front[w]
        delta = dw - dv
        if delta >= 0 or rnd.random() < math.exp(delta/T):
            S.add(w); E += delta
            if E > best: best = E; bestS = set(S)
        else:
            S.add(v)
    assert count_edges(list(bestS), sh) == best
    return best, bestS
Ds = [1, 3, 7, 13, 49, 91]
shells = {D: shell(D) for D in Ds}
print("shell sizes r(D):", {D: len(shells[D]) for D in Ds})
res = {}
for N in range(2, 31):
    row = {}
    for D in Ds:
        b = 0
        for seed in range(3):
            e, _ = anneal(N, shells[D], iters=2500 if N <= 21 else 3500, seed=seed)
            b = max(b, e)
        row[D] = b
    best = max(row.values())
    ref = uN.get(N, ubounds.get(N))
    tag = ""
    if N in uN: tag = "= u(N) OPTIMAL" if best == uN[N] else f"< u(N)={uN[N]}"
    else: tag = f"AMP bounds {ref}" + (" (meets best-known lower bound)" if best >= ref[0] else "")
    print(f"N={N:2d}: best section edges {best:3d} by D: {row}  {tag}"); sys.stdout.flush()
