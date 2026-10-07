# Stronger heuristic (tabu swap search with restarts, product seeds) for N = 13..20 sections of Z[w].
import random, math, sys, time
from itertools import product
def nrm(a, b): return a*a + a*b + b*b
def shell(D):
    R = int(math.isqrt(4*D)) + 2
    return [(a, b) for a in range(-R, R+1) for b in range(-R, R+1) if nrm(a, b) == D]
uN = {13:30,14:33,15:37,16:41,17:43,18:46,19:50,20:54}
def edges_of(S, sh):
    Ss = set(S); return sum(1 for (a, b) in S for (s, t) in sh if (a+s, b+t) in Ss)//2
def tabu(N, sh, seedset, iters, rnd):
    S = set(seedset)
    while len(S) > N:  # trim lowest degree
        v = min(S, key=lambda p: (sum((p[0]+s, p[1]+t) in S for s, t in sh), rnd.random())); S.remove(v)
    while len(S) < N:
        front = {}
        for (a, b) in S:
            for (s, t) in sh:
                p = (a+s, b+t)
                if p not in S: front[p] = front.get(p, 0) + 1
        mx = max(front.values()); S.add(rnd.choice([p for p, d in front.items() if d == mx]))
    E = edges_of(S, sh); best = E; tabu_list = {}
    for it in range(iters):
        degs = {p: sum((p[0]+s, p[1]+t) in S for s, t in sh) for p in S}
        front = {}
        for (a, b) in S:
            for (s, t) in sh:
                p = (a+s, b+t)
                if p not in S: front[p] = front.get(p, 0) + 1
        bestmove = None; bestdelta = -10**9
        for v, dv in degs.items():
            if tabu_list.get(v, -1) > it: continue
            for w, dw in front.items():
                if tabu_list.get(w, -1) > it: continue
                adj = 1 if (w[0]-v[0], w[1]-v[1]) in shset else 0
                delta = dw - adj - dv
                if delta > bestdelta or (delta == bestdelta and rnd.random() < 0.3):
                    bestdelta = delta; bestmove = (v, w)
        if bestmove is None: break
        v, w = bestmove; S.remove(v); S.add(w); E += bestdelta
        tabu_list[v] = it + rnd.randint(3, 8); tabu_list[w] = it + rnd.randint(1, 3)
        if E > best: best = E
    return best
t0 = time.time()
for N in range(13, 21):
    out = {}
    for D in (7, 49, 91, 133, 637):
        sh = shell(D); shset = set(sh); globals()["shset"] = shset
        rnd = random.Random(N*1000 + D)
        b = 0
        for r in range(6):
            # seed: random walk cluster
            seed = {(0, 0)}
            while len(seed) < N + 3:
                p = rnd.choice(list(seed)); s = rnd.choice(sh); seed.add((p[0]+s[0], p[1]+s[1]))
            b = max(b, tabu(N, sh, seed, 300, rnd))
        out[D] = b
    best = max(out.values())
    print(f"N={N}: tabu best {best} by D {out}  u(N)={uN[N]}  {'OPTIMAL' if best == uN[N] else 'short by ' + str(uN[N]-best)}  [{time.time()-t0:.0f}s]"); sys.stdout.flush()
