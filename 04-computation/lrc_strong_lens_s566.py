#!/usr/bin/env python3
"""The strong lens on the LRC. opus-2026-06-02-S566.
THM-381: LRC <=> observer is a source at some time. S525/HYP-2000: at a lonely time
the runner sub-tournament has #SCC in {1, m} -- transitive (runners in a semicircle)
or a single STRONG block (runners wrap fully around the observer). We classify how
configs realize loneliness through strong-connectivity."""
from fractions import Fraction
from math import gcd
import random
n = 14

def safe(V, t):
    return all(min((v*t) % 1, 1-(v*t) % 1) >= Fraction(1, n) for v in V)

def lonely_time(V):
    eps = set()
    for v in V:
        for k in range(v+1):
            for s in (-1, 1):
                eps.add(Fraction(k*n+s, n*v) % 1)
    pts = sorted(eps); L = len(pts)
    for i in range(L):                          # prefer open/interior (generic)
        a, b = pts[i], pts[(i+1) % L]; ln = (b-a) if b > a else (b-a+1)
        mid = (a + ln/2) % 1
        if all(min((v*mid) % 1, 1-(v*mid) % 1) > Fraction(1, n) for v in V):
            return mid, True
    for t in pts:                               # boundary (tight)
        if safe(V, t):
            return t, False
    return None, None

def runner_scc(V, t):
    m = len(V); pos = [(v*t) % 1 for v in V]
    adj = [[(i != j and 0 < (pos[i]-pos[j]) % 1 < Fraction(1, 2))
            for j in range(m)] for i in range(m)]
    def reach(s):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for w in range(m):
                if adj[u][w] and w not in seen:
                    seen.add(w); st.append(w)
        return seen
    R = [reach(i) for i in range(m)]
    comp = [-1]*m; c = 0
    for i in range(m):
        if comp[i] == -1:
            for j in (j for j in range(m) if j in R[i] and i in R[j]):
                comp[j] = c
            c += 1
    return c

def in_semicircle(V, t):
    pos = sorted(float((v*t) % 1) for v in V); pos += [pos[0]+1]
    return max(pos[i+1]-pos[i] for i in range(len(V))) > 0.5

def prim(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v//g for v in V))

if __name__ == "__main__":
    rng = random.Random(2)
    samp = [prim(tuple(rng.sample(range(1, 60), 13))) for _ in range(150)]
    samp += [tuple(range(1, 14)), (1,2,3,4,5,6,7,8,9,10,11,13,24)]
    strong = transitive = other = semi = 0
    for V in samp:
        t, isopen = lonely_time(V)
        if t is None:
            continue
        sc = runner_scc(V, t)
        strong += sc == 1; transitive += sc == 13; other += 1 < sc < 13
        semi += in_semicircle(V, t)
    print(f"{len(samp)} configs, runner-block at a lonely time:")
    print(f"  STRONG (#SCC=1): {strong}   transitive (#SCC=13): {transitive}   "
          f"intermediate: {other}")
    print(f"  runners fit in a semicircle (transitive-realizable) at that time: {semi}")
    # AP and V* specifically
    for name, V in [("AP", tuple(range(1, 14))), ("V*", (1,2,3,4,5,6,7,8,9,10,11,13,24))]:
        t, op = lonely_time(V)
        print(f"  {name}: lonely t={t} ({'open' if op else 'boundary'}), "
              f"runner #SCC={runner_scc(V,t)}, semicircle={in_semicircle(V,t)}")
