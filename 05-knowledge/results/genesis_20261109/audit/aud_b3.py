"""AUDIT G1: Lemma B3 structural test.  T on 0..n-1 whose back arcs (later -> earlier) have span <= s.
For every Hamiltonian path P and every v: the vertices <= v-1 form at most s+1 maximal runs of P (heads of runs not
containing P's first vertex lie in [v-s, v-1]); the encoding (pred choice, succ choice) per v, with runs indexed
canonically, is injective on paths; and H <= (s+3)^(2n)."""
import random, itertools
random.seed(11)
def paths(n, out):
    res = []
    def dfs(p, used):
        if len(p) == n:
            res.append(tuple(p)); return
        for w in range(n):
            if not used >> w & 1 and out[p[-1]] >> w & 1:
                p.append(w); dfs(p, used | 1 << w); p.pop()
    for v in range(n):
        dfs([v], 1 << v)
    return res
bad = 0; tests = 0; maxruns_seen = {}
for t in range(300):
    n = random.randint(4, 9); s = random.randint(1, 3)
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if j - i <= s and random.random() < 0.5: out[j] |= 1 << i
            else: out[i] |= 1 << j
    P = paths(n, out)
    tests += 1
    if len(P) > (s + 3) ** (2 * n): bad += 1
    codes = set()
    for p in P:
        pos = {v: k for k, v in enumerate(p)}
        enc = []
        for v in range(n):
            # runs of P among vertices < v
            runs = []
            k = 0
            while k < n:
                if p[k] < v:
                    st = k
                    while k < n and p[k] < v: k += 1
                    runs.append((p[st], p[k - 1], st))   # (head, tail, start index)
                else:
                    k += 1
            if len(runs) > s + 1: bad += 1
            maxruns_seen[s] = max(maxruns_seen.get(s, 0), len(runs))
            for (h, tl, st) in runs:
                if st != 0 and not (v - s <= h <= v - 1): bad += 1
            runs.sort()                       # canonical index: by head label
            k = pos[v]
            pred = 'none' if k == 0 else ('later' if p[k - 1] > v else [i for i, r in enumerate(runs) if r[1] == p[k - 1]][0])
            succ = 'none' if k == n - 1 else ('later' if p[k + 1] > v else [i for i, r in enumerate(runs) if r[0] == p[k + 1]][0])
            enc.append((pred, succ))
        codes.add(tuple(enc))
    if len(codes) != len(P): bad += 1
print('B3 structural tests:', tests, 'failures:', bad, ' max runs seen per s:', maxruns_seen)
