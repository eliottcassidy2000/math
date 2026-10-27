"""Analyses of the class move graphs dumped by dipole_exactness --dump (orders 4-8).
(1) P_7's transitive distance: the L-descending DAG below P_7, level by level, with failures and local minima.
(2) For every L-local minimum (target TT): is it an endpoint obstruction or a colour obstruction?
    For score potentials the change under a move depends only on the endpoints (a, b) of the path (dipole law
    d -> -d + e_a - e_b), so descent at T <=> some Hamiltonian path joins a descending start to a descending end.
(3) The endpoint relation Pi(T) = {(a, b): T has a Hamiltonian path from a to b} versus arc-strong connectivity.
Usage: python3 analyze_graph.py n
"""
import sys
from collections import deque, Counter
from itertools import permutations

n = int(sys.argv[1])
cls = []
for line in open(f"graph_n{n}.txt"):
    head, anyn, blk, arcs = line.rstrip("\n").split("|")
    f = head.split()
    c = dict(id=int(f[0]), scores=f[2], comps=f[3], aut=int(f[4]), sc=int(f[5]), lam=int(f[6]), d0=int(f[7]),
             L=int(f[8]) // 2, fas=int(f[9]), R2=int(f[10]), H=int(f[11]),
             any=[int(x) for x in anyn.split()], black=[int(x) for x in blk.split()], arcs=arcs.split()[1])
    cls.append(c)
C = len(cls)

def out_masks(arcs):
    out = [0] * n; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if arcs[k] == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    return out

def ham_endpoints(out):
    """set of (start, end) of Hamiltonian paths, by DP over subsets"""
    full = (1 << n) - 1
    res = set()
    for a in range(n):
        reach = {1 << a: {a}}
        # dp[S] = set of possible ends of paths from a covering S
        dp = [0] * (1 << n); dp[1 << a] = 1 << a
        for S in range(1 << n):
            e = dp[S]
            if not e or not (S >> a & 1): continue
            for v in range(n):
                if e >> v & 1:
                    nx = out[v] & ~S
                    while nx:
                        w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                        dp[S | 1 << w] |= 1 << w
        e = dp[full]
        for b in range(n):
            if e >> b & 1: res.add((a, b))
    return res

def bfs(targets, key):
    d = [-1] * C; q = deque()
    for c in range(C):
        if targets(cls[c]): d[c] = 0; q.append(c)
    while q:
        c = q.popleft()
        for t in cls[c][key]:
            if d[t] < 0: d[t] = d[c] + 1; q.append(t)
    return d

dTT = bfs(lambda c: c['L'] == 0, 'any')
fail = [c for c in range(C) if dTT[c] > cls[c]['L']]
lmin = [c for c in range(C) if cls[c]['L'] > 0 and not any(cls[t]['L'] == cls[c]['L'] - 1 for t in cls[c]['any'])]
print(f"n={n}: classes {C}; TT/L failures {len(fail)}; L-local minima {len(lmin)}")

# (1) P_7 (or any vertex-transitive regular class at odd n): descending DAG
if n % 2 == 1:
    regs = [c for c in range(C) if cls[c]['scores'] == str((n - 1) // 2) * n]
    for r in regs:
        print(f"regular class {r}: |Aut| {cls[r]['aut']} SC {cls[r]['sc']} H {cls[r]['H']} L {cls[r]['L']} d(TT) {dTT[r]} moves to {len(cls[r]['any'])} classes")
        level = {r}
        while level:
            Lv = cls[next(iter(level))]['L']
            nf = sum(1 for c in level if dTT[c] > cls[c]['L'])
            nl = sum(1 for c in level if c in set(lmin))
            print(f"   L={Lv}: {len(level)} classes reachable by L-descending words, {nf} of them fail (d > L), {nl} local minima;"
                  f" symmetric among failing: {sum(1 for c in level if dTT[c] > cls[c]['L'] and (cls[c]['aut'] > 1 or cls[c]['sc']))}")
            nxt = set()
            for c in level:
                for t in cls[c]['any']:
                    if cls[t]['L'] == cls[c]['L'] - 1: nxt.add(t)
            level = nxt
        # which local minima are hit
        seen = {r}; st = [r]; hit = Counter()
        while st:
            c = st.pop()
            for t in cls[c]['any']:
                if cls[t]['L'] == cls[c]['L'] - 1 and t not in seen:
                    seen.add(t); st.append(t)
        mins = [c for c in seen if c in set(lmin)]
        print(f"   local minima below it: {len(mins)}; by components: {Counter(cls[c]['comps'] for c in mins).most_common(8)}")
        print(f"   successful descents: {'yes' if any(cls[c]['L'] == 0 for c in seen) else 'none'}")

# (2) classify L-local minima: endpoint obstruction (no Hamiltonian path joins a descending pair) vs colour
def descending_pairs(scores_vec):
    """pairs (a, b) such that the move with endpoints a -> b lowers L by one"""
    L0 = Lof(scores_vec); res = set()
    for a in range(n):
        for b in range(n):
            if a == b: continue
            s2 = [n - 1 - x for x in scores_vec]; s2[a] += 1; s2[b] -= 1
            if Lof(s2) == L0 - 1: res.add((a, b))
    return res

def Lof(s):
    t = sorted(s); return sum(abs(t[i] - i) for i in range(n)) // 2

kinds = Counter(); lam_hist = Counter()
for c in lmin:
    out = out_masks(cls[c]['arcs'])
    sv = [bin(out[v]).count("1") for v in range(n)]
    D = descending_pairs(sv); P = ham_endpoints(out)
    if not D: kinds['no descending pair at all'] += 1
    elif not (D & P): kinds['endpoint obstruction (Pi misses D)'] += 1
    else: kinds['descending path exists but no descending neighbour?!'] += 1
    lam_hist[(cls[c]['lam'], 'strong' if cls[c]['lam'] > 0 else 'nonstrong')] += 1
print("L-local minima kinds:", dict(kinds))
print("L-local minima by lambda:", dict(lam_hist))

# (3) endpoint relation versus arc-strong connectivity: for which classes is Pi(T) 'semi-connected'
#     (every pair joined in at least one direction) / 'connected' (both directions)?
semi = Counter(); full = Counter()
for c in range(C):
    if cls[c]['lam'] == 0: continue
    out = out_masks(cls[c]['arcs']); P = ham_endpoints(out)
    allsemi = all(((a, b) in P) or ((b, a) in P) for a in range(n) for b in range(a + 1, n))
    allfull = all(((a, b) in P) and ((b, a) in P) for a in range(n) for b in range(a + 1, n))
    semi[(cls[c]['lam'], allsemi)] += 1; full[(cls[c]['lam'], allfull)] += 1
print("strong classes by (lambda, every pair joined in some direction):", dict(semi))
print("strong classes by (lambda, Hamiltonian-connected):", dict(full))
