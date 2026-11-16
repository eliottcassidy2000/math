#!/usr/bin/env python3
"""converse_census.py: (A) near-join values h + 2^k B (k >= 1) from bases of order <= 8, up to 10000; the strong digit-0
pairs (S, z) of orders 4-8 whose base S - z has no sink (terminal component D strong): is H(S) a near-join value, and from
which least base order?  (B) the target chooses the core: representations v = h(A)h(D) + B_t(A) P_D(2) of odd v <= 125 with
cores A of order <= 6 (every root t) and right sources D in {M_1..M_8} or strong of order <= 5.  (C) one longest dead-end
descent chain at order 8 with the escapes H(S_i) - h(S_{i+1})."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
import functools
from collections import Counter, defaultdict
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def Bvec(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T))
    return H[N-1], [sum(dpr[m][t] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1) for t in range(n)]
def classes1(n): return [[[0]]] if n == 1 else classes(n)
def P2(D):
    n = len(D); N = 1 << n; dp, H = ham_all(D); return sum(H[S] * H[(N-1) & ~S] for S in range(N)), H[N-1]
def term_comp(B):
    n = len(B)
    def reach(s):
        r = {s}; f = [s]
        while f:
            a = f.pop()
            for b in range(n):
                if B[a][b] and b not in r: r.add(b); f.append(b)
        return r
    for v in range(n):
        c = scc_of(B, v)
        if reach(v) == c: return sorted(c)
def digit0(S):
    H, B = Bvec(S); return H, [z for z in range(len(S)) if B[z] == H]
LIMIT = 10000
# (A)
NJ = {}; maxH = Counter()
for n in range(1, 9):
    for T in classes1(n):
        h, B = Bvec(T); maxH[n] = max(maxH[n], h)
        for b in B:
            k = 1
            while h + (1 << k) * b <= LIMIT:
                v = h + (1 << k) * b
                if v not in NJ or NJ[v] > n: NJ[v] = n
                k += 1
print(f"(A) largest counts by order: {dict(sorted(maxH.items()))}; near-join values <= {LIMIT} from bases <= 8: {len(NJ)}; odd values <= 1000 missed: {[v for v in range(1, 1001, 2) if v not in NJ]}")
comp_cases = []; single = 0
for n in range(4, 9):
    for S in classes(n):
        if not strong(S): continue
        H, d0 = digit0(S)
        for z in d0:
            rest = [v for v in range(n) if v != z]; base = sub(S, rest); comp = term_comp(base)
            if len(comp) == 1: single += 1; continue
            t0 = rest.index([x for x in range(n) if S[z][x]][0])
            A_idx = [v for v in range(len(base)) if v not in comp]; A = sub(base, A_idx); D = sub(base, comp)
            hA, BA = Bvec(A); pD, hD = P2(D)
            comp_cases.append((H, n, hA, BA[A_idx.index(t0)], hD, pD, NJ.get(H)))
print(f"    strong digit-0 pairs of orders 4-8: {single} with a sink in the base (literal near-join, k >= 1), {len(comp_cases)} with a strong right source D")
by_order = Counter(c[6] for c in comp_cases)
print(f"    least near-join base order for the strong-D counts (None = not a near-join value from bases <= 8): {dict(sorted(by_order.items(), key=lambda x: (x[0] is None, x[0])))}")
print("    the strong-D cases (H, order, h(A), B_t(A), h(D), P_D(2), least near-join base order):")
for c in sorted(set(comp_cases)): print("     ", c)
# (B)
readers = [(1, 1 << k, f"M_{k}") for k in range(1, 9)]
for n in range(3, 6):
    for D in classes(n):
        if strong(D):
            pD, hD = P2(D); readers.append((hD, pD, f"strong{n}"))
print(f"(B) readers (h(D), P_D(2)): {[(a, b, c) for a, b, c in readers]}")
pairs = Counter()
for n in range(1, 7):
    for A in classes1(n):
        h, B = Bvec(A)
        for b in B: pairs[(h, b)] += 1
reps = defaultdict(lambda: Counter())
for (h, b), mult in pairs.items():
    for hD, pD, lab in readers:
        v = h * hD + b * pD
        if v <= 125: reps[v][lab.startswith('strong')] += mult
none = [v for v in range(1, 126, 2) if not reps[v]]
only_comp = [v for v in range(1, 126, 2) if reps[v] and reps[v][False] == 0]
only_nj = [v for v in range(1, 126, 2) if reps[v] and reps[v][True] == 0]
print(f"    odd v <= 125 with no representation: {none}")
print(f"    odd v <= 125 with only composite (strong D) representations: {only_comp}")
print(f"    odd v <= 125 with only near-join (D = M_k) representations: {only_nj}")
print("    v: (near-join reps, composite reps) for odd v <= 125:")
print("    " + "; ".join(f"{v}:({reps[v][False]},{reps[v][True]})" for v in range(1, 126, 2)))
# (C)
def key(adj): return tuple(tuple(r) for r in adj)
@functools.lru_cache(None)
def descent(k):
    adj = [list(r) for r in k]; n = len(adj); H, d0 = digit0(adj)
    if not d0 or n == 1: return 0
    return 1 + max(descent(key(sub(adj, [v for v in range(n) if v != z]))) for z in d0)
done = False
for S in classes(8):
    if not strong(S) or descent(key(S)) != 7: continue
    cur = S; print("(C) a longest dead-end descent at order 8 (order, H, deleted vertex, out-degree, escape H(S) - h(S - z)):")
    while True:
        n = len(cur); H, d0 = digit0(cur)
        if not d0 or n == 1: print(f"    sterile endpoint: order {n}, H = {H}"); break
        z = max(d0, key=lambda z: descent(key(sub(cur, [v for v in range(n) if v != z]))))
        nxt = sub(cur, [v for v in range(n) if v != z]); hn = ham_all(nxt)[1][-1] if n > 1 else 1
        print(f"    order {n}, H = {H}, z = {z}, out-degree {sum(cur[z])}, escape {H - hn}")
        cur = nxt
    break
