"""Wave 59.  Count spectrum of the class D(P): closure of K1 under ordinal sums and substitution Q[S_1..S_k] for prime quotients Q in P.
Profiles as in dgen.py; substitution chi(t) = sum_r N^Q_t(r) prod_i p_{S_i}(r_i) with N^Q_t(r) = # t-tuples of possibly empty
walks in Q with visit vector r (THM-3121).  Usage: pgen.py NMAX X LETTERS   (LETTERS: comma list from C3,Q9,Q13,Q15,ALL6)"""
import sys, itertools
from math import comb
from functools import lru_cache
NMAX = int(sys.argv[1]); X = int(sys.argv[2]); LET = sys.argv[3].split(',')
def modules(n, out):
    for m in range(1, (1 << n) - 1):
        if bin(m).count('1') < 2: continue
        if all(len({out[x] >> y & 1 for y in range(n) if m >> y & 1}) == 1 for x in range(n) if not m >> x & 1): return True
    return False
def Hq(n, out):
    from itertools import permutations
    return sum(1 for p in permutations(range(n)) if all(out[p[i]] >> p[i + 1] & 1 for i in range(n - 1)))
primes = {}
for line in open('w57_reps.txt'):
    t = list(map(int, line.split())); n = t[0]; out = tuple(t[1:])
    if n in (3, 5, 6) and not modules(n, out):
        h = Hq(n, out)
        if n == 3: primes['C3'] = out
        elif n == 5: primes['Q%d' % h] = out
        else: primes.setdefault('P6', []).append(out)
quots = []
for L in LET:
    if L == 'ALL6': quots += [(6, o) for o in primes['P6']]
    else: quots.append((len(primes[L]), primes[L]))
def chi_of(p, t): return sum(comb(t, d) * p[d - 1] for d in range(1, len(p) + 1))
def p_of_chi(ch, n): return tuple(sum((-1) ** (d - t) * comb(d, t) * ch[t] for t in range(d + 1)) for d in range(1, n + 1))
def join(pa, pb):
    n = len(pa) + len(pb); return p_of_chi([chi_of(pa, t) * chi_of(pb, t) for t in range(n + 1)], n)
@lru_cache(None)
def walkpoly(q, out, sizes):
    """W: dict visit-vector -> number of nonempty walks in Q (vertex sequences along arcs) with that visit vector, visits <= sizes"""
    W = {}
    frontier = {}
    for v in range(q):
        r = [0] * q; r[v] = 1
        if r[v] <= sizes[v]: frontier[(v, tuple(r))] = 1
    while frontier:
        nf = {}
        for (v, r), c in frontier.items():
            W[r] = W.get(r, 0) + c
            for w in range(q):
                if w != v and out[v] >> w & 1 and r[w] < sizes[w]:
                    r2 = list(r); r2[w] += 1; key = (w, tuple(r2)); nf[key] = nf.get(key, 0) + c
        frontier = nf
    return tuple(W.items())
@lru_cache(None)
def Npoly(q, out, sizes):
    n = sum(sizes); W = walkpoly(q, out, sizes)
    res = []; cur = {tuple([0] * q): 1}; res.append(tuple(cur.items()))
    for t in range(1, n + 1):
        nxt = dict(cur)
        for r, c in cur.items():
            for v, w in W:
                s = tuple(a + b for a, b in zip(r, v))
                if all(s[i] <= sizes[i] for i in range(q)): nxt[s] = nxt.get(s, 0) + c * w
        cur = nxt; res.append(tuple(cur.items()))
    return tuple(res)
def subst(q, out, ps):
    sizes = tuple(len(p) for p in ps); n = sum(sizes); N = Npoly(q, out, sizes)
    ch = []
    for t in range(n + 1):
        s = 0
        for r, c in N[t]:
            if all(r):
                pr = c
                for i in range(q): pr *= ps[i][r[i] - 1]
                s += pr
        ch.append(s)
    return p_of_chi(ch, n)
def compositions(n, k):
    if k == 1:
        yield (n,); return
    for a in range(1, n - k + 2):
        for rest in compositions(n - a, k - 1): yield (a,) + rest
prof = {1: {(1,)}}; strong = {1: {(1,)}}
for n in range(2, NMAX + 1):
    S = set()
    for (q, out) in quots:
        if q > n: continue
        for sz in compositions(n, q):
            lists = [prof[s] for s in sz]
            for ps in itertools.product(*lists):
                p = subst(q, out, ps)
                if p[0] <= X * 50: S.add(p)          # prune huge counts (cannot feed counts <= X)
    A = set(S)
    for a in range(1, n):
        for pa in strong[a]:
            for pb in prof[n - a]:
                if pa[0] * pb[0] <= X * 50: A.add(join(pa, pb))
    strong[n] = S; prof[n] = A
    print(f"order {n}: profiles {len(A)} (strong {len(S)})", flush=True)
spec = set(p[0] for n in prof for p in prof[n])
miss = [h for h in range(1, X + 1, 2) if h not in spec]
print(f"letters {LET}: odd h <= {X} missed at orders <= {NMAX}: {len(miss)}; first: {miss[:50]}")
def isprime(q): return q > 1 and all(q % d for d in range(2, int(q ** .5) + 1))
def isatom(h): return (isprime(h) and h != 7) or (h % 7 == 0 and isprime(h // 7) and h // 7 >= 5) or h in (63, 343)
print(f"atoms < 700 realized: {sorted(h for h in spec if h < 700 and isatom(h))}")
