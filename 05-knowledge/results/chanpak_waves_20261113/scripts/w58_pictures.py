"""Matching pictures X -> Q[..., X at slot s, ...] for prime quotients Q (others singletons).
Count = sum_v W_v p_X(v), p_X(v) = ordered v-path covers = Delta^v chi_X(0)  (THM-3121, THM-3166).
Express as sum_j lam_j chi_X(j).  Check against brute force on random X.  Primes of order 3 and 5 from reps.txt."""
import itertools, random, sys
sys.path.insert(0, '.')
from w45_hk import H_of
from math import comb
def modules(n, out):
    for m in range(1, (1 << n) - 1):
        if bin(m).count('1') < 2: continue
        if all(len({out[x] >> y & 1 for y in range(n) if m >> y & 1}) == 1 for x in range(n) if not m >> x & 1):
            return True
    return False
def walkW(q, out, s, vmax=8):
    """W_v = # walks in Q visiting s exactly v times and every other vertex exactly once (consecutive distinct, arcs)."""
    from functools import lru_cache
    full = ((1 << q) - 1) & ~(1 << s)
    W = [0] * (vmax + 1)
    @lru_cache(None)
    def go(cur, mask, v):
        res = {}
        if mask == full and v >= 1: res[v] = res.get(v, 0) + 1
        for nx in range(q):
            if nx == cur or not out[cur] >> nx & 1: continue
            if nx == s:
                if v < vmax:
                    for k, c in go(nx, mask, v + 1).items(): res[k] = res.get(k, 0) + c
            elif not mask >> nx & 1:
                for k, c in go(nx, mask | 1 << nx, v).items(): res[k] = res.get(k, 0) + c
        return res
    tot = {}
    for st in range(q):
        r = go(st, 0, 1) if st == s else go(st, 1 << st, 0)
        for k, c in r.items(): tot[k] = tot.get(k, 0) + c
    for k, c in tot.items(): W[k] = c
    return W
def substitute(q, out, s, xn, xout):
    """blow up slot s of Q into tournament X"""
    others = [v for v in range(q) if v != s]; N = q - 1 + xn
    pos = {}; idx = 0
    for v in range(q):
        if v == s: pos[v] = list(range(idx, idx + xn)); idx += xn
        else: pos[v] = [idx]; idx += 1
    o = [0] * N
    for v in range(q):
        for w in range(q):
            if v != w and out[v] >> w & 1:
                for a in pos[v]:
                    for b in pos[w]: o[a] |= 1 << b
    for i in range(xn):
        for j in range(xn):
            if xout[i] >> j & 1: o[pos[s][i]] |= 1 << pos[s][j]
    return N, o
def chi(xn, xout, t):
    """chi_X(t) = sum over maps V -> [t] of prod H(fibres) = ordered t-tuples of possibly empty paths"""
    from functools import lru_cache
    Hs = {}
    def Hm(mask):
        if mask not in Hs:
            vs = [v for v in range(xn) if mask >> v & 1]
            if not vs: Hs[mask] = 1
            else:
                idx = {v: i for i, v in enumerate(vs)}
                Hs[mask] = H_of(len(vs), [sum(1 << idx[w] for w in vs if xout[v] >> w & 1) for v in vs])
        return Hs[mask]
    @lru_cache(None)
    def f(mask, t):
        if t == 1: return Hm(mask)
        tot = 0; sub = mask
        while True:
            tot += Hm(sub) * f(mask & ~sub, t - 1)
            if sub == 0: break
            sub = (sub - 1) & mask
        return tot
    return f((1 << xn) - 1, t)
primes = []
for line in open('w57_reps.txt'):
    t = list(map(int, line.split())); n = t[0]; out = t[1:]
    if n in (3, 5) and not modules(n, out): primes.append((n, out))
print("prime quotients of order 3 and 5:", len(primes))
random.seed(1)
for (q, out) in primes:
    seenW = set()
    for s in range(q):
        W = walkW(q, out, s)
        if tuple(W) in seenW: continue
        seenW.add(tuple(W)); vm = max(k for k in range(len(W)) if W[k])
        # lam_j = sum_v W_v (-1)^(v-j) C(v,j)
        lam = [sum(W[v] * (-1) ** (v - j) * comb(v, j) for v in range(j, vm + 1)) for j in range(vm + 1)]
        HQ = H_of(q, out)
        # check on random X
        ok = True
        for _ in range(12):
            xn = random.randint(1, 5); xo = [0] * xn
            for i in range(xn):
                for j in range(i + 1, xn):
                    if random.random() < .5: xo[i] |= 1 << j
                    else: xo[j] |= 1 << i
            N, o = substitute(q, out, s, xn, xo)
            pred = sum(lam[j] * chi(xn, tuple(xo), j) for j in range(1, vm + 1)) if vm >= 1 else 0
            if H_of(N, o) != pred: ok = False
        terms = " + ".join(f"{lam[j]}*chi({j})" for j in range(1, vm + 1) if lam[j])
        print(f"Q order {q}, H(Q)={HQ}, slot {s}: walk counts W_v={W[1:vm+1]}  ->  H(Q[X]) = {terms}   (checked: {ok})")
