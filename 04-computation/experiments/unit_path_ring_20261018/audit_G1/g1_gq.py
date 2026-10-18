#!/usr/bin/env python3
"""audit G1 (independent): finite unit-circle graphs G_q = Cay(F_(q^2), mu_(q+1)) for THM-4701.
Own implementation: F_(p^k) via a PRIMITIVE polynomial found by search, elements encoded as base-p integers,
mu_(q+1) = <x^(q-1)> (x primitive), so no norm formula is used. Petersen copies are counted on the Kneser model K(5,2)
with a different vertex order from the session script. Usage: python3 g1_gq.py [part]"""
import sys, itertools, time
from sympy import factorint, primerange

def prim_poly(p, k):
    """Return monic primitive polynomial coefficients c[0..k-1] (x^k = sum c_i x^i) and the power table of x."""
    Q = p ** k
    fac = list(factorint(Q - 1))
    for coeffs in itertools.product(range(p), repeat=k):
        if coeffs[0] == 0: continue
        # powers of x: vectors of length k
        cur = [0] * k; cur[0] = 1
        pw = []
        ok = True
        for i in range(Q - 1):
            pw.append(tuple(cur))
            # multiply by x
            top = cur[-1]
            cur = [0] + cur[:-1]
            if top:
                cur = [(cur[t] + top * coeffs[t]) % p for t in range(k)]
            if i < Q - 2 and cur == [1] + [0] * (k - 1):
                ok = False; break
        if ok and cur == [1] + [0] * (k - 1):
            return coeffs, pw
    raise RuntimeError('no primitive polynomial')

class GF:
    def __init__(self, p, k):
        self.p, self.k, self.Q = p, k, p ** k
        self.coeffs, pw = prim_poly(p, k)
        enc = lambda v: sum(c * p ** t for t, c in enumerate(v))
        self.exp = [enc(v) for v in pw]                 # exp[i] = code of x^i
        self.log = {c: i for i, c in enumerate(self.exp)}
        self.digits = [[(n // p ** t) % p for t in range(k)] for n in range(self.Q)]
    def add(self, a, b):
        p = self.p; da, db = self.digits[a], self.digits[b]
        return sum(((x + y) % p) * p ** t for t, (x, y) in enumerate(zip(da, db)))
    def neg(self, a):
        p = self.p
        return sum(((-x) % p) * p ** t for t, x in enumerate(self.digits[a]))

def gq(q):
    (p, e), = factorint(q).items()
    F = GF(p, 2 * e)
    Q = F.Q
    mu = [F.exp[(q - 1) * i % (Q - 1)] for i in range(q + 1)]
    assert len(set(mu)) == q + 1
    # translation table: add[a][m] for m in mu
    addm = [[F.add(a, m) for m in mu] for a in range(Q)]
    nb = [frozenset(r) for r in addm]
    return F, mu, nb

def local_data(q):
    F, mu, nb = gq(q)
    n0 = nb[0]
    mx = max(len(n0 & nb[d]) for d in range(1, F.Q))
    tri = sorted({len(n0 & nb[m]) for m in mu})
    k4 = any(all(b in nb[a] for a, b in itertools.combinations(c, 2)) for c in itertools.combinations(sorted(n0), 3))
    return mx, tri, k4, F, mu, nb

def srg_params(vs, nb):
    vs = list(vs); S = set(vs)
    degs = {len(nb[v] & S) for v in vs}
    lam = {len(nb[a] & nb[b] & S) for a in vs for b in vs if a != b and b in nb[a]}
    mu_ = {len(nb[a] & nb[b] & S) for a in vs for b in vs if a != b and b not in nb[a]}
    return len(vs), degs, lam, mu_

# Petersen graph as Kneser K(5,2)
PV = list(itertools.combinations(range(5), 2))
PE = {(a, b) for a in PV for b in PV if not set(a) & set(b)}
padj = {v: [w for w in PV if (v, w) in PE] for v in PV}
# order: fixed edge (0,1)-(2,3), then 5-cycle (0,1)-(2,3)-(0,4)-(1,2)-(3,4), then the rest by max placed neighbours
order = [(0, 1), (2, 3), (0, 4), (1, 2), (3, 4)]
rest = [v for v in PV if v not in order]
while rest:
    rest.sort(key=lambda v: -sum(w in order for w in padj[v]))
    order.append(rest.pop(0))
assert ((0, 4), (1, 2)) in PE and ((1, 2), (3, 4)) in PE and ((3, 4), (0, 1)) in PE and ((2, 3), (0, 4)) in PE
back = [[order.index(w) for w in padj[v] if order.index(w) < i] for i, v in enumerate(order)]

def petersen_count(q, exist_only=False, F=None, mu=None, nb=None):
    if F is None: F, mu, nb = gq(q)
    one = F.exp[0]
    img = [None] * 10; img[0] = 0; img[1] = one
    used = {0, one}
    cnt = 0
    def rec(i):
        nonlocal cnt
        if i == 10:
            cnt += 1
            return exist_only
        cands = None
        for j in back[i]:
            s = nb[img[j]]
            cands = s if cands is None else cands & s
        for c in cands:
            if c in used: continue
            img[i] = c; used.add(c)
            stop = rec(i + 1)
            used.discard(c)
            if stop: return True
        return False
    rec(2)
    return cnt

if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if part in ('all', 'local'):
        print('(1) local data: q, max common neighbours of two distinct vertices, triangles per edge (set over edges at 0), K4 present, predicted')
        for q in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 64]:
            mx, tri, k4, F, mu, nb = local_data(q)
            p = list(factorint(q))[0]; e = factorint(q)[p]
            pred_tri = 2 if (q % 6 == 5 or (p == 2 and e % 2 == 1)) else (1 if p == 3 else 0)
            pred_k4 = (p == 2 and e % 2 == 1)
            print('  q = %2d: max CN %d ; triangles per edge %s (pred %d) ; K4 %s (pred %s) %s' % (q, mx, tri, pred_tri, k4, pred_k4,
                  'OK' if (mx <= 2 and tri == [pred_tri] and k4 == pred_k4) else 'MISMATCH'), flush=True)
            if q == 4:
                print('     G_4 srg parameters (v, k, lambda, mu):', srg_params(range(F.Q), nb), '(Clebsch is the unique srg(16,5,0,2))')
                non = [v for v in range(1, F.Q) if v not in nb[0]]
                sums = {F.add(mu[a], mu[b]) for a in range(5) for b in range(a + 1, 5)}
                print('     non-neighbours of 0: %d ; equal to the ten two-sums of mu_5: %s ; induced srg parameters %s (Petersen is the unique srg(10,3,0,1))'
                      % (len(non), set(non) == sums, srg_params(non, nb)))
                lab = {F.add(mu[a], mu[b]): (a, b) for a in range(5) for b in range(a + 1, 5)}
                kn = all((b in nb[a]) == (not set(lab[a]) & set(lab[b])) for a in non for b in non if a != b)
                print('     adjacency of two-sums = disjointness of index pairs (Kneser K(5,2)):', kn)
                print('     sum of the two-sums of disjoint pairs is the fifth root:', all(F.add(F.add(mu[a], mu[b]), F.add(mu[c], mu[d])) == mu[({0,1,2,3,4} - {a,b,c,d}).pop()]
                      for (a, b) in itertools.combinations(range(5), 2) for (c, d) in itertools.combinations(range(5), 2) if not {a, b} & {c, d}))
            if q == 3:
                print('     G_3 srg parameters:', srg_params(range(F.Q), nb), '(K3 x K3 = srg(9,4,1,2))')
    if part in ('all', 'census'):
        print('(2) Petersen census: injective homomorphisms with a fixed directed edge -> (0, 1)')
        for q in list(primerange(3, 48)) + [9, 25, 27]:
            t0 = time.time(); c = petersen_count(q)
            print('  q = %2d: %d  (%.1f s)' % (q, c, time.time() - t0), flush=True)
    if part in ('all', 'exist'):
        print('(3) existence spot-checks and exhaustive q = 73')
        for q in [53, 59, 61, 67, 71, 79, 83, 89, 97, 101, 103, 107, 109, 113]:
            t0 = time.time(); c = petersen_count(q, exist_only=True)
            print('  q = %3d: exists %s  (%.1f s)' % (q, c > 0, time.time() - t0), flush=True)
    if part in ('all', 'q73'):
        t0 = time.time(); c = petersen_count(73)
        print('  q = 73 exhaustive count: %d  (%.1f s)' % (c, time.time() - t0), flush=True)
    if part == 'count113':
        print('(4) exhaustive counts for all primes 53 <= q <= 113 (and prime powers 49, 121 is skipped)')
        for q in list(primerange(53, 114)) + [49]:
            t0 = time.time(); c = petersen_count(q)
            print('  q = %3d (q mod 20 = %2d): %d  (%.1f s)' % (q, q % 20, c, time.time() - t0), flush=True)
