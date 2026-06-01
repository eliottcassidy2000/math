from math import gcd
from functools import reduce
from itertools import combinations
import sys

def norm_q(num, q):
    # ||num/q|| as a fraction value num'/q in [0, q/2], returns integer numerator over q
    r = num % q
    return min(r, q - r)

def M_fast(vs):
    """Exact max of D over t in [0,1). Returns (best_num, best_q, t_num, t_q).
    D max occurs at t = a/q with q in denominators from crossings/peaks.
    We collect candidate t as (num,q) reduced, evaluate D exactly via integer mod.
    Value at t=a/q for runner v: ||v*a/q|| = min(m, q-m)/q where m=(v*a)%q.
    D(t) = min_i that /q  => we compare fractions with possibly different q via cross-mult.
    """
    vs = list(vs)
    # collect candidate denominators q and the t-values a/q
    cand = set()  # store (num, q) in lowest terms
    def add(num, q):
        if q <= 0: return
        num %= q
        g = gcd(num, q)
        cand.add((num//g, q//g))
    # tent peaks t=(2k+1)/(2v)
    for v in vs:
        q = 2*v
        for k in range(q):
            num = 2*k+1
            if num < q:
                add(num, q)
    # crossings t=k/(vi+vj), k/(vi-vj)
    n = len(vs)
    for i in range(n):
        for j in range(i+1, n):
            for s in (vs[i]+vs[j], vs[i]-vs[j]):
                if s == 0: continue
                s = abs(s)
                for k in range(s):
                    add(k, s)
    # evaluate
    best_num, best_den = 0, 1
    bt = (0,1)
    for (a, q) in cand:
        # D = min over v of norm_q(v*a, q) ; this value = mnum/q
        mnum = q  # large
        for v in vs:
            r = (v*a) % q
            nv = r if r <= q-r else q-r
            if nv < mnum:
                mnum = nv
                if mnum == 0:
                    break
        # compare mnum/q vs best_num/best_den
        if mnum * best_den > best_num * q:
            best_num, best_den = mnum, q
            bt = (a, q)
    g = gcd(best_num, best_den)
    return (best_num//g, best_den//g, bt)

def search(m, cap):
    res = []
    for vs in combinations(range(1, cap+1), m):
        if reduce(gcd, vs) != 1: continue
        bn, bd, bt = M_fast(vs)
        # M*(m+1) = bn*(m+1)/bd
        res.append((bn*(m+1), bd, vs, (bn,bd), bt))
    # sort by value bn*(m+1)/bd ascending
    res.sort(key=lambda r: (r[0]/r[1], r[2]))
    return res

if __name__ == '__main__':
    m = int(sys.argv[1]); cap = int(sys.argv[2])
    res = search(m, cap)
    print(f'=== FAST m={m}, speeds<= {cap}, primitive sets={len(res)} ===')
    v0 = res[0][0]/res[0][1]
    print(f'min M*(m+1)={v0:.6f}')
    print(f'{"M*(m+1)":>10} {"speeds":<26} {"M":>10} {"t":>10} {"den(t)":>7} {"gap":>14}')
    from fractions import Fraction as Fr
    for num_m1, den, vs, (bn,bd), (ta,tq) in res[:12]:
        Mval = Fr(bn,bd)
        val = Fr(num_m1, den)
        gap = Mval - Fr(1, m+1)
        print(f'{str(val):>10} {str(list(vs)):<26} {str(Mval):>10} {str(Fr(ta,tq)):>10} {tq:>7} {str(gap):>14}')
