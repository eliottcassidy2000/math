#!/usr/bin/env python3
"""reader_ab1: intermediate claims of the written proof of THM-5141 (A), half split, (B1)
that the earlier verifiers did not test directly.
 T1  (A)1 two-line proof: for delta in {-2..2}^w, delta != 0, D(delta) = 0  ==>  the lowest nonzero digit is +-2
     and the partial reversed sums u_i = sum_{l<=i} delta_l 2^(i-l) satisfy |u_i| >= 2 for all i >= m
     (so D*(delta) = u_(w-1) != 0).  Exhaustive w <= 8.
 T2  (D1) per-direction multiplicity: for fixed P_u and fixed direction (k,k') at most 2 points P_v;
     for fixed P_v and fixed direction at most 4 points P_u.  Exhaustive over ALL integer points
     |a|,|b| < 2^s1 (b odd), |a'|,|b'| < 2^s2 (b' odd), w = 2..10 (not only the support of r).
 T3  (D2) intermediate: Deg_K <= 3 (2^(s2+1) (1+ln K) + K), Deg'_K <= 3 (2^(s1+1) (1+ln K) + K) for K <= 2^s2
     (the unsimplified bounds), same ranges, all K = 1..2^s2.
"""
import sys, itertools, math
from math import gcd
from collections import defaultdict

def t1(wmax=8):
    ok = True
    for w in range(1, wmax + 1):
        cnt = 0
        for d in itertools.product(range(-2, 3), repeat=w):
            if not any(d): continue
            if sum(di << i for i, di in enumerate(d)) != 0: continue
            cnt += 1
            m = next(i for i, di in enumerate(d) if di)
            if abs(d[m]) != 2: ok = False
            u = 0
            for i in range(w):
                u = 2 * u + d[i]
                if i >= m and abs(u) < 2: ok = False
        print(f"T1 w={w}: nonzero delta with D=0: {cnt}; lowest digit +-2 and |u_i|>=2 for i>=m: {ok}")
    return ok

def t23(wmin=2, wmax=10):
    allok = True
    for w in range(wmin, wmax + 1):
        s1 = (w + 1) // 2; s2 = w // 2
        A1 = 1 << s1; A2 = 1 << s2
        perdir_u = defaultdict(int); perdir_v = defaultdict(int)
        # degree tables by height
        degu = defaultdict(lambda: defaultdict(int)); degv = defaultdict(lambda: defaultdict(int))
        for a in range(-A1 + 1, A1):
            for b in range(-A1 + 1, A1):
                if b % 2 == 0: continue
                for ap in range(-A2 + 1, A2):
                    X0 = ap * A1 + b
                    for bp in range(-A2 + 1, A2):
                        if bp % 2 == 0: continue
                        Xi = a * A2 + bp
                        g = gcd(X0, Xi)
                        if Xi < 0: g = -g
                        k = X0 // g; kp = Xi // g
                        M = max(abs(k), kp)
                        perdir_u[(a, b, k, kp)] += 1
                        perdir_v[(ap, bp, k, kp)] += 1
                        if M <= A2:
                            degu[(a, b)][M] += 1; degv[(ap, bp)][M] += 1
        mu = max(perdir_u.values()); mv = max(perdir_v.values())
        # (D2) unsimplified, all K <= 2^s2
        worst_u = worst_v = 0.0
        for tab, s_other, store in ((degu, s2, 'u'), (degv, s1, 'v')):
            for P, hist in tab.items():
                c = 0
                for K in range(1, A2 + 1):
                    c += hist.get(K, 0)
                    bound = 3 * ((1 << (s_other + 1)) * (1 + math.log(K)) + K)
                    r = c / bound
                    if store == 'u': worst_u = max(worst_u, r)
                    else: worst_v = max(worst_v, r)
        ok = (mu <= 2) and (mv <= 4) and worst_u <= 1 and worst_v <= 1
        allok &= ok
        print(f"T2/T3 w={w} (s1={s1}, s2={s2}): max #P_v per (P_u,direction) = {mu} (<=2), "
              f"max #P_u per (P_v,direction) = {mv} (<=4); (D2) unsimplified worst ratios Deg {worst_u:.4f}, Deg' {worst_v:.4f}; ok={ok}")
        sys.stdout.flush()
    return allok

if __name__ == "__main__":
    wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    a = t1(8)
    b = t23(2, wmax)
    print("ALL OK" if (a and b) else "FAILURE")
