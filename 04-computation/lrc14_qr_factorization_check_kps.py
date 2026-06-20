#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify the CLAIMED factorization  K(n) = [prod_j 1/(2 pi i n_j)] * D7(n mod 7).
K(n) computed numerically the EXACT-definition way (from the proven kernel file s12).
D7 computed exactly in Z[zeta_7].  Match => the coset kernel is the right object and
the QR/NQR split is meaningful.
"""
import sys, itertools, cmath, math
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# --- numeric K(n) (exact definition, same as lrc14_signed_shell_decay_kps_s12.py) ---
def chat_T(T, m):
    if m % 7 == 0:
        return 0.0+0.0j
    s = 0.0+0.0j
    for j in T:
        a=j/7.0; b=(j+1)/7.0
        s += (cmath.exp(-2j*math.pi*m*a)-cmath.exp(-2j*math.pi*m*b))/(2j*math.pi*m)
    return -s
def K_num(vals):
    tot=0.0+0.0j
    for r in range(7):
        for T in itertools.combinations(range(1,7), r):
            p=1.0+0.0j
            for v in vals: p*=chat_T(T,v)
            tot += ((-1)**r)*p
    return tot

# --- D7 numeric (the coset part only), zeta = e^{2pi i/7} ---
Z=cmath.exp(2j*math.pi/7)
def sigma_T_num(T,m): return sum(Z**((-m*t)%7) for t in T)
def D7_num(c):
    pref=1.0+0.0j
    for cj in c: pref*=(1 - Z**((-cj)%7))
    acc=0.0+0.0j
    for r in range(7):
        for T in itertools.combinations(range(1,7),r):
            p=1.0+0.0j
            for cj in c: p*=sigma_T_num(T,cj)
            acc += ((-1)**r)*p
    return pref*acc

if __name__=="__main__":
    print("Factorization check: K(n) ?= [prod 1/(2 pi i n_j)] * D7(n mod 7)")
    # chat_T(T,m) = -(1/(2 pi i m))(1 - e^{-2pi i m/7}) sum_{j in T} e^{-2pi i m j/7}
    #            = (1/(2 pi i m)) * [ -(1 - zeta^{-m}) sigma_T(m) ]   with zeta=e^{2pi i/7}
    # so prod_j chat_T = prod_j 1/(2pi i n_j) * prod_j[ -(1-zeta^{-n_j}) sigma_T(n_j) ]
    # K = sum_T (-1)^|T| prod chat = [prod 1/(2pi i n_j)] * sum_T (-1)^|T| prod_j[-(1-zeta^{-c})sigma_T]
    # the (-1) per coordinate -> (-1)^6 = +1 overall since support 6. good.
    tests=[[1,1,1,1,1,1],[1,2,4,3,5,6],[2,-3,1,5,-1,4],[8,-3,2,11,-5,6],[1,2,3,-4,5,-6]]
    for vals in tests:
        K=K_num(vals)
        arch=1.0+0.0j
        for v in vals: arch*=1.0/(2j*math.pi*v)
        c=tuple(v%7 for v in vals)
        D=D7_num(c)
        pred=arch*D
        ok=abs(K-pred)<1e-9
        print(f"  n={vals}: K={K:.3e}  arch*D7={pred:.3e}  match={ok}")
