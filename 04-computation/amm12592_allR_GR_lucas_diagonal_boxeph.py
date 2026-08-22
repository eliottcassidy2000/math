#!/usr/bin/env python3
"""ANGLE B3 step 2 -- explicit algebraic expansions of G_R (routes (b),(c)).

G_R := (q^{R-1} - E_{R-1})/2 at dyadic R (integer coefficients exactly there).
Derived exactly here, then machine-verified:

(1) STEP LAW      2(G_{R+1} - G_R) = -p L_{R-1},  L_m := p^m + q^m  (Lucas in e = pq).
(2) DOUBLING LAW  G_{2R} = 2 q G_R^2 + (p-q)(2 G_R - 1) - p^R q^{R-1}.
    [derivation: 2 q G_R = L_R - (p-q); L_{2R} = L_R^2 - 2 e^R; (p-q)^2 = 1-4e.]
(3) LUCAS DIAGONAL EXPANSION (the staircase basis {p^k q^{k-1}} is triangular,
    so these coefficients are UNIQUE):
        2 G_R = 2 + sum_{k=1}^{floor(R/2)} a_k(R) p^k q^{k-1},
        a_k(R) = (-1)^k (R/(R-k)) C(R-k,k) = (-1)^k [C(R-k,k) + C(R-k-1,k-1)],
    i.e. L_R = 1 + sum_k a_k e^k and division of L_R - (p-q) by 2q is EXPLICIT:
    (1 - (p-q))/q = 2 and e^k/q = p^k q^{k-1}.
(4) DYADIC EVENNESS: all a_k(R), k >= 1, are even  <=>  R = 2^t
    (L_R == 1 mod 2 iff R dyadic; L_{2m} == L_m^2 mod 2, L_1 = 1, L_3 = 1+e,...).
    So m_k := a_k/2 in Z exactly at dyadic R: G_R = 1 + sum_k m_k p^k q^{k-1}.
(5) closed forms: G_2 = q, G_4 = q^2 - p^3 = 1 - 2p + p^2 q;  the R = 4 epoch
    has the closed-form witness gamma_0 = q^2, gamma_1 = -p^2, gamma_2 = 0.
(6) magnitude profile: max_k |m_k| ~ phi^R-scale vs row capacity 2^{d_i};
    exact peak location and bit sizes tabulated for R = 8..1024.
(7) x=1 ledger: G_R(1) = -(R-2)/2 carried ENTIRELY by the head 1 + m_1 p
    (all k >= 2 atoms vanish at x = 1); m_1 = -R/2.

Exact int arithmetic only.
"""
import sys
from math import comb

WT = "/tmp/math-wt-boxeph-multifront"
OUT = WT + "/05-knowledge/results/amm12592_allR_GR_lucas_diagonal_boxeph.out"

log_lines = []
def log(s=""):
    print(s); log_lines.append(s)
def flush():
    with open(OUT, "w") as f: f.write("\n".join(log_lines) + "\n")
FAIL = 0
def check(name, ok, detail=""):
    global FAIL
    if not ok: FAIL += 1
    log("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- "+detail) if detail else ""))
    return ok

def ptrim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a
def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i,v in enumerate(a): r[i] += v
    for i,v in enumerate(b): r[i] += v
    return ptrim(r)
def psub(a, b): return padd(a, [-v for v in b])
def pmul(a, b):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        if u:
            for j,v in enumerate(b): r[i+j] += u*v
    return ptrim(r)
def pshift(a, s): return ([0]*s + list(a)) if a else []
def pscale(a, c): return [c*v for v in a] if c else []
def qpow(m): return [((-1)**k)*comb(m,k) for k in range(m+1)]
def Em(m): return ptrim([-1] + [1]*m)
def twoG(R): return psub(qpow(R-1), Em(R-1))
def G(R):
    t = twoG(R)
    assert all(v % 2 == 0 for v in t), "G_%d not integer (R not dyadic)" % R
    return [v//2 for v in t]
def Lpoly(m):
    """L_m = p^m + q^m as poly in x."""
    return padd(pshift([1], m), qpow(m))
def a_coef(R, k):
    return (-1)**k * (comb(R-k, k) + (comb(R-k-1, k-1) if k >= 1 else 0)) if k >= 1 else 1

log("="*78); log("(1) step law: 2(G_{R+1} - G_R) = -p L_{R-1}"); log("="*78)
ok = True
for R in range(2, 400):
    lhs = psub(twoG(R+1), twoG(R))
    rhs = pscale(pshift(Lpoly(R-1), 1), -1)
    if lhs != rhs: ok = False; break
check("2(G_{R+1}-G_R) = -p L_{R-1}, R = 2..399", ok)

log(""); log("="*78); log("(2) doubling law"); log("="*78)
ok = True
for R in (2,3,4,5,6,8,12,16,24,32,48,64,100,128,256):
    # law: G_{2R} = 2q G_R^2 + (p-q)(2G_R - 1) - p^R q^{R-1}.  With g := 2G_R
    # (integer for every R), multiply by 2:
    #   twoG(2R) = q g^2 + 2(p-q)(g - 1)/... : x4 then /2:
    #   2*twoG(2R) = 2q g^2 + 4(p-q)(g-1) - 4 p^R q^{R-1}
    g = twoG(R)
    q1 = [1, -1]; pmq = [-1, 2]
    lhs = pscale(twoG(2*R), 2)
    rhs = padd(padd(pscale(pmul(q1, pmul(g, g)), 2),
                    pscale(pmul(pmq, psub(g, [1])), 4)),
               pscale(pmul(pshift([1], R), qpow(R-1)), -4))
    if lhs != rhs: ok = False; break
check("doubling law G_{2R} = 2qG_R^2 + (p-q)(2G_R-1) - p^R q^{R-1} (x4 form), sample R<=256", ok)

log(""); log("="*78); log("(3) Lucas diagonal expansion (staircase basis, unique)"); log("="*78)
ok = True; okL = True
for R in list(range(2, 200)) + [256, 512, 1024]:
    # L_R = 1 + sum a_k e^k
    e = [0,1,-1]  # pq = x - x^2
    acc = [1]; epow = [1]
    for k in range(1, R//2 + 1):
        epow = pmul(epow, e)
        acc = padd(acc, pscale(epow, a_coef(R,k)))
    if acc != Lpoly(R): okL = False; break
    # 2G_R = 2 + sum a_k p^k q^{k-1}
    acc2 = [2]
    for k in range(1, R//2 + 1):
        acc2 = padd(acc2, pscale(pshift(qpow(k-1), k), a_coef(R,k)))
    if acc2 != twoG(R): ok = False; break
check("L_R = 1 + sum_k a_k e^k with a_k = (-1)^k[C(R-k,k)+C(R-k-1,k-1)], R<=1024 sample", okL)
check("2G_R = 2 + sum_k a_k p^k q^{k-1}  (staircase expansion), R<=1024 sample", ok)
# uniqueness: staircase atom p^k q^{k-1} has lowest degree term x^k -> triangular
ok = True
for k in range(1, 40):
    at = pshift(qpow(k-1), k)
    if any(at[j] != 0 for j in range(k)) or at[k] != 1: ok = False
check("staircase triangularity: p^k q^{k-1} = x^k + higher", ok)

log(""); log("="*78); log("(4) dyadic evenness of a_k"); log("="*78)
alleven = [R for R in range(2, 1025)
           if all(a_coef(R,k) % 2 == 0 for k in range(1, R//2 + 1))]
check("all a_k(R) even (k>=1) iff R = 2^t, R = 2..1024",
      alleven == [2,4,8,16,32,64,128,256,512,1024], str(alleven[:12]))

log(""); log("="*78); log("(5) closed forms at R = 2, 4 and the R = 4 closed-form witness"); log("="*78)
check("G_2 = q", G(2) == [1,-1])
check("G_4 = q^2 - p^3", G(4) == psub(pmul([1,-1],[1,-1]), pshift([1],3)))
check("G_4 = 1 - 2p + p^2 q  (diagonal form)",
      G(4) == padd(padd([1], pscale([0,1],-2)), pmul(pshift([1],2), [1,-1])))
# closed-form witness R=4: gamma_0 = q^2, gamma_1 = -p^2, gamma_2 = 0
w = padd(pmul([1,-1],[1,-1]), pshift(pscale(pshift([1],2), -1), 1))
check("R=4 witness: q^2 - x*p^2 = G_4 exactly", w == G(4))

log(""); log("="*78); log("(6) magnitude profile: |m_k| = |a_k|/2 vs capacities (dyadic R)"); log("="*78)
def floorgs(m):
    # exact floor(m gamma*), gamma* = log_5 phi^2 via Fib/Lucas
    def fp(n):
        if n == 0: return (0,1)
        a,b = fp(n>>1); c = a*(2*b-a); d2 = a*a+b*b
        return (d2, c+d2) if n & 1 else (c, d2)
    def le(d, m):
        if d < 0: return True
        F, F1 = fp(2*m); L = 2*F1 - F
        A = 2*5**d - L
        return A <= 0 or A*A < 5*F*F
    d = int(m*0.5979874356654402)
    while le(d+1, m): d += 1
    while not le(d, m): d -= 1
    return d
log("  R | argmax k | k/R | bits|m_k| | bits at cap of row k: d_k-1 | head m_1")
for t in range(3, 11):
    R = 2**t
    best, bk = 0, 0
    for k in range(1, R//2+1):
        v = abs(a_coef(R,k))//2
        if v > best: best, bk = v, k
    dk = floorgs(R + bk)
    log("  %5d | %6d | %.4f | %9d | %9d | m_1 = %d" %
        (R, bk, bk/R, best.bit_length(), dk-1, -(R//2)))
flush()

log(""); log("="*78); log("(7) x = 1 ledger"); log("="*78)
ok = True
for t in range(1, 11):
    R = 2**t
    g = G(R)
    if sum(g) != -(R-2)//2: ok = False
check("G_R(1) = -(R-2)/2 (all mass in head; atoms k>=2 vanish at 1), R dyadic <= 1024", ok)

log(""); log("SUMMARY: %s" % ("ALL CHECKS PASSED" if FAIL == 0 else "%d FAILURES" % FAIL))
flush()
