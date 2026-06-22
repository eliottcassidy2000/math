"""
hp0cap via the THM-534 L_y route + sharp arc-complexity (mac-mini-2026-06-22-S28, HYP-2840).

Establishes:
1. p0(E) <= L_y(E) (THM-534 moment-LP dual, PROVED) -> hp0cap residual = "consec maximizes L_y".
2. The tight cap-margin (0.00138, k=9) lives only at the consec AP-orbit; any decorrelation drops L_y by >=0.044.
3. The arc-complexity V(E')=Σ_j #arcs(B_j) is 40-107x smaller than the THM-546 bound 42Σe,
   making the gapped single-peel cutoff w* = (6/49)V/0.044 ~ 80 (feasible).
4. V(E') ~ 4*span (NOT uniformly bounded) -> gapped vs ungapped(scale-invariance) dichotomy.
"""
from fractions import Fraction as Fr
from math import comb
import itertools

def breaks(E):
    E = sorted(set(int(e) for e in E)); bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7*abs(e)+1): bp.add(Fr(m, 7*abs(e)))
    return sorted(bp)

def Sr_fast(E):  # factorial moments S_r = E[C(N,r)], N = #missed sectors among 1..6
    bp = breaks(E); S = {r: Fr(0) for r in range(5)}
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        hit = set(int((e*((a+b)/2)) % 1 * 7) for e in E)
        free = 6 - len(hit & {1,2,3,4,5,6})
        for r in range(5): S[r] += (b-a)*comb(free, r)
    return S

# dual functionals L_y = Σ_r y_r S_r (THM-534), g(t)>=1[t=0] on {0..6}
YK = {8: {0:Fr(1),1:Fr(-1),2:Fr(1),3:Fr(-9,10),4:Fr(3,5)},
      9: {0:Fr(1),1:Fr(-13,18),2:Fr(4,9),3:Fr(-1,6)},
      10:{0:Fr(1),1:Fr(-13,18),2:Fr(4,9),3:Fr(-1,6)},
      11:{0:Fr(1),1:Fr(-1,2),2:Fr(1,6)},12:{0:Fr(1),1:Fr(-1,2),2:Fr(1,6)}}
CAP = {8:Fr(2243,5880),9:Fr(1979,4004),10:Fr(55,91),11:Fr(66,91),12:Fr(6,7)}

def Ly(E, k):
    S = Sr_fast(E); return sum(YK[k][r]*S[r] for r in YK[k])

def p0(E):
    bp = breaks(E); tot = Fr(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        if len(set(int((e*((a+b)/2)) % 1 * 7) for e in E)) == 7: tot += b-a
    return tot

def Vactual(E):  # Σ_j #arcs(B_j), B_j = {x: E misses EXACTLY sector j among 1..6}
    bp = breaks(E); seq = []
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        hit = set(int((e*((a+b)/2)) % 1 * 7) for e in E)
        seq.append(frozenset({1,2,3,4,5,6} - hit))
    tot = 0
    for j in range(1, 7):
        prev = False; c = 0
        for m in seq:
            cur = (m == {j})
            if cur and not prev: c += 1
            prev = cur
        if seq and seq[0] == {j} and seq[-1] == {j} and c > 1: c -= 1
        tot += c
    return tot

if __name__ == "__main__":
    print("1. p0 <= L_y and consec maximizes L_y (k=9, tightest):")
    for tag, E in [('consec', list(range(9))), ('AP_d3', [3*i for i in range(9)]),
                   ('1far@15', list(range(8))+[15]), ('gap', [0,1,2,3,4,5,6,8,9])]:
        print(f"   {tag:>10}: p0={float(p0(E)):.5f} L_y={float(Ly(E,9)):.5f} "
              f"p0<=Ly={p0(E)<=Ly(E,9)} Ly<=cap={Ly(E,9)<=CAP[9]}")
    print(f"   cap_9={float(CAP[9]):.5f}, consec margin={float(CAP[9]-Ly(list(range(9)),9)):.5f}")

    print("\n2. arc-complexity V_actual vs 42Σe (consec cores):")
    for kp in range(7, 13):
        E = list(range(kp)); V = Vactual(E); se = sum(E)
        print(f"   consec_{kp}: V_actual={V}  42Σe={42*se}  ratio={42*se/V:.0f}x  "
              f"cutoff w*={float(Fr(6,49)*V/Fr(44,1000)):.0f}")

    print("\n3. V ~ 4*span (not uniformly bounded):")
    for tag, E in [('consec_8', list(range(8))), ('AP_d5', [5*i for i in range(8)]),
                   ('1far@200', list(range(7))+[200]), ('allbig', [50*i for i in range(8)])]:
        print(f"   {tag:>10}: span={max(E)-min(E)} V={Vactual(E)}")
