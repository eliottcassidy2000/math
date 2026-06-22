"""
hp0cap residual analysis via the L_y route (mac-mini-2026-06-22-S29, HYP-2852).

Findings (k=9 dual unless noted):
1. ADDITIVE-ENERGY FRAME (approximate): consec maximizes both L_y and additive energy
   A(E)=Σ_s r_E(s)², but L_y is NOT a function of A (same-A configs differ; envelope not
   monotone). APs maximize A (classical) -- the rigorous half; the L_y↑A bridge is only a trend.
2. MOMENT STRUCTURE: consec MIN's S_3,S_4 but is mid for S_1,S_2; the dual L_y is intrinsically
   COUPLED (g=[1,0,0,0.1,0,0,1] for k=8 => L_y≈p_0, the cover extremality, genuinely hard).
3. CLOSURE: bounded span<=14 DONE (HYP-2830); WIDE decomposes by #far, ALL CLOSE
   (1far margin 0.035 binding @ consec_8+far21; 2far 0.112).
4. SQRT-CANCELLATION (the key): peel deviation |Δ_w·w| ~ C·√V (NOT (6/49)V), via the
   mean-zero sawtooth F_j; avg_w|Δ_w·w|² ≤ c·V (Parseval). Explains THM-546's 5-8x looseness.
   Generic w => √-cancellation (Parseval, mine); resonant w (max outliers) => kps HYP-2842.
"""
from fractions import Fraction as Fr
from math import comb
from collections import Counter
import random

def cells_frac(E):
    E = sorted(set(int(e) for e in E)); bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7*abs(e)+1): bp.add(Fr(m, 7*abs(e)))
    bp = sorted(bp); out = []
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        hit = set(int((e*((a+b)/2)) % 1 * 7) for e in E)
        out.append((b-a, 6 - len(hit & {1,2,3,4,5,6})))
    return out

YK = {8: [Fr(1),Fr(-1),Fr(1),Fr(-9,10),Fr(3,5)], 9: [Fr(1),Fr(-13,18),Fr(4,9),Fr(-1,6)]}
def gval(k, t):
    return sum(YK[k][r]*comb(t, r) for r in range(len(YK[k])))
def Ly(E, k): return sum(L*gval(k, N) for L, N in cells_frac(E))
def energy(E):
    c = Counter()
    for a in E:
        for b in E: c[a+b] += 1
    return sum(v*v for v in c.values())

# float versions for fast search
def Ly_f(E, k):
    g = [float(gval(k, t)) for t in range(7)]
    return sum(float(L)*g[N] for L, N in cells_frac(E))

if __name__ == "__main__":
    k = 9; cap = Fr(1979, 4004)
    print("consec_9: A=%d L_y=%.5f cap=%.5f (margin %.5f)" %
          (energy(list(range(9))), float(Ly(list(range(9)), 9)), float(cap),
           float(cap - Ly(list(range(9)), 9))))
    # 1far / 2far closure
    random.seed(5)
    best1 = 0
    for _ in range(800):
        sp = random.randint(7, 14); C = [0]+sorted(random.sample(range(1, sp+1), 7))
        for w in range(15, 60):
            best1 = max(best1, Ly_f(C+[w], 9))
    print("1far max L_y = %.4f (margin %.4f) -> closes: %s" %
          (best1, float(cap)-best1, best1 < float(cap)))
    # sqrt-cancellation RMS over w
    import statistics
    for name, C in [("consec_8", list(range(8))), ("AP_d4", [4*i for i in range(8)])]:
        cs = cells_frac(C)
        Ld = sum(float(L)*((N/7)*([1,5/18,0,0,1/9,1/6,0][N-1] if N>=1 else 0)
                 + (1-N/7)*[1,5/18,0,0,1/9,1/6,0][N]) for L, N in cs)
        vals = [((Ly_f(C+[w], 9)-Ld)*w)**2 for w in range(max(C)+1, max(C)+200)]
        print("  %s: #cells=%d RMS|Δw·w|=%.3f (~√cells => sqrt-cancellation)" %
              (name, len(cs), statistics.mean(vals)**0.5))
