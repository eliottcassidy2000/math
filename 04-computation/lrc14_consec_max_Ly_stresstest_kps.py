"""
lrc14_consec_max_Ly_stresstest_kps.py  (kind-pasteur-2026-06-27-S31ag)

Stress-test the load-bearing assumption of the binding-row cover bound:
  "consec maximizes L_y" for k=9,10 (degree-3 dual) and k=8 (degree-4).
The valid Bonferroni duals (w_t >= 1[t=0]) make meas(S7)=q_0 <= L_y(E) ALWAYS;
the bound L_y(E) <= cap_k then needs consec to be the MAXIMIZER of L_y.
consec-max is FALSE at k=12 (HYP-2780). Does it hold at the BINDING k=8,9,10?

Search a diverse pool (bounded, wide single/double-far, dilated, cap-minimizers,
random) for any config with L_y(E) > L_y(consec_k). A counterexample breaks the
route; robust consec-max confirms it.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb

def sector_of(p): return int((p % 1) * 7)

def N_dist(E):
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        cov = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        t = 7 - len(cov)
        if 0 <= t <= 6: q[t] += x1 - x0
    return q

# valid Bonferroni duals in q-readout form: w_t weights, t=0..6
WEIGHTS = {
    8:  [10, 0, 0, 1, 0, 0, 10],          # L_yK8 = 10q0 + q3 + 10q6
    9:  [18, 5, 0, 0, 2, 3, 0],           # L_yK9 = 18q0+5q1+2q4+3q5
    10: [18, 5, 0, 0, 2, 3, 0],
}
CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def Ly(E, k):
    q = N_dist(E); w = WEIGHTS[k]
    return sum(w[t] * q[t] for t in range(7))

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(910)
    for k in (8, 9, 10):
        consec = tuple(range(k))
        Lc = Ly(consec, k); cap = CAPS[k]
        pool = [("consec", consec)]
        # structured competitors -- ALL must contain 0 (the stationary anchor; sector 0 auto-hit)
        pool.append(("even-AP", tuple(2*i for i in range(k))))
        pool.append(("single-far", consec[:-1] + (21,)))
        pool.append(("single-far2", consec[:-1] + (40,)))
        pool.append(("doublet", consec[:-2] + (20, 21)))
        pool.append(("cap-min", tuple(sorted({0,1,5,7,8,9,10,11,12,13})) [:k]))
        # random bounded, wide, far
        for _ in range(2000):
            mode = random.random()
            if mode < 0.4:
                cfg = tuple(sorted([0] + random.sample(range(1, 16), k-1)))
            elif mode < 0.7:
                cfg = tuple(sorted([0] + random.sample(range(1, 30), k-1)))
            else:
                base = sorted(random.sample(range(1, 12), k-2))
                cfg = tuple(sorted([0] + base + [random.randint(20, 60)]))
            if len(set(cfg)) == k and 0 in cfg:
                pool.append(("rand", cfg))
        best = Lc; bestcfg = consec; nbeat = 0; beaters = []
        for name, E in pool:
            L = Ly(E, k)
            if L > Lc + F(1, 10**12):
                nbeat += 1
                if L > best: best = L; bestcfg = E
                if len(beaters) < 4: beaters.append((float(L), name, E))
        print(f"k={k}: L_y(consec)={float(Lc):.5f}  cap={float(cap):.5f}  "
              f"consec<=cap? {Lc <= cap}")
        print(f"   configs BEATING consec's L_y: {nbeat}/{len(pool)}  "
              f"{'consec is MAX (route holds)' if nbeat==0 else 'CONSEC NOT MAX -> '+str(beaters[:3])}")
        if nbeat:
            print(f"   max L_y = {float(best):.5f} at {bestcfg}  (> cap? {best > cap})")
