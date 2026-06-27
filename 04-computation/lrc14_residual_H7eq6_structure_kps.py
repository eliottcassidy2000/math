"""
lrc14_residual_H7eq6_structure_kps.py  (kind-pasteur-2026-06-27-S31ag)

LEAD 8: THM-573/574 close |H_7| >= 7. The residual is |H_7| <= 6. The boundary
|H_7| = 6 is where the c=7 lift JUST fails the survivor count (7 coprime speeds,
7 lifts). It fails only if those 7 coprime speeds kill 7 DISTINCT lifts (a perfect
kill-system); any overlap leaves a survivor.

Test: over primitive covering 13-sets with EXACTLY |H_7| = 6:
  - how often is M > 1/14 anyway (lift survives / overlap)?
  - the min M, and whether tight/near-tight ones are dilated-AP / GW;
  - directly: count the "perfect kill-system" failures.
"""
from math import gcd
import random

def nrm(x):
    f = x - int(x)
    if f < 0: f += 1.0
    return min(f, 1 - f)

def M_value(S):
    S = sorted(set(S)); denoms = set()
    for i in range(len(S)):
        a = S[i]; denoms.add(2 * a)
        for j in range(len(S)):
            if i == j: continue
            b = S[j]; denoms.add(a + b)
            if a != b: denoms.add(abs(a - b))
    best = -1.0
    for D in denoms:
        if D == 0: continue
        for n in range(1, D):
            t = n / D; m = 1.0
            for s in S:
                v = nrm(s * t)
                if v < m: m = v
                if m <= best: break
            if m > best: best = m
    return best

def gcdl(S):
    g = 0
    for s in S: g = gcd(g, s)
    return g

def kill_system_perfect(S):
    """For |H_7|=6: is there a P-safe phase v* whose 7 lifts t_j=(v*+j)/7 are ALL killed
    by the 7 coprime speeds (each killing a distinct lift)? If NO safe v* exists or some
    lift always survives, return False (residual is benign). Heuristic over v* candidates."""
    H = [s for s in S if s % 7 == 0]
    if len(H) != 6: return None
    P = [h // 7 for h in H]
    coprime = [s for s in S if s % 7 != 0]   # 7 of them
    # P-safe phase: need ||p v*|| >= 1/14 for all p in P. Search rational v* = n/D.
    # use the same denom candidates restricted to P
    Pden = set()
    for i in range(len(P)):
        a = P[i]; Pden.add(2*a)
        for j in range(len(P)):
            if i==j: continue
            b=P[j]; Pden.add(a+b)
            if a!=b: Pden.add(abs(a-b))
    # find any v* maximizing min_p ||p v*||, then test whether a lift survives the coprime kills
    best_survive = False
    for D in Pden:
        if D==0: continue
        for n in range(1, D):
            vstar = n / D
            if min(nrm(p*vstar) for p in P) < 1/14 - 1e-12:
                continue
            # 7 lifts; a lift t_j survives if all coprime speeds safe there
            for j in range(7):
                tj = (vstar + j)/7
                if all(nrm(w*tj) >= 1/14 - 1e-9 for w in coprime):
                    best_survive = True; break
            if best_survive: break
        if best_survive: break
    return best_survive  # True = a lift survives (benign); False = no surviving lift found

if __name__ == "__main__":
    random.seed(706)
    thr = 1/14
    n_total = 0; n_tight = 0; n_gt = 0; minM = 1.0
    tight_examples = []; lowM = []
    survive_true = 0; survive_false = 0
    for trial in range(3000):
        # 6 multiples of 7 (at least one multiple of 14), 7 coprime-to-7
        mult14 = random.choice([14, 28, 42, 56, 70])
        m7 = set([mult14])
        while len(m7) < 6:
            m7.add(7 * random.randint(1, 12))
        cop = set()
        while len(cop) < 7:
            x = random.randint(1, 40)
            if x % 7 != 0: cop.add(x)
        S = sorted(m7 | cop)
        if len(S) != 13: continue
        if gcdl(S) != 1: continue
        if sum(1 for s in S if s % 7 == 0) != 6: continue
        n_total += 1
        M = M_value(S)
        if M < minM: minM = M
        if M > thr + 1e-9: n_gt += 1
        if abs(M - thr) < 1e-6:
            n_tight += 1
            if len(tight_examples) < 8: tight_examples.append(S)
        if M < thr + 0.01:
            lowM.append((round(M,5), S))
    print(f"Primitive covering 13-sets with |H_7|=6: tested {n_total}")
    print(f"  M > 1/14 (survives, benign): {n_gt}/{n_total} = {100*n_gt/max(n_total,1):.1f}%")
    print(f"  tight (M=1/14): {n_tight}   minM = {minM:.5f}")
    if tight_examples:
        print("  tight examples (check dilated-AP/GW):")
        for S in tight_examples:
            g = gcdl(S)
            print(f"    {S}  gcd={g}  reduced={[s//g for s in S]}")
    lowM.sort()
    print("  lowest-M sets:")
    for M, S in lowM[:6]:
        print(f"    M={M}  {S}")
