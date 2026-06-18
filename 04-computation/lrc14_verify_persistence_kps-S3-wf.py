"""
Persistence test (FIXED covering): does C(S) keep holding as the large cluster scale grows?

Construction that STAYS covering for all scales t:
  small part P subset {1..13} chosen to cover {2,3,...,13} (all moduli except 14),
  large cluster L = {14*t + d : d in Delta} with 0 in Delta so 14*t covers modulus 14.
As t -> infinity, Vmax = 14*t + max(Delta) -> infinity, but covering holds for ALL t.

For each (P, Delta) we sweep t and record best-C-margin via Wwidth (cheap-ish; cost ~ sum speeds,
so we cap t). We test whether the margin stays >1 (C holds) and whether the witness runner
and margin stabilize (the author's 'persistence' claim: #good ruler periods ~ linear in w0,
margin bounded below by a positive constant).

This directly probes GAP G1 (universality / uniform positive floor) in the dangerous direction:
arbitrarily large Vmax with a FIXED offset pattern.

kind-pasteur-S3-wf
"""
from fractions import Fraction as F
from math import gcd
import itertools, sys, time

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g

def case_of(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13*Vmin: return 'S2'
    return 'S3'

def best_C(S):
    best = F(0); arg = None
    for v in S:
        A = [x for x in S if x != v]
        r = Wwidth(A) * 7 * v
        if r > best: best = r; arg = v
    return best, arg

def covers_2_to_13(P):
    return all(any(v % q == 0 for v in P) for q in range(2, 14))

def main():
    print("PERSISTENCE (fixed-covering) — C(S) margin vs large-cluster scale t")
    # small parts that cover 2..13. e.g. {4,6,7,9,10,11,12,13} covers? 2|4,6,..;3|6,9,12;5|10;... check.
    # We'll search small subsets of size 8..11 of {1..13} covering 2..13.
    smalls = []
    for sz in (8, 9, 10, 11):
        for P in itertools.combinations(range(1, 14), sz):
            if covers_2_to_13(P):
                smalls.append(list(P))
        if len(smalls) >= 6: break
    # Delta offset patterns (include 0). nlarge = 13 - |P|. We pick a few.
    patterns = {
        'AP1':   lambda k: [i for i in range(k)],          # 0,1,2,...
        'AP2':   lambda k: [2*i for i in range(k)],        # 0,2,4,...
        'AP3':   lambda k: [3*i for i in range(k)],
        'mix':   lambda k: [0,3,9,12,15][:k],
        'wide':  lambda k: [0,5,17,31,45][:k],
    }
    overall_min = F(100); overall_arg = None
    for P in smalls[:4]:
        nlarge = 13 - len(P)
        if nlarge < 2: continue
        for pname, pf in patterns.items():
            Delta = pf(nlarge)
            if len(set(Delta)) != nlarge: continue
            margins = []
            for t in [1, 2, 3, 5, 10, 20, 40, 80, 160]:
                anchor = 14*t
                large = [anchor + d for d in Delta]
                S = sorted(set(P) | set(large))
                if len(S) != 13: continue
                if gcd_all(S) != 1: continue
                if not is_covering(S): continue
                if case_of(S) != 'S3': continue
                bc, arg = best_C(S)
                margins.append((t, max(S), bc, arg))
                if bc < overall_min: overall_min = bc; overall_arg = (P, Delta, t, S)
            if margins:
                tail = margins[-1]
                print(f"  P={P} {pname} Delta={Delta}: "
                      f"margins(t,Vmax,marg,v)= " +
                      ", ".join(f"({t},{vm},{float(bc):.2f},v{arg})" for t,vm,bc,arg in margins))
                # does C hold for ALL tested t?
                allhold = all(bc > 1 for _,_,bc,_ in margins)
                print(f"      C holds for all tested t: {allhold}; min margin = {float(min(bc for _,_,bc,_ in margins)):.3f}")
                sys.stdout.flush()
    print(f"\n  OVERALL min best-C-margin over persistence sweep = {float(overall_min):.4f}")
    if overall_arg:
        P, Delta, t, S = overall_arg
        print(f"   at P={P} Delta={Delta} t={t}  S={S}")

if __name__ == '__main__':
    main()
