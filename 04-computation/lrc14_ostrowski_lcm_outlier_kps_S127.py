# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont55: WHY the covering-min outlier is (N-1)N -- the lcm-outlier mechanism.
#
# The covering-min family {1,..,N-2, (N-1)N} (klein's deep well {1..12,182} at N=14) has ONE large outlier.
# WHY (N-1)N exactly? Because for the base {1,..,N-2} (which already covers d=2..N-2), the two MISSING
# divisors are d=N-1 and d=N, and the SMALLEST single speed carrying BOTH is lcm(N-1,N) = (N-1)N
# (gcd(N-1,N)=1). So {1..N-2,(N-1)N} is the minimal-outlier single-outlier covering family, and on the
# ladder M_k=k/((N-1)k+1) the outlier (N-1)N is rung k=N -- the FIRST covering rung, the ladder's min.
# This explains klein's 182 = 13*14 = lcm(13,14) structurally, and why the covering-min sits at N/Phi6(N).
from math import gcd
from fractions import Fraction as F

def lcm(a, b): return a * b // gcd(a, b)
def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mexact(v):
    qcap = 3 * max(v) + 2; best = F(0)
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best = m
    return best
def is_cov(v, N): return all(any(x % d == 0 for x in v) for d in range(2, N + 1))
def missing_div(v, N): return [d for d in range(2, N + 1) if not any(x % d == 0 for x in v)]

def main():
    print("THE lcm-OUTLIER MECHANISM: base {1..N-2} covers d=2..N-2; the two missing divisors are d=N-1, d=N;")
    print("the minimal single speed carrying BOTH is lcm(N-1,N) = (N-1)N. So the min single-outlier covering")
    print("family is {1..N-2, (N-1)N}, and on the ladder that outlier is rung k=N = the first covering rung.\n")
    print(f"{'N':>2} | base {{1..N-2}} missing | lcm(N-1,N)=(N-1)N | family {{1..N-2,(N-1)N}} covering? | M | =N/Phi6(N)?")
    for N in range(3, 15):
        base = list(range(1, N - 1))
        miss = missing_div(base, N)
        outlier = lcm(N - 1, N)
        fam = base + [outlier]
        cov = is_cov(fam, N)
        formula = F(N, N * N - N + 1)
        M = Mexact(fam) if outlier < 400 else None
        okM = (M == formula) if M is not None else "skip(big)"
        print(f"{N:>2} | {str(miss):>18} | {str((N-1)*N)+'='+str(outlier):>16} | {str(cov):>28} | {str(M):>7} | {okM}")
    print()
    print("CHECK the minimality: is (N-1)N the SMALLEST outlier making {1..N-2, w} covering? (need mult of N-1 AND N)")
    print(f"{'N':>2} | smallest w with {{1..N-2,w}} covering | = (N-1)N? | (any single w<(N-1)N covering?)")
    for N in range(3, 15):
        base = list(range(1, N - 1))
        found = None
        for w in range(1, (N - 1) * N + 1):
            if is_cov(base + [w], N): found = w; break
        print(f"{N:>2} | {str(found):>36} | {str(found==(N-1)*N):>9} | {'no single w below (N-1)N covers both N-1 and N' if found==(N-1)*N else 'SMALLER OUTLIER EXISTS: '+str(found)}")
    print()
    print("=> The covering-min single-outlier family is {1..N-2, lcm(N-1,N)} = {1..N-2,(N-1)N}, the first")
    print("   covering rung k=N of the ladder, M=N/Phi6(N). klein's 182=13*14=lcm(13,14) is the N=14 instance.")
    print("   (Whether a MULTI-outlier / compressed family beats it: small N yes {2,3}/{1,3,4}; N=14 no, klein ILP.)")

if __name__ == "__main__":
    main()
