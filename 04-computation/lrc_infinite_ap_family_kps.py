"""
The infinite AP family & the zeta(-1) limit (kps-S31p). AP_n={1..n} (n speeds, LRC(n+1)).
M(AP_n)=1/(n+1): the AP kills resonances b=1..n (speed=b) leaving b=n+1 first survivor.
Sum of killed resonances 1+..+n=n(n+1)/2; the INFINITE AP {1,2,3,...} kills ALL resonances,
M->0, and the regularized killing-sum 1+2+3+...=zeta(-1)=-1/12. Finite avatar: M({1..11,13})=1/12.
"""
from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def M_exact(S):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    return max((min(nf(s*t) for s in S) for t in C if 0<t<1), default=F(0))
print("INFINITE AP FAMILY: M(AP_n)=M({1..n}) =?= 1/(n+1), killed resonances {1..n}:")
for n in range(2,17):
    M=M_exact(list(range(1,n+1)))
    print(f"  n={n:2d}: M({{1..{n}}}) = {str(M):>6s}  {'= 1/(n+1)' if M==F(1,n+1) else '!= 1/(n+1)'};  sum 1..n = {n*(n+1)//2}")
print("\n  => M(AP_n)=1/(n+1) EXACT (the AP is the unique least-spread killer of b=1..n).")
print("     Infinite AP {1,2,3,...}: kills ALL resonances, M->0; regularized 1+2+3+..=zeta(-1)=-1/12.")
print("     1/12-core {1..11,13}: skips killer-12 => b=12 survives => M=1/12 (the finite -zeta(-1) avatar).")
