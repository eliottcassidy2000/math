"""
Verify klein THM-587 (per-level SIGNED CYCLE INDEX) and EXTEND the SC sequence past n=8.
P_n(x) = (1/n!) sum_{sigma in S_n} prod_{pair-cycles c} (1 + s_c x^{ell_c}),
 s_c = -1 if sigma^{ell_c} REVERSES the pair (orientation-reversing), else +1.
 P_n(1) = A000568 (tournaments);  P_n(-1) = SC (self-converse / antipodal Euler number).
Also: NS = (A000568 - SC)/2 = R-odd block (the Borsuk-Ulam obstruction dimension).
"""
from itertools import permutations, combinations
from fractions import Fraction as F

def Pn_eval(n):
    pairs=list(combinations(range(n),2)); idx={p:k for k,p in enumerate(pairs)}; m=len(pairs)
    tot1=0; totm1=0   # accumulate n!*P_n(1), n!*P_n(-1)
    for s in permutations(range(n)):
        # pair-permutation
        perm=[0]*m
        for k,(i,j) in enumerate(pairs):
            a,b=s[i],s[j]; perm[k]=idx[(a,b) if a<b else (b,a)]
        seen=[False]*m; prod1=1; prodm1=1
        for start in range(m):
            if seen[start]: continue
            # follow cycle
            cyc=[]; x=start
            while not seen[x]:
                seen[x]=True; cyc.append(x); x=perm[x]
            L=len(cyc)
            # s_c: does sigma^L reverse pair 'start'? track element order
            i,j=pairs[start]
            # apply s, L times to (i,j)
            ii,jj=i,j
            for _ in range(L): ii,jj=s[ii],s[jj]
            sc = 1 if (ii,jj)==(i,j) else -1   # else it's (j,i) (reversed)
            # factor (1 + sc * x^L) at x=1 and x=-1
            prod1 *= (1 + sc*(1**L))
            prodm1*= (1 + sc*((-1)**L))
        tot1+=prod1; totm1+=prodm1
    import math
    nf=math.factorial(n)
    return tot1//nf, totm1//nf   # both integers

print(f"{'n':>2} {'P_n(1)=A000568':>15} {'P_n(-1)=SC':>11} {'NS=(A-SC)/2':>12} {'V_merged=(A+SC)/2':>17}")
A000568=[1,1,2,4,12,56,456,6880,191536]
for n in range(3,10):
    a,sc=Pn_eval(n)
    ns=(a-sc)//2; vm=(a+sc)//2
    chk = "OK" if a==A000568[n-1] else f"!=A000568({A000568[n-1]})"
    print(f"{n:>2} {a:>15} {sc:>11} {ns:>12} {vm:>17}   {chk}")
print()
print("klein THM-587: SC = 2,2,8,12,88,176 (n=3..8). Extension below if n=9 computed.")
print("NS = R-odd block dim = where the Borsuk-Ulam/Ky-Fan obstruction lives (klein: the LRC odd index).")
