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
    best=F(0); arg=None
    for t in C:
        if 0<t<1:
            m=min(nf(s*t) for s in S)
            if m>best: best=m; arg=t
    return best,arg
print("SMALL covering sets (all speeds bounded) -- lonely (M>1/14)? 1/14=%.5f"%(1/14))
tests = {
 "{2..14}": list(range(2,15)),
 "{2,3,4,5,6,7,8,9,10,11,12,13,14}": list(range(2,15)),
 "even+odd cover {2,4,6,8,10,12,14,9,5,3,11,13,7}": [2,4,6,8,10,12,14,9,5,3,11,13,7],
 "{4,6,8,9,10,12,14,5,7,11,13,2,3}": [4,6,8,9,10,12,14,5,7,11,13,2,3],
}
for name,S in tests.items():
    cov=all(any(s%q==0 for s in S) for q in range(2,15))
    M,arg=M_exact(S)
    print(f"  {name:42s}: cov={cov}, M={str(M):>8s}={float(M):.5f} {'LONELY' if M>F(1,14) else 'TIGHT/below'} t={arg}")
print("\n=> If small covering sets are all LONELY (M>1/14 with margin), the covering crux splits:")
print("   BOUNDED covering (finite/compactness, THM-527) + UNBOUNDED (equidistribution). Neither is tight.")
