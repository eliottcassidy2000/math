from fractions import Fraction as F
P=[1,2,4,7,9,11,12]
def danger(v):
    w=F(1,14*v);o=[]
    for j in range(v):
        c=F(j,v);lo=(c-w)%1;hi=(c+w)%1
        if lo<hi:o.append((lo,hi))
        else:o.append((lo,F(1)));o.append((F(0),hi))
    return o
def subf(S,a):
    for clo,chi in sorted(a):
        n=[]
        for lo,hi in S:
            if chi<=lo or clo>=hi:n.append((lo,hi));continue
            if clo>lo:n.append((lo,clo))
            if chi<hi:n.append((chi,hi))
        S=n
    return S
S0=[(F(0),F(1))]
for v in P:S0=subf(S0,danger(v))
def L(K):
    S=S0
    for k in K:S=subf(S,danger(k))
    return max(hi-lo for lo,hi in S) if S else F(0)
def R(K):
    l=L(K);return F(1,7)/(l*max(K))
# EXACT per-k5 counterexample to "consecutive is extremal": k5=165
consec=[157,159,161,163,165]; nonc=[159,161,162,163,165]
print("k5=165  consecutive %s: R=%s=%.5f"%(consec,R(consec),float(R(consec))))
print("k5=165  NON-consec   %s: R=%s=%.5f  <-- BEATS consecutive"%(nonc,R(nonc),float(R(nonc))))
print("  => for fixed k5=165, a non-consecutive shape has strictly larger R_sharp (exact).")
print()
# and the global maximizer is consecutive:
gmax=[171,173,175,177,179]
print("global maximizer %s: R=%s=%.6f (the record)"%(gmax,R(gmax),float(R(gmax))))
