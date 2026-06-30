"""
THREAD B: Ky Fan forced-existence for the LRC lonely point.
Danger labeling d(t)=(1[||s_i t||<1/14])_i traces a walk on Q_{n-1}; lonely = {d=0}. t->-t is even
(d(-t)=d(t)) => lonely set is MIRROR-symmetric => lonely intervals come in pairs => COUNT EVEN
(Borsuk-Ulam regime; Redei's odd-forcing FAILS -- the crux of LRC hardness).
FORCING TEST: factor out the '2' (the mirror, 14=2*7). Is the HALF-circle (0,1/2) lonely-interval
count ODD (=> total = 2*odd > 0, forced)? Compute parity for covering sets.
"""
from fractions import Fraction as F
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def lonely_intervals(S,n=14,lo=F(0),hi=F(1,2)):
    # boundaries where some ||s_i t||=1/n: t=(n k +-1)/(n s_i)
    B=set([lo,hi])
    for s in S:
        k=0
        while True:
            for t in (F(n*k+1,n*s), F(n*k-1,n*s)):
                if lo<t<hi: B.add(t)
            if F(n*k-1,n*s) >= hi and k>0: break
            k+=1
            if k> n*int(hi)+ s+5: break
    B=sorted(B)
    cnt=0; intervals=[]
    for i in range(len(B)-1):
        mid=(B[i]+B[i+1])/2
        if min(nrm(s*mid) for s in S)>=F(1,n):
            cnt+=1; intervals.append((float(B[i]),float(B[i+1])))
    return cnt,intervals
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
import random
random.seed(3)
sets={"{2..14}":list(range(2,15)),
      "{1..11,13,84}":[1,2,3,4,5,6,7,8,9,10,11,13,84],
      "{1..11,13,28}":[1,2,3,4,5,6,7,8,9,10,11,13,28],
      "{1..12,182}":[1,2,3,4,5,6,7,8,9,10,11,12,182]}
# add random covering
cnt=0
while len(sets)<12 and cnt<50000:
    cnt+=1; S=tuple(sorted(random.sample(range(1,30),13)))
    if is_cov(S): sets["rand "+str(S[:3])+"..."]=list(S)
print("LRC half-circle (0,1/2) lonely-interval count parity (is it ODD => total=2*odd>0 forced?):")
odd=0;tot=0
for nm,S in sets.items():
    h,iv=lonely_intervals(S)
    full,_=lonely_intervals(S,hi=F(1))
    par="ODD" if h%2==1 else "even"
    print(f"  {nm:28s}: half-count={h} ({par}), full-count={full}  full mod4={full%4}")
    tot+=1; odd+= (h%2)
print(f"\n  {odd}/{tot} covering sets have ODD half-count.  If ALWAYS odd => Ky-Fan forces LRC on the 7-core.")
