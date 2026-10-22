# Test the higher Kummer congruences for strong tournaments:
#   sum_{j=0}^{2k-1} C(2k-1,j) (-eps)^(2k-1-j) s(m + j(p-1)) = 0 (mod p^k)  for m >= 2k-1,
# and report the failing m below the threshold.  Exact s(m) mod p^k via Wright's recurrence.
from math import comb
def smod(M,q):
    T=[pow(2,j*(j-1)//2,q) for j in range(M+1)]; s=[0]*(M+1); row=[1]
    for n in range(1,M+1):
        row=[1]+[(row[i-1]+row[i])%q for i in range(1,n)]+[1]
        s[n]=(T[n]-sum(row[k]*s[k]*T[n-k] for k in range(1,n)))%q
    return s
for p,k in [(5,1),(5,2),(7,2),(11,2),(13,2),(5,3),(7,3),(11,3),(13,3),(17,3),(7,4),(11,4),(5,5)]:
    q=p**k; e=1 if p%8 in (1,7) else -1; r=2*k-1; M=40+r*(p-1)+60
    s=smod(M,q)
    fails=[m for m in range(1,M-r*(p-1)+1) if sum(comb(r,j)*(-e)**(r-j)*s[m+j*(p-1)] for j in range(r+1))%q]
    print(f'p={p} k={k}: order-{r} difference with step p-1: failures at m = {fails}  (threshold 2k-1 = {r})')
