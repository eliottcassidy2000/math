# D(m) = s(m) - 2^((m-1)/2) mod p for odd m: zero sets and periods.  s = strong labeled tournaments.
from math import comb
def s_mod(M,p):
    T=[pow(2,j*(j-1)//2,p) for j in range(M+1)]
    s=[0]*(M+1)
    # binomials mod p by rows
    row=[1]
    for k in range(1,M+1):
        row=[1]+[(row[i-1]+row[i])%p for i in range(1,k)]+[1]
        acc=0
        for i in range(1,k): acc=(acc+row[i]*s[i]*T[k-i])%p
        s[k]=(T[k]-acc)%p
    return s
M=400
for p in [3,5,7,11,13,17,19,23,29,31,37,41,43]:
    s=s_mod(M,p)
    D=[None]*(M+1)
    for m in range(1,M+1,2): D[m]=(s[m]-pow(2,(m-1)//2,p))%p
    zeros=[m for m in range(1,M+1,2) if D[m]==0]
    # find minimal period L (in m, step over odd m) with D[m+L]==D[m] for all odd m in range
    per=None
    for L in range(2,2*p*p+2,2):
        if all(D[m]==D[m+L] for m in range(1,M+1-L,2)): per=L; break
    # even m: s(m) mod p
    ev=[s[m] for m in range(2,40,2)]
    print(f'p={p}: period(odd m)={per}  p-1={p-1}  zeros(odd m<=80)={[m for m in zeros if m<=80]}  even s mod p (m=2..38)={ev}')
