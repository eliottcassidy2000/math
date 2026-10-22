# Independent (pure Python) recomputation of failure sets of s(m + phi(p^k)) = (2/p) s(m) (mod p^k), 1 <= m <= 2k+2,
# for a sample of (p,k); compares with the C table in HYP-9358(e).
def s_upto(N,q):
    T=[pow(2,j*(j-1)//2,q) for j in range(N+1)]; s=[0]*(N+1); row=[1]
    for n in range(1,N+1):
        row=[1]+[(row[i-1]+row[i])%q for i in range(1,n)]+[1]
        s[n]=(T[n]-sum(row[k]*s[k]*T[n-k] for k in range(1,n)))%q
    return s
expected={(3,4):5,(3,5):7,(3,6):9,(5,3):4,(5,4):5,(7,2):2,(7,3):5,(11,3):5,(13,2):2,(13,3):5}
for (p,k),tau in expected.items():
    q=p**k; L=p**(k-1)*(p-1); eps=1 if p%8 in (1,7) else -1; M=L+2*k+3
    s=s_upto(M,q)
    fails=[m for m in range(1,2*k+3) if (s[m+L]-eps*s[m])%q]
    t=max(fails)+1 if fails else 1
    print(f'p={p} k={k}: failures {fails} -> tau={t} (table {tau}) {"OK" if t==tau else "MISMATCH"}')
