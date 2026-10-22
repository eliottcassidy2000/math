# Independent check (pure Python, forward EGF-style recurrence written differently from the C code):
# v_p(D(p^a)) computed modulo p^(a+1); freezing lemma prediction: once v_p(D(p^b)) < b+1, it stays constant.
# D(n) = s(n) - c(n), s = strong labelled tournaments, c(n) = 2^((n-1)/2) for odd n.
def s_upto(N, q):
    # s(n) = T(n) - sum_{k=1}^{n-1} C(n,k) s(k) T(n-k), T(j) = 2^C(j,2) mod q; binomials via multiplicative formula mod q is
    # unsafe (division), so use Pascal rows kept as a list (O(N^2) memory-light: one row at a time).
    T=[pow(2,j*(j-1)//2,q) for j in range(N+1)]
    s=[0]*(N+1); row=[1]
    for n in range(1,N+1):
        row=[1]+[(row[i-1]+row[i])%q for i in range(1,n)]+[1]
        tot=0
        for k in range(1,n):
            tot+=row[k]*s[k]%q*T[n-k]
        s[n]=(T[n]-tot)%q
    return s
def v(x,p,cap):
    if x==0: return cap
    c=0
    while x%p==0: x//=p; c+=1
    return c
cases=[(3,6),(5,4),(7,3),(11,2),(13,3)]
for p,A in cases:
    q=p**(A+1); N=p**A; s=s_upto(N,q)
    vals=[]
    for a in range(1,A+1):
        n=p**a; D=(s[n]-pow(2,(n-1)//2,q))%q; vals.append(v(D,p,A+1))
    print(f'p={p}: v_p(D(p^a)) for a=1..{A} (computed mod p^{A+1}; value {A+1} means >= {A+1}): {vals}')
