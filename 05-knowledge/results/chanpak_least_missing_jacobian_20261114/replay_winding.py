from functools import cache
from itertools import permutations

def check(ok):
    if not ok:
        raise AssertionError('physical relative-winding check')

def profile(T):
    n=len(T); M=(1<<n)-1
    @cache
    def A(s,i):
        t=s^(1<<i)
        if not t:return tuple(int(i==j) for j in range(n))
        ls=[A(t,j) for j in range(n) if ((t&T[i])>>j)&1]
        return tuple(sum(a[j] for a in ls) for j in range(n))
    @cache
    def H(s):
        return sum(sum(A(s,i)) for i in range(n) if (s>>i)&1) if s else 1
    data=[[[0]*n for _ in range(n)] for _ in range(2)]
    Q=0
    for s in range(M+1):
        h=H(M^s);Q+=h*H(s)
        for i in range(n):
            if (s>>i)&1:
                for j,v in enumerate(A(s,i)):
                    data[s.bit_count()%2][i][j]+=h*v
    return H(M),Q,*data

def apex(T,U):
    n=len(T)
    return tuple(x|((1<<n) if not (U>>i)&1 else 0) for i,x in enumerate(T))+(U,)

def counts(T,U):
    n=len(T)
    P=apex(T,U);N=n+1
    @cache
    def literal(s,i):
        m=s^(1<<i)
        return sum(literal(m,j) for j in range(N) if (m&P[i])>>j&1) if m else 1
    return sum(literal((1<<N)-1,i) for i in range(N))

plus=(50,44,25,17,34,12)
minus=(50,44,25,33,10,20)
for T in (plus,minus):
    check(sum(all(((T[i]>>j)&1)==((T[p[i]]>>p[j])&1)
                  for i in range(6) for j in range(i+1,6))
              for p in permutations(range(6)))==3)
    h,Q,K,O=profile(T)
    check((h,Q)==(45,432))
    for U in range(64):
        H=counts(T,U)
        check(H==h+2*sum(K[i][j] for i in range(6) for j in range(6)
                         if (U>>i)&1 and not (U>>j)&1))
        if T==minus:
            check(H%3==0)
    if T==minus:
        for mat in (K,O):
            check(all(sum(row)%3==0 for row in mat))
            check(all(sum(row[i] for row in mat)%3==0 for i in range(6)))
        check(all((K[i][j]+K[j][i])%3==0 for i in range(6) for j in range(6)))
        check([sum(row) for row in O]==[39]*3+[33]*3)
        check([sum(row[i] for row in O) for i in range(6)]==[33]*3+[39]*3)
check((counts(plus,1),counts(plus,8))==(109,95))
print('128 literal apex counts; common Hamilton count, cover count and symmetry order; endpoint lock controls passed.')
