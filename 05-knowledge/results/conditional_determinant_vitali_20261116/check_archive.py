from functools import cache
from itertools import product
A = (2,12,17,37,75,23,47)
B = (2,12,17,37,43,71,31)
assert tuple(a.bit_count() for a in A)==tuple(b.bit_count() for b in B)
perm=[0,1,3,2,5,4,6]
Bprime=tuple(sum(((B[perm[i]]>>perm[j])&1)<<j for j in range(7)) for i in range(7))
assert {(i,j) for i in range(7) for j in range(i+1,7) if (Bprime[i]^A[i])>>j&1}=={(2,3)}

def lift(T,n,m,r):
    blocks=[i for i,s in enumerate((n,m,1,1,1,1,r)) for _ in range(s)]
    return tuple(sum((int(i<j) if a==b else (T[a]>>b)&1)<<j for j,b in enumerate(blocks)) for i,a in enumerate(blocks))

def H(T):
    @cache
    def f(M,i):
        if M==1<<i:return 1
        S=M^(1<<i)
        return sum(f(S,j) for j in range(len(T)) if (S&T[i])>>j&1)
    return sum(f((1<<len(T))-1,j) for j in range(len(T)))

def recover(d):
    h=1-d
    n=(h&-h).bit_length()-1
    Q=h>>n
    if d>=0 or n<1:return None
    if (Q+1)&Q==0:
        return n,1,(Q+1).bit_length()-1
    L=(Q-1).bit_length()
    v=(1<<L)+1-Q
    s=[i for i in range(v.bit_length()) if v>>i&1]
    if len(s)==1:s=[s[0]-1]*2
    if len(s)==2 and min(s)>=2 and sum(s)==L:return (n,*s)

assert (H(A),H(B))==(67,69)
assert (H(lift(A,1,1,2)),H(lift(B,1,1,2)))==(135,133)
for n,m,r in product(range(1,4),repeat=3):
    d=1-(1<<n)*((1<<m)-1)*((1<<(r-1))-1)
    assert H(lift(B,n,m,r))-H(lift(A,n,m,r))==2*d
for n,x,y in product(range(1,32),repeat=3):
    d=1-(1<<n)*((1<<x)-1)*((1<<y)-1)
    assert recover(d)==(n,*sorted((x,y)))
print('Two literal count fibres verified, 27 inflated signs, and 29,791 Mersenne inverse checks.')
