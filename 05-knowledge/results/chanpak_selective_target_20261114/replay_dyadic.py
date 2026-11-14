"""Exact prescribed-target decoders with literal tournament controls."""
from functools import cache

def check(flag):
    if not flag:
        raise AssertionError('literal or inverse disagreement')

def p2(x):
    return x.bit_length()-1 if x>0 and x&(x-1)==0 else None

def divpow(N, D):
    if D==0 or N%D: return None
    return p2(N//D)

def scan(N, a, b, c, offsets, ties):
    for j,d in enumerate(offsets):
        M=N-d
        if M<=0:
            continue
        t=(M&-M).bit_length()-1
        # These families have v2(b)=0 and v2(c)=1.
        k=t
        if k>=1:
            m=divpow(M-b*2**k,a*2**k+c)
            if m is not None and m>=k:
                yield j,k,m
        m=t-1
        if m>=1:
            k=divpow(M-c*2**m,a*2**m+b)
            if k is not None and k>=m+2:
                yield j,k,m
        for m in ties(t):
            if m>=1 and M==a*2**(2*m+1)+b*2**(m+1)+c*2**m:
                yield j,m+1,m

def decode_pair(N):
    # Apply the shared decoder after dividing both sides by 3.
    for j,d in enumerate((71,73)):
        if (N-d)%3==0:
            for _,k,m in scan((N-d)//3,40,-19,-26,[0],
                              lambda t: ([1] if t==6 else [])+
                                        ([2] if t==10 else [])+
                                        ([t-6] if t>=9 else [])):
                yield j,k,m

def decode_triple(N):
    yield from scan(N,432,-303,-174,(225,227,229),lambda t:[t-2])

def transitive(n):
    return tuple((1<<n)-(2<<i) for i in range(n))

def ordinal(A,B):
    return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)

def apex(B,U):
    return tuple(x|((1<<len(B)) if not (U>>i)&1 else 0) for i,x in enumerate(B))+(U,)

def from_code(n,c):
    z=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if c&1:
                z[i]|=1<<j
            else:
                z[j]|=1<<i
            c>>=1
    return tuple(z)

def hp(B):
    n=len(B)
    @cache
    def starts(s,i):
        q=s^(1<<i)
        return sum(starts(q,j) for j in range(n) if (B[i]&q)>>j&1) if q else 1
    full=(1<<n)-1
    return sum(starts(full,i) for i in range(n))

R5=tuple((1<<((i+3)%5))|(1<<((i+4)%5)) for i in range(5))
cases=[
 (120,-57,-78,(71,73),[(R5,17),(R5,9)],decode_pair),
 (432,-303,-174,(225,227,229),[(from_code(6,code),mask)
     for code,mask in ((1332,23),(408,15),(408,27))],decode_triple),
]
for a,b,c,offsets,cores,decode in cases:
    cap=40000
    table={}
    for j,d in enumerate(offsets):
        for k in range(1,50):
            for m in range(1,50):
                v=a*2**(k+m)+b*2**k+c*2**m+d
                check((j,k,m) in set(decode(v)))
                if v<=cap:
                    table.setdefault(v,set()).add((j,k,m))
                if max(k,m)<=4:
                    R,U=cores[j]
                    S=ordinal(ordinal(transitive(k),R),transitive(m))
                    mask=((1<<k)-1)|(U<<k)
                    check(hp(apex(S,mask))==v)
    for N in range(1,cap+1,2):
        check(set(decode(N))==table.get(N,set()))

check(set(decode_pair(13897))=={(1,4,3)})
print('80 literal two-height controls; exact inverse below 40001 and at forward heights through 49; pair/triple witnesses pass.')
