from itertools import product

def tt(n):
    return tuple(((1<<n)-1) ^ ((2<<i)-1) for i in range(n))

def ordinal(a,b):
    mask=((1<<len(b))-1)<<len(a)
    return tuple(u|mask for u in a)+tuple(v<<len(a) for v in b)

def covers(a):
    n=len(a); F=(1<<n)-1
    dp=[[0]*n for _ in range(F+1)]
    hs=[1]+[0]*F
    for i in range(n):
        dp[1<<i][i]=1
    for m in range(1,F+1):
        for i,w in enumerate(dp[m]):
            if not w: continue
            hs[m]+=w
            nxt=a[i]&~m
            while nxt:
                l=nxt&-nxt; nxt-=l
                dp[m|l][l.bit_length()-1]+=w
    return hs[F],sum(hs[m]*hs[F^m] for m in range(F+1))

def wrap(b):
    n=len(b)
    return (2, ((1<<(n+2))-1)^3)+tuple((x<<2)|1 for x in b)

def vertices(y):
    width=max(1,y.bit_length())
    return tuple((y>>(width-i-1))&1 for i in range(width))

def reader(y):
    bits=vertices(y)
    n=len(bits); a=[0]*(n+1)
    for i,b in enumerate(bits,1):
        for j in range(i+1,n+1):
            a[j]|=1<<i
        if b: a[i]|=1
        else: a[0]|=1<<i
    return tuple(a)

def seam(src,rdr):
    a=list(ordinal(src,rdr)); s=len(src)
    a[0]^=1<<s; a[s]^=1
    return tuple(a)

def g(x):
    a=1
    while x:
        u=x//2
        if not x%2: a+=u+u%2
        x=u
    return a

def independent_reader(y):
    # Counts by actual apex-cycle leg choices, not by the recurrence.
    bits=vertices(y)
    h=1; e=0; o=1
    for i,b in enumerate(bits):
        if not b: continue
        later=len(bits)-1-i
        if not later: e+=1
        else: e+=1<<(later-1); o+=1<<(later-1)
        for j in range(i+1,len(bits)):
            if not bits[j]:
                h+=1<<max(1,j-i-1)
    return h,e,o

def phase(j,v):
    beta=int(f'{v:0{j}b}'[::-1],2)^((1<<j)-1)
    M=(1<<j)+beta
    return M,g(v)+v-M+1

def guarded_inverse(N,h,Q,j,v):
    W=h+Q
    if (N+h)%W: return []
    M,J=phase(j,v)
    L=(N+h)//W-J
    out=[]
    for k in range(1,max(1,L.bit_length())):
        A=M<<k
        if A>=L: break
        u=(L-1)//A
        if 1<=u<=A and L-A*u==g(u-1):
            out.append((u,k))
    return out

BASES=((tt(2),1,0,3,45,9), (tt(1),1,0,5,45,11),
    (tt(3),4,6,8505,1913625,67913),
    ((36,33,26,99,75,84,7),3,0,75,1913625,67915),
    (ordinal((24,73,11,32,78,87,13),tt(1)),3,0,45,1913625,67917))

def run():
    literals=0
    for b,y in product((tt(1),tt(2),tt(3),(2,4,1)),range(32)):
        h,Q=covers(b)
        eta,e,o=independent_reader(y)
        assert eta==g(y) and (e,o)==((y+1)//2,1+y//2)
        assert covers(seam(wrap(b),reader(y)))[0]==(h+Q)*(eta+e+o)-h*(o-e)
        literals+=1
    for b,j,v,u,K,O in BASES:
        h,Q=covers(b); W=h+Q; M,J=phase(j,v)
        assert (W*M*u,W*(g(u-1)+J)-h)==(K,O)
        for k in range(1,30):
            y=((u<<k)-1)*(1<<j)+v
            eta,e,o=independent_reader(y)
            N=(K<<k)+O
            assert W*(eta+e+o)-h*(o-e)==N
            if u<=M<<k: assert (u,k) in guarded_inverse(N,h,Q,j,v)
        if K==45:
            y=((u<<1)-1)*(1<<j)+v
            assert covers(seam(wrap(b),reader(y)))[0]==2*K+O
    x=(50,92,105,17,68,26,41); z=(12,73,34,84,7,27,53)
    assert (covers(x),covers(z))==((139,1664),(143,1660))
    for k in (1,2,3):
        for b in (x,z):
            h,Q=covers(b)
            assert covers(seam(wrap(b),reader(k)))[0]==1803*(g(k)+k+1)-h*(1-k%2)
    print(f'PASS: {literals} literal seam comparisons, 2 literal first endpoints, '
          'separate apex-leg recount at heights 1..29 for every rail, '
          'guarded inverses, and exact gap-4 sources.')
if __name__=='__main__': run()
