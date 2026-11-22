from collections import Counter,defaultdict
from itertools import product
# left port L, left neighbor, x, y, right neighbor, right port R, apex A, apex B.
def rows(orientation=1, leftbit=0,rightbit=0):
 A=[0]*8
 for i in range(8):
  for j in range(i+1,8):
   if j<6: forward=True
   elif i==6:forward=orientation
   else:
    tags=((1,1),(leftbit,leftbit),(1,0),(0,1),(rightbit,rightbit),(0,0))
    forward=not tags[i][j-6]
   if forward:A[i]|=1<<j
   else:A[j]|=1<<i
 return tuple(A)

def sig(A,ports=(0,5)):
 n=len(A);U=((1<<n)-1)^sum(1<<j for j in ports)
 D={(0,(0,0),-1):1};R=Counter()
 for steps in range(n+4+1):
  E=defaultdict(int)
  for (mask,r,i),v in D.items():
   if mask==U and min(r)>=1:R[r]+=v
   for j in range(n):
    if i>=0 and (A[i]>>j&1)==0:continue
    if j in ports:
     t=ports.index(j);s=list(r);s[t]+=1
     if s[t]<=3:E[mask,tuple(s),j]+=v
    elif not (mask>>j&1):
     E[mask|1<<j,r,j]+=v
  D=E
 return R
from itertools import combinations
from random import Random

def flip(a,i,j):
    r=list(a);r[i]^=1<<j;r[j]^=1<<i
    return tuple(r)

def h(a):
    n=len(a)
    if not n:return 1
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for m in range(1,1<<n):
        for i,x in enumerate(d[m]):
            q=a[i]&~m
            while x and q:
                b=q&-q
                d[m|b][b.bit_length()-1]+=x
                q^=b
    return sum(d[-1])

def strong(a):
    full=(1<<len(a))-1
    for i in range(len(a)):
        s=1<<i
        while True:
            t=s
            for j,r in enumerate(a):
                if s>>j&1:t|=r
            if s==t:break
            s=t
        if s!=full:return False
    return True

def inflate(a,w):
    word=[i for i in range(len(a)) for _ in range(w.get(i,1))]
    return tuple(sum(1<<j for j,x in enumerate(word)
              if (x!=y and a[y]>>x&1) or (x==y and i<j))
              for i,y in enumerate(word))

def paired(a,u,v):
    assert a[u]>>v&1
    d=[j for j in range(len(a)) if j not in (u,v)]
    full=(1<<len(d))-1
    B=lambda M: tuple(sum(1<<k for k,j in enumerate(M) if a[i]>>j&1) for i in M)
    from itertools import permutations
    out=0
    for mask in range(1,full+1):
        S=[j for i,j in enumerate(d) if mask>>i&1]
        if len(S)%2==0:continue
        rest=[j for i,j in enumerate(d) if not mask>>i&1]
        mass=h(B(rest))
        for q in permutations(S):
            if not all(a[i]>>j&1 for i,j in zip(q,q[1:])):continue
            one=a[u]>>q[0]&1 and a[q[-1]]>>v&1
            two=a[v]>>q[0]&1 and a[q[-1]]>>u&1
            out+=mass*(int(one)-int(two))
    return 2*out

def allsmall():
    seen=set()
    for n in (3,4,5):
        pairs=list(combinations(range(n),2))
        for flags in range(1<<len(pairs)):
            a=[0]*n
            for k,(i,j) in enumerate(pairs):
                if flags>>k&1:a[i]|=1<<j
                else:a[j]|=1<<i
            if not(a[0]&2) or ((a[0]^a[1])>>2).bit_count()!=1:continue
            b=flip(a,0,1)
            if not strong(a) or not strong(b):continue
            x,y=h(a),h(b)
            if y-x==2:seen.add((x,y))
    assert seen=={(9,11),(13,15)}

def phase_rail():
    a=(1022,4092,1016,4080,3040,1984,3968,3840,512,3072,2325,293)
    b=flip(a,10,11)
    S,T=sig(a,(0,9)),sig(b,(0,9))
    assert tuple(tuple(S[i,j] for j in (1,2,3)) for i in (1,2,3))==((13919,20804,9250),(28900,45776,21626),(16090,26666,13122))
    def mid(X,Z):
        return 18998*X*X-30440*X*Z-10988*X+1982+7488*Z+13122*Z*Z
    for k in range(1,15):
        X,Z=2**k,3**k
        F={1:1,2:X-2,3:Z-3*X+3}
        M=mid(X,Z)
        assert M==1+sum(v*F[i]*F[j] for (i,j),v in S.items())
        assert M==-1+sum(v*F[i]*F[j] for (i,j),v in T.items())
        p,q=(inflate(d,{0:k,9:k}) for d in (a,b))
        assert strong(p) and strong(q)
        if k<=2:assert (h(p),h(q))==(M-1,M+1)
        if k%2==0:assert M%3==2
        elif k>=5:assert M%32==16
    assert [mid(2**k,3**k) for k in (1,3)]==[13920,4323024]
    for s in (5,6,7,17,24,33,49):
        mod=1<<s
        for R in (0,32,64,2**(s-1),mod-32):
            if R%32:continue
            u=0
            for i in range(s-5):
                if ((706+8433*u+26244*u*u-R//32)>>i)&1:u+=1<<i
            Zwant=1+8*u
            m=0
            for i in range(3,s):
                if ((pow(9,m,1<<(i+1))-Zwant)>>i)&1:m+=1<<(i-3)
            m+=(1<<(s-3))*s
            k=2*m
            assert mid(pow(2,k,mod),pow(3,k,mod))%mod==R%mod

def main():
    allsmall()
    phase_rail()
    r=Random(8643)
    for n in range(3,8):
        for _ in range(16):
            a=[0]*n
            for i,j in combinations(range(n),2):
                if (i,j)==(0,1) or r.randrange(2):a[i]|=1<<j
                else:a[j]|=1<<i
            assert paired(a,0,1)==h(flip(a,0,1))-h(a)
    lo=rows(1);hi=rows(0)
    assert lo==(62,252,184,112,224,192,133,9)
    assert hi==flip(lo,6,7)
    S=sig(lo)
    assert tuple(tuple(S[i,j] for j in (1,2,3)) for i in (1,2,3))==((259,236,66),(514,592,222),(258,362,162))
    for k in range(1,21):
        fills={1:1,2:2**k-2,3:3**k-3*2**k+3}
        pred=162*9**k-388*6**k+298*4**k+128*3**k-250*2**k+53
        assert pred==sum(c*fills[i]*fills[j] for (i,j),c in S.items())
        assert pred+2==sum(c*fills[i]*fills[j] for (i,j),c in sig(hi).items())
        x,y=(inflate(a,{0:k,5:k}) for a in (lo,hi))
        assert strong(x) and strong(y)
        if k<=4:assert h(x)==pred and h(y)==pred+2
    print("paired contraction, sharp nested-contact cases, nested +2 family PASS")
if __name__=='__main__':main()
