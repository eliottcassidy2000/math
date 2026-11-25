from collections import Counter
from itertools import combinations, permutations, product
from math import factorial

def word_h(w):
    if not w: return 1
    ans=w[0]+1-w[-1]
    b=0
    for i,x in enumerate(w):
        if not x: ans+=b
        if i and not w[i-1] and x: ans+=1
        b=2*b+x
    return ans

def value(w):
    b=0
    for x in w:b=2*b+x
    return b

def mate_h(w):
    assert w and w[-1]==0
    empty=(0,0,0,1) # clipped run count, top run singleton, pending bad descent, endpoint permitted
    D={(empty,empty,2):1}
    for i,a in enumerate(w):
        out=1-a if i+1<len(w) else 0 # flipped apex outcontact a -> i
        N=Counter()
        for (l,r,last), num in D.items():
            for j,old in enumerate((l,r)):
                c,s,b,ok=old
                # For left lane keep the low run's last vertex, for right the high run's first.
                good=(1-out) if j==0 else out
                z=[]
                if not c:
                    z.append((1,1,0,good))
                else:
                    if j==last:
                        z.append((c,0,0,good if j==0 and c==1 else ok))
                    if not b:
                        z.append((2,1,int(j==last and s),good if j==1 else ok))
                for new in z:
                    N[(new,r,j) if j==0 else (l,new,j)]+=num
        D=N
    return sum(n for (l,r,j),n in D.items() if not l[2] and not r[2] and l[3] and r[3])

def tournament(n, rev=0):
    A=[[0]*n for i in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[j][i]=(rev>>k)&1
        A[i][j]=1-A[j][i]
    return A

def hp(A):
    n=len(A); M=1<<n
    dp=[[0]*n for _ in range(M)]
    for s in range(1,M):
        for i in range(n):
            if s>>i&1:
                t=s^(1<<i)
                dp[s][i]=sum(dp[t][j]*A[j][i] for j in range(n)) if t else 1
    return sum(dp[-1])

def direct(w,Y=(),contacts=(),K=0,flip=False):
    ny=len(Y); L=len(w); n=L+ny+K+(1 if K else 0)+1
    A=tournament(n)
    if K: word=tuple(w)+tuple(contacts)+(1,)*K+(0,)
    else: word=tuple(w)
    for i in range(n-1):A[-1][i]=word[i];A[i][-1]=1-word[i]
    for i in range(ny):
        for j in range(ny):
            A[L+i][L+j]=Y[i][j]
    if flip:
        for i,j in combinations(range(n),2):
            if j>i+1:A[i][j],A[j][i]=A[j][i],A[i][j]
    return A

def build_graph(d):
    if d in (2,3):
        G=[set() for _ in range(1+d*d)]
        def add(a,b):G[a].add(b);G[b].add(a)
        for i in range(5):
            add(i,(i+1)%5)
            if d==3:add(i,i+5);add(i+5,5+(i+2)%5)
    else:
        assert d==7
        G=[set() for _ in range(50)]
        def add(a,b):G[a].add(b);G[b].add(a)
        for h in range(5):
            for x in range(5):
                add(5*h+x,5*h+(x+1)%5)
                add(25+5*h+x,25+5*h+(x+2)%5)
                for i in range(5):add(5*h+x,25+5*i+(x+h*i)%5)
    return G


from math import factorial
from math import isqrt
def bits(a,l):return tuple(map(int,format(a,f'0{l}b')))
def val(w):
 b=0
 for x in w:b=2*b+x
 return b
def word(W,k):return (1,0,1,0)+bits(k,W)+(0,)
def transform(v,L):return tuple(1-x for x in bits(v,L)[::-1])
def half(w):
 L=len(w)//2
 return val(w[:L]),val(1-x for x in w[L+1:][::-1])
def K(x):return (x+1)//2
def lower(u,v):
 return 2*u*v+word_h(bits(u,u.bit_length()))+word_h(bits(v,v.bit_length()))-1+2*min(K(u),K(v))
def inverse(H):
 s=isqrt(H//2)
 for a in range(max(1,s-2),s+1):
  for b in range(a+1,a+4):
   if abs(K(a)-K(b))!=1 or lower(a,b)!=H:continue
   for u,v in [(a,b),(b,a)]:
    L=u.bit_length()
    if L>=5 and v.bit_length()==L and u>>(L-4)==10:
     w=bits(u,L)+(0,)+transform(v,L)
     yield w

for W in range(3,13):
 H=[word_h(word(W,k)) for k in range(1<<W)]
 S=[]
 for k in range(1<<W):
  for i in range(W):
   l=k^(1<<i)
   if k<l and abs(H[k]-H[l])==2:
    S.append((k,l))
    if W>=5:
     assert i==W//2+1
     u,v=half(word(W,k))
     assert lower(u,v)==min(H[k],H[l])
     assert word(W,k) in set(inverse(min(H[k],H[l])))
 assert len(S)==({3:1,4:3}.get(W,0 if W%2 else 2**(W//2)))
 print('width',W,'seams',len(S))
 if W==8:
  sources=[(k,l) for k,l in S if l<144]
  mates=[]
  for k in range(144):
   for i in range(8):
    l=k^(1<<i)
    if k<l<144 and abs(mate_h(word(W,k))-mate_h(word(W,l)))==2:mates.append((k,l))
  assert len(sources)==8 and mates==[(6,38)]
  for k,l in sources+mates:
   got=[]
   for z in (k,l):
    w=word(W,z);hh=hp(direct(w));cc=hp(direct(w,flip=True))
    assert hh==H[z] and cc==mate_h(w)
    got.append((hh,cc))
   print(k,l,got,half(word(W,k)))
