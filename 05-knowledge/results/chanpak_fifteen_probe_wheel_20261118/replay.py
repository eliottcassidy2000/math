from itertools import product, combinations, permutations
from collections import deque, Counter

HOSTS={'P':((10,12,9,0),0),'Q':((2,4,9,3),0),'R':((2,4,9,3),1)}
def substitute(A,letter):
    K,u=HOSTS[letter];m=len(A)
    spec=[(v,j) for v in range(4) for j in (range(m) if v==u else [-1])]
    return [sum(1<<l for l,(w,k) in enumerate(spec)
                if (v==w==u and A[j]>>k&1) or (v!=w and K[v]>>w&1)) for v,j in spec]

def reads(A,J):
    n=len(A);N=1<<n
    E=[[0]*n for _ in range(N)]
    for i in range(n):E[1<<i][i]=1
    for s in range(1,N):
        for i,b in enumerate(E[s]):
            if not b:continue
            t=A[i]&~s
            while t:
                c=t&-t;t^=c;E[s|c][c.bit_length()-1]+=b
    H=[sum(e) for e in E];H[0]=1
    ans=[0]
    F=[0]*N;F[0]=1
    for j in range(1,J+1):
        C=[0]*N
        for s in range(N):
            a=s
            while True:
                C[s]+=H[a]*F[s^a]
                if not a:break
                a=(a-1)&s
        F=C;ans.append(F[-1])
    return ans
def act(F,c,t):
    f,plus,minus=F[t],F[t+1],F[t-1]
    z=t**3*f
    if c=='P':return z+t*t*(plus-minus)
    if c=='Q':return z+t*((t+1)*plus+(1-t)*minus)
    return z+2*t*t*(plus-minus)+t*(plus+minus-2*f)
def all_rows(n):
    E=list(combinations(range(n),2))
    for b in range(1<<len(E)):
        R=[0]*n
        for k,(i,j) in enumerate(E):
            a,c=(j,i) if b>>k&1 else (i,j)
            R[a]|=1<<c
        yield R
checks=0
for n in range(1,5):
    for A in all_rows(n):
        F=reads(A,6)
        for c in HOSTS:
            G=reads(substitute(A,c),5)
            assert G[1:]==[act(F,c,t) for t in range(1,6)]
            checks+=1
print("literal five-read checks",checks)

# 0 indicates even order. Only transitions that preserve 3-units are present.
M={(0,'P'):((1,1),(1,2)),(1,'Q'):((1,2),(3,2)),
   (0,'R'):((4,3),(4,4)),(1,'R'):((4,3),(4,4))}
for (e,c),F in M.items():
    for i,(u,v) in enumerate(((1,0),(0,1))):
        z=(-1)**e
        A=[0,u,v,z*v,z*u,0]
        assert tuple(act(A,c,t)%5 for t in (1,2))==tuple(row[i] for row in F)
states=list(product(range(2),range(1,5),range(5)))
def step(s,c):
    p,u,v=s
    return (1-p,*[sum(a*b for a,b in zip(r,(u,v)))%5 for r in M[p,c]])
def neighbors(s):
    return [(c,out) for c in ('P' if s[0]==0 else 'Q','R') if (out:=step(s,c))[1]]
live={s for s in states if neighbors(s)}
dead=set(states)-live
assert dead=={(1,u,2*u%5) for u in range(1,5)}
def bfs(s):
    D={s:''};Q=deque([s])
    while Q:
        z=Q.popleft()
        for c,w in neighbors(z):
            if w in live and w not in D:
                D[w]=D[z]+c;Q.append(w)
    return D
diam=0
for s in live:
    F=bfs(s)
    assert set(F)==live
    diam=max(diam,*map(len,F.values()))
print('live, traps, diameter:',len(live),len(dead),diam)
A=(2,0);seed=(0,1,4);B=bfs(seed)
for s in sorted(set((p,u) for p,u,v in live)):
    path=min((w for (p,u,v),w in B.items() if (p,u)==s),key=lambda w:(len(w),w))
    # Exact integer reads using only the necessary initial values.
    arr=reads(A,len(path)+2)
    for c in path:
        arr=[0]+[act(arr,c,t) for t in range(1,len(arr)-1)]
    h=arr[1]
    assert h%5==s[1] and h%3==(-1)**s[0]%3 and h%2
    print('residue',h%30,'word',path or '(empty)','H',h)

# The 15 path signs for the cyclic five-tournament, in cyclic label order.
paths=[p for p in permutations(range(5)) if all((v-u)%5 in (1,2) for u,v in zip(p,p[1:]))]
sgn=lambda p:(-1)**sum(p[i]>p[j] for i,j in combinations(range(5),2))
signs=Counter(map(sgn,paths))
assert signs=={-1:10,1:5}
print('cyclic five-tournament path permutation signs:',sorted(signs.items()))
