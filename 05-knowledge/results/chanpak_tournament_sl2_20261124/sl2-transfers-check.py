from itertools import permutations
from collections import Counter

def paths(g,start=None,end=None):
    size=len(g)
    def run(v,seen,word,weight):
        if len(word)==size:
            if end is None or end==v:
                yield word,weight
            return
        for j,a in enumerate(g[v]):
            if a and not (seen>>j&1):
                yield from run(j,seen|1<<j,word+(j,),weight*a)
    for i in range(size) if start is None else [start]:
        yield from run(i,1<<i,(i,),1)

def mul(A,B,m=None):
    a,b,c,d=A; e,f,g,h=B
    v=(a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
    return v if m is None else tuple(x%m for x in v)

U=(1,1,0,1); L=(1,0,1,1); I=(1,0,0,1)
Ms=[(1,3,-1,2),L,(4,2,3,1)]
exits=[2,0]
g=[[0]*12 for _ in range(12)]
for k in range(4):
    for i in range(3):
        g[3*k+i][3*k+(i+1)%3]=1
    if k<3:
        for i in range(2):
            for j in range(2):
                g[3*k+exits[i]][3*(k+1)+j]=Ms[k][2*i+j]
assert tuple(sum(w for _,w in paths(g,i,9+exits[j]))
             for i in range(2) for j in range(2))==mul(mul(Ms[0],Ms[1]),Ms[2])

for i,j,M in [(0,3,U),(2,4,L)]:
    g=[[int((a//3==b//3 and (b-a)%3==1) or a//3<b//3)
        for b in range(6)] for a in range(6)]
    g[i][j],g[j][i]=0,1
    pp=[w for w,z in paths(g)]
    assert len(pp)==17
    framed=[w for w in pp if max(w[:3])<3]
    assert len(framed)==8
    assert tuple(sum(w[0]==a and w[3]==b+3 for w in framed)
                 for a in range(2) for b in range(2))==M

for m,order in [(2,6),(3,24),(5,120),(7,336),(8,384),(12,1152)]:
    Q=[I]; seen=set(Q)
    for a in Q:
        for b in [U,L]:
            p=mul(a,b,m)
            if p not in seen:
                Q.append(p); seen.add(p)
    assert len(Q)==order

for m in [3,4,5]:
    g=[[0]*(2*m) for _ in range(2*m)]
    for i in range(1,m):
        g[i][m+i]=1
    for i in range(m):
        for j in range(m):
            if i!=j: g[m+i][j]=1
    got=Counter(tuple(sorted((a-m,b) for a,b in zip(w,w[1:]) if a>=m))
                for w,z in paths(g))
    target=Counter()
    for pp in permutations(range(1,m)):
        w=(0,)+pp+(0,)
        target[tuple(sorted(zip(w,w[1:])))]+=1
    assert got==target
print("checked weighted integer transfers, unweighted leakage, six congruence orders, and exact polynomial supports")
