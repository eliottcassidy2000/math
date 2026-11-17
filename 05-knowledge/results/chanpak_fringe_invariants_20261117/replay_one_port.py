from itertools import product, permutations, combinations
from collections import Counter
from functools import lru_cache
from random import Random
from math import factorial

def require(value):
    if not value: raise ArithmeticError('control failed')

def literal_h(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i,w in enumerate(dp[s]):
            if not w:continue
            b=rows[i]&~s
            while b:
                bit=b&-b;b-=bit;dp[s|bit][bit.bit_length()-1]+=w
    return sum(dp[-1])

def profile(D,out):
    n=len(D);ans=Counter()
    for perm in permutations(range(n)):
        for r in range(1,n+1):
            for cuts in combinations(range(1,n),r-1):
                if any(j not in cuts and not (D[perm[j-1]]>>perm[j]&1) for j in range(1,n)):continue
                ans['R',r]+=1
                e=not ((out>>perm[(cuts+(n,))[0]-1])&1)
                s=(out>>perm[0])&1
                if e:ans['E',r]+=1
                if s:ans['S',r]+=1
                if r>=2 and e and out>>perm[cuts[0]]&1:ans['B',r]+=1
    return ans

def masks(D,out,a,b,x,y,one_side):
    n=len(D);v=n; A=list(range(n+1,n+1+a));F=list(range(n+1+a,n+1+a+b))
    size=n+1+a+b; rows=list(D)+[out]+[0]*(a+b)
    def edge(i,j):rows[i]|=1<<j
    for i in range(n):
        if not (out>>i&1):edge(i,v)
        for j in A:edge(i,j)
        for j in F:
            if one_side:edge(i,j)
            else:edge(j,i)
    for i,j in combinations(A+F,2):edge(i,j)
    for j in A:edge(v,j) if (x if one_side else not y) else edge(j,v)
    for j in F:edge(v,j) if (y if one_side else x) else edge(j,v)
    return tuple(rows)

@lru_cache(None)
def coeff(n,a,b,x,y,one_side=False):
    base=masks((0,),0,a,b,x,y,one_side);total=len(base)
    def visit(part,left,r):
        if not left and not r:
            if part.count(0)==0:return
            v=part.index(1)
            e=v>0 and part[v-1]==0; s=v+1<len(part) and part[v+1]==0
            typ='B' if e and s else 'E' if e else 'S' if s else 'R'
            ans[typ,part.count(0)]+=1;return
        for j in ([0] if r else [])+list(left):
            if part and not(0 in (part[-1],j) and 1 in (part[-1],j)):
                if not (base[part[-1]]>>j&1):continue
            if j==0:visit(part+(j,),left,r-1)
            else:visit(part+(j,),left-{j},r)
    ans=Counter()
    for r in range(1,n+1):visit((),set(range(1,total)),r)
    return ans

def signed(n,*probes):
    ans=Counter()
    for sign,args in probes:
        for k,v in coeff(n,*args).items():ans[k]+=sign*v
    return {k:v for k,v in ans.items() if v}

def diff(n,a,b,x,y,side,axis):
    if axis==0:return signed(n,(1,(a,b,1,y,side)),(-1,(a,b,0,y,side)))
    return signed(n,(1,(a,b,x,1,side)),(-1,(a,b,x,0,side)))
def mix(n,a,b,side):
    return signed(n,*[((1 if x==y else -1),(a,b,x,y,side)) for x,y in product(range(2),repeat=2)])

n=4
require(mix(n,1,1,True)=={('R',1):-2})
for r in range(2,n+1):
    z=mix(n,r-1,r-1,False)
    require(all(k=='R' and j<=r for k,j in z) and z['R',r]==(r-1)*factorial(r-1)**2)
    for axis,a,b,t in [(0,r-2,r-1,'E'),(1,r-1,r-2,'S')]:
        z=diff(n,a,b,0,0,False,axis)
        require(set(k for k in z if k[1]>=r)=={(t,r)} and all(k[0]!='B' for k in z))
        require(z[t,r]==(r-1)*factorial(r-2)*factorial(r-1))
    z=coeff(n,r-2,r-2,0,0,False)
    require(set(k for k in z if k[1]>=r)=={('B',r)} and z['B',r]==(r-1)*factorial(r-2)**2)
z=diff(n,1,1,0,0,True,1)
require(z['S',1]==-2 and all(k[0] in 'RS' for k in z))
z=diff(n,1,1,0,0,True,0)
require(z['E',1]==1 and z['S',1]==-1 and all(k[0] in 'RES' for k in z))

rng=Random(250731)
tests=0
for n in range(1,5):
    for trial in range(7):
        D=[0]*n
        for i,j in combinations(range(n),2):
            if rng.randrange(2):D[i]|=1<<j
            else:D[j]|=1<<i
        out=rng.randrange(1<<n)
        data=profile(D,out)
        for a,b,one in [(0,0,False),(1,1,True),(0,1,False),(1,0,False),(1,1,False),(2,1,False),(1,2,False),(2,2,False)]:
            for x,y in product(range(2),repeat=2):
                val=sum(data[k]*v for k,v in coeff(n,a,b,x,y,one).items())
                require(val==literal_h(masks(D,out,a,b,x,y,one)))
                tests+=1
print(f'PASS: triangular coefficients through r=4 and {tests} literal DP profiles')

D=(14,28,8,16,5)
U,V=[masks(D,z,0,0,0,0,False) for z in (10,5)]
common={k:x for name,vals in [('R',(9,66,186,240,120)),('E',(5,38,110,144,72)),('S',(4,28,76,96,48)),('B',(0,20,56,72,36))] for k,x in [((name,j+1),v) for j,v in enumerate(vals)] if x}
require(profile(D,10)==profile(D,5)==common)
require(literal_h(U)==literal_h(V)==29)
def expose0(rows):
    return tuple(r+((i!=0)<<6) for i,r in enumerate(rows))+(1,)
require((literal_h(expose0(U)),literal_h(expose0(V)))==(115,75))
print('PASS: active six-vertex collision, split by second-port read')
def add(rows,bits):
    return tuple(r|((1-(bits>>i&1))<<len(rows)) for i,r in enumerate(rows))+(bits,)
require((literal_h(add(U,51)),literal_h(add(V,51)))==(135,133))
print('PASS: second read exposing two additional D exceptions gives 135/133')
