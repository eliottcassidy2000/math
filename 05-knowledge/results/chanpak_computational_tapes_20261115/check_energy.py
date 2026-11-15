from math import gcd
from itertools import combinations
from check_head import tournament,block,lift,path_dp

def path_ends(a):
    n=len(a); f=(1<<n)-1
    d=[[0]*n for _ in range(f+1)]
    for u in range(n):d[1<<u][u]=1
    for s,row in enumerate(d):
        for u,w in enumerate(row):
            if w:
                for v in range(n):
                    if a[u][v] and not s>>v&1:d[s|1<<v][v]+=w
    return d

def ring(n):
    return [[int(1<=(v-u)%n<=n//2) for v in range(n)] for u in range(n)]

def ring_energies(a):
    d=path_ends(a); e=path_ends(list(zip(*a))); n=len(a); f=(1<<n)-1
    w=[0]+[sum(d[s][0]*e[f^s][v]+d[s][v]*e[f^s][0] for s in range(f+1)) for v in range(1,n//2+1)]
    strides=[]; longest=0
    for k in range(2,n-1):
        values=sorted({sum(w[min(v-u,n-v+u)] for u,v in combinations(V,2)) for V in combinations(range(n),k)})
        stride=gcd(*(v-values[0] for v in values)); strides.append(stride)
        b=a_run=1
        for x,y in zip(values,values[1:]):
            a_run=a_run+1 if y-x==2 else 1
            b=max(a_run,b)
        longest=max(longest,b)
    return w[1:],strides,longest

if __name__=='__main__':
    D=tournament(5,76)
    assert tuple(block(D,{0,2}))==(15,120,42,8,42,6,6)
    assert tuple(block(D,{0,1}))==(15,120,42,10,42,6,6)
    for a,b in [(1,1),(1,4),(3,2),(2,2)]:
        tape=[tournament(a,(1<<(a*(a-1)//2))-1),D,tournament(b,(1<<(b*(b-1)//2))-1)]
        for V,J in [({0,1},10),({0,2},8)]:
            assert path_dp(lift(tape,[set(),V,set()]))==15+120*2**(a+b)+42*(2**a-2**b)-J
    for n in (5,7,9,11,13):
        print(n, ring_energies(ring(n)))
