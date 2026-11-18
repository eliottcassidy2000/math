"""All-order restricted inverse move tests from strong component compression."""
from itertools import combinations
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moonshine_prime_blocks'))
from replay import dp,profile

def H(T):return sum(dp(T)[-1]) if T else 1
def canon_lane(n):
    s=1 if n%4==1 else -1
    z=3*n+s
    return z//(z&-z)
def restrict(T,Q):
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(Q)) for i in Q)
def reps(q):
    E=set()
    for line in q.open():
        T=tuple(map(int,line.split()))
        E.add(T)
        for n in range(2,len(T)):
            for J in combinations(range(len(T)),n): E.add(restrict(T,J))
    return E
def rep8(path):
    import csv
    with path.open() as f:
        for r in csv.DictReader(f):
            z=int(r['code'])
            yield tuple((z>>(8*i))&255 for i in range(8))
def flip(A,i,j):
    U=list(A);U[i]^=1<<j;U[j]^=1<<i
    return tuple(U)

def run():
    E=reps(Path(sys.argv[1]))
    cap=43
    if len(sys.argv)>2:
        E.update(rep8(Path(sys.argv[2])));cap=73
    all_edges={}
    target={tuple(sorted((canon_lane(N),N))) for N in range(3,cap+2,2) if N not in (7,21) and canon_lane(N) not in (7,21)}
    other={tuple(sorted((N-2,N))) for N in range(3,cap+2,2) if N not in (7,21) and N-2 not in (7,21)}
    need=target|other
    start={h for a,b in need for h in (a,b)}
    for A in E:
        v=H(A)
        if v not in start:continue
        for i,j in combinations(range(len(A)),2):
            B=flip(A,i,j);w=H(B)
            if w not in start:continue
            X=tuple(sorted((v,w)))
            if X not in need:continue
            z=(len(A),A if v<w else B,(i,j),B if v<w else A)
            if X not in all_edges or z<all_edges[X]:all_edges[X]=z
    for v in sorted(target,key=lambda r:r[1]):
        print(v,all_edges.get(v))
    missing=target-all_edges.keys()
    print('same-order one-arc failures',sorted(missing),'relevant edges',len(all_edges))
    if cap==73:
        assert missing == {(1,11),(1,43),(5,27),(13,17),(13,35),(17,45),(23,61),
                           (25,67),(29,39),(35,47),(47,63)}
    print('unit-step failures' ,sorted(other-all_edges.keys()))
if __name__=='__main__':run()
