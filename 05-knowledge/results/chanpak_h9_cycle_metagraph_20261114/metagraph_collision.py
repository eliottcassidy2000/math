"""Compare rooted line-neighbor spectra of two rigid typed-spectrum twins."""
from itertools import permutations
from collections import Counter
A=(34,68,33,87,71,90,5)
B=(70,64,98,39,47,3,56)
P=list(permutations(range(7)))

def check(ok):
    if not ok:raise AssertionError('collision line control')

def ham(T):
    return [p for p in P if all(T[i]>>j&1 for i,j in zip(p,p[1:]))]

def dual(T,p):
    # Nonpath complement; flipping only path arcs gives its converse.
    U=list(T)
    for i,j in zip(p,p[1:]):
        U[i]^=1<<j;U[j]^=1<<i
    return tuple(U)

def neighbor(T):
    ans=Counter()
    for p in ham(T):
        X=dual(T,p)
        hp=ham(X)
        nC=sum(p[0]==0 and X[p[-1]]&1!=0 for p in hp)
        ans[bool(T[p[-1]]>>p[0]&1),len(hp),nC]+=1
    return ans

def direct_spectra(T):
    ways=[]
    cycles=[]
    for n in range(1,8):
        orders=[p for p in permutations(range(7),n)
                if all(T[i]>>j&1 for i,j in zip(p,p[1:]))]
        ways.append(len(orders))
        if n>=3:
            cycles.append(sum(p[0]==min(p) and bool(T[p[-1]]>>p[0]&1)
                              for p in orders))
    return tuple(ways),tuple(cycles)
def group_count(T, anti):
    return sum(all((((T[i]>>j)&1)^((T[p[i]]>>p[j])&1))==anti
                   for i in range(7) for j in range(i+1,7)) for p in P)
for T, expected in [(A,[33,37,37]),(B,[23,25,33,33,37,37])]:
    check(direct_spectra(T)==((7,21,53,107,161,151,65),(9,15,16,11,3)))
    check((group_count(T,0),group_count(T,1))==(1,0))
    N=neighbor(T)
    check(sum(N.values())==65)
    check(sum(q for (closed,h,c),q in N.items() if closed)==21)
    check(sorted(h for (closed,h,c),q in N.items() if closed and c==0
                 for _ in range(q))==expected)
a,b=sorted(x.bit_count() for x in A),sorted(x.bit_count() for x in B)
check(a!=b and sorted(6-x for x in a)!=b)
print('Equal all-length path/cycle totals, distinct rigid merged nodes, and different nonstrong line-neighborhoods.')
