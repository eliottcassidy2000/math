from collections import Counter
from itertools import permutations, combinations, product

Q={1,2,4}
arc=lambda a,b:(b-a)%7 in Q
fano={}
for a in range(7):
    b,c=(a+1)%7,(a+3)%7
    for x,y,z in ((a,b,c),(b,c,a),(c,a,b)):
        assert (x,y) not in fano
        fano[x,y]=(1,z)
        fano[y,x]=(-1,z)
fib=(1,3,5,2,6,7,4)
back={c:i for i,c in enumerate(fib)}
assert all(fib[z]==fib[x]^fib[y] for (x,y),(_,z) in fano.items())
def pmul(x,y):
    if x==7:return 1,y
    if y==7:return 1,x
    if x==y:return -1,7
    return fano[x,y]
def lp(path):
    val,u=1,7
    for v in path:
        s,u=pmul(u,v);val*=s
    return val,u

allp=[w for w in permutations(range(7)) if all(arc(a,b) for a,b in zip(w,w[1:]))]
assert len(allp)==189
assert Counter(lp(w) for w in allp)=={(-1,7):21,(1,7):168}
assert all((lp(w)[0]==-1)==(len({(b-a)%7 for a,b in zip(w,w[1:])})==1) for w in allp)
N={w for w in allp if w[:2]==(0,1)}
assert len(N)==9
def rev(w):
    x,y=w[-1],w[-2]
    a=pow((y-x)%7,-1,7)
    return tuple(a*(k-x)%7 for k in reversed(w))
assert all(rev(w) in N and rev(rev(w))==w and lp(rev(w))==lp(w) for w in N)
assert sum(rev(w)==w for w in N)==3
for w in sorted(N):
    print("".join(map(str,w)),lp(w)[0],"".join(map(str,rev(w))))

def lin(p):
    return all(fib[p[z]]==fib[p[x]]^fib[p[y]] for (x,y),(_,z) in fano.items())

def lifts(p):
    for s in product((1,-1),repeat=7):
        if all(t*s[z] == s[x]*s[y]*fano[p[x],p[y]][0] for (x,y),(t,z) in fano.items()):
            yield s
P={p:list(lifts(p)) for p in permutations(range(7)) if lin(p)}
assert len(P)==168 and all(len(l)==8 for l in P.values())
assert sum((1,)*7 in l for l in P.values())==21
assert {p for p,l in P.items() if (1,)*7 in l} == {
  tuple((a*x+b)%7 for x in range(7)) for a in Q for b in range(7)
}
# V4 displacement c=001. Every nonidentity involution lift negates e_c.
c=1
for f in (2,4,6):
    g=tuple(back[x^c if (x&f).bit_count()%2 else x] for x in fib)
    L=[s for s in P[g] if all(s[i]*s[g[i]]==1 for i in range(7))]
    assert len(L)==4 and all(s[back[c]]==-1 for s in L)
print('monomial lifts:',sum(map(len,P.values())),'pure:',sum((1,)*7 in l for l in P.values()),'local V4 verified')
for K in range(1,51):
    for L in range(1,K+1):
        D=abs(2**K-3**L)
        assert D%3!=0
        if D%24==7 and D<10000 and all(D%d for d in range(2,int(D**.5)+1)):
            assert (K,L,D)==(4,2,7)
assert [3+2**k for k in (1,2,3)] == [5,7,11]
print('raw-gap controls pass')
