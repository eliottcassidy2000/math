from fractions import Fraction as Q
from itertools import combinations
import math

def cage(k):
    if k == 3:
        labels = tuple(combinations(range(5), 2))
        N = [{j for j,B in enumerate(labels) if not set(A)&set(B)} for A in labels]
    elif k == 7:
        N = [set() for _ in range(50)]
        def add(u,v):
            N[u].add(v);N[v].add(u)
        for i in range(5):
            for x in range(5):
                add(5*i+x,5*i+(x+1)%5)
                add(25+5*i+x,25+5*i+(x+2)%5)
                for j in range(5):
                    add(5*i+x,25+5*j+(x+i*j)%5)
    else:
        raise ValueError(k)
    return N

def mask(S):
    return sum(1<<i for i in S)

def cliques(rows, j, P, prefix=()):
    if j == 0:
        yield prefix
    while j and P.bit_count() >= j:
        lo=P&-P; P^=lo
        i=lo.bit_length()-1
        yield from cliques(rows,j-1,P&rows[i],prefix+(i,))

def test(k,q):
    N=cage(k);n=k*k+1
    assert all(len(R)==k for R in N)
    assert all(len(N[i]&N[j])+(j in N[i])==1 for i in range(n) for j in range(i))
    full=(1<<n)-1
    def test_parity(G):
        rows=[mask(r) for r in G]
        return all(((rows[i]&rows[j]).bit_count()+(j in G[i]))%2==1
                   for i in range(n) for j in range(n))
    assert test_parity(N)
    found=False
    for a,b,c,d in combinations(range(n),4):
        for (x,y,z,w) in ((a,b,c,d),(a,c,b,d),(a,d,b,c)):
            if y not in N[x] or w not in N[z] or z in N[x] or w in N[y]:
                continue
            NN=[row.copy() for row in N]
            for a0,b0 in ((x,y),(z,w),(x,z),(y,w)):
                NN[a0].symmetric_difference_update({b0})
                NN[b0].symmetric_difference_update({a0})
            def is_moore(G):
                return all(len(G[i]&G[j])+(j in G[i])==1 for i in range(n) for j in range(i))
            if not is_moore(NN):
                assert not test_parity(NN)
                found=True
                break
        if found:
            break
    assert found
    rows=[full^((1<<i)|mask(N[i])) for i in range(n)]
    m=(q-1)*k+1
    S=set(next(cliques(rows,m,full)))
    V=sorted(set(range(n))-S)
    D=[N[i]&S for i in V]
    assert all(len(B)==q for B in D)
    assert len(D)==len({frozenset(B) for B in D})
    for x,y in combinations(S,2):
        assert sum({x,y} <= B for B in D)==1
    R=[(N[v]|{v})&set(V) for v in V]
    for v,W in zip(V,R):
        assert len(W)==m//q
        assert set.union(*(N[z]&S for z in W))==S
    for i,j in combinations(range(len(V)),2):
        v,w=V[i],V[j]
        assert len(R[i]&R[j])==(2 if w in N[v] else 0 if D[i]&D[j] else 1)
    rows=[mask(j for j,W in enumerate(D) if j!=i and B&W) for i,B in enumerate(D)]
    families=list(cliques(rows,k,(1<<len(V))-1))
    nstar=0;nplane=0;trades=set()
    for fam in families:
        assert set.union(*(R[i] for i in fam))==set(V)
        if set.intersection(*(D[i] for i in fam)):
            nstar+=1
            continue
        nplane+=1
        W=set.union(*(D[i] for i in fam))
        assert len(W)==k
        new=frozenset((S-W)|{V[j] for j in fam})
        assert len(new)==m and all(not (N[i]&new) for i in new)
        trades.add(new)
    if k==7:
        pos={x:i for i,x in enumerate(sorted(S))}
        M=[mask(pos[x] for x in B) for B in D]
        ker=[x for x in range(1,1<<len(S)) if all((x & b).bit_count()%2==0 for b in M)]
        assert len(ker)==15
        # Each point gives a distinct nonzero evaluation functional on that 4-space.
        labels={mask(j for j,w in enumerate(ker) if w>>i&1) for i in range(len(S))}
        assert len(labels)==15 and 0 not in labels
        # All displayed block equations x+y+z=0, so this is PG(3,2).
    count_other=sum(len(set(C)&S)==m-k for C in cliques(
      [full^((1<<i)|mask(N[i])) for i in range(n)],m,full))
    assert count_other==len(trades)
    assert (nstar,nplane)==({3:(4,4),7:(15,15)}[k])
    print(k,"maximum coclique",sorted(S),"stars",nstar,"plane trades",nplane)
    if trades:
        print("one translated coclique",sorted(min(trades,key=lambda C: sorted(C))))
    return nstar,nplane

def variance_checks():
    k,q,m,b=57,8,400,2850
    assert k==q*q-q+1 and (q-1)*k+1==m and b*q==k*m
    assert n_check(k,n=m+b) and (k-1)*m//4==5600
    assert math.ceil(Q(5600-2850,16-1))==184
    assert math.ceil(Q(5600,16))==350
def n_check(k,n):
    return n==k*k+1

if __name__=="__main__":
    variance_checks()
    test(3,2)
    test(7,3)
