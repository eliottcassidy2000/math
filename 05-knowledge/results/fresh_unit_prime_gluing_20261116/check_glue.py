#!/usr/bin/env python3
from itertools import permutations
from math import factorial
def graph(n, code):
    A = [[0]*n for _ in range(n)]
    bit=0
    for j in range(n):
        for i in range(j):
            A[i][j] = code>>bit & 1
            A[j][i] = 1-A[i][j]
            bit += 1
    return A

def count(G, V):
    return sum(all(G[u][v] for u,v in zip(pi,pi[1:])) for pi in permutations(V))

def runs(G,V,w):
    counts=[0]*(len(V)+2)
    V=tuple(v for v in V if v!=w)
    for p in permutations(V):
        optional=[]
        required=[]
        for i,(x,y) in enumerate(zip((-1,)+p,p+(-1,))):
            insert=(x==-1 or G[x][w]) and (y==-1 or G[w][y])
            plain = x==-1 or y==-1 or G[x][y]
            if insert and plain:
                optional.append(i)
            elif insert:
                required.append(i)
            elif not plain:
                break
        else:
            from math import comb
            for i in range(len(optional)+1):
                counts[len(required)+i]+=comb(len(optional),i)
    return tuple(counts[1:])

if __name__=='__main__':
    for n,m,u,v,a,b,expected_r in [
            (6,1105,4,5,13,13,5),
            (7,470924,1,4,23,29,13),
            (7,393231,1,3,17,23,9),
            (7,1053296,0,6,17,19,9),
            (7,1615491,2,5,23,31,13)]:
        G,vs=graph(n,m),set(range(n))
        A,B,W=vs-{u},vs-{v},vs-{u,v}
        h,p,q,r=(count(G,z) for z in (vs,A,B,W))
        assert sorted([p,q])==[a,b] and h==p+q-1 and r==expected_r
        print(f'{n=} {m=} {p=} {q=} {h=} {r=} mixed={h-p-q+r}')
        for w in sorted(W):
            lens=len(vs)+1
            cs,ca,cb=[(s+(0,)*lens)[:lens] for s in (runs(G,V,w) for V in (vs,A,B))]
            assert cs[0]==h and ca[0]==p and cb[0]==q
            if tuple(x-y-z for x,y,z in zip(cs,ca,cb))==(-1,)+(0,)*(lens-1):
                print('all lengths',w,cs,ca,cb)

def paired_runs(G,V,x,y):
    from collections import Counter
    V=tuple(v for v in V if v not in (x,y))
    xy=(x,y) if G[x][y] else (y,x)
    ret=Counter()
    for p in permutations(V):
        out=Counter({(0,0):1})
        for a,b in zip((-1,)+p,p+(-1,)):
            choices=[]
            for mid in [(),(x,),(y,),xy]:
                q=tuple(z for z in (a,)+mid+(b,) if z!=-1)
                if all(G[u][v] for u,v in zip(q,q[1:])):
                    choices.append((int(x in mid),int(y in mid)))
            nxt=Counter()
            for f,k in out.items():
                for i,j in choices:
                    nxt[(f[0]+i,f[1]+j)]+=k
            out=nxt
        ret.update({k:v for k,v in out.items() if min(k)>0})
    return ret

if __name__=='__main__':
    from itertools import combinations
    for n,m,u,v in [(6,1105,4,5),(7,470924,1,4),(7,393231,1,3),(7,1053296,0,6),(7,1615491,2,5)]:
        G,vs=graph(n,m),set(range(n))
        A,B,W=vs-{u},vs-{v},vs-{u,v}
        for x,y in combinations(sorted(W),2):
            ps,pa,pb=[paired_runs(G,V,x,y) for V in (vs,A,B)]
            if all(ps[k]-pa[k]-pb[k]==-int(k==(1,1)) for k in ps.keys()|pa.keys()|pb.keys()|{(1,1)}):
                print('paired',m,x,y,dict(ps),dict(pa),dict(pb))
            if (m,x,y)==(393231,0,2):
                keys=[(1,1),(2,1),(1,2),(2,2)]
                assert set(ps)==set(pa)==set(pb)==set(keys)
                assert [[g[k] for k in keys] for g in (ps,pa,pb)]==[[39,20,29,16],[23,12,14,8],[17,8,15,8]]
