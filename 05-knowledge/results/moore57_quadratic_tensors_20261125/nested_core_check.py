"""Exact small Moore controls and core-filter replay on the fourteen-seat chart."""
from collections import Counter
from itertools import combinations
import sys
from planted_stars import orbit_chart, full_graph
from three_stars import SLOTS


def hs():
    A=[set() for v in range(50)]
    def add(x,y):
        A[x].add(y); A[y].add(x)
    for i in range(5):
        for x in range(5):
            add(5*i+x,5*i+(x+1)%5)
            add(25+5*i+x,25+5*i+(x+2)%5)
            for j in range(5):
                add(5*i+x,25+5*j+(x+i*j)%5)
    return A


def mask(xs):
    return sum(1 << x for x in xs)


def pentagons(A):
    def dfs(p):
        if len(p)==5:
            if p[1]<p[-1] and p[0] in A[p[-1]]:
                yield mask(p)
            return
        for v in A[p[-1]]:
            if v>p[0] and v not in p:
                yield from dfs(p+[v])
    X=sorted({z for i in range(len(A)) for z in dfs([i])})
    assert len(X)==1260
    return X


def petersen_copies(A):
    L=pentagons(A)
    adj=[mask(xs) for xs in A]
    roots={P:tuple(i for i in range(50) if P>>i&1) for P in L}
    decompositions=Counter()
    for P,Q in combinations(L,2):
        if not P&Q and all((adj[i]&Q).bit_count()==1 for i in roots[P]) and \
                       all((adj[i]&P).bit_count()==1 for i in roots[Q]):
            decompositions[P|Q]+=1
    assert len(decompositions)==525 and set(decompositions.values())=={6}
    at_pent=Counter()
    for S in decompositions:
        for P in L:
            if P&S==P:
                at_pent[P]+=1
    assert set(at_pent.values())=={5} and len(at_pent)==1260
    S=set(range(5))|set(range(25,30))
    assert mask(S) in decompositions
    ports={x:{v for v in range(50) if v not in S and {x}==A[v]&S} for x in S}
    assert all(len(P)==4 for P in ports.values())
    assert sum(map(len,ports.values()))==40
    for x,y in combinations(S,2):
        assert all(len(A[u]&ports[y])==(0 if y in A[x] else 1) for u in ports[x])
    print("Hoffman–Singleton has 525 Petersen supports and 5 copies per pentagon.")


def orbits(rows,n):
    unseen=set(range(n)); ans=[]
    while unseen:
        x=min(unseen); frontier=[x]; seen={x}; unseen.remove(x)
        for u in frontier:
            for f in rows:
                v=f[u]
                if v not in seen:
                    seen.add(v);unseen.remove(v);frontier.append(v)
        ans.append(seen)
    return ans


def triangle_type(A):
    count=0
    for w in range(50):
        heads=sorted(A[w])
        fibers={i:sorted(A[i]-{w}) for i in heads}
        match={}
        for i,j in combinations(heads,2):
            match[i,j]={x:next(iter(A[x]&set(fibers[j]))) for x in fibers[i]}
            match[j,i]={v:u for u,v in match[i,j].items()}
        for i,j,h in combinations(heads,3):
            X=fibers[i]
            f=[X.index(match[j,i][match[h,j][match[i,h][x]]]) for x in X]
            assert sorted(map(len,orbits([f],6)))==[2,2,2]
            count+=1
    assert count==1750
    print("All 1750 rooted branch triangles have type 2^3 in Hoffman–Singleton.")


def chart_filter():
    M,stars=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
    assert full_graph(M,stars)==(13231,240)
    census=Counter(tuple(sorted(map(len,orbits([M[0,28],M[h,0],M[h,28]],56)))) for h in SLOTS)
    assert census=={(56,):42,(14,42):4,(14,14,28):4,(14,14,14,14):3,(28,28):1}
    print("Orbit-size census with all three completed branches",sorted(census.items()))
    expect=[(14,19,20,29,37,45),(3,10,19,48),(5,6,10,14,20,45)]
    for e,(i,j) in zip(expect,combinations((56,0,28),2)):
        L={}
        for h in SLOTS:
            f=[M[j,i][M[h,j][t]] for t in M[i,h]]
            if sum(len(o)==2 for o in orbits([f],56))>=3:
                L[h]=f
        assert tuple(L)==e
        num=0
        for hs0 in combinations(L,5):
            assert len(orbits([L[h] for h in hs0],56))==1
            num+=1
        print((i,j),"extra heads",tuple(L),"five-head trials",num)
    print("Any root-containing Hoffman–Singleton copy uses at most one completed branch.")


from collections import defaultdict
from pathlib import Path
import json

def graph(C,stars):
    A=[set() for _ in range(3250)]
    def j(i,k):
        A[i].add(k);A[k].add(i)
    def L(b,x):return 58+56*b+x
    for b in range(57):
        j(0,1+b)
        for x in range(56):j(1+b,L(b,x))
    for (i,k),f in C.items():
        if i<k:
            for x,y in enumerate(f):j(L(i,x),L(k,y))
    for (i,x),E in stars.items():
        for k,y in E:j(L(i,x),L(k,y))
    return A
def color_classes(A,mark=()):
    c=[int(i in mark) for i in range(len(A))]
    for r in range(len(A)):
        lab={};d=[]
        for i,R in enumerate(A):
            key=c[i],tuple(sorted(c[j] for j in R))
            d.append(lab.setdefault(key,len(lab)))
        assert len(lab)>=len(set(c))
        if len(set(c))==len(lab):
            P=defaultdict(list)
            for i,k in enumerate(d):P[k].append(i)
            return d,dict(P)
        c=d
def autgroups(A):
    lab,P=color_classes(A)
    print('WL',len(P),Counter(map(len,P.values())))
    if max(map(len,P.values()))>2:return
    K={k for k,v in P.items() if len(v)==2}; comp=[]
    while K:
        k=min(K);K.remove(k);Q=[k]
        for g in Q:
            v=P[g][0]
            N=Counter(lab[u] for u in A[v])
            for h,m in N.items():
                if h in K and m==1:K.remove(h);Q.append(h)
        comp.append(Q)
    print('automorphism group rank',len(comp),'component sizes',sorted(map(len,comp),reverse=True)[:20])
    assert len(comp)==1
    sig=list(range(len(A)))
    if len(comp)==1:
        for g in comp[0]:
            x,y=P[g];sig[x]=y;sig[y]=x
        assert all({sig[w] for w in R}==A[sig[v]] for v,R in enumerate(A))
        print("root?",sig[0],"heads",[(i-1,sig[i]-1) for i in range(1,58) if sig[i]!=i])
    return sig
def legal(A,v,w):
    if v==w:return False
    reached=set(A[v])|{v}
    for u in tuple(reached):reached.update(A[u])
    if w in reached:return False
    # New edge would close a cycle of length <=4 iff there's an old route of length <=3.
    return not (A[w]&reached)


def rigidity_and_near_patch(C,S):
    A=graph(C,S)
    sig=autgroups(A)
    assert sig is not None and sum(i!=x for i,x in enumerate(sig))==3194
    assert [(i-1,sig[i]-1) for i in range(1,58) if sig[i]!=i]==[(27,49),(49,27)]
    G=[R.copy() for R in A]
    assert legal(G,115,177)
    G[115].add(177);G[177].add(115)
    assert len(color_classes(G)[1])==3250
    E=json.loads(Path(__file__).with_name('rigid_51_contacts.json').read_text())
    assert len(E)==51
    assert len({h for h,y in E})==51 and {h for h,y in E}==set(SLOTS)-{1,13,33}
    for h,y in E:
        w=58+56*h+y
        assert legal(A,115,w)
        A[115].add(w);A[w].add(115)
    assert len(A[115])==55 and sum(map(len,A))//2==13282
    assert max(map(len,A))==57 and sum(len(R)==57 for R in A)==240
    assert len(color_classes(A)[1])==3250
    print("Retained every old edge and all 14 full neighborhoods; added 51 contacts at (1,1); nested chart is rigid.")


def board_packet():
    C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
    ref=dict(json.loads(Path(__file__).with_name('rigid_51_contacts.json').read_text()))
    b,x=1,1
    centers=(56,0,28)
    R=[h for h in SLOTS if h!=b]
    used={h:{y for choices in S.values() for k,y in choices if k==h} for h in R}
    avail=[set(range(56))-({C[b,c][x]}|{C[d,c][C[b,d][x]] for d in centers if d!=c}) for c in centers]
    print(len(R))
    for h in R:
        row=[]
        for y in sorted(set(range(56))-used[h]):
            if h==2 and y!=7:continue
            z=[C[h,c][y] for c in centers]
            if all(z[i] in avail[i] for i in range(3)):row.append([y,*z])
        assert h not in ref or any(a[0]==ref[h] for a in row)
        print(ref.get(h,-1),len(row))
        for a in row:print(*a)


if __name__=='__main__':
    if '--board' in sys.argv:
        board_packet()
        raise SystemExit
    A=hs()
    assert all(len(A[i]&A[j])+(j in A[i])==1 for i in range(50) for j in range(i))
    petersen_copies(A)
    triangle_type(A)
    chart_filter()


    C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
    rigidity_and_near_patch(C,S)
