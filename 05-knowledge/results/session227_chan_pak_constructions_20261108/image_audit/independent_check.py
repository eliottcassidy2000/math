#!/usr/bin/env python3
"""Independent Session227 graft-image audit. No producer imports.
Simple paths are enumerated by recursive extension; no subset DP is used.
All exact checks use explicit exceptions and survive Python optimization.
"""
import itertools as it
import json
from math import gcd
from fractions import Fraction as F

TOTAL=0

def check(ok, msg):
    global TOTAL
    TOTAL+=1
    if not ok:
        raise RuntimeError(msg)

def tournaments(n):
    pairs=list(it.combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        arcs=set()
        for k,(i,j) in enumerate(pairs):
            arcs.add((i,j) if bits>>k&1 else (j,i))
        yield n,arcs

def from_rows(rows):
    n=len(rows)
    check(all(not (rows[i]>>i&1) for i in range(n)), 'loops')
    check(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i,j in it.combinations(range(n),2)), 'tournament')
    return n,{(i,j) for i in range(n) for j in range(n) if rows[i]>>j&1}

def path_tables(graph):
    n,arcs=graph
    h=[0]*(1<<n); h[0]=1
    starts=[[0]*(1<<n) for _ in range(n)]
    ends=[[0]*(1<<n) for _ in range(n)]
    def visit(first,last,mask):
        h[mask]+=1; starts[first][mask]+=1; ends[last][mask]+=1
        for v in range(n):
            if not mask>>v&1 and (last,v) in arcs:
                visit(first,v,mask|1<<v)
    for v in range(n): visit(v,v,1<<v)
    return h,starts,ends

def subs(mask):
    k=mask
    while True:
        yield k
        if k==0: return
        k=(k-1)&mask

def all_stats(graph):
    n,arcs=graph; h,ss,ee=path_tables(graph); top=(1<<n)-1
    pairs=[sum(h[a]*h[m^a] for a in subs(m)) for m in range(1<<n)]
    output=[]
    for r in range(n):
        alpha=sum(ss[r][m]*h[top^m] for m in range(1<<n))
        delta=sum(ee[r][m]*h[top^m] for m in range(1<<n))
        T=sum(ss[r][m]*pairs[top^m] for m in range(1<<n))
        output.append(dict(D=h[top],E=ee[r][top],U=alpha,d=h[top^(1<<r)],
                           H=h[top],Q=pairs[top],s=ss[r][top],e=ee[r][top],
                           alpha=alpha,delta=delta,T=T,q=h[top^(1<<r)],b=pairs[top^(1<<r)]))
    return output

def state(s): return tuple(s[k] for k in ('D','E','U','d'))
def slacks(x):
    D,E,U,d=x
    return E,d-E,D-d,U-D-d+E

def mv(M,v): return tuple(sum(a*b for a,b in zip(row,v)) for row in M)
def det(M):
    n=len(M); ans=0
    for p in it.permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        prod=sign
        for i in range(n):prod*=M[i][p[i]]
        ans+=prod
    return ans

def matrix(s,n):
    H,Q,S,e,a,de,T,q,b=(s[k] for k in ('H','Q','s','e','alpha','delta','T','q','b'))
    return ((H,-H,Q-H,2*H-Q),(e,-e,de-e,e-de),
            (a,S-a,2*T-a,3*a-S-2*T),
            (q,-q,b-q,2*q-b) if n>1 else (1,0,0,0))

def strong(g):
    n,a=g
    for r in range(n):
        seen={r}
        while True:
            nxt=seen|{v for u,v in a if u in seen}
            if nxt==seen:break
            seen=nxt
        if len(seen)!=n:return False
    return True

def cyclic(m):
    n=m+2; u=m;v=m+1
    a={(i,j) for i in range(m) for j in range(i+1,m)}
    a|={(i,u) for i in range(m)}|{(u,v)}|{(v,i) for i in range(m)}
    return n,a

def transitive(n):return n,{(i,j) for i in range(n) for j in range(i+1,n)}

def graft(a,r,b,br):
    n,aa=a;m,bb=b
    arcs=aa|{(i+n,j+n) for i,j in bb}
    arcs|={(j+n,r) for j in range(m)}|{(i,j+n) for i in range(n) if i!=r for j in range(m)}
    return (n+m,arcs),n+br

def module_status(g,M):
    n,arcs=g
    losses=[]
    for w in range(n):
        if w in M:continue
        relation=[(w,v) in arcs for v in M]
        if any(relation) and not all(relation):return False,[]
        if not any(relation):losses.append(w)
    return True,losses

def certs(g,r,b,br):
    n,a=g;m,ba=b;out=[]
    if n<=m:return out
    for tail in it.combinations([v for v in range(n) if v!=r],m-1):
        M={r,*tail};mod,loss=module_status(g,M)
        if not mod or len(loss)!=1:continue
        nonroots=[v for v in range(m) if v!=br]
        for image in it.permutations(tail):
            f={br:r,**dict(zip(nonroots,image))}
            if all(((i,j) in ba)==((f[i],f[j]) in a) for i,j in it.permutations(range(m),2)):
                out.append((M,loss[0]));break
    return out

def inv4(y):
    A,B,C,G=map(F,y)
    return C+5*B-(A+13*G)/2,C+4*B-(5*A+11*G)/4,G-B,(A+3*G-8*B)/8

def main():
    roots=0;small_graphs=[]
    for n in range(1,6):
        for g in tournaments(n):
            statistics=all_stats(g)
            for r,s in enumerate(statistics):
                x=state(s);D,E,U,d=x
                check(min(slacks(x))>=0,'physical cone')
                check(D%2==d%2==1 and (U-E)%2==0,'physical affine lattice')
                if D==1:
                    check(U>0 and U&(U-1)==0,'H=1 rooted cover power of two')
                    check(not any((i,j) in g[1] and (j,k) in g[1] and (k,i) in g[1]
                                  for i,j,k in it.permutations(range(n),3)),'H=1 transitive')
                M=matrix(s,n);Delta=s['q']*s['Q']-s['H']*s['b']
                generators=((2,0,0,0),(0,-1,1,0),(0,0,2,0),(0,0,0,2))
                increments=[sum(a*b for a,b in zip(M[0],v)) for v in generators]
                check(gcd(*increments)==gcd(2*s['H'],s['Q']),'strongest affine scalar congruence')
                check(mv(M,(1,1,1,1))[0]==s['H'],'scalar affine basepoint')
                if n>1:check(det(M)==s['s']*s['e']*Delta,'det factorization')
                if n<=4:small_graphs.append((g,r,s,M))
                roots+=1
    check(roots==5405,'root universe')
    # Independent exact ray coordinate transformation.
    rays=((1,1,1,1),(1,0,2,1),(1,0,1,0),(0,0,1,0))
    check(abs(det(rays))==1,'unimodular cone rays')
    for j,r in enumerate(rays):check(slacks(r)==tuple(int(i==j) for i in range(4)),'ray slacks')
    for n in range(1,11):
        s=all_stats(transitive(n))
        for r in range(n):
            check(state(s[r])==(1,int(r==n-1),2**(n-1-r),1),'all transitive root states')
    family=[]
    for m in range(1,10):
        g=cyclic(m);stats=all_stats(g)
        H=2**m+1;Q=2*3**m+2**(m+2)-2
        check(strong(g),'cyclic family strong')
        if m>=3 and m%2:
            check(H%3==Q%3==0,'full rank family scalar divisibility barrier')
        for r,s in enumerate(stats):
            check((s['H'],s['Q'])==(H,Q),'cyclic H,Q')
            if r>=m:
                expected=2*(3**m+2**m-1-4**m)
                check((s['s'],s['e'],s['q'],s['b'])==(1,1,1,2**(m+1)),'singleton stats')
                check(expected==0 if m==1 else expected<0,'negative Delta')
            elif m>=2:
                expected=(2**m+4)*3**(m-1)+3*2**m
                check((s['s'],s['e'],s['q'],s['b'])==(2**(m-1-r),2**r,2**(m-1)+1,2*3**(m-1)+2**(m+1)-2),'block stats')
                check(expected>0,'positive Delta')
            else:expected=0
            check(s['q']*s['Q']-s['H']*s['b']==expected,'family Delta formula')
            check(det(matrix(s,m+2))==s['s']*s['e']*expected,'family determinant')
        check(state(stats[m])==(H,1,H,1),'cone-ray family state')
        family.append({'m':m,'H':H,'Q':Q,'singleton_delta':2*(3**m+2**m-1-4**m),
                       'block_delta':((2**m+4)*3**(m-1)+3*2**m) if m>=2 else 0})
    graphs=[(transitive(1),0),(transitive(2),0),(transitive(3),0),(cyclic(1),0),(cyclic(2),0)]
    witnesses=[state(all_stats(g)[r]) for g,r in graphs]
    check(witnesses==[(1,1,1,1),(1,0,2,1),(1,0,4,1),(3,1,3,1),(5,1,9,3)],'lattice states')
    differences=[tuple(x-y for x,y in zip(w,witnesses[0])) for w in witnesses[1:]]
    check(abs(det(differences))==8,'affine difference determinant')
    # Literal graphs for every rooted gadget through order four on four spanning sources.
    sources=[graphs[i] for i in (0,1,3,2)]
    grafts=0
    for g,r,s,M in small_graphs:
        for a,ar in sources:
            joined,jr=graft(a,ar,g,r)
            check(state(all_stats(joined)[jr])==mv(M,state(all_stats(a)[ar])),'literal graft state')
            check(bool(certs(joined,jr,g,r)),'literal graph certificate')
            grafts+=1
    check(grafts==1140,'graft universe')
    B=from_rows([8,1,3,6]);K=matrix(all_stats(B)[0],4)
    check(det(K)==-8,'B4 det')
    for x in it.product(range(-2,3),repeat=4):
        check(inv4(mv(K,x))==x,'B4 exact inverse')
        xx=inv4(x)
        check(all(v.denominator==1 for v in xx)==((x[0]+3*x[3])%8==0),'B4 complete integer image')
    z=(1,0,6,1)
    check(min(slacks(z))>=0 and z[0]%2==z[3]%2==1 and (z[2]-z[1])%2==0,'hostile hull gates')
    check(mv(K,z)==(145,31,159,37),'hostile B4 image')
    good=from_rows([18,32,3,7,14,29]);bad=from_rows([18,40,3,5,14,29])
    check(strong(good) and strong(bad),'collision strong')
    check(state(all_stats(good)[0])==state(all_stats(bad)[0])==(37,7,43,9),'collision states')
    good_certs=certs(good,0,B,0);bad_certs=certs(bad,0,B,0)
    check(good_certs==[({0,2,3,4},1)],'good exact certificates')
    check(bad_certs==[],'bad absence')
    obstruction=[]
    for tail in it.combinations(range(1,6),3):
        M={0,*tail};mod,loss=module_status(bad,M)
        check(not mod,'bad has no root-containing four-module')
        witnesses=[]
        for w in range(6):
            if w in M:continue
            beaten=sorted(v for v in M if (w,v) in bad[1]);beats=sorted(M-set(beaten))
            if beaten and beats:witnesses.append({'exterior':w,'beats':beaten,'loses_to':beats})
        check(bool(witnesses),'mixed-edge obstruction witness')
        obstruction.append({'module_candidate':sorted(M),'obstruction':witnesses[0]})
    check(inv4((37,7,43,9))==(1,0,2,1),'physical inverse collision')
    print(json.dumps({'status':'PASS','checks':TOTAL,'universe':'All 5405 labeled rooted tournaments at orders 1..5; no filters',
                      'independent_method':'Recursive enumeration of all simple directed paths; ordered subset partitions; Leibniz determinants; no producer imports',
                      'literal_graft_checks':grafts,'cyclic_family_range':'m=1..9, every root',
                      'family':family,'lattice_witnesses':[list(state(all_stats(g)[r])) for g,r in graphs],
                      'B4_input_output_grid':625,'bad_graph_module_obstructions':obstruction},indent=2,sort_keys=True))

if __name__=='__main__':main()
