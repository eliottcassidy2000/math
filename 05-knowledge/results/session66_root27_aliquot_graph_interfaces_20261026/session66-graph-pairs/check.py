#!/usr/bin/env python3
"""Standalone, deterministic standard-library checks; no removable assertions."""
import itertools as it
import json, random
from functools import lru_cache

def need(x,m):
    if not x: raise RuntimeError(m)

def tournament(n,b):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(it.combinations(range(n),2)):
        a[i][j]=(b>>k)&1;a[j][i]=1-a[i][j]
    return a

def key(a): return tuple(a[i][j] for i,j in it.combinations(range(len(a)),2))
def relabel(a,p):
    b=[[0]*len(a) for _ in a]
    for i in range(len(a)):
        for j in range(len(a)): b[p[i]][p[j]]=a[i][j]
    return b

def canon(a):
    return min((key(relabel(a,p)) for p in it.permutations(range(len(a)))),default=())

def add(g,u,v,c=1):g[u][v]|=c;g[v][u]|=c

def joint(a):
    n=len(a);g=[[0]*(7*n) for _ in range(7*n)]
    # Tree center0, arms [1], [2,3], [4,5,6]. Tail1, head3.
    for i in range(n):
        for x,y in ((0,1),(0,2),(2,3),(0,4),(4,5),(5,6)):add(g,7*i+x,7*i+y,1)
    for i in range(n):
        for j in range(n):
            if a[i][j]:add(g,7*i+1,7*j+3,2)
    return g

def graph(a):
    n=len(a);g=[[0]*(9*n) for _ in range(9*n)]
    # Core C0,T1,H2, their leaf banks [3],[4,5],[6,7,8].
    for i in range(n):
        for x,y in ((0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),(2,8)):add(g,9*i+x,9*i+y)
    for i in range(n):
        for j in range(n):
            if a[i][j]:add(g,9*i+1,9*j+2)
    return g

def decode_graph(g):
    N=len(g);adj=[{v for v in range(N) if g[u][v]} for u in range(N)]
    leaves={u for u in range(N) if len(adj[u])==1};cores=set(range(N))-leaves
    typ={u:len(adj[u]&leaves) for u in cores}
    centers=sorted(u for u in cores if typ[u]==1)
    tails=[];heads=[]
    for c in centers:
        tails.append(next(u for u in adj[c]&cores if typ[u]==2))
        heads.append(next(u for u in adj[c]&cores if typ[u]==3))
    return [[int(g[t][h]>0) for h in heads] for t in tails]

def decode_joint(g):
    N=len(g);adj=[{v for v in range(N) if g[u][v]&1} for u in range(N)]
    centers=sorted(u for u in range(N) if len(adj[u])==3)
    tails=[];heads=[]
    for c in centers:
        arms=[]
        for x in adj[c]:
            arm=[x];prev=c
            while len(adj[x])==2:
                y=next(v for v in adj[x] if v!=prev);prev,x=x,y;arm.append(x)
            arms.append(arm)
        bylen={len(x):x for x in arms}
        need(set(bylen)=={1,2,3},'bad tree');tails.append(bylen[1][-1]);heads.append(bylen[2][-1])
    return [[int(bool(g[t][h]&2)) for h in heads] for t in tails]

def graph_iso(g,h):
    """Independent exact colored adjacency backtracker; no codec decoder used."""
    n=len(g)
    if n!=len(h):return False
    if n==0:return True
    def degrees(a):return [tuple(row.count(c) for c in (1,2,3)) for row in a]
    ds=degrees(g)+degrees(h);palette={x:i for i,x in enumerate(sorted(set(ds)))}
    c=[palette[x] for x in ds]
    for _ in range(n):
        sig=[]
        for k,a in enumerate((g,h)):
            cs=c[k*n:(k+1)*n]
            sig.extend((cs[u],tuple(sorted((a[u][v],cs[v]) for v in range(n) if a[u][v]))) for u in range(n))
        pal={x:i for i,x in enumerate(sorted(set(sig)))};d=[pal[x] for x in sig]
        if d==c:break
        c=d
    c1,c2=c[:n],c[n:]
    if sorted(c1)!=sorted(c2):return False
    mapped={};used=set()
    def rec():
        if len(mapped)==n:return True
        choices=[]
        for u in range(n):
            if u in mapped:continue
            cand=[v for v in range(n) if v not in used and c1[u]==c2[v] and all(g[u][x]==h[v][y] for x,y in mapped.items())]
            if not cand:return False
            choices.append((len(cand),-sum(bool(g[u][x]) for x in mapped),u,cand))
        _,_,u,cand=min(choices)
        for v in cand:
            mapped[u]=v;used.add(v)
            if rec():return True
            used.remove(v);del mapped[u]
        return False
    return rec()

def hamilton(a):
    n=len(a)
    if not n:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,1<<n):
        for v in range(n):
            if s>>v&1:
                t=s^(1<<v)
                if t:dp[s][v]=sum(dp[t][u] for u in range(n) if a[u][v])
    return sum(dp[-1])

def cycles(a):
    n=len(a);ans=[]
    for k in range(3,n+1,2):
        for sub in it.combinations(range(n),k):
            for rest in it.permutations(sub[1:]):
                c=(sub[0],)+rest
                if all(a[c[i]][c[(i+1)%k]] for i in range(k)):ans.append(c)
    return ans

def omega(a):
    cs=cycles(a);g=[[0]*len(cs) for _ in cs]
    for i,j in it.combinations(range(len(cs)),2):
        if set(cs[i])&set(cs[j]):add(g,i,j)
    return g

def independence(g):
    adj=[sum((1<<j) for j in range(len(g)) if g[i][j]) for i in range(len(g))]
    @lru_cache(None)
    def f(s):
        if not s:return 1
        v=(s&-s).bit_length()-1;t=s&~(1<<v)
        return f(t)+2*f(t&~adj[v])
    return f((1<<len(g))-1)

def homogeneous(a):
    n=len(a);g=[[0]*n for _ in a]
    for i,j in it.combinations(range(n),2):
        if all(a[i][k]==a[j][k] for k in range(n) if k not in (i,j)):add(g,i,j)
    return g

def clone_compress(a):
    j=homogeneous(a);remaining=set(range(len(a)));comps=[]
    while remaining:
        c={min(remaining)};todo=list(c)
        while todo:
            u=todo.pop()
            for v in remaining:
                if j[u][v] and v not in c:c.add(v);todo.append(v)
        remaining-=c;comps.append(sorted(c,key=lambda u:-sum(a[u][v] for v in c)))
    q=[[a[c[0]][d[0]] if c!=d else 0 for d in comps] for c in comps]
    order=[v for c in comps for v in c];expanded=[]
    for i,c in enumerate(comps):
        for r,u in enumerate(c):
            row=[]
            for k,d in enumerate(comps):
                for ss,v in enumerate(d):row.append(int(r<ss) if i==k else q[i][k])
            expanded.append(row)
    need(expanded==[[a[u][v] for v in order] for u in order],'weighted clone reconstruction')
    return q,[len(c) for c in comps]

def strong(a):
    n=len(a)
    for u in range(n):
        seen={u};todo=[u]
        while todo:
            x=todo.pop()
            for y in range(n):
                if a[x][y] and y not in seen:seen.add(y);todo.append(y)
        if len(seen)<n:return False
    return True

def main():
    out={'scope':'Authored elementary proofs plus explicitly finite exact tests; no all-order census or complexity-class claim.'}
    reps={};tot=0;counts={};hist={}
    for n in range(6):
        reps[n]={};totn=0
        for b in range(1<<(n*(n-1)//2)):
            a=tournament(n,b);ca=canon(a);reps[n].setdefault(ca,a)
            need(decode_joint(joint(a))==a,'joint roundtrip')
            need(decode_graph(graph(a))==a,'single roundtrip')
            clone_compress(a)
            hp=hamilton(a)
            brute=sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in it.permutations(range(n)))
            need(hp==brute==independence(omega(a)),'OCF/direct/DP')
            totn+=1
        counts[n]=len(reps[n]);tot+=totn
    need(list(counts.values())==[1,1,1,2,4,12],'known small class counts')
    out['all_labelled_through_order5']={'tournaments':tot,'classes_by_n':counts,'roundtrip_codecs':2,'weighted_clone_reconstruction_checks':tot,'three_way_path_count_checks':tot}
    pairtests=0;positive=0;rng=random.Random(660721)
    for n,dd in reps.items():
        aa=list(dd.values())
        for i,a in enumerate(aa):
            for codec in (joint,graph):
                g=codec(a);p=list(range(len(g)));rng.shuffle(p)
                need(graph_iso(g,relabel(g,p)),'positive scrambled graph iso');positive+=1
                for j,b in enumerate(aa):
                    if j<i:continue
                    need(graph_iso(g,codec(b))==(i==j),'codec iso census');pairtests+=1
    out['independent_graph_backtracking']={'all_class_pairs_same_order_through5_two_codecs':pairtests,'scrambled_positive_cases':positive}
    # Strong spectra through six, dynamic program over every labelled tournament.
    for n in range(3,7):
        vals=set()
        for b in range(1<<(n*(n-1)//2)):
            a=tournament(n,b);hp=hamilton(a)
            need(hp not in (7,21),'gap hostile')
            if strong(a):vals.add(hp)
        hist[n]=sorted(vals)
    need([min(hist[n]) for n in range(3,7)]==[3,5,9,15],'strong minima')
    out['finite_strong_spectra']=hist
    defect_checks=0
    for n in range(5):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits)
            for p in it.permutations(range(n)):
                defect=sum(a[i][j]!=a[p[i]][p[j]] for i,j in it.combinations(range(n),2))
                for codec,m in ((joint,7),(graph,9)):
                    g=codec(a);lp=[m*p[u//m]+u%m for u in range(m*n)]
                    delta=sum(g[i][j]!=g[lp[i]][lp[j]] for i,j in it.combinations(range(m*n),2))
                    need(delta==2*defect,'exact lifted defect')
                    defect_checks+=1
    out['lifted_graph_defect_checks']=defect_checks
    # Forward/backward pair loses its order at the smallest possible n=3.
    trans=[[int(i<j) for j in range(3)] for i in range(3)]
    t1=[r[:] for r in trans];t1[0][1],t1[1][0]=0,1
    t2=[r[:] for r in trans];t2[0][2],t2[2][0]=0,1
    def fb(a):
        g=[[0]*len(a) for _ in a]
        for i,j in it.combinations(range(len(a)),2):add(g,i,j,1 if a[i][j] else 2)
        return g
    need(hamilton(t1)==1 and hamilton(t2)==3 and graph_iso(fb(t1),fb(t2)),'order-loss hostile')
    # First non-self-converse pair, order4 source plus directed triangle.
    a=[[0]*4 for _ in range(4)]
    for j in (1,2,3):a[0][j]=1
    for i,j in ((1,2),(2,3),(3,1)):a[i][j]=1
    for i,j in it.combinations(range(1,4),2):a[j][i]=1-a[i][j]
    b=[list(r) for r in zip(*a)]
    need(canon(a)!=canon(b),'converse hostile noniso')
    need(hamilton(a)==hamilton(b)==3,'converse H')
    need(graph_iso(omega(a),omega(b)) and graph_iso(homogeneous(a),homogeneous(b)),'natural pair hostile')
    ga=joint(a);gb=joint(b)
    la=lambda g,c:[[int(bool(x&c)) for x in row] for row in g]
    need(graph_iso(la(ga,1),la(gb,1)) and graph_iso(la(ga,2),la(gb,2)) and not graph_iso(ga,gb),'independent alignment hostile')
    need(not graph_iso(graph(a),graph(b)),'repaired independent codec')
    k3=[[int(i!=j) for j in range(3)] for i in range(3)]
    p4=[[int(abs(i-j)==1) for j in range(4)] for i in range(4)]
    need(independence(k3)==7 and independence(p4)==21,'arbitrary-graph gap hostile')
    # Every simple graph invariant under C3 is empty or complete.
    invariant=[]
    for bits in range(8):
        g=[[0]*3 for _ in range(3)]
        for k,(i,j) in enumerate(it.combinations(range(3),2)):
            if bits>>k&1:add(g,i,j)
        if relabel(g,[1,2,0])==g:invariant.append(bits)
    need(invariant==[0,7],'same-vertex naturality obstruction')
    # Canonical clone paths preserve weights; J alone loses all quotient directions.
    tt3=trans;j=homogeneous(tt3)
    need(j[0][1] and j[1][2] and not j[0][2],'homogeneity not equivalence')
    def involution_defect(a,p):return sum(a[i][j]!=a[p[i]][p[j]] for i,j in it.combinations(range(len(a)),2))
    need(min(involution_defect(tournament(4,b),[1,0,3,2]) for b in range(64))==2,'two-transposition minimum')
    out['hostiles']={'forward_backward_min_order':3,'natural_pair_and_alignment_min_order':4,'converse_scores':[[sum(r) for r in a],[sum(r) for r in b]],'arbitrary_graph_I_at2':[7,21],'C3_invariant_graph_masks':invariant,'two_transposition_min_defect':2,'TT3_homogeneity_is_not_equivalence':True}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':main()
