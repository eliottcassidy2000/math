#!/usr/bin/env python3
"""Independent exact census and signed-Walsh coupling audit. Standard library only.
Masks: lexicographic pairs i<j; graph 1 means present, tournament 1 means i->j.
No repository imports. No asserts, so all checks remain active under python -O.
"""
import argparse, json, math, itertools, hashlib
from fractions import Fraction as Q
from collections import Counter, defaultdict
from pathlib import Path

def require(ok, message):
    if not ok: raise RuntimeError(message)

def sign_perm(p):
    return -1 if sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2 else 1

def maps(n):
    es=list(itertools.combinations(range(n),2)); ix={e:i for i,e in enumerate(es)}
    out=[]
    for p in itertools.permutations(range(n)):
        dest=[]; inv=0; inv_dest=0
        for i,(u,v) in enumerate(es):
            a,b=p[u],p[v]; j=ix[tuple(sorted((a,b)))]; dest.append(1<<j)
            if a>b: inv |=1<<i; inv_dest|=1<<j
        out.append((p,dest,inv,inv_dest))
    return es,out

def move(mask,dest):
    r=0
    while mask:
        bit=mask&-mask; mask-=bit; r |=dest[bit.bit_length()-1]
    return r

def orbit_data(n):
    es,ps=maps(n); N=len(es); size=1<<N
    tc=[-1]*size; gc=[-1]*size; ts=[]; gs=[]; copies={}; aut={}
    for a in range(size):
        if tc[a]<0:
            vals={move(a,d)^r for p,d,i,r in ps}; t=len(ts); ts.append(a)
            for b in vals: tc[b]=t
        if gc[a]<0:
            vals={}; autos=[]
            for p,d,iv,r in ps:
                b=move(a,d); sg=-1 if (a&iv).bit_count()%2 else 1
                vals.setdefault(b,set()).add(sg)
                if b==a: autos.append((p,d,iv,r))
            g=len(gs); gs.append(a)
            for b in vals:gc[b]=g
            aut[a]=autos
            if all(len(s)==1 for s in vals.values()): copies[a]={b:next(iter(s)) for b,s in vals.items()}
    return es,ps,tc,gc,ts,gs,copies,aut

def rank(matrix):
    a=[[Q(x) for x in row] for row in matrix]; r=0
    for c in range(len(a[0]) if a else 0):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]; z=a[r][c]; a[r]=[v/z for v in a[r]]
        for i in range(r+1,len(a)):
            z=a[i][c]
            if z:a[i]=[x-z*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

def matmul(A,B):
    bt=list(zip(*B)); return [[sum(a*b for a,b in zip(row,col)) for col in bt] for row in A]

def matching(support):
    match={}
    def augment(i,seen):
        for j in support[i]:
            if j in seen:continue
            seen.add(j)
            if j not in match or augment(match[j],seen):match[j]=i;return True
        return False
    for i in range(len(support)):require(augment(i,set()),'Hall matching failed')
    return [next(j for j,k in match.items() if i==k) for i in range(len(support))]

def birkhoff(B):
    R=[row[:] for row in B]; parts=[]
    while any(R[0]):
        match=matching([[j for j,x in enumerate(row) if x] for row in R])
        q=min(R[i][j] for i,j in enumerate(match));parts.append((q,match))
        for i,j in enumerate(match):R[i][j]-=q
    require(all(all(x==0 for x in row) for row in R),'Birkhoff remainder')
    require(sum(q for q,p in parts)==1,'Birkhoff weights')
    return parts

def comps(n,mask,es):
    adj=[set() for _ in range(n)]
    for i,(u,v) in enumerate(es):
        if mask>>i&1:adj[u].add(v);adj[v].add(u)
    left=set(range(n));cs=[]
    while left:
        seen={min(left)};todo=list(seen)
        while todo:
            u=todo.pop()
            for v in adj[u]-seen:seen.add(v);todo.append(v)
        cs.append(seen);left-=seen
    return cs,adj

def homology_character_check(n,g,es,autos):
    """Independent H1 action via spanning-forest fundamental cycles and chord coords."""
    cs,adj=comps(n,g,es); edges=[i for i in range(len(es)) if g>>i&1]
    parent=list(range(n))
    def root(u):
        while parent[u]!=u:u=parent[u]
        return u
    tree=[];chords=[]
    for e in edges:
        u,v=es[e];a,b=root(u),root(v)
        if a!=b:parent[a]=b;tree.append(e)
        else:chords.append(e)
    tadj=[[] for _ in range(n)]
    for e in tree:
        u,v=es[e];tadj[u].append((v,e,1));tadj[v].append((u,e,-1))
    cycles=[]
    for ch in chords:
        u,v=es[ch]; stack=[(v,{ch:1},set([v]))]; found=None
        while stack:
            a,vec,seen=stack.pop()
            if a==u:found=vec;break
            for b,e,sg in tadj[a]:
                if b not in seen:stack.append((b,dict(vec,**{})|{e:sg},seen|{b}))
        require(found is not None,'missing fundamental cycle');cycles.append(found)
    def determinant(A):
        if not A:return Q(1)
        a=[[Q(v) for v in row] for row in A];det=Q(1)
        for j in range(len(a)):
            k=next((k for k in range(j,len(a)) if a[k][j]),None)
            if k is None:return Q(0)
            if j!=k:a[k],a[j]=a[j],a[k];det=-det
            z=a[j][j];det*=z
            for k in range(j+1,len(a)):
                t=a[k][j]/z
                for l in range(j+1,len(a)):a[k][l]-=t*a[j][l]
        return det
    ex={e:i for i,e in enumerate(edges)};co={v:i for i,c in enumerate(cs) for v in c};checked=0
    for p,d,iv,rd in autos:
        mat=[]
        for cyc in cycles:
            image={}
            for e,c in cyc.items():image[d[e].bit_length()-1]=c*(-1 if iv>>e&1 else 1)
            mat.append([image.get(e,0) for e in chords])
        dh=determinant(mat); sv=sign_perm(p)
        se=sign_perm([ex[d[e].bit_length()-1] for e in edges])
        sc=sign_perm([co[p[min(c)]] for c in cs]); eps=-1 if (g&iv).bit_count()%2 else 1
        require(eps==sv*se*sc*dh,('homology character',n,g,p,eps,sv,se,sc,dh));checked+=1
    return checked

def run(n,do_homology):
    es,ps,tc,gc,ts,gs,copies,aut=orbit_data(n);N=len(es);size=1<<N;nf=math.factorial(n)
    eg=sorted(copies);rT=[tc.count(i) for i in range(len(ts))];rG=[len(copies[g]) for g in eg]
    require(len(ts)==len(eg),'equinumerosity')
    M=[[sum(sg*(-1 if (g&(~t)).bit_count()%2 else 1) for g,sg in copies[e].items()) for e in eg] for t in ts]
    # Column weighted orthogonality. This includes every zero/off-diagonal check.
    for j in range(len(eg)):
        for k in range(j,len(eg)):
            require(sum(r*M[i][j]*M[i][k] for i,r in enumerate(rT))==(size*rG[j] if j==k else 0),'Walsh orthogonality')
    B=[[Q(rT[i]*M[i][j]**2,size*rG[j]) for j in range(len(eg))] for i in range(len(ts))]
    require(all(sum(row)==1 for row in B),'row stochastic')
    require(all(sum(col)==1 for col in zip(*B)),'column stochastic')
    A=[[0]*len(ts) for t in ts]
    for i,t in enumerate(ts):
        for k in range(N):A[i][tc[t^(1<<k)]]+=1
    require(all(rT[i]*A[i][j]==rT[j]*A[j][i] for i in range(len(ts)) for j in range(len(ts))),'detailed balance')
    AM=matmul(A,M)
    require(all(AM[i][j]==(N-2*eg[j].bit_count())*M[i][j] for i in range(len(ts)) for j in range(len(eg))),'spectral eigenvectors')
    # Full transition matrices for powers 0..4, not merely diagonal tests.
    Ap=[[int(i==j) for j in range(len(ts))] for i in range(len(ts))]
    for q in range(5):
        for i in range(len(ts)):
            for j in range(len(ts)):
                spectral=Q(rT[j],size)*sum(Q((N-2*g.bit_count())**q*M[i][k]*M[j][k],rG[k]) for k,g in enumerate(eg))
                require(Ap[i][j]==spectral,('full transition',n,q,i,j))
            require(Ap[i][i]==sum(B[i][k]*(N-2*g.bit_count())**q for k,g in enumerate(eg)),'return spectral measure')
        Ap=matmul(Ap,A)
    converse=[tc[(size-1)^t] for t in ts]
    require(all(B[i]==B[converse[i]] for i in range(len(ts))),'converse probability loss')
    require(all(M[converse[i]][j]==(-1)**g.bit_count()*M[i][j] for i in range(len(ts)) for j,g in enumerate(eg)),'converse amplitudes')
    # Rational noise t=1/3: each coordinate correlation z=1/3, exact independent-flip law.
    z=Q(1,3);pf=(1-z)/2
    for i,t in enumerate(ts):
        direct=[Q(0)]*len(ts)
        for mask in range(size):
            h=(t^mask).bit_count();direct[tc[mask]]+=pf**h*(1-pf)**(N-h)
        for j in range(len(ts)):
            spectral=Q(rT[j],size)*sum(Q(M[i][k]*M[j][k],rG[k])*z**g.bit_count() for k,g in enumerate(eg))
            require(direct[j]==spectral,('noise kernel',n,i,j))
    # Canonical support matching and exact Birkhoff decomposition.
    parts=birkhoff(B)
    # Probability equality fibres and degree-aggregated spectral-measure equality fibres.
    fibres=defaultdict(list);rf=defaultdict(list)
    for i,row in enumerate(B):
        fibres[tuple(row)].append(i)
        rr=tuple(sum(row[k] for k,g in enumerate(eg) if g.bit_count()==d) for d in range(N+1));rf[rr].append(i)
    counts=Counter(g.bit_count() for g in eg)
    homchecks=0
    if do_homology:
        for g in gs:homchecks+=homology_character_check(n,g,es,aut[g])
    out={'n':n,'tournament_classes':len(ts),'all_graph_classes':len(gs),'automorphism_even_graph_classes':len(eg),'even_order_automorphism_graph_classes':sum(len(aut[g])%2==0 for g in gs),'labelled_even_graphs':sum(rG),'even_graph_edge_counts':dict(sorted(counts.items())),'self_converse_tournament_classes':sum(i==j for i,j in enumerate(converse)),'coupling_rank':rank(B),'coupling_distinct_rows':len(fibres),'return_spectral_distinct_rows':len(rf),'converse_merged_classes':len({min(i,j) for i,j in enumerate(converse)}),'birkhoff_terms':len(parts),'homology_automorphisms_checked':homchecks,'checks':'all exact checks passed','tournament_masks':ts,'tournament_orbit_sizes':rT,'even_graph_masks':eg,'even_graph_orbit_sizes':rG,'converse_indices':converse,'Walsh_matrix':M,'coupling':[[str(x) for x in row] for row in B],'full_flip_adjacency':A,'birkhoff_decomposition':[{'weight':str(q),'permutation':p} for q,p in parts],'coupling_nontrivial_fibres':[[ts[i] for i in v] for v in fibres.values() if len(v)>1],'return_nontrivial_fibres':[[ts[i] for i in v] for v in rf.values() if len(v)>1]}
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--max-n',type=int,default=6);ap.add_argument('--out',default='results.json');args=ap.parse_args()
    result={'convention':'lex pairs; tournament bit1 means lower beats higher; graph bit1 means present','scope':'exhaustive unlabelled classes, all labels and all permutations through max_n; exact rational identities','orders':[]}
    for n in range(1,args.max_n+1):
        row=run(n,do_homology=True);result['orders'].append(row)
        print(json.dumps({k:row[k] for k in ['n','tournament_classes','even_order_automorphism_graph_classes','coupling_rank','coupling_distinct_rows','return_spectral_distinct_rows','converse_merged_classes','birkhoff_terms','homology_automorphisms_checked']}),flush=True)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
