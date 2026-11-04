#!/usr/bin/env python3
"""Exact path-complement operator in Walsh coordinates. Standard library only."""
import itertools,json,math
from fractions import Fraction as Q
from pathlib import Path
from check import maps,orbit_data,matmul,require,comps

def path_data(n,es):
    ix={e:i for i,e in enumerate(es)};out=[]
    for p in itertools.permutations(range(n)):
        support=0;positive=0
        for u,v in zip(p,p[1:]):
            e=ix[tuple(sorted((u,v)))];support|=1<<e
            if u<v:positive|=1<<e
        out.append((p,support,positive))
    return out

def directed_paths(mask,data):
    return [(p,s,up) for p,s,up in data if mask&s==up]

def all_support_check(n):
    es,ps=maps(n);N=len(es);M=1<<N;paths=path_data(n,es);undirected=[row for row in paths if row[0]<row[0][::-1]]
    # For every Fourier input h, form D chi_h by direct labelled path complements,
    # then Walsh-transform it by the elementary butterfly algorithm.
    for h in range(M):
        vals=[]
        for t in range(M):
            z=0
            for p,s,up in directed_paths(t,paths):
                u=t^((M-1)^s);z+=-1 if (h&(~u)).bit_count()%2 else 1
            vals.append(z)
        # Our characters are product x_e=(-1)^(|h|-h dot mask), not (-1)^(h dot mask).
        span=1
        while span<M:
            for start in range(0,M,2*span):
                for j in range(start,start+span):a,b=vals[j],vals[j+span];vals[j],vals[j+span]=a+b,a-b
            span*=2
        for g in range(M):
            observed=Q((-1)**g.bit_count()*vals[g],M);f=g^h
            expected=Q(0)
            if f.bit_count()%2==0:
                for p,s,up in undirected:
                    if f&~s:continue
                    sg=(-1)**((h&~s).bit_count()+(f&~up).bit_count())
                    expected+=Q(sg,2**(n-2))
            require(observed==expected,('path-Walsh',n,g,h,observed,expected))
    return M*M

def run(n):
    es,ps,tc,gc,ts,gs,copies,aut=orbit_data(n);N=len(es);M=1<<N;eg=sorted(copies);paths=path_data(n,es)
    rt=[tc.count(i) for i in range(len(ts))];rg=[len(copies[g]) for g in eg]
    W=[[sum(sg*(-1 if (g&~t).bit_count()%2 else 1) for g,sg in copies[e].items()) for e in eg] for t in ts]
    D=[[0]*len(ts) for t in ts];checked_pairs=0
    for i,t in enumerate(ts):
        a=math.factorial(n)//rt[i]
        for p,s,up in directed_paths(t,paths):
            u=t^((M-1)^s);j=tc[u];D[i][j]+=1
            require(((~(t^u))&(M-1))==s,'common path support')
            # On each oriented representative pair, check common automorphisms directly.
            common=[]
            for q,d,iv,rd in ps:
                def move(x):
                    z=rd
                    for e in range(N):
                        if x>>e&1:z^=d[e]
                    return z
                if move(t)==t and move(u)==u:common.append(q)
            require(common==[tuple(range(n))],'pair stabilizer is not trivial');checked_pairs+=1
        require(all(x%a==0 for x in D[i]),'incidence multiplicity not divisible by Aut(T)')
        if n>=3:require(D[i][i]%(2*a)==0,'self incidence parity')
    require(all(rt[i]*D[i][j]==rt[j]*D[j][i] for i in range(len(ts)) for j in range(len(ts))),'path adjacency reversibility')
    DW=matmul(D,W)
    C=[[sum(Q(rt[i]*W[i][j]*DW[i][k],M*rg[j]) for i in range(len(ts))) for k in range(len(eg))] for j in range(len(eg))]
    require(matmul(W,C)==DW,'invariant Walsh operator')
    require(all(C[j][k]==0 for j,g in enumerate(eg) for k,h in enumerate(eg) if (g.bit_count()+h.bit_count())%2),'parity mixing')
    return {'n':n,'ordered_representative_pairs_checked':checked_pairs,'labelled_Walsh_entries_checked':all_support_check(n) if n<=4 else 0,'tournament_masks':ts,'even_graph_masks':eg,'path_complement_adjacency':D,'path_counts':[sum(row) for row in D],'Walsh_basis_operator':[[str(v) for v in row] for row in C],'diagonal':all(C[i][j]==0 for i in range(len(C)) for j in range(len(C)) if i!=j),'checks':'passed'}

if __name__=='__main__':
    rows=[]
    for n in range(2,7):
        r=run(n);rows.append(r);print({k:r[k] for k in ['n','ordered_representative_pairs_checked','labelled_Walsh_entries_checked','diagonal']},flush=True)
    Path(__file__).with_name('path_results.json').write_text(json.dumps(rows,indent=2)+'\n')
