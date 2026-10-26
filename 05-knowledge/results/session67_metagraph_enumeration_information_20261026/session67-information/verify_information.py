#!/usr/bin/env python3
"""Exact independent census: all labelled tournaments and all C3 inputs n<=5.
No inherited code imports, no third-party libraries. Canonical actions are brute
permutations, positive constraints use independent parity propagation, and
triangle tests are checked by explicit directed cycles.
"""
from collections import Counter,defaultdict
from fractions import Fraction
from itertools import combinations,permutations,product
import json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def require(ok,why):
    if not ok: raise RuntimeError(why)

def graph(n,b):
    A=[[0]*n for _ in range(n)]
    for q,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=(b>>q)&1; A[j][i]=1-A[i][j]
    return A

def cycle_h(A,triples):
    return sum((A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i])<<q for q,(i,j,k) in enumerate(triples))

def compile_solutions(n,h,edges,triples):
    idx={e:q for q,e in enumerate(edges)}; m=len(edges)
    adj=[[] for _ in edges]
    for q,(i,j,k) in enumerate(triples):
        if (h>>q)&1:
            a,b,c=idx[i,j],idx[j,k],idx[i,k]
            for u,v,d in [(a,b,0),(a,c,1)]:
                adj[u].append((v,d));adj[v].append((u,d))
    comp=[None]*m; offset=[0]*m; p=0
    for e in range(m):
        if comp[e] is not None: continue
        comp[e]=p; stack=[e]
        while stack:
            u=stack.pop()
            for v,d in adj[u]:
                if comp[v] is None:
                    comp[v]=p;offset[v]=offset[u]^d;stack.append(v)
                elif offset[v]!=(offset[u]^d): return set(),False,p
        p+=1
    out=set()
    for z in range(1<<p):
        xs=[((z>>comp[e])&1)^offset[e] for e in range(m)]
        good=True
        for q,(i,j,k) in enumerate(triples):
            if not (h>>q)&1:
                aa,bb,cc=xs[idx[i,j]],xs[idx[j,k]],1-xs[idx[i,k]]
                if aa==bb==cc: good=False;break
        if good: out.add(sum(x<<e for e,x in enumerate(xs)))
    return out,True,p

def census(n):
    edges=list(combinations(range(n),2));tri=list(combinations(range(n),3));m=len(edges)
    perms=list(permutations(range(n)));eidx={e:i for i,e in enumerate(edges)};tidx={t:i for i,t in enumerate(tri)}
    # Each permutation is a new-label -> old-label coordinate pullback.
    ep=[];tp=[]
    for p in perms:
        ep.append([(eidx[tuple(sorted((p[i],p[j])))],int(p[i]>p[j])) for i,j in edges])
        tp.append([tidx[tuple(sorted(p[i] for i in t))] for t in tri])
    def act(b,emap):return sum((((b>>old)&1)^flip)<<new for new,(old,flip) in enumerate(emap))
    def acth(h,tmap):return sum(((h>>old)&1)<<new for new,old in enumerate(tmap))
    hcanon={h:min(acth(h,p) for p in tp) for h in range(1<<len(tri))}
    class_mult=Counter();fibres=defaultdict(set);rows={};hclasses=defaultdict(set);curl_checks=0
    for b in range(1<<m):
        A=graph(n,b);h=cycle_h(A,tri);fibres[h].add(b)
        images=[act(b,p) for p in ep];key=min(images);aut=images.count(b)
        class_mult[key]+=1;hclasses[hcanon[h]].add(key)
        s=[sum(row) for row in A];r=[2*x-(n-1) for x in s]
        S=[[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
        c=lambda i,j,k:S[i][j]+S[j][k]+S[k][i]
        for q,(i,j,k) in enumerate(tri):
            require((c(i,j,k)**2-1)//8==((h>>q)&1),'cycle/curl identity')
        for i,j in edges:
            require(n*S[i][j]==r[i]-r[j]+sum(c(i,j,k) for k in range(n)),'signed curl + score reconstruction')
            curl_checks+=1
        c3=h.bit_count()
        require(c3==math.comb(n,3)-sum(math.comb(x,2) for x in s),'score/c3 identity')
        rows[b]=(key,aut)
    input_count=0;positive_consistent=0;parity_only_false=[]
    for h in range(1<<len(tri)):
        sols,ok,p=compile_solutions(n,h,edges,tri);input_count+=1;positive_consistent+=ok
        require(sols==fibres.get(h,set()),'phase solution fibre not exact')
        if ok and not sols:parity_only_false.append(h)
    orbit_fibres=[]
    for h in sorted(hclasses):
        sols=fibres[h];auts=[i for i,p in enumerate(tp) if acth(h,p)==h]
        groups={min(act(b,ep[i]) for i in auts) for b in sols}
        fixsum=sum(sum(act(b,ep[i])==b for b in sols) for i in auts)
        require(fixsum==len(auts)*len(groups),'Burnside fibre')
        require(len(groups)==len(hclasses[h]),'C3 fibre orbit count')
        orbit_fibres.append({'hypergraph_code':h,'positive_triples':[list(t) for q,t in enumerate(tri) if (h>>q)&1], 'labelled_realizations':len(sols),'aut_h_order':len(auts),'unlabelled_tournament_fibre':len(groups),'fixed_point_sum':fixsum})
    total=1<<m;entropy=-sum((v/total)*math.log2(v/total) for v in class_mult.values())
    rhs=m-math.log2(math.factorial(n))+sum(math.log2(aut) for key,aut in rows.values())/total
    require(abs(entropy-rhs)<1e-12,'class Shannon identity')
    for key,mult in class_mult.items():require(mult*rows[key][1]==math.factorial(n),'orbit stabilizer')
    require(sum(x['unlabelled_tournament_fibre'] for x in orbit_fibres)==len(class_mult),'fibre family sum')
    return {'n':n,'labelled_tournaments':total,'A000568':len(class_mult),'unlabelled_realizable_C3_hypergraphs':len(hclasses),'all_C3_inputs_checked':input_count,'positive_consistent_inputs':positive_consistent,'positive_consistent_unrealizable_inputs':len(parity_only_false),'curl_edge_reconstruction_checks':curl_checks,'class_entropy_uniform_labelled':entropy,'class_entropy_formula':rhs,'uniform_unlabelled_entropy':math.log2(len(class_mult)),'fibre_counts':orbit_fibres}

def hostile_checks():
    # Canonical sorted pair convention: 111 is transitive; cut at vertex1 ->010 cyclic.
    a=graph(3,7);b=graph(3,2)
    require(cycle_h(a,[(0,1,2)])==0 and cycle_h(b,[(0,1,2)])==1,'cut-switch hostile')
    require(((7).bit_count()%2)==((2).bit_count()%2),'GF2 parity invariant')
    n=5;edges=list(combinations(range(n),2));tri=list(combinations(range(n),3))
    h=sum(1<<tri.index(t) for t in [(0,1,4),(0,2,3)])
    sols,ok,p=compile_solutions(n,h,edges,tri)
    require(ok and p==6 and not sols,'negative NAE obstruction')
    # For a uniform latent ranking and p=1/2 every T has identical likelihood.
    # True mutual information exactly zero; old log n!-E log H surrogate positive at n3.
    wrong=math.log2(6)-math.log2(3)/4
    require(wrong>2,'ranking-entropy hostile')
    return {'transitive_bits':7,'cut_switched_cyclic_bits':2,'common_GF2_triangle_parity':1,'NAE_hostile_hypergraph':[[0,1,4],[0,2,3]],'NAE_hostile_phase_count':p,'independent_noisy_channel_true_mutual_information_bits':0,'old_path_entropy_mutual_information_surrogate_n3':wrong}

if __name__=='__main__':
    out={'status':'PASS','scope':'All labelled tournaments and every possible C3-hypergraph input for each n=0,...,5; no random cases. No claim above n=5 from finite tests.','implementation':'stdlib-only direct permutation oracle + independent parity graph solver; require checks survive python -O','census':[census(n) for n in range(6)],'hostiles':hostile_checks()}
    expected=[1,1,1,2,4,12]
    require([r['A000568'] for r in out['census']]==expected,'A000568 control')
    dst=ROOT/'verification.json';dst.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','rows':[{k:v for k,v in r.items() if k!='fibre_counts'} for r in out['census']],'hostiles':out['hostiles']},indent=2,sort_keys=True))
