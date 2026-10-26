#!/usr/bin/env python3
"""Portable exact controls for the reset/involution response-rank theorem.
No repository imports, network calls, dependencies, or removable assertions.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from itertools import combinations, permutations, product
import argparse, hashlib, json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rank(rows, modulus=None):
    basis = {}
    for row in rows:
        v = [F(x) for x in row] if modulus is None else [int(x) % modulus for x in row]
        for p, b in sorted(basis.items()):
            if v[p]:
                a = v[p]
                v = [x-a*y for x,y in zip(v,b)] if modulus is None else [(x-a*y)%modulus for x,y in zip(v,b)]
        p = next((i for i,x in enumerate(v) if x), None)
        if p is not None:
            a = 1/v[p] if modulus is None else pow(v[p], -1, modulus)
            basis[p] = [a*x for x in v] if modulus is None else [(a*x)%modulus for x in v]
    return len(basis)


def rational_determinant(rows):
    A=[[F(x) for x in row] for row in rows];n=len(A);det=F(1)
    need(all(len(row)==n for row in A),'square determinant')
    for k in range(n):
        i=next((i for i in range(k,n) if A[i][k]),None)
        if i is None:return F(0)
        if i!=k:A[k],A[i]=A[i],A[k];det=-det
        pivot=A[k][k];det*=pivot
        for j in range(k+1,n):
            if A[j][k]:
                c=A[j][k]/pivot;A[j]=[x-c*y for x,y in zip(A[j],A[k])]
    return det


def span_basis(rows):
    out=[]; r=0
    for row in rows:
        s=rank(out+[row])
        if s>r: out.append(row); r=s
    return out


def pullback(M, row):
    return [sum(a*b for a,b in zip(m,row)) for m in M]


def close_rows(rows, matrices):
    b=span_basis(rows); dims=[len(b)]
    while True:
        c=span_basis(b+[pullback(M,r) for M in matrices for r in b])
        if len(c)==len(b): return b,dims
        b=c; dims.append(len(b))


def incidence(E,q):
    return [[int((a if side==0 else b)==v) for a,b in E] for side in (0,1) for v in range(q)]


def graph_data(E,q):
    adj=defaultdict(list)
    for i,(a,b) in enumerate(E):
        adj[a].append((q+b,i)); adj[q+b].append((a,i))
    seen=set(); components=0; tree=set(); chord=[]
    tree_adj=defaultdict(list)
    for root in range(2*q):
        need(root in adj,'every observed class must occur on each shore')
        if root in seen: continue
        components+=1; seen.add(root); queue=deque([root])
        while queue:
            u=queue.popleft()
            for v,i in adj[u]:
                if v not in seen:
                    seen.add(v); queue.append(v); tree.add(i)
                    tree_adj[u].append((v,i)); tree_adj[v].append((u,i))
    chord=[i for i in range(len(E)) if i not in tree]
    cycles=[]
    for i in chord:
        a,b=E[i]; target=q+b; queue=deque([a]); prev={a:None}
        while target not in prev:
            u=queue.popleft()
            for v,j in tree_adj[u]:
                if v not in prev: prev[v]=(u,j); queue.append(v)
        path=[]; v=target
        while prev[v] is not None:
            u,j=prev[v]; path.append(j); v=u
        v=[0]*len(E); v[i]=1
        for h,j in enumerate(path): v[j]=(-1)**(h+1)
        cycles.append(v)
    B=incidence(E,q)
    need(all(all(sum(a*b for a,b in zip(row,c))==0 for row in B) for c in cycles),'alternating forest cycle in kernel')
    need(len(cycles)==len(E)-2*q+components,'cycle dimension formula')
    need(all(c[i]==1 and all(c[j]==0 for j in chord if j!=i) for i,c in zip(chord,cycles)),'chord coordinates independently repair all cycle directions')
    return components, cycles, chord


def finite_lumping(obs, matrices):
    def ids(sigs):
        d={}; return [d.setdefault(s,len(d)) for s in sigs]
    cells=ids(obs)
    while True:
        k=max(cells)+1
        sig=[]
        for i,o in enumerate(obs):
            sig.append((o,tuple(tuple(sum(M[i][j] for j,c in enumerate(cells) if c==h) for h in range(k)) for M in matrices)))
        new=ids(sig)
        if len(set(new))==len(set(cells)): return cells
        cells=new


def abstract_controls(maxq=4):
    counts=[]
    for q in range(1,maxq+1):
        orbits=[[(a,b)] if a==b else [(a,b),(b,a)] for a in range(q) for b in range(a,q)]
        count=0
        for mask in range(1,1<<len(orbits)):
            E=sorted(e for i,o in enumerate(orbits) if (mask>>i)&1 for e in o)
            if set(a for a,b in E)!=set(range(q)): continue
            count+=1; N=len(E); index={e:i for i,e in enumerate(E)}
            T=[[int(j==index[(b,a)]) for j in range(N)] for a,b in E]
            deg=Counter(a for a,b in E)
            R=[[F(int(c==a),deg[a]) for c,d in E] for a,b in E]
            K=[[int(c==a) for c,d in E] for a,b in E]
            B=incidence(E,q); base=B[:q]; c,cycles,chord=graph_data(E,q)
            expected=2*q-c
            for field in (None,2,3,5):
                need(rank(B,field)==expected,'incidence rank over Q/F2/F3/F5')
            closure, dims=close_rows(base,[T,R,K])
            need(len(closure)==expected,'all terminal traces close at marginals')
            need(dims[-1]==expected and len(dims)<=2,'at most one orbit-observation expansion')
            need(len(set(finite_lumping([a for a,b in E],[T,R])))==N,'state-lumping cells are all realized pairs')
            products=[[x*y for x,y in zip(l,r)] for l in B[:q] for r in B[q:]]
            need(rank(products)==N,'same-sample two-time products recover full pair law')
            if cycles:
                z=cycles[0]; pos=[max(x,0) for x in z]; neg=[max(-x,0) for x in z]
                need(sum(pos)==sum(neg)>0,'cycle gives equal positive masses')
                need([sum(a*b for a,b in zip(row,pos)) for row in B]==[sum(a*b for a,b in zip(row,neg)) for row in B],'cycle probability marginals agree')
        counts.append({'q':q,'symmetric_supports_with_full_shores':count})
    return counts


def presentation(n,mask):
    rows=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            forward=True if j==i+1 else bool((mask>>k)&1)
            if j!=i+1:k+=1
            if forward: rows[i]|=1<<j
            else: rows[j]|=1<<i
    return tuple(rows)


def code_at(rows,p):
    n=len(rows)
    return sum(((rows[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(combinations(range(n),2)))


def canonical(rows):
    groups=defaultdict(list)
    for v,row in enumerate(rows):groups[row.bit_count()].append(v)
    choices=[list(permutations(groups[d])) for d in sorted(groups)]
    return min(code_at(rows,tuple(v for block in blocks for v in block)) for blocks in product(*choices))


def all_canonical(rows):
    return min(code_at(rows,p) for p in permutations(range(len(rows))))


def tournament_controls(maxn=6):
    result=[]; witnesses=[]
    known=[None,1,1,2,4,12,56]
    for n in range(1,maxn+1):
        M=(n-1)*(n-2)//2; N=1<<M
        ts=[presentation(n,m) for m in range(N)]
        qs=[canonical(t) for t in ts]; labels=sorted(set(qs)); ix={q:i for i,q in enumerate(labels)}
        E=sorted(set((ix[qs[m]],ix[qs[m^(N-1)]]) for m in range(N))); q=len(labels)
        need(q==known[n],'ordinary tournament class census')
        c,cycles,chord=graph_data(E,q); B=incidence(E,q)
        expected=2*q-c
        for field in (None,2,3):need(rank(B,field)==expected,'native tournament incidence rank')
        if n<=5:
            alt=[all_canonical(t) for t in ts]; fwd=defaultdict(set); rev=defaultdict(set)
            for a,b in zip(qs,alt): fwd[a].add(b);rev[b].add(a)
            need(all(len(s)==1 for s in fwd.values()) and all(len(s)==1 for s in rev.values()),'degree-stratified versus all-permutation canonical equivalence')
        row={'n':n,'literal_states':N,'ordinary_classes':q,'pair_cells':len(E),'bipartite_components':c,'terminal_linear_rank':expected,'hidden_pair_law_rank':len(cycles),'normalized_terminal_affine_dimension':expected-1}
        if cycles and not witnesses:
            z=cycles[0]; support=[{'pair':[labels[a],labels[b]],'coefficient':z[i], 'literal_mask':next(m for m in range(N) if (ix[qs[m]],ix[qs[m^(N-1)]] )==(a,b))} for i,(a,b) in enumerate(E) if z[i]]
            witnesses.append({'n':n,'encoding':'upper-triangle forward bits after nondecreasing degree block order, minimum within blocks','signed_cycle':support})
        result.append(row)
    return result,witnesses


def native_trace_hostile():
    n=5; N=64; ts=[presentation(n,m) for m in range(N)]
    q=[all_canonical(t) for t in ts]; fibres=defaultdict(list)
    for m,a in enumerate(q):fibres[a].append(m)
    plus=[F(int(m in (15,21)),2) for m in range(N)]
    minus=[F(int(m in (47,5)),2) for m in range(N)]
    def marginal(law):
        return {a:sum(law[m] for m in ms) for a,ms in fibres.items()}
    def step(law,op):
        if op=='T':return [law[m^63] for m in range(N)]
        masses=marginal(law)
        return [masses[q[m]]/len(fibres[q[m]]) for m in range(N)]
    checked=0
    stack=[('',plus,minus)]
    while stack:
        word,a,b=stack.pop(); checked+=1
        need(marginal(a)==marginal(b),'actual tournament terminal law equality')
        if len(word)<8:
            for op in 'TR':stack.append((word+op,step(a,op),step(b,op)))
    event=lambda law:sum(law[m] for m in range(N) if (q[m],q[m^63])==(4,10))
    need(event(plus)==F(1,2) and event(minus)==0,'actual tournament two-time event distinguishes')
    return {'word_alphabet':{'T':'literal complement','R':'uniform fresh unrestricted re-marking, implemented on presentation fibres'},'all_words_max_length':8,'words_checked':checked,'positive_law_masks':[15,21],'negative_law_masks':[47,5],'each_mass':'1/2','joint_event_full_canonical_pair':[4,10],'joint_event_probabilities':[str(event(plus)),str(event(minus))],'presentations':[{'mask':m,'rows':ts[m],'target_rows':ts[m^63],'full_min_pair':[q[m],q[m^63]]} for m in (15,21,47,5)]}


def test_return_controls():
    out=[]
    for n in range(1,6):
        N=1<<((n-1)*(n-2)//2); ts=[presentation(n,m) for m in range(N)]
        ps=list(permutations(range(n))); guards=[]; test_count=0
        def remark(T,p):
            if T is None:return None
            if not all((T[p[i]]>>p[i+1])&1 for i in range(n-1)):return None
            return tuple(sum(((T[p[i]]>>p[j])&1)<<j for j in range(n)) for i in range(n))
        for p in ps:
            inverse=tuple(p.index(i) for i in range(n)); row=[]
            for T in ts:
                moved=remark(T,p); restored=remark(moved,inverse); legal=moved is not None
                need(restored==(T if legal else None),'native inverse remark restores state iff guard passes')
                row.append(int(legal)); test_count+=1
            guards.append(row)
        need(len({tuple(row[x] for row in guards) for x in range(N)})==N,'native guard signatures separate')
        fields=[None,2,3,5,7]
        raw={('Q' if f is None else 'F'+str(f)):rank(guards,f) for f in fields}
        need(all(rank([[1]*N]+guards,f)==raw['Q' if f is None else 'F'+str(f)] for f in fields),'constant does not increase tested raw ranks')
        gm=[sum(bit<<j for j,bit in enumerate(row)) for row in guards]
        basis={};selected=[];max_tests=0
        for depth, candidates in [(1,gm),(2,[a&b for a in gm for b in gm])]:
            for orig in candidates:
                v=orig
                while v:
                    pivot=v.bit_length()-1
                    if pivot in basis:v^=basis[pivot]
                    else:basis[pivot]=v;selected.append(orig);max_tests=depth;break
        selected_rows=[[(bits>>j)&1 for j in range(N)] for bits in selected]
        pair={('Q' if f is None else 'F'+str(f)):rank(selected_rows,f) for f in fields}
        determinant=rational_determinant(selected_rows)
        need(abs(determinant)==1,'native at-most-two-product unimodular certificate')
        need(all(r==N for r in pair.values()),'at most two native guard products give full rank in tested fields')
        out.append({'n':n,'valid_states':N,'permutations':len(ps),'inverse_return_instances':test_count,'raw_guard_rank':raw,'at_most_two_guard_product_rank':pair,'selected_basis_max_test_returns':max_tests,'selected_actual_product_determinant':str(determinant),'selected_actual_product_row_bits':[hex(v) for v in selected]})
    need(out[-1]['raw_guard_rank']=={'Q':44,'F2':43,'F3':44,'F5':44,'F7':44},'exact field-dependent n5 raw guard rank')
    return out


def degree_six_and_fourier_controls():
    n=6;N=1024;ts=[presentation(n,m) for m in range(N)];ps=list(permutations(range(n)))
    guards=[sum(int(all(T[p[i]]>>p[i+1]&1 for i in range(n-1)))<<m for m,T in enumerate(ts)) for p in ps]
    basis={};chosen=[]
    def add(v,word):
        original=v
        while v:
            pivot=v.bit_length()-1
            if pivot in basis:v^=basis[pivot]
            else:basis[pivot]=v;chosen.append(word);return
    add((1<<N)-1,[])
    for i,g in enumerate(guards):add(g,[i])
    raw=len(basis); examined=0
    for i,g in enumerate(guards):
        for j in range(i,len(guards)):
            add(g&guards[j],[i,j]);examined+=1
            if len(basis)==N:break
        if len(basis)==N:break
    need(raw==272 and len(basis)==1024,'n6 exact F2 rank and depth-two recovery')
    # Rebuild the selected ACTUAL product rows (not reduced XOR rows) independently.
    selected=[]
    for word in chosen:
        v=(1<<N)-1
        for i in word:v&=guards[i]
        selected.append(v)
    pivots={}
    for v in selected:
        while v:
            p=(v&-v).bit_length()-1
            if p in pivots:v^=pivots[p]
            else:pivots[p]=v;break
    need(len(pivots)==N,'opposite-pivot replay of actual n6 certificate')
    n=5;N=64;ts=[presentation(n,m) for m in range(N)];ps=list(permutations(range(n)))
    G=[[int(all(T[p[i]]>>p[i+1]&1 for i in range(n-1))) for T in ts] for p in ps]
    H=[[sum(v*(-1)**((m&S).bit_count()) for m,v in enumerate(row)) for S in range(N)] for row in G]
    visible=[j for j in range(N) if any(row[j] for row in H)]
    edges=[(i,j) for i in range(n) for j in range(i+2,n)]; forests=[]
    for mask in range(N):
        degrees=[0]*n;parent=list(range(n));cyclic=False
        def root(x):
            while parent[x]!=x:x=parent[x]
            return x
        for j,(a,b) in enumerate(edges):
            if mask>>j&1:
                degrees[a]+=1;degrees[b]+=1;u,v=root(a),root(b)
                if u==v:cyclic=True
                parent[u]=v
        if max(degrees)<=2 and not cyclic:forests.append(mask)
    need(visible==forests and len(visible)==46,'n5 Fourier support equals46 linear forests')
    need(rank(H)==44,'n5 two extra dependencies among46 visible Fourier columns')
    return {'n6':{'states':1024,'raw_rank_F2':raw,'two_guard_product_rank_F2':len(basis),'rational_rank_lower_bound':1024,'unordered_pairs_examined':examined,'selected_lexicographic_permutation_index_words':chosen,'certificate':'Actual0/1 product rows, independent opposite-pivot Gaussian elimination overF2; determinant odd, hence nonzero overQ.'},'n5_fourier':{'free_edges':edges,'visible_character_masks':visible,'linear_forest_count':46,'absent_character_count':18,'visible_column_rank_Q':44,'additional_visible_dependencies_Q':2}}


def hostiles():
    E=list(product(range(2),repeat=2)); B=incidence(E,2)
    T=[[int(e==(b,a)) for e in E] for a,b in E]
    Rbad=[[int(e==((1,1) if (a,b)==(1,1) else (0,0))) for e in E] for a,b in E]
    closure,_=close_rows(B[:2],[T,Rbad])
    need(rank(B)==3 and len(closure)==4,'reset not q-dependent defeats marginal closure')
    even=[F(1,2),0,0,F(1,2)]; odd=[0,F(1,2),F(1,2),0]
    need([sum(a*b for a,b in zip(row,even)) for row in B]==[sum(a*b for a,b in zip(row,odd)) for row in B],'square cycle same terminal response')
    need(even[3]!=odd[3],'same-sample trace separates parity laws')
    # A symmetric forest support: compatible mass equations need not admit positivity.
    forest=[(0,0),(0,1),(1,0)]; signed=[-1,1,1]
    marg=[sum(a*b for a,b in zip(row,signed)) for row in incidence(forest,2)]
    need(marg==[0,1,0,1],'nonnegative marginals with only negative edge realization')
    # Partial-affine domain hostile (inherited S71, independent tiny control).
    def split(word,x):
        for L in word:
            if not (0<x<L): return None
            x=L-x
        return x
    need(split([2,2],F(5,2)) is None and split([3,3],F(5,2))==F(5,2),'formal affine identity does not preserve domain')
    return {'four_cycle_terminal_rank':3,'four_cycle_pair_law_rank':4,'nonreset_terminal_rank':len(closure),'parity_law_full_trace_probability':[str(even[3]),str(odd[3])],'forest_feasibility_hostile':{'edges':forest,'unique_signed_edge_law':signed,'marginals':marg},'domain_hostile':{'words':[[2,2],[3,3]],'input':'5/2','outputs':['failure','5/2']}}


def main():
    p=argparse.ArgumentParser();p.add_argument('--corrupt',choices=['rank','trace','cone','return']);p.add_argument('--max-n',type=int,default=6);args=p.parse_args()
    abstract=abstract_controls(); native,witnesses=tournament_controls(args.max_n); hostile=hostiles(); returns=test_return_controls()
    if args.corrupt=='rank': need(native[-1]['terminal_linear_rank']==native[-1]['pair_cells'],'intentional corruption: pair states are not response rank')
    if args.corrupt=='trace': need(hostile['four_cycle_terminal_rank']==hostile['four_cycle_pair_law_rank'],'intentional corruption: terminal responses are not joint trace responses')
    if args.corrupt=='cone': need(all(x>=0 for x in hostile['forest_feasibility_hostile']['unique_signed_edge_law']),'intentional corruption: nonnegative marginals do not guarantee feasibility')
    if args.corrupt=='return': need(returns[-1]['raw_guard_rank']['Q']==returns[-1]['valid_states'],'intentional corruption: raw guards are not full test-return response')
    print(json.dumps({'abstract':abstract,'native_tournaments':native,'native_cycle_witness':witnesses,'native_trace_hostile':native_trace_hostile(),'test_return_controls':returns,'degree_and_fourier_controls':degree_six_and_fourier_controls(),'hostiles':hostile},sort_keys=True,indent=2))

if __name__=='__main__':main()
