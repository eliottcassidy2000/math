#!/usr/bin/env python3
"""Exact ordinary/merged incidence, independent symmetry and signed-kernel audit.
No repository producer imported. Native move complements nonpath arcs.
"""
from collections import Counter,defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,permutations,product
from math import factorial
from pathlib import Path
import argparse,json


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def direct(n):
    pairs=list(combinations(range(n),2));loc={p:i for i,p in enumerate(pairs)}
    top=(1<<len(pairs))-1;ps=list(permutations(range(n)))
    def arc(code,u,v):return bool(code>>loc[tuple(sorted((u,v)))]&1)==(u<v)
    def move(code,p):return sum(int(arc(code,p[u],p[v]))<<i for i,(u,v) in enumerate(pairs))
    owner={};aut={}
    for code in range(top+1):
        if code in owner:continue
        orbit={move(code,p) for p in ps}
        need(not any(c in owner for c in orbit),'disjoint ordinary orbits')
        owner.update({c:code for c in orbit});aut[code]=factorial(n)//len(orbit)
    ordinary=sorted(aut);R={c:owner[top^c] for c in ordinary}
    C={c:Counter() for c in ordinary};blue={c:Counter() for c in ordinary};H={}
    for c in ordinary:
        hp=0
        for p in ps:
            if not all(arc(c,p[i],p[i+1]) for i in range(n-1)):continue
            hp+=1
            pathmask=sum(1<<loc[tuple(sorted((p[i],p[i+1])))] for i in range(n-1))
            target=owner[c^top^pathmask]
            C[c][target]+=1
            rev=[None]*n
            for u,v in zip(p,p[::-1]):rev[u]=v
            if move(c,rev)==top^c:blue[c][target]+=1
        H[c]=hp
        for counts in (C[c],blue[c]):
            for t,v in list(counts.items()):
                need(v%aut[c]==0,'free source path orbit')
                counts[t]=v//aut[c]
    merged=sorted({min(c,R[c]) for c in ordinary})
    A={c:Counter() for c in merged};B={c:Counter() for c in merged}
    for c in ordinary:
        for raw,out in ((C,A),(blue,B)):
            for t,v in raw[c].items():out[min(c,R[c])][min(t,R[t])]+=v
    for c in ordinary:
        need(sum(C[c].values())==H[c]//aut[c],'ordinary degree')
        for t in ordinary:
            need(C[c][t]==C[t][c] and blue[c][t]==blue[t][c],'symmetric ordinary matrix')
            need(C[R[c]][R[t]]==C[c][t],'R commutation')
        if R[c]!=c:need(not blue[c],'blue confined to SC')
    nodes=[{'code':c,'converse':R[c],'H':H[c],'aut':aut[c],'SC':c==R[c],
            'degree':sum(C[c].values())} for c in ordinary]
    return {'n':n,'ordinary_nodes':nodes,'merged_nodes':[dict(next(x for x in nodes if x['code']==c),
             degree=sum(A[c].values())) for c in merged],
             'C':[[c,t,v] for c in ordinary for t,v in sorted(C[c].items())],
             'blueC':[[c,t,v] for c in ordinary for t,v in sorted(blue[c].items())],
             'A':[[c,t,v] for c in merged for t,v in sorted(A[c].items())],
             'blueA':[[c,t,v] for c in merged for t,v in sorted(B[c].items())]}


def dense(codes,triples):
    loc={c:i for i,c in enumerate(codes)};M=[[0]*len(codes) for _ in codes]
    for c,t,v in triples:M[loc[c]][loc[t]]=v
    return M


def refine(matrices,labels):
    def pack(signatures):
        vals={s:i for i,s in enumerate(sorted(set(signatures)))}
        return [vals[s] for s in signatures]
    n=len(labels)
    colors=pack([(labels[i],tuple((sum(M[i]),M[i][i]) for M in matrices)) for i in range(n)])
    history=[len(set(colors))]
    while True:
        signatures=[]
        for i in range(n):
            cells=Counter((colors[j],tuple(M[i][j] for M in matrices)) for j in range(n))
            signatures.append((colors[i],tuple(sorted(cells.items()))))
        new=pack(signatures)
        history.append(len(set(new)))
        if history[-1]==history[-2]:return new,history
        colors=new


def automorphisms(matrices,labels):
    colors,history=refine(matrices,labels);n=len(colors)
    options=[[j for j in range(n) if colors[i]==colors[j]] for i in range(n)]
    assigned={};used=set();answer=[];visits=0
    def search():
        nonlocal visits
        visits+=1
        if len(assigned)==n:
            p=tuple(assigned[i] for i in range(n))
            need(all(M[i][j]==M[p[i]][p[j]] for M in matrices for i in range(n) for j in range(n)),
                 'complete automorphism check')
            answer.append(p);return
        best=None
        for i in range(n):
            if i in assigned:continue
            valid=[j for j in options[i] if j not in used and
                   all(M[i][u]==M[j][v] for u,v in assigned.items() for M in matrices)]
            if not valid:return
            if best is None or len(valid)<len(best[1]):best=(i,valid)
        i,valid=best
        for j in valid:
            assigned[i]=j;used.add(j);search();used.remove(j);del assigned[i]
    search()
    return answer,history,visits


def gauge_solutions(D,g):
    n=len(D);assigned={};components=[]
    for i in range(n):
        if D[i][i]!=D[g[i]][g[i]]:return []
        for j in range(n):
            if abs(D[i][j])!=abs(D[g[i]][g[j]]):return []
    for root in range(n):
        if root in assigned:continue
        assigned[root]=1;members=[root];queue=[root]
        while queue:
            i=queue.pop()
            for j in range(n):
                if i==j or not D[i][j]:continue
                ratio=D[g[i]][g[j]]//D[i][j]
                target=assigned[i]*ratio
                if j in assigned:
                    if assigned[j]!=target:return []
                else:assigned[j]=target;members.append(j);queue.append(j)
        components.append(members)
    solutions=[]
    for signs in product([-1,1],repeat=len(components)):
        e=[0]*n
        for sign,component in zip(signs,components):
            for i in component:e[i]=sign*assigned[i]
        need(all(D[g[i]][g[j]]==e[i]*e[j]*D[i][j] for i in range(n) for j in range(n)), 'gauge solution')
        solutions.append(e)
    return solutions


def gauge_controls():
    K=[[0,1,1,1],[1,0,1,1],[1,1,0,-1],[1,1,-1,0]]
    checks=[]
    for twist in [(1,1,1,1),(-1,1,1,1)]:
        D=[[K[i][j]*twist[i]*twist[j] for j in range(4)] for i in range(4)]
        lifts=set();good=0
        for g in permutations(range(4)):
            solutions=gauge_solutions(D,g)
            brute=[list(e) for e in product([-1,1],repeat=4)
                   if all(D[g[i]][g[j]]==e[i]*e[j]*D[i][j] for i in range(4) for j in range(4))]
            need(sorted(solutions)==sorted(brute),'gauge BFS versus all signs')
            good+=bool(solutions)
            for e in solutions:lifts.add(tuple(2*g[i//2]+((i%2)^(e[i//2]<0)) for i in range(8)))
        C=[[0]*8 for _ in range(8)]
        for i in range(8):
            for j in range(8):
                if i//2!=j//2:C[i][j]=(3+(1 if i%2==j%2 else -1)*D[i//2][j//2])//2
        literal=set()
        for p in permutations(range(8)):
            if any(p[i^1]!=(p[i]^1) for i in range(8)):continue
            if all(C[i][j]==C[p[i]][p[j]] for i in range(8) for j in range(8)):literal.add(p)
        need(lifts==literal and good==4 and len(lifts)==8,'complete synthetic ordinary symmetry')
        checks.append({'gauge':twist,'base_automorphisms':24,'liftable_base':good,'ordinary_group_order':len(lifts)})
    return checks


def signed_analysis(row):
    oc=[x['code'] for x in row['ordinary_nodes']];mc=[x['code'] for x in row['merged_nodes']]
    info={x['code']:x for x in row['ordinary_nodes']};R={c:info[c]['converse'] for c in oc}
    C=dense(oc,row['C']);A=dense(mc,row['A']);B=dense(mc,row['blueA']);oi={c:i for i,c in enumerate(oc)}
    ch=[c for c in mc if R[c]!=c];sc=[c for c in mc if R[c]==c]
    D=[[C[oi[c]][oi[t]]-C[oi[c]][oi[R[t]]] for t in ch] for c in ch]
    need(all(D[i][j]==D[j][i] for i in range(len(ch)) for j in range(len(ch))),'signed symmetry')
    # Reconstruct every entry of C from the merged matrix and signed current.
    mi={c:i for i,c in enumerate(mc)};ci={c:i for i,c in enumerate(ch)}
    for c in oc:
        for t in oc:
            u,v=min(c,R[c]),min(t,R[t]);a=A[mi[u]][mi[v]]
            if u==R[u] and v==R[v]:expected=Fraction(a)
            elif u==R[u] or v==R[v]:expected=Fraction(a,2)
            else:
                s=1 if (c==u)==(t==v) else -1
                expected=Fraction(a,4)+s*Fraction(D[ci[u]][ci[v]],2)
            need(expected.denominator==1 and expected==C[oi[c]][oi[t]],'integral reconstruction')
    todo=set(range(len(ch)));components=[]; spanning=[]
    while todo:
        seen={min(todo)};queue=list(seen)
        while queue:
            i=queue.pop()
            for j in range(len(ch)):
                if i!=j and D[i][j] and j not in seen:
                    seen.add(j);queue.append(j)
                    spanning.append([ch[i],ch[j],C[oi[ch[i]]][oi[ch[j]]],C[oi[ch[i]]][oi[R[ch[j]]]],D[i][j]])
        todo-=seen;components.append([ch[i] for i in sorted(seen)])
    regimes={}
    for name,mat,labels in [('A',[A],[()]*len(mc)),('A_type',[A],[(c==R[c],) for c in mc]),
                          ('A_type_H',[A],[(c==R[c],info[c]['H']) for c in mc]),
                          ('colored_type_H',[B,[[A[i][j]-B[i][j] for j in range(len(mc))] for i in range(len(mc))]],
                           [(c==R[c],info[c]['H']) for c in mc])]:
        G,history,visits=automorphisms(mat,labels)
        inv=[g for g in G if all(g[g[i]]==i for i in range(len(g)))]
        regimes[name]={'group_order':len(G),'refinement_cell_counts':history,'search_nodes':visits,
                       'quotient_involutions':len(inv),'minimum_fixed_singletons':min(sum(c==R[c] and g[i]==i for i,c in enumerate(mc)) for g in inv),
                       'permutations_codes':[[mc[g[i]] for i in range(len(mc))] for g in G]}
    matrix_groups={}
    CB=dense(oc,row['blueC'])
    for regime in ['A_type','A_type_H','colored_type_H']:
        permutations_found=[];fixed=Counter();admissible=0
        for g_codes in regimes[regime]['permutations_codes']:
            gm=dict(zip(mc,g_codes));gch=[ci[gm[c]] for c in ch]
            solutions=gauge_solutions(D,gch)
            admissible+=bool(solutions)
            for epsilon in solutions:
                lifted={}
                for c in oc:
                    representative=min(c,R[c]);target=gm[representative]
                    if representative==R[representative]:lifted[c]=target
                    else:
                        flip=(c!=representative)^(epsilon[ci[representative]]<0)
                        lifted[c]=R[target] if flip else target
                need(all(C[oi[c]][oi[t]]==C[oi[lifted[c]]][oi[lifted[t]]] for c in oc for t in oc),'full ordinary lift')
                if regime=='colored_type_H':
                    need(all(CB[oi[c]][oi[t]]==CB[oi[lifted[c]]][oi[lifted[t]]] for c in oc for t in oc),'colored ordinary lift')
                if regime.endswith('_H'):need(all(info[c]['H']==info[lifted[c]]['H'] for c in oc),'H ordinary lift')
                if all(lifted[lifted[c]]==c for c in oc):fixed[sum(lifted[c]==c for c in oc)]+=1
                permutations_found.append([lifted[c] for c in oc])
        matrix_groups[regime]={'liftable_base_order':admissible,'ordinary_group_order':len(permutations_found),
                              'involution_fixed_polynomial':dict(sorted(fixed.items())),
                              'minimum_fixed':min(fixed),'permutations_codes':permutations_found}
    # Check every signed component flip as an actual ordinary matrix symmetry.
    component_checks=0
    for comp in components:
        flip=set(comp);p={c:(R[c] if min(c,R[c]) in flip else c) for c in oc}
        need(all(C[oi[c]][oi[t]]==C[oi[p[c]]][oi[p[t]]] for c in oc for t in oc),'component flip physical C')
        component_checks+=1
    bad_edges=[]
    for i,c in enumerate(ch):
        for j,t in enumerate(ch):
            if i<j and D[i][j]:
                p={x:R[x] if min(x,R[x])==c else x for x in oc}
                need(C[oi[c]][oi[t]]!=C[oi[p[c]]][oi[p[t]]],'individual flip hostile')
                bad_edges.append([c,t,D[i][j]])
    poly=Counter({len(sc):1})
    for comp in components:
        out=Counter()
        for k,v in poly.items():out[k]+=v;out[k+2*len(comp)]+=v
        poly=out
    return {'n':row['n'],'ordinary':len(oc),'merged':len(mc),'self_converse':len(sc),'chiral_pairs':len(ch),
            'regimes':regimes,'ordinary_matrix_groups':matrix_groups,'signed_codes':ch,'signed_matrix':D,'signed_edges':bad_edges,
            'signed_spanning_forest_C_and_D':spanning,'signed_support_components':components,'physical_kernel_rank':len(components),
            'merged_kernel_rank':len(ch),'component_sizes':[len(c) for c in components],
            'physical_kernel_fixed_polynomial':dict(sorted(poly.items())),
            'component_generator_checks':component_checks,
            'ordinary_nodes':row['ordinary_nodes'],'source_data':row}


def producer_row(row):
    nodes=[{'code':x['code'],'converse':x['converse'],'H':x['H'],'aut':x['aut'],
            'SC':x['self_converse'],'degree':x['degree']} for x in row['ordinary_nodes']]
    reverse={x['code']:x['converse'] for x in nodes}
    merged=[{'code':x['code'],'converse':reverse[x['code']],'H':x['H'],'aut':x['aut'],
             'SC':x['self_converse'],'degree':x['degree']} for x in row['merged_nodes']]
    return {'n':row['n'],'ordinary_nodes':nodes,'merged_nodes':merged,
            'C':row['ordinary_A'],'blueC':row['ordinary_blue'],
            'A':row['merged_A'],'blueA':row['merged_blue']}


def selected_check(row):
    n=row['n'];pairs=list(combinations(range(n),2));loc={p:i for i,p in enumerate(pairs)}
    top=(1<<len(pairs))-1;ps=list(permutations(range(n)))
    def arc(c,u,v):return bool(c>>loc[tuple(sorted((u,v)))]&1)==(u<v)
    def transport(c,p):return sum(int(arc(c,p[u],p[v]))<<i for i,(u,v) in enumerate(pairs))
    cache={}
    def canon(c):
        if c not in cache:cache[c]=min(transport(c,p) for p in ps)
        return cache[c]
    chosen=sorted(row['ordinary_nodes'],key=lambda x:(x['H'],x['code']))[:6]
    C={(c,t):v for c,t,v in row['C']};BC={(c,t):v for c,t,v in row['blueC']}
    result=[]
    for node in chosen:
        c=node['code'];counts=Counter();blue=Counter();aut=0;hp=0
        for p in ps:
            aut+=transport(c,p)==c
            if not all(arc(c,p[i],p[i+1]) for i in range(n-1)):continue
            hp+=1
            pathmask=sum(1<<loc[tuple(sorted((p[i],p[i+1])))] for i in range(n-1))
            target=canon(c^top^pathmask);counts[target]+=1
            rev=[None]*n
            for u,v in zip(p,p[::-1]):rev[u]=v
            if transport(c,rev)==top^c:blue[target]+=1
        need(aut==node['aut'] and hp==node['H'],'selected independent metadata')
        need(canon(top^c)==node['converse'],'selected independent converse')
        need({t:v//aut for t,v in counts.items()}=={t:v for (a,t),v in C.items() if a==c},'selected independent C row')
        need({t:v//aut for t,v in blue.items()}=={t:v for (a,t),v in BC.items() if a==c},'selected independent blue row')
        result.append({'code':c,'H':hp,'aut':aut,'ordinary_C_row':sorted(counts.items())})
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-n',type=int,default=7);args=parser.parse_args()
    folder=Path(__file__).resolve().parent
    producer=folder.parent/'closure/data.json'
    frozen=producer.read_bytes();incoming={row['n']:producer_row(row) for row in json.loads(frozen)['rows']}
    rows=[];checks={}
    for n in range(3,args.max_n+1):
        if n<=6:
            row=direct(n)
            for key in ['ordinary_nodes','merged_nodes','C','blueC','A','blueA']:
                lhs=row[key];rhs=incoming[n][key]
                if key in ['C','blueC','A','blueA']:lhs=sorted(lhs);rhs=sorted(rhs)
                need(lhs==rhs,('complete independent producer comparison',n,key))
            checks[n]='complete independent ordinary orbit and physical-path census'
        else:
            row=incoming[n];checks[n]=selected_check(row)
        rows.append(signed_analysis(row))
    return {'status':'RESERVED proof candidate with FINITE-EXACT controls',
            'universe':'all ordinary/merged classes n3..max_n; independent full census through6; frozen producer plus selected exact rows at7',
            'native_operation':'complement all nonpath arcs; preserve the displayed Hamilton path',
            'producer_sha256':sha256(frozen).hexdigest(),'independent_checks':checks,'gauge_controls':gauge_controls(),
            'rows':rows,'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
