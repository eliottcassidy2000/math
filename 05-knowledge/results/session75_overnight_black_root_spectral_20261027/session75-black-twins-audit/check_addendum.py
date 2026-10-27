#!/usr/bin/env python3
"""Bounded independent audit: exactly two n10 matrices, no new census.
Cycle collections are enumerated as directed permutation covers with fixed
points unused and all nontrivial cycles odd. No producer imports.
"""
from collections import Counter
import json,pathlib,itertools,hashlib
from audit_checks import *
H=pathlib.Path(__file__).resolve().parent
OLD=(966,1004,1016,1009,995,193,384,768,544,96)
NEW=(1008,1001,995,997,494,448,640,768,576,48)

def permutation_profile(a):
    n=len(a);profile=Counter();cycles_seen=set();accepted_permutations=[];all_covers=0;lengths=(3,5,7,9)
    permitted=[[j for j in range(n) if j==i or a[i][j]] for i in range(n)]
    def finish(p):
        nonlocal all_covers
        all_covers+=1;visited=set();cycles=[];multiplicities=[0]*4
        for v in range(n):
            if v in visited:continue
            word=[];u=v
            while u not in visited:
                visited.add(u);word.append(u);u=p[u]
            need(u==v,'permutation cycle decomposition failed')
            if len(word)==1:continue
            if len(word)%2==0:return
            need(word[0]==min(word),'cycle root is not least label')
            multiplicities[lengths.index(len(word))]+=1;cycles.append(tuple(word))
        profile[tuple(multiplicities)]+=1;cycles_seen.update(cycles);accepted_permutations.append(tuple(p))
    # Enumerate all bijections satisfying arcs or fixed points; no cycle-word
    # or disjoint-cycle packing generator is used.
    def match(row,available,p):
        if row==n:finish(p);return
        for v in permitted[row]:
            if available>>v&1:match(row+1,available^(1<<v),p+[v])
    match(0,(1<<n)-1,[])
    q=Counter()
    for mult,c in profile.items():q[sum(mult)]+=c
    return profile,q,cycles_seen,all_covers,accepted_permutations

def incoming_recipe():
    # Read pinned stack/C3sub/doubled definitions, then reconstruct only M=C3.
    D=[[False]*5 for _ in range(5)]
    for u,v in [(0,1),(0,2),(0,3),(1,2),(2,3),(3,1),(1,4),(2,4),(3,4),(4,0)]:D[u][v]=True
    T=[[False]*10 for _ in range(10)]
    for u in range(5):
        for v in range(5):T[u][v]=D[u][v];T[u+5][v+5]=D[v][u];T[u][v+5]=True
    T[0][5]=False;T[5][0]=True
    return tuple(tuple(r) for r in T)

def main():
    producer=json.loads((H/'producer/normal.json').read_text());results=[];all_profiles=[];triangles=[]
    expected_profiles=[{(0,0,0,0):1,(1,0,0,0):14,(0,1,0,0):30,(0,0,1,0):52,(0,0,0,1):30,(2,0,0,0):35,(1,1,0,0):48,(1,0,1,0):6,(0,2,0,0):4,(3,0,0,0):4},{(0,0,0,0):1,(1,0,0,0):14,(0,1,0,0):54,(0,0,1,0):90,(0,0,0,1):36,(2,0,0,0):22,(1,1,0,0):36,(0,2,0,0):9}]
    for name,r,expected in zip(('old','new'),(OLD,NEW),expected_profiles):
        a=from_rows(r);F,Q,cycles,cover_count,valid=permutation_profile(a);need(dict(F)==expected,'cycle polynomial differs')
        src=producer[name];need(list(r)==src['rows'],'producer named source mismatch')
        reported={tuple(x['multiplicities']):x['count'] for x in src['odd_cycle_collection_polynomial']['collection_coefficients']}
        need(dict(F)==reported,'producer profile mismatch');need(cycles==set(map(tuple,src['independently_enumerated_directed_odd_cycle_words'])),'producer literal cycle list incomplete or wrong')
        pp=list(paths(a));blue=[p for p in pp if anti_failure(a,p) is None];black=[p for p in pp if anti_failure(a,p) is not None];exits=[p for p in black if not strong(reverse(a,p))]
        need((len(pp),len(blue),len(black),len(exits))==((657,33,624,6) if name=='old' else (657,99,558,0)),'Hamilton/color/exit mismatch')
        at2=sum(count*2**sum(mult) for mult,count in F.items());need(at2==len(pp),'F at2 differs from direct Hamilton count')
        terminal,moves=replay(a,[m['R_path'] for m in src['core_route']],require_core=True,detail=True)
        for m,p in zip(moves,src['core_route']):need(m['K_source_rows']==p['K_source_rows'] and m['K_target_rows']==p['K_target_rows'],'core native edge mismatch')
        need(cut_certificate(a)[:3]==(1,2,[992]),'exact invariants mismatch')
        tri={tuple(sorted(c)) for c in cycles if len(c)==3};triangles.append(tri);all_profiles.append(Q)
        results.append({'name':name,'rows':r,'H':len(pp),'blue':len(blue),'black':len(black),'black_one_step_exits':len(exits),'all_directed_permutation_covers_including_even_cycles':cover_count,'odd_cycle_collections':len(valid),'F':[{'mult':k,'coefficient':v} for k,v in sorted(F.items())],'Q':dict(sorted(Q.items())),'F_at_all_variables_2':at2,'distinct_odd_cycle_words':len(cycles),'odd_collection_permutations_sha256':hashlib.sha256(json.dumps(valid,separators=(',',':')).encode()).hexdigest(),'named_core_ports':core(terminal),'core_moves':moves})
    diff={k:all_profiles[1][k]-all_profiles[0][k] for k in range(4)};need(diff=={0:0,1:68,2:-26,3:-4},'Q difference factor mismatch')
    a,b=from_rows(OLD),from_rows(NEW);need(degrees(a)==degrees(b),'named score mismatch')
    values=[]
    for mask in range(1,1023):
        S={v for v in range(10) if mask>>v&1};ca,cb=cut(a,S),cut(b,S);need(ca==cb,'named cut mismatch');values.append(ca)
    need(len(triangles[0])==len(triangles[1])==14 and triangles[0]!=triangles[1],'triangle scalar/set distinction failed')
    inc=incoming_recipe();phi=(4,1,3,2,0,9,6,8,7,5);need(sorted(phi)==list(range(10)) and all(b[u][v]==inc[phi[u]][phi[v]] for u in range(10) for v in range(10)),'incoming isomorphism failed')
    need(rows(inc)==producer['incoming_THM4811']['rows_from_source_recipe'],'incoming recipe named rows differ');need(list(phi)==producer['incoming_THM4811']['isomorphism_new_label_to_incoming_label'],'incoming named permutation differs')
    print(json.dumps({'status':'PASS','scope':'two fixed order10 sources and one pinned incoming E10 recipe only','cycle_convention':'unordered directed odd cycles modulo rotation but not reversal; unused vertices are fixed points','independent_method':'all allowed directed permutation covers, reject even nontrivial cycles; no producer imports or cycle-packing generator','sources':results,'Q_E_minus_Q_old_coefficients':diff,'Q_difference_factorization':'2*t*(2-t)*(2*t+17)','all_1022_named_cuts_equal':True,'labelled_triangle_sets_equal':False,'scalar_c3_each':14,'only_old_cyclic_triples':sorted(triangles[0]-triangles[1]),'only_E10_cyclic_triples':sorted(triangles[1]-triangles[0]),'incoming':{'rows':rows(inc),'source_commit':'87d61f67673a6b10cbafff1bf9e1377f3c41e56c','E10_to_incoming_permutation':phi,'all_100_ordered_pair_equations_checked':True},'not_certified':['incoming full n10 census','general determination of black distance by F','any new search range']},sort_keys=True,indent=2))
if __name__=='__main__':main()
