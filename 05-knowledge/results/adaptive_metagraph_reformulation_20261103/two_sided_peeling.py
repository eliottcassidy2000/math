#!/usr/bin/env python3
"""Independent two-sided deletion critical-pair audit; exact small universe.
Uses only the separately authored literal ancestry auditor, no producer imports.
"""
from functools import lru_cache
from collections import Counter
from itertools import permutations
import importlib.util
from pathlib import Path
import json

# Resolve the independent local auditor explicitly, including under python -I.
_spec=importlib.util.spec_from_file_location('literal_ancestry_audit',Path(__file__).with_name('ancestry_audit.py'))
_audit=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_audit)
decode,strong,canonical,add,need,automorphisms=(getattr(_audit,name) for name in
    ('decode','strong','canonical','add','need','automorphisms'))


def induced(T,V):
    return tuple(sum(((T[u]>>v)&1)<<j for j,v in enumerate(V)) for u in V)


def moves(T):
    n=len(T);out=[]
    for x,row in enumerate(T):
        d=n-1-row.bit_count();o=row.bit_count()
        U=induced(T,tuple(v for v in range(n) if v!=x))
        if min(d,o)==1 and strong(U,tuple(range(n-1))):out.append((x,'in' if d==1 else 'out',U))
    return out


@lru_cache(None)
def terminals(T):
    children=moves(T)
    if not children:return frozenset(((len(T),canonical(T)),))
    return frozenset().union(*(terminals(U) for _,_,U in children))


def construction(K,a):
    """x indegree1 at a, b indegree1 at x, y outdegree1 at b."""
    k=len(K);T=add(add(K,a),k);n=len(T)
    # New y beats only b=k+1.
    return tuple(row|((1<<n) if i!=k+1 else 0) for i,row in enumerate(T))+(1<<(k+1),)


def named_branch(T,sequence):
    V=tuple(range(len(T)));steps=[]
    for x,role in sequence:
        rest=tuple(v for v in V if v!=x)
        d=sum(T[v]>>x&1 for v in rest);o=sum(T[x]>>v&1 for v in rest)
        need((d if role=='in' else o)==1 and strong(T,rest),'named legal deletion')
        steps.append({'delete':x,'role':role,'remaining':rest});V=rest
    U=induced(T,V)
    return {'steps':steps,'terminal_names':V,'terminal_rows':U,'terminal_code':canonical(U),
            'terminal':not moves(U)}


def main():
    universe=[]
    for n in range(3,6):
        total=0
        for code in range(1<<(n*(n-1)//2)):
            T=decode(n,code)
            if not strong(T,tuple(range(n))):continue
            total+=1;need(len(terminals(T))==1,'no smaller ordinary nonconfluence')
        universe.append({'n':n,'all_labelled_strong':total})
    T=(16,9,3,37,14,23)
    need(strong(T,tuple(range(6))),'hostile strong source')
    left=named_branch(T,((5,'in'),));right=named_branch(T,((0,'out'),(4,'in'),(5,'in')))
    need(left['terminal'] and right['terminal'] and len(left['terminal_names'])==5 and len(right['terminal_names'])==3,'different terminal orders')
    need(terminals(T)==frozenset(((3,2),(5,40))),'entire terminal ordinary set')
    examples=[];wrappers=[]
    for name,K in (('C3',(2,4,1)),('S4',(2,4,9,3)),('R5',(6,12,24,17,3))):
        k=len(K)
        def legal_word(A,p):return all(A[u]>>v&1 for u,v in zip(p,p[1:]))
        h=sum(legal_word(K,p) for p in permutations(range(k)))
        c2=sum(legal_word(K,p[:j]) and legal_word(K,p[j:]) for p in permutations(range(k)) for j in range(1,k))
        for a in range(len(K)):
            k=len(K);U=construction(K,a)
            need(strong(U,tuple(range(k+3))),'general family strong source')
            primitive=induced(U,tuple(i for i in range(k+3) if i!=k))
            need(strong(primitive,tuple(range(k+2))) and not moves(primitive),'cyclic K11 core irreducible')
            if a==0:
                wh=sum(legal_word(primitive,p) for p in permutations(range(k+2)))
                wa=len(automorphisms(primitive));ka=len(automorphisms(K))
                need(wh==3*h+c2 and wa==ka,'wrapper cover transfer and exact automorphism group')
                need(sum(k+1-r.bit_count()==1 for r in primitive)==1 and sum(r.bit_count()==1 for r in primitive)==1,'characteristic wrapper ports')
                wrappers.append({'K':name,'H_K':h,'C2_nonempty':c2,'H_wrapper':wh,'Aut_K_and_wrapper':ka})
            seq=named_branch(U,((k+2,'out'),(k+1,'in'),(k,'in')))
            need(tuple(seq['terminal_names'])==tuple(range(k)) and tuple(seq['terminal_rows'])==K,'return exactly named K')
            examples.append({'K':name,'root':a,'source_order':k+3,'irreducible_branch_order':k+2,'other_branch_reaches':k})
    return {'status':'REFUTED two-sided confluence; PROVED family at every order n>=6',
            'minimality_universe':universe,'smallest_hostile':{'rows':T,'branches':[left,right]},
            'literal_family_controls':examples,
            'wrapper_controls':wrappers,
            'scope':'Global all-order nonconfluence follows the written construction; finite checks validate named critical pairs, not every larger tournament.'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
