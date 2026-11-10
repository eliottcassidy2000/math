#!/usr/bin/env python3
"""Recount missing extension pairs from a separate audited strong inventory.

No producer imports; literal directed-word enumeration, not subset DP.
Completeness of the imported inventory remains an inherited audited input.
"""
from collections import Counter
from pathlib import Path
from hashlib import sha256
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=ROOT/'04-computation/experiments/h_spectrum_multiplicities_20261021/audit_S50/inventory_375.json'
TARGETS=((1,7),(3,11),(15,47),(35,107),(93,281),(17,55),(35,109),(41,127),(87,265),(13,27))

def require(x,msg):
    if not x: raise RuntimeError(msg)

def decode(n,b):
    a=[set() for _ in range(n)];k=0
    for i in range(n):
        for j in range(i+1,n):
            x,y=(i,j) if b>>k&1 else (j,i)
            a[x].add(y);k+=1
    return a

def count(a,remove=None):
    vertices=set(range(len(a))) - ({remove} if remove is not None else set())
    def visit(last,left):
        if not left:return 1
        return sum(visit(v,left-{v}) for v in (left if last is None else left & a[last]))
    return visit(None,vertices)

def strong(a):
    for v in range(len(a)):
        seen={v};todo=[v]
        while todo:
            x=todo.pop()
            for y in a[x]-seen: seen.add(y);todo.append(y)
        if len(seen)!=len(a):return False
    return True

def floor(n):
    if n%3==0:return 3*5**((n-3)//3)
    if n%3==1:return 5**((n-1)//3)
    return 9*5**((n-5)//3)

def main():
    d=json.loads(SOURCE.read_text());require(d['limit']==375,'input cap')
    needed={N for h,N in TARGETS}|{13,17,19,29,31,35,37,41,43}
    byN=Counter();deletions=Counter();actual={};controls={}
    for ns,bank in d['orders'].items():
        n=int(ns)
        for code,H,_,_ in bank:
            if H not in needed:continue
            a=decode(n,int(code))
            require(strong(a),'strong input')
            require(count(a)==H,'literal H disagrees with input')
            require(H>=floor(n),'Busch bound')
            byN[H]+=1
            for v in range(n):
                h=count(a,v);deletions[H]+=1
                actual[h,H]=actual.get((h,H),0)+1
                if (h,H) in ((3,13),(5,11),(5,13),(5,17),(5,19)) and (h,H) not in controls:
                    controls[h,H]={'order':n,'code':code,'root':v}
    for h,N in TARGETS:
        require(actual.get((h,N),0)==0,('false absence',h,N))
        require(N<floor(12),'target outside complete input bound')
    require(len(controls)==5,'positive controls')
    print(json.dumps({
        'status':'PASS','input_status':'inherited independent audited complete strong inventory H<=375',
        'count_method':'literal directed-word enumeration with no producer imports',
        'missing_pairs':[list(t) for t in TARGETS],
        'selected_strong_classes_by_child_count':dict(sorted(byN.items())),
        'class_recounts':sum(byN.values()),'vertex_deletion_recounts':sum(deletions.values()),
        'positive_controls':{f'{h},{N}':v for (h,N),v in controls.items()},
        'strong_order12_floor':floor(12),
        'source_sha256':{str(SOURCE.relative_to(ROOT)):sha256(SOURCE.read_bytes()).hexdigest()},
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()
    },indent=2,sort_keys=True))

if __name__=='__main__':main()
