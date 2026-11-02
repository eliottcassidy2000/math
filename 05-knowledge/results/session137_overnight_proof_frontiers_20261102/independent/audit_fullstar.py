"""Independent all-class order7 full-star graph and exact H147 barrier."""
from tournament_kernel import *
from collections import defaultdict
from math import gcd
import json

def run():
    reps=list(class_bank(7))[-1][1];bank={c:hp(matrix(7,c)) for c in reps if strong(matrix(7,c))}
    require(len(bank)==353,'strong7 classes')
    parent={c:c for c in bank};cards=defaultdict(set)
    def find(c):
        while parent[c]!=c:parent[c]=parent[parent[c]];c=parent[c]
        return c
    for c in bank:
        a=matrix(7,c)
        for x in range(7):cards[canon(induced(a,[i for i in range(7) if i!=x]))].add(c)
    for members in cards.values():
        for c,d in combinations(members,2):
            if abs(bank[c]-bank[d])<=2:
                p,q=find(c),find(d)
                if p!=q:parent[max(p,q)]=min(p,q)
    groups=defaultdict(list)
    for c in bank:groups[find(c)].append(c)
    require(len(groups)==11,'full-star component count')
    rows=[6,112,74,35,13,21,57];a=tuple(tuple(rows[i]>>j&1 for j in range(7)) for i in range(7));cc=canon(a)
    require(hp(a)==147 and strong(a),'H147 source')
    require(sum(h==147 for h in bank.values())==1,'H147 class uniqueness')
    require(len(groups[find(cc)])==1,'H147 isolated full-star class')
    automorphisms=sum(induced(a,p)==a for p in permutations(range(7)))
    require(automorphisms==1,'H147 rigidity')
    images=set();equal=0;deck_H=[];deck_gcd=[]
    for x in range(7):
        old=induced(a,[i for i in range(7) if i!=x])
        require(strong(old),'all H147 deletion cards strong')
        old_H=hp(old);deck_H.append(old_H);response_gcd=0
        for cut in range(64):
            b=append(old,cut)
            h=hp(b);response_gcd=gcd(response_gcd,h-old_H)
            if not strong(b):continue
            images.add(h)
            require(h not in (145,149),'star target absence')
            if h==147:require(canon(b)==cc,'equal-H star isomorphism');equal+=1
        deck_gcd.append(response_gcd)
    require(deck_H==[43,29,29,41,43,43,43] and deck_gcd==[2]*7,'deletion sidecars')
    code=encoding(a);edges=list(combinations(range(7),2));repairs=[];distance_counts={}
    for distance in (1,2,3):
        found=[]
        for ks in combinations(range(21),distance):
            childcode=code
            for k in ks:childcode^=1<<k
            b=matrix(7,childcode)
            if strong(b) and hp(b)==145:found.append({'code':childcode,'flips':[edges[k] for k in ks]})
        distance_counts[distance]=len(found)
        if distance==3:repairs=found
    require(distance_counts=={1:0,2:0,3:4},'minimal Hamming repair')
    outputgroups=sorted([{'classes':len(cs),'H':sorted({bank[c] for c in cs})} for cs in groups.values()],key=lambda r:r['H'])
    print(json.dumps({'scope':'FINITE-EXACT all353 strong order7 classes and arbitrary six-vertex deletion cards; all448 star presentations of explicit H147 source',
                      'components':outputgroups,'H147_code':code,'strong_star_H':sorted(images),
                      'H147_automorphisms':automorphisms,'deletion_H':deck_H,'deletion_response_gcd':deck_gcd,
                      'H147_equal_star_presentations':equal,'H145_Hamming_distance_counts':distance_counts,
                      'distance3_repairs':repairs},sort_keys=True,indent=2))
if __name__=='__main__':run()
