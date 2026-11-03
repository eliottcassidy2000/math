#!/usr/bin/env python3
"""Independent literal audit of triangle two-birth determinant and rank.
Uses sets of directed arcs, all vertex permutations, and rational elimination.
No producer/repository code or external dependency is imported.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product, combinations
import json


def need(ok,why):
    if not ok: raise RuntimeError(why)


def h(vertices,E):
    return sum(all((a,b) in E for a,b in zip(p,p[1:])) for p in permutations(vertices))


def old(n,mask):
    return {(j,i) if mask>>k&1 else (i,j) for k,(i,j) in enumerate(combinations(range(n),2))}


def add(n,E,mask):
    return E|{(i,n) if mask>>i&1 else (n,i) for i in range(n)}


def canon(n,E):
    es=tuple(combinations(range(n),2))
    return min(sum((1<<k) for k,(i,j) in enumerate(es) if (p[i],p[j]) in E) for p in permutations(range(n)))


def det(A):
    M=[list(map(Fraction,row)) for row in A];sign=1;d=Fraction(1)
    for k in range(len(M)):
        p=next((i for i in range(k,len(M)) if M[i][k]),None)
        if p is None:return 0
        if p!=k:M[p],M[k]=M[k],M[p];sign=-sign
        pivot=M[k][k];d*=pivot
        for i in range(k+1,len(M)):
            c=M[i][k]/pivot
            for j in range(k,len(M)):M[i][j]-=c*M[k][j]
    return int(sign*d)


def rank(A):
    M=[list(map(Fraction,row)) for row in A];r=0
    for j in range(len(M[0])):
        p=next((i for i in range(r,len(M)) if M[i][j]),None)
        if p is None:continue
        M[p],M[r]=M[r],M[p];pivot=M[r][j]
        for i in range(r+1,len(M)):
            c=M[i][j]/pivot
            for k in range(j,len(M[0])):M[i][k]-=c*M[r][k]
        r+=1
        if r==len(M):break
    return r


Tplus={(0,1),(1,2),(2,0)}
Tminus={(1,0),(2,1),(0,2)}

def children(a,b,c):
    am=sum(x<<i for i,x in enumerate(a))
    bm=sum(x<<i for i,x in enumerate(b))|(c<<3)
    return add(4,add(3,Tplus,am),bm),add(4,add(3,Tminus,am),bm)


def main():
    histogram=Counter();partition=Counter();canonical_rows=[]
    for a in product((0,1),repeat=3):
        for b in product((0,1),repeat=3):
            D=det([[1,a[i],b[i]] for i in range(3)])
            for c in (0,1):
                U,V=children(a,b,c)
                hu,hv=h(range(5),U),h(range(5),V)
                need(hu-hv==2*D*(1-2*c+sum(b)-sum(a)),'detector formula')
                repeated=len(set(zip(a,b)))<3
                iso=canon(5,U)==canon(5,V)
                need(iso==repeated,'ordinary unpointed class iff')
                histogram[hu-hv]+=1
                key='same_class' if iso else 'different_same_H' if hu==hv else 'different_H'
                partition[key]+=1
                if not repeated and hu==hv:need(hu==9,'equal scalar distinct class')
    expected_hist={0:92,2:12,-2:12,4:6,-4:6}
    need(dict(histogram)==expected_hist,'difference histogram')
    corners=tuple(product((0,1),repeat=2))
    for missing in corners:
        signatures=[x for x in corners if x!=missing]
        a,b=zip(*signatures)
        for c in (0,1):
            U,V=children(a,b,c)
            hu,hv=h(range(5),U),h(range(5),V)
            decks=[sorted(h(tuple(x for x in range(5) if x!=v),G) for v in range(5)) for G in (U,V)]
            if hu==hv:need(sorted(decks)==[[1,1,5,5,5],[1,3,3,5,5]],'equal-H deck certificate')
            canonical_rows.append({'missing':missing,'c':c,'H':[hu,hv],
                                   'Delta':det([[1,a[i],b[i]] for i in range(3)]),'decks':decks})
    contexts=((0,0),(0,1),(0,2),(0,3),(0,9),(0,10),(1,2))
    matrix=[[h(range(5),add(4,add(3,old(3,m),a),b)) for m in range(8)] for a,b in contexts]
    expected_matrix=[[1,1,3,1,1,3,1,1],[3,5,9,9,3,9,5,9],[5,3,9,3,9,9,9,5],
                     [5,5,13,9,9,13,11,11],[1,3,5,5,1,5,3,5],[3,1,5,1,5,5,5,3],
                     [5,5,13,9,9,11,15,13]]
    need(matrix==expected_matrix,'seven context matrix')
    minor=det([row[:7] for row in matrix])
    need(minor==-128 and rank(matrix)==7,'seven query certificate')
    need(all(sum(x*(-1)**m.bit_count() for m,x in enumerate(row))==0 for row in matrix),'cubic Walsh kernel')
    one_birth=[[h(range(4),add(3,old(3,m),cut)) for m in range(8)] for cut in range(8)]
    need(rank(one_birth)==4,'one-birth rank four')
    # Named first children are ordinarily isomorphic, but the same next cut separates them.
    Xp=add(3,Tplus,1);Xm=add(3,Tminus,1)
    need(h(range(4),Xp)==h(range(4),Xm)==5 and canon(4,Xp)==canon(4,Xm),'same first child class')
    Up=add(4,Xp,2);Um=add(4,Xm,2)
    rows=lambda G:tuple(sum(1<<v for u,v in G if u==i) for i in range(5))
    need(rows(Up)==(10,20,1,6,13) and rows(Um)==(12,17,2,6,13),'named final rows')
    need((h(range(5),Up),h(range(5),Um))==(13,11),'prime-target hostile')
    two_birth= [[h(range(5),add(4,add(3,old(3,m),a),b)) for m in range(8)]
                for a in range(8) for b in range(16)]
    need(rank(two_birth)==7,'all two-birth query rank')
    # At order two the complete one-birth response determines the named orientation.
    small=[tuple(h(range(3),add(2,old(2,m),a)) for a in range(4)) for m in range(2)]
    need(small[0]!=small[1],'minimum collision order')
    return {'status':'INDEPENDENT PASS, literal full vertex permutations; no producer imports',
            'canonical_eight_rows':canonical_rows,'all_128_partition':dict(partition),
            'difference_histogram':dict(sorted(histogram.items())),
            'seven_query_matrix':matrix,'first_seven_columns_determinant':minor,
            'one_birth_rank':4,'all_two_birth_rank':7,
            'first_children_H':[5,5],'final_named_H':[13,11],
            'order_two_one_birth_profiles':small}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
