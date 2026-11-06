#!/usr/bin/env python3
"""Independent sidecar audit: exact rank and residue-state tests."""
from fractions import Fraction
from collections import deque
import json
from independent_check import ROOT, bank, bound, D, U, raw, raw_steps, need, word, target, section

def rank(rows):
    a=[[Fraction(x) for x in row] for row in rows]
    pivot=0
    for col in range(len(a[0])):
        found=next((i for i in range(pivot,len(a)) if a[i][col]),None)
        if found is None:continue
        a[pivot],a[found]=a[found],a[pivot]
        c=a[pivot][col];a[pivot]=[x/c for x in a[pivot]]
        for i in range(len(a)):
            if i!=pivot:
                c=a[i][col];a[i]=[x-c*y for x,y in zip(a[i],a[pivot])]
        pivot+=1
    return pivot

def collision_matrices():
    basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
    cases=[('O=J',0,2,[(1,1,1,1)]),('Z=J',1,2,[(1,0,2,1)]),
           ('O=Z',0,1,[(1,0,0,0),(0,1,1,0)])]
    out={}
    for label,i,j,kernel_basis in cases:
        columns=[tuple(x-y for x,y in zip(raw_steps(v)[i],raw_steps(v)[j])) for v in basis]
        rows=[tuple(columns[c][r] for c in range(4)) for r in range(4)]
        r=rank(rows)
        need(r+len(kernel_basis)==4,('nullity',label))
        need(rank(kernel_basis)==len(kernel_basis),('kernel independence',label))
        for v in kernel_basis:
            need(all(sum(x*y for x,y in zip(row,v))==0 for row in rows),('kernel containment',label))
        out[label]={'difference_rank':r,'kernel_basis':kernel_basis}
    return out

def residue_closure():
    states=[(1,0,0,1),(1,1,1,1)]
    need(tuple(x%2 for x in (-1,-7,-1,-7)) in states,'initial parity')
    for s in states:
        for q in raw_steps(s):
            need(tuple(x%2 for x in q) in states,('residue closure',s,q))
    return {'invariant_residue_states':states,'transitions_checked':6,'conclusion':'All integer descendants have odd d, so cannot enter the O=Z kernel.'}

def first_gap():
    seeds=bank(ROOT,128)
    count=0
    for s in seeds:
        Ds,E,Us,d=raw(s);F=Ds-E+Us
        need(F>Ds,'strict growth')
        for k in range(11):
            q=word(s,'O'*k)
            for letter in 'ZJ':
                need(D(word(q,letter))==(F if k==0 else Ds+(2**k)*Us),('first growth formula',s,k))
        values=sorted({D(q) for q in bank(s,8*Us)})
        need(values[:2]==[Ds,F],('first two D values',s,values[:2]))
        if s[1]>0:
            values=sorted({target(q) for q in bank(s,bound(s,F)) if section(q)})
            need(values[:2]==[Ds,F],('first two section targets',s,values[:2]))
            count+=1
    return {'source_height_cap':128,'sources':len(seeds),'first_growth_O_lengths':[0,10],
            'positive_f_first_two_section_checks':count,'D_tail_bank_height':'8 U_s'}

def hostiles():
    for r in range(51):
        lam=3**r
        s=word(ROOT,'Z'+'OJ'*r+'OOZ')
        need(raw(s)==tuple(lam*x for x in (15,1,23,1)),('hostile scaling',r))
        Ds,E,Us,d=raw(s)
        need(s[1]==8*lam and Ds<target(s)<Ds-E+Us,('hostile gap',r))
        v=word(ROOT,'OOOOJ'+'OJ'*r)
        need(D(v)==U(v)==33*lam,('represented balanced target',r))
        q=word(v,'OJJJJ')
        need(section(q) and target(q)==target(s),('section representation',r))
        if r<3:
            vals={target(q) for q in bank(s,bound(s,target(s))) if section(q)}
            need(target(s) not in vals,('finite no-repair check',r))
    v=word(ROOT,'OOOOZ')
    need(raw(v)==(33,1,63,1),'separate initial representative')
    return {'scaled_formula_r':[0,50],'complete_target_obstruction_search_r':[0,2]}

print(json.dumps({'status':'PASS','scope':'Separate first-growth and formal parity sidecars only; no target coverage.',
                  'exact_collision_nullspaces':collision_matrices(),'all_order_parity':residue_closure(),
                  'first_growth':first_gap(),'hostile_family':hostiles()},indent=2,sort_keys=True))
