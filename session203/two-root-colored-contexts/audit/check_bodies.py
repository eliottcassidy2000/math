#!/usr/bin/env python3
"""Literal colored cut capacities and direct labelled-bin convolution checks."""
from itertools import permutations,combinations,product
from collections import Counter
import random,json
from check_tensors import need,valid,CASES,ZERO,H,W,FACTORS

VISITS=0
def hp(a):
    global VISITS
    if not a:return 1
    full=(1<<len(a))-1;total=0
    def visit(last,mask):
        global VISITS
        nonlocal total
        VISITS+=1
        if mask==full:total+=1;return
        todo=a[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            visit(bit.bit_length()-1,mask|bit)
    for v in range(len(a)):visit(v,1<<v)
    return total

def tournament(n,code):
    a=[0]*n
    for bit,(i,j) in enumerate(combinations(range(n),2)):
        if code&(1<<bit):a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def color_functions(a,col):
    n=len(a);f=[[0]*(1<<n) for _ in range(4)]
    def visit(last,mask,first):
        f[2*col[first]+col[last]][mask]+=1
        todo=a[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            visit(bit.bit_length()-1,mask|bit,first)
    for v in range(n):visit(v,1<<v,v)
    return f

def evaluate(poly,f,n):
    total=0
    for powers,coefficient in poly.items():
        slots=[kind for kind,count in enumerate(powers) for _ in range(count)]
        if not slots:
            total+=coefficient*int(n==0);continue
        value=0
        for labels in product(range(len(slots)),repeat=n):
            masks=[0]*len(slots)
            for v,label in enumerate(labels):masks[label]|=1<<v
            term=1
            for kind,mask in zip(slots,masks):term*=f[kind][mask]
            value+=term
        total+=coefficient*value
    return total

def capacities(a):
    n=len(a);Q=[[0]*n for _ in range(n)]
    def path(p):return all(a[x]&(1<<y) for x,y in zip(p,p[1:]))
    for word in permutations(range(n)):
        for cut in range(1,n):
            if path(word[:cut]) and path(word[cut:]):Q[word[cut-1]][word[cut]]+=1
    for i,j in combinations(range(n),2):need(Q[i][j]+Q[j][i]>0,'cut capacity not positive')
    return Q

def cut_check(a,col,Q=None):
    n=len(a);Q=capacities(a) if Q is None else Q
    apex=a[:]+[0]
    for v,color in enumerate(col):
        if color==0:apex[n]|=1<<v
        else:apex[v]|=1<<n
    valid(apex);hb=hp(a);defect=hp(apex)-hb
    weighted=sum(Q[i][j]+Q[j][i] for i in range(n) for j in range(n)
                 if a[i]&(1<<j) and col[i]==0 and col[j]==1)
    bad=any(a[i]&(1<<j) and col[i]==0 and col[j]==1 for i in range(n) for j in range(n))
    need(defect==weighted and defect>=0 and (defect>0)==bad,'cut formula/equality')
    f=color_functions(a,col)
    need(evaluate(W,f,n)==defect and evaluate(H,f,n)==hb,'colored cut polynomial')
    return hb,defect,bad,f

def context(core,masks,body,col):
    n=len(core);a=core[:]+[row<<n for row in body]
    for j,color in enumerate(col):
        for c in range(n):
            if masks[color]&(1<<c):a[n+j]|=1<<c
            else:a[c]|=1<<(n+j)
    valid(a);return a

def canon(a):
    return min(tuple(sum(1<<j for j,v in enumerate(p) if a[u]&(1<<v)) for u in p)
               for p in permutations(range(len(a))))

def main():
    cut_cases=0;representatives={};capacity_pairs=0
    for n in range(5):
        reps=set()
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code);valid(a);Q=capacities(a);capacity_pairs+=n*(n-1)//2
            reps.add(canon(a))
            for col in product((0,1),repeat=n):cut_check(a,col,Q);cut_cases+=1
        representatives[n]=sorted(reps)
    need(cut_cases==1099,'colored cut universe')
    inputs=[(list(a),col) for n,reps in representatives.items() for a in reps for col in product((0,1),repeat=n)]
    need(len(inputs)==87,'complete representative-colored universe')
    rng=random.Random(203221)
    for _ in range(12):inputs.append((tournament(5,rng.randrange(1<<10)),tuple(rng.randrange(2) for _ in range(5))))
    block_body=[58,60,57,16,32,8] # C3 -> C3; each color class can be nontransitive.
    inputs.append((block_body,(1,1,1,0,0,0)))
    inputs.append((block_body,(0,0,0,1,1,1)))
    inputs.append(([2,0],(0,1)))
    pairs=strict=0;hostile=[];block_controls=[];minimum_positive={}
    for index,(body,col) in enumerate(inputs):
        hb,defect,bad,f=cut_check(body,col)
        for name,p,pm,c,cm,flip in CASES:
            parent=context(p,pm,body,col);child=context(c,cm,body,col)
            oldmap=list(range(len(p)))+list(range(len(c),len(child)))
            changed=[(i,j) for i,j in combinations(range(len(parent)),2)
                     if bool(parent[i]&(1<<j))!=bool(child[oldmap[i]]&(1<<oldmap[j]))]
            need(changed==[flip],'physical old-edge edits')
            ph=hp(parent);ch=hp(child);excess=ch-2*ph-hb
            need(excess==evaluate(FACTORS[name],f,len(body)),'literal context versus convolution')
            need(excess>=0 and (excess>0)==bad,'context equality iff ordered colors')
            need(excess>=(4 if name=='A' else 21)*defect,'positive cut lower bound')
            if excess>0:minimum_positive[name]=min(excess,minimum_positive.get(name,excess))
            pairs+=1;strict+=excess>0
            if index==len(inputs)-1:
                expected=(133,279,12) if name=='A' else(115,295,64)
                need((ph,ch,excess)==expected,'specified hostile')
                hostile.append({'family':name,'parent_rows':parent,'child_rows':child,'parent_H':ph,'child_H':ch,'body_H':hb,'excess':excess})
            if index in (len(inputs)-3,len(inputs)-2):
                need(hb==9,'two-cycle body count')
                if index==len(inputs)-3:
                    expected=(1413,2835) if name=='A' else (765,1539)
                    need((ph,ch)==expected and excess==0,'nontransitive ordered blocks')
                block_controls.append({'family':name,'colors':col,'body_H':hb,'parent_H':ph,'child_H':ch,'cut_defect':defect,'excess':excess})
    need(pairs==204,'context pair universe')
    print(json.dumps({'status':'PASS','all_labelled_colored_cut_cases_through4':cut_cases,
                      'capacity_pairs_checked':capacity_pairs,'uncolored_class_counts':{n:len(r) for n,r in representatives.items()},
                      'complete_representative_colorings_through4':87,'additional_seeded_order5_cases':12,
                      'additional_nontransitive_block_order6_cases':2,'nontransitive_block_controls':block_controls,
                      'explicit_unsorted_spine_case':1,'context_pairs':pairs,'strict_cases':strict,
                      'minimum_positive_excess_observed':minimum_positive,'literal_DFS_visits':VISITS,'hostiles':hostile,
                      'scope':'Only the two specified physical repair pairs; empty body and arbitrary internal body orientations included.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
