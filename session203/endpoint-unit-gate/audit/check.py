#!/usr/bin/env python3
"""Independent literal DFS and permutation/cut checks for unit-root curvature."""
from itertools import product, permutations, combinations
import json

def need(b,s):
    if not b:raise RuntimeError(s)

def tournaments(n):
    pairs=list(combinations(range(n),2))
    for code in range(1<<len(pairs)):
        a=[0]*n
        for b,(i,j) in enumerate(pairs):
            if code&(1<<b):a[i]|=1<<j
            else:a[j]|=1<<i
        yield a

def valid(a):
    n=len(a)
    for i,row in enumerate(a):
        need(0<=row<1<<n and not row&(1<<i),'row')
        for j in range(i):need(((row>>j)&1)+((a[j]>>i)&1)==1,'pair')

def literal(a):
    full=(1<<len(a))-1; h=s=c=calls=0
    def dfs(last,mask,start):
        nonlocal h,s,c,calls
        calls+=1
        if mask==full:
            h+=1
            if start==0:
                s+=1;c+=bool(a[last]&1)
            return
        todo=a[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            dfs(bit.bit_length()-1,mask|bit,start)
    for v in range(len(a)):dfs(v,1<<v,v)
    return h,s,c,calls

def input_graph(a,b,reverse_twins=False):
    # Labels v,x1,x2,x3,x4; construct the private chain first.
    rows=[0]*5
    for i in range(1,5):
        for j in range(i+1,5):rows[i]|=1<<j
    for v,bit in enumerate((1,a,b,0),1):
        if bit:rows[0]|=1<<v
        else:rows[v]|=1
    if reverse_twins:rows[2]^=8;rows[3]^=4
    return rows

def context(K,mask,U,complement=False):
    n=len(K); size=n+len(U)-1; a=[0]*size
    for i,row in enumerate(K):a[i]=row
    image=[0]+list(range(n,size))
    for i,j in combinations(range(len(U)),2):
        forward=bool(U[i]&(1<<j))
        if i==0 and complement:forward=not forward
        x,y=image[i],image[j]
        if forward:a[x]|=1<<y
        else:a[y]|=1<<x
    for x in range(n,size):
        for y in range(1,n):
            if mask&(1<<y):a[x]|=1<<y
            else:a[y]|=1<<x
    valid(a);return a

def induced(a,keep):
    return [sum(1<<j for j,v in enumerate(keep) if a[w]&(1<<v)) for w in keep]

def rooted_profile_by_permutations(U):
    h=s=L=0;n=len(U)
    def ispath(p):return all(U[x]&(1<<y) for x,y in zip(p,p[1:]))
    for word in permutations(range(n)):
        if ispath(word):h+=1;s+=word[0]==0
        # The first path must end at the root, fixing the only possible cut.
        cut=word.index(0)+1
        if ispath(word[:cut]) and ispath(word[cut:]):L+=1
    return h,s,L

def main():
    sources=[input_graph(a,b) for a,b in product((0,1),repeat=2)]
    source_data=[];calls=0
    for U in sources:
        valid(U);h,s,c,w=literal(U);calls+=w
        need(s==1 and c>0,'physical strong unit-start source')
        source_data.append({'rows':U,'H':h,'start':s,'cycles':c})
    need([x['H'] for x in source_data]==[9,9,11,9],'source count')
    contexts=curvature=reverse_checks=literal_graphs=0
    quotient_values=set();delta_mod4=set();zero_root_controls=0
    for n in range(1,5):
        for K in tournaments(n):
            valid(K)
            for mask0 in range(1<<(n-1)):
                contexts+=1;mask=mask0<<1
                for complement in (False,True):
                    for reverse in (False,True) if n<=3 else (False,):
                        four=[context(K,mask,input_graph(a,b,reverse),complement)
                              for a,b in product((0,1),repeat=2)]
                        counts=[]
                        for T in four:
                            h,_,_,w=literal(T);counts.append(h);calls+=w;literal_graphs+=1
                        x,y=n+1,n+2
                        T=four[0]
                        for z in range(len(T)):
                            if z not in (0,x,y):need(bool(T[x]&(1<<z))==bool(T[y]&(1<<z)),'not twins')
                        Q=induced(T,[i for i in range(len(T)) if i not in (0,y)])
                        valid(Q);qh,_,_,w=literal(Q);calls+=w;literal_graphs+=1
                        delta=counts[3]-counts[2]-counts[1]+counts[0]
                        need(qh%2==1 and delta==-2*qh,'twin curvature')
                        quotient_values.add(qh);delta_mod4.add(delta%4)
                        if reverse:reverse_checks+=1
                        else:curvature+=1
                # Discarding the root erases both free bits; the same obstruction does not apply.
                deleted=[induced(U,[1,2,3,4]) for U in sources]
                need(all(a==deleted[0] for a in deleted),'zero-root boundary')
                zero_root_controls+=1
    need(contexts==549 and curvature==1098 and delta_mod4=={2},'context census')
    positive=0;positive_samples=[]
    for n in range(1,6):
        for U in tournaments(n):
            h,s,L=rooted_profile_by_permutations(U)
            P=context([2,0],0,U)
            C=context([4,1,2],0,U)
            ph,_,_,wp=literal(P);ch,_,_,wc=literal(C)
            calls+=wp+wc;literal_graphs+=2;positive+=1
            need(ph==h+L-s and ch==h+2*L and ch-2*ph==2*s-h,'positive identity')
            if n==1 or n==3 and h==3:positive_samples.append({'U':U,'H':h,'start':s,'L':L,'parent_H':ph,'child_H':ch})
    need(positive==1099,'positive census')
    print(json.dumps({'status':'PASS','source_data':source_data,'contexts':contexts,
                      'curvature_checks':curvature,'extra_reversed_twin_checks':reverse_checks,
                      'quotient_H_values':sorted(quotient_values),'curvature_mod4':sorted(delta_mod4),
                      'zero_copied_root_scope_controls':zero_root_controls,
                      'positive_root_zero_inputs':positive,'positive_samples':positive_samples,
                      'literal_graphs':literal_graphs,'literal_DFS_calls':calls},indent=2,sort_keys=True))

if __name__=='__main__':main()
