#!/usr/bin/env python3
"""Independent literal DFS audit of all frozen seed and profile graphs."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json

HERE=Path(__file__).resolve().parent

def need(b,s):
    if not b:raise RuntimeError(s)

def valid(a):
    n=len(a);need(n>0,'empty graph')
    for i,row in enumerate(a):
        need(0<=row<1<<n and not row&(1<<i),'bad row')
        for j in range(i):need(((row>>j)&1)+((a[j]>>i)&1)==1,'bad pair')

def literal(a,root=None):
    n=len(a);full=(1<<n)-1;h=[0]*(1<<n);h[0]=1
    end=[0]*(1<<n);start=0;calls=0;total=0
    def visit(last,mask,first):
        nonlocal calls,total,start
        calls+=1;h[mask]+=1
        if last==root:end[mask]+=1
        if mask==full:
            total+=1
            if first==root:start+=1
            return
        todo=a[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            visit(bit.bit_length()-1,mask|bit,first)
    for v in range(n):visit(v,1<<v,v)
    L=sum(end[mask]*h[full^mask] for mask in range(1,full+1)) if root is not None else None
    return total,start,L,calls

def child(a,root):
    core=[2,8,3,5];n=len(a);out=core[:]+[0]*(n-1)
    image={root:0}
    for i,v in enumerate(w for w in range(n) if w!=root):image[v]=i+4
    for i,j in combinations(range(n),2):
        x,y=image[i],image[j]
        if a[i]&(1<<j):out[x]|=1<<y
        else:out[y]|=1<<x
    for c in (1,2,3):
        for v in range(4,len(out)):out[c]|=1<<v
    valid(out)
    for i,j in combinations(range(n),2):
        need(bool(out[image[i]]&(1<<image[j]))==bool(a[i]&(1<<j)),'input not retained')
    return out

def main():
    seeds={};seed_orders=Counter();source_orders=Counter();j_hist=Counter()
    calls=graphs=0;keys=set();max_target=0;examples=[]
    seed_file=HERE/'data/SEEDS.tsv';profile_file=HERE/'data/PROFILE_CERTIFICATES.tsv'
    for line in seed_file.read_text().splitlines():
        target,n,*a=map(int,line.split());need(n==len(a) and target not in seeds,'seed record')
        valid(a);h,_,_,work=literal(a);calls+=work;graphs+=1
        need(h==target,'seed count')
        seeds[target]=a;seed_orders[n]+=1
    need(set(seeds)==set(range(23,410,2)) and len(seeds)==194,'seed interval')
    for line in profile_file.read_text().splitlines():
        B,r,j,root,n,*a=map(int,line.split())
        need(n==len(a) and 0<=root<n,'source size/root')
        need(B in range(137,998,2) and r in (1,2) and 0<=j<=8,'profile range')
        need((B,r) not in keys,'duplicate profile');keys.add((B,r));valid(a)
        h,s,L,work=literal(a,root);calls+=work;graphs+=1
        need(h+L-s==B-2*j and s==r+3*j,'physical profile equations')
        need(h+L-s>0 and (h+L-s)%2==1,'positive odd physical base')
        out=child(a,root);hc,_,_,work=literal(out);calls+=work;graphs+=1
        need(hc==3*B+2*r and hc==3*h+3*L-s,'literal child count')
        source_orders[n]+=1;j_hist[j]+=1;max_target=max(max_target,hc)
        if B in (137,139,995,997):examples.append({'B':B,'r':r,'j':j,'root':root,'source_order':n,'h':h,'s':s,'L':L,'child_order':len(out),'child_H':hc})
    need(keys=={(B,r) for B in range(137,998,2) for r in (1,2)},'profile interval')
    need(graphs==1918 and max_target==2995,'graph totals')
    boundaries=[]
    for n in (409,411,413,415):
        if n<=409:boundaries.append({'n':n,'route':'seed'})
        elif n%3==0:boundaries.append({'n':n,'route':'triangle','factor':n//3})
        else:
            r=1 if n%3==2 else 2;B=(n-2*r)//3
            need(B==137,'first residue base');boundaries.append({'n':n,'route':'profile','B':B,'r':r})
    print(json.dumps({'status':'PASS','seed_graphs':len(seeds),'seed_interval':[23,409],
                      'seed_order_histogram':dict(seed_orders),'profile_certificates':len(keys),
                      'profile_base_interval':[137,997],'both_residues':True,'J':8,
                      'source_order_histogram':dict(source_orders),'j_histogram':dict(j_hist),
                      'literal_graphs':graphs,'literal_DFS_visits':calls,'largest_child_target':max_target,
                      'seed_sha256':hashlib.sha256(seed_file.read_bytes()).hexdigest(),
                      'profile_sha256':hashlib.sha256(profile_file.read_bytes()).hexdigest(),
                      'induction_boundary_controls':boundaries,'examples':examples,
                      'all_order_profile_hypothesis':'OPEN','eventual_threshold':'CONDITIONAL ONLY'},indent=2,sort_keys=True))

if __name__=='__main__':main()
