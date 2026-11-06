#!/usr/bin/env python3
"""Literal permutation/cut path covers and literal DFS, no producer imports."""
from itertools import permutations, combinations
from collections import Counter
import json

CASES = [
    ('A', [2,24,3,37,13,23], 15, [2,8,67,101,15,87,19], 79, (1,4)),
    ('B', [32,5,1,39,47,70,31], 5, [0,5,129,39,47,199,159,27], 5, (0,5)),
    ('C', [2,4,41,3,15,91,31], 27, [0,5,41,131,15,219,159,23], 155, (0,1)),
]

def need(b,s):
    if not b:
        raise RuntimeError(s)

def valid(a):
    n=len(a)
    for i,row in enumerate(a):
        need(0<=row<1<<n and not row&(1<<i),'bad row')
        for j in range(i):
            need(((row>>j)&1)+((a[j]>>i)&1)==1,'bad pair')

def all_modules(n):
    edges=list(combinations(range(n),2))
    for code in range(1<<len(edges)):
        a=[0]*n
        for bit,(i,j) in enumerate(edges):
            if code&(1<<bit):
                a[i]|=1<<j
            else:
                a[j]|=1<<i
        yield code,a

def weak_lengths(n,q):
    if q==0:
        if n==0:
            yield ()
        return
    for first in range(n+1):
        for rest in weak_lengths(n-first,q-1):
            yield (first,)+rest

def cover_count(a,q):
    # Each ordered tuple of paths has one concatenated permutation and cut vector.
    n=len(a); count=0
    for lengths in weak_lengths(n,q):
        for word in permutations(range(n)):
            pos=0; good=True
            for size in lengths:
                part=word[pos:pos+size];pos+=size
                if any(not a[x]&(1<<y) for x,y in zip(part,part[1:])):
                    good=False;break
            if good:
                count+=1
    return count

def join_sink(a):
    n=len(a)
    return [row|(1<<n) for row in a]+[0]

def inflate(core,mask,u):
    n=len(core); a=core[:]+[row<<n for row in u]
    for v in range(n,n+len(u)):
        for w in range(n):
            if mask&(1<<w):a[v]|=1<<w
            else:a[w]|=1<<v
    valid(a);return a

def restrict(a,order):
    return [sum(1<<j for j,v in enumerate(order) if a[w]&(1<<v)) for w in order]

def literal(a):
    full=(1<<len(a))-1; calls=0
    def visit(last,used,root):
        nonlocal calls
        calls+=1
        if used==full:
            return 1,int(root==0 and bool(a[last]&1))
        h=c=0;todo=a[last]&~used
        while todo:
            bit=todo&-todo;todo-=bit
            nh,nc=visit(bit.bit_length()-1,used|bit,root)
            h+=nh;c+=nc
        return h,c
    h=c=0
    for root in range(len(a)):
        nh,nc=visit(root,1<<root,root);h+=nh;c+=nc
    return h,c,calls

def cover_profile_by_assignments(a):
    # Independent higher-coordinate hostile: enumerate every simple path,
    # then assign each labelled vertex to a bin and multiply path choices.
    from itertools import product
    n=len(a);h=[0]*(1<<n);h[0]=1
    def visit(last,used):
        h[used]+=1
        todo=a[last]&~used
        while todo:
            bit=todo&-todo;todo-=bit
            visit(bit.bit_length()-1,used|bit)
    for v in range(n):visit(v,1<<v)
    result=[int(n==0)]
    for q in range(1,4):
        total=0
        for assignment in product(range(q),repeat=n):
            bins=[0]*q
            for v,label in enumerate(assignment):bins[label]|=1<<v
            ways=1
            for mask in bins:ways*=h[mask]
            total+=ways
        result.append(total)
    return result

def main():
    modules=graphs=rewrites=dfs_calls=0
    order_counts=Counter(); profiles=Counter(); samples=[]
    empty_rebase=None
    for n in range(5):
        for code,u in all_modules(n):
            valid(u);modules+=1;order_counts[n]+=1
            b=[cover_count(u,q) for q in range(4)]
            b0,h,b2,b3=b
            need(h>0,'Hamilton path existence control')
            need(b0==int(n==0),'empty boundary')
            profiles[tuple(b)]+=1
            up=join_sink(u)
            padded=[cover_count(up,q) for q in range(4)]
            need(padded==[q*b[q] for q in range(4)],'all-q single-sink padding law')
            up2=join_sink(up)
            for q in range(4):
                need(cover_count(up2,q)==q*q*b[q],'all-q two-sink padding law')
            if n:
                need(b2%2==0 and b2>=(n+1)*h,'coupling bound')
            observed=[]
            for name,p,pm,t,tm,flip in CASES:
                pp=inflate(p,pm,u);tt=inflate(t,tm,u)
                core_n=len(p)
                old_order=list(range(core_n))+list(range(core_n+1,len(tt)))
                old=restrict(tt,old_order)
                changed=[(i,j) for i,j in combinations(range(len(pp)),2)
                         if bool(pp[i]&(1<<j))!=bool(old[i]&(1<<j))]
                need(changed==[flip] and pp[flip[0]]&(1<<flip[1]),'physical edit')
                ph,pc,work=literal(pp);dfs_calls+=work;graphs+=1
                th,tc,work=literal(tt);dfs_calls+=work;graphs+=1
                if name=='A':
                    expected=(h+17*b2+8*b3-b0,3*h+34*b2+16*b3-2*b0,h+b2,h+2*b2)
                    rh,rc=th,tc
                else:
                    need(tt[0]==0,'universal sink')
                    rr=restrict(tt,list(range(1,len(tt))))
                    rh,rc,work=literal(rr);dfs_calls+=work;graphs+=1
                    need((rh,tc)==(th,0),'sink deletion')
                    # Retained old vertices agree, including all edges inside arbitrary U.
                    retained_parent=restrict(pp,list(range(1,len(pp))))
                    retained_child=restrict(rr,list(range(core_n-1))+list(range(core_n,len(rr))))
                    need(retained_parent==retained_child,'retained old graph')
                    need(len(rr)==len(pp),'same-order replacement')
                    expected=(5*h+20*b2,11*h+40*b2,h,3*h) if name=='B' else (27*b2-h-b0,54*b2-h-2*b0,h,3*h)
                need((ph,rh,pc,rc)==expected,'literal formula mismatch')
                need(rh==2*ph+h and pc>0 and rc>0,'correction/strongness')
                if name=='A':need(rc==2*pc-h,'A cycle formula')
                rewrites+=1
                observed.append({'family':name,'parent_H':ph,'child_H':rh,'parent_C':pc,'child_C':rc})
            cparent=inflate(CASES[2][1],CASES[2][2],up)
            rebased_h,rebased_c,work=literal(cparent);dfs_calls+=work;graphs+=1
            if n:
                need((rebased_h,rebased_c)==(observed[2]['child_H'],h),'C rebase')
            else:
                empty_rebase={'child_empty_H':observed[2]['child_H'],'parent_singleton_H':rebased_h}
                need(empty_rebase=={'child_empty_H':51,'parent_singleton_H':53},'empty rebase hostile')
            if n==0 or h>1 and not any(s['b']==b for s in samples):
                samples.append({'module_order':n,'module_code':code,'module_rows':u,'b':b,'rewrites':observed})
    need(modules==76 and graphs==684 and rewrites==228,'universe counts')
    # Equal H is not permission to treat the two graphs as identical modules.
    ccore=CASES[2]
    cchild=inflate(ccore[3],ccore[4],[0])
    replacement=restrict(cchild,list(range(1,len(cchild))))
    rebased=inflate(ccore[1],ccore[2],join_sink([0]))
    replacement_profile=cover_profile_by_assignments(replacement)
    rebased_profile=cover_profile_by_assignments(rebased)
    need(replacement_profile==[0,107,2600,26235],'replacement profile hostile')
    need(rebased_profile==[0,107,4804,73041],'rebased profile hostile')
    print(json.dumps({'status':'PASS','method':'literal permutation/cut path covers and nonmemoized directed-path DFS',
                      'modules':modules,'module_order_counts':dict(order_counts),'literal_graphs':graphs,
                      'rewrite_pairs':rewrites,'literal_DFS_calls':dfs_calls,
                      'profile_histogram':{str(k):v for k,v in sorted(profiles.items())},
                      'empty_rebase_hostile':empty_rebase,
                      'equal_H_unequal_profile_hostile':{'replacement_b':replacement_profile,'rebased_parent_b':rebased_profile},
                      'samples':samples},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
