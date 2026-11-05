#!/usr/bin/env python3
"""Exact endpoint-typed path-cover gluing checks; no external dependencies."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import factorial
import hashlib,json,pathlib

HERE=pathlib.Path(__file__).resolve().parent

def require(condition,message):
    if not condition: raise RuntimeError(message)

def tournament(n,mask):
    a=[0]*n;p=0
    for i in range(n):
        for j in range(i+1,n):
            if mask>>p&1:a[i]|=1<<j
            else:a[j]|=1<<i
            p+=1
    return tuple(a)

def induced(a,vertices):
    return tuple(sum(((a[v]>>u)&1)<<j for j,u in enumerate(vertices)) for v in vertices)

def path_profile(a,colors,k):
    """Unordered spanning nonempty directed path covers, indexed by k*k type counts.
    Works for arbitrary loopless digraphs. Components are chosen by smallest remaining vertex.
    """
    n=len(a);zero=(0,)*(k*k)
    if not n:return {zero:1}
    d=[{} for _ in range(1<<n)]
    for v in range(n):d[1<<v][(colors[v],v)]=1
    h=[Counter() for _ in range(1<<n)]
    for s in range(1,1<<n):
        for (start,v),count in d[s].items():
            h[s][start*k+colors[v]]+=count
            todo=a[v]&~s
            while todo:
                b=todo&-todo;todo-=b;w=b.bit_length()-1
                key=(start,w);d[s|b][key]=d[s|b].get(key,0)+count
    @lru_cache(None)
    def covers(s):
        if not s:return {zero:1}
        root=s&-s;rest=s^root;sub=rest;out=Counter()
        while True:
            part=sub|root
            for t,c in h[part].items():
                for sig,m in covers(s^part).items():
                    ss=list(sig);ss[t]+=1;out[tuple(ss)]+=c*m
            if not sub:break
            sub=(sub-1)&rest
        return dict(out)
    return covers((1<<n)-1)

def hamilton(a):
    n=len(a)
    if not n:return 1
    f=[[0]*n for _ in range(1<<n)]
    for v in range(n):f[1<<v][v]=1
    for s in range(1,1<<n):
        for v,c in enumerate(f[s]):
            if not c:continue
            todo=a[v]&~s
            while todo:
                b=todo&-todo;todo-=b;f[s|b][b.bit_length()-1]+=c
    return sum(f[-1])

def profile_h(p):return sum(c for s,c in p.items() if sum(s)==1)

def cross_join(a,b,ca,cb,k,C):
    na=len(a);out=list(a)+[x<<na for x in b]
    for i,ci in enumerate(ca):
        for j,cj in enumerate(cb):
            if C[ci*k+cj]:out[i]|=1<<(na+j)
            else:out[na+j]|=1<<i
    return tuple(out)

@lru_cache(None)
def kernel(ma,mb,k,C,maps=None,l=None):
    # A fragment has a side and two endpoint colors. Its identity is retained.
    if maps is None:maps=(tuple(range(k)),tuple(range(k)))
    if l is None:l=k
    types=[]
    for side,m in enumerate((ma,mb)):
        for t,count in enumerate(m):types += [(side,t//k,t%k)]*count
    a=[0]*len(types);starts=[];ends=[]
    for side,s,t in types:starts.append(maps[side][s]);ends.append(maps[side][t])
    for i,(side,s,t) in enumerate(types):
        for j,(other,u,v) in enumerate(types):
            if side==other:continue
            edge=C[t*k+u] if side==0 else not C[u*k+t]
            if edge:a[i]|=1<<j
    # Need start color of first fragment and END color of last fragment,
    # not one color assigned to the whole fragment. Use the same subset recurrence,
    # with custom start/end tags. This is a digraph of fragment identities.
    n=len(a);zero=(0,)*(l*l);d=[{} for _ in range(1<<n)]
    for v in range(n):d[1<<v][(starts[v],v)]=1
    h=[Counter() for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for (s,v),count in d[mask].items():
            h[mask][s*l+ends[v]]+=count
            todo=a[v]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;w=bit.bit_length()-1
                key=(s,w);d[mask|bit][key]=d[mask|bit].get(key,0)+count
    @lru_cache(None)
    def rec(mask):
        if not mask:return {zero:1}
        root=mask&-mask;rest=mask^root;sub=rest;out=Counter()
        while True:
            part=sub|root
            for t,count in h[part].items():
                for sig,m in rec(mask^part).items():
                    ss=list(sig);ss[t]+=1;out[tuple(ss)]+=count*m
            if not sub:break
            sub=(sub-1)&rest
        return dict(out)
    return rec((1<<n)-1)

def compose(pa,pb,k,C,maps=None,l=None):
    out=Counter()
    for ma,ca in pa.items():
        for mb,cb in pb.items():
            for sig,c in kernel(ma,mb,k,tuple(C),maps,l).items():out[sig]+=ca*cb*c
    return dict(out)

def relabel(profile,k,mapping,l):
    out=Counter()
    for sig,c in profile.items():
        ns=[0]*(l*l)
        for t,m in enumerate(sig):ns[mapping[t//k]*l+mapping[t%k]]+=m
        out[tuple(ns)]+=c
    return dict(out)

def uniform_rule(a,L,R,colors,k):
    rule=[None]*(k*k)
    for x in L:
        for y in R:
            t=colors[x]*k+colors[y];v=(a[x]>>y)&1
            if rule[t] is not None and rule[t]!=v:return None
            rule[t]=v
    return tuple(0 if x is None else x for x in rule)

def prime(a):
    n=len(a)
    for s in range(1,1<<n):
        if not 2<=s.bit_count()<n:continue
        if all((a[v]&s) in (0,s) for v in range(n) if not(s>>v&1)):return False
    return True

def alternating(r):
    n=2*r+1;a=[0]*n
    for i in range(1,n):
        for j in range(i+1,n):a[i]|=1<<j
        if i%2:a[0]|=1<<i
        else:a[i]|=1
    return tuple(a)

def two_apex(word):
    n=2+len(word);a=[0]*n;a[0]|=2
    for i,c in enumerate(word,2):
        for j in (0,1):
            if int(c)>>j&1:a[j]|=1<<i
            else:a[i]|=1<<j
        for j in range(2,i):a[j]|=1<<i
    return tuple(a)

def run():
    census=[]
    for n in range(1,5):
        tests=0;profiles=0
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask)
            for colors in product(range(2),repeat=n):
                full=path_profile(a,colors,2);profiles+=1
                require(profile_h(full)==hamilton(a),'H coefficient disagrees')
                require(sum(c for sig,c in full.items() if sum(sig)==n)==1,'singleton cover normalization')
                require(relabel(full,2,(0,0),1)==path_profile(a,(0,)*n,1),'relabel failed')
                # Every unordered proper cut once, using vertex zero on the left.
                for cut in range(1,(1<<n)-1,2):
                    L=[v for v in range(n) if cut>>v&1];R=[v for v in range(n) if not(cut>>v&1)]
                    C=uniform_rule(a,L,R,colors,2)
                    if C is None:continue
                    pa=path_profile(induced(a,L),tuple(colors[v] for v in L),2)
                    pb=path_profile(induced(a,R),tuple(colors[v] for v in R),2)
                    require(compose(pa,pb,2,C)==full,f'glue failed n={n}, mask={mask}, colors={colors}, cut={cut}')
                    tests+=1
        census.append({'order':n,'profiles':profiles,'uniform_colored_cuts':tests})
    # Side-specific output recoloring, including splitting same input color by operand.
    recolor_tests=0
    for C in product((0,1),repeat=4):
        a=tournament(3,2);b=tournament(2,1);ca=(0,1,0);cb=(0,1)
        glued=cross_join(a,b,ca,cb,2,C)
        for lm in product(range(2),repeat=2):
            for rm in product(range(2),repeat=2):
                result=compose(path_profile(a,ca,2),path_profile(b,cb,2),2,C,(lm,rm),2)
                expected=path_profile(glued,tuple(lm[c] for c in ca)+tuple(rm[c] for c in cb),2)
                require(result==expected,'side-specific recolor mismatch')
                recolor_tests+=1
    # Full-profile recursive prime family, not merely scalar checks.
    p=path_profile((0,),(0,),2);a=(0,);colors=(0,)
    recursive=[]
    singleton=path_profile((0,),(1,),2)
    for m in range(1,8):
        C=(0,1 if m%2 else 0,0,1)
        p=compose(p,singleton,2,C)
        a=cross_join(a,(0,),colors,(1,),2,C);colors+=(1,)
        require(p==path_profile(a,colors,2),'recursive prime-profile mismatch')
        if m%2==0:
            r=m//2;h=(4**(r+1)+6*r+5)//9
            require(a==alternating(r) and prime(a) and profile_h(p)==h,'prime family control')
            recursive.append({'r':r,'order':len(a),'H':h,'profile_terms':len(p)})
    # Scalar and module checks on larger family members, independently of cover algebra.
    prime_controls=[]
    for r in range(1,7):
        a=alternating(r);h=hamilton(a)
        require(prime(a),'all-order prime family finite control')
        require(9*h==4**(r+1)+6*r+5,'closed form failed')
        prime_controls.append({'r':r,'order':len(a),'H':h})
    # H13 prime and non-module input retained in a physical growth operation.
    a=two_apex('120');colors=(0,1,2,2,2);C=(0,0,1,0,0,0,0,0,1)
    pa=path_profile(a,colors,3);pb=path_profile((0,),(2,),3)
    child=cross_join(a,(0,),colors,(2,),3,C)
    require(child==two_apex('1201'),'physical append convention')
    pchild=compose(pa,pb,3,C)
    require(pchild==path_profile(child,colors+(2,),3),'H13-child all-cover profile')
    require(prime(a) and prime(child) and profile_h(pa)==13 and profile_h(pchild)==23,'H13/H23 values')
    # Heterogeneous operands: a prime alternating bank, and a prime two-apex bank.
    a=alternating(2);ca=(0,1,1,1,1);b=two_apex('120');cb=(0,1,2,2,2)
    C=(0,1,1,1,1,1,1,1,1) # only root-root arc reversed relative to A -> B
    joined=cross_join(a,b,ca,cb,3,C)
    # For size ten use the Hamilton slice from summing complete covers on small operands.
    mixed=compose(path_profile(a,ca,3),path_profile(b,cb,3),3,C)
    mixed_h=hamilton(joined)
    require(profile_h(mixed)==mixed_h,'mixed-bank gluing H')
    receipt={'status':'PASS','census':census,'prime_recursive_profiles':recursive,'prime_scalar_controls':prime_controls,
        'H13_H23':{'parent_rows':list(two_apex('120')),'child_rows':list(child),'H':[13,23]},
        'heterogeneous_join':{'left_rows':list(a),'right_rows':list(b),'cross_rule':list(C),'H':mixed_h,'profile_terms':len(mixed)},
        'side_specific_recolor_tests':recolor_tests,'kernel_cache':str(kernel.cache_info()),'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__':run()
