"""Exact nonnegative forest transfers with deferred active-active edges."""
from collections import Counter
from itertools import product, combinations, combinations_with_replacement

def patterns(width,p):
    pairs=list(product(range(-1,width),repeat=2));out=[]
    for size in range(width+p+1):
        for ps in combinations_with_replacement(pairs,size):
            if sum(a==-1 for a,b in ps)>p:continue
            if all(sum((a==u)+(b==u) for a,b in ps)<=1 or
                   ((u,u) in ps and sum(a==u or b==u for a,b in ps)==1)
                   for u in range(width)):out.append(ps)
    return out

def introduce(states,u):
    out=Counter()
    for ps,count in states.items():
        out[tuple(sorted(ps+((u,u),)))]+=count
        for i,(start,end) in enumerate(ps):
            if end==-1:
                new=list(ps);new[i]=(start,u)
                out[tuple(sorted(new))]+=count
    return out

def process_edge(states,a,b):
    out=Counter(states)
    for ps,count in states.items():
        for i,(start,end) in enumerate(ps):
            if end!=a:continue
            for j,(start2,end2) in enumerate(ps):
                if start2!=b or i==j:continue
                new=[part for k,part in enumerate(ps) if k not in (i,j)]
                new.append((start,end2));out[tuple(sorted(new))]+=count
    return out

def update(states,width,p,drop,mask,permutation=None):
    out=introduce(states,width)
    remain=[u for u in range(width+1) if u!=drop]
    for bit,u in enumerate(remain):
        out=process_edge(out,drop,u) if mask>>bit&1 else process_edge(out,u,drop)
    if permutation is not None:remain=[remain[i] for i in permutation]
    labels={-1:-1,drop:-1,**{u:i for i,u in enumerate(remain)}}
    result=Counter()
    for ps,count in out.items():
        q=tuple(sorted((labels[a],labels[b]) for a,b in ps))
        if sum(a==-1 for a,b in q)<=p:result[q]+=count
    return result

def complete(states,width,bits):
    for bit,(a,b) in enumerate(combinations(range(width),2)):
        states=process_edge(states,a,b) if bits>>bit&1 else process_edge(states,b,a)
    return states

def read(states,root=0):
    """h,s are valid at p>=1; L,B require p>=2 and an active root."""
    h=s=L=0
    for ps,count in states.items():
        if len(ps)==1:
            h+=count
            if ps[0][0]==root:s+=count
        if len(ps) in (1,2) and any(b==root for a,b in ps):L+=count
    return h,s,L,h+L-s
