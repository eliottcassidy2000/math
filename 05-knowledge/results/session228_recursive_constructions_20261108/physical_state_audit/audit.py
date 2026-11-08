#!/usr/bin/env python3
"""Independent literal path-list verification. No producer or repository imports."""
import itertools,json,hashlib,random
from collections import defaultdict

def require(x,label):
    if not x: raise RuntimeError(label)
def graph(n,c):
    a=[set() for _ in range(n)]
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i if (c>>k)&1 else j].add(j if (c>>k)&1 else i)
    return a
def paths(a,vertices):
    out=defaultdict(list);out[frozenset()].append(())
    def extend(p,unused):
        if p:out[frozenset(p)].append(p)
        for v in unused:
            if not p or v in a[p[-1]]:extend(p+(v,),unused-{v})
    extend((),set(vertices));return out
def state(a,vs):
    ps=paths(a,vs);full=frozenset(vs)
    return (len(ps[full]),sum(p[-1]==0 for p in ps[full]),sum(len(ps[full-s])*sum(p[0]==0 for p in pl) for s,pl in list(ps.items()) if 0 in s),len(ps[full-{0}]))
def blocks(a):
    remaining=set(range(1,len(a)));reach={v:{v}|(a[v]&remaining) for v in remaining}
    changed=True
    while changed:
        changed=False
        for v in remaining:
            r=set.union(*(reach[w] for w in reach[v]))
            if r!=reach[v]:reach[v]=r;changed=True
    bs=[]
    while remaining:
        v=min(remaining);b={w for w in remaining if w in reach[v] and v in reach[w]};bs.append(b);remaining-=b
    old=bs[:]
    bs.sort(key=lambda b:sum(next(iter(c)) in a[next(iter(b))] for c in old if c!=b),reverse=True)
    return bs
def audit(a):
    x=state(a,range(len(a)));bs=blocks(a);d=1;v=0;activated=False;rm=set();prefix={0};rows=[]
    for b in bs:
        ps=paths(a,b);h=len(ps[frozenset(b)])
        q=sum(len(pl)*len(ps[frozenset(b)-s]) for s,pl in list(ps.items()))
        alpha=sum(sum(p[0] in a[0] for p in pl)*len(ps[frozenset(b)-s]) for s,pl in list(ps.items()) if s)
        require(q>=2**len(b),'Q>=2^size')
        require(bool(alpha)==bool(b&a[0]),'activation')
        if not activated and not alpha and len(b)==1:rm|=b
        activated|=bool(alpha);d,v=d*h,v*q+alpha*d;prefix|=b
        require(state(a,prefix)[2:]==(d+v,d),'prefix recurrence')
        if len(b)>1:require(len(b)<=h,'Camion')
        rows.append((len(b),h,q,alpha))
    require(x==state(a,set(range(len(a)))-rm),'four-coordinate compression')
    require(x[2:]==(d+v,d),'final recurrence')
    if v==0:
        require(x==(d,d,d,d),'sink diagonal');bound=1 if d==1 else d+1
    else:bound=(0 if d==1 else d)+v.bit_length()+1
    require(len(a)-len(rm)<=bound,'order bound')
    require(sum(len(b) for b in bs if len(b)>1)<=(d if d>1 else 0),'strong budget')
    return [x,rows,sorted(rm)]
def main():
    digest=hashlib.sha256();totals={};rng=random.Random(228)
    for n in range(1,8):
        codes=range(1<<(n*(n-1)//2)) if n<=5 else sorted({0,(1<<(n*(n-1)//2))-1}|{rng.randrange(1<<(n*(n-1)//2)) for _ in range(128)})
        count=0
        for code in codes:
            result=audit(graph(n,code));count+=1
            digest.update(json.dumps([n,code,result],separators=(',',':')).encode())
        totals[n]=count
    # Long inactive prefix: each new source dominates existing rooted tournament.
    a=graph(4,17);baseline=state(a,range(4))
    for _ in range(8):a.append(set(range(len(a))))
    result=audit(a);require(state(a,range(len(a)))==baseline,'long source prefix')
    print(json.dumps({'status':'PASS','method':'literal enumeration of every directed path; independently computed prefix states','all_root0_labeled_through_n':5,'cases_by_order':totals,'long_prefix_order':len(a),'digest_sha256':digest.hexdigest()},indent=2))
if __name__=='__main__':main()
