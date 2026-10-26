#!/usr/bin/env python3
"""Separate signed-component verification of order-ten old-ledger dimensions.

The producer's n10 diagnostic reused the earlier independent Gaussian.
This second implementation uses no Gaussian methods for n10 equations.
"""
from core import old_words,endpair_subcells,path,need,literal_source,literal_ledger
from itertools import combinations,product
from collections import Counter
import json

class Signs:
    def __init__(self,n):
        self.n=n
        self.pairs=list(combinations(range(n),2))
        self.index={p:i+1 for i,p in enumerate(self.pairs)}
        self.parent=list(range(len(self.pairs)+1));self.sign=[0]*len(self.parent)
        self.size=[1]*len(self.parent);self.ok=True
    def copy(self):
        t=Signs(self.n);t.parent=self.parent[:];t.sign=self.sign[:];t.size=self.size[:];t.ok=self.ok
        return t
    def root(self,v):
        if self.parent[v]!=v:
            p,d=self.root(self.parent[v]);self.sign[v]^=d;self.parent[v]=p
        return self.parent[v],self.sign[v]
    def equal(self,a,b,d):
        a,x=self.root(a);b,y=self.root(b);d^=x^y
        if a==b:
            if d:self.ok=False
            return
        if self.size[a]<self.size[b]:a,b=b,a
        self.parent[b]=a;self.sign[b]=d;self.size[a]+=self.size[b]
    def literal(self,a,b,flips=frozenset()):
        p=(min(a,b),max(a,b))
        return self.index[p],int(a>b)^int(p in flips)
    def force(self,a,b,value):
        v,d=self.literal(a,b);self.equal(v,0,value^d)
    def anti(self,q,flips=frozenset()):
        r=dict(zip(q,reversed(q)))
        for a,b in self.pairs:
            v,d=self.literal(a,b,flips);w,e=self.literal(r[b],r[a],flips)
            self.equal(v,w,d^e)
    def forced(self,a,b,flips):
        v,d=self.literal(a,b,flips);r,e=self.root(v);z,f=self.root(0)
        return r==z and (d^e^f)==1
    def dimension(self):
        need(self.ok,'dimension of inconsistent signed family')
        return len({self.root(v)[0] for v in self.parent})-1
    def singleton(self):
        need(self.dimension()==0,'singleton decoder only')
        z,c=self.root(0);rows=[0]*self.n
        for (a,b),v in self.index.items():
            r,d=self.root(v);need(r==z,'all variables forced')
            u,w=(a,b) if d^c else (b,a);rows[u]|=1<<w
        return tuple(rows)

def source(n,I,opened,chord):
    t=Signs(n);t.anti(list(range(n)))
    for v in range(n-1):t.force(v,v+1,1)
    t.force(n-1,0,1)
    for v in range(1,n-1):t.force(0,v,v not in I)
    t.force(1,n-2,opened);t.force(2,n-2,chord)
    need(t.ok,'nonempty raw source')
    return t

def main():
    n=10;counts=Counter();dimensions=Counter();failed=[];models=set()
    for ports in combinations(range(2,n-1),3):
        I=set(ports)
        for opened,chord in product((False,True),repeat=2):
            s=source(n,I,opened,chord);dimensions[s.dimension()]+=1;counts['all_source_cells']+=1
            if not opened and {n-3,n-2}<=I:
                h=next(iter(I-{n-3,n-2}));leaves=endpair_subcells(s,h,chord);counts['endpair_cells']+=1
            else:
                label,qs=old_words(n,I,opened,chord);leaves=[(label,[],qs)];counts['base_cells']+=1
            refined=[]
            for label,guards,qs in leaves:
                t=s.copy()
                for a,b,value in guards:t.force(a,b,value)
                need(t.ok,'nonempty refined leaf')
                refined.append((guards,t.dimension()))
                counts['leaves']+=1;flips=set();blues=[]
                for step,q in enumerate(qs,1):
                    need(sorted(q)==list(range(n)),'complete old word')
                    need(all(t.forced(a,b,flips) for a,b in zip(q,q[1:])),'every old join uniformly directed')
                    trial=t.copy();trial.anti(q,flips)
                    if trial.ok:
                        need(trial.dimension()==0,'exact blue locus is singleton')
                        rows=trial.singleton();literal_source(rows,I,opened,chord)
                        actual,_,_=literal_ledger(rows,qs,False)
                        need(actual[step-1],'conditioned blue source literally blue')
                        models.add(rows)
                        blues.append({'step':step,'dimension':trial.dimension(),'source_rows':rows,'actual_old_blue_flags':actual})
                    flips^=set(t.pairs)^path(q)
                need(flips==set(t.pairs)^path(range(n)),'complete literal symbolic target')
                if blues:
                    counts['failed_leaves']+=1
                    failed.append({'ports':ports,'middle_open':opened,'x_to_c':chord,'label':label,
                                   'guards':guards,'dimension':t.dimension(),'blue_loci':blues})
                else:counts['passing_leaves']+=1
            need(sum(1<<d for _,d in refined)==1<<s.dimension(),'complete leaf dimension mass')
            for (a,_),(b,_) in combinations(refined,2):
                t=s.copy()
                for u,v,value in a+b:t.force(u,v,value)
                need(not t.ok,'disjoint leaf guards')
    need(dict(counts)=={'all_source_cells':140,'base_cells':130,'endpair_cells':10,'leaves':151,'passing_leaves':142,'failed_leaves':9},'exact n10 counts')
    need(dimensions=={10:140} and len(models)==9,'full source dimensions and nine distinct singleton failures')
    print(json.dumps({'status':'PASS / DIAGNOSTIC ONLY','method':'independent signed components, not the producer Gaussian','counts':dict(counts),'raw_dimensions':dict(dimensions),
                      'literal_marked_sources':140*(1<<10),'failed_leaves':failed,'all_joins_and_targets_valid':True,
                      'scope':'Complete dimensions-only old-menu diagnostic. Nine singleton old-blue sources remain; no replacement search and no connectivity obstruction claimed.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
