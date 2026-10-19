"""Exhaustive small arithmetic tail recognition and routing, with actual edges."""
from itertools import product
from random import Random
from verify_partition_routes import route,states,check,random_state
import json

class DSU:
    def __init__(self,n):self.p=list(range(n))
    def find(self,x):
        while self.p[x]!=x:self.p[x]=self.p[self.p[x]];x=self.p[x]
        return x
    def join(self,x,y):
        x=self.find(x);y=self.find(y)
        if x==y:return False
        self.p[x]=y;return True

def main():
    fibers=routes=0;nonstar=[];noncomplete=[];counts=[];by_r={};family_routes=0
    for m in range(1,11):
        row=0
        for r in range(1,min(m,5)+1):
            for tail in product(*(range(m+1-d) for d in range(r+1,m+1))):
                dsu=DSU(m+1)
                if not all(dsu.join(a,a+d) for d,a in zip(range(r+1,m+1),tail)):continue
                contact=set();edges=[]
                for d in range(1,r+1):
                    color=[]
                    for a in range(m+1-d):
                        x,y=dsu.find(a),dsu.find(a+d)
                        pair=tuple(sorted((x,y))) if x!=y else None
                        color.append(pair)
                        if pair is not None:contact.add(pair)
                    edges.append(color)
                if len(contact)!=r:continue
                # Color one connects the quotient, so r distinct links on r+1
                # vertices certify that the simplification is a tree.
                ids={p:i for i,p in enumerate(sorted(contact))}
                A=[tuple(None if p is None else ids[p] for p in col) for col in edges]
                st=states(A);check(bool(st),'arithmetic staircase completion')
                endpoints=[(st[0],st[-1])]
                if len(st)>2:endpoints.append((st[len(st)//2],st[0]))
                for x,y in endpoints:route(A,x,y);routes+=1
                row+=1;fibers+=1;by_r[r]=by_r.get(r,0)+1
                degrees={v:sum(v in e for e in contact) for e in contact for v in e}
                if r>=3 and max(degrees.values())<r and len(nonstar)<3:
                    nonstar.append({'m':m,'r':r,'tail_offsets':tail,'contact_edges':sorted(contact)})
                if any(len({p for p in a if p is not None})<r for a in A) and len(st)>1 and len(noncomplete)<3:
                    noncomplete.append({'m':m,'r':r,'tail_offsets':tail,'states':len(st),'contact_supports':[sorted({p for p in a if p is not None}) for a in A]})
        counts.append({'m':m,'partition_contact_fibers':row})
    # Explicit higher-rank arithmetic positive controls from the proved
    # factorial family. These are constructed individually, not enumerated.
    rng=Random(20261019)
    for k in range(1,13):
        m=k*k+2*k;L=k+1;centers={i*L:i for i in range(1,k+1)}
        dsu=DSU(m+1)
        for d in range(k+1,m+1):
            e=(0,d) if d not in centers else (centers[d],d+centers[d])
            check(dsu.join(*e),'factorial fixed forest')
        contact=set();edges=[]
        for d in range(1,k+1):
            color=[]
            for a in range(m+1-d):
                x,y=dsu.find(a),dsu.find(a+d)
                pair=tuple(sorted((x,y))) if x!=y else None
                color.append(pair)
                if pair is not None:contact.add(pair)
            edges.append(color)
        check(len(contact)==k,'factorial contact tree')
        ids={p:i for i,p in enumerate(sorted(contact))}
        A=[tuple(None if p is None else ids[p] for p in col) for col in edges]
        for trial in range(10):
            route(A,random_state(A,rng),random_state(A,rng));family_routes+=1
    print(json.dumps({'status':'PASS' ,'all_m_at_most':10,'free_colors_at_most':5,'partition_contact_fibers':fibers,'endpoint_routes':routes,'by_m':counts,'by_free_rank':by_r,'constructed_factorial_max_rank':12,'constructed_factorial_routes':family_routes,'first_nonstar_contact_trees':nonstar,'first_noncomplete_assignment_fibers':noncomplete,'scope':'Exhaustive fixed upper-tail choices only in the displayed bounded range; loops retained as unavailable actual choices and parallel realizations remain distinct.'},indent=2))

if __name__=='__main__':main()
