"""Reflection-symmetric matching reversals in odd round tournaments.

Select an antipodal reflection-fixed triple {0,d,d/2}; reverse an invariant
perfect matching on its complement; add the two extreme ports. Their
uniform two-run bound makes four exact corner readings decisive at all
lengths. This is a declared structured universe, not an all-tournament scan.
"""
from itertools import combinations
from fractions import Fraction
from collections import Counter
import sys,time,json
import numpy as np
from census_nine import need,save,strong
from stable_census import selected_defects
from joint_nine import cones,key,quotient,kernel

def matchings(vertices):
    if not vertices:
        yield ();return
    a=vertices[0]
    for i,b in enumerate(vertices[1:],1):
        rest=vertices[1:i]+vertices[i+1:]
        for tail in matchings(rest):yield ((a,b),)+tail

def run(max_n=13):
    need(max_n in [7,9,11,13],'exact arithmetic universe through literal order17')
    records=[];levels=[];start=time.time()
    for n in range(7,max_n+1,2):
        stats=Counter();seen=set();base=[sum(1<<((i+j)%n) for j in range(1,(n+1)//2)) for i in range(n)]
        for d in range(1,n):
            center=d*pow(2,-1,n)%n;O={0,d,center};remaining=[v for v in range(n) if v not in O]
            for matching in matchings(remaining):
                M={tuple(sorted(e)) for e in matching};ref={tuple(sorted(((d-a)%n,(d-b)%n))) for a,b in M}
                if M!=ref:continue
                stats['reflection_matching_presentations']+=1;R=list(base)
                for a,b in M:R[a]^=1<<b;R[b]^=1<<a
                # Port p has unique outgoing neighbour0; q has unique
                # incoming neighbourd; q->p. All selected vertices are fixed.
                Q=[r for r in R]+[1,0];p=n;q=n+1
                for v in range(n):
                    if v!=0:Q[v]|=1<<p
                    if v==d:Q[v]|=1<<q
                    else:Q[q]|=1<<v
                Q[q]|=1<<p
                if not strong(Q):continue
                stats['strong_presentations']+=1;readings={}
                for a,b in [(1,1),(1,2),(2,1),(2,2)]:
                    rr,sel,_=cones.expand(Q,sorted(O),[p,q],[a,b]);m=sum(1<<v for v in sel)
                    readings[a,b]=int(selected_defects(np.array(rr,dtype=np.int64),np.array([m],dtype=np.int64))[0])
                need(readings[1,2]==readings[2,1],'physical marked converse symmetry')
                k=Fraction(readings[2,2]-2*readings[1,2]+readings[1,1],4)
                b=Fraction(readings[1,2]-readings[1,1],2)-2*k
                c=readings[1,1]-4*k-4*b
                need(k.denominator==b.denominator==c.denominator==1,'integer two-run coefficients')
                if not k:continue
                root=-b/k
                if root.denominator!=1 or root<2:continue
                root=int(root)
                if root&(root-1) or c!=1+k*root*root:continue
                A=root.bit_length()-1;stats['crossing_presentations']+=1;mask=sum(1<<v for v in O);canonical=key(Q,mask)
                if canonical in seen:continue
                seen.add(canonical);_,_,comps=quotient(Q,mask)
                # The independently compiled full kernel is inexpensive for
                # positive hits and certifies the two-run reading directly.
                K=kernel(Q,mask,[p,q])
                from tensor import expo
                E=expo(K);expected={(2,2):int(k),(2,1):-int(k)*root,(1,2):-int(k)*root,(1,1):1+int(k)*root*root}
                need(E==expected,'independent crossing kernel')
                r={'interior_order':n,'gate_distance':d,'selected':sorted(O),'reversed_matching':sorted(M),'interior_rows':R,'quotient_rows':Q,'quotient_selected':mask,'ports':[p,q],'crossing_length':A,'coefficient':int(k),'quotient_twin_blocks':[c for c in comps if len(c)>1],'corner_defects':[[a,b,v] for (a,b),v in readings.items()]}
                records.append(r);print('REFLECTION_CROSSING '+json.dumps(r),flush=True)
        levels.append({'interior_order':n,'statistics':dict(stats),'distinct_marked_quotients':len(seen)})
        print(json.dumps({'level':levels[-1],'seconds':round(time.time()-start,1)}),flush=True)
    out={'scope':f'Every reflection-invariant perfect-matching reversal on the complement of {{0,d,d/2}} in the odd round tournaments of orders7..{max_n}; all d=1..n-1. Two specified extreme unselected ports. Exact uniform two-run interpolation, confirmed by full kernels on positive hits.','levels':levels,'records':records}
    save('reflection_crossings.json.gz',out)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 13)
