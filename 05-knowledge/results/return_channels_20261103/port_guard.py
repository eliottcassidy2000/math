"""The independent port-start coordinate and its minimal Boolean legality guard."""
from itertools import permutations,combinations
import json,sympy as sp
from common import *

def starting(rows,port):
    rest=[i for i in range(len(rows)) if i!=port]
    return sum(all(rows[u]>>v&1 for u,v in zip((port,)+p,((port,)+p)[1:])) for p in permutations(rest))

def cover_factor(block):
    q=len(block)-1;b=0
    total=0
    for rest in permutations([i for i in range(len(block)) if i!=q]):
        p=(q,)+rest;j=p.index(b)
        total+=all(block[u]>>v&1 for k,(u,v) in enumerate(zip(p,p[1:])) if k!=j)
    return total

def join(old,block,port):
    n=len(old);m=len(block)
    rows=[r|(((1<<m)-1)<<n) for r in old]+[r<<n for r in block]
    rows[port]&=~(1<<n);rows[n]|=1<<port
    return tuple(rows)

def subset_start(rows,port):
    n=len(rows);dp=[{} for _ in range(1<<n)];dp[1<<port][port]=1
    for mask in range(1<<n):
        for u,c in dp[mask].items():
            opts=rows[u]&~mask
            while opts:
                bit=opts&-opts;opts-=bit;v=bit.bit_length()-1
                d=dp[mask|bit];d[v]=d.get(v,0)+c
    return sum(dp[-1].values())

def main():
    factors={};tests=0
    for l in MATS:
        block=actual_target((l,));k=cover_factor(block);psi=starting(block,len(block)-1)
        need(k>=psi>=1,('positive actual block port',l))
        need(starting(join((0,),block,0),len(block))==k,'independent singleton join factor')
        factors[l]={'kappa':k,'first_atom_start_count':psi}
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for code in range(1<<len(pairs)):
            rows=[0]*n
            for bit,(u,v) in enumerate(pairs):
                if code>>bit&1:rows[u]|=1<<v
                else:rows[v]|=1<<u
            for port in range(n):
                s=starting(rows,port)
                for l in ('a','b','q1','q4'):
                    block=actual_target((l,));child=join(rows,block,port)
                    need(subset_start(child,len(child)-1)==s*factors[l]['kappa'],'all small pointed graph joins')
                    tests+=1
    # All65 blocks on an actual start-zero source, and on a nontrivial start-positive source.
    for rows,port in (((10,12,9,0),3),((2,4,1),0)):
        s=starting(rows,port)
        for l in MATS:
            child=join(rows,actual_target((l,)),port)
            need(subset_start(child,len(child)-1)==s*factors[l]['kappa'],'all block guard controls');tests+=1
    extended=(26,28,25,0,8);representative=actual_target(('q1',))
    need(starting(extended,4)==0 and starting(representative,3)==1,'same hxy but different start feasibility')
    U=('r5','a','a');rows=[vector(U+w)+(subset_start(actual_target(U+w),sum(COSTS[l] for l in U+w)-1),) for w in ((),('b',),('q2',))]
    # The equal-hxy pair supplies the independent fourth coordinate.
    witnesses=rows+[(3,3,3,0),(3,3,3,1)];rank=sp.Matrix(witnesses).rank();need(rank==4,'sharp four-coordinate linear lift')
    result={'status':'PASS','factors':factors,'literal_pointed_joins':tests,
            'same_observer_control':{'extended_rows':extended,'bank_rows':representative,'start_counts':[0,1]},
            'linear_rank_witness':witnesses,'linear_rank':rank,
            'scope':'Coarsest future (H,port-start-feasibility) observer is hxy plus one bit for these joins; full source reconstruction is not claimed.'}
    (HERE/'port_guard.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','literal_joins':tests,'kappa_values':sorted({d['kappa'] for d in factors.values()}),'linear_rank':rank}))
if __name__=='__main__':main()
