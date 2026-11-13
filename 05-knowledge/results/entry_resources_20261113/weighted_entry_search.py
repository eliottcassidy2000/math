"""Exact unit-marking search as a nonnegative quadratic subset-sum problem."""
from common import *
from math import gcd
from functools import reduce

def search(H,W,use_gcd=True):
    n=len(W);target=int(H)-1;W=[[int(x) for x in row] for row in W]
    adjacent=[sum(1<<j for j in range(n) if j!=i and W[i][j]<=target) for i in range(n)]
    result=[];visits=0;gcd_cuts=0;budget_cuts=0
    def visit(chosen,candidates,cost):
        nonlocal visits,gcd_cuts,budget_cuts
        visits+=1
        if chosen and cost==target:result.append(chosen)
        options=[];filtered=0
        while candidates:
            bit=candidates&-candidates;candidates^=bit;v=bit.bit_length()-1
            marginal=sum(W[v][j] for j in range(n) if chosen>>j&1)
            if cost+marginal<=target:options.append((v,marginal));filtered|=bit
            else:budget_cuts+=1
        if use_gcd and cost<target and options:
            g=0
            for v,m in options:g=gcd(g,m)
            for i,(v,_) in enumerate(options):
                for u,_ in options[i+1:]:
                    if adjacent[v]>>u&1:g=gcd(g,W[v][u])
            if g and (target-cost)%g:gcd_cuts+=1;return
        for v,m in options:
            filtered^=1<<v
            visit(chosen|(1<<v),filtered&adjacent[v],cost+m)
    visit(0,(1<<n)-1,0)
    return sorted(result),{'visited_subsets':visits,'gcd_cuts':gcd_cuts,'budget_cuts':budget_cuts,'all_subsets':1<<n}

def run():
    bank=load(HERE.parent/'endpoint_coverage_compiler_20261111/inflation_bank_7.json')['records']
    unique={tuple(r['rows']) for r in bank}
    chain=load(HERE/'moving_root.json.gz')['records'];sources=sorted(unique,key=lambda t:(len(t),t))+[tuple(r['rows']) for r in chain]
    rows=[]
    for rr in sources:
        h,u,w,wanted=tables(np.array(rr,dtype=np.int64))
        found,stats=search(h,w);plain,plain_stats=search(h,w,False)
        need(found==plain==sorted(map(int,wanted)),'full subset recurrence and both exact pruned searches agree')
        rows.append({'rows':rr,'H':int(h),'unit_masks':found,'pruned':stats,'without_gcd':plain_stats})
    # Zero weights are allowed by the search contract, even if a producer's
    # actual class imposes stricter inequalities. Do not stop at equality.
    z,stats=search(1,[[0]*4 for _ in range(4)]);need(z==list(range(1,16)),'all zero-cost extensions retained')
    save('weighted_entry_search',{'scope':'Exact algorithm for symmetric nonnegative pair weights; unit sets satisfy sum_{i<j in O}W_ij=H-1. Complete agreement on all distinct inherited stable-bank tournaments through7 and eight moving-root cores12..19. This accelerates selected-marking tests after endpoint tables exist; it is not a faster complete tournament generator.','records':rows,'zero_weight_hostile_control':len(z)})
    print('cores',len(rows),'chain stats',[(len(r['rows']),r['pruned']) for r in rows[-8:]],flush=True)
if __name__=='__main__':run()
