"""Intrinsic selected-budget invariants across every tournament through7."""
from common import *
from collections import Counter
from math import factorial
from weighted_entry_search import search
sys.path.insert(0,str(HERE.parent/'modular_centers_collatz_pell_20261111'))
import content_checks as cc
import tournament_checks as tc

def run(max_n=7):
    codes={0};levels=[];records=[];expected=[1,1,2,4,12,56,456,6880]
    need(1<=max_n<=8,'declared complete universe through8')
    for n in range(1,max_n+1):
        if n>1:
            codes={cc.canonical(cc.child(cc.decode(n-1,c),m)) for c in codes for m in range(1<<(n-1))}
        need(len(codes)==expected[n-1],'inherited A000568 universe control')
        totals=Counter();max_arity=0;nonflag=0;vectors=set();mass=0
        for code in sorted(codes):
            rr=cc.decode(n,code);h,u,w,masks=tables(np.array(rr,dtype=np.int64));target=int(h)-1
            costs=[0]*(1<<n);f=[0]*(n+1);minimal=[]
            for O in range(1,1<<n):
                b=O&-O;i=b.bit_length()-1;smaller=O^b
                costs[O]=costs[smaller]+sum(int(w[i,j]) for j in range(n) if smaller>>j&1)
            for O,cost in enumerate(costs):
                if cost<=target:f[O.bit_count()]+=1
                elif all(costs[O^(1<<j)]<=target for j in range(n) if O>>j&1):minimal.append(O)
            arity=max((O.bit_count() for O in minimal),default=0);max_arity=max(max_arity,arity);nonflag+=int(arity>2)
            for i,v in enumerate(f):totals[i]+=v
            vectors.add(tuple(f));aut=tc.automorphism_count(rr);mass+=factorial(n)//aut
            unit=list(map(int,masks));need(len(unit)%aut==0,'marked unit rigidity orbit divisor')
            found,stats=search(h,w);need(found==sorted(unit),'pruned unit search against full table')
            records.append({'n':n,'code':code,'H':int(h),'budget_f_vector':f,'minimal_nonfaces':minimal,'maximum_nonface_size':arity,'nonempty_unit_markings':len(unit),'marked_unit_orbits':len(unit)//aut,'automorphisms':aut,'search_visits':stats['visited_subsets']})
        need(mass==1<<(n*(n-1)//2),'complete labelled orbit mass')
        group=[r for r in records if r['n']==n]
        level={'n':n,'ordinary_types':len(codes),'budget_face_totals_by_size':[totals[i] for i in range(n+1)],'distinct_budget_f_vectors':len(vectors),'nonflag_budget_complexes':nonflag,'maximum_minimal_nonface_size':max_arity,'unit_marked_orbits':sum(r['marked_unit_orbits'] for r in group),'cores_with_unit_marking':sum(r['nonempty_unit_markings']>0 for r in group)}
        levels.append(level);print(level,flush=True)
    save('budget_complexes',{'scope':f'Every unmarked tournament class of orders1..{max_n}, no converse merging. Budget complex includes empty selection; unit-orbit counts exclude it. Face totals count vertex subsets on one representative of each ordinary class, not marked orbits. Counts are intrinsic under relabeling. A000568 and labelled orbit masses are controls, not new sequence terms.','levels':levels,'records':records})
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 7)
