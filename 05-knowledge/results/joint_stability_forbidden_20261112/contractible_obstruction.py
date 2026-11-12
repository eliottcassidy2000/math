"""Independent all-length and terminal-root audit of nonconfluence."""
from itertools import product,combinations
from pathlib import Path
import numpy as np
from census_nine import need,save,tables,defect,strong,axis_values
from joint_nine import cones,kernel,expo,evaluate,key,quotient
from obstruction import permutation_counts
HERE=Path(__file__).resolve().parent

def describe(rows,mask):
    h,u,w,_=tables(np.array(rows,dtype=np.int64));sel=[i for i in range(len(rows)) if mask>>i&1];d=defect(h,w,sel)
    ports=[];axes={}
    for v in range(len(rows)):
        if mask>>v&1:continue
        R=1+min(rows[v].bit_count(),len(rows)-1-rows[v].bit_count());vals=axis_values(rows,mask,v,R)
        axes[v]={'R':R,'values':vals}
        if all(x==1 for x in vals):ports.append(v)
    q,o,comps=quotient(rows,mask);eligible=[c for c in comps if len(c)>1 and set(c)&set(ports)]
    return {'rows':rows,'selected':mask,'strong':strong(rows),'H':int(h),'W':int(h)-d,'defect':d,'U':sum(int(u[i]) for i in sel),'stable_ports':ports,'axis_certificates':axes,'eligible_blocks':eligible,'canonical_key':key(rows,mask)}

def run():
    Q=[2,100,49,7,75,409,173,31,223];O=[1,4,5];P=[0,8];mask=50
    need(permutation_counts(Q,mask)==(665,660),'independent quotient permutation counts')
    D,U=cones.kernel(Q,P,O);need(D=={(1,1):5,(1,2):-2,(2,1):-2,(2,2):1},'independent full kernel')
    need(D==kernel(Q,mask,P),'signed-word kernel control');E=expo(D)
    grid=[]
    for a,b in product(range(1,5),repeat=2):
        rr,oo,_=cones.expand(Q,O,P,[a,b]);h,u,w,_=tables(np.array(rr,dtype=np.int64));d=defect(h,w,oo)
        need(d==1+(2**a-4)*(2**b-4)==evaluate(E,[a,b]),'literal full grid')
        grid.append({'lengths':[a,b],'H':int(h),'W':int(h)-d,'defect':d,'U':sum(int(u[i]) for i in oo)})
    states=[]
    for sizes in [[1,2],[2,1],[2,2]]:
        rr,oo,_=cones.expand(Q,O,P,sizes);r=describe(rr,sum(1<<v for v in oo));r['lengths']=sizes;states.append(r)
    a,b,c=states;need(a['strong'] and b['strong'] and c['strong'],'strong graph control')
    need(a['defect']==b['defect']==c['defect']==1,'unit states')
    need(not a['eligible_blocks'] and not b['eligible_blocks'],'both terminal births')
    need(a['canonical_key']!=b['canonical_key'],'nonisomorphic terminal roots')
    reverse=[((1<<10)-1)^a['rows'][i]^(1<<i) for i in range(10)]
    need(key(reverse,a['selected'])==tuple(b['canonical_key']),'converse terminal roots')
    need(c['eligible_blocks']==[[0,1],[9,10]],'exact competing contraction blocks')
    children=set()
    for block in c['eligible_blocks']:
        kept=[i for i in range(11) if i!=block[0]];rows=[sum(1<<j for j,v in enumerate(kept) if c['rows'][u]>>v&1) for u in kept];m=sum(1<<j for j,v in enumerate(kept) if c['selected']>>v&1);children.add(key(rows,m))
    need(children=={tuple(a['canonical_key']),tuple(b['canonical_key'])},'literal deletion edges')
    result={'scope':'All-length physical counterexample to joint stability of stable ports in different contractible unselected blocks and to ordinary stable-twin confluence. Minimality not asserted by this local proof.','quotient_rows':Q,'quotient_selected':mask,'quotient_ports':P,'quotient_H_W':[665,660],'defect_formula':'1+(2^a-4)*(2^b-4)','exact_positive_unit_locus':'a=2 or b=2','Stirling_terms':[[*k,v] for k,v in sorted(D.items())],'literal_grid':grid,'states':states,'distinct_terminal_births':2,'terminal_births_are_converses':True}
    save('contractible_obstruction.json.gz',result);print('Verified contractible obstruction, two nonisomorphic terminal births, full 4x4 literal grid.')
if __name__=='__main__':run()
