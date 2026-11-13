"""Cross-certificate checks, selected-gate hostile, and exact spectral audit."""
from context import *
from math import comb,factorial
from itertools import product
from fractions import Fraction
from collections import Counter
from descent_grammar import count_capped
from independent_descent import tail
from endpoint_profiles import canonical_observer,advance
import sympy as sy
import time

def run():
    start=time.time();obs=load('endpoint_profiles.json.gz',HERE)['records'];split=load('split_endpoint_profiles.json.gz',HERE)['records']
    need(len({r['family'] for r in obs})==70 and len(obs)==78,'70 families,78 presentations')
    need(Counter(len(r['observer'][0]) for r in obs)=={2:39,4:39},'complete observer mode boundary')
    need(Counter(len(r['observer'][0]) for r in split)=={3:78},'split-clock mode boundary')
    quotient_checks=0
    for r in obs:
        X=tuple(map(tuple,r['observer']));padded=tuple(p+(0,0) for p in X)
        need(canonical_observer(padded)==X,'zero mode padding is not a new observer')
        for letter in range(12):
            need(canonical_observer(advance(padded,letter))==canonical_observer(advance(X,letter)),'canonical quotient commutes with every command')
            quotient_checks+=1
    a=next(r for r in obs if r['family']=='even:0');b=next(r for r in obs if r['family']=='even:2')
    need(a['observer']==b['observer'],'equal complete future ordinary-H observers')
    old=load('birth_sources_synchrony_12.json.gz')['records'];hostile=[]
    for idx,expectedU in [(0,11292),(2,17025)]:
        r=old[idx];h,u,w,_=tables(np.array(r['rows'],dtype=np.int64));O=[v for v in range(len(r['rows'])) if r['selected']>>v&1]
        U=sum(int(u[v]) for v in O);need(h==2663 and U==expectedU and defect(h,w,O)==1,'literal selected gate hostile')
        need(int(count_capped(np.array(tail(r['rows'],r['selected']),dtype=np.int64),2*U+2))==2*U+1,'literal hostile tail distinguishes selected gate')
        hostile.append({'family':f'even:{idx}','H':int(h),'U':U,'selected':O,'tail_H':1+2*U})
    m=17
    cell=max(factorial(m)*comb(m,r-1) for r in range(1,m+2))
    cell=max(cell,max(factorial(m-1)*comb(m,r) for r in range(m+1)))
    total=max(factorial(m)*comb(m+1,r) for r in range(m+2))
    need(cell<2**63 and total<2**64,'run-DP exact arithmetic bound')
    # A separate symbolic characteristic-polynomial computation checks every
    # word through length3, for every stride in the recorded modular search.
    env=load('two_port_envelope.json',HERE.parent/'adaptive_state_induction_20261103')
    M={a:sy.Matrix([[sy.Rational(s) for s in row] for row in env['rank']['query_matrices'][a]]) for a in [3,4,5,6,7,8]}
    periods=load('endpoint_periods.json.gz',HERE);lookup={(r['stride'],r['length']):{tuple(w) for w in r['modular_survivors']} for r in periods['records']};checks=0
    for length in range(1,4):
        for word in product(M,repeat=length):
            P=sy.eye(17)
            for a in word:P=P*M[a]
            cp=P.charpoly()
            for stride in range(1,4):
                if (stride,length) not in lookup:continue
                exact=all(cp.eval(x)==0 for x in [1,2**stride,3**stride])
                mod=all(int(cp.eval(x).p)%101==0 for x in [1,2**stride,3**stride])
                need(mod==(word in lookup[stride,length]),'independent symbolic polynomial/modular agreement')
                need(not exact or mod,'exact survivor never lost mod101');checks+=1
    for r in periods['records']:
        need(r['tested']==6**r['length'],'complete finite word universe')
        expected=[list(w) for w in product([5,8],repeat=r['stride'])] if r['length']==r['stride'] else []
        need(r['modular_survivors']==expected,'only static-root survivors')
    out={'scope':'Cross-checks and hostile controls; symbolic independent polynomial audit through period length3, not an independent rerun of every long word.','selected_gate_hostile':hostile,'exact_polynomial_tests':checks,'bounded_periods_total':sum(r['tested'] for r in periods['records']),'arithmetic_bounds':{'largest_run_cell':cell,'largest_run_sum':total},'ordinary_observer_classes':load('observer_rewrites.json.gz',HERE)['classes'],'old_remaining_repairs':473}
    out['zero_padding_quotient_controls']=quotient_checks
    save('verification.json.gz',out);print(out,'seconds',round(time.time()-start,2),flush=True)
if __name__=='__main__':run()
