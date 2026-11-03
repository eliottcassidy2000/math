"""Positive search through small two-entry cores and arbitrary declared tails."""
import sys
from itertools import combinations,product
from functools import lru_cache
from support import HERE,need,save,junction
sys.path.insert(0,str(HERE.parent/'ready_return_pentagon_20261103'))
import ready_inverse as inv
from research import core_data,tail_coefficients

@lru_cache(None)
def coefficients(t):
    words=[()] if t==0 else [w+(1,) for w in product((0,1),repeat=t-1)]
    return [(w,tail_coefficients(w,2)[0],tail_coefficients(w,1)[1:]) for w in words]

def main():
    bank=[]
    for n in range(3,9):
        cores=[(T,(0,)) for T in inv.two_start(n)]
        for T in inv.unique_cycle_classes(n):
            starts=[r for r in range(n) if inv.unique_start_path(T,r) is not None]
            cores.extend((T,J) for J in combinations(starts,2))
        for T,J in cores:
            data=core_data(T,J);bank.append((T,J,data));inv.tables.cache_clear()
        print('small core bank through',n,len(bank),flush=True)
    searched=[];hit=None
    for order in range(15,20):
        tests=0
        for T,J,d in bank:
            t=order-len(T)-2
            if t<0:continue
            for word,A,(D,B) in coefficients(t):
                tests+=1;h=A*(d['G']+d['p2'])+B*d['h']-D*d['e']-d['C']-d['F']-d['b']
                if h!=18979:continue
                P,pr=inv.funnel(T,set(J),word);U,r=inv.birth(P,pr);z=junction.profile(U,r)
                need(z[0]==h and z[3]-h==2 and junction.strong(U),'literal ready two-entry target')
                hit={'rows':U,'root':r,'profile':z,'parent_order':order-1,'budget':2,'core':T,'core_ports':J,'word':word,'core_data':d};break
            if hit:break
        searched.append({'child_order':order,'tests':tests});print('two-entry small-core child',order,tests,'hit',hit is not None,flush=True)
        if hit:break
    save('funnel_search.json',{'status':'FINITE_EXACT_POSITIVE_WITNESS' if hit else 'FINITE_EXACT_BOUNDED_LANGUAGE_ABSENCE','core_presentations':len(bank),'banks':searched,'hit':hit,
      'scope':'two-start one-port cores or unique-cycle two-port cores through order8; one additional funnel then newborn; child orders15..19; duplicates allowed'})
    if hit:save('ready18979.json',{'status':'FINITE_EXACT_POSITIVE_WITNESS','hit':hit,'scope':'small-core funnel search; no minimality claim'})
    print('PASS funnel search',None if hit is None else hit['profile'],flush=True)
if __name__=='__main__':main()
