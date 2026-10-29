"""Weighted stabilizer recurrences; literals only in the declared small universes."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
import importlib.util
import json

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('oldcanon',HERE.parent/'weighted_crt_clocks_20261028/check_canonical.py')
old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

# Bivariate monomial: selected-vertex weight, exponent of stabilizer order3.
def mul(p,q):
    out=Counter()
    for (a,e),v in p.items():
        for (b,f),w in q.items(): out[a+b,e+f]+=v*w
    return dict(out)

def step(p):
    cube=Counter(mul(mul(p,p),p)); diagonal={(3*a,3*e):v for (a,e),v in p.items()}
    for key,v in diagonal.items(): cube[key]-=v
    need(all(v%3==0 for v in cube.values()),'nonconstant cyclic orbits')
    out=Counter({key:v//3 for key,v in cube.items() if v})
    for (a,e),v in diagonal.items(): out[a,e+1]+=v
    return dict(out)

def eval_stabilizer(p,t): return sum(v*t**e for (a,e),v in p.items())
def output_poly(p): return [{'weight':a,'stabilizer_exponent':e,'orbits':v} for (a,e),v in sorted(p.items())]

@lru_cache(None)
def poly_a(k,palette):
    p={(a,0):v for a,v in enumerate(palette) if v}
    for _ in range(k): p=step(p)
    return p

def poly_c(k,palette):
    seed={(a,0):v for a,v in enumerate(palette) if v}; pair=mul(seed,seed)
    q=seed
    for _ in range(k): q=mul(q,pair); pair=step(pair)
    return q

def literal(k,q,centralizer):
    n=3**k; ga=old.w.permutations_a(k); t=old.w.negation(k)
    group=[p for p in ga if all(p[t[i]]==t[p[i]] for i in range(n))] if centralizer else ga
    canonical=old.cc if centralizer else old.ca
    seen=set(); out=Counter(); stabs={}
    for values in product(range(q),repeat=n):
        key=canonical(values,k)
        if key in seen: continue
        seen.add(key)
        stabilizer=[p for p in group if old.moved(values,p)==values]
        size=len(stabilizer); e=0
        while size>1:
            need(size%3==0,'ternary stabilizer order'); size//=3; e+=1
        # q2 has literal selected-vertex weight; q3 uses colour-value sum.
        out[sum(values),e]+=1
        stabs[key]=stabilizer
    return dict(out),stabs

def free_numbers(q,kmax):
    out=[q]
    for _ in range(kmax): out.append((out[-1]**3-out[-1])//3)
    return out

def main():
    rows=[]
    for k in range(5):
        pa=poly_a(k,(1,1)); pc=poly_c(k,(1,1)); n=3**k; ea=(n-1)//2; ec=(n-2*k-1)//4
        need(3**ea*eval_stabilizer(pa,Fraction(1,3))==2**n,'A labelled mass')
        need(3**ec*eval_stabilizer(pc,Fraction(1,3))==2**n,'C labelled mass')
        fa={a:v for (a,e),v in pa.items() if e==0}
        need(fa=={(n-1)//2:1,(n+1)//2:1},'inherited two rigid palettes')
        fc={a:v for (a,e),v in pc.items() if e==0}
        minimum=0 if k<=1 else (3**(k-1)-1)//2
        need(min(fc)==minimum,'minimum regular C subset')
        need(fc[minimum]==(1 if k==0 else 2**(k-1)),'minimum regular-subset orbit multiplicity')
        free=free_numbers(4,max(k-1,0)); predicted=2
        for j in range(k): predicted*=free[j]
        need(sum(fc.values())==predicted,'C free-orbit recurrence')
        if k<=2:
            for centralizer,poly in [(False,pa),(True,pc)]:
                actual,_=literal(k,2,centralizer)
                need(actual==poly,'literal binary stabilizers')
        rows.append({'depth':k,'A_stabilizer_polynomial':output_poly(pa),'C_stabilizer_polynomial':output_poly(pc),
                     'A_orbits':sum(pa.values()),'C_orbits':sum(pc.values()),'C_regular_orbits':sum(fc.values()),
                     'minimum_C_regular_subset':minimum,'minimum_C_regular_subset_orbits':fc[minimum]})
    # A nonbinary palette independently tests the generalization beyond S68.
    for centralizer,poly in [(False,poly_a(2,(1,1,1))),(True,poly_c(2,(1,1,1)))]:
        actual,_=literal(2,3,centralizer); need(actual==poly,'literal three-color stabilizers')
    # Shared-source orbit-fibre identity: sum over first colouring's stabilizer.
    joint=[]
    for k in range(3):
        for centralizer in [False,True]:
            _,stabs=literal(k,2,centralizer)
            total=0
            for group in stabs.values():
                total+=sum(old.w.color_poly(old.w.actual_index(group)))
            expected=sum((poly_c if centralizer else poly_a)(k,(1,2,1)).values())
            need(total==expected,'joint fibre/stabilizer recurrence')
            joint.append({'k':k,'group':'C' if centralizer else 'A','two_binary_joint_orbits':total})
    result={'status':'FINITE-EXACT checks of scoped recurrences; A binary specialization inherited from Session68',
            'rows':rows,'joint_fibres':joint,'checks':CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'rows':[{key:v for key,v in row.items() if 'polynomial' not in key} for row in rows],
                      'joint_fibres':joint,'checks':CHECKS},indent=2))

if __name__=='__main__': main()
