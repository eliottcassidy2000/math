"""Top-weight zero-heavy boundary; exact inherited full compiler plus local DP."""
from collections import Counter
from fractions import Fraction
from itertools import product, combinations
from math import comb
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)
def load(name,p):
    s=spec_from_file_location(name,p);m=module_from_spec(s);s.loader.exec_module(m);return m
P=Path(__file__).resolve().parents[1]
old=load('inherited_sparse',P/'defect_boundary_20261020/sparse_count.py')

def matchings(xs):
    if not xs:yield ();return
    a=xs[0]
    for b in xs[1:]:
        for rest in matchings(tuple(x for x in xs[1:] if x!=b)):yield ((a,b),)+rest
def compositions(total,parts):
    if parts==1:yield(total,);return
    for x in range(total+1):
        for q in compositions(total-x,parts-1):yield(x,)+q
def model(k,matching,gaps,strong=False):
    if strong:gaps=(0,)+gaps+(0,)
    slots={};weights=[]
    for i,g in enumerate(gaps):
        if g:weights.append(g)
        if i<2*k:slots[i]=len(weights);weights.append(1)
    n=len(weights);t=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for a,b in matching:
        a,b=slots[a],slots[b];t[a]^=1<<b;t[b]^=1<<a
    return tuple(t),tuple(weights)
def strong(t):
    for i in range(len(t)):
        seen=1<<i
        while True:
            nxt=seen
            for j in range(len(t)):
                if seen>>j&1:nxt|=t[j]
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(t))-1:return False
    return True
def indecomp(m,k):return all(any(a<c<=b for a,b in m) for c in range(1,2*k))
def fixed_cycle(n,mask):
    ds=[i for i in range(1,(n+1)//2)]
    steps=[d if mask>>(d-1)&1 else n-d for d in ds]
    return tuple(sum(1<<((i+d)%n) for d in steps) for i in range(n))

def main():
    complete=[]
    for k in (1,2):
        states,_=old.template_states(k); mass=Counter();smass=Counter()
        for q,code,w in states:
            if max(w)>k:continue
            t=old.E.decode(q,code);mass[sum(w)]+=1
            if q>1 and strong(t):smass[sum(w)]+=1
        need(max(mass)==2*k*k+3*k,'sharp zero-heavy mass')
        need(max(smass)==2*k*k+k,'sharp strong mass')
        complete.append({'k':k,'all_zero_heavy_by_mass':dict(sorted(mass.items())),
                         'strong_zero_heavy_by_mass':dict(sorted(smass.items()))})
    generated=[];dp_checks=0
    for k in range(2,6):
        ms=list(matchings(tuple(range(2*k)))); a=sum(indecomp(m,k) for m in ms)
        for delta in range(min(k-2,2)+1):
            for isstrong in (False,True):
                g=2*k-1 if isstrong else 2*k+1; count=0
                # Full structure/clone/strong checks for every generated presentation.
                for m in ms:
                    if isstrong and not indecomp(m,k):continue
                    for ds in compositions(delta,g):
                        t,w=model(k,m,tuple(k-d for d in ds),isstrong)
                        need(len(old.G.clone_chains(t))==len(t),'no hidden clone merger')
                        if isstrong:need(strong(t),'indecomposable matching is strong')
                        if k<=2 or (k==3 and count<15):
                            need(old.fas(t,w)==k,'independent weighted optimum');dp_checks+=1
                        count+=1
                expected=(a if isstrong else len(ms))*comb(g+delta-1,delta)
                need(count==expected,'top layer count')
                generated.append({'k':k,'delta':delta,'strong':isstrong,'count':count})
    # Same prime automorphism cycle and score data need different feedback budgets.
    cyc=Counter();witness={}
    for mask in range(8):
        t=fixed_cycle(7,mask);f=old.fas(t);cyc[f]+=1;witness.setdefault(f,t)
        need(len(old.G.clone_chains(t))==7,'single odd cycle pair free')
    need(len(cyc)>1,'feedback is not fixed-permutation count')
    print(json.dumps({'status':'FINITE-EXACT controls of top-mass theorem',
        'complete_budget_compiler':complete,'matching_layers':generated,
        'weighted_DP_controls':dp_checks,'single_seven_cycle_feedback':dict(cyc),
        'cycle_witnesses':witness},indent=2))
if __name__=='__main__':main()
