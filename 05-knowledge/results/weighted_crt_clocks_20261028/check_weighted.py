"""Cycle-index recurrences, weighted attachments, and A000568 operator recurrence."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import gcd, factorial
from pathlib import Path
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

# A monomial is the sorted list of its cycle lengths, including multiplicity.
def mul(a,b):
    c=defaultdict(Q)
    for u,x in a.items():
        for v,y in b.items(): c[tuple(sorted(u+v))]+=x*y
    return dict(c)

def step_a(z):
    out={p:c/3 for p,c in mul(mul(z,z),z).items()}
    for p,c in z.items():
        key=tuple(3*j for j in p); out[key]=out.get(key,Q(0))+2*c/3
    return out

def step_c(c,a): return mul(c,{tuple(sorted(p+p)):v for p,v in a.items()})

def color_poly(z):
    out=defaultdict(Q)
    for p,c in z.items():
        poly={0:1}
        for j in p:
            new=Counter(poly)
            for d,v in poly.items(): new[d+j]+=v
            poly=dict(new)
        for d,v in poly.items(): out[d]+=c*v
    need(all(v.denominator==1 for v in out.values()),"integral orbit polynomial")
    return [int(out[i]) for i in range(max(out)+1)]

def tournament_eval(z):
    total=Q(0)
    for p,c in z.items():
        e=sum((j-1)//2 for j in p)+sum(gcd(j,l) for i,j in enumerate(p) for l in p[i+1:])
        total+=c*2**e
    need(total.denominator==1,"integral tournament evaluation")
    return int(total)

def inv(p):
    out=[0]*len(p)
    for i,x in enumerate(p): out[x]=i
    return tuple(out)

@lru_cache(None)
def permutations_a(k):
    if not k: return [(0,)]
    old=permutations_a(k-1); d=3**(k-1); out=[]
    for shift in range(3):
        for ps in product(old,repeat=3):
            out.append(tuple(((b+shift)%3)*d+ps[b][i] for b in range(3) for i in range(d)))
    return out

def negation(k):
    out=[]
    for x in range(3**k):
        y=0; d=1
        for _ in range(k): y+=((-x%3)%3)*d; x//=3; d*=3
        out.append(y)
    return tuple(out)

def cycles(p):
    seen=set(); out=[]
    for x in range(len(p)):
        if x in seen: continue
        y=x; d=0
        while y not in seen: seen.add(y); d+=1; y=p[y]
        out.append(d)
    return tuple(sorted(out))

def actual_index(group):
    count=Counter(cycles(p) for p in group)
    return {p:Q(v,len(group)) for p,v in count.items()}

def literal_mask_poly(group):
    n=len(group[0]); seen=set(); out=[0]*(n+1)
    for mask in range(1<<n):
        if mask in seen: continue
        orbit={sum(1<<p[i] for i in range(n) if mask>>i&1) for p in group}
        seen.update(orbit); out[mask.bit_count()]+=1
    return out

def operator_counts(limit):
    js=tuple(range(1,limit+1,2)); visits=0
    @lru_cache(None)
    def rec(n,e):
        nonlocal visits
        visits+=1
        if n==0: return Q(1)
        return sum((2**((j-1)//2+e[t])*rec(n-j,tuple(e[s]+gcd(j,l) for s,l in enumerate(js)))
                    for t,j in enumerate(js) if j<=n),Q(0))/n
    values=[rec(n,(0,)*len(js)) for n in range(limit+1)]
    need(all(x.denominator==1 for x in values),"operator integrality")
    inherited=module('old_burnside',ROOT/'04-computation/a000568_exact_burnside.py')
    for n,x in enumerate(values): need(x==inherited.a000568(n),"independent odd-partition sum")
    strong=[0]*(limit+1)
    for n in range(1,limit+1): strong[n]=int(values[n])-sum(strong[i]*int(values[n-i]) for i in range(1,n))
    return {"A000568":[int(x) for x in values],"strong_tournament_counts":strong,"operator_cache_states":visits}

def main():
    a=c={(1,):Q(1)}; rows=[]
    for k in range(4):
        need(sum(a.values())==sum(c.values())==1,"cycle-index normalization")
        moments_a=[sum(v*p.count(1)**r for p,v in a.items()) for r in range(4)]
        moments_c=[sum(v*p.count(1)**r for p,v in c.items()) for r in range(4)]
        need(moments_a[1:]==[1,2*k+1,6*k*k+2*k+1],"A moments")
        need(moments_c[1:]==[2*k+1,8*k*k+1,48*k**3-40*k*k+18*k+1],"C moments")
        pa,pc=color_poly(a),color_poly(c)
        if k<=2:
            ga=permutations_a(k); t=negation(k)
            gc=[p for p in ga if all(p[t[i]]==t[p[i]] for i in range(len(p)))]
            need(a==actual_index(ga) and c==actual_index(gc),"literal centralizer cycle indices")
            need(pa==literal_mask_poly(ga) and pc==literal_mask_poly(gc),"literal subset orbits")
        rows.append({"k":k,"vertices":3**k,"A_cycle_types":len(a),"C_cycle_types":len(c),
                     "A_binary_weight_polynomial":pa,"C_binary_weight_polynomial":pc,
                     "A_binary_orbits":sum(pa),"C_binary_orbits":sum(pc),
                     "A_tournament_orbits":tournament_eval(a),"C_tournament_orbits":tournament_eval(c)})
        a,c=step_a(a),step_c(c,a)
    oracle=module('root_oracle',HERE.parent/'relation_tower_frobenius_ports_20261028/check_root_switch.py')
    hostile=oracle.hostile(); polys=[]
    for name in ['rows_AAA','rows_BAA']:
        group=oracle.isomorphisms(hostile[name],hostile[name],True)
        poly=color_poly(actual_index(group))
        need(poly==literal_mask_poly(group),"81/27 literal masks")
        polys.append(poly)
    need([sum(p) for p in polys]==[176,512],"81/27 attachment distinction")
    hostile['binary_weight_polynomials']=polys
    result={"status":"PROVED recurrences; FINITE-EXACT checks","towers":rows,
            "a000568_operator":operator_counts(20),"hostile_81_27":hostile,"checks":CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
