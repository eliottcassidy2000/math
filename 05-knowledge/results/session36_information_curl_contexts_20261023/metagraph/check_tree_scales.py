#!/usr/bin/env python3
"""Two-color spanning-tree polynomial: sharp exponent support and contraction endpoints."""
import json,itertools as it
from pathlib import Path
import sympy as sp
from check_dyadic_network import lap,need,components,tree,val2
from check_dyadic_lift import dyadic_smith
HERE=Path(__file__).resolve().parent

def contract(n,edges,parts):
    where={v:i for i,c in enumerate(parts) for v in c}
    return [(where[i],where[j],w) for i,j,w in edges if where[i]!=where[j]]

def component_trees(n,edges,parts):
    L=lap(n,edges)
    return [tree(L.extract(c,c)) for c in parts]

def run(row):
    n=row['nodes'];blue=row['blue_edges'];black=row['paired_black_edges'];cb=components(n,blue);ck=components(n,black)
    Lb=lap(n,blue);Lk=lap(n,black);lo=len(ck)-1;hi=n-len(cb);z=sp.symbols('z')
    # Degree hi is the blue graphic rank, sharper than n-1.
    values=[(t,tree(t*Lb+2*Lk)) for t in range(hi+1)]
    p=sp.Poly(sp.interpolate(values,z),z);coeff=[int(p.nth(k)) for k in range(hi+1)]
    need(all(p.nth(k)==coeff[k] for k in range(hi+1)))
    need([k for k,c in enumerate(coeff) if c]==list(range(lo,hi+1)))
    need(all(c>0 for c in coeff[lo:]))
    need(p.eval(-2)==tree(-2*Lb+2*Lk))
    need(p.eval(hi+2)==tree((hi+2)*Lb+2*Lk))
    tblue=component_trees(n,blue,cb);tblack=component_trees(n,black,ck)
    qblue=tree(lap(len(ck),contract(n,blue,ck)));qblack=tree(lap(len(cb),contract(n,black,cb)))
    low=2**(n-len(ck))*int(sp.prod(tblack))*qblue
    high=2**(len(cb)-1)*int(sp.prod(tblue))*qblack
    need(coeff[lo]==low);need(coeff[hi]==high)
    for k,c in enumerate(coeff):need(c%2**(n-1-k)==0)
    return dict(n=row['n'],minimum_blue_edges=lo,maximum_blue_edges=hi,coefficients_by_blue_count=coeff,blue_component_tree_counts=tblue,black_component_tree_counts=tblack,contract_black_blue_tree_count=qblue,contract_blue_paired_black_tree_count=qblack)

if __name__=='__main__':
 rows=json.loads((HERE/'dyadic_network_receipt.json').read_text());out=[];bluegroups=[]
 for row in rows:
    n=row['nodes'];Lb=lap(n,row['blue_edges']);vals=[]
    for c in components(n,row['blue_edges']):
        if len(c)>1:vals += [v for v in dyadic_smith(Lb.extract(c[:-1],c[:-1]).tolist()) if v]
    bluegroups.append(dict(n=row['n'],blue_two_primary_exponents=sorted(vals)))
    if row['n']<=6:
        a=run(row);out.append(a);print(a,flush=True)
 (HERE/'tree_scales_receipt.json').write_text(json.dumps(dict(tree_polynomials=out,blue_critical_groups=bluegroups),indent=2)+'\n')
 print('BLUE GROUPS',bluegroups,flush=True);print('ALL CHECKS PASSED')
