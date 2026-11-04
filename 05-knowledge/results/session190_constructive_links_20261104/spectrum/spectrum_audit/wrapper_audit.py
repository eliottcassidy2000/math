#!/usr/bin/env python3
"""Independent literal wrapper controls and negative boundary controls."""
from independent_audit import profile, valid, strong, A, ROOT_A, flip, need
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def wrap(t,v,k):
    n=len(t);w=n+k
    # Build the arc set explicitly, without using producer row formulas.
    arcs={(i,j) for i in range(n) for j in range(n) if (t[i]>>j)&1}
    z=range(n,n+k)
    arcs.update((i,j) for i in z for j in z if i<j)
    arcs.update((i,j) for i in z for j in range(n))
    arcs.update((w,j) for j in range(w) if j!=v)
    arcs.add((v,w))
    return tuple(sum(1<<j for i2,j in arcs if i2==i) for i in range(w+1))

def root_is_reachable(t,v):
    n=len(t)
    for i in range(n):
        seen={i};todo=[i]
        while todo:
            a=todo.pop()
            for b in range(n):
                if (t[a]>>b)&1 and b not in seen:seen.add(b);todo.append(b)
        if v not in seen:return False
    return True

def main():
    states=[A]
    for edge in [(1,2),(0,1),(2,3),(3,5)]:states.append(flip(states[-1],edge))
    controls=[]
    for k in range(5):
        counts=[]
        for t in states:
            p=profile(t,ROOT_A);r=wrap(t,ROOT_A,k);valid(r);got=profile(r,ROOT_A)
            expected=p['H']+(1<<k)*p['F']-(p['s'] if k==0 else 0)
            need(got['H']==expected)
            if k>=1:need(got['strong'])
            counts.append(got['H'])
        if k>=1:need(counts[-1]==counts[0]+2 and max(counts[1:-1])<counts[0])
        controls.append(dict(k=k,counts=counts,endpoint_difference=counts[-1]-counts[0]))
    all_controls=0
    for n in range(1,5):
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        for mask in range(1<<len(pairs)):
            r=[0]*n
            for b,(i,j) in enumerate(pairs):
                if mask>>b&1:r[i]|=1<<j
                else:r[j]|=1<<i
            for v in range(n):
                p=profile(r,v)
                for k in [0,1,2]:
                    w=wrap(r,v,k);got=profile(w,v)
                    need(got['H']==p['H']+(1<<k)*p['F']-(p['s'] if k==0 else 0))
                    if k:need(got['strong']==root_is_reachable(r,v))
                    all_controls+=1
    out=dict(status='PASS',method='Explicit enumeration of simple directed paths in literal wrappers; no producer code imported.',
             selected_patch_controls=controls,small_universe='Every labeled tournament of orders 1 through 4, every root, k in {0,1,2}',
             small_universe_cases=all_controls,
             general_formula_k_ge_1='H(W_k(T,v)) = H(T) + 2^k F_T(v)',
             boundary_formula_k_eq_0='H(W_0(T,v)) = H(T) + F_T(v) - s_v(T)',
             strong_criterion_k_ge_1='Every vertex of T reaches v by a directed path in T (length zero allowed).',
             all_k_nonovershoot_differences=['12 - 23*2^k','-68*2^k','-4 - 36*2^k'],
             all_k_endpoint_counts=['29+157*2^k','31+157*2^k'],
             excluded_conclusion='No universal successor operation, arbitrary wrapped-order minimality, or global spectrum closure is established.')
    (HERE/'wrapper_audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
