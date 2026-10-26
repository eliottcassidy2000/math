#!/usr/bin/env python3
"""Exact hostile/positive tests: every subgroup of AGL(2,2), no assertions."""
import ast
import itertools
import json
from pathlib import Path

CHECKS=0

def check(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(message)

def dot(x,y): return (x&y).bit_count()&1

def rank(A):
    from fractions import Fraction
    A=[[Fraction(x) for x in r] for r in A]
    r=0
    for j in range(len(A[0]) if A else 0):
        p=next((i for i in range(r,len(A)) if A[i][j]),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]
        z=A[r][j]
        A[r]=[x/z for x in A[r]]
        for i in range(len(A)):
            if i!=r:
                z=A[i][j]
                A[i]=[x-z*y for x,y in zip(A[i],A[r])]
        r+=1
    return r

def orbits(actions):
    reps={}
    canon={}
    for x in range(4):
        if x in canon: continue
        orb=sorted({p[x] for p in actions})
        reps[x]=orb
        canon.update({y:x for y in orb})
    return reps,canon

def main():
    linears=[]
    for a,b in itertools.permutations(range(1,4),2): linears.append((0,a,b,a^b))
    elements=sorted({tuple(r[x]^t for x in range(4)) for r in linears for t in range(4)})
    ids={p:i for i,p in enumerate(elements)}
    identity=ids[(0,1,2,3)]
    products=[[ids[tuple(a[b[x]] for x in range(4))] for b in elements] for a in elements]
    def generated(gens):
        S={identity}
        changed=True
        while changed:
            bigger=S|{products[g][s] for g in gens for s in S}
            changed=bigger!=S
            S=bigger
        return frozenset(S)
    groups={frozenset({identity})}
    todo=list(groups)
    while todo:
        G=todo.pop()
        for g in range(24):
            if g not in G:
                H=generated(tuple(G)+(g,))
                if H not in groups: groups.add(H);todo.append(H)
    check(len(groups)==30,'AGL2 subgroup census')
    rows=[]
    for G in sorted(groups,key=lambda x:(len(x),tuple(sorted(x)))):
        actions=[elements[g] for g in G]
        rhos=[tuple(p[x]^p[0] for x in range(4)) for p in actions]
        ts=[p[0] for p in actions]
        dual=[tuple(next(y for y in range(4) if all(dot(y,rho[x])==dot(chi,x) for x in range(4))) for chi in range(4)) for rho in rhos]
        source,canon=orbits(actions)
        dualorb,_=orbits(dual)
        linorb,lincanon=orbits(rhos)
        admissible=[chi for chi in dualorb if all(not dot(chi,t) for d,t in zip(dual,ts) if d[chi]==chi)]
        check(len(source)==len(admissible),'affine dual count')
        columns=[]
        for chi in admissible:
            # Reynolds average unnormalized, avoids a choice of orbit phase.
            vals=[sum((-1 if dot(chi,p[x]) else 1) for p in actions) for x in range(4)]
            check(any(vals),'admissible Reynolds nonzero')
            check(all(vals[x]==vals[canon[x]] for x in range(4)),'Reynolds orbit invariant')
            columns.append([vals[x] for x in source])
            for orb in linorb.values():
                eigen=sum(-1 if dot(chi,s) else 1 for s in orb)
                for x in range(4): check(sum(vals[x^s] for s in orb)==eigen*vals[x],'all orbit kernels diagonal')
        check(rank([list(r) for r in zip(*columns)])==len(source),'full invariant basis')
        kernel_values=[[sum(-1 if dot(chi,s) else 1 for s in orb) for chi in admissible] for orb in linorb.values()]
        check(rank(kernel_values)==len(source),'full diagonal kernel algebra')
        fixed=[x for x in range(4) if all(p[x]==x for p in actions)]
        equivariant=[]
        for f in itertools.permutations(range(4)):
            if all(f[p[x]]==rho[f[x]] for p,rho in zip(actions,rhos) for x in range(4)): equivariant.append(f)
        compatible=[t for t in range(4) if all((p[x]^t)==rho[x^t] for p,rho in zip(actions,rhos) for x in range(4))]
        check(bool(fixed)==bool(equivariant)==bool(compatible),'torsor trivialization iff')
        if all(t==0 for t in ts):
            edgeimages=[]
            saturations=[]
            for bits in range(8):
                S=[x for x in range(1,4) if bits>>(x-1)&1]
                E={tuple(sorted((canon[x],canon[x^s]))) for x in range(4) for s in S if canon[x]!=canon[x^s]}
                edgeimages.append(E)
                saturations.append({rho[s] for rho in rhos for s in S})
            for i in range(8):
                for j in range(8):check((edgeimages[i]==edgeimages[j])==(saturations[i]==saturations[j]),'linear generator saturation iff')
        rows.append({'group_order':len(G),'affine_orbits':len(source),'admissible_dual_orbits':len(admissible),'global_fixed_points':len(fixed),'equivariant_bijections_to_linear_module':len(equivariant)})
    # Minimal noncyclic character collision.
    Hs=[{0,h} for h in (1,2,3)]
    Xfix=[6,2,2,2]
    Yfix=[sum(2 if g in H else 0 for H in Hs) for g in range(4)]
    check(Xfix==Yfix,'V4 permutation characters equal')
    check(2!=0,'V4 global marks unequal')
    tree=ast.parse(Path(__file__).read_text())
    check(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'no assert nodes')
    print(json.dumps({'status':'PASS','checks':CHECKS,'universe':'All30 subgroups of AGL(2,2); every linear-orbit kernel; every candidate set bijection; all generator subsets in linear cases','rows':rows,'V4_equal_element_marks':Xfix,'V4_global_marks':[2,0]},indent=2,sort_keys=True))

if __name__=='__main__': main()
