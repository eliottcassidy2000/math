#!/usr/bin/env python3
"""Exact portable controls, no external imports, assertions or floating point.
The finite checks validate implementation, not the all-order proof.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
import ast
import hashlib
import json

CHECKS = Counter()
def check(condition, tag):
    CHECKS[tag] += 1
    if not condition:
        raise RuntimeError('FAILED: '+tag)

def parity(x): return x.bit_count() & 1
def chi(x, y): return 1-2*parity(x & y)
def span(gens):
    ans={0}
    for g in gens: ans |= {x^g for x in ans}
    return sorted(ans)
def rank(a):
    if not a: return 0
    a=[list(map(Fraction,r)) for r in a]; r=0
    for j in range(len(a[0])):
        i=next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None: continue
        a[r],a[i]=a[i],a[r]
        pivot=a[r][j]; a[r]=[v/pivot for v in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:
                t=a[i][j]; a[i]=[v-t*w for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r

def matmul(a,b):
    return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]

def transform(x, mapping):
    y=0
    while x:
        b=x & -x; y ^= mapping[b.bit_length()-1]; x ^= b
    return y

class Model:
    def __init__(self,n,kind):
        self.n=n; self.kind=kind
        self.edges=list(combinations(range(n),2)); self.m=len(self.edges)
        index={e:i for i,e in enumerate(self.edges)}; self.index=index
        self.generators=[]
        for i in range(n-1):
            p=list(range(n)); p[i],p[i+1]=p[i+1],p[i]
            self.generators.append(self.permdata(p))
        cuts=[sum(1<<i for i,e in enumerate(self.edges) if v in e) for v in range(n)]
        if kind=='zero': U=[0]
        elif kind=='cut': U=span(cuts)
        elif kind=='reversal': U=span([(1<<self.m)-1])
        elif kind=='full': U=list(range(1<<self.m))
        else: raise ValueError(kind)
        self.U=U
        self.q=[min(x^u for u in U) for x in range(1<<self.m)]
        self.states=sorted(set(self.q))
        self.dual=[x for x in range(1<<self.m) if all(not parity(x&u) for u in U)]
        self.aorbs,self.aclass=self.orbits(self.states,True)
        self.dorbs,self.dclass=self.orbits(self.states,False)
        self.good=[]; self.bad=[]
        todo=set(self.dual)
        while todo:
            x=min(todo); coeff={x:1}; stack=[x]; consistent=True
            while stack:
                y=stack.pop()
                for mapping,inv in self.generators:
                    z=transform(y,mapping); s=coeff[y]*chi(inv,y)
                    if z in coeff:
                        if coeff[z]!=s: consistent=False
                    else: coeff[z]=s; stack.append(z)
            todo.difference_update(coeff)
            (self.good if consistent else self.bad).append(coeff)
        self.reps=[min(o) for o in self.aorbs]
        self.basis=[[sum(s*chi(x,a) for x,s in orb.items()) for orb in self.good] for a in self.reps]
        check(len(self.good)==len(self.aorbs),'orbit_equinumerosity')
        check(rank(self.basis)==len(self.aorbs),'signed_basis_full_rank')
        check(len(self.dual)*len(U)==1<<self.m,'dual_quotient_cardinality')
        for u in U:
            check(all(transform(u,g) in U for g,_ in self.generators),'submodule')

    def permdata(self,p):
        mapping=[]; inv=0
        for i,(u,v) in enumerate(self.edges):
            a,b=p[u],p[v]
            if a>b: inv |= 1<<i
            mapping.append(1<<self.index[tuple(sorted((a,b)))])
        return mapping,inv

    def action(self,a,data,affine):
        mapping,inv=data
        y=transform(a,mapping)
        if affine: y ^= transform(inv,mapping)
        return self.q[y]

    def orbits(self,states,affine):
        todo=set(states); orbs=[]; cls={}
        while todo:
            start=min(todo); orb={start}; stack=[start]
            while stack:
                a=stack.pop()
                for g in self.generators:
                    b=self.action(a,g,affine)
                    if b not in orb: orb.add(b); stack.append(b)
            for a in orb: cls[a]=len(orbs)
            orbs.append(sorted(orb)); todo.difference_update(orb)
        return orbs,cls

    def operator(self,kernel):
        A=[[0]*len(self.aorbs) for _ in self.aorbs]
        for i,a in enumerate(self.reps):
            for d,k in kernel.items(): A[i][self.aclass[self.q[a^d]]] += k
        return A

    def verify_operators(self):
        lambdas=[]; matrices=[]
        for orb in self.dorbs:
            A=self.operator(dict.fromkeys(orb,1)); matrices.append(A)
            lam=[sum(chi(next(iter(g)),d) for d in orb) for g in self.good]; lambdas.append(lam)
            AB=matmul(A,self.basis)
            check(AB==[[v*l for v,l in zip(row,lam)] for row in self.basis],'full_kernel_diagonalization')
            check(all(sum(r)==len(orb) for r in A),'multiplicity_row_sum')
            for i in range(len(A)):
                for j in range(len(A)):
                    check(len(self.aorbs[i])*A[i][j]==len(self.aorbs[j])*A[j][i],'detailed_balance')
        check(rank(lambdas)==len(self.good),'full_diagonal_algebra')
        return matrices

    def verify_arc(self):
        kernel=Counter(self.q[1<<i] for i in range(self.m))
        A=self.operator(kernel)
        vals=[self.m-2*next(iter(g)).bit_count() for g in self.good]
        check(matmul(A,self.basis)==[[v*l for v,l in zip(r,vals)] for r in self.basis],'arc_spectrum')
        spectrum=Counter(self.m-v for v in vals)
        if self.kind=='cut':
            if self.n>=4:
                check(min(v for v in spectrum if v)>0,'switching_positive_spectrum')
                check(min(v for v in spectrum if v)==8 and spectrum[8]==1,'switching_gap_eight_simple')
            else: check(len(spectrum)==1 and spectrum[0]==1,'switching_small_boundary')
        return A,dict(sorted(spectrum.items()))

    def independent_even_test(self):
        if self.kind!='zero' or self.n>5: return
        ps=[self.permdata(p) for p in permutations(range(self.n))]
        for orb in self.good+self.bad:
            x=min(orb)
            direct=all(chi(inv,x)==1 for mapping,inv in ps if transform(x,mapping)==x)
            check(direct==(orb in self.good),'independent_all_permutations_sign')

def coarsening(source,target):
    check(source.n==target.n and set(source.U)<=set(target.U),'coarsening_scope')
    fibre=sorted({source.q[u] for u in target.U})
    targetdual=set(target.dual)
    for orb in source.good:
        representative=next(iter(orb)); survives=representative in targetdual
        for a in source.reps:
            lhs=sum(sum(s*chi(x,a^u) for x,s in orb.items()) for u in fibre)
            rhs=len(fibre)*sum(s*chi(x,a) for x,s in orb.items()) if survives else 0
            check(lhs==rhs,'coarsening_signed_mode_selection')
    # All invariant increment-orbit kernels, literal counts before Fourier.
    for orb in source.dorbs:
        pushed=Counter(target.q[d] for d in orb)
        for a in source.states:
            observed=Counter(target.q[a^d] for d in orb)
            expected=Counter()
            for d,k in pushed.items(): expected[target.q[a^d]]+=k
            check(observed==expected,'coarsening_pushed_kernel')


def main():
    rows=[]; models={}; hostiles={}
    for n in range(6):
        for kind in ['zero','cut','reversal','full']:
            model=Model(n,kind); models[n,kind]=model
            matrices=model.verify_operators(); A,spec=model.verify_arc(); model.independent_even_test()
            rows.append({'n':n,'subspace':kind,'labelled_quotient_size':len(model.states),'affine_types':len(model.aorbs),'admissible_dual_types':len(model.good),'all_dual_types':len(model.good)+len(model.bad),'difference_kernel_types':len(model.dorbs),'arc_laplacian_multiplicities':spec})
            if n==4 and kind=='zero':
                # Sums of complete Hamming layers commute before support erasure.
                radial=[]
                for k in range(model.m+1):
                    radial.append(model.operator({d:1 for d in model.states if d.bit_count()==k}))
                for i in range(len(radial)):
                    for j in range(i+1,len(radial)):
                        check(matmul(radial[i],radial[j])==matmul(radial[j],radial[i]),'radial_commutation')
                        b=[[int(v>0 and a!=c) for c,v in enumerate(r)] for a,r in enumerate(radial[i])]
                        c=[[int(v>0 and a!=d) for d,v in enumerate(r)] for a,r in enumerate(radial[j])]
                        bc,cb=matmul(b,c),matmul(c,b)
                        if bc!=cb and 'boolean_noncommutation' not in hostiles:
                            witness=next((a,d,bc[a][d],cb[a][d]) for a in range(len(bc)) for d in range(len(bc)) if bc[a][d]!=cb[a][d])
                            hostiles['boolean_noncommutation']={'n':n,'layers':[i,j],'type_representatives':model.reps,'B_i':b,'B_j':c,'witness_row_col_products':witness}
        coarsening(models[n,'zero'],models[n,'cut'])
        coarsening(models[n,'zero'],models[n,'reversal'])
        coarsening(models[n,'cut'],models[n,'full'])
    model=Model(6,'cut'); model.verify_operators(); A,spec=model.verify_arc()
    rows.append({'n':6,'subspace':'cut','labelled_quotient_size':len(model.states),'affine_types':len(model.aorbs),'admissible_dual_types':len(model.good),'all_dual_types':len(model.good)+len(model.bad),'difference_kernel_types':len(model.dorbs),'arc_laplacian_multiplicities':spec,'weighted_arc_operator':A,'type_representatives':model.reps})
    # Real hostile: same total count is not the same graph predicate.
    m=models[3,'zero']; p3=(1<<m.index[0,1])|(1<<m.index[0,2]); c3=7
    goodmasks=set().union(*(set(o) for o in m.good))
    check(p3 in goodmasks and c3 not in goodmasks,'sign_even_vs_degree_even')
    check(len(models[4,'zero'].aorbs)==4 and len(models[4,'zero'].dorbs)==11,'omitted_cocycle_hostile')
    # Distinct invariant kernels may become equal after orbit descent.
    a=m.operator({d:1 for d in m.states if d.bit_count()==1})
    b=m.operator({d:1 for d in m.states if d.bit_count()==2})
    check(a==b,'distinct_kernels_equal_on_affine_invariants')
    hostiles['kernel_invisibility_n3']={'one_edge':a,'two_edge':b}
    # Source-dependent directed-triangle-reversal is not translation-invariant.
    def guarded(x): return x^7 if x in (2,5) else x
    check(guarded(0)^2 != guarded(0^2),'source_dependent_translation_failure')
    hostiles['source_dependent_guard']={'translation':2,'F(0)^2':guarded(0)^2,'F(0^2)':guarded(2),'warning':'On the two isomorphism types this particular map is identity; hostile is explicitly LABELLED.'}
    check('boolean_noncommutation' in hostiles,'boolean_noncommutation_detected')
    # Enforce normal/-O reliability and exact arithmetic in this verifier.
    tree=ast.parse(Path(__file__).read_text())
    check(not any(isinstance(n,ast.Assert) for n in ast.walk(tree)),'no_asserts')
    check(not any(isinstance(n,ast.Constant) and isinstance(n.value,float) for n in ast.walk(tree)),'no_floats')
    print(json.dumps({'status':'PASS','checks':dict(sorted(CHECKS.items())),'total_checks':sum(CHECKS.values()),'models':rows,'hostiles':hostiles},sort_keys=True,indent=2))

if __name__=='__main__': main()
