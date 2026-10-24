#!/usr/bin/env python3
"""Independent finite audit of Session47's modular carrier-transfer claims."""
from itertools import product
from math import gcd
from fractions import Fraction
from random import Random
import json


def mv(M, x, m):
    return tuple(sum(a*b for a,b in zip(row,x)) % m for row in M)

def tr(M):
    return list(map(list, zip(*M)))

def kernel(M,m):
    return {x for x in product(range(m),repeat=len(M[0])) if all(t == 0 for t in mv(M,x,m))}

def image(M,m):
    return {mv(M,x,m) for x in product(range(m),repeat=len(M[0]))}

def block(B,C,D):
    return [list(b)+list(c) for b,c in zip(B,C)] + [list(c)+list(d) for c,d in zip(tr(C),D)]

def test_equivalences(A,B,C,D,m):
    n=len(A); k=len(D)
    Ap=block(B,C,D); E=[list(c) for c in C]+[list(d) for d in D]
    F=[list(b) for b in B]+[list(c) for c in tr(C)]
    old=kernel(A,m); new=kernel(Ap,m)
    extension=all(any(z[:n] == x for z in new) for x in old)
    primal=all(mv(F,x,m) in image(E,m) for x in old)
    dual=all(mv(tr(F),z,m) in image(tr(A),m) for z in kernel(tr(E),m))
    old_image=image(A,m)
    current_dual=all(d in old_image for d in product(range(m),repeat=n) if d+(0,)*k in image(Ap,m))
    if len({extension,primal,dual,current_dual}) != 1:
        raise ValueError(('criterion disagreement',A,B,C,D,m,extension,primal,dual,current_dual))
    if gcd(det(D),m)==1:
        # Avoid a rational inverse: old x extends uniquely at the new rows;
        # test its remaining old-row residual directly.
        schur=True
        for x in old:
            ys=[y for y in product(range(m),repeat=k) if mv(D,y,m)==tuple(-t % m for t in mv(tr(C),x,m))]
            if len(ys)!=1: raise ValueError('unit D is not invertible')
            if any((a+b)%m for a,b in zip(mv(B,x,m),mv(C,ys[0],m))): schur=False
        if extension!=schur: raise ValueError('Schur disagreement')
    return extension

def det(M):
    if len(M)==1:return M[0][0]
    if len(M)==2:return M[0][0]*M[1][1]-M[0][1]*M[1][0]
    return sum((-1)**j*M[0][j]*det([row[:j]+row[j+1:] for row in M[1:]]) for j in range(len(M)))

def v2(x):
    a=0
    while x%2==0:a+=1;x//=2
    return a

def diag_extension(ds,B,C,D,m):
    E=C+D; F=B+tr(C); im=image(E,m)
    return all(mv(F,tuple((m//gcd(d,m))*(i==j) for j in range(len(ds))),m) in im for i,d in enumerate(ds))

def diagonal_section(ds,B,C,D,m):
    # An ambient P=[I;T] is determined independently column by column.
    E=C+D; F=B+tr(C); k=len(D)
    for i,d in enumerate(ds):
        h=m//gcd(d,m)
        if not any(all(h*(F[j][i]+v)%m==0 for j,v in enumerate(mv(E,t,m)))
                   for t in product(range(m),repeat=k)):
            return False
    return True

def require(condition, context):
    if not condition:
        raise ValueError(context)

def main():
    rng=Random(470024)
    result={"status":"PASS", "checks":{}}
    # Full kernels and images, rather than just one algebraic formulation.
    count=0
    for m in (2,3,4,6):
        for a,b,c,d in product(range(4),repeat=4):
            test_equivalences([[a]],[[b]],[[c]],[[d]],m); count+=1
    result['checks']['full_scalar_modular_equivalence_cases']=count
    count=0
    for _ in range(120):
        a,b,c=[rng.randrange(-2,4) for _ in range(3)]
        A=[[a,b],[b,c]]
        x,y,z=[rng.randrange(-2,4) for _ in range(3)]
        B=[[x,y],[y,z]]; C=[[rng.randrange(-2,4)] for _ in range(2)]; D=[[rng.randrange(-2,4)]]
        for m in (2,4):test_equivalences(A,B,C,D,m);count+=1
    result['checks']['full_nondiagonal_modular_equivalence_cases']=count
    # Independent hostile enumeration.
    A=[[2]]; B=[[4]]; C=[[-1,-1]]; D=[[2,0],[0,2]]
    Ap=block(B,C,D)
    ext4=sorted(z for z in kernel(Ap,4) if z[0]==2)
    require(ext4==[(2,1,3),(2,3,1)], ("hostile extensions",ext4))
    require(all(any(2*t%4 for t in z) for z in ext4), "hostile has an order-two lift")
    require(test_equivalences(A,B,C,D,4), "hostile mod4 failure")
    require(not test_equivalences(A,B,C,D,2), "hostile mod2 success")
    require(not diagonal_section([2],B,C,D,4), "hostile has a section")
    result['hostile']={'matrix':Ap,'old_mod4_nonzero_extensions':ext4,'mod4_extension':True,'mod4_additive_section':False,'mod2_extension':False,'determinant':det(Ap)}
    ds=[2,4];B=[[4,0],[0,4]];C=[[-1,-1],[0,0]];D=[[2,0],[0,2]]
    require(test_equivalences([[2,0],[0,4]],B,C,D,4), "true-top hostile fails at top")
    require(not test_equivalences([[2,0],[0,4]],B,C,D,2), "true-top hostile passes below top")
    require(not diagonal_section(ds,B,C,D,4), "true-top hostile has a section")
    result['saturation_top_hostile']={'old_smith':[2,4],'true_saturation_level':2,'mod4_extension':True,'mod2_extension':False,'matrix':block(B,C,D)}
    count=0; positive=0; higher_pass_lower_fail=0
    for _ in range(1800):
        ds=[rng.choice((1,2,3,4,6,8,12)) for _ in range(rng.choice((1,2)))]
        n=len(ds); k=1
        B=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i,n):B[i][j]=B[j][i]=rng.randrange(-4,9)
        C=[[rng.randrange(-4,5)] for _ in range(n)];D=[[rng.randrange(-4,9)]]
        s=max(map(v2,ds)); levels=sorted({v2(d) for d in ds if v2(d)>0})
        vals={r:diag_extension(ds,B,C,D,2**r) for r in range(1,s+3)}
        finite=all(vals[r] for r in range(1,s+1))
        smith=all(vals[r] for r in levels)
        section=diagonal_section(ds,B,C,D,2**s)
        require(finite==smith==section, ("finite equivalence",ds,B,C,D,vals,section))
        require(not finite or all(vals.values()), ("upward saturation",ds,B,C,D,vals))
        if s:require(not vals[s] or all(vals[r] for r in range(s,s+3)), "top-to-higher failure")
        else:require(all(vals.values()), "odd old determinant nontrivial kernel")
        if finite:positive+=1
        if s and vals[s] and not finite:higher_pass_lower_fail+=1
        count+=1
    result['checks']['smith_level_and_uniform_section_random_cases']=count
    result['checks']['smith_positive_cases']=positive
    result['checks']['true_top_pass_lower_fail_random_cases']=higher_pass_lower_fail
    # Explicit kernel scaling with a nontrivial unimodular change of basis.
    count=0
    for ds in ([2,4],[3,8],[1,12],[6,8]):
        for t in (-2,-1,1,2):
            # A=W diag(ds) W^T, W=[[1,t],[0,1]], right Smith basis V=W^-T.
            A=[[ds[0]+t*t*ds[1],t*ds[1]],[t*ds[1],ds[1]]]
            s=max(map(v2,ds)); Hs=kernel(A,2**s)
            for r in (s,s+1,s+2):
                scaled={tuple(2**(r-s)*z for z in x) for x in Hs}
                require(scaled==kernel(A,2**r), ("kernel scaling",A,s,r))
                count+=1
    result['checks']['nontrivial_basis_kernel_scaling_cases']=count
    # General readouts: one integer Z preserves the equations on harmonic
    # inputs, without asserting an ambient right inverse of the readout.
    count=0; positive=0
    for _ in range(400):
        ds=[rng.choice((1,2,3,4,6,8,12)) for _ in range(rng.choice((1,2)))]
        n=len(ds); p=rng.choice((1,2)); q=3
        M=[[rng.randrange(-3,4) for _ in range(p)] for _ in range(q)]
        J=[[rng.randrange(-3,4) for _ in range(n)] for _ in range(q)]
        s=max(map(v2,ds)); levels=sorted({v2(d) for d in ds if v2(d)>0})
        def covers(r):
            m=2**r; im=image(M,m)
            return all(mv(J,tuple((m//gcd(d,m))*(i==j) for j in range(n)),m) in im
                       for i,d in enumerate(ds))
        vals={r:covers(r) for r in range(1,s+3)}
        finite=all(vals[r] for r in levels)
        require(finite==all(vals.values()), ("general-readout level criterion",ds,M,J,vals))
        if finite:
            cols=[]
            for i,d in enumerate(ds):
                m=2**v2(d)
                z=next(z for z in product(range(m),repeat=p) if mv(M,z,m)==tuple(row[i]%m for row in J))
                cols.append(z)
            Z=tr(cols)
            defect=[[sum(M[i][h]*Z[h][j] for h in range(p))-J[i][j] for j in range(n)] for i in range(q)]
            require(all(Fraction(defect[i][j],ds[j]).denominator%2==1 for i in range(q) for j in range(n)),
                    ("general-readout odd denominators",ds,M,J,Z))
            for r in range(1,s+3):
                m=2**r
                require(all(not any(mv(defect,tuple((m//gcd(d,m))*(i==j) for j in range(n)),m))
                            for i,d in enumerate(ds)), ("general-readout common lift",ds,M,J,Z,r))
            positive+=1
        count+=1
    result['checks']['general_readout_common_lift_cases']=count
    result['checks']['general_readout_positive_cases']=positive
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
