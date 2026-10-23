#!/usr/bin/env python3
"""Independent exact rational tests for the continuum-comb claims.
Finite-grid rank tests are witnesses on an explicitly finite-dimensional subspace,
not a replacement for the continuum Fourier proof in AUDIT.md.
"""
from fractions import Fraction as F
from pathlib import Path
import json, itertools, math
import sympy as S
HERE=Path(__file__).resolve().parent

def check(ok,msg):
    if not ok:raise RuntimeError(msg)

def split_arc(a,b):
    a,b=F(a),F(b);length=b-a
    check(0<=length<=1,'bad arc')
    a%=1;b=a+length
    return [(a,b)] if b<=1 else [(a,F(1)),(F(0),b-1)]

def comb(b,n,phase=F(0)):
    return [piece for i in range(n) for piece in split_arc(F(i,n)+phase-F(1,2*b*n),F(i,n)+phase+F(1,2*b*n))]

def overlap(A,B):
    return sum((max(F(0),min(b,d)-max(a,c)) for a,b in A for c,d in B),F(0))

def regular_hostile(b,J):
    q=b**J
    K=[(F(4*r,4*q),F(4*r+1,4*q)) for r in range(q)]
    L=[(F(4*r+2,4*q),F(4*r+3,4*q)) for r in range(q)]
    check(overlap(K,L)==0,'hostile is not disjoint')
    number=0
    for j in range(J):
        n=b**j;D=comb(b,n)
        # Breakpoints plus one point per intervening phase cell certify equality
        # of the piecewise-linear overlap responses for this interval example.
        endpoints=[z for I in K+L for z in I];walls=[z for I in D for z in I]
        cuts=sorted({(x-y)%1 for x in endpoints for y in walls}|{F(0),F(1)})
        phases=cuts[:-1]+[(x+y)/2 for x,y in zip(cuts,cuts[1:])]
        for s in phases:
            a=overlap(K,comb(b,n,s));d=overlap(L,comb(b,n,s))
            check(a==d==F(1,4*b),'retained hostile response incorrect')
        number+=len(phases)
    nextK=overlap(K,comb(b,q));nextL=overlap(L,comb(b,q))
    check(nextK==F(1,2*b) and nextL==0,'next-comb separation failed')
    return dict(b=b,J=J,q=q,exact_phase_tests=number,next_response_K=str(nextK),next_response_L=str(nextL))

def finite_grid(b,J):
    # Grid-cell indicator functions, with M chosen to align every comb wall.
    q=b**J;M=2*q;rows=[]
    for j in range(J):
        n=b**j
        mask=[int(min((n*F(2*r+1,2*M))%1,1-(n*F(2*r+1,2*M))%1)<F(1,2*b)) for r in range(M)]
        for s in range(M):rows.append([mask[(r-s)%M] for r in range(M)])
    mat=S.Matrix(rows);rank=mat.rank();null=mat.nullspace()
    check(M-rank==M//q-1,'unexpected sampled-grid nullity')
    for x in null:
        check(sum(x)==0 and all(x[r]==x[(r+M//q)%M] for r in range(M)),'grid kernel not periodic mean zero')
    return dict(b=b,J=J,cells=M,rank=rank,kernel_dimension=M-rank,predicted=M//q-1)

def phase_hostile():
    probes=[(1,F(0)),(2,F(1,11)),(7,F(2,13)),(13,F(3,17)),(23,F(5,19))]
    points=[F(1,101),F(4,37),F(9,19)]
    walls=sorted({F(0),F(1),*points,*[x%1 for n,s in probes for I in comb(7,n,s) for x in I]})
    a,z=max(zip(walls,walls[1:]),key=lambda p:p[1]-p[0]);d=(z-a)/8
    K=[(a+d,a+2*d)];L=[(a+4*d,a+5*d)]
    values=[]
    for n,s in probes:
        u=overlap(K,comb(7,n,s));v=overlap(L,comb(7,n,s));check(u==v,'finite bank separates hostile');values.append(str(u))
    check(all(not any(a<=p<=b for a,b in K+L) for p in points),'membership sample hits hostile')
    return dict(probes=[(n,str(s)) for n,s in probes],K=[[str(x) for x in I] for I in K],L=[[str(x) for x in I] for I in L],equal_overlap_values=values)

def AP_skeleton():
    # Every strict-danger wall is rational. Safe cells are decided by midpoints;
    # the remaining weak-safe endpoints are tested individually, exactly.
    walls=sorted({F(0),F(1),*[x%1 for j in range(1,14) for I in comb(7,j) for x in I]})
    def safe(t):return all(min((j*t)%1,1-(j*t)%1)>=F(1,14) for j in range(1,14))
    check(not any(safe((a+b)/2) for a,b in zip(walls,walls[1:])),'AP has positive-length safe interval')
    isolated=[t for t in walls[:-1] if safe(t)]
    expected=[F(a,14) for a in range(1,14) if math.gcd(a,14)==1]
    check(isolated==expected,'AP endpoint set differs')
    return dict(safe_arcs=0,safe_isolated_points=list(map(str,isolated)),wall_count=len(walls)-1)

def run():
    out={}
    out['regular_hostiles']=[regular_hostile(b,J) for b,J in [(2,1),(2,2),(3,2),(4,2),(6,1),(7,1),(7,2)]]
    out['finite_grid_restrictions']=[finite_grid(b,J) for b,J in [(2,1),(2,2),(3,1),(3,2),(4,1),(7,1)]]
    out['finite_phase_hostile']=phase_hostile();out['AP_skeleton']=AP_skeleton()
    out['status']='PASS'
    (HERE/'runner_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2));print('ALL INDEPENDENT RUNNER CHECKS PASSED')
if __name__=='__main__':run()
