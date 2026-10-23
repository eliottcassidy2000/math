#!/usr/bin/env python3
"""Small exhaustive and randomized exact Smith/half-form controls."""
from check_metagraph_independent import binary_kernel,check
from pathlib import Path
import random,itertools,json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent

def inspect(A):
    M=S.Matrix(A)
    if M.det()==0:return None
    rank,basis=binary_kernel(A)
    if basis:
        V=S.Matrix(basis);raw=V*M*V.T
        check(all(x%2==0 for x in raw),'nonintegral half form')
        beta=raw.applyfunc(lambda x:(int(x)//2)%2)
        brank,_=binary_kernel(beta.tolist());rad=len(basis)-brank
    else:rad=0
    D=smith_normal_form(M,domain=S.ZZ);diag=[abs(int(D[i,i])) for i in range(M.rows)]
    check(len(basis)==sum(x%2==0 for x in diag),'kernel vs Smith mismatch')
    check(rad==sum(x%4==0 for x in diag),'half-form vs Smith mismatch')
    if M.rows<=3:
        images={tuple(sum(A[i][j]*z[j] for j in range(M.rows))%2 for i in range(M.rows)) for z in itertools.product(range(2),repeat=M.rows)}
        liftable=0
        for x in itertools.product(range(2),repeat=M.rows):
            Ax=[sum(A[i][j]*x[j] for j in range(M.rows)) for i in range(M.rows)]
            if any(y%2 for y in Ax):continue
            if tuple((y//2)%2 for y in Ax) in images:liftable+=1
        check(liftable==2**rad,'direct mod4 lifts mismatch')
    return (M.rows,len(basis),rad)

def run():
    results=[]
    for a,b,c in itertools.product(range(-4,5),repeat=3):
        z=inspect([[a,b],[b,c]])
        if z:results.append(z)
    random.seed(360071)
    for n in range(1,8):
        for _ in range(40):
            A=[[0]*n for i in range(n)]
            for i in range(n):
                for j in range(i,n):A[i][j]=A[j][i]=random.randrange(-16,17)
            z=inspect(A)
            if z:results.append(z)
    out=dict(status='PASS',nonsingular_matrices_checked=len(results),sizes=sorted({x[0] for x in results}),nonzero_half_form_radicals=sum(x[2]>0 for x in results),direct_mod4_lift_tests=sum(x[0]<=3 for x in results))
    (HERE/'algebra_controls_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
