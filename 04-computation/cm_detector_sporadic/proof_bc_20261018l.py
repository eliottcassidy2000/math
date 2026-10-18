#!/usr/bin/env python3
"""Exact checks for THM-4711; uses existing mechanism() on integer sequences.
Usage: python3 .../proof_bc_20261018l.py --bound 199 --out FILE
Proof is in the theorem. Finite checks include entire polynomials, carry
indices, the characteristic-zero parameter tangent and the Euler prefactor.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sympy as s
from ordinary_detector_and_mechanism import mechanism
from sporadic_cm_detector import primes_between, zagier_rec

z=s.Symbol('z')
checks=0

def require(ok, msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)

def series(a,n):
    h=[Q(1)];j=[Q(0)];t=Q(0)
    for k in range(1,n):
        h.append(h[-1]*(a+k-1)**2/k**2)
        t+=2*(1/(a+k-1)-Q(1,k))
        j.append(h[-1]*t)
    return h,j

def convolve(a,b,n):
    return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(n)]

def modq(q,p):
    q=Q(q)
    require(q.denominator%p!=0,'nonintegral coefficient')
    return q.numerator*pow(q.denominator,-1,p)%p

def poly(a,p):return s.Poly.from_list(a[::-1],z,modulus=p)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=199);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    n=40;H=[Q(0)]
    for k in range(1,2*n):H.append(H[-1]+Q(1,k))
    fc=[Q(comb(2*k,k)**3,64**k) for k in range(n)]
    gc=[fc[k]*6*(H[2*k]-H[k]) for k in range(n)]
    h,j=series(Q(1,4),n);b,c=series(Q(3,4),n)
    hh=convolve(h,h,n);hj=convolve(h,j,n);bb=convolve(b,b,n);bc=convolve(b,c,n)
    require(fc==hh,'Clausen rational coefficients')
    require(gc==hj,'Clausen tangent rational coefficients')
    require(fc==[bb[k]-(bb[k-1] if k else 0) for k in range(n)],'Euler square')
    require(gc==[bc[k]-(bc[k-1] if k else 0) for k in range(n)],'Euler tangent')
    # Derivative slopes are three numerator shifts and denominator slopes 2+1.
    for k in range(n):
        tangent=3*sum((1/(Q(1,2)+j) for j in range(k)),Q(0))-3*H[k]
        require(fc[k]*tangent==gc[k],'deformation tangent mismatch')
    rows=[]
    for p in primes_between(7,a.bound):
        m=(p-1)//2;hs=[0]
        for k in range(1,p):hs.append((hs[-1]+pow(k,-1,p))%p)
        cases={}
        for r in (2,3):
            u=[comb(2*k,k)**r for k in range(2*p)]
            require(all((u[p+k]-u[1]*u[k])%p==0 for k in range(p)),f'Lucas {p,r}')
            fd=[u[k]*pow(4**r,-k,p)%p for k in range(p)]
            gd=[((u[p+k]-u[1]*u[k])//p)*pow(u[1],-1,p)*pow(4**r,-k,p)%p for k in range(p)]
            for k in range(p):
                expected=fd[k]*2*r*(hs[2*k]-hs[k])%p if k<=m else 0
                require(gd[k]==expected,f'shift/carry {p,r,k}')
                if k>m:require(fd[k]==gd[k]==0,'tail')
            f=poly(fd,p);g=poly(gd,p)
            if r==2:
                R=f*f+s.Poly(z,z,modulus=p)*(f*g.diff()-f.diff()*g)
                require(R==s.Poly((1-z)**(p-1),z,modulus=p),'Legendre Wronskian')
                # Explicit hostile to deleting the Frobenius term.
                require(s.Poly(1-z,z,modulus=p)*R==s.Poly(1-z**p,z,modulus=p),'degree-p boundary')
            else:
                aa=Q(1,4) if p%4==1 else Q(3,4);d=(p-int(4*aa))//4
                bh,cj=series(aa,p)
                B=poly([modq(t,p) for t in bh[:d+1]],p);C=poly([modq(t,p) for t in cj[:d+1]],p)
                for k in range(d+1,p):
                    require(modq(bh[k],p)==modq(cj[k],p)==0,'quarter omitted tail')
                S=s.Poly(1 if p%4==1 else 1-z,z,modulus=p)
                require(f==S*B*B and g==S*B*C,'cube product identities')
                R=B*B+s.Poly(z,z,modulus=p)*(B*C.diff()-B.diff()*C)
                require(R==s.Poly((1-z)**(2*d),z,modulus=p),'quarter Wronskian')
                require(s.gcd(f,g).monic()==(S*B).monic(),'exact gcd is radical')
            cases[str(r)]=mechanism(u,p,r==3)
            require(cases[str(r)]['rad|V'] if r==3 else cases[str(r)]['squarefree'] and cases[str(r)]['gcd(F,V)=1'],'existing mechanism')
        # Open-route controls, tested rather than promoted: Franel and Apery.
        if p<=199:
            fr=zagier_rec(7,2,-8,2*p)
            au=[1,5]
            for k in range(1,2*p-1):
                num=(2*k+1)*(17*k*k+17*k+5)*au[-1]-k**3*au[-2]
                require(num%(k+1)**3==0,'Apery recurrence integrality')
                au.append(num//(k+1)**3)
            cases['Franel_OPEN']=mechanism(fr,p,False)
            cases['Apery_OPEN']=mechanism(au,p,True)
            require(cases['Franel_OPEN']['gcd(F,V)=1'] and cases['Franel_OPEN']['squarefree'],'Franel control')
            require(cases['Apery_OPEN']['rad|V'],'Apery control')
        rows.append({'p':p,'mechanism':cases})
    out={'status':'FINITE-EXACT; all-prime proof in THM-4711','prime_bound':a.bound,'prime_count':len(rows),'rational_coefficients_checked':n,'live_gates':checks,'rows':rows}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='rows'})
if __name__=='__main__':main()
