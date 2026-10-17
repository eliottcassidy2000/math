#!/usr/bin/env python3
"""Exact cyclotomic coefficient check of the closed Keller Fourier formula.

No FFT or floating-point test: coefficient vectors agree modulo Phi_p.
"""
import json
from itertools import product
import numpy as np
import sympy as s


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def symbolic():
    x,y,w=s.symbols('x y w');u=1+x*y;z=(2-3*x*y)/x**2
    f1=u**3*z+y*y*u*(4+3*x*y)
    f2=y+3*x*u*u*z+3*x*y*y*(4+3*x*y)
    need(s.cancel(f1-(y*y+3*y/x+2/x**2))==0,'Gamma=0 first coordinate')
    need(s.cancel(f2-(4*y+6/x))==0,'Gamma=0 second coordinate')
    a,b,c,t=s.symbols('a b c t')
    phase=a*(-w**3+w*w)/c**2+b*(-3*w*w+4*w)/c
    need(s.expand(phase.subs(w,c*t)-(a*t*t+4*b*t-c*(a*t**3+3*b*t*t)))==0,'E root-counter substitution')
    z=s.symbols('z');u=1+x*y
    F=s.Matrix([u**3*z+y*y*u*(4+3*x*y),y+3*x*u*u*z+3*x*y*y*(4+3*x*y),2*x-3*x*x*y-x**3*z])
    need(s.factor(F.jacobian([x,y,z]).det())==-2,'literal constant Jacobian')
    phase=(a*F[0]+b*F[1]+c*F[2]).subs(y,t-1/x)
    expanded=s.expand(phase)
    need(s.factor(expanded.coeff(z)-x**3*(a*t**3+3*b*t*t-c))==0,'linear-z cubic selector')
    reduced=s.cancel(expanded.subs(c,a*t**3+3*b*t*t))
    need(s.cancel(reduced-(a*t*t+4*b*t+(a*t+2*b)/x))==0,'retained reciprocal-source phase')
    return True


def run(p):
    pts=np.array(list(product(range(p),repeat=3)),dtype=np.int64)
    x,y,z=pts.T;u=(1+x*y)%p
    F=np.column_stack(((u**3*z+y*y*u*(4+3*x*y))%p,
                      (y+3*x*u*u*z+3*x*y*y*(4+3*x*y))%p,
                      (2*x-3*x*x*y-x**3*z)%p))
    fibre=np.bincount((F[:,0]*p+F[:,1])*p+F[:,2],minlength=p**3)
    total=p**3*sum(int(v)**2 for v in fibre)-p**6
    checked=0;cone_points=0
    # A character coefficient vector has zero value iff all p entries agree.
    for a,b,c in product(range(p),repeat=3):
        if a==b==c==0:continue
        actual=np.bincount((F[:,0]*a+F[:,1]*b+F[:,2]*c)%p,minlength=p)
        expected=np.zeros(p,dtype=np.int64)
        if a:
            inv=pow(a,-1,p);phase=-4*b*b*inv%p
            if (a*a*c-4*b**3)%p==0:expected[phase]+=p*p
            for t in range(p):
                if (a*t**3+3*b*t*t-c)%p==0:expected[(a*t*t+4*b*t)%p]-=p
        elif b:
            for t in range(p):
                if (3*b*t*t-c)%p==0:expected[4*b*t%p]-=p
        else:expected[0]=p*p
        delta=actual-expected
        need(np.all(delta==delta[0]),'exact cyclotomic Fourier identity')
        cone_points+=(a*a*c-4*b**3)%p==0
        checked+=1
    need(cone_points==p*p-1,'complete dual-cone count')
    cone_energy=None
    if p>3:
        cone_energy=p*p*(p-1)*(p**3-p*p+4*p-3)
        need(0<cone_energy<=total,'analytic cone energy and direct Parseval total')
        need(total==p**3*(p-1)*(p*p+2),'exact all-prime collision energy')
    return dict(p=p,frequencies=checked,cone_points=cone_points,
                total_fourier_energy=total,cone_energy=cone_energy,
                fibre_histogram={str(k):int(np.count_nonzero(fibre==k)) for k in np.unique(fibre)})


def prime_power_25():
    p=5;q=25
    add=np.empty((q,q),dtype=np.int64);mul=add.copy()
    for a,b in product(range(q),repeat=2):
        x,y=a%p,a//p;z,w=b%p,b//p
        add[a,b]=(x+z)%p+p*((y+w)%p)
        mul[a,b]=(x*z+2*y*w)%p+p*((x*w+y*z)%p)
    neg=np.array([next(b for b in range(q) if add[a,b]==0) for a in range(q)])
    inv=np.array([0]+[next(b for b in range(1,q) if mul[a,b]==1) for a in range(1,q)])
    x,y,z=np.array(list(product(range(q),repeat=3)),dtype=np.int64).T
    xy=mul[x,y];u=add[1,xy];u2=mul[u,u];u3=mul[u,u2]
    y2=mul[y,y];x2=mul[x,x];x3=mul[x,x2];v=add[4,mul[3,xy]]
    F=np.column_stack((add[mul[u3,z],mul[mul[y2,u],v]],
        add[y,add[mul[3,mul[mul[x,u2],z]],mul[3,mul[mul[x,y2],v]]]],
        add[mul[2,x],neg[add[mul[3,mul[x2,y]],mul[x3,z]]]]))
    tt=np.arange(q);tt2=mul[tt,tt];tt3=mul[tt,tt2]
    trace=np.array([2*(a%p)%p for a in range(q)])
    checked=0;cone=0
    for a,b,c in product(range(q),repeat=3):
        if a==b==c==0:continue
        phase=add[add[mul[a,F[:,0]],mul[b,F[:,1]]],mul[c,F[:,2]]]
        actual=np.bincount(trace[phase],minlength=p);expected=np.zeros(p,dtype=np.int64)
        is_cone=mul[mul[a,a],c]==mul[4,mul[mul[b,b],b]]
        if a:
            marked=neg[mul[4,mul[mul[b,b],inv[a]]]]
            if is_cone:expected[trace[marked]]+=q*q
            roots=tt[add[mul[a,tt3],mul[3,mul[b,tt2]]]==c]
            for t in roots:expected[trace[add[mul[a,mul[t,t]],mul[4,mul[b,t]]]]]-=q
        elif b:
            roots=tt[mul[3,mul[b,tt2]]==c]
            for t in roots:expected[trace[mul[4,mul[b,t]]]]-=q
        else:expected[0]=q*q
        delta=actual-expected
        need(np.all(delta==delta[0]),'exact prime-power trace-character formula')
        checked+=1;cone+=is_cone
    fibres=np.bincount((F[:,0]*q+F[:,1])*q+F[:,2],minlength=q**3)
    total=q**3*sum(int(v)**2 for v in fibres)-q**6
    need(total==q**3*(q-1)*(q*q+2) and cone==q*q-1,'prime-power collision and cone controls')
    return dict(q=q,field='F5[t]/(t^2-2)',character='exp(2pi i Tr(z)/5)',frequencies=checked,
                cone_points=int(cone),total_fourier_energy=total,cone_energy=q*q*(q-1)*(q**3-q*q+4*q-3))


if __name__=='__main__':
    symbolic()
    print(json.dumps(dict(status='FINITE-EXACT all frequencies; universal formula proved algebraically in companion note',
                         runs=[run(p) for p in (3,5,7,11,13,17,19)],prime_power=prime_power_25()),sort_keys=True,separators=(',',':')))
