#!/usr/bin/env python3
"""Exact root-box and resultant certificates for the four-anchor domain."""
from fractions import Fraction as Q
import hashlib
import json
import sympy as S

v,f,w=S.symbols('v f w')
F=v**5-13*v**4+55*v**3-84*v**2+35*v
B=F-f
C=v**4-12*v**3+45*v**2-56*v+15
D=v**4-11*v**3+36*v**2-35*v+5
GATES=0
TRACE=hashlib.sha256()


def gate(ok,label):
    global GATES
    if ok != True:
        raise RuntimeError(label+': '+str(ok))
    GATES+=1
    TRACE.update((label+'\n').encode())


def zero(expr,label):
    gate(S.simplify(expr)==0,label)


def enclosure(expr,lo,hi):
    a=b=Q(0)
    for c in S.Poly(expr,v).all_coeffs():
        values=[a*lo,a*hi,b*lo,b*hi]
        a,b=min(values)+Q(c),max(values)+Q(c)
    return a,b


def boxes(poly,intervals,targets,name):
    rows=[]
    previous=Q(0)
    for j,((lo,hi),(lower,upper)) in enumerate(zip(intervals,targets)):
        gate(previous<lo<hi,name+' disjoint positive box '+str(j))
        gate(poly.subs(v,S.Rational(lo))*poly.subs(v,S.Rational(hi))<0,
             name+' exact root sign change '+str(j))
        bounds=enclosure(F,lo,hi)
        gate(lower<bounds[0]<=bounds[1]<upper,name+' full F image enclosure '+str(j))
        rows.append(dict(box=[str(lo),str(hi)],F_bounds=list(map(str,bounds))))
        previous=hi
    return rows


def main():
    cb=[(Q(136,373),Q(101,277)),(Q(461,262),Q(783,445)),
        (Q(2114,537),Q(2425,616)),(Q(2245,378),Q(3807,641))]
    ctargets=[(Q(4),Q(41,10)),(Q(-7),Q(-6)),(Q(14),Q(15)),(Q(-19),Q(-18))]
    crows=boxes(C,cb,ctargets,'C')
    gate(enclosure(S.diff(F,v),*cb[0])[1]<-6,'simple C collision derivative')
    rc=f**4+6*f**3-285*f**2-778*f+7125
    zero(S.resultant(B,C,v)-rc,'C full resultant identity')
    zero(D-(v*v-6*v+1)*(v*v-5*v+5),'D two-quadratic factorization')
    dd=[3-2*S.sqrt(2),(5-S.sqrt(5))/2,(5+S.sqrt(5))/2,3+2*S.sqrt(2)]
    values=[-16+14*S.sqrt(2),(15-15*S.sqrt(5))/2,
            (15+15*S.sqrt(5))/2,-16-14*S.sqrt(2)]
    for j,(root,value) in enumerate(zip(dd,values)):
        zero(D.subs(v,root),'D exact root '+str(j))
        zero(F.subs(v,root)-value,'D exact F value '+str(j))
    gate(Q(7,5)**2<2<Q(3,2)**2 and Q(11,5)**2<5<Q(23,10)**2,
         'D roots positive and ordered by disjoint radical intervals')
    gate(Q(8,7)**2<2 and 9800<9801,'D upper boundary 0<L<19/5')
    rd=(f*f-15*f-225)*(f*f+32*f-136)
    zero(S.resultant(B,D,v)-rd,'D full resultant identity')
    den=S.Poly(S.expand(w**5*B.subs(v,1/w)),w)
    num=S.Poly(S.expand(w**4*D.subs(v,1/w)),w)
    nu=[]
    for j in range(9):
        nu.append(S.expand(num.nth(j)-sum(den.nth(k)*nu[j-k] for k in range(1,min(j,5)+1))))
    zero(S.Matrix(5,5,lambda i,j:nu[i+j]).det(method='domain-ge')-rd,
         'fifth moment determinant equals resultant')
    eb=[(Q(131,472),Q(68,245)),(Q(564,401),Q(737,524)),
        (Q(1521,457),Q(1957,588)),(Q(4391,815),Q(4655,864))]
    etargets=[(Q(43,10),Q(22,5)),(Q(-10),Q(-9)),(Q(26),Q(28)),(Q(-63),Q(-62))]
    erows=boxes(S.diff(F,v),eb,etargets,'critical')
    rb=3125*f**4+125758*f**3-4825473*f**2-30354226*f+211485225
    zero(S.resultant(B,S.diff(F,v),v)-rb,'B full discriminant identity')
    gate(Q(19,5)<4<Q(41,10)<Q(43,10)<Q(22,5)<5,'strict interval hierarchy')
    zero(S.diff(F,v).subs(v,0)-35,'simple zero-root boundary')
    manifest=dict(C_boxes=crows,critical_boxes=erows,
                  B_domain='[0,kappa], unique R_B root in (43/10,22/5)',
                  C_domain='[0,rho], unique R_C root in (4,41/10)',
                  D_domain='[0,14sqrt(2)-16]',
                  simultaneous_domain='[0,14sqrt(2)-16]')
    print('UNIVERSE exact eight root boxes, formal resultant/moment identities, radical boundary values; no shape census')
    print(json.dumps(manifest,sort_keys=True,separators=(',',':')))
    print('CLAIM exact B/C/D admissible hierarchy; B-only domain contained in [0,5]')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__=='__main__':
    main()
