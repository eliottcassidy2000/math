#!/usr/bin/env python3
"""Exact roof integration in two independent coordinate representations.

Named carrier hostiles, not another height atlas. Every gate stays enabled
under python -O. No mathematical library is required.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import importlib.util

CHECKS=0


def need(value, detail):
    global CHECKS
    CHECKS+=1
    if not value:
        raise RuntimeError(detail)


def clip(poly, A, B, D):
    out=[]
    for u,v in zip(poly,poly[1:]+poly[:1]):
        fu=A*u[0]+B*u[1]-D
        fv=A*v[0]+B*v[1]-D
        if fu<=0:
            out.append(u)
        if (fu<0<fv) or (fv<0<fu):
            t=fu/(fu-fv)
            out.append((u[0]+t*(v[0]-u[0]),u[1]+t*(v[1]-u[1])))
    return out


def integrate_affine(poly, f):
    if len(poly)<3:
        return Q(0)
    u=poly[0]
    answer=Q(0)
    for v,z in zip(poly[1:-1],poly[2:]):
        twice=(v[0]-u[0])*(z[1]-u[1])-(v[1]-u[1])*(z[0]-u[0])
        answer+=abs(twice)*(f(u)+f(v)+f(z))/6
    return answer


def polygon_roof(w,i):
    a,b,c=map(Q,w)
    widths=[3*(a+b+c-x)/14 for x in (a,b,c)]
    X,Y=widths[:2]
    poly=[(-X,-Y),(X,-Y),(X,Y),(-X,Y)]
    poly=clip(clip(poly,a,b,c*widths[2]),-a,-b,c*widths[2])
    A,B=((Q(1),Q(0)),(Q(0),Q(1)),(-a/c,-b/c))[i]
    j,k=[d for d in range(3) if d!=i]
    denom=w[j]*w[k]
    q=3/(7*c)
    threshold=widths[i]-q*denom
    need(threshold>0,('positive roof plateau width',w,i))
    plateau=clip(clip(poly,A,B,threshold),-A,-B,threshold)
    plus=clip(poly,-A,-B,-threshold)
    minus=clip(poly,A,B,-threshold)
    raw=integrate_affine(plateau,lambda _:q)
    raw+=integrate_affine(plus,lambda v:(widths[i]-A*v[0]-B*v[1])/denom)
    raw+=integrate_affine(minus,lambda v:(widths[i]+A*v[0]+B*v[1])/denom)
    return 2*raw/(9*c)


def formula(w,i):
    a,b,c=map(Q,w)
    j,k=[d for d in range(3) if d!=i]
    return (6*c*(a+b+c)-6*w[j]*w[k]-2*a*b)/(343*c*c)


def main():
    count=0
    # The continuous formula does not need parity, residues, or primitivity.
    for w in combinations(range(1,13),3):
        js=tuple(formula(w,i) for i in range(3))
        need(js[0]<js[1]<js[2]<Q(12,343),'continuous order and largest bound')
        need(js[0]<Q(10,343),'smallest continuous bound')
        need(js==tuple(polygon_roof(w,i) for i in range(3)),('independent polygon integrals',w))
        count+=1
    print('CONTINUOUS POLYGON CONTROLS',count,'all rational integrals match')
    path=Path(__file__).with_name('lrc14_two_ray_overnight_hexagon_sep05.py')
    spec=importlib.util.spec_from_file_location('rows',path)
    rows=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rows)
    for w in ((1,5,11),(1,19,79),(19,23,29),(5,23,47),(41,47,49),(1,137,277),(5,191,199),(7,611,613)):
        live=rows.row_carriers(w)
        es=rows.projections(w,live)
        js=tuple(formula(w,i) for i in range(3))
        need(js==tuple(polygon_roof(w,i) for i in range(3)),('hostile polygon replay',w))
        ratios=tuple(e/j for e,j in zip(es,js))
        if w==(19,23,29):
            need(all(e>j for e,j in zip(es,js)),'density alone fails all three coordinates')
        if w==(1,19,79):
            need(es[1]<es[0] and js[0]<js[1],'actual selector reversal is arithmetic discrepancy')
        if w==(5,23,47):
            need(ratios[2]==Q(335251,166175)>2,'factor-two discrete/continuous hostile')
            need(len({rows.primitive(x) for x in live})==3,'hostile is truly multi-direction')
        print('NAMED',w,'N',len(live),'E',es,'J',js,'E/J',ratios)
    for T in (3,5,11,101):
        eps=Q(1,4)
        live=[(x,y) for x in (-1,1) for y in range(1-T,T)]
        E=sum(min(Q(1),T-abs(y)) for x,y in live)
        J=Q(4,3)*(1+eps)*(2*T-1)
        need(E==4*T-2 and E/J==Q(6,5),'abstract transverse boundary discrepancy')
        print('ABSTRACT STRIP',T,'E',E,'J',J,'E-J',E-J)
    print('CHECKS',CHECKS,'ROW HELPER CHECKS',rows.CHECKS)
    print('PROVED: continuous roof formula/order/bounds; named discrete hostiles')
    print('OPEN: owner-coset boundary control for arbitrary support; no universal network conclusion')


if __name__=='__main__':
    main()
