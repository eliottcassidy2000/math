#!/usr/bin/env python3
"""Finite proof head and exact determinant controls for two-ray closure.

Primary lattice enumeration solves one congruence per row. The separately
written one-ray verifier supplies only literal shifted-sheet controls.
"""

import importlib.util
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def primitive(C):
    g = gcd(gcd(abs(C[0]), abs(C[1])), abs(C[2]))
    v = tuple(x//g for x in C)
    return v if next(x for x in v if x) > 0 else tuple(-x for x in v)


def cross(u,v):
    return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])


def row_carriers(w):
    a,b,c = w
    bx,by,bz = [(3*(sum(w)-v)-1)//14 for v in w]
    g = gcd(b,c)
    modulus = c//g
    inverse = pow(b//g,-1,modulus) if modulus != 1 else 0
    live = set()
    for x in range(-bx,bx+1):
        if x % 3 == 0 or (a*x) % g:
            continue
        residue = (-(a*x//g)*inverse) % modulus
        first = -by+(residue+by) % modulus
        for y in range(first,by+1,modulus):
            z = -(a*x+b*y)//c
            if y % 3 and z % 3 and abs(z)<=bz:
                live.add((x,y,z))
    return live


def projections(w,live):
    out = [Q(0)]*3
    for C in live:
        for i in range(3):
            j,k = [v for v in range(3) if v != i]
            value = min(Q(3,7*w[2]),Q(3*(w[j]+w[k])-14*abs(C[i]),14*w[j]*w[k]))
            need(value>0,('live projection',w,C))
            out[i] += value
    return tuple(out)


def audit(w,live):
    vv = sorted({primitive(C) for C in live})
    need(len(vv)==2,('exactly two unoriented rays',w,vv))
    u,v = vv
    ms = [max(map(abs,t)) for t in vv]
    for direction,M in zip(vv,ms):
        need(sum(map(abs,direction))>=16 and M>=7,('short-relation obstruction',w,direction))
        B = min(Q(3*(sum(w)-w[i]),14*abs(direction[i])) for i in range(3))
        K = (B.numerator-1)//B.denominator
        actual = {C for C in live if primitive(C)==direction}
        predicted = {tuple(sign*k*x for x in direction) for k in range(1,K+1)
                     if k%3 for sign in (-1,1)}
        need(actual==predicted,('complete per-ray address',w,direction))
        need(B<Q(3*w[2],7*M),('narrow coordinate',w,direction))
        need(len(actual)<Q(4,3)*(B+1),('residue deletion',w,direction))
    product = cross(u,v)
    k = Q(product[2],w[2])
    need(k.denominator==1 and k and k.numerator%3==0,('integral determinant multiple3',w,vv,product))
    need(product==tuple(k*x for x in w),('exact oriented determinant',w,vv,product))
    need(2*ms[0]*ms[1]>=3*w[2],('product lower bound',w,ms))
    need(Q(1,ms[0])+Q(1,ms[1])<=Q(1,7)+Q(14,3*w[2]),('reciprocal envelope',w,ms))
    values = projections(w,live)
    need(max(values)<Q(12,343)+Q(16,7*w[2]),('every-projection analytic envelope',w,values))
    if w[2]>=55:
        need(max(values)<Q(6,77),('analytic tail',w,values))
    return values,vv,k


def main():
    path = Path(__file__).with_name('lrc14_one_ray_overnight_hexagon_sep05.py')
    spec = importlib.util.spec_from_file_location('literal_one_ray',path)
    literal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(literal)
    need(Q(12,343)+Q(16,7*55)==Q(1444,18865)<Q(6,77),'tail threshold')
    digest = sha256()
    rows = selected = 0
    leader = (Q(0),None,None)
    for w in combinations([v for v in range(1,55,2) if v%3],3):
        if gcd(gcd(w[0],w[1]),w[2])!=1:
            continue
        rows += 1
        live = row_carriers(w)
        # Independent complete boxes verify the enumeration before filtering.
        need(live==literal.carriers(w),('row/box equality',w))
        if len({primitive(C) for C in live})!=2:
            continue
        selected += 1
        values,vv,k = audit(w,live)
        need(values==literal.literal_projection_data(w)[0],('raw/literal sheets',w))
        need(max(values)<Q(6,77),('complete two-ray finite head',w,values))
        if max(values)>leader[0]:
            leader=max(values),w,values
        digest.update(repr((w,tuple(sorted(live)),values,tuple(vv),k)).encode())
    need(selected==192,('finite head count',selected))
    need(leader[0]==Q(114,2233) and leader[1]==(11,23,29),('finite head leader',leader))
    print('PROVED: complete TWO-RAY support implies EVERY network projection<6/77')
    print('ALL-HEIGHT ENVELOPE E_i<12/343+16/(7c); c>=55 closes analytically')
    print('PROOF HEAD: primitive sorted odd ternary-unit triples c<55;',rows,'total;',selected,'two-ray')
    print('FINITE HEAD LEADER',leader)
    for w in ((1,11,55),(1,5,101),(5,49,251)):
        live=row_carriers(w)
        values,vv,k=audit(w,live)
        need(live==literal.carriers(w),('wide row/box equality',w))
        need(values==literal.literal_projection_data(w)[0],('wide literal projection',w))
        print('TAIL CONTROL',w,'N',len(live),'directions',vv,'determinant multiplier',k,
              'E',tuple(map(str,values)))
    print('SEMANTIC SHA256',digest.hexdigest())
    print('CHECKS',CHECKS,'LITERAL CHECKS',literal.CHECKS)
    print('OPEN: three or more primitive live directions; entry; synchronization; LRC14')


if __name__=='__main__':
    main()
