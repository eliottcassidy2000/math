#!/usr/bin/env python3
"""Parity-free physical 6/55 ceiling: complete coefficient box and finite head.

Standalone exact verifier, no repository mathematics imports. The continuous
knapsack width implementation is inherited from our independent THM4434 audit.
All gates survive Python optimization.
"""
import argparse
import csv
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path
import sys

R = Q(3,14)
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def cross(u,v):
    return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])


def speed_vertices(v):
    vertices = set()
    for free in range(3):
        if v[free] == 0:
            continue
        other = [j for j in range(3) if j != free]
        for values in product((Q(0),Q(1)),repeat=2):
            w = [Q(0)]*3
            for j,x in zip(other,values):
                w[j] = x
            w[free] = -sum(v[j]*w[j] for j in other)/v[free]
            if 0 <= w[free] <= 1:
                vertices.add(tuple(w))
    need(bool(vertices), ("nonempty normalized speed polytope",v))
    return sorted(vertices)


def fill_value(items, amount, reverse):
    total = Q(0)
    for ratio, capacity in sorted(items, reverse=reverse):
        take = min(capacity,amount)
        total += take*ratio
        amount -= take
    need(amount == 0, "continuous resource is filled exactly")
    return total


def slice_width(v,w,delta):
    """LP width via fractional knapsack; no error polygon vertices or clipping."""
    pivot = next(j for j in range(3) if v[j])
    j,k = (pivot+1)%3,(pivot+2)%3
    coefficients = [Q(0)]*3
    coefficients[j],coefficients[k] = -Q(w[k],v[pivot]),Q(w[j],v[pivot])
    items,free_width = [],Q(0)
    for p,ell in zip(v,coefficients):
        if p:
            # Change e to sign(p)e, then shift [-R,R] to [0,2R].
            # The constant objective shift cancels between max and min.
            slope = ell*(1 if p>0 else -1)/abs(p)
            items.append((slope,2*R*abs(p)))
        else:
            free_width += 2*R*abs(ell)
    amount = Q(delta)+R*sum(map(abs,v))
    need(0 <= amount <= 2*R*sum(map(abs,v)), "error slice resource range")
    return free_width+fill_value(items,amount,True)-fill_value(items,amount,False)


def pattern_slope(p):
    S = sum(p)
    defects = tuple(d for d in range(-3*S,3*S+1)
                    if 14*abs(d)<3*S and ((d%3 == 0) == all(x%3 for x in p)))
    unit = all(x%3 for x in p)
    rho = Q(2 if unit else 1,3)
    intercept = Q(4,3)*len(defects) if unit else Q(len(defects))
    maximum = Q(0)
    for negative in range(3):
        if p[negative] == 0:
            continue
        v = tuple(-x if j == negative else x for j,x in enumerate(p))
        for w in speed_vertices(v):
            value = rho*sum((slice_width(v,w,d) for d in defects),Q(0))
            maximum = max(maximum,value)
    return defects,maximum,intercept


def residue_fibers():
    cases=0
    for w in product((1,2),repeat=3):
        kernel=[v for v in product(range(3),repeat=3) if sum(a*b for a,b in zip(v,w))%3==0]
        for v in kernel:
            if v==(0,0,0):
                continue
            unit=all(v)
            need(sum(x==0 for x in v)==(0 if unit else 1), "primitive residue type dichotomy")
            for delta in range(3):
                fiber=[C for C in kernel if all((x-delta*y)%3==0 for x,y in zip(cross(v,C),w))]
                need(len(fiber)==3, "complete affine residue line")
                live=sum(all(C) for C in fiber)
                expected=(2 if delta==0 else 0) if unit else (0 if delta==0 else 1)
                need(live==expected, "exact owner word density and allowed defect")
                cases+=1
    return cases


def parity_free_box():
    exceptions={(0,1,1),(1,1,1),(1,1,2),(1,2,2)}
    patterns=sorted({tuple(sorted(p)) for p in product(range(19),repeat=3)
                     if sum(x!=0 for x in p)>=2
                     and gcd(gcd(p[0],p[1]),p[2])==1
                     and sum(x%3==0 for x in p)<=1
                     and tuple(sorted(p)) not in exceptions})
    need(len(patterns)==747, 'complete parity-free M<=18 coefficient universe')
    digest=sha256()
    equalities=[]
    for p in patterns:
        ds,slope,B=pattern_slope(p)
        if not ds:
            need((p,slope,B)==((0,1,2),Q(0),Q(0)),
                 'empty defect list is discharged as N=0, never 0<0')
        need(slope<=Q(15,98), ('complete coefficient slope',p,slope))
        need(B<=Q(2*sum(p),7)+Q(4,3), ('all exact defect intercepts',p,B))
        if slope==Q(15,98):
            equalities.append(p)
        digest.update(repr((p,ds,slope,B)).encode())
    need(equalities==[(1,7,8)], 'unique normalized slope equality')
    need(pattern_slope((1,2,2))==((0,),Q(3,14),Q(4,3)),
         'the only additional count exception has fixed short intercept')
    need(pattern_slope((1,1,1))[1]==Q(2,7), 'norm-three count hostile')
    need(pattern_slope((1,1,2))[1]==Q(2,7), 'norm-four count hostile')
    return len(patterns),Counter(sum(x!=0 for x in p) for p in patterns),equalities,digest.hexdigest()


def cutoff_controls():
    gap=Q(14,55)-Q(15,98)
    need(gap==Q(547,5390), 'physical count threshold gap')
    A,B=Q(2,7)/gap,Q(4,3)/gap
    need((A,B)==(Q(1540,547),Q(21560,1641)), 'count cutoff constants')
    need(A*18+B==Q(104720,1641)<64, 'short relation S<=18 tail')
    g19=Q(3*19*19,16)-A*19-B
    need(g19==Q(27763,26256)>0 and Q(3*19,8)>A, 'quadratic S>=19 tail')
    need(Q(6,49)+Q(4,7*19)==Q(142,931)<Q(15,98), 'M>=19 slope')
    need(Q(3,14)+Q(4,3*34)<Q(14,55), 'norm-five count strict tail c>=34')
    need(Q(9,98)+Q(6,7*50)<Q(6,55), 'additive physical strict tail c>=50')
    need(Q(3,49)+Q(6,7*18)<Q(6,55), 'norm-four physical strict tail c>=18')
    # Rational identity controls for the projected l1-ball area proof.
    for a,b in product((Q(1,17),Q(2,9),Q(1,2),Q(4,5),Q(1)),repeat=2):
        lhs=8*(a*b+a+b)-3*(a+b)*(a+1)*(b+1)
        rhs=3*a*(1-a)*(b+1)+3*b*(1-b)*(a+1)+2*a*(1-b)+2*b*(1-a)
        need(lhs==rhs and lhs>=0, 'positive l1-area coefficient identity')
        if a<1 or b<1:
            need(lhs>0, 'strict area for distinct positive speeds')
    # The cap/one-double-roof half-integral identity is uniform in p,q<=c.
    for p,q,c in ((1,5,11),(5,11,11),(2,20,20),(7,19,23),(1,2,3)):
        H=Q(3,7*c); intercept=Q(3*(p+q),14*p*q); slope=Q(2,p*q)
        need(intercept>=H, 'cap lies below double-roof intercept')
        half=H*(intercept-H/2)/slope
        need(half==Q(9*(c*(p+q)-p*q),196*c*c)<=Q(9,196),
             'norm-four half-integral bound')
    return A,B,g19


def raw_mass(w):
    """Enumerate first and third carrier coordinates; solve the middle one."""
    a,b,c=w
    P=[3*(sum(w)-x) for x in w]
    limits=[(x-1)//14 for x in P]
    mass=N=0
    E=[0,0,0]
    D=42*a*b*c
    for x in range(-limits[0],limits[0]+1):
        if x%3==0:
            continue
        for z in range(-limits[2],limits[2]+1):
            if z%3==0 or (-a*x-c*z)%b:
                continue
            y=(-a*x-c*z)//b
            if y%3==0 or abs(y)>limits[1]:
                continue
            C=(x,y,z)
            terms=[min(18*a*b,3*w[i]*(P[i]-14*abs(C[i]))) for i in range(3)]
            need(min(terms)>0, 'strict complete raw carrier')
            N+=1
            mass+=min(terms)
            for i in range(3):
                E[i]+=terms[i]
    return D,tuple(E),mass,N


def complete_head(height, tsv):
    universe=[w for w in combinations([n for n in range(1,height+1) if n%3],3)
              if gcd(gcd(w[0],w[1]),w[2])==1]
    maximum=Q(0)
    equalities=[]
    network_maximum=Q(0)
    network_equalities=[]
    rows=[]
    old_bad=[]
    for w in universe:
        D,E,mass_numerator,N=raw_mass(w)
        mass=Q(mass_numerator,D)
        network=Q(min(E),D)
        need(mass<=Q(6,55), ('complete physical head ceiling',w,mass))
        need(network<=Q(6,55), ('complete selected-network head ceiling',w,network))
        need(mass<=network, ('physical is bounded by every projection',w))
        if mass==Q(6,55):
            equalities.append(w)
        if network==Q(6,55):
            network_equalities.append(w)
        if mass>Q(6,77):
            old_bad.append((w,mass))
        maximum=max(maximum,mass)
        network_maximum=max(network_maximum,network)
        rows.append((w,mass,N,D,E,mass_numerator))
    need(equalities==[(1,10,11)], 'complete head unique physical equality')
    need(network_equalities==[(1,10,11)], 'complete head unique selected-network equality')
    controls={(2,5,7):Q(22,245),(1,7,8):Q(31,392),(1,5,11):Q(6,77),
              (1,10,11):Q(6,55),(2,11,20):Q(11,140),(1,2,4):Q(0)}
    observed={row[0]:row[1] for row in rows}
    for w,expected in controls.items():
        need(observed[w]==expected, ('positive and hostile physical control',w))
    if tsv:
        with Path(tsv).open('w',newline='',encoding='utf-8') as f:
            writer=csv.writer(f,delimiter='\t',lineterminator='\n')
            writer.writerow(('a','b','c','denominator','E0_numerator','E1_numerator','E2_numerator','mass_numerator','raw_carriers'))
            for w,mass,N,D,E,mass_numerator in rows:
                writer.writerow((*w,D,*E,mass_numerator,N))
    top=[(row[0],row[1],row[2]) for row in sorted(rows,key=lambda t:t[1],reverse=True)[:10]]
    return len(rows),maximum,equalities,network_maximum,network_equalities,len(old_bad),top


def main():
    sys.stdout.reconfigure(newline='\n')
    parser=argparse.ArgumentParser()
    parser.add_argument('--head-tsv')
    args=parser.parse_args()
    box=parity_free_box()
    fibers=residue_fibers()
    cutoff=cutoff_controls()
    head=complete_head(63,args.head_tsv)
    print('universe=primitive distinct positive sorted ternary-unit triples; no parity filter')
    print('coefficient_box='+str(box))
    print('exact_owner_fibers='+str(fibers))
    print('empty_defect_rule=N=0 directly; strict N<F+B only for nonempty defect lists')
    print('cutoff_constants='+str(cutoff))
    print('general_tail=c>=64; norm5_tail=c>=34; norm3_tail=c>=50; norm4_tail=c>=18')
    print('complete_head63='+str(head))
    print('all_height_candidate=physical mass<=6/55, equality iff(1,10,11); proof audit required')
    print('selected_network_extension=complete head PASS; additive fixed-projection and norm4 fixed-roof inputs under audit')
    print('optimization_live_gates='+str(CHECKS))


if __name__=='__main__':
    main()
