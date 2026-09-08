#!/usr/bin/env python3
"""Exact controls for THM-4457; standard library only, no numerical PDE claim.

Universe: 3^8 integer matrices with the first eight entries in {-1,0,1}
and the ninth chosen for trace zero, against all six coordinate shears with
amplitudes -3,...,3. Additional rational rotation, sharpness, pressure-profile,
and orbit-fixing controls are independent of that finite universe.
Run: python3 04-computation/euler_shear_distance_budget_20260908.py
All checks remain active under python -O. Stdout is deterministic.
"""
from fractions import Fraction as Q
from itertools import product
import hashlib
import json


def require(value, label):
    if not value:
        raise RuntimeError(label)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mv(m, v):
    return tuple(dot(m[3*i:3*i+3], v) for i in range(3))


def transpose(m):
    return tuple(m[3*j+i] for i in range(3) for j in range(3))


def mm(a, b):
    bt = transpose(b)
    return tuple(dot(a[3*i:3*i+3], bt[3*j:3*j+3])
                 for i in range(3) for j in range(3))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def omega(m):
    return (m[7]-m[5], m[2]-m[6], m[3]-m[1])


def check_matrix(m, s):
    w = omega(m)
    e = sub(m, s)
    a = dot(w, mv(m, w))
    w2 = dot(w, w)
    # Squared universal sharp inequality, including the zero-vorticity boundary.
    require(6*a*a <= 7*dot(e, e)*w2*w2, 'sharp Frobenius inequality')
    ws, we = omega(s), omega(e)
    require(mv(s, ws) == (0, 0, 0), 'right self-vorticity annihilation')
    require(mv(transpose(s), ws) == (0, 0, 0), 'left annihilation')
    require(dot(w, mv(s, w)) == dot(we, mv(s, we)), 'scalar cancellation')


def run():
    identity = (1,0,0,0,1,0,0,0,1)
    count = 0
    coordinate_shears = []
    for i in range(3):
        for j in range(3):
            if i != j:
                for a in range(-3, 4):
                    s = [0]*9
                    s[3*i+j] = a
                    coordinate_shears.append(tuple(s))
    for x in product((-1,0,1), repeat=8):
        m = x + (-x[0]-x[4],)
        for s in coordinate_shears:
            check_matrix(m, s)
            count += 1
    # Dense rational rotation: conjugation preserves all consequence objects.
    rotation = (Q(-2,3),Q(2,15),Q(11,15), Q(2,3),Q(-1,3),Q(2,3),
                Q(1,3),Q(14,15),Q(2,15))
    require(mm(rotation, transpose(rotation)) == identity, 'orthogonal control')
    rotated = 0
    for a in (-7,-1,0,2,9):
        s = (0,a,0,0,0,0,0,0,0)
        for z in (-3,-1,0,2):
            m = (z,2,-1,3,-2*z,1,-2,4,z)
            mr = mm(mm(rotation,m),transpose(rotation))
            sr = mm(mm(rotation,s),transpose(rotation))
            check_matrix(mr,sr)
            require(dot(omega(mr),mv(mr,omega(mr))) == dot(omega(m),mv(m,omega(m))),
                    'orthogonal stretching covariance')
            require(dot(sub(mr,sr),sub(mr,sr)) == dot(sub(m,s),sub(m,s)),
                    'orthogonal residual norm covariance')
            rotated += 1
    # Rational sharpness family; alpha=7/6 while e^2=7/6+epsilon^2/2.
    s = (Q(1,2),Q(1,2),0,Q(-1,2),Q(-1,2),0,0,0,0)
    sharp = []
    for n in (1,2,10,100,1000):
        eps = Q(1,n)
        m = (Q(7,6),0,0,0,Q(-5,6),-eps/2,0,eps/2,Q(-1,3))
        check_matrix(m,s)
        w = omega(m)
        require(w == (eps,0,0),'sharp vorticity')
        alpha = dot(w,mv(m,w))/dot(w,w)
        e2 = dot(sub(m,s),sub(m,s))
        require(alpha == Q(7,6) and e2 == Q(7,6)+eps*eps/2,'sharp family')
        ratio2 = alpha*alpha/e2
        require(ratio2 < Q(7,6),'limit approached from below')
        sharp.append([n,str(ratio2)])
    # A pure shear has arbitrarily large norm and zero self-stretching.
    for a in (0,1,10,10**12):
        s = (0,a,0,0,0,0,0,0,0)
        require(mv(s,omega(s)) == (0,0,0),'pure shear hostile')
    # Poisson identity for the explicit arctangent profile, rational cos(theta).
    pressure_cases = 0
    for delta in (Q(1,100),Q(1,7),Q(1),Q(3)):
        r = 1/(1+delta)
        for c in (Q(-1),Q(-3,5),Q(0),Q(3,5),Q(1)):
            derivative = ((1+delta)*c-1)/((1+delta)**2-2*(1+delta)*c+1)
            poisson = (1-r*r)/(1-2*r*c+r*r)
            require(derivative == (poisson-1)/2,'Poisson identity')
            require(-1/(2+delta) <= derivative <= 1/delta,'sharp slope range')
            # a=coupling=1, m=e1: adverse positive eigenvalue <= 2/(2+delta).
            require(-2*derivative <= 2/(2+delta),'one-sided pressure bound')
            pressure_cases += 1
        require(2/delta > 2/(2+delta),'negative coupling hostile at theta zero')
    # An integer transvection can have huge norm yet fix the entire speed line.
    # v=(1,...,13), m=(2,-1,0,...), a=e3 gives m.v=m.a=0.
    v = tuple(range(1,14))
    covector = (2,-1)+(0,)*11
    direction = (0,0,1)+(0,)*10
    require(dot(covector,v)==dot(covector,direction)==0,'orbit-fixing hypotheses')
    for k in (1,13,10**6):
        image = tuple(v[i]+k*direction[i]*dot(covector,v) for i in range(13))
        require(image==v,'pointwise speed-line fixing')
    report = {
        'status':'FINITE-EXACT controls; universal proof is in THM-4457',
        'integer_matrices':3**8,
        'matrix_shear_pairs':count,
        'rational_rotation_cases':rotated,
        'sharp_constant_squared':'7/6',
        'sharpness_squared_ratios':sharp,
        'profile_identity_cases':pressure_cases,
        'hostiles':['pure shear has zero self-stretching',
                    'negative pressure coupling destroys one-sided bound',
                    'unimodular norm growth can fix every LRC phase'],
        'scope':'No PDE solver, no exact optimization of the shear distance, no LRC row exclusion.'
    }
    serialized = json.dumps(report,sort_keys=True,separators=(',',':'))
    print(json.dumps(report,sort_keys=True,indent=2))
    print('semantic_sha256='+hashlib.sha256(serialized.encode()).hexdigest())
    print('PASS')


if __name__ == '__main__':
    run()
