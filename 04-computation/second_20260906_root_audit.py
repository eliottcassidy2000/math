#!/usr/bin/env python3
"""Independent SymPy Smith and literal balanced-entry audit, no producer imports."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, prod
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def val(x,p):
    need(x != 0,'finite valuation argument')
    x = abs(int(x))
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def predicted(a,b,p):
    A,B,C = val(a,p),val(b,p),val(b-a,p)
    e,f = min(A,B,C),max(A,B,C)
    if A == e:
        r = min(2*e,e+int(p==2))
        return (0,0,r,4*e-r,2*e+2*f)
    if p != 2:
        r = min(f,2*e)
        return (0,0,r,f+3*e-r,3*f+e)
    d = f-e
    if d >= 2:
        r = min(2*e,e+d+1)
        return (0,0,r,4*e+d+1-r,4*e+3*d-1)
    if e == 0:
        return (0,0,0,2,2)
    if e == 1:
        return (0,0,2,5,5)
    u,v = a//2**(e+1),b//2**e
    epsilon = int((u*pow(v,-1,4)) % 4 == 1)
    return (0,0,e+2,3*e+2-epsilon,4*e+epsilon)


def smith_audit():
    bank = {(a,b) for a in range(-12,13) for b in range(-12,13)
            if a and b and a != b}
    for p in (2,3,5,7):
        for e in range(0,6):
            for d in range(1,5):
                for sign in (-1,1):
                    a,b = p**(e+d),sign*p**e
                    bank.update(((a,b),(b,a),(a,a-b)))
    count = 0
    for a,b in sorted(bank):
        # Direct ordinary monomial values/derivatives; not the residual matrix.
        rows = []
        for x,m in ((0,2),(a,2),(b,1)):
            rows.append([x**k for k in range(5)])
            if m == 2:
                rows.append([0]+[k*x**(k-1) for k in range(1,5)])
        matrix = Matrix(rows)
        diagonal = smith_normal_form(matrix,domain=ZZ)
        full = [abs(int(diagonal[i,i])) for i in range(5)]
        need(abs(matrix.det()) == a**4*b**2*(b-a)**2,'mixed determinant')
        for p in (2,3,5,7):
            need(tuple(val(x,p) for x in full) == predicted(a,b,p),'entire mixed partition')
            count += 1
    for e in range(2,16):
        zero,one = predicted(2**(e+1),-2**e,2),predicted(2**(e+1),2**e,2)
        for N in range(1,4*e+4):
            diff = sum(min(N,z) for z in zero)-sum(min(N,z) for z in one)
            need(diff == int(3*e+2 <= N <= 4*e),'whole kernel precision window')
    return len(bank),count


def norm(x):
    return min(x % 1,(-x) % 1)


def actual_atlas(s):
    if not 3 <= s <= 356:
        return False
    p = 2
    while p*p <= s:
        e = 0
        while s % p == 0:
            s //= p
            e += 1
        if e and (p % 3 != 2 or e > 2):
            return False
        p += 1
    return s == 1 or s % 3 == 2


def balanced_audit():
    Q = 91**6
    qs = (179,181,183,185,187)
    P = prod(qs)
    V = sorted([P]+[(356-q)*P//q for q in qs])
    U = [1,49,109,331,331**2,331**3,331**4]
    need(reduce(gcd,V)==reduce(gcd,U)==1,'primitive shapes')
    for t in (36883259177,36883259192):
        row = [t*v for v in V]+U
        need(len(set(row)) == 13 and sum(row) <= Q**2,'actual box')
        need(gcd(t,49*109*331)==1,'coprime scale')
        relation_rows = []
        for i,j in combinations(range(13),2):
            D = gcd(row[i],row[j])
            if actual_atlas((row[i]+row[j])//D):
                need((i<6)==(j<6),'no cross decoder edge')
                edge = [0]*13
                edge[i],edge[j] = row[j]//D,-row[i]//D
                relation_rows.append(edge)
        need(Matrix(relation_rows).rank() == 11,'actual decoder rank')
        need(min(t*gcd(v,w)//gcd(t*gcd(v,w),u)
                 for v,w in combinations(V,2) for u in U) > Q,'first mixed orientation')
        need(t*min(V) > Q*(U[-1]+U[-2]),'second mixed orientation')
        for sub in combinations(row,7):
            need(reduce(gcd,sub)==1,'every seven-subset gcd')
        need(min(V)>3*Q//28,'old balanced minimum comparison fails')
        need(56*U[-1]*min(V)>6*Q*(U[-1]+1),'old native comparison fails')
        for u in U[:-1]:
            D = gcd(u,U[-1])
            A,B = u//D,U[-1]//D
            delta = gcd(D,t*min(V))
            R = Q*(A+B)-(A-1)*(B-1)
            need(delta == 1 and D//delta <= Q,'native relation gate is paid')
            need(6*delta*R < 56*U[-1]*min(V),'every maximum-endpoint native gate fails')
        need(6*t<56*U[-1] and 7<49*max(V),'both simple arc comparisons fail')
        phase = F(1,2) if t % 2 else F(11,23)
        clearance = min(norm(x*phase) for x in row)
        need(clearance == (F(1,2) if t % 2 else F(3,23)),'actual safe phase')
        if t % 2 == 0:
            for q in range(2,23):
                for r in range(1,q):
                    need(min(norm(x*F(r,q)) for x in row)<F(1,14),'no weak phase through q22')
            need(all(min(norm(x*F(r,16)) for x in row)==F(1,16)
                     for r in range(1,16,2)),'odd-sixteenth gate hostile')
    return 2,2*1716


def main():
    matrices,partitions = smith_audit()
    rows,subsets = balanced_audit()
    print('SECOND_SESSION_INDEPENDENT_ROOT_AUDIT')
    print(f'mixed_matrices={matrices} full_prime_partitions={partitions} primes=2,3,5,7')
    print('whole_kernel_precision_windows=e2..15 PASS')
    print(f'balanced_actual_rows={rows} independently_checked_seven_subsets={subsets}')
    print('actual_graph_rank_both_crossing_orientations_box_phase_and_failed_gates=PASS')
    print('all_six_maximum_endpoint_native_gates_fail_for_both_balanced_rows=PASS')
    print('PASS')


if __name__ == '__main__':
    main()
