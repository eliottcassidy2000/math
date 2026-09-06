#!/usr/bin/env python3
"""Independent rational certificate audit for the cubic first-return family.

Reads only literal candidate factor arrays using ast.literal_eval. Imports no
producer or CAS. Twenty-five rational specializations certify identities of
proved degree at most 24; see the companion for the degree argument.
"""
import ast
import hashlib
import json
from fractions import Fraction as F
from math import factorial, gcd
from pathlib import Path

CHECKS = 0
ROOT = Path(__file__).resolve().parents[1]
PRODUCER = ROOT/'04-computation/trinomial_cubic_empty_core_sep06.py'
PIN = '7ebdffa7cd529cc1d641dc319954c23a594cff5b4e98ae538e21ed18849b77f5'


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def product(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def falling(n, k):
    return product(range(n-k+1, n+1))


def horner(coefficients, value):
    answer = 0
    for coefficient in coefficients:
        answer = answer*value+coefficient
    return answer


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, s):
    return tuple(x*s for x in a)


def times_tau(a, u):
    # Independent companion recurrence, all coordinates in basis (1,tau,tau^2).
    return (a[2]*F(-u*(u-1)*(u-2), 72),
            a[0]-a[2]*F(7*u*(u-1), 6),
            a[1]-7*u*a[2])


def trace(a, u):
    # Trace of the multiplication operator, by its three diagonal entries.
    t_a = times_tau(a, u)
    tt_a = times_tau(t_a, u)
    return a[0]+t_a[1]+tt_a[2]


def determinant(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    return sum((-1)**j*matrix[0][j]*determinant(
        [row[:j]+row[j+1:] for row in matrix[1:]]) for j in range(n))


def row(a, b, c, mass):
    answer = []
    for z in range(mass+1):
        for y in range(mass-z+1):
            x = mass-y-z
            if -a*x+b*y+c*z == 0:
                answer.append(((x,y,z), factorial(mass)//
                               (factorial(x)*factorial(y)*factorial(z))))
    return answer


def main():
    source = PRODUCER.read_bytes()
    need(hashlib.sha256(source).hexdigest() == PIN, 'literal coefficient bank source pin')
    data = {}
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in ('SHIFTED','DENOMINATORS'):
                    data[target.id] = ast.literal_eval(node.value)
    shifted, denominators = data['SHIFTED'], data['DENOMINATORS']
    need({key:len(values)-1 for key,values in shifted.items()} ==
         {'F5':5,'F11':11,'F3':3,'F15':15}, 'complete four factors')
    need(sum(map(len,shifted.values())) == 38 and
         all(v > 0 for values in shifted.values() for v in values), '38 positive coefficients')
    need(len(denominators) == 3 and all(d > 0 for d in denominators), 'positive denominators')
    manifest = []
    for g in range(11,36):
        u = g-7
        K = F(falling(2*g,14),128)
        q = [F(falling(2*g,21-j), factorial(21-3*j)*factorial(2*j))/K
             for j in range(8)]
        canceled = F(1024*product(2*g-k for k in (15,17,19,20)),factorial(21))
        need(q[0]/(u*(u-1)*(u-2)) == canceled, 'inverse denominator cancellation')
        need(all(q[j] == F(128*product(2*g-k for k in range(14,21-j)),
                            factorial(21-3*j)*factorial(2*j)) for j in range(1,8)),
             'normalized complete higher channels')
        inverse = (F(-84,u-2), F(-504,(u-1)*(u-2)), F(-72,u*(u-1)*(u-2)))
        need(times_tau(inverse,u) == (1,0,0), 'inverse carry identity')
        r = scale(inverse,q[0])
        power = (F(1),F(0),F(0))
        for j in range(1,8):
            r = add(r,scale(power,q[j]))
            power = times_tau(power,u)
        expected = (
            F(-(g-8)*(g-7)*horner([981739413,-29190405089,323857306851,-1589986674295,2916051235200],g),28510570408320000),
            F(-(g-7)*horner([1715049346,-49273046168,530071519257,-2530936512265,4525872382400],g),593970216840000),
            F(-horner([70405697381,-1950797880338,20268905917687,-93593992728910,162061152680400],g),4157791517880000))
        need(r == expected, 'complete rational recurrence residue')
        weighted_powers = [r]
        for _ in range(4):
            weighted_powers.append(times_tau(weighted_powers[-1],u))
        H = [[trace(weighted_powers[i+j],u) for j in range(3)] for i in range(3)]
        minors = [(-1)**k*determinant([r[:k] for r in H[:k]]) for k in range(1,4)]
        factors = {name:horner(values,g-11) for name,values in shifted.items()}
        targets = [F((g-7)*factors['F5'],denominators[0]),
                   F((g-8)*(g-7)**2*factors['F11'],denominators[1]),
                   F((g-8)**2*(g-7)**4*factors['F3']*factors['F15'],denominators[2])]
        need(minors == targets, 'all three exact leading-minor identities')
        need(all(value > 0 for value in minors), 'signed minors positive')
        a,b,c,d = 72,504*u,84*u*(u-1),u*(u-1)*(u-2)
        discriminant = b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
        need(discriminant == 15552*(g-8)*(g-7)**2*factors['F3'], 'discriminant polynomial identity')
        if g == 11:
            need([minors[k-1]*K**k for k in range(1,4)] ==
                 [F(752350432547692),F(21236578848830718128306804,3),
                  F(942083106929885721679784497035400,27)], 'inherited unnormalized minors')
        manifest.append([g,[str(value) for value in r],[str(value) for value in minors]])
    for g in (11,13,16,22):
        need(gcd(g,21)==1,'primitive named row')
        first,second = row(21,2*g-21,3*g-21,g),row(21,2*g-21,3*g-21,2*g)
        need([v for v,_ in first] == [(g-10+j,9-3*j,1+2*j) for j in range(4)],'all first channels')
        need([v for v,_ in second] == [(2*g-21+j,21-3*j,2*j) for j in range(8)],'all second channels')
        u = g-7
        content = F(falling(g,7),factorial(9))
        need([F(weight)/content for _,weight in first] ==
             [u*(u-1)*(u-2),84*u*(u-1),504*u,72], 'literal first normalization')
        need([weight for _,weight in second] ==
             [F(falling(2*g,21-j),factorial(21-3*j)*factorial(2*j)) for j in range(8)],
             'literal complete second coefficients')
        need(all(not row(21,2*g-21,3*g-21,m) for m in range(1,g)), 'no earlier named support returns')
    need(next(m for m in range(1,13) if row(21,3,15,m)) == 4, 'nonprimitive first-return hostile')
    print('PASS: independent stdlib Fraction recurrence and companion traces')
    print('identity_parameters 11..35 inclusive; 25 distinct exact specializations')
    print('proved_degree_bounds residue_coefficients [6,5,4]; signed_minors [6,14,24]')
    print('certificate_bank four factors, 38 positive coefficients; original producer hash',PIN)
    print('named_literal_rows g=11,13,16,22; nonprimitive g=12 has first return 4')
    print('PROVED via degree argument: three signed Hermite minors positive for every real g>=11')
    print('first-return interpretation: integer g>=11, gcd(g,21)=1 only')
    print('semantic_sha256',hashlib.sha256(json.dumps(manifest,separators=(',',':')).encode()).hexdigest())
    print('explicit_checks',CHECKS)


if __name__ == '__main__':
    main()
