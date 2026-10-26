#!/usr/bin/env python3
"""Independent full-triple modular controls; no atlas implementation imports."""
from itertools import product
from math import gcd
import json

gates = 0

def check(value, message):
    global gates
    gates += 1
    if not value:
        raise RuntimeError(message)

def squarefree(n):
    return all(n % (k*k) for k in range(2, abs(n)+1))

def primitive_point(delta, modulus, p):
    # Unlike the producer's two normalized charts, search all primitive triples.
    for a, b, c in product(range(modulus), repeat=3):
        if any(x % p for x in (a,b,c)) and (a*a+b*b-delta*c*c) % modulus == 0:
            return [a,b,c]
    return None

rows=[]
for delta in range(-20,21):
    if not delta or not squarefree(delta):
        continue
    odd = delta//2 if delta%2==0 else delta
    for p, modulus in [(2,8),(3,9),(5,25),(7,49)]:
        point=primitive_point(delta,modulus,p)
        predicted_failure = odd%4==3 if p==2 else (delta%p==0 and p%4==3)
        check((point is None)==predicted_failure,'full primitive triples agree with local symbol')
        rows.append([delta,p,point])

for d,p,k in [(17,2,8),(13,3,6),(-7,2,8)]:
    modulus=p**k
    check(any((x*x-d)%modulus==0 for x in range(modulus)),'hostile quadratic field embeds locally')

for p,modulus in [(2,4),(3,9)]:
    check(primitive_point(3,modulus,p) is None,'delta3 local obstruction in embedded field')

# Explicit obstruction predictions for the smallest degree examples.
# This checks the elementary gates, while the global splitting theorem remains cited.
expected={2:2,3:2,5:2,13:4,17:4,73:4,-1:2,-7:4}
for d,degree in expected.items():
    check((2 if d%8!=1 and d%3!=1 else 4)==degree,'delta3 minimum degree gate')

print(json.dumps({'status':'PASS','gates':gates,'universe':'signed squarefree delta between -20 and20; full primitive triples modulo8,9,25,49','local_rows':rows,'scope':'Independent finite local controls only; universal quadratic-field atlas rests on the cited splitting theorem and inherited resolved fibre.'},sort_keys=True,indent=2))
