"""Prime-digit and compatible-CRT compiler for the actual <2,3> slope group.

Lattice (h,c,v) has basis (h,0),(c,v), with h,v>0 and 0<=c<h.
No floating point, external packages, or search over Collatz starts.
"""
from math import gcd, lcm
from pathlib import Path
import json


def egcd(a, b):
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    d, x, y = egcd(b, a % b)
    return d, y, x - (a // b) * y


def normalize(a, b):
    det = abs(a[0]*b[1] - a[1]*b[0])
    v, s, t = egcd(a[1], b[1])
    if not det or not v:
        raise ValueError('rank-two lattice required')
    h = det // v
    return h, (s*a[0] + t*b[0]) % h, v


def contains(lat, u, w):
    h, c, v = lat
    return w % v == 0 and (u - (w//v)*c) % h == 0


def prime_layer(p):
    powers = {}
    a = 1
    while a not in powers:
        powers[a] = len(powers)
        a = 2*a % p
    h = len(powers)
    w, a = 1, 3 % p
    while a not in powers:
        w += 1
        a = 3*a % p
    return h, (-powers[a]) % h, w


def lift(lat, p, precision):
    h, c, v = lat
    a, b = (h, 0), (c, v)
    modulus = precision*p
    values = [(pow(2,x,modulus)*pow(3,y,modulus)) % modulus for x,y in (a,b)]
    if any((x-1) % precision for x in values):
        raise RuntimeError('old lattice did not reduce to identity')
    lam = [((x-1)//precision) % p for x in values]
    if lam == [0,0]:
        return lat, lam
    if lam[0] == 0:
        a, b = b, a
        lam = lam[::-1]
    k = lam[1]*pow(lam[0],-1,p) % p
    new = normalize((p*a[0],p*a[1]), (b[0]-k*a[0],b[1]-k*a[1]))
    return new, [((x-1)//precision) % p for x in values]


def prime_power(p, e):
    lat = prime_layer(p)
    trace = [{'modulus':p,'lattice':lat}]
    modulus = p
    for _ in range(1,e):
        lat, lam = lift(lat,p,modulus)
        modulus *= p
        trace.append({'modulus':modulus,'lattice':lat,'digit_functional':lam})
    return lat, trace


def intersect(a,b):
    h1,c1,v1 = a
    h2,c2,v2 = b
    base_y = lcm(v1,v2)
    d = gcd(h1,h2)
    delta = c1*(base_y//v1)-c2*(base_y//v2)
    y = base_y*(d//gcd(d,delta))
    x1,x2 = c1*(y//v1),c2*(y//v2)
    k = ((x2-x1)//d)*pow(h1//d,-1,h2//d) % (h2//d) if h2//d>1 else 0
    h = lcm(h1,h2)
    return h,(x1+h1*k)%h,y


def sum_invariants(a,b):
    h1,c1,v1 = a
    h2,c2,v2 = b
    vectors = [(h1,0),(c1,v1),(h2,0),(c2,v2)]
    content = gcd(*(x for row in vectors for x in row))
    index = gcd(*(abs(x[0]*y[1]-x[1]*y[0]) for i,x in enumerate(vectors) for y in vectors[i+1:]))
    return content,index//content


def factor(n):
    factors=[]
    p=2
    while p*p<=n:
        if n%p==0:
            e=0
            while n%p==0:
                n//=p;e+=1
            factors.append((p,e))
        p+=1
    if n>1:factors.append((n,1))
    return factors


def compile_modulus(n):
    if n<1 or gcd(n,6)!=1:
        raise ValueError('use the ternary layer separately')
    lat=(1,0,1)
    trace=[]
    for p,e in factor(n):
        local,steps=prime_power(p,e)
        compatibility=sum_invariants(lat,local)
        lat=intersect(lat,local)
        trace.append(dict(prime=p,exponent=e,local_steps=steps,
                          compatibility=compatibility,merged_lattice=lat))
    return lat,trace


def subgroup(n):
    seen={1%n};stack=list(seen)
    while stack:
        x=stack.pop()
        for a in(2,3):
            y=a*x%n
            if y not in seen:
                seen.add(y);stack.append(y)
    return seen


def slope_image(n,lat):
    h,c,v=lat
    return {pow(2,u,n)*pow(3,w,n)%n for u in range(h) for w in range(v)}


def main():
    count=0
    for n in range(5,200,2):
        if gcd(n,6)!=1:continue
        lat,_=compile_modulus(n)
        actual=subgroup(n)
        if len(actual)!=lat[0]*lat[2] or slope_image(n,lat)!=actual:
            raise RuntimeError(('subgroup',n,lat))
        for u in range(-8,9):
            for v in range(-8,9):
                if contains(lat,u,v)!=(pow(2,u,n)*pow(3,v,n)%n==1):
                    raise RuntimeError(('relation',n,u,v,lat))
                count+=1
    towers=[]
    for p,e in((5,4),(7,3),(11,3),(13,3),(1093,2)):
        lat,trace=prime_power(p,e)
        n=p**e
        for u in range(-12,13):
            for v in range(-12,13):
                if contains(lat,u,v)!=(pow(2,u,n)*pow(3,v,n)%n==1):
                    raise RuntimeError(('tower relation',p,e,u,v))
        towers.append(dict(prime=p,steps=trace))
    examples=[]
    for n in(35,65,95,119,143,145,1225,1729):
        lat,trace=compile_modulus(n)
        examples.append(dict(modulus=n,lattice=lat,order=lat[0]*lat[2],trace=trace))
    # All fixed affine shortcut words through length six, both signs, on
    # every odd modulus through45: ternary image/fibre rule versus literal map.
    word_checks=0
    for n in range(3,46,2):
        for length in range(1,7):
            for mask in range(1<<length):
                ones=mask.bit_count();carry=0
                for j in range(length):
                    if mask>>j&1:carry=3*carry+(1<<j)
                d=gcd(3**ones,n)
                for sign in(-1,1):
                    outputs={}
                    for x in range(n):
                        y=x
                        for j in range(length):
                            y=((3*y+sign) if mask>>j&1 else y)*pow(2,-1,n)%n
                        outputs[y]=outputs.get(y,0)+1
                    target={y for y in range(n) if ((1<<length)*y-sign*carry)%d==0}
                    if set(outputs)!=target or set(outputs.values())!={d}:
                        raise RuntimeError(('ternary fibre',n,length,mask,sign))
                    word_checks+=1
    report=dict(status='FINITE-EXACT controls; all-modulus proofs in PROOF.md',
                universe='all admissible odd moduli5..199, exponent pairs[-8,8]^2; five named prime-power towers; all binary words lengths1..6 on odd moduli3..45 and both signs',
                relation_checks=count,word_fibre_checks=word_checks,towers=towers,examples=examples)
    Path(__file__).with_name('receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(relation_checks=count,word_fibre_checks=word_checks,
                         examples=[{k:v for k,v in x.items() if k!='trace'} for x in examples]),indent=2))


if __name__=='__main__':main()
