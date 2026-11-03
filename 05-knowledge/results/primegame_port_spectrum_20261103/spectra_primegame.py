"""Reconstruct the actual blue metagraph and compile quartet execution cells."""
from collections import Counter
from fractions import Fraction
from itertools import combinations,permutations
from math import gcd,lcm
from pathlib import Path
import json
import sympy as s
from port_scout import classes,converse,encode,rows
HERE=Path(__file__).resolve().parent
FRACS=[(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
PRIMES=[2,3,5,7,11,13,17,19,23,29]
BLOCKS=['JK','AIMN','CEFH','BDGL']

def chart(mask):
    a=[0]*5
    for i in range(4):a[i]|=1<<(i+1)
    for k,(i,j) in enumerate((i,j) for i in range(5) for j in range(i+2,5)):
        if mask>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def canon_table():
    result={}
    for m,a,_ in classes(5):
        for p in permutations(range(5)):
            b=[sum(((a[p[i]]>>p[j])&1)<<j for j in range(5)) for i in range(5)]
            result[encode(b)]=m
    return result

def next_rule(n):
    for i,(a,b) in enumerate(FRACS):
        if n%b==0:return i,n//b*a
    raise RuntimeError('denominator-one instruction missing')

def run(n,steps):
    out=[]
    for _ in range(steps):i,n=next_rule(n);out.append(chr(65+i))
    return ''.join(out),n

def cell(word):
    prefix=Fraction(1);terms=[];D=1
    for c in word:
        j=ord(c)-65
        ds=[prefix.denominator*b//gcd(prefix.numerator,b) for _,b in FRACS]
        D=lcm(D,ds[j]);terms.extend(ds[:j]);prefix*=Fraction(*FRACS[j])
    forb=sorted({f//gcd(f,D) for f in terms})
    minimal=[f for f in forb if not any(f%g==0 for g in forb if g<f)]
    return {'word':word,'minimum':D,'forbidden_multiplier_divisors':minimal,'nonempty':1 not in minimal,'product':str(prefix)}

def main():
    x=s.symbols('x');ct=canon_table()
    def node(a):return min(ct[encode(a)],ct[encode(converse(a))])
    nodes=sorted({node(chart(m)) for m in range(64)});idx={v:i for i,v in enumerate(nodes)}
    incidence={v:[] for v in nodes};B=s.zeros(10);blue=[];black=[]
    for m in range(64):
        a=chart(m);rho=[sum((1-((a[4-i]>>(4-j))&1))<<j for j in range(5) if i!=j) for i in range(5)]
        fixed=a==rho;incidence[node(a)].append(fixed)
        if m<(m^63):
            u,v=node(a),node(chart(m^63));(blue if fixed else black).append([u,v])
            if fixed:B[idx[u],idx[v]]+=1;B[idx[v],idx[u]]+=1
    pure=[v for v in nodes if all(incidence[v])];mixed=[v for v in nodes if any(incidence[v]) and not all(incidence[v])]
    active=sorted(pure+mixed)
    char=lambda inds:str(s.factor(B.extract([idx[v] for v in inds],[idx[v] for v in inds]).charpoly(x).as_expr()))
    polya=s.expand(x*(x*x+2)*(x*x+x-s.Rational(1,2)))
    polyb=s.expand((x-2)*(x*x+1)*(x*x+x-1))
    cells={b:[cell(''.join(p)) for p in permutations(b)] for b in BLOCKS}
    guards=0
    for records in cells.values():
        for c in records:
            for t in range(1,101):
                expected=all(t%d for d in c['forbidden_multiplier_divisors'])
                observed=run(c['minimum']*t,len(c['word']))[0]==c['word']
                if expected!=observed:raise RuntimeError('priority cylinder')
                guards+=1
    inc=[]
    for a,b in FRACS:
        vec=[]
        for p in PRIMES:
            z=0
            while a%p==0:a//=p;z+=1
            while b%p==0:b//=p;z-=1
            vec.append(z)
        if a!=1 or b!=1:raise RuntimeError('prime support')
        inc.append(vec)
    matrix=[]
    for i,v in enumerate(inc):matrix.append(v+[-int(chr(65+i) in BLOCKS[1]),-int(chr(65+i) in BLOCKS[2])])
    null=s.Matrix(matrix).nullspace()
    if len(null)!=1:raise RuntimeError('grading uniqueness')
    grading=[int(v) for v in null[0]]
    result={'blue_metagraph':{'nodes':nodes,'pure_blue':pure,'mixed':mixed,'pure_black':sorted(set(nodes)-set(active)),'blue_lines':blue,'black_lines':black,'characteristic':char(nodes),'self_converse_characteristic':char(active),'mixed_characteristic':char(mixed),'pure_characteristic':char(pure)},
            'proposed_spectra':{'first_charpoly':str(polya),'second_charpoly':str(polyb),'traces':[-1,1],'first_nonintegral_coefficients':True},
            'quartet_cells':cells,'guard_checks':guards,'grading':grading}
    (HERE/'spectra_primegame.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='quartet_cells'},indent=2))
    print({b:[(v['word'],v['minimum'],v['forbidden_multiplier_divisors']) for v in cs if v['nonempty']] for b,cs in cells.items()})

if __name__=='__main__':main()
