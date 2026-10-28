"""Independent small controls: no imports from the producer."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb, isqrt, prod
import argparse
import json
from pathlib import Path


def need(p, label):
    if not p: raise RuntimeError(label)


def parity_queries():
    def target(x):
        v=0
        for j in range(7):
            if x>>j&1: v ^= j+1
        return v
    depths={}
    for kind in ('one_check','square'):
        value=lambda x: (target(x)&1) if kind=='one_check' else (target(x)==0)
        @lru_cache(None)
        def depth(states):
            if len({value(x) for x in states})<=1: return 0
            choices=[]
            for j in range(7):
                a=tuple(x for x in states if not (x>>j&1))
                b=tuple(x for x in states if x>>j&1)
                if a and b: choices.append(1+max(depth(a),depth(b)))
            return min(choices)
        plain=depth(tuple(range(128)))
        known=max(depth(tuple(x for x in range(128) if x.bit_count()%2==p)) for p in (0,1))
        need((plain,known)==((4,3) if kind=='one_check' else (7,6)), 'adaptive query optimum')
        depths[kind]=dict(no_parity=plain, known_parity=known)
    return depths


def dual_enumerator(r):
    n=2**r-1; half=2**(r-1)
    # Character averaging, independent of the subset producer.
    result={}
    for k in range(n+1):
        signed=sum((-1)**j*comb(half,j)*comb(n-half,k-j)
                   for j in range(max(0,k-(n-half)),min(half,k)+1))
        numerator=comb(n,k)+(2**r-1)*signed
        need(numerator%(2**r)==0, 'MacWilliams divisibility')
        if numerator: result[k]=numerator//(2**r)
    return result


def six_labels():
    def matchings(vertices):
        if not vertices: return [()]
        a=vertices[0]; out=[]
        for b in vertices[1:]:
            rest=[v for v in vertices if v!=a and v!=b]
            out.extend((((a,b),)+tail) for tail in matchings(rest))
        return out
    match=matchings(list(range(6)))
    need(len(match)==15, 'six-label matchings')
    facts=[m for m in combinations(match,5) if len({e for row in m for e in row})==15]
    need(len(facts)==6, 'six-label factorizations')
    # Four even-subset classes e_i={i,5}; quotient by complement.
    chart={}
    for x in range(16):
        mask=0
        for i in range(4):
            if x>>i&1: mask ^= (1<<i)|(1<<5)
        chart[min(mask,mask^63)]=x
    need(len(chart)==16,'even-subset basis')
    triads=[]
    for m in match:
        row=[chart[min((1<<a)|(1<<b), ((1<<a)|(1<<b))^63)] for a,b in m]
        need(row[0]^row[1]^row[2]==0, 'matching zero sum')
        triads.append(row)
    return dict(matchings=len(match), factorizations=len(facts), triples=triads)


def trace_lift():
    out={}
    for m,poly in [(1,7),(2,19),(3,67)]:
        r=2*m; size=2**r
        def mult(x,y):
            z=0
            for i in range(r):
                if y>>i&1: z ^= x<<i
            for i in range(2*r-2,r-1,-1):
                if z>>i&1: z ^= poly<<(i-r)
            return z
        def power(a,n):
            z=1
            for _ in range(n): z=mult(z,a)
            return z
        def tr(a):
            z=0
            for _ in range(r): z ^= a; a=mult(a,a)
            need(z in (0,1), 'absolute trace lands in F2')
            return z
        sub=[a for a in range(size) if power(a,2**m)==a]
        need(len(sub)==2**m,'subfield size')
        compatible=[]
        for c in sub[1:]:
            def b(x,y): return tr(mult(c,mult(x,power(y,2**m))))
            need(all(b(x,x)==0 for x in range(size)), 'alternating')
            need(all(any(b(x,1<<j) for j in range(r)) for x in range(1,size)), 'nondegenerate')
            need(all(b(x,mult(a,x))==0 for x in range(size) for a in sub),'isotropic subfield lines')
            if all(b(mult(x,x),mult(y,y))==b(x,y) for x in range(size) for y in range(size)):
                compatible.append(c)
        need(compatible==[1], 'unique Frobenius invariant trace form')
        out[m]=dict(subfield_size=len(sub), candidate_forms=len(sub)-1, Frobenius_invariant_coefficients=compatible)
    return out


def main():
    here=Path(__file__).resolve().parent
    receipt=json.loads((here/'normal.json').read_text())
    for r in (3,4):
        row=receipt['arithmetic_codes'][str(r)]
        need({int(k):v for k,v in row['weight_enumerator'].items()}==dual_enumerator(r),'dual weight enumerator')
        for vector,d,p,t in row['bank']:
            # Sieve on [0,p], instead of the producer's per-candidate trial division.
            sieve=bytearray(b'\x01')*(p+1); sieve[:2]=b'\x00\x00'
            for q in range(2,isqrt(p)+1):
                if sieve[q]:
                    for multiple in range(q*q,p+1,q): sieve[multiple]=0
            need(bool(sieve[p]) and p+1==d*t*t, 'prime squareclass witness')
            need(all(not sieve[d*u*u-1] for u in range(1,t)), 'least prime for fixed class')
    result=dict(adaptive_queries=parity_queries(), matchings=six_labels(), trace_reconstruction=trace_lift(),
                Hamming_weight_enumerators={r:dual_enumerator(r) for r in (3,4)}, prime_bank_independent_sieve='PASS')
    parser=argparse.ArgumentParser(); parser.add_argument('--output'); args=parser.parse_args()
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')


if __name__=='__main__': main()
