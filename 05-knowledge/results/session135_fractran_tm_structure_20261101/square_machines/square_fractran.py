"""Literal 14-fraction square recognizer; exponent execution is a separate path."""
from math import gcd, isqrt

FRACTIONS = [(43,14),(31,7),(66,43),(65,33),(17,11),(19,26),(37,13),
             (11,19),(23,34),(31,17),(58,23),(123,145),(99,29),(29,41)]
PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43]

def factor(x):
    out={}
    for p in PRIMES:
        while x%p == 0:
            out[p]=out.get(p,0)+1
            x//=p
    if x != 1:
        raise ValueError('undeclared prime')
    return out

FACTORED=[(factor(a),factor(b)) for a,b in FRACTIONS]

def run(n, literal=False, limit=10000000):
    if not isinstance(n,int) or n<0:
        raise ValueError('input exponent must be a nonnegative integer')
    steps=0
    if literal:
        value=7*(1<<n)
        while steps<limit:
            for a,b in FRACTIONS:
                if value%b==0:
                    value=value//b*a
                    steps+=1
                    break
            else:
                return factor(value),steps
    else:
        val={2:n,7:1}
        while steps<limit:
            for a,b in FACTORED:
                if all(val.get(p,0)>=e for p,e in b.items()):
                    for p,e in b.items(): val[p]=val.get(p,0)-e
                    for p,e in a.items(): val[p]=val.get(p,0)+e
                    steps+=1
                    break
            else:
                return {p:e for p,e in val.items() if e},steps
    raise RuntimeError('step limit reached')

def expected_steps(n):
    if n==0: return 1
    k=isqrt(n)
    if k*k==n: return 5*n+2
    return 3*n+2*k*k+4*k+4

if __name__ == '__main__':
    from pathlib import Path
    Path(__file__).with_name('square.fractions.txt').write_text(', '.join(f'{a}/{b}' for a,b in FRACTIONS)+'\n')
