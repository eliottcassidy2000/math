"""A fixed 16-fraction FRACTRAN program for the literal-line count L_n (n>=3).

Input: 306 * 13**(n-3).
Halt: 59 * 2**L_n * 3**(n-1), where L_n=2**(binom(n-1,2)-1).
The output is the exponent of 2, not the whole final integer.
"""
import math

FRACTIONS = (
    (19,221), (59,17),
    (1265,57), (29,19),
    (19,23),
    (93,319), (37,29),
    (29,31),
    (41,185), (51,37),
    (2107,82), (47,41),
    (41,43),
    (106,329), (37,47),
    (47,53),
)
PRIMES = (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59)


def expected(n):
    if type(n) is not int or n<3: raise ValueError('n must be an integer >=3')
    return 2 ** ((n-1)*(n-2)//2-1)


def exact_steps(n):
    t=n-3
    return 6*expected(n)-5+(7*t*t+29*t)//2


def run_integer(n,step_limit=1_000_000):
    expected(n)
    value=306*13**(n-3)
    steps=0; maxbits=value.bit_length()
    while steps<step_limit:
        for numerator,denominator in FRACTIONS:
            if value%denominator==0:
                value=value//denominator*numerator
                break
        else:
            out=0; residual=value
            while residual%2==0: out+=1; residual//=2
            if residual!=59*3**(n-1): raise RuntimeError('bad terminal shape')
            return {'n':n,'output':out,'steps':steps,'max_integer_bits':maxbits}
        steps+=1; maxbits=max(maxbits,value.bit_length())
    raise RuntimeError('literal step limit reached')


def run_exponents(n,step_limit=10_000_000):
    expected(n)
    def factor(x):
        e=[0]*len(PRIMES)
        for j,p in enumerate(PRIMES):
            while x%p==0: e[j]+=1; x//=p
        if x!=1: raise RuntimeError('unknown prime')
        return e
    rows=[(factor(a),factor(b)) for a,b in FRACTIONS]
    exps=[0]*len(PRIMES)
    exps[0]=1;exps[1]=2;exps[5]=n-3;exps[6]=1
    steps=0
    while steps<step_limit:
        for pos,neg in rows:
            if all(a>=b for a,b in zip(exps,neg)):
                exps=[a+b-c for a,b,c in zip(exps,pos,neg)]
                break
        else:
            if exps[1]!=n-1 or exps[-1]!=1 or any(exps[2:-1]):
                raise RuntimeError('bad terminal exponent shape')
            return {'n':n,'output':exps[0],'steps':steps}
        steps+=1
    raise RuntimeError('exponent step limit reached')
