"""Direct 14-row tail / 20-row full reduced-core-floor counters.

The recurrence and first-eligible program equivalence are unconditional. The
identification with minimum tournament H uses the separately audited core theorem.
"""
from compiler import require

TAIL14 = (
    (13,55), (9728,11),
    (17,65), (622592,13),
    (77,85), (637534208,17),
    (23,133), (41,19),
    (783,46), (31,23),
    (23,29),
    (74,93), (19,31),
    (31,37),
)
FULL20 = ((47,215), (328,43),
          (53,235), (1312,47),
          (11,265), (20992,53)) + TAIL14
PRIMES = (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53)


def expected(n):
    require(type(n) is int and n >= 3, 'core order must be an integer >=3')
    if n <= 5:
        return (3,5,9)[n-3]
    k, r = divmod(n-6,3)
    return (9,15,25)[r] * 3**k


def exact_steps(n, tail=False):
    require(type(n) is int and n >= (6 if tail else 3), 'invalid core input')
    if n <= 5:
        return n-2
    k, r = divmod(n-6,3)
    return 4*(expected(n)-(9,15,25)[r])+6*k+r+2+(0 if tail else 3)


def input_value(n, tail=False):
    require(type(n) is int and n >= (6 if tail else 3), 'invalid core input')
    return (11 if tail else 43)*5**(n-(6 if tail else 3))


def run_integer(n, tail=False, limit=1_000_000):
    rows = TAIL14 if tail else FULL20
    value = input_value(n,tail)
    steps = 0
    while steps < limit:
        for a,b in rows:
            if value % b == 0:
                value = value//b*a
                break
        else:
            out = 0
            while value % 2 == 0:
                out += 1
                value //= 2
            require(value == 41, 'core terminal must be exactly 41*2^output')
            return dict(n=n, output=out, steps=steps)
        steps += 1
    raise RuntimeError('core integer step limit')


def run_exponents(n, tail=False, limit=2_000_000):
    rows = TAIL14 if tail else FULL20
    require(type(n) is int and n >= (6 if tail else 3), 'invalid core input')
    def factor(value):
        out = {}
        for p in PRIMES:
            while value%p == 0:
                out[p] = out.get(p,0)+1
                value //= p
        require(value==1,'unknown core prime')
        return out
    factored = [(factor(a),factor(b)) for a,b in rows]
    exps = {p:0 for p in PRIMES}
    exps[11 if tail else 43] = 1
    exps[5] = n-(6 if tail else 3)
    steps = 0
    while steps < limit:
        for pos,neg in factored:
            if all(exps[p]>=e for p,e in neg.items()):
                for p,e in neg.items(): exps[p]-=e
                for p,e in pos.items(): exps[p]+=e
                break
        else:
            require(exps[41]==1 and all(e==0 for p,e in exps.items() if p not in (2,41)),
                    'pure core terminal invariant')
            return dict(n=n, output=exps[2], steps=steps)
        steps += 1
    raise RuntimeError('core exponent step limit')
