"""Standalone literal checks and one-scalar codec; no repository imports."""
from math import factorial
from functools import cache
from itertools import product
from random import Random


def tt(k):
    return tuple(((1 << k) - 1) ^ ((2 << i) - 1) for i in range(k))


def ordinal(a, b):
    mask = ((1 << len(b)) - 1) << len(a)
    return tuple(x | mask for x in a) + tuple(x << len(a) for x in b)


def wrap(b, k=1):
    # Blocks r -> TT_k -> B -> r, with r at 0.
    return (((1 << (k+1)) - 2),) + tuple(
        x << 1 | (((1 << len(b))-1) << (k+1)) for x in tt(k)
    ) + tuple((x << (k+1)) | 1 for x in b)


def literal_eval(rows):
    n = len(rows)
    allmask = (1 << n) - 1

    @cache
    def starts(s, v):
        r = s ^ (1 << v)
        return sum(starts(r, w) for w in range(n)
                   if ((rows[v] & r) >> w) & 1) if r else 1

    @cache
    def h(s):
        return sum(starts(s, v) for v in range(n) if (s >> v) & 1) if s else 1

    @cache
    def p(t, s=allmask):
        if not t:
            return int(s == 0)
        if t == 1:
            return h(s)
        z, total = s, 0
        while True:
            total += h(z)*p(t-1, s^z)
            if z == 0:
                return total
            z = (z-1) & s
    return p


def wide_value(f, t, k):
    return t**(k+1)*f(t) + t*(
        ((t+1)**k-t**k)*(f(t+1)-f(t)) +
        ((t-1)**k-t**k)*(f(t-1)-f(t)))


def compiler(n):
    # c represents the functional of the outer context, built outside in.
    d = max(1, (n-1)//2)
    c, pads = [0, 1], []
    for s in range(1, d):
        for p in range(2*s+1):
            z = [0]*(s+2)
            for t in range(1, s+1):
                a = t**p*c[t]
                z[t] += t*t*a
                z[t+1] += t*a
                z[t-1] -= t*a
            z[0] = 0
            if all(z[1:]):
                break
        else:
            raise ArithmeticError('no pad found')
        c, pads = z, pads+[p]
    bound = factorial(n)*d**n*sum(map(abs,c))
    big = d*(2*bound).bit_length()
    return c, pads, big


def decode(n, count):
    c, _, k = compiler(n)
    values = []
    for t in range(len(c)-1, 0, -1):
        b = t**k
        q = (2*count+b)//(2*b)
        value, r = divmod(q,c[t])
        if r or value < 1:
            raise ValueError('not an encoded evaluation')
        count -= q*b
        values.append(value)
    if count:
        raise ValueError('nonzero remainder')
    return values[::-1]


def evaluate_circuit(f, pads, big):
    # Does not count the constructed tournament; evaluates its context.
    g = lambda t: t**big*f(t)
    for p in reversed(pads):
        g = (lambda g,p: lambda t: t**p*wide_value(g,t,1) if t else 0)(g,p)
    return g(1)


def literal_circuit(b, pads, big):
    b = ordinal(tt(big),b)
    for p in reversed(pads):
        b = ordinal(tt(p),wrap(b))
    return b


def sources(n):
    e = [(i,j) for i in range(n) for j in range(i+1,n)]
    for flips in product((0,1),repeat=len(e)):
        b = [0]*n
        for (i,j),flip in zip(e,flips):
            if flip:
                i,j = j,i
            b[i] |= 1<<j
        yield tuple(b)


def run():
    samples = wraps = 0
    for n in range(1,5):
        c,pads,big = compiler(n)
        for b in sources(n):
            p = literal_eval(b)
            for k in (1,2,3):
                w = literal_eval(wrap(b,k))
                for t in (1,2,3):
                    assert w(t) == wide_value(p,t,k)
                    wraps += 1
            # A literal small-pad, scalar read, distinct from the general
            # safety bound used for the large codec.
            for tiny in (0,1):
                assert literal_eval(literal_circuit(b,[0],tiny))(1) == evaluate_circuit(p,[0],tiny)
            samples += 1
    random = Random(67*69)
    for n in range(1,51):
        c,pads,big = compiler(n)
        # Synthetic odd-cycle-shaped polynomials, not asserted realizable.
        coeff = {j:random.randrange(1,100) for j in range(n%2 or 2,n,2)}
        while sum(coeff.values())+1 > factorial(n):
            coeff[max(coeff)] //= 2
        p = lambda t: t**n + sum(a*t**j for j,a in coeff.items())
        value = sum(c[t]*t**big*p(t) for t in range(1,len(c)))
        assert decode(n,value) == [p(t) for t in range(1,len(c))]
    # A nontrivial literal circuit, with smaller-than-proof padding.
    b = wrap(tt(3))
    p = literal_eval(b)
    assert literal_eval(literal_circuit(b,[1,0],2))(1) == evaluate_circuit(p,[1,0],2)
    for rows,q in (((22,24,10,1,12),108), ((6,12,24,17,3),120)):
        p = literal_eval(rows)
        assert (p(1),p(2)) == (15,q)
        for k in (1,2,3):
            assert literal_eval(wrap(rows,k))(1) == 30+(2**k-1)*(q-15)
    assert sorted(s*literal_eval(tt(2))(2)+j for s in (72,73,74)
                  for j in (67,69)) == list(range(355,366,2))
    print(f'{samples} literal sources; {wraps} literal wide-wrap comparisons; '
          'small literal assembled circuits; exact codec controls through n=50; '
          'six-count window checked.')


if __name__ == '__main__':
    run()
