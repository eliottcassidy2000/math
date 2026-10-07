"""Exact frieze/Collatz interfaces: marked coordinates, guards, and obstructions.

Standard library only. No ROOT lookup, bounded orbit guess, or imported tests.
Run with Python normally and with -B -O; checks do not use assert.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product, combinations

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def word(w):
    need(type(w) is tuple and all(type(a) is int and a > 0 for a in w),
         'exact tuple of positive valuations')


def natural(n):
    need(type(n) is int and n >= 0, 'exact natural integer')


def odd(n):
    need(type(n) is int and n > 0 and n % 2, 'exact positive odd source')


def carrier(w):
    word(w)
    p = q = 1
    b = 0
    for a in w:
        p, q, b = 3*p, q << a, 3*b+q
    return p, q, b


def prefixes(w):
    word(w)
    out = [(1, 1, 0)]
    p = q = 1
    b = 0
    for a in w:
        p, q, b = 3*p, q << a, 3*b+q
        out.append((p, q, b))
    return tuple(out)


def step(n):
    odd(n)
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def replay(n, w):
    odd(n)
    word(w)
    out = [n]
    for a in w:
        need(n != 1, 'no formal ROOT padding in exported receipts')
        n, actual = step(n)
        need(actual == a, 'wrong source valuation')
        out.append(n)
    return tuple(out)


def z_chart(w):
    return tuple(F(b, p) for p, q, b in prefixes(w))


def log2_power(n):
    need(type(n) is int and n > 0 and n & (n-1) == 0, 'exact power of two')
    return n.bit_length()-1


def decode_chart(z, terminal_q):
    need(type(z) is tuple and z and
         all(type(x) in (int, F) for x in z), 'exact rational marked chart')
    need(z[0] == 0, 'fixed translation gauge')
    total = log2_power(terminal_q)
    if len(z) == 1:
        need(total == 0, 'empty chart terminal cost')
        return ()
    gaps = tuple(F(z[i+1]-z[i]) for i in range(len(z)-1))
    need(gaps[0] == F(1, 3) and all(d > 0 for d in gaps), 'initial gap and order')
    w = []
    for i in range(1, len(gaps)):
        ratio = 3*gaps[i]/gaps[i-1]
        need(ratio.denominator == 1, 'integral gap ratio')
        a = log2_power(ratio.numerator)
        need(a >= 1, 'positive valuation')
        w.append(a)
    last = total-sum(w)
    need(last >= 1, 'terminal valuation is retained')
    w.append(last)
    w = tuple(w)
    need(z_chart(w) == z, 'full marked chart authentication')
    return w


def matmul(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in (0, 1))
                       for j in (0, 1)) for i in (0, 1))


def frieze_matrix(w):
    word(w)
    m = ((1, 0), (0, 1))
    for a in w:
        m = matmul(m, ((a, -1), (1, 0)))
    return m


def insert_ear(w, i):
    word(w)
    need(type(i) is int and 0 <= i < len(w)-1, 'interior edge index')
    return w[:i]+(w[i]+1, 1, w[i+1]+1)+w[i+2:]


def is_quiddity(w):
    """Exact constructive recognizer: delete polygon ears down to the triangle."""
    word(w)
    q = list(w)
    while len(q) > 3:
        if 1 not in q:
            return False
        i = q.index(1)
        left, right = (i-1) % len(q), (i+1) % len(q)
        if q[left] <= 1 or q[right] <= 1:
            return False
        q[left] -= 1
        q[right] -= 1
        q.pop(i)
    return q == [1, 1, 1]


@dataclass(frozen=True)
class Descent:
    word: tuple
    first: int
    period: int
    endpoint: int
    increment: int


def descent_row(w):
    word(w)
    need(is_quiddity(w) and len(w) >= 5, 'contracting positive quiddity')
    p, q, b = carrier(w)
    need(p < q, 'strict slope contraction')
    cut = max([1, b//(q-p)]+[(qi-bi)//pi for pi, qi, bi in prefixes(w)[:-1]])
    residue = ((q-b)*pow(p, -1, 2*q)) % (2*q)
    first = residue + 2*q*max(0, (cut-residue)//(2*q)+1)
    endpoint = (p*first+b)//q
    return Descent(w, first, 2*q, endpoint, 2*p)


def apply_descent(row, t):
    need(type(row) is Descent and type(row.word) is tuple and
         all(type(x) is int for x in (row.first, row.period, row.endpoint, row.increment)),
         'exact descent record')
    natural(t)
    need(row == descent_row(row.word), 'authenticated descent record')
    n, y = row.first+row.period*t, row.endpoint+row.increment*t
    need(replay(n, row.word)[-1] == y and y < n, 'actual paid endpoint')
    return n, y


def endpoint_label(w):
    p, q, b = carrier(w)
    return p, (b*pow(q, -1, p)) % p


def compatible(w, v):
    p, c = endpoint_label(w)
    q, d = endpoint_label(v)
    return (c-d) % min(p, q) == 0


def join_row(w, v):
    """Forward t>=0 part of an exact simultaneous endpoint progression.

    Source w must grow faster than child v as the common endpoint increases.
    A strict finite cut pays positivity, no early ROOT, and child<source.
    """
    word(w)
    word(v)
    need(compatible(w, v), 'incompatible ternary endpoint labels')
    p, q, b = carrier(w)
    pp, qq, bb = carrier(v)
    slope = F(q, p)-F(qq, pp)
    need(slope > 0, 'designated child orientation')
    pm = max(p, pp)
    label = endpoint_label(w if p >= pp else v)[1]
    residue = label if label % 2 else label+pm
    cuts = [F(0), (F(b, p)-F(bb, pp))/slope]
    for ww, p0, q0, b0 in ((w, p, q, b), (v, pp, qq, bb)):
        for pi, qi, bi in prefixes(ww)[:-1]:
            cuts.append(F(p0*(qi-bi)+pi*b0, pi*q0))
        cuts.append(F(b0, q0))
    cut = max(x.numerator//x.denominator for x in cuts)
    end = residue+2*pm*max(0, (cut-residue)//(2*pm)+1)
    n, h = (q*end-b)//p, (qq*end-bb)//pp
    need(n > h > 0, 'positive smaller child')
    return n, 2*q*(pm//p), h, 2*qq*(pm//pp), end, 2*pm


def join_instance(w, v, t):
    natural(t)
    n, dn, h, dh, y, dy = join_row(w, v)
    n, h, y = n+dn*t, h+dh*t, y+dy*t
    need(replay(n, w)[-1] == replay(h, v)[-1] == y and h < n,
         'source-authenticated actual common future')
    return n, h, y


def pentagon(z):
    need(len(z) == 5, 'five marked vertices')
    edge = lambda i, j: z[max(i, j)]-z[min(i, j)]
    data = {(i, i+1): edge(i, i+1) for i in range(4)}
    data[(0, 4)] = edge(0, 4)
    data[(0, 2)], data[(0, 3)] = edge(0, 2), edge(0, 3)
    for i, j, k, l, old, new in (
        (0, 1, 2, 3, (0, 2), (1, 3)),
        (0, 1, 3, 4, (0, 3), (1, 4)),
        (1, 2, 3, 4, (1, 3), (2, 4)),
        (0, 1, 2, 4, (1, 4), (0, 2)),
        (0, 2, 3, 4, (2, 4), (0, 3)),
    ):
        value = (data[(i, j)]*data[(k, l)]+data[(i, l)]*data[(j, k)])/data[old]
        del data[old]
        data[new] = value
        need(value == edge(*new), 'exact Pluecker flip')
    return data


def rejected(f):
    try:
        f()
    except (ValueError, TypeError):
        return
    raise ValueError('hostile was accepted')


def main():
    charts = flips = 0
    for r in range(6):
        for w in product(range(1, 5), repeat=r):
            z = z_chart(w)
            p, q, b = carrier(w)
            need(decode_chart(z, q) == w, 'lossless marked chart')
            for i, j, k, l in combinations(range(r+1), 4):
                need((z[k]-z[i])*(z[l]-z[j]) ==
                     (z[j]-z[i])*(z[l]-z[k])+(z[l]-z[i])*(z[k]-z[j]), 'Pluecker')
            if r == 4:
                pentagon(z)
                flips += 5
            charts += 1
    need(z_chart((1,)) == z_chart((2,)), 'terminal cost hostile')
    need(replay(3, (1,)) == (3, 5) and replay(9, (2,)) == (9, 7),
         'same chart opposite drift with the exact declared valuations')
    rejected(lambda: decode_chart((F(0), F(1, 3), F(2, 3)), 8))
    need(frieze_matrix((1,)*9) == ((-1, 0), (0, -1)) and
         not is_quiddity((1,)*9), 'negative identity alone is not a positive polygon frieze')

    ear_controls = 0
    for a, b in product(range(1, 17), repeat=2):
        u, v = (a, b), (a+1, 1, b+1)
        need(frieze_matrix(u) == frieze_matrix(v), 'ear transfer identity')
        need(not compatible(u, v), 'single ear ternary wall')
        ear_controls += 1
    descendants = {(1, 1)}
    rows = []
    first_compatible = []
    first_reset = []
    quiddities = 0
    for block_len in range(2, 10):
        matches = []
        reset_matches = []
        for block in sorted(descendants):
            need(frieze_matrix(block) == frieze_matrix((1, 1)), 'ear-class invariant')
            w = (1,)+block
            need(is_quiddity(w) and frieze_matrix(w) == ((-1, 0), (0, -1)), 'closed quiddity')
            need(sum(w) == 3*len(w)-6, 'triangle incidence budget')
            if len(w) > 3:
                need(all(not (w[i] == w[(i+1) % len(w)] == 1) for i in range(len(w))), 'no adjacent ears')
            if len(w) >= 5:
                row = descent_row(w)
                for t in (0, 1, 7):
                    apply_descent(row, t)
            if block != (1, 1) and compatible(block, (1, 1)):
                matches.append(block)
            if block != (1, 1) and compatible(w, (1, 1, 1)):
                reset_matches.append(w)
            quiddities += 1
        if matches and not first_compatible:
            first_compatible = matches
        if reset_matches and not first_reset:
            first_reset = reset_matches
        rows.append((block_len, len(descendants), len(matches), len(reset_matches)))
        descendants = {insert_ear(w, i) for w in descendants for i in range(len(w)-1)}
    need(first_compatible == [(2, 2, 2, 3, 1, 2, 5)], 'first compatible block in finite universe')
    w = (1, 2, 2, 2, 2, 2, 2, 1, 7)
    need(first_reset == [w], 'first compatible marked reset2 in finite universe')
    need(join_row(w, (1, 1, 1)) == (2583211, 4194304, 7183, 11664, 24245, 39366), 'join constants')
    for t in range(128):
        n, h, y = join_instance(w, (1, 1, 1), t)
        need(replay(n, w)[3] < n, 'join is already ordinary three-step descent')
    u = first_compatible[0]
    need(join_row(u, (1, 1)) == (66817, 262144, 495, 1944, 1115, 4374), 'short join constants')
    for t in range(32):
        join_instance(u, (1, 1), t)

    for k in range(512):
        n = 8*k+7
        y, a = step(n)
        z, b = step(y)
        need((a, b) == (1, 1), 'whole excluded dyadic cell')
    for r in range(1, 33):
        source = (1 << (r+1))-1
        need(replay(source, (1,)*r)[-1] == 2*3**r-1, 'variable one-run boundary')

    base = descent_row(w)
    for bad in (True, 1.0, 0, 2, -1):
        rejected(lambda bad=bad: replay(bad, ()))
    for bad in (True, 1.0, -1):
        rejected(lambda bad=bad: apply_descent(base, bad))
    rejected(lambda: apply_descent(replace(base, first=True), 0))
    rejected(lambda: apply_descent(replace(base, endpoint=base.endpoint+2), 0))
    rejected(lambda: decode_chart((F(0), F(1, 3)), True))
    rejected(lambda: descent_row((1, 1, 1)))
    rejected(lambda: join_row((2, 1, 2), (1, 1)))
    rejected(lambda: replay(1, (2,)))

    print('Frieze / first-reset2 exact guarded interfaces')
    print('Marked charts: lengths 0..5, alphabet 1..4:', charts)
    print('Pentagon coordinate flips:', flips)
    print('Single-ear endpoint obstructions: a,b=1..16:', ear_controls)
    print('Ear descendants of block11; columns: length,total,compatible11,compatible111 after prefix1')
    for row in rows:
        print(*row)
    print('Constructively recognized marked quiddities:', quiddities)
    print('Compatible reset2 join:', join_row(w, (1, 1, 1)))
    print('No added first-descent coverage: its third actual endpoint is already smaller.')
    print('All-height direct-quiddity obstruction: every n=7 mod8.')
    print('No assertion that a smaller-child certificate is a completed ROOT certificate.')
    print('Exact checks:', CHECKS)


if __name__ == '__main__':
    main()
