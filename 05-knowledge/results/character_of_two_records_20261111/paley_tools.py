"""Paley tournaments on GF(q), q = 3 (mod 4) prime or q = 27: out-neighbour bitmasks (a -> b iff b - a is a non-zero square)."""
def field(q):
    if q == 27:
        def dig(a): return [a % 3, (a // 3) % 3, (a // 9) % 3]
        def und(d): return d[0] + 3 * d[1] + 9 * d[2]
        add = lambda a, b: und([(x + y) % 3 for x, y in zip(dig(a), dig(b))])
        neg = lambda a: und([(-x) % 3 for x in dig(a)])
        def mul(a, b):
            A, B = dig(a), dig(b); c = [0] * 5
            for i in range(3):
                for j in range(3): c[i + j] = (c[i + j] + A[i] * B[j]) % 3
            c[2] = (c[2] + c[4]) % 3; c[1] = (c[1] + c[4]) % 3     # x^4 = x^2 + x
            c[1] = (c[1] + c[3]) % 3; c[0] = (c[0] + c[3]) % 3     # x^3 = x + 1
            return und(c[:3])
    else:
        add = lambda a, b: (a + b) % q
        neg = lambda a: (-a) % q
        mul = lambda a, b: (a * b) % q
    return add, neg, mul
def paley_out(q):
    add, neg, mul = field(q)
    sq = {mul(a, a) for a in range(1, q)}
    assert len(sq) == (q - 1) // 2
    out = [0] * q
    for a in range(q):
        for b in range(q):
            if a != b and add(b, neg(a)) in sq: out[a] |= 1 << b
    return out
