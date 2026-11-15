"""Check the explicit short-loop generators and the dyadic phase kicks."""
from math import gcd


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def run(x, q, w, m):
    for b in w:
        if b == 0:
            q = (q + (x + (x & 1)) // 2) % m
        x = (2*x + b) % (2*m)
    return x, q


def digits(x, l):
    return [int(i) for i in f'{x:0{l}b}']


def loops(m):
    s = 2*m
    l = (s-1).bit_length()
    fwd = [digits(v, l) for v in range(s)]
    back = [digits((-2**l*v) % s, l) for v in range(s)]
    for v in range(s):
        yield fwd[v] + back[v]
        for b in (0, 1):
            w = (2*v+b) % s
            yield fwd[v] + [b] + back[w]


def main():
    tested = 0
    for m in range(1, 257):
        g = m
        for w in loops(m):
            x, drift = run(0, 0, w, m)
            need(x == 0, (m, w))
            g = gcd(g, drift)
            tested += 1
        need(g == 1, m)
    count = 0
    for a in range(1, 31):
        m = 2**a
        for j in range(a):
            w = [0] + [1]*(j+1) + [0] + [1]*(a+1)
            need(run(2*m-1, 0, w, m) == (2*m-1, (-2**j) % m), (a, j))
            count += 1
    print(f'{tested} short loops, all moduli 1..256; {count} dyadic kicks, exponents 1..30')


if __name__ == '__main__':
    main()
