"""Certificate: counts generate (Z/M)^x for EVERY M whose prime factors are <= P (any exponents).
Reduction: the kernel of (Z/M)^x -> (Z/M')^x, M' = 8 * prod_{odd p|M} p^2, lies in the Frattini subgroup, so it is
enough that the counts generate the Frattini quotient of (Z/8)^x x prod_{3<=p<=P} (Z/p^2)^x, i.e. for every prime l the
matrix of l-coordinates has full rank.  (Z/p^2)^x is cyclic of order p(p-1): its l-coordinate (l | p(p-1)) is the
discrete log mod l, read off from x^(n/l).  (Z/8)^x = {1,3,5,7} = F_2^2 via 7 -> (1,0), 5 -> (0,1).
Counts used: primes q in [lo, hi] (every odd number <= 20,405 except 7, 21 is a count: wave 46, n <= 11)."""
import sys
def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
def factor(n):
    f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f
def primroot_mod_p2(p):
    n = p * (p - 1); fs = factor(n)
    for g in range(2, p * p):
        if g % p == 0: continue
        if all(pow(g, n // f, p * p) != 1 for f in fs): return g
def rank_mod(rows, l):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % l), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], l - 2, l) if l > 2 else 1
        rows[rk] = [(v * inv) % l for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % l:
                f = rows[i][c]; rows[i] = [(a - f * b) % l for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def certificate(P, lo, hi):
    odd = [p for p in primes_upto(P) if p > 2]
    gens = [q for q in primes_upto(hi) if q >= lo]
    # columns per l
    cols = {}            # l -> list of (p, kind)
    cols.setdefault(2, []).extend([('8', 0), ('8', 1)])
    data = {}
    for p in odd:
        n = p * (p - 1); g = primroot_mod_p2(p)
        for l in factor(n):
            h = pow(g, n // l, p * p)
            table = {pow(h, k, p * p): k for k in range(l)}
            data[(p, l)] = (n, table)
            cols.setdefault(l, []).append((p, 0))
    ok = True; report = []
    for l, cl in sorted(cols.items()):
        rows = []
        for q in gens:
            row = []
            for (p, kind) in cl:
                if p == '8':
                    r = q % 8; v = {1: (0, 0), 3: (1, 1), 5: (0, 1), 7: (1, 0)}[r]; row.append(v[kind])
                else:
                    n, table = data[(p, l)]; row.append(table[pow(q, n // l, p * p)])
            rows.append(row)
        rk = rank_mod(rows, l)
        report.append((l, len(cl), rk)); ok &= (rk == len(cl))
    return ok, report, len(gens)
if __name__ == '__main__':
    P = int(sys.argv[1]); lo = int(sys.argv[2]); hi = int(sys.argv[3])
    ok, rep, ng = certificate(P, lo, hi)
    bad = [r for r in rep if r[1] != r[2]]
    print(f"P={P}: counts = primes in [{lo},{hi}] ({ng} of them); Frattini quotient blocks (l, columns, rank): "
          f"{len(rep)} primes l; full rank for all: {ok}; deficient: {bad[:5]}")
    print("   largest blocks:", sorted(rep, key=lambda r: -r[1])[:4])
