"""Audit B (2026-10-12) of THM-4624 / HYP-9256: independent recomputation of |Aut G_q|,
G_q = Cay(F_{q^2}, mu_{q+1}).  Independent of gq_automorphisms.py: own field arithmetic
(F_{q^2} = F_p[x]/(m) with m a RANDOM irreducible monic polynomial of degree 2e, seed 20261012;
multiplication through exp/log tables of a primitive element g), mu_{q+1} = {g^{(q-1) i}},
own dreadnaut driver (sparse nauty), explicit check that the generators of Aff(q)
(translation by 1, multiplication by g^{q-1}, Frobenius z -> z^p) are graph automorphisms,
stabiliser-of-0 orbit count compared with 1 + #(orbits of j -> p j on Z/(q-1)).
Usage: python3 audit_B_gq_census.py q1 q2 ...
"""
import sys, random, subprocess, re, time
random.seed(20261012)

def factor_pp(q):
    for p in range(2, q + 1):
        if q % p == 0:
            e = 0
            while q % p == 0: q //= p; e += 1
            return (p, e) if q == 1 else None

def prime_factors(n):
    out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: out.append(n)
    return out

# polynomials over F_p as lists of coefficients, low degree first
def pmod(f, m, p):
    f = f[:]
    d = len(m) - 1
    while len(f) > d:
        c = f[-1]
        if c:
            sh = len(f) - 1 - d
            for i in range(d + 1): f[sh + i] = (f[sh + i] - c * m[i]) % p
        f.pop()
    return f + [0] * (d - len(f))

def pmul(f, g, m, p):
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g): r[i + j] = (r[i + j] + a * b) % p
    return pmod(r, m, p)

def ppow(f, n, m, p):
    d = len(m) - 1
    r = [1] + [0] * (d - 1)
    while n:
        if n & 1: r = pmul(r, f, m, p)
        f = pmul(f, f, m, p); n >>= 1
    return r

def pgcd(a, b, p):
    def trim(x):
        while x and x[-1] == 0: x.pop()
        return x
    a, b = trim(a[:]), trim(b[:])
    while b:
        # a mod b
        while len(a) >= len(b) and a:
            c = a[-1] * pow(b[-1], -1, p) % p
            sh = len(a) - len(b)
            for i in range(len(b)): a[sh + i] = (a[sh + i] - c * b[i]) % p
            a = trim(a)
        a, b = b, a
    return a

def irreducible(m, p):
    """Rabin test for monic m of degree n."""
    n = len(m) - 1
    x = [0, 1] + [0] * (n - 2) if n >= 2 else [0, 1]
    x = pmod(x, m, p)
    xp = ppow(x, p ** n, m, p)
    if xp != x: return False
    for r in prime_factors(n):
        y = ppow(x, p ** (n // r), m, p)
        diff = [(a - b) % p for a, b in zip(y, x)]
        g = pgcd(diff, m[:], p)
        if len(g) > 1: return False
    return True

def build_field(p, n):
    while True:
        m = [random.randrange(p) for _ in range(n)] + [1]
        if m[0] and irreducible(m, p): break
    N = p ** n
    def enc(f):
        v = 0
        for c in reversed(f): v = v * p + c
        return v
    def dec(v):
        f = []
        for _ in range(n): f.append(v % p); v //= p
        return f
    # primitive element
    order = N - 1
    pf = prime_factors(order)
    one = [1] + [0] * (n - 1)
    while True:
        g = [random.randrange(p) for _ in range(n)]
        if not any(g): continue
        if all(ppow(g, order // r, m, p) != one for r in pf): break
    exp = [0] * (N - 1); cur = one
    for i in range(N - 1):
        exp[i] = enc(cur); cur = pmul(cur, g, m, p)
    assert enc(cur) == enc(one)
    log = {v: i for i, v in enumerate(exp)}
    assert len(log) == N - 1
    return m, enc, dec, exp, log

def add(u, v, p, n):
    if p == 2: return u ^ v
    r = 0; mult = 1
    for _ in range(n):
        r += ((u % p + v % p) % p) * mult
        u //= p; v //= p; mult *= p
    return r

def run_dreadnaut(adj, fix=None):
    n = len(adj)
    lines = ["n=%d $=0 As g" % n]
    for i in range(n):
        lines.append(" ".join(str(j) for j in adj[i] if j > i) + (";" if i < n - 1 else "."))
    if fix is not None: lines.append("f=[%d]" % fix)
    lines.append("x"); lines.append("q")
    out = subprocess.run(["dreadnaut"], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout
    g = re.findall(r"grpsize=([0-9.e+]+)", out)
    o = re.findall(r"(\d+) orbits?;", out)
    return g[-1] if g else None, int(o[-1]) if o else None, out

def main(qs):
    print("q   p^e  |V|   |Aut G_q|   q^2(q+1)2e  ratio  orbits  |Stab_0|  Stab orbits  1+FrobOrb  gens-ok  time(s)")
    for q in qs:
        t0 = time.time()
        p, e = factor_pp(q); n = 2 * e
        m, enc, dec, exp, log = build_field(p, n)
        N = p ** n
        a = exp[q - 1]                       # generator of mu_{q+1}
        mu = [exp[(q - 1) * i % (N - 1)] for i in range(q + 1)]
        assert len(set(mu)) == q + 1
        # sanity: every z in mu has z^{q+1} = 1 and nothing else does (norm = z^{q+1})
        assert all(((q - 1) * i * (q + 1)) % (N - 1) == 0 for i in range(q + 1))
        adj = [sorted(add(z, w, p, n) for w in mu) for z in range(N)]
        adjs = [set(x) for x in adj]
        # Aff generators: translation by 1, multiplication by a, Frobenius
        def mul(z, w):
            if z == 0 or w == 0: return 0
            return exp[(log[z] + log[w]) % (N - 1)]
        def frob(z):
            return 0 if z == 0 else exp[(log[z] * p) % (N - 1)]
        gens_ok = True
        for name, f in (("transl", lambda z: add(z, 1, p, n)), ("mul", lambda z: mul(z, a)), ("frob", frob)):
            img = [f(z) for z in range(N)]
            if len(set(img)) != N: gens_ok = False
            for z in range(N):
                fz = img[z]
                if any(img[w] not in adjs[fz] for w in adj[z]): gens_ok = False; break
        grp, orb, out = run_dreadnaut(adj)
        sg, sorb, _ = run_dreadnaut(adj, fix=0)
        aff = q * q * (q + 1) * 2 * e
        # Frobenius orbits on F_q^x = orbits of j -> p j on Z/(q-1)
        seen = set(); forb = 0
        for j in range(q - 1):
            if j in seen: continue
            forb += 1; k = j
            while k not in seen: seen.add(k); k = (k * p) % (q - 1)
        G = int(float(grp))
        print("%-3d %d^%-2d %-5d %-11d %-11d %-6s %-7s %-9s %-12s %-10d %-8s %.1f" % (
            q, p, e, N, G, aff, ("%g" % (G / aff)), orb, sg, sorb, 1 + forb, gens_ok, time.time() - t0))
        sys.stdout.flush()

if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]])
