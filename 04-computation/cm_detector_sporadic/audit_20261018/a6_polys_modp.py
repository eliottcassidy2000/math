"""Audit check 6 (HYP-9330 (B), (C)) with hand-written polynomial arithmetic over F_p.

F_p(x) = sum_{k<p} u_k x^k,  V_p(x) = sum_{j<p} v_j x^j,  v_j = (u_{p+j} - u_1 u_j)/p   (all mod p).
Weight 1: is F_p squarefree?  gcd(F_p, V_p) = 1?
Weight 2: does rad(F_p) divide V_p?  simple-root part S of F_p divides the singular-fibre polynomial?
          multiplicities of the other irreducible factors (is F_p/S a square?)  which primes have S != 1?
Also records the Lucas congruence u_{p+j} == u_1 u_j (mod p), j < p.

usage: python a6_polys_modp.py PLO PHI [families...]
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, primes_in, FAMILIES

# singular-fibre polynomials in x = 1/m (coefficients low degree first)
SING = {
    'Apery': [1, -34, 1],          # x^2 - 34x + 1
    'Domb': [1, -20, 64],          # (1-4x)(1-16x)
    's7': [1, -26, -27],           # (1+x)(1-27x)
    's10': [1, -12, -64],          # (1+4x)(1-16x)
    'A274786': [1, -44, -16],      # 1 - 44x - 16x^2
    'A': [1, -7, -8],              # 1 - 7x - 8x^2 = (1+x)(1-8x)
    'B': [1, -9, 27], 'C': [1, -10, 9], 'D': [1, -11, -1], 'E': [1, -12, 32], 'F': [1, -17, 72],
}


# ---------- polynomial arithmetic over F_p, lists low-degree-first, normalised (no trailing zeros) ----------
def norm(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def pmod(a, b, p):
    a = a[:]
    norm(a)
    b = norm(b[:])
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        c = a[-1] * inv % p
        s = len(a) - len(b)
        for i, bi in enumerate(b):
            a[s + i] = (a[s + i] - c * bi) % p
        norm(a)
    return a


def pdivexact(a, b, p):
    a = norm([x % p for x in a])
    b = norm([x % p for x in b])
    inv = pow(b[-1], -1, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        s = len(a) - len(b)
        q[s] = c
        for i, bi in enumerate(b):
            a[s + i] = (a[s + i] - c * bi) % p
        norm(a)
    assert not a, 'not exact'
    return norm(q)


def pgcd(a, b, p):
    a = norm([x % p for x in a])
    b = norm([x % p for x in b])
    while b:
        a, b = b, pmod(a, b, p)
    if not a:
        return a
    inv = pow(a[-1], -1, p)
    return [x * inv % p for x in a]


def pderiv(a, p):
    return norm([(i * a[i]) % p for i in range(1, len(a))])


def deg(a):
    return len(norm(a[:])) - 1


def radical(F, p):
    # valid when every multiplicity is < p (true here: deg F < p)
    g = pgcd(F, pderiv(F, p), p)
    return pdivexact(F, g, p) if deg(g) > 0 else norm([x % p for x in F])


def sqfree_decomposition(F, p):
    """Yun-type: returns dict e -> product of irreducible factors of multiplicity exactly e (requires all e < p)."""
    out = {}
    e = 1
    G = norm([x % p for x in F])
    R = radical(G, p)
    while deg(G) > 0:
        # factors of multiplicity >= e are exactly rad(G); after dividing G by rad(G) multiplicities drop by 1
        Gn = pdivexact(G, R, p)
        Rn = radical(Gn, p) if deg(Gn) > 0 else [1]
        Ee = pdivexact(R, Rn, p) if deg(Rn) > 0 else R
        if deg(Ee) > 0:
            out[e] = Ee
        G, R, e = Gn, Rn, e + 1
    return out


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def kron(D, p):  # Kronecker symbol (D/p) for odd prime p
    return legendre(D, p)


def analyse(name, primes):
    w = FAMILIES[name][1]
    u = seq(name, 2 * max(primes))
    rows = []
    for p in primes:
        lucas = all((u[p + j] - u[1] * u[j]) % p == 0 for j in range(p))
        F = norm([u[k] % p for k in range(p)])
        V = norm([((u[p + j] - u[1] * u[j]) // p) % p for j in range(p)]) if lucas else None
        row = dict(p=p, lucas=lucas, degF=deg(F))
        if not lucas:
            rows.append(row)
            continue
        g = pgcd(F, pderiv(F, p), p)
        row['squarefree'] = deg(g) == 0
        row['gcdFV_deg'] = deg(pgcd(F, V, p)) if V else deg(F)
        R = radical(F, p)
        row['radF_divides_V'] = (not V) or (not pmod(V, R, p))
        dec = sqfree_decomposition(F, p)
        row['mult'] = {e: deg(f) for e, f in dec.items()}       # multiplicity -> total degree of such factors
        S = dec.get(1, [1])
        row['simple_deg'] = deg(S)
        sing = SING.get(name)
        if sing is not None and deg(S) > 0:
            row['simple_divides_sing'] = not pmod([c % p for c in sing], S, p)
        else:
            row['simple_divides_sing'] = True
        # odd multiplicities other than the simple part?
        row['odd_mult_ge3'] = [e for e in dec if e % 2 == 1 and e >= 3]
        rows.append(row)
    return rows


if __name__ == '__main__':
    plo, phi = int(sys.argv[1]), int(sys.argv[2])
    fams = sys.argv[3:] or ['A', 'B', 'C', 'D', 'E', 'F', 'Apery', 'Domb', 's7', 's10', 'A274786']
    primes = primes_in(plo, phi)
    summary = {}
    for name in fams:
        t0 = time.time()
        rows = analyse(name, primes)
        w = FAMILIES[name][1]
        nl = [r['p'] for r in rows if not r['lucas']]
        if w == 1:
            bad = [r['p'] for r in rows if r['lucas'] and not (r['squarefree'] and r['gcdFV_deg'] == 0)]
            print(f'{name:8s} w=1, {len(primes)} primes in [{plo},{phi}]: Lucas fails at {nl}; '
                  f'(B) [F_p squarefree and gcd(F_p,V_p)=1] fails at {bad}   [{time.time() - t0:.1f}s]', flush=True)
        else:
            bad_rad = [r['p'] for r in rows if r['lucas'] and not r['radF_divides_V']]
            bad_sing = [r['p'] for r in rows if r['lucas'] and not r['simple_divides_sing']]
            simple_primes = [r['p'] for r in rows if r['lucas'] and r['simple_deg'] > 0]
            odd3 = [(r['p'], r['odd_mult_ge3']) for r in rows if r['lucas'] and r['odd_mult_ge3']]
            # which discriminant's inert primes match the simple-root primes?
            cand = []
            for D in [-3, -4, -7, -8, -11, -12, -15, -16, -19, -20, -24, -27, -28, -35, -40, -51, -52, -60, -84, -88, -120]:
                inert = [p for p in primes if kron(D, p) == -1]
                if inert == simple_primes:
                    cand.append(D)
            mults = sorted({e for r in rows if r['lucas'] for e in r['mult']})
            print(f'{name:8s} w=2, {len(primes)} primes in [{plo},{phi}]: Lucas fails at {nl}; rad(F_p)|V_p fails at {bad_rad}; '
                  f'simple part divides sing-poly fails at {bad_sing}; simple roots at {len(simple_primes)} primes '
                  f'= inert primes of D in {cand}; odd multiplicities >=3 at {odd3}; multiplicities seen {mults}   '
                  f'[{time.time() - t0:.1f}s]', flush=True)
            print(f'          simple-root primes: {simple_primes}', flush=True)
        summary[name] = rows
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'a6_polys_modp_{plo}_{phi}.json'), 'w') as f:
        json.dump(summary, f, default=str)
