# Residue-colouring bounds for the unit-distance graph on a multiquadratic field L = Q(sqrt a_1,...,sqrt a_k) in C.
# Lemma (PROVED in report): if P is a prime of L with c(P)=P (c = complex conjugation), every unit vector u
# (u*c(u)=1) is a P-unit; reducing mod P gives a graph homomorphism  L -> Cay(k_P, U_P)  (cosetwise),
#   U_P = {+-1}            if c in I_P  (P ramified over L+)      -> chi(L) <= 2 (p=2) or 3 (p odd)
#   U_P = mu_{q+1}, q=|k_p| if c in D_P \ I_P (P inert over L+)  -> chi(L) <= kappa(q)
import itertools, sys
from sympy import factorint, primerange, legendre_symbol

KAPPA = {2: 4, 3: 3, 4: 4, 5: 4, 7: 4, 8: 4, 9: 3, 11: 5}   # from s2_kappa.py (SAT-exact)

def sqfree(n):
    s = -1 if n < 0 else 1
    f = factorint(abs(n))
    r = 1
    for p, e in f.items():
        if e % 2: r *= p
    return s*r

def qdisc(D):
    return D if D % 4 == 1 else 4*D

def local_type(D, p):
    """behaviour of p in Q(sqrt D), D squarefree != 1: 'R','S','I'"""
    d = qdisc(D)
    if d % p == 0: return 'R'
    if p == 2:
        return 'S' if D % 8 == 1 else 'I'
    return 'S' if legendre_symbol(D % p, p) == 1 else 'I'

def analyse(avec, pmax=200, verbose=False):
    k = len(avec)
    vecs = [t for t in itertools.product((0, 1), repeat=k) if any(t)]
    Dt = {}
    for t in vecs:
        prod = 1
        for ti, a in zip(t, avec):
            if ti: prod *= a
        Dt[t] = sqfree(prod)
    assert all(D != 1 for D in Dt.values()), "a_i not independent mod squares"
    cvec = tuple(1 if a < 0 else 0 for a in avec)
    dot = lambda e, t: sum(x*y for x, y in zip(e, t)) % 2
    group = list(itertools.product((0, 1), repeat=k))
    res = []
    for p in primerange(2, pmax+1):
        T_unr = [t for t in vecs if local_type(Dt[t], p) != 'R']
        T_spl = [t for t in vecs if local_type(Dt[t], p) == 'S']
        I = [e for e in group if all(dot(e, t) == 0 for t in T_unr)]
        Dg = [e for e in group if all(dot(e, t) == 0 for t in T_spl)]
        f = len(Dg) // len(I)
        cinI = cvec in I
        cinD = cvec in Dg
        if cinI:
            bound = 2 if p == 2 else 3
            kind = 'ramified in L/L+'
        elif cinD:
            q = p**(f//2)
            bound = KAPPA.get(q, None)
            kind = f'inert in L/L+, k(L+)=F_{q}'
        else:
            continue
        res.append((p, kind, bound, f, len(I)))
    return res, Dt

def best(avec, pmax=200):
    res, _ = analyse(avec, pmax)
    bs = [r for r in res if r[2] is not None]
    return (min(bs, key=lambda r: r[2]) if bs else None), res

if __name__ == "__main__":
    fields = {
        'Q(sqrt-1)=Q^2-plane': [-1],
        'Q(sqrt-2)': [-2],
        'Q(sqrt-3) Eisenstein': [-3],
        'Q(sqrt-7)': [-7],
        'Q(zeta8)=Q(sqrt2)^2': [-1, -2],
        'Q(zeta12)=Q(sqrt3)^2': [-1, -3],
        'Moser field Q(sqrt-3,sqrt-11) [Z[w1,w3]]': [-3, -11],
        'Polymath Q(sqrt-3,sqrt-11,sqrt-15) [Z[w1,w3,w4]]': [-3, -11, -15],
        'Heegner rung Q(sqrt-3,sqrt-7)': [-3, -7],
        'Heegner rung Q(sqrt-3,sqrt-19)': [-3, -19],
        'Heegner rung Q(sqrt-3,sqrt-43)': [-3, -43],
        'Heegner rung Q(sqrt-3,sqrt-67)': [-3, -67],
        'Heegner rung Q(sqrt-3,sqrt-163)': [-3, -163],
        'Q(sqrt-3,sqrt-11,sqrt-19)': [-3, -11, -19],
        'Heegner compositum d=3,11,19,43,67,163': [-3, -11, -19, -43, -67, -163],
        'Heegner compositum + sqrt-7': [-3, -7, -11, -19, -43, -67, -163],
        'Polymath + sqrt-7': [-3, -11, -15, -7],
        'Polymath + sqrt-19': [-3, -11, -15, -19],
        'Polymath + sqrt-23 (w6)': [-3, -11, -15, -23],
        'Q(sqrt-3,sqrt-11,sqrt-7)': [-3, -11, -7],
        'Q(sqrt-3,sqrt-15)': [-3, -15],
        'Q(sqrt-11,sqrt-15)': [-11, -15],
    }
    for name, av in fields.items():
        b, res = best(av)
        unram = all(r[1] != 'ramified in L/L+' for r in res)
        lst = ', '.join(f"p={r[0]}:{'<=' + str(r[2]) if r[2] else '?(q=' + r[1].split('F_')[1] + ')'}" for r in res[:8])
        print(f"{name:55s} best={'chi<=' + str(b[2]) + ' via p=' + str(b[0]) + ' (' + b[1] + ')' if b else 'none<=200'} | L/L+ unram(fin)<=200: {unram} | {lst}")
