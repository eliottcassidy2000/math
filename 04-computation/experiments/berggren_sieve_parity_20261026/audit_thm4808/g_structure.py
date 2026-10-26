#!/usr/bin/env python3
"""Independent audit of THM-4808 (G) -- written from the theorem text only, no lane code imported.

On X = P^1(Z/N), N odd (a field when N is prime), with branches given by their MATRICES
    U = [2 -1; 1 0],  A = [2 1; 1 0],  D = [1 2; 0 1],  tau = [-1 2; 0 1]  (t -> 2 - t),
    iota = [-1 0; 0 1] (t -> -t),  S = [0 -1; 1 0] (t -> -1/t),  V = {1, iota, S, S iota},
we check, exactly and for every dart e in X:
  (T1) U e = D(S e), A e = D(S(iota e)), tau e = D(iota e)            (theta-group factorisation)
  (T2) M(e) = [Ue, Ae, De, tau e] equals the multiset D(V.e) with every element counted 4/|V.e| times;
       no V-orbit has size 1; census of the small orbits and their stabilisers
  (T3) tau^2 = 1, tau e in D(V.e), tail(tau e) = head(e), head(tau e) = tail(e)
       (head(e) = V.e, tail(e) = V.D^-1 e)
  (T4) B = P_U + P_A + P_D satisfies B[e,f] = (4/|head e|) [head e = tail f] - [f = tau e]  (>= 0),
       i.e. B is the Hashimoto operator of the edge-indexed graph whose dart f is repeated 4/|tail f| times
  (T5) Ihara-Bass, EXACT polynomial identity (python-flint, integer arithmetic):
          charpoly(B)(x) * (x - 1) = (x - 1)^alpha (x + 1)^beta * x^nV F_A(x + 3/x),
       F_A = charpoly of the weighted adjacency A[O,O'] = #{v in V : D v s_O in O'} (row sums 4)
  (T6) the same for the branch-signed B_chi = P_U - P_A + P_D, with the twisted adjacency A_chi on
       chi-equivariant functions, chi(1) = chi(S) = 1, chi(iota) = chi(S iota) = -1 (chi(A) = -1):
          charpoly(B_chi)(x) = (x - 1)^alpha' (x + 1)^beta' * x^n' F_{A_chi}(x + 3/x)
  (T7) the C-involution C(t) = (t+1)/(t-1) (matrix [1 1; 1 -1]): C V C = V, C D C = U, C A C = A, C tau C = tau,
       so C commutes with B and B_chi (an automorphism of G_N swapping the U and D labels).
Prints alpha, beta, alpha', beta' and the vertex census, for comparison with berggren_ihara.out.
"""
import sys
from math import gcd
from collections import Counter
import flint


class PrimeCanon(dict):
    """canonical form on P^1(F_q): (t, 1) or (1, 0)."""
    def __init__(self, q):
        super().__init__(); self.q = q

    def __missing__(self, key):
        x, y = key; q = self.q
        return ((x * pow(y, -1, q)) % q, 1) if y % q else (1, 0)


def make_space(N):
    """canonical representatives of P^1(Z/N): primitive vectors mod N up to units."""
    if is_prime(N):
        canon = PrimeCanon(N)
        pts = [(t, 1) for t in range(N)] + [(1, 0)]
        return pts, {p: i for i, p in enumerate(pts)}, canon
    units = [u for u in range(1, N) if gcd(u, N) == 1]
    canon = {}
    pts = []
    for x in range(N):
        for y in range(N):
            if gcd(gcd(x, y), N) != 1:
                continue
            key = min(((u * x) % N, (u * y) % N) for u in units)
            canon[(x, y)] = key
    for key in sorted(set(canon.values())):
        pts.append(key)
    idx = {p: i for i, p in enumerate(pts)}
    return pts, idx, canon


def matact(N, canon, M, p):
    a, b, c, d = M
    x, y = p
    return canon[((a * x + b * y) % N, (c * x + d * y) % N)]


U_ = (2, -1, 1, 0); A_ = (2, 1, 1, 0); D_ = (1, 2, 0, 1); Dinv = (1, -2, 0, 1)
TAU = (-1, 2, 0, 1); IOTA = (-1, 0, 0, 1); S_ = (0, -1, 1, 0); SIOTA = (0, 1, 1, 0)  # S o iota : t -> 1/t
C_ = (1, 1, 1, -1)
VLIST = [(1, 0, 0, 1), IOTA, S_, SIOTA]
CHI = [1, -1, 1, -1]


def x_plus_3_over_x(F, n):
    """x^n F(x + 3/x) as an fmpz_poly, F of degree n."""
    X = flint.fmpz_poly([0, 1]); Q = flint.fmpz_poly([3, 0, 1])
    coeffs = F.coeffs()
    G = flint.fmpz_poly(0)
    for k, ck in enumerate(coeffs):
        G += int(ck) * Q ** k * X ** (n - k)
    return G


def pm1_split(P):
    """if P = c (x-1)^a (x+1)^b return (a, b) else None."""
    a = b = 0
    xm1 = flint.fmpz_poly([-1, 1]); xp1 = flint.fmpz_poly([1, 1])
    while P.degree() > 0:
        q, r = divmod(P, xm1)
        if r == 0:
            P = q; a += 1; continue
        q, r = divmod(P, xp1)
        if r == 0:
            P = q; b += 1; continue
        return None
    return (a, b)


def audit(N):
    pts, idx, canon = make_space(N)
    n = len(pts)
    act = lambda M, p: matact(N, canon, M, p)
    ok = {k: True for k in ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7')}
    # orbits
    orb_of = {}
    orbits = []
    for p in pts:
        if p in orb_of:
            continue
        O = frozenset(act(v, p) for v in VLIST)
        for r in O:
            orb_of[r] = len(orbits)
        orbits.append(O)
    sizes = Counter(len(O) for O in orbits)
    if 1 in sizes:
        ok['T2'] = False
    small = []
    for O in orbits:
        if len(O) == 2:
            s = min(O)
            stab = [i for i, v in enumerate(VLIST) if act(v, s) == s]
            small.append((tuple(sorted(O)), tuple(['1', 'iota', 'S', 'S.iota'][i] for i in stab)))
    head = lambda e: orb_of[e]
    tail = lambda e: orb_of[act(Dinv, e)]
    B = [[0] * n for _ in range(n)]
    Bc = [[0] * n for _ in range(n)]
    for e in pts:
        Ue, Ae, De, Te = act(U_, e), act(A_, e), act(D_, e), act(TAU, e)
        # T1
        if Ue != act(D_, act(S_, e)) or Ae != act(D_, act(S_, act(IOTA, e))) or Te != act(D_, act(IOTA, e)):
            ok['T1'] = False
        # T2
        O = orbits[orb_of[e]]
        want = Counter({act(D_, p): 4 // len(O) for p in O})
        if Counter([Ue, Ae, De, Te]) != want:
            ok['T2'] = False
        # T3
        if act(TAU, Te) != e or Te not in want or tail(Te) != head(e) or head(Te) != tail(e):
            ok['T3'] = False
        i = idx[e]
        for f, s in ((Ue, 1), (Ae, 1), (De, 1)):
            B[i][idx[f]] += 1
        for f, s in ((Ue, 1), (Ae, -1), (De, 1)):
            Bc[i][idx[f]] += s
        # T7
        if act(C_, act(D_, act(C_, e))) != act(U_, e) or act(C_, act(A_, act(C_, e))) != act(A_, e) \
                or act(C_, act(TAU, act(C_, e))) != act(TAU, e):
            ok['T7'] = False
        if any(orb_of[act(C_, act(v, act(C_, e)))] != orb_of[e] for v in VLIST):
            ok['T7'] = False
    # T4
    for e in pts:
        for f in pts:
            w = (4 // len(orbits[head(e)])) * (head(e) == tail(f)) - (f == act(TAU, e))
            if w < 0 or B[idx[e]][idx[f]] != w:
                ok['T4'] = False
    # T5: weighted adjacency on orbits
    nV = len(orbits)
    rep = [min(O) for O in orbits]
    Aw = [[0] * nV for _ in range(nV)]
    for a, s in enumerate(rep):
        for v in VLIST:
            Aw[a][orb_of[act(D_, act(v, s))]] += 1
    FA = flint.fmpz_mat(Aw).charpoly()
    PB = flint.fmpz_mat(B).charpoly()
    G = x_plus_3_over_x(FA, nV)
    q, r = divmod(PB * flint.fmpz_poly([-1, 1]), G)
    split = pm1_split(q) if r == 0 else None
    if split is None:
        ok['T5'] = False
    # T6: twisted adjacency on chi-equivariant functions
    alive = []
    for a, s in enumerate(rep):
        stab = [i for i, v in enumerate(VLIST) if act(v, s) == s]
        if all(CHI[i] == 1 for i in stab):
            alive.append(a)
    pos = {a: k for k, a in enumerate(alive)}
    nC = len(alive)
    Ac = [[0] * nC for _ in range(nC)]
    for a in alive:
        s = rep[a]
        for vi, v in enumerate(VLIST):
            p = act(D_, act(v, s))
            b = orb_of[p]
            if b not in pos:
                continue
            # p = u s_b, find u
            us = [ui for ui, u in enumerate(VLIST) if act(u, rep[b]) == p]
            vals = set(CHI[ui] for ui in us)
            assert len(vals) == 1
            Ac[pos[a]][pos[b]] += CHI[vi] * vals.pop()
    FC = flint.fmpz_mat(Ac).charpoly() if nC else flint.fmpz_poly([1])
    PC = flint.fmpz_mat(Bc).charpoly()
    GC = x_plus_3_over_x(FC, nC)
    qc, rc = divmod(PC, GC)
    splitc = pm1_split(qc) if rc == 0 else None
    if splitc is None:
        ok['T6'] = False
    return dict(N=N, npts=n, nV=nV, sizes=dict(sizes), small=small, ok=ok, split=split, splitc=splitc,
                nC=nC, FA=FA, FC=FC, Aw=Aw, Ac=Ac)


def is_prime(N):
    return N > 1 and all(N % d for d in range(2, int(N ** 0.5) + 1))


if __name__ == '__main__':
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 212
    mods = [q for q in range(3, QMAX) if is_prime(q)] + [9, 15, 21, 25, 27, 33, 35, 45, 49, 63, 75, 105]
    allok = True
    pred_bad = 0
    for N in mods:
        r = audit(N)
        allok &= all(r['ok'].values())
        a, b = r['split'] if r['split'] else (None, None)
        ac, bc = r['splitc'] if r['splitc'] else (None, None)
        line = (f"N={N:4d} |P1|={r['npts']:4d} vertices={r['nV']:4d} orbit sizes={r['sizes']} "
                f"checks={''.join('1' if r['ok'][k] else '0' for k in sorted(r['ok']))} "
                f"B: mult(x-1)={a} mult(x+1)={b} | B_chi: mult(x-1)={ac} mult(x+1)={bc} twisted dim={r['nC']}")
        if is_prime(N) and N > 3:
            # predicted multiplicities (derived in REPORT.md): B: q=3 mod 4: (q-3)/4, (q+1)/4; q=1 mod 4: (q-5)/4, (q-1)/4
            # B_chi: q=3 mod 4: (q+9)/4, (q+1)/4; q=1 mod 4: (q+7)/4, (q-1)/4
            if N % 4 == 3:
                pred = ((N - 3) // 4, (N + 1) // 4, (N + 9) // 4, (N + 1) // 4)
            else:
                pred = ((N - 5) // 4, (N - 1) // 4, (N + 7) // 4, (N - 1) // 4)
            match = pred == (a, b, ac, bc)
            pred_bad += not match
            line += f" predicted={pred} match={match}"
        if N <= 13 or not all(r['ok'].values()) or not is_prime(N):
            line += f" small orbits={r['small']}"
        print(line)
    print(f"ALL STRUCTURAL CHECKS T1-T7 PASS on every modulus: {allok}; multiplicity-formula mismatches: {pred_bad}")
