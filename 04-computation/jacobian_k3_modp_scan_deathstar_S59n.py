#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- k=3 emptiness: mod-p projective scan.
For each prime p and each delta in F_p^*: build the c1-linear system mod p,
compute its kernel basis, then scan PROJECTIVE lambda-space for common zeros
of the c0-quadrics.  Nonempty C-variety of dim >= 0 => F_p-points for a
positive density of p (splitting) or all large p (Lang-Weil, dim >= 1).
Zero points across several primes => strong emptiness evidence.
VALIDATION: the same scanner at k=2 MUST find the witness point mod every p.
"""
import itertools, sys

def build_eqs(k, degB, degC, degE0, delta, p):
    """returns (basis over F_p of c1-kernel, c0-quadrics as dicts over basis idx)."""
    NB, NC, NE = degB+1, degC+1, degE0+1
    NU = NB + NC + NE
    # polynomial helpers over F_p, univariate lists
    def pmul(a, b):
        out = [0]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y: out[i+j] = (out[i+j] + x*y) % p
        return out
    def padd(*ps):
        n = max(len(q) for q in ps)
        out = [0]*n
        for q in ps:
            for i, x in enumerate(q): out[i] = (out[i] + x) % p
        return out
    def pscale(q, s): return [(x*s) % p for x in q]
    def pd(q): return [(q[i]*i) % p for i in range(1, len(q))]
    def shift(q, m): return [0]*m + list(q)
    A = [1]
    for _ in range(k+1): A = pmul(A, [1, 1])
    D = pscale(pmul([1,1],[1,1]), delta)
    Ap, Dp = pd(A), pd(D)
    # c1 rows: linear in unknowns (B coeffs, C coeffs, E0 coeffs)
    # build by unit vectors
    def c1_of(B, C, E0):
        Bp, Cp, E0p = pd(B), pd(C), pd(E0)
        M1 = padd(pmul(Ap, D), pscale(pmul(A, Dp), p-1))
        return padd(
            pmul(E0, M1), pscale(pmul(E0p, pmul(A, D)), k-1),
            pscale(shift(pmul(Bp, D), k), p-2), pscale(shift(pmul(B, D), k-1), (p - 2*k) % p),
            pscale(shift(pmul(B, Dp), k), k),
            pscale(pmul(A, C), k+1), pscale(shift(pmul(A, Cp), 1), k+1),
            pscale(shift(pmul(Ap, C), 1), p-1))
    rows = {}
    for j in range(NU):
        B = [0]*NB; C = [0]*NC; E0 = [0]*NE
        if j < NB: B[j] = 1
        elif j < NB+NC: C[j-NB] = 1
        else: E0[j-NB-NC] = 1
        c1 = c1_of(B, C, E0)
        for td, c in enumerate(c1):
            if c:
                rows.setdefault(td, [0]*NU)
                rows[td][j] = c
    M = [rows[td] for td in sorted(rows)]
    # kernel mod p
    m = len(M); r = 0; piv = []
    M = [row[:] for row in M]
    for c in range(NU):
        pv = next((i for i in range(r, m) if M[i][c] % p), None)
        if pv is None: continue
        M[r], M[pv] = M[pv], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(x*inv) % p for x in M[r]]
        for i in range(m):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f*b) % p for a, b in zip(M[i], M[r])]
        r += 1; piv.append(c)
    free = [c for c in range(NU) if c not in piv]
    basis = []
    for fc in free:
        v = [0]*NU; v[fc] = 1
        for i, c in enumerate(piv): v[c] = (-M[i][fc]) % p
        basis.append(v)
    if not basis: return [], []
    NL = len(basis)
    # c0 quadrics in lambda: evaluate c0 for lambda unit/pair combos to extract
    def funcs_of_lam(lam):
        vec = [0]*NU
        for bi, l in enumerate(lam):
            if l:
                for j in range(NU): vec[j] = (vec[j] + l*basis[bi][j]) % p
        B = vec[:NB]; C = vec[NB:NB+NC]; E0 = vec[NB+NC:]
        return B, C, E0
    def c0_of(B, C, E0):
        Bp, Cp, E0p = pd(B), pd(C), pd(E0)
        K0 = padd(E0, shift(E0p, 1))
        M0 = padd(shift(pmul(Bp, D), k), pscale(shift(pmul(B, D), k-1), k),
                  pscale(pmul(A, C), p-1), pscale(shift(pmul(A, Cp), 1), p-1))
        return padd(
            pmul(K0, M0),
            pscale(pmul(E0p, padd(shift(pmul(Bp, D), k+1),
                                  pscale(shift(pmul(A, Cp), 2), p-1))), p-1),
            pscale(shift(padd(pmul(Bp, C), pscale(pmul(B, Cp), p-k)), k+1), p-1))
    # extract quadric coefficients: Q(l) = sum q_ij l_i l_j via polarization
    unit = lambda i: [1 if j == i else 0 for j in range(NL)]
    c0_cache = {}
    def c0_lam(lam):
        key = tuple(lam)
        if key not in c0_cache:
            c0_cache[key] = c0_of(*funcs_of_lam(lam))
        return c0_cache[key]
    qii = [c0_lam(unit(i)) for i in range(NL)]
    qij = {}
    for i in range(NL):
        for j in range(i+1, NL):
            lam = unit(i); lam[j] = 1
            full = c0_lam(lam)
            n = max(len(full), len(qii[i]), len(qii[j]))
            ext = lambda q: q + [0]*(n-len(q))
            qij[(i,j)] = [(a - b - c) % p for a, b, c in zip(ext(full), ext(qii[i]), ext(qii[j]))]
    def eval_c0(lam):
        n = 12
        out = [0]*n
        for i in range(NL):
            if lam[i]:
                li2 = (lam[i]*lam[i]) % p
                for td, c in enumerate(qii[i]):
                    if c: out[td] = (out[td] + li2*c) % p
        for (i, j), q in qij.items():
            if lam[i] and lam[j]:
                lij = (lam[i]*lam[j]) % p
                for td, c in enumerate(q):
                    if c: out[td] = (out[td] + lij*c) % p
        return out
    return basis, eval_c0

def scan(k, degB, degC, degE0, p, need_nondeg=True):
    total_pts = 0
    for delta in range(1, p):
        basis, eval_c0 = build_eqs(k, degB, degC, degE0, delta, p)
        if not basis: continue
        NL = len(basis)
        # projective scan: first nonzero coord = 1
        for lead in range(NL):
            for rest in itertools.product(range(p), repeat=NL-lead-1):
                lam = [0]*lead + [1] + list(rest)
                c0 = eval_c0(lam)
                # constancy: all t^>=1 coeffs zero; t^0 (const) must be NONZERO (Keller)
                if any(c0[1:]):
                    continue
                if need_nondeg and c0[0] == 0:
                    continue
                total_pts += 1
                if total_pts <= 3:
                    print(f"    p={p} delta={delta}: POINT lam={lam} c0={c0[0]}")
    return total_pts

print("=== VALIDATION: k=2 scanner must find the witness mod p ===")
for p in [5, 7, 11]:
    n = scan(2, 2, 2, 1, p)
    print(f"  k=2, p={p}: {n} projective Keller points (expect > 0)")

print("\n=== k=3 scan ===")
for p in [5, 7, 11]:
    n = scan(3, 4, 4, 3, p)
    print(f"  k=3, p={p}: {n} projective Keller points")
