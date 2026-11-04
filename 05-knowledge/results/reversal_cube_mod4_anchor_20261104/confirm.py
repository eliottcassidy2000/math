"""Confirmatory runs P1-P5 of PREREGISTRATION.md (kind-pasteur-2026-11-04-S9).

Exact integer / F2 computations only.  Run:  g++ -O2 -std=c++17 fcoef.cpp -o fcoef.exe ; python confirm.py
Writes confirm.json; prints a summary.  Every check raises on failure only via explicit counters (no asserts
are relied on; python -O gives identical output).
"""
import itertools, json, os, random, subprocess, sys
from wmodel import gamma, contrast, Lfun

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'fcoef.exe' if os.name == 'nt' else 'fcoef')

# ---------------------------------------------------------------- tournaments
def random_tournament(m, rng):
    Q = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < 0.5: Q[i][j] = 1
            else: Q[j][i] = 1
    return Q

def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        Q = [[0] * n for _ in range(n)]
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: Q[i][j] = 1
            else: Q[j][i] = 1
        yield Q

def conv(M): return [list(r) for r in zip(*M)]
def TT(k): return [[1 if p < q else 0 for q in range(k)] for p in range(k)]

def subst(Q, children):
    blocks, n = [], 0
    for C in children:
        blocks.append(list(range(n, n + len(C)))); n += len(C)
    A = [[0] * n for _ in range(n)]
    for i, Bi in enumerate(blocks):
        for p, u in enumerate(Bi):
            for q, v in enumerate(Bi): A[u][v] = children[i][p][q]
        for j, Bj in enumerate(blocks):
            if i != j and Q[i][j]:
                for u in Bi:
                    for v in Bj: A[u][v] = 1
    return A, blocks

def is_module(A, S):
    S = set(S)
    return all(len({A[v][u] for u in S}) == 1 for v in range(len(A)) if v not in S)

def is_prime(Q):
    m = len(Q)
    return not any(is_module(Q, S) for r in range(2, m) for S in itertools.combinations(range(m), r))

_EP = {}
def endpoint_parity(Q):
    key = tuple(map(tuple, Q))
    if key in _EP: return _EP[key]
    _EP[key] = _endpoint_parity(Q); return _EP[key]

def _endpoint_parity(Q):
    n = len(Q); N = 1 << n; S = [0] * n; E = [0] * n
    for u in range(n):
        dp = [[0] * n for _ in range(N)]; dp[1 << u][u] = 1
        for m in range(N):
            for v in range(n):
                x = dp[m][v]
                if not x: continue
                for w in range(n):
                    if not m >> w & 1 and Q[v][w]: dp[m | 1 << w][w] += x
        for v in range(n):
            S[u] += dp[N - 1][v]; E[v] += dp[N - 1][v]
    return [(S[i] + E[i]) % 2 for i in range(n)]

# ---------------------------------------------------------------- kernel model
def kernel(Q, j, orient, wbeta, wo):
    """Weighted tournament Q[pair at j]; u=j, beta=j+1.  orient 0: u->beta, 1: beta->u."""
    m = len(Q); pair = [[0, 1], [0, 0]] if orient == 0 else [[0, 0], [1, 0]]
    Z, bl = subst(Q, [pair if i == j else [[0]] for i in range(m)])
    w = []
    for i in range(m):
        if i == j: w += [1, wbeta]
        else: w.append(wo[i])
    return Z, w

def kernel_checks(Q, j, kk, orient, wo, counters, tag):
    """kk in {2,3}: kernel type (3 stands for any k>=3).  Returns nothing; updates counters."""
    m = len(Q); wbeta = 1 if kk == 2 else 0
    Z, w = kernel(Q, j, orient, wbeta, wo)
    Zu, _ = kernel(conv(Q), j, orient, wbeta, wo)          # unchanged child, conversed quotient
    Zr, _ = kernel(conv(Q), j, 1 - orient, wbeta, wo)      # reversed child = whole converse
    g = gamma(Z, w); cu = contrast(g, gamma(Zu, w)); cr = contrast(g, gamma(Zr, w))
    u, be = j, j + 1
    verts = list(range(m + 1))
    arc = lambda a, b: Z[a][b] == 1
    hard = not any(wo[i] for i in range(m) if i != j)
    def bump(key, ok):
        counters.setdefault(tag + ':' + key, [0, 0]); counters[tag + ':' + key][0] += 1
        if not ok: counters[tag + ':' + key][1] += 1
    bump('(i) contrast nonzero', bool(cu) and bool(cr))
    bump('(ii) reversed L_u=1', Lfun(cr, arc, u, verts) == 1)
    if hard and kk == 3:
        lam = Lfun(cu, arc, u, verts) ^ Lfun(cu, arc, be, verts)
        missing = (frozenset([u, be]), frozenset([be])) if orient == 0 else (frozenset([u, be]), frozenset([u]))
        lam ^= cu.get(missing, 0)
        bump("(iii) hard unchanged k>=3 Lambda'=1+p_j", lam == (1 + endpoint_parity(Q)[j]) % 2)
    if hard and kk == 2:
        bump('(iv) hard unchanged k=2 L_u=1', Lfun(cu, arc, u, verts) == 1)
    if not hard:
        for i in range(m):
            if i != j and wo[i]:
                z = i if i < j else i + 1
                bump('(v) nonhard L_z=1 both targets', Lfun(cu, arc, z, verts) == 1 and Lfun(cr, arc, z, verts) == 1)

# ---------------------------------------------------------------- exact named coefficients (C++)
def coeffs(jobs):
    """jobs: list of (A, a0, c0, varlist). Returns list of dict (S,R)->int over variable names."""
    inp = [str(len(jobs))]
    for A, a0, c0, var in jobs:
        inp.append(str(len(A))); inp += [''.join(str(x) for x in r) for r in A]
        inp.append(''.join(str(x) for x in a0)); inp.append(''.join(str(x) for x in c0))
        inp.append(str(len(var))); inp.append(' '.join(map(str, var)))
    out = subprocess.run([EXE], input='\n'.join(inp) + '\n', capture_output=True, text=True, check=True).stdout.split('\n')
    res = []
    for k, (A, a0, c0, var) in enumerate(jobs):
        nv, d = map(int, out[2 * k].split()); vals = list(map(int, out[2 * k + 1].split()))
        D = [()] + [(var[i],) for i in range(nv)] + [(var[i], var[j]) for i, j in itertools.combinations(range(nv), 2)]
        res.append({(D[i], D[j]): vals[i * d + j] for i in range(d) for j in range(d)})
    return res

def restrict_half(CA, CB, labels):
    """Half contrast mod 2 after collapsing named variables to labels (coarse restriction)."""
    ra, rb = {}, {}
    for C, r in ((CA, ra), (CB, rb)):
        for (S, R), v in C.items():
            key = (frozenset(labels[x] for x in S), frozenset(labels[x] for x in R))
            r[key] = r.get(key, 0) + v
    out, odd = {}, 0
    for k in ra:
        x = ra[k] - rb[k]
        if x % 2: odd += 1
        if (x // 2) % 2 and (k[0] or k[1]): out[k] = 1
    return out, odd

def rand_child(k, rng): return random_tournament(k, rng) if k > 1 else [[0]]

def main():
    rng = random.Random(20261104)
    C = {}
    # ---------------- P1
    for n in (3, 4, 5):
        for Z in all_tournaments(n):
            Zo = conv(Z)
            for wbits in range(1 << n):
                w = [wbits >> i & 1 for i in range(n)]
                d = contrast(gamma(Z, w), gamma(Zo, w))
                for v in range(n):
                    key = 'P1 exhaustive n=%d' % n; C.setdefault(key, [0, 0]); C[key][0] += 1
                    if Lfun(d, lambda a, b: Z[a][b] == 1, v, range(n)) != w[v]: C[key][1] += 1
    for n in (6, 7, 8, 9):
        for t in range(300):
            Z = random_tournament(n, rng); w = [rng.randint(0, 1) for _ in range(n)]
            d = contrast(gamma(Z, w), gamma(conv(Z), w))
            for v in range(n):
                key = 'P1 random n=%d' % n; C.setdefault(key, [0, 0]); C[key][0] += 1
                if Lfun(d, lambda a, b: Z[a][b] == 1, v, range(n)) != w[v]: C[key][1] += 1
    print('P1 done', flush=True)
    # ---------------- P2
    for m in (2, 3, 4, 5):
        for Q in all_tournaments(m):
            for j in range(m):
                if m <= 4:
                    wos = [[b >> i & 1 for i in range(m)] for b in range(1 << m) if not (b >> j & 1)]
                else:
                    wos = [[0] * m] + [[1 if i == z else 0 for i in range(m)] for z in range(m) if z != j]
                for wo in wos:
                    for kk in (2, 3):
                        for orient in (0, 1):
                            kernel_checks(Q, j, kk, orient, wo, C, 'P2 m=%d' % m)
        print('P2 m=%d done' % m, flush=True)
    # ---------------- P3
    for t in range(400):
        m = rng.randint(6, 9); Q = random_tournament(m, rng); j = rng.randrange(m)
        wo = [0] * m if t % 2 == 0 else [rng.randint(0, 1) if i != j else 0 for i in range(m)]
        kernel_checks(Q, j, rng.choice((2, 3)), rng.choice((0, 1)), wo, C, 'P3')
    print('P3 done', flush=True)
    # ---------------- P4: S188 prime-eight endpoint obstruction quotient
    rows = (238, 4, 240, 118, 195, 146, 34, 74)
    Q8 = [[rows[i] >> j & 1 for j in range(8)] for i in range(8)]
    p8 = endpoint_parity(Q8)
    C['P4 Q8 prime'] = [1, 0 if is_prime(Q8) else 1]
    C['P4 Q8 endpoint parity zero'] = [1, 0 if not any(p8) else 1]
    for j in range(8):
        for kk in (2, 3):
            for orient in (0, 1):
                kernel_checks(Q8, j, kk, orient, [0] * 8, C, 'P4 Q8 kernel')
        for k in (2, 3, 4, 5):     # full named TT_k child, other children nontrivial (weight 0)
            for eps in (0, 1):
                kids = [(TT(k) if i == j else [[0]]) for i in range(8)]
                Z, bl = subst(Q8, kids); w = [1 if v in bl[j] else 0 for v in range(len(Z))]
                kid2 = [((conv(TT(k)) if eps else TT(k)) if i == j else [[0]]) for i in range(8)]
                Zp, _ = subst(conv(Q8), kid2)
                d = contrast(gamma(Z, w), gamma(Zp, w))
                key = 'P4 Q8 named TT_k both targets nonzero'; C.setdefault(key, [0, 0]); C[key][0] += 1
                if not d: C[key][1] += 1
        # coarse-only (child j contracted to weight 0 too): predicted zero (S188 blindness)
        d = contrast(gamma(Q8, [0] * 8), gamma(conv(Q8), [0] * 8))
        key = 'P4 Q8 coarse contrast zero'; C.setdefault(key, [0, 0]); C[key][0] += 1
        if d: C[key][1] += 1
    print('P4 done', flush=True)
    # ---------------- P5 actual graphs (exact C++)
    def rand_prime(m):
        while True:
            Q = random_tournament(m, rng)
            if is_prime(Q): return Q
    # (a) singleton child, coarse
    for t in range(12):
        m = rng.choice((5, 6)); Q = rand_prime(m)
        while True:
            sizes = [1] + [rng.choice((1, 2, 3)) for _ in range(m - 1)]
            rng.shuffle(sizes)
            if sum(sizes) <= 11: break
        W, bl = subst(Q, [rand_child(k, rng) for k in sizes]); Wp, _ = subst(conv(Q), [rand_child(k, rng) for k in sizes])
        n = len(W); var = list(range(n))
        CW, CWp = coeffs([(W, [0] * n, [0] * n, var), (Wp, [0] * n, [0] * n, var)])
        lab = {v: i for i, b in enumerate(bl) for v in b}
        d, odd = restrict_half(CW, CWp, lab)
        for z in range(m):
            if sizes[z] == 1:
                key = 'P5a singleton coarse L_z=1'; C.setdefault(key, [0, 0]); C[key][0] += 1
                if Lfun(d, lambda a, b: Q[a][b] == 1, z, range(m)) != 1: C[key][1] += 1
        key = 'P5a contrasts even'; C.setdefault(key, [0, 0]); C[key][0] += 1
        if odd: C[key][1] += 1
    # (b) transitive child kernel on actual graphs, (c) ambient controls
    amb_done = 0
    for t in range(12):
        m = rng.choice((4, 5)); Q = random_tournament(m, rng); j = rng.randrange(m); k = rng.choice((3, 4))
        while True:
            sizes = [k if i == j else rng.choice((2, 3)) for i in range(m)]
            if sum(sizes) <= 11: break
        kids = [TT(k) if i == j else rand_child(sizes[i], rng) for i in range(m)]
        W, bl = subst(Q, kids); n = len(W)
        pj = endpoint_parity(Q)[j]
        for eps in (0, 1):
            kids2 = [((conv(TT(k)) if eps else TT(k)) if i == j else rand_child(sizes[i], rng)) for i in range(m)]
            Wp, _ = subst(conv(Q), kids2)
            var = list(range(n))
            CW, CWp = coeffs([(W, [0] * n, [0] * n, var), (Wp, [0] * n, [0] * n, var)])
            lab = {}
            for i, b in enumerate(bl):
                for v in b: lab[v] = i
            lab[bl[j][0]] = 'u'
            for v in bl[j][1:]: lab[v] = 'beta'
            d, odd = restrict_half(CW, CWp, lab)
            verts = [i for i in range(m) if i != j] + ['u', 'beta']
            def arcZ(a, b):
                qa = j if a in ('u', 'beta') else a; qb = j if b in ('u', 'beta') else b
                if qa == qb: return (a, b) == ('u', 'beta')
                return Q[qa][qb] == 1
            key = 'P5b nonzero'; C.setdefault(key, [0, 0]); C[key][0] += 1
            if not d: C[key][1] += 1
            if eps == 1:
                key = 'P5b reversed L_u=1'; C.setdefault(key, [0, 0]); C[key][0] += 1
                if Lfun(d, arcZ, 'u', verts) != 1: C[key][1] += 1
            else:
                lam = Lfun(d, arcZ, 'u', verts) ^ Lfun(d, arcZ, 'beta', verts) ^ d.get((frozenset(['u', 'beta']), frozenset(['beta'])), 0)
                key = "P5b unchanged Lambda'=1+p_j"; C.setdefault(key, [0, 0]); C[key][0] += 1
                if lam != (1 + pj) % 2: C[key][1] += 1
            # (c) ambient control: embed the core as a module in an exterior of 2-3 vertices
            if amb_done < 4 and n <= 8:
                amb_done += 1 if eps == 1 else 0
                r = rng.choice((2, 3)); R = random_tournament(r + 1, rng)   # vertex 0 of R = the core
                TA, blA = subst(R, [W] + [[[0]]] * r); TB, _ = subst(R, [Wp] + [[[0]]] * r)
                N = len(TA); core = blA[0]
                a0 = [0] * N; c0 = [0] * N
                for e in range(1, r + 1):
                    for v in blA[e]:
                        a0[v] = R[e][0]; c0[v] = R[e][0]    # ports imitate the core's exterior incidences
                CA, CB = coeffs([(TA, a0, c0, core), (TB, a0, c0, core)])
                dA, _ = restrict_half(CA, CB, lab)
                key = 'P5c ambient restricted mod-4 contrast == bare core'; C.setdefault(key, [0, 0]); C[key][0] += 1
                if dA != d: C[key][1] += 1
    print('P5 done', flush=True)
    json.dump(C, open(os.path.join(HERE, 'confirm.json'), 'w'), indent=1, sort_keys=True)
    bad = sum(v[1] for v in C.values())
    for k in sorted(C): print('%-60s checks=%7d failures=%d' % (k, C[k][0], C[k][1]))
    print('TOTAL FAILURES', bad)

if __name__ == '__main__':
    main()
