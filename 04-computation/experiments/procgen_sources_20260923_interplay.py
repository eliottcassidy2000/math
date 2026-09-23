#!/usr/bin/env python3
"""procgen_sources_20260923_interplay.py

Small tests behind the typed interplay map of
05-knowledge/results/procgen_sources_20260923_interplay.md.

  T1  Erdos 1062(ii) (Kruer-Kohlmeyer exposition): the finite formula (4) against
      an exact ILP brute force of f(n) for n <= 70; f(n)/n versus the explicit
      density L (7); the coefficient words a_k, b_k are codings of the rotations by
      log_3 4 and log_4 3 (the critical Collatz slope, doubled); repetition estimates.
  T2  Theorem Y in q-Pochhammer and Tate-curve form (2-adic, exact mod 2^N):
      Jacobi triple product for theta_3(rho); the square-swap and pronic-swap
      Bernstein numbers from their parity words; j(rho^2) computed twice
      (q-expansion with c(1) = 196884, and 256(1-l+l^2)^3/(l^2(1-l)^2) with the
      2-adic theta constants), rho = 2^10/3^9.
  T3  Radchenko-Wheeler: the N = 5 finite quantum dilogarithm E_a satisfies the
      defining identity (29) with a nonzero constant; Theorem 5 (i), (ii).
  T4  Erdos 859: Ford's delta against the session's large-deviation rates; PSLQ
      search for small integer relations.
  T5  Bounded-discrepancy words: Terras bijection (every finite parity word is
      realized by exactly one residue class mod 2^s) and strip-word counts.
  T6  The Out(S_6) twist on degree-6 symmetric functions applied to the
      Redei-Berge functions U_T of all 2^15 tournaments on 6 vertices.
  T7  Hexagon numerology: centered hexagonal numbers that are Phi_6 values
      (Pell equation X^2 - 3Y^2 = -2).

Deterministic; one process; peak memory < 200 MB.
"""
import math
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations

import mpmath
import numpy as np


def say(*a):
    print(*a)
    sys.stdout.flush()


# ------------------------------------------------------------------ T1
def ilog(b, t):
    k = 0
    p = b
    while p <= t:
        p *= b
        k += 1
    return k


def R_count(t):
    if t <= 0:
        return 0
    l2, l3 = ilog(2, t), ilog(3, t)
    eta = 1 if (l2 % 2 == 1 and 2 * 3 ** l3 <= t) else 0
    return l3 + 1 + l2 // 2 + eta


def B_corr(t):
    if t <= 0:
        return 0
    j = 0
    ok4 = False
    while 4 ** j <= t:
        if 4 * t < 5 * 4 ** j:
            ok4 = True
        j += 1
    if not ok4:
        return 0
    b = 0
    while 10 * 3 ** b <= t:
        if t < 12 * 3 ** b:
            return 1
        b += 1
    return 0


def f_formula(n, gcache={}):
    tot = 0
    for q in range(1, n + 1):
        if q % 2 and q % 3:
            t = n // q
            if t not in gcache:
                gcache[t] = R_count(t) - B_corr(t)
            tot += gcache[t]
    return tot


def f_bruteforce(n):
    """max |A|, A subset of [1..n], no a in A with two distinct proper multiples in A (exact ILP)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    if n <= 2:
        return n
    rows = []
    for a in range(1, n + 1):
        mult = list(range(2 * a, n + 1, a))
        for b, c in combinations(mult, 2):
            r = np.zeros(n)
            r[a - 1] = r[b - 1] = r[c - 1] = 1
            rows.append(r)
    A = np.array(rows) if rows else np.zeros((0, n))
    res = milp(c=-np.ones(n), constraints=[LinearConstraint(A, -np.inf, 2)] if rows else [],
               integrality=np.ones(n), bounds=Bounds(0, 1))
    return int(round(-res.fun))


def coeff_a(k):
    x = 4 ** k
    y = 3 ** ilog(3, x)
    for cond, val in ((15 * x < 16 * y, 78), (9 * x < 10 * y, 30), (3 * x < 4 * y, -30),
                      (2 * x < 3 * y, 30), (x < 2 * y, 0), (3 * x < 8 * y, -60)):
        if cond:
            return val
    return -12


def coeff_b(k):
    u = 3 ** k
    v = 4 ** ilog(4, u)
    for cond, val in ((3 * u < 4 * v, 30), (5 * u < 8 * v, 35), (3 * u < 5 * v, 29), (u < 2 * v, 24)):
        if cond:
            return val
    return -60


def section_T1():
    say('=' * 78)
    say('T1. Erdos 1062(ii): formula (4), brute force, the density, and the rotation words')
    mism = []
    for n in range(1, 71):
        fb = f_bruteforce(n)
        ff = f_formula(n)
        if fb != ff:
            mism.append((n, fb, ff))
    say(f'  exact ILP brute force of f(n) = formula (4) for all n <= 70: {not mism}' +
        (f'  mismatches: {mism[:5]}' if mism else ''))
    say(f'  f(n), n = 10..20 (brute force): {[f_bruteforce(n) for n in range(10, 21)]}')
    mpmath.mp.dps = 50
    K = 130
    A = [coeff_a(k) for k in range(K)]
    B = [coeff_b(k) for k in range(K)]
    L = (mpmath.mpf(32) + sum(mpmath.mpf(A[k]) / mpmath.mpf(4) ** k for k in range(K))
         + sum(mpmath.mpf(B[k]) / mpmath.mpf(3) ** k for k in range(K))) / 180
    say(f'  L_explicit = (32 + sum a_k/4^k + sum b_k/3^k)/180 = {mpmath.nstr(L, 40)} '
        f'(truncation error < 1e-58)')
    for n in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
        fn = f_formula(n)
        say(f'  n = {n:>8}: f(n)/n = {fn / n:.10f}   (f(n)/n - L = {fn / n - float(L):+.2e})')
    # rotation codings
    al = mpmath.log(4) / mpmath.log(3)
    be = mpmath.log(3) / mpmath.log(4)

    def arcA(fr):
        r = mpmath.mpf(3) ** fr
        for lim, val in ((mpmath.mpf(16) / 15, 78), (mpmath.mpf(10) / 9, 30), (mpmath.mpf(4) / 3, -30),
                         (mpmath.mpf(3) / 2, 30), (mpmath.mpf(2), 0), (mpmath.mpf(8) / 3, -60)):
            if r < lim - mpmath.mpf(10) ** -40:  # an exact endpoint belongs to the next row
                return val
        return -12

    def arcB(fr):
        r = mpmath.mpf(4) ** fr
        for lim, val in ((mpmath.mpf(4) / 3, 30), (mpmath.mpf(8) / 5, 35), (mpmath.mpf(5) / 3, 29),
                         (mpmath.mpf(2), 24)):
            if r < lim - mpmath.mpf(10) ** -40:
                return val
        return -60
    KK = 3000
    okA = all(coeff_a(k) == arcA(k * al - mpmath.floor(k * al)) for k in range(KK))
    okB = all(coeff_b(k) == arcB(k * be - mpmath.floor(k * be)) for k in range(KK))
    say(f'  a_k = A(frac(k log_3 4)) (arc endpoints log_3 of 16/15, 10/9, 4/3, 3/2, 2, 8/3) for k < {KK}: {okA}'
        f' (k = 1 hits the endpoint 4/3 exactly; the tables send it to the next row)')
    say(f'  b_k = B(frac(k log_4 3)) (arc endpoints log_4 of 4/3, 8/5, 5/3, 2) for k < {KK}: {okB}')
    say(f'  log_3 4 = 2 log_3 2 = {mpmath.nstr(al, 12)}: twice the critical 3x+1 slope log_3 2')
    wa = [coeff_a(k) for k in range(4000)]
    wb = [coeff_b(k) for k in range(4000)]
    for name, w in (('a', wa), ('b', wb)):
        comp = [len({tuple(w[i:i + m]) for i in range(len(w) - m)}) for m in (1, 2, 4, 8, 16, 32)]
        say(f'  word {name}: factor complexity p(1,2,4,8,16,32) = {comp} (linear: rotation coding)')


# ------------------------------------------------------------------ T2
def section_T2():
    say('=' * 78)
    say('T2. Theorem Y as a 2-adic theta / q-Pochhammer value; the Tate curve with q = rho^2')
    N = 4000
    M = 1 << N
    inv3 = pow(3, -1, M)
    rho = (1 << 10) * pow(inv3, 9, M) % M
    # theta_3(rho) = sum_{k in Z} rho^{k^2}
    th = 1
    k = 1
    while 10 * k * k < N:
        th = (th + 2 * pow(rho, k * k, M)) % M
        k += 1
    # Jacobi triple product: prod_{m>=1} (1 - rho^{2m}) (1 + rho^{2m-1})^2
    pr = 1
    m = 1
    while 10 * (2 * m - 1) < N:
        a = (1 - pow(rho, 2 * m, M)) % M
        b = (1 + pow(rho, 2 * m - 1, M)) % M
        pr = pr * a % M * b % M * b % M
        m += 1
    say(f'  Jacobi triple product: sum_(k in Z) rho^(k^2) = prod (1-rho^(2m))(1+rho^(2m-1))^2 mod 2^{N}: {th == pr}')

    def bernstein_from_swapset(S_pred, nbits):
        """-sum_l 2^{d_l} 3^{-l} over the ones of the block word (mod 2^nbits)."""
        Mb = 1 << nbits
        i3 = pow(3, -1, Mb)
        tot = 0
        pos = 0
        l = 0
        blk = 0
        p3 = 1
        while pos < nbits:
            block = [1] * 8 + [0, 1] if (blk > 0 and S_pred(blk)) else [1] * 9 + [0]
            for bit in block:
                if bit:
                    l += 1
                    p3 = p3 * i3 % Mb
                    if pos < nbits:
                        tot = (tot + pow(2, pos, Mb) * p3) % Mb
                pos += 1
            blk += 1
        return (-tot) % Mb

    NB = 3000
    Mb = 1 << NB
    rb = rho % Mb
    c0 = (-1 - 512 * pow(3 ** 9 - 2 ** 10, -1, Mb)) % Mb
    c1 = 256 * pow(3 ** 9, -1, Mb) % Mb
    issq = lambda n: int(math.isqrt(n)) ** 2 == n
    ispron = lambda n: (lambda r: r * (r + 1) == n)(int(math.isqrt(n)))
    sq = sum(pow(rb, k * k, Mb) for k in range(1, int(math.isqrt(NB // 10)) + 2)) % Mb
    prn = sum(pow(rb, k * (k + 1), Mb) for k in range(1, int(math.isqrt(NB // 10)) + 2)) % Mb
    phiY = bernstein_from_swapset(issq, NB)
    phiP = bernstein_from_swapset(ispron, NB)
    say(f'  square-swap word:  Phi_2(Y)    = -1 - 512/(3^9-2^10) - (256/3^9) sum_(k>=1) rho^(k^2)    mod 2^{NB}: '
        f'{phiY == (c0 - c1 * sq) % Mb}')
    say(f'  pronic-swap word:  Phi_2(Y_pr) = -1 - 512/(3^9-2^10) - (256/3^9) sum_(k>=1) rho^(k(k+1)) mod 2^{NB}: '
        f'{phiP == (c0 - c1 * prn) % Mb}')
    # Tate curve: q = rho^2, j(q) q = E4^3 / prod (1-q^n)^24
    nterms = N // 20 + 3
    sig3 = [0] + [sum(d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, nterms + 1)]
    E4 = [1] + [240 * sig3[n] for n in range(1, nterms + 1)]

    def smul(a, b, L):
        r = [0] * (L + 1)
        for i, x in enumerate(a[:L + 1]):
            if x:
                for jdx, y in enumerate(b[:L + 1 - i]):
                    r[i + jdx] += x * y
        return r
    P = [0] * (nterms + 1)
    P[0] = 1
    for n in range(1, nterms + 1):
        for _ in range(24):
            for kk in range(nterms, n - 1, -1):
                P[kk] -= P[kk - n]
    invP = [0] * (nterms + 1)
    invP[0] = 1
    for kk in range(1, nterms + 1):
        invP[kk] = -sum(P[i] * invP[kk - i] for i in range(1, kk + 1))
    jq_coeffs = smul(smul(smul(E4, E4, nterms), E4, nterms), invP, nterms)
    say(f'  j q = 1 + {jq_coeffs[1]} q + {jq_coeffs[2]} q^2 + {jq_coeffs[3]} q^3 + ... (integer coefficients)')
    q = rho * rho % M
    jq1 = 0
    qp = 1
    for n in range(nterms + 1):
        jq1 = (jq1 + jq_coeffs[n] * qp) % M
        qp = qp * q % M
    # lambda route
    psi = 0
    n = 0
    while 10 * n * (n + 1) < N + 40:
        psi = (psi + pow(rho, n * (n + 1), M)) % M
        n += 1
    U = pow(inv3, 9, M) * pow(psi, 4, M) % M * pow(pow(th, 4, M), -1, M) % M   # lambda = 2^14 U
    lam = (1 << 14) * U % M
    num = pow((1 - lam + lam * lam) % M, 3, M)
    den = U * U % M * pow((1 - lam) % M, 2, M) % M
    jq2 = pow(inv3, 18, M) * num % M * pow(den, -1, M) % M
    Nc = N - 60
    say(f'  j(rho^2) * rho^2 from the q-expansion  ==  from lambda = theta_2^4/theta_3^4 (mod 2^{Nc}): '
        f'{(jq1 - jq2) % (1 << Nc) == 0}')
    say(f'  j(rho^2) * rho^2 is a 2-adic unit: {jq1 % 2 == 1}  =>  v_2(j(rho^2)) = -20; '
        f'v_2(196884 q) = {(196884 & -196884).bit_length() - 1} + 20')
    say('  BDGP 1996 (p-adic Mahler-Manin, CITED): j(rho^2) is transcendental; lambda is a rational function')
    say('  of theta_3(rho) and psi(rho^2) over Q, so theta_3(rho) and psi(rho^2) are not both algebraic,')
    say('  i.e. the square-swap and pronic-swap Bernstein numbers are not both algebraic over Q.')
    # q-hypergeometric test: c_{k+1}/c_k for squares vs cubes
    say('  coefficient ratios: squares rho^((k+1)^2-k^2) = rho * (rho^2)^k (q-hypergeometric, q = rho^2);')
    say('  cubes rho^(3k^2+3k+1): exponent quadratic in k, not of the form R(q^k) (no first-order q-difference eq.)')


# ------------------------------------------------------------------ T3
def section_T3():
    say('=' * 78)
    say('T3. Radchenko-Wheeler: the N = 5 finite quantum dilogarithm E_a')
    mpmath.mp.dps = 40
    Nn = 5
    z10 = mpmath.exp(2j * mpmath.pi / 10)
    E = {0: (3 + mpmath.sqrt(5)) / 2, 1: z10, 4: z10, 2: z10 ** -1, 3: z10 ** -1}
    gauss = {x: z10 ** ((3 * x * (x + 5)) % 10) for x in range(5)}
    pair = lambda x, y: gauss[(x + y) % 5] / (gauss[x] * gauss[y])
    # bicharacter check
    bich = all(abs(pair(x, (y + z) % 5) - pair(x, y) * pair(x, z)) < 1e-30
               for x in range(5) for y in range(5) for z in range(5))
    Cs = []
    for x in range(5):
        for y in range(5):
            lhs = pair(x, y) * E[x] * E[y]
            integ = sum(E[(y - z) % 5] * E[z] * gauss[z] * E[(x - z) % 5] for z in range(5)) / mpmath.sqrt(5)
            Cs.append(lhs - integ)
    spread = max(abs(c - Cs[0]) for c in Cs)
    say(f'  <x;y> = <x+y>/(<x><y>) is a bicharacter: {bich}')
    say(f'  (29): <x;y>E(x)E(y) - (1/sqrt5) sum_z E(y-z)E(z)<z>E(x-z) is constant over all (x,y): '
        f'spread {mpmath.nstr(spread, 3)}; C = {mpmath.nstr(Cs[0], 15)} (nonzero)')
    e0 = E[0]
    say(f'  Theorem 5(i): E(0)^2 = sqrt5 E(0) + 1: {abs(e0 ** 2 - mpmath.sqrt(5) * e0 - 1) < 1e-30}; '
        f'E(0) = phi^2 = eps^(1/2), eps = phi^4 (trace N+2 = 7)')
    ok2 = all(abs(E[u] * E[(-u) % 5] - 1 / gauss[u]) < 1e-30 for u in range(1, 5))
    say(f'  Theorem 5(ii): E(u)E(-u) = <u>^-1 for u != 0: {ok2}')
    say(f'  numerology: the AMM 12592 constant 1 + 2 log(phi)/log 5 = 1 + log_5 E_a(0) = '
        f'{mpmath.nstr(1 + mpmath.log(e0) / mpmath.log(5), 15)}')


# ------------------------------------------------------------------ T4
def section_T4():
    say('=' * 78)
    say("T4. Erdos 859 (Ford's delta) against the session's large-deviation rates")
    mpmath.mp.dps = 40
    ln2, ln3 = mpmath.log(2), mpmath.log(3)
    delta = 1 - (1 + mpmath.log(ln2)) / ln2
    lam = 1 / ln2
    Q = lam * mpmath.log(lam) - lam + 1
    p = ln2 / ln3
    Hn = -p * mpmath.log(p) - (1 - p) * mpmath.log(1 - p)
    D = ln2 - Hn
    cstar = mpmath.log(mpmath.mpf(128) / 81) / 4
    IG = -mpmath.log(mpmath.mpf('0.758751'))
    say(f'  delta = 1 - (1+ln ln 2)/ln 2 = {mpmath.nstr(delta, 20)} = Q(1/ln 2), Q(l) = l ln l - l + 1: '
        f'{abs(delta - Q) < 1e-35}')
    say(f'  Collatz critical rate D(log_3 2 || 1/2) = ln 2 - H(log_3 2) = {mpmath.nstr(D, 12)} nats '
        f'= {mpmath.nstr(D / ln2, 12)} bits (= 1 - h(log_3 2), h = {mpmath.nstr(Hn / ln2, 8)})')
    say(f'  greedy-G rate -ln 0.758751 = {mpmath.nstr(IG, 10)}; budget c* = ln(128/81)/4 = {mpmath.nstr(cstar, 12)}')
    say(f'  ratios: delta/D = {mpmath.nstr(delta / D, 8)}, delta/c* = {mpmath.nstr(delta / cstar, 8)}')
    rel = mpmath.pslq([delta, D, cstar, 1], maxcoeff=10 ** 4, maxsteps=10 ** 6)
    say(f'  PSLQ integer relation among (delta, D, c*, 1) with |coeff| <= 10^4: {rel}')
    rel2 = mpmath.pslq([delta * ln2, 1, ln2, mpmath.log(ln2)], maxcoeff=100, maxsteps=10 ** 5)
    say(f'  sanity: PSLQ on (delta ln2, 1, ln2, ln ln2) finds {rel2} (delta ln 2 = ln 2 - 1 - ln ln 2)')
    # both are relative entropies at a 'log 2' threshold
    say('  structure: delta is the Poisson rate at mean-multiple 1/ln 2 (2^omega = log scale);')
    say('  D is the Bernoulli(1/2) rate at frequency ln2/ln3 (3^a = 2^s). Same shape, different laws.')


# ------------------------------------------------------------------ T5
def section_T5():
    say('=' * 78)
    say('T5. Bounded-discrepancy words: every finite word is realized (Terras), strips have positive entropy')
    for s in (8, 12, 16):
        words = set()
        for n in range(1 << s):
            x = n
            w = 0
            for i in range(s):
                bit = x & 1
                w |= bit << i
                x = (3 * x + 1) >> 1 if bit else x >> 1
            words.add(w)
        say(f'  s = {s:2d}: residues mod 2^s -> first s parity bits of T(x) = x/2, (3x+1)/2 is a bijection: '
            f'{len(words) == 1 << s}')
    for alpha, W in ((Fraction(9, 10), Fraction(1)), (Fraction(9, 10), Fraction(3, 2)), (Fraction(3, 4), Fraction(2))):
        # count words with c <= a_t - alpha t <= c + W for all t <= s, c = -W/2 (centred strip)
        c = -W / 2
        states = Counter({0: 1})  # key: a_t (number of ones)
        counts = []
        for t in range(1, 61):
            new = Counter()
            for a, m in states.items():
                for bit in (0, 1):
                    a2 = a + bit
                    dev = a2 - alpha * t
                    if c <= dev <= c + W:
                        new[a2] += m
            states = new
            counts.append(sum(states.values()))
        rate = math.log2(counts[59]) - math.log2(counts[39])
        say(f'  strip slope {alpha}, width {W}: #words of length 60 = {counts[59]}, growth '
            f'{rate / 20:.3f} bits/letter; each is realized by exactly one class mod 2^60 (Terras bijection)')


# ------------------------------------------------------------------ T6
def section_T6():
    say('=' * 78)
    say('T6. The Out(S_6) twist on Lambda^6 applied to Redei-Berge functions of 6-vertex tournaments')
    V = range(6)
    pairs = list(combinations(V, 2))
    triples = Counter()
    fixed = 0
    realized = set()
    stats = []
    for mask in range(1 << 15):
        adj = [[False] * 6 for _ in V]
        for bitpos, (i, j) in enumerate(pairs):
            if mask >> bitpos & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        # directed 3-cycles (as vertex sets; each 3-set carries at most one)
        c3sets = [S for S in combinations(V, 3)
                  if (adj[S[0]][S[1]] and adj[S[1]][S[2]] and adj[S[2]][S[0]])
                  or (adj[S[0]][S[2]] and adj[S[2]][S[1]] and adj[S[1]][S[0]])]
        c3 = len(c3sets)
        d33 = sum(1 for A, B in combinations(c3sets, 2) if not set(A) & set(B))
        c5 = 0
        for S in combinations(V, 5):
            first = S[0]
            for perm in permutations(S[1:]):
                cyc = (first,) + perm
                if all(adj[cyc[i]][cyc[(i + 1) % 5]] for i in range(5)):
                    c5 += 1
        # Hamiltonian paths by DP
        dp = [[0] * 6 for _ in range(1 << 6)]
        for v in V:
            dp[1 << v][v] = 1
        for S in range(1 << 6):
            for v in V:
                if dp[S][v]:
                    for u in V:
                        if not S >> u & 1 and adj[v][u]:
                            dp[S | 1 << u][u] += dp[S][v]
        H = sum(dp[63])
        stats.append((c3, c5, d33, H))
        realized.add((c3, c5, d33))
    ocf_ok = all(H == 1 + 2 * c3 + 2 * c5 + 4 * d33 for c3, c5, d33, H in stats)
    say(f'  OCF check on all 32768 tournaments: H(T) = 1 + 2 c3 + 2 c5 + 4 d33: {ocf_ok}')
    say(f'  distinct (c3, c5, d33) triples: {len(realized)}: {sorted(realized)}')
    img_ok = [t for t in realized if t[0] % 2 == 0 and (2 * t[2], t[1], t[0] // 2) in realized]
    fixed = [t for t in realized if t[0] == 2 * t[2]]
    say('  U_T = p1^6 + 2 c3 p3 p1^3 + 2 c5 p5 p1 + 4 d33 p3^2 (Grinberg-Stanley, n = 6);')
    say('  tau (from Out(S_6)) swaps p_{3,1,1,1} <-> p_{3,3} and fixes p_{1^6}, p_{5,1}:')
    say('  tau(U_T) = U_T\' requires (c3\', c5\', d33\') = (2 d33, c5, c3/2).')
    say(f'  triples whose tau-image is again a Redei-Berge function: {len(img_ok)} of {len(realized)}: {sorted(img_ok)}')
    say(f'  tau-fixed triples (c3 = 2 d33): {sorted(fixed)}')
    cnt_fixed = sum(1 for c3, c5, d33, H in stats if c3 == 2 * d33)
    say(f'  labelled tournaments with tau(U_T) = U_T: {cnt_fixed} of 32768; tau preserves H (zeta(p_lambda) = 1)')


# ------------------------------------------------------------------ T7
def section_T7():
    say('=' * 78)
    say('T7. Hexagon numerology: 3k(k-1)+1 = n^2 - n + 1 = Phi_6(n)  <=>  (2n-1)^2 - 3(2k-1)^2 = -2')
    sols = []
    for k in range(1, 10 ** 6):
        H = 3 * k * (k - 1) + 1
        n = (1 + math.isqrt(4 * H - 3)) // 2
        if n * n - n + 1 == H:
            sols.append((k, n, H))
    say(f'  (k, n, N) with N = 3k(k-1)+1 = Phi_6(n), k < 10^6: {sols}')
    X = [(2 * n - 1, 2 * k - 1) for k, n, H in sols]
    say(f'  Pell pairs (X, Y) = (2n-1, 2k-1): {X}; all satisfy X^2 - 3Y^2 = -2: '
        f'{all(x * x - 3 * y * y == -2 for x, y in X)}')
    say('  the 7-flower (k = 2) is Phi_6(3) = 7; LRC14 witness modulus Phi_6(14) = 183 is not centered hexagonal: '
        f'{not any(H == 183 for _, _, H in sols)}')
    # the flower centres {0, +-1, +-w, +-w^2} as residues of Z[w]/(3+w) = Z/7 (w -> -3)
    wmap = {'0': 0, '1': 1, '-1': -1, 'w': -3, '-w': 3, 'w2': 9, '-w2': -9}
    res = sorted(v % 7 for v in wmap.values())
    say(f'  flower centres mod (3 + w) (w = e^(2 pi i/3) -> -3 in Z/7): residues {res}; complete system: {res == list(range(7))}')
    petals = [1, -9, -3, -1, 9, 3]  # 1, -w^2, w, -1, w^2, -w in counterclockwise order (angles 0, 60, ..., 300)
    rot = all((3 * petals[i]) % 7 == petals[(i + 1) % 6] % 7 or (3 * petals[i]) % 7 == petals[(i - 1) % 6] % 7
              for i in range(6))
    say(f'  multiplication by 3 = -w on Z/7 rotates the six petals by 60 degrees (a zeta_6 action): {rot}; '
        f'order of 3 mod 7 = {next(k for k in range(1, 7) if pow(3, k, 7) == 1)}')


def main():
    say('procgen_sources_20260923_interplay.py')
    section_T1()
    section_T2()
    section_T3()
    section_T4()
    section_T5()
    section_T6()
    section_T7()
    say('=' * 78)
    say('interplay: done')


if __name__ == '__main__':
    main()
