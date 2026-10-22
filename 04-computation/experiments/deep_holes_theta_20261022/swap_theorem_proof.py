r"""PROOF of the swap theorem and of HYP-9363 for all 23 deep-hole types (session opus-2026-10-22-S53).

Setting.
- II = II_{25,1}. rho is the Leech cusp (Weyl vector). z is the cusp of a Niemeier lattice N, with z . rho = -h.
- M = {rho, z}^perp has discriminant group D = M^*/M = (Z/h)^2. D is hyperbolic: it contains two isotropic lines,
  H_L (from Lambda) and H_N (from N).
- sigma is the swap of the two lines. It lies in O(D).
- Theta_M = sum_g theta_(M+g) e_g is a modular form of weight 12 for the Weil representation rho_M (Borcherds 1998,
  Thm 4.1; the signature is 24, so the eighth-root phase is 1).
- Since O(D) commutes with rho_M, sigma^* Theta_M = sum_g theta_(M + sigma g) e_g is also one.

Lemma (vector-valued Sturm bound). Let rho be unitary with finite image and let F be in M_k(rho). If every component
of F vanishes to order > k/12 at infinity, then F = 0.
Proof. Each component of F|g is a linear combination of components of F. So for each j the product N_j over
Gamma_rho \ SL_2(Z) of F_j|_k g is a scalar modular form of weight k [SL_2(Z) : Gamma_rho] for SL_2(Z). Its order at
infinity exceeds its weight divided by 12, so by the valence formula N_j = 0. Hence some F_j|g = 0, and F_j = 0.

Proof of the swap theorem for one type: a finite check, carried out exactly below.
1. Enumerate the vectors of M^* with |x|^2 <= 2 and record each with its class in D. (The number of nonzero
   vectors is cross-checked against PARI's qfminim.)
2. Put E = Theta_M - sigma^* Theta_M. Check that E has no terms of exponent < 1, i.e. no differences at norms < 2.
3. Let v_g be the q^1 coefficient of E_g, i.e. the norm-2 difference. Check that v is rho_M-invariant:
   - v is supported on isotropic classes, so it is T-invariant;
   - (1/h) sum_d e(-b(g, d)) v_d = v_g holds exactly in Z[zeta_h], so it is S-invariant.
4. Then F = E - Delta v lies in M_12(rho_M) and vanishes to order > 1, so F = 0. That is,

        Theta_M - sigma^* Theta_M = Delta . v          (the swap theorem).

Consequence (HYP-9363). theta_(Lambda + c) is the sum of theta_(M+g) over the line A, which consists of the classes of
projections of {x : x . rho = -1}; theta_(N + rho_N/h) is the sum over B = sigma(A). Hence

        theta_(Lambda + c) - theta_(N + rho_N/h) = (sum_{g in A} v_g) Delta,

and the script checks that sum_{g in A} v_g = 24 + k, the number of hole vertices, read from Borcherds' table.
It also reports the decomposition v = sum_d c_d 1_{S_d} over the maximal isotropic subgroups S_d.
Needs cypari, numba and sympy. Prints ALL CHECKS PASSED."""
import sys, time, math
import numpy as np
from fractions import Fraction
from cypari import pari
from numba import njit
from sympy import Poly, symbols, cyclotomic_poly, rem, ZZ

ONLY = sys.argv[1:]
pari.allocatemem(2 * 10 ** 9)
T0 = time.time()
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('swap_theorem_proof.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = []
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows.append((t[4], int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]]))

gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
Am = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(Am)
basis = [[int(Hb[i, j]) for i in range(n + 1)] for j in range(n + 1)]
rho2 = [2 * i for i in range(n)] + [140]
r0 = [2, -2] + [0] * (n - 1)                      # e_0 - e_1 (doubled): a Leech root, r0 . rho = -1


def kernel(forms, dim):
    F = pari.matrix(len(forms), n + 1, [lor2(b, f) // 4 for f in forms for b in basis])
    K = pari.matkerint(F)
    return [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(dim)]


@njit(cache=True)
def enum_dual(Q, Bnd):
    d = Q.shape[0]
    q = Q.astype(np.float64).copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    out = []
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    rad = np.sqrt(Bnd / q[i, i])
    a[i] = np.int64(np.ceil(-rad - 1e-9))
    hi[i] = np.int64(np.floor(rad + 1e-9))
    while True:
        if a[i] > hi[i]:
            i += 1
            if i == d:
                break
            a[i] += 1
            continue
        t = a[i] - cen[i]
        T[i] = T[i + 1] + q[i, i] * t * t
        if T[i] > Bnd * (1 + 1e-9) + 1e-9:
            a[i] += 1
            continue
        if i == 0:
            out.append(a.copy())
            a[i] += 1
            continue
        j = i - 1
        s = 0.0
        for l in range(j + 1, d):
            s += q[j, l] * a[l]
        cen[j] = -s
        rad = np.sqrt(max(Bnd * (1 + 1e-9) + 1e-9 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return out


X = symbols('X')
INTER = {}
for name, h, n0, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    t1 = time.time()
    Mv = kernel([rho2, z2], 24)
    GM = [[Fraction(lor2(a, b), 4) for b in Mv] for a in Mv]
    GMp = pari.matrix(24, 24, [int(x) for row in GM for x in row])
    Gi = pari.matsolve(GMp, pari.matid(24))
    Ginv = [[Fraction(str(Gi[i, j])) for j in range(24)] for i in range(24)]
    den = 1
    for row in Ginv:
        for x in row:
            den = den * x.denominator // math.gcd(den, x.denominator)
    QD = np.array([[int(x * den) for x in row] for row in Ginv], dtype=np.int64)    # den * Gram of the dual basis

    def projM(x2):               # M-component of x (doubled coords), as M-basis coordinates
        g = [Fraction(lor2(a, x2), 4) for a in Mv]
        return [sum(Ginv[i][j] * g[j] for j in range(24)) for i in range(24)]

    def order(c):
        o = 1
        for t in c:
            o = o * t.denominator // math.gcd(o, t.denominator)
        return o
    gL = max((projM(x) for x in kernel([rho2], 25)), key=order)
    gN = max((projM(x) for x in kernel([z2], 25)), key=order)
    eL = sum(gL[i] * sum(GM[i][j] * gN[j] for j in range(24)) for i in range(24))   # g_L . g_N
    e = int((eL * h) % h)
    ok_gen = order(gL) == h and order(gN) == h and math.gcd(e, h) == 1 and int(pari.matdet(GMp)) == h * h

    def cls_dual(y):             # class of x = sum y_i e_i^* : (h x.g_N, h x.g_L) mod h
        return (int((sum(Fraction(int(y[i])) * gN[i] for i in range(24)) * h) % h),
                int((sum(Fraction(int(y[i])) * gL[i] for i in range(24)) * h) % h))

    def cls_M(c):                # class of an M (x) Q vector given in M-basis coordinates
        x_dot_gN = sum(c[i] * sum(GM[i][j] * gN[j] for j in range(24)) for i in range(24))
        x_dot_gL = sum(c[i] * sum(GM[i][j] * gL[j] for j in range(24)) for i in range(24))
        return (int((x_dot_gN * h) % h), int((x_dot_gL * h) % h))
    einv = pow(e, -1, h)
    qD = lambda g: Fraction((einv * g[0] * g[1]) % h, h)                  # q in these coordinates
    bD = lambda g, d: ((einv * (g[0] * d[1] + g[1] * d[0])) % h)          # h * b(g, d) mod h
    # consequence: the Leech-side line A and the N-side line B
    rp = kernel([rho2], 25)
    HL = {cls_M(projM(x)) for x in rp}
    HLsub = set()
    frontier = {(0, 0)}
    gens_L = list(HL)
    while frontier:
        nxt = set()
        for a_ in frontier:
            for b_ in gens_L:
                c_ = ((a_[0] + b_[0]) % h, (a_[1] + b_[1]) % h)
                if c_ not in HLsub:
                    HLsub.add(c_); nxt.add(c_)
        frontier = nxt
    a0 = cls_M(projM(r0))
    Aline = {((a0[0] + s[0]) % h, (a0[1] + s[1]) % h) for s in HLsub}
    # the N-side line B: projections of {x : x . z = -1} = y0 + (z^perp n II)
    zrow = [lor2(b, z2) // 4 for b in basis]
    gcoef, coefv = 0, [0] * len(zrow)
    for j, val in enumerate(zrow):
        if val == 0:
            continue
        if gcoef == 0:
            gcoef, coefv = val, [0] * len(zrow); coefv[j] = 1; continue
        a_0, b_0, a_1, b_1, r_0, r_1 = 1, 0, 0, 1, gcoef, val
        while r_1:
            qq = r_0 // r_1
            r_0, r_1 = r_1, r_0 - qq * r_1
            a_0, a_1 = a_1, a_0 - qq * a_1
            b_0, b_1 = b_1, b_0 - qq * b_1
        coefv = [a_0 * c for c in coefv]; coefv[j] += b_0; gcoef = r_0
    coefv = [-c * gcoef for c in coefv]
    y0 = [sum(coefv[j] * basis[j][t] for j in range(n + 1)) for t in range(n + 1)]
    assert lor2(y0, z2) // 4 == -1
    HN = {cls_M(projM(x)) for x in kernel([z2], 25)}
    HNsub, frontier = set(), {(0, 0)}
    while frontier:
        nxt = set()
        for a_ in frontier:
            for b_ in HN:
                c_ = ((a_[0] + b_[0]) % h, (a_[1] + b_[1]) % h)
                if c_ not in HNsub:
                    HNsub.add(c_); nxt.add(c_)
        frontier = nxt
    b0 = cls_M(projM(y0))
    Bline = {((b0[0] + s[0]) % h, (b0[1] + s[1]) % h) for s in HNsub}
    # A is the horizontal line {(*, a0[1])}, B the vertical line {(b0[0], *)}; choose the swap sigma_t(a, b) = (lam b, a/lam)
    lam = (b0[0] * pow(a0[1], -1, h)) % h
    laminv = pow(lam, -1, h)
    sigma = lambda g: ((lam * g[1]) % h, (laminv * g[0]) % h)
    lines_ok = (Bline == {sigma(g) for g in Aline} and len(HLsub) == h and len(HNsub) == h and not (HLsub & HNsub) - {(0, 0)}
                and all(qD(sigma(g)) == qD(g) for g in [(1, 0), (0, 1), (1, 1)]) and {sigma(g) for g in HLsub} == HNsub)
    # 1. enumerate M^* vectors of norm <= 2
    vecs = enum_dual(QD, 2.0 * den)
    # independent count of the nonzero vectors of M^* of norm <= 2 (PARI qfminim, counts +-v), against the enumerator
    n_pari = int(pari.qfminim(pari.matrix(24, 24, [int(x) for row in QD for x in row]), 2 * den, 0, 0)[0])
    n_enum = sum(1 for y in vecs if 0 < int(y @ QD @ y) <= 2 * den)
    count_ok = n_pari == n_enum
    table = {}
    qok = True
    for y in vecs:
        nrm = Fraction(int(y @ QD @ y), den)
        if nrm > 2:
            continue
        g = cls_dual(y)
        qok &= (nrm / 2 - qD(g)).denominator == 1
        table.setdefault(g, {}).setdefault(nrm, 0)
        table[g][nrm] += 1
    classes = [(a, b) for a in range(h) for b in range(h)]
    # 2. E below norm 2 vanishes
    below_ok = True
    v = {}
    for g in classes:
        tg, ts = table.get(g, {}), table.get(sigma(g), {})
        for nr in set(tg) | set(ts):
            if nr < 2 and tg.get(nr, 0) != ts.get(nr, 0):
                below_ok = False
        v[g] = tg.get(Fraction(2), 0) - ts.get(Fraction(2), 0)
    # 3. invariance of v: T (isotropic support) and S (exact in Z[zeta_h])
    T_ok = all(v[g] == 0 or qD(g) == 0 for g in classes)
    Phi = Poly(cyclotomic_poly(h, X), X, domain=ZZ)
    S_ok = True
    supp = [d for d in classes if v[d] != 0]
    for g in classes:
        w = [0] * h
        for d in supp:
            w[(-bD(g, d)) % h] += v[d]
        w[0] -= h * v[g]
        if any(w):
            P = Poly(sum(c * X ** j for j, c in enumerate(w)), X, domain=ZZ)
            if not rem(P, Phi).is_zero:
                S_ok = False
                break
    G = sum(v[g] for g in Aline)
    # decomposition v = sum_d c_d 1_{S_d} (S_d = maximal isotropic subgroups in these coordinates)
    divs = [d for d in range(1, h + 1) if h % d == 0]
    Ssets = {}
    for d in divs:
        Ssets[d] = {(a, b) for a in range(0, h, d) for b in range(0, h, h // d)}
    # solve least squares over the indicator vectors restricted to supp
    Mx = np.array([[1.0 if g in Ssets[d] else 0.0 for d in divs] for g in classes])
    vv = np.array([v[g] for g in classes], dtype=float)
    cd, *_ = np.linalg.lstsq(Mx, vv, rcond=None)
    recon_ok = np.allclose(Mx @ cd, vv, atol=1e-8)
    cd_int = {d: int(round(c)) for d, c in zip(divs, cd)}
    check(ok_gen and qok and count_ok and below_ok and T_ok and S_ok and len(Aline) == h and lines_ok and G == n0,
          f'{name:9s} h={h:2d}: M^* vectors |x|^2<=2: {sum(sum(t.values()) for t in table.values())} (PARI qfminim agrees: {count_ok}); E = 0 below norm 2; '
          f'v (norm-2 swap defect) is rho_M-invariant (T, S exact) => Theta_M - sigma^*Theta_M = Delta v (PROVED); '
          f'B = sigma(A); sum over the Leech line A of v = {G} = 24 + k = {n0} => theta(Lambda+c) - theta(N+rho/h) = {n0} Delta [{time.time() - t1:.0f}s]')
    # intermediate unimodular lattices L_d (from the maximal isotropic subgroups S_d): root counts and Coxeter numbers
    roots_d = {}
    for d in divs:
        Sd = Ssets[d]
        roots_d[d] = sum(table.get(g, {}).get(Fraction(2), 0) for g in Sd)
    inter = {d: (roots_d[d], roots_d[d] // 24 if roots_d[d] % 24 == 0 else None) for d in divs}
    print(f'       v = sum_d c_d 1_(S_d): {recon_ok}; c_d = {cd_int}', flush=True)
    print(f'       intermediate lattices L_d: (roots, Coxeter number) = {inter}', flush=True)
    INTER[name] = (h, inter, cd_int)

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
