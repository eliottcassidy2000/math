"""Independent audit of the session collatz-necklace-20260929 (note collatz_necklace_20260929_fair_splits_power_clocks_basins.md,
THM-4515 / THM-4516 / THM-4517, HYP-9165). Auditor: Claude (Fable 5.1) session, 2026-09-29. All code here is the auditor's own;
nothing from the session's scripts is imported or re-run (the session used gmpy2; this audit uses a residue sieve + exact
integer roots in pure Python/numpy, exact Z[zeta_j] arithmetic by reduction modulo Phi_j, and its own C sieves).

Run from the worktree root:   python3 04-computation/experiments/collatz_necklace_20260929_audit.py [output-path]
Default output: 05-knowledge/results/collatz_necklace_20260929_audit.out.  Needs python3 with numpy and sympy; gcc optional
(without it the basin sieves fall back to a tiny pure-python run).  Runtime about 6 minutes with gcc (2 GB RAM for the sieves).

Sections: A = THM-4515 (fair splits are circulant), B = THM-4516 (perfect-power clocks), C = THM-4517 (basins, sheet cap),
D = HYP-9165 and the status-label audit, E = summary of verdicts."""
import sys, os, time, platform, subprocess

# ================================================================================================
# SECTION A code
# ================================================================================================
"""Section A of the audit: THM-4515 (fair consecutive splits are circulant)."""
import sys, random, time
from math import gcd, comb, log, sqrt, pi, lgamma, exp
from fractions import Fraction
from itertools import combinations
from sympy import cyclotomic_poly, Poly, Symbol
_z = Symbol('z')

def phi_coeffs(j):
    """integer coefficients of Phi_j(z), low degree first"""
    P = Poly(cyclotomic_poly(j, _z), _z)
    return [int(c) for c in reversed(P.all_coeffs())]

class Zzeta:
    """exact arithmetic in Z[zeta_j] = Z[z]/(Phi_j(z)); elements = int lists of length deg Phi_j"""
    def __init__(self, j):
        self.j = j; self.phi = phi_coeffs(j); self.deg = len(self.phi) - 1
    def reduce(self, a):
        a = list(a)
        # divide by monic Phi_j
        for i in range(len(a) - 1, self.deg - 1, -1):
            c = a[i]
            if c:
                for t in range(self.deg + 1):
                    a[i - self.deg + t] -= c * self.phi[t]
        a = a[:self.deg] + [0] * max(0, self.deg - len(a))
        return a
    def zpow(self, e):
        e %= self.j
        a = [0] * (e + 1); a[e] = 1
        return self.reduce(a)
    def const(self, c):
        return self.reduce([c])
    def add(self, a, b):
        n = max(len(a), len(b))
        return self.reduce([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
    def scal(self, c, a):
        return self.reduce([c * t for t in a])
    def mul(self, a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for k, y in enumerate(b):
                    out[i + k] += x * y
        return self.reduce(out)
    def iszero(self, a):
        return all(t == 0 for t in self.reduce(a))

def carry(w, q):
    """c_w = sum_{w_j = 1} q^{(ones after j)} 2^j  (THM-4484 accumulation c <- q c + 2^j)"""
    c = 0
    for j, b in enumerate(w):
        if b:
            c = q * c + (1 << j)
    return c

def rot(w, r):
    return w[r:] + w[:r]

def T(y, q, d):
    return y // 2 if y % 2 == 0 else (q * y + d) // 2

def is_primitive(w):
    K = len(w)
    return all(rot(w, r) != w for r in range(1, K))

def fair_cuts(w, j):
    """cut positions r in [0, m) at which all j consecutive windows of length m have x = X/j ones"""
    K = len(w); m = K // j; X = sum(w); x = X // j
    ww = w + w
    return [r for r in range(m) if all(sum(ww[r + i * m: r + (i + 1) * m]) == x for i in range(j))]

def canon(w):
    return min(tuple(rot(list(w), r)) for r in range(len(w)))

def prim_necklaces(K, X):
    seen = set(); out = []
    for ones in combinations(range(K), X):
        w = [0] * K
        for p in ones: w[p] = 1
        c = canon(w)
        if c in seen: continue
        seen.add(c)
        if is_primitive(list(c)): out.append(list(c))
    return out

def sectionA(out):
    P = out
    P("=" * 78)
    P("SECTION A: THM-4515 fair consecutive splits are circulant (independent re-derivation)")
    P("=" * 78)
    rng = random.Random(4515)
    # ---------- A1: identities on random necklaces (rational cycles), several q, d, shapes ----------
    P("\nA1. Share equation, circulant DFT identities, clock factorisation, c_w = sum q^{(j-1-i)x} 2^{im} c_i,")
    P("    all checked EXACTLY (numerators over the common denominator D = 2^K - q^X) on random words.")
    n_tests = 0; n_fair = 0; n_bad = 0
    ZZ = {j: Zzeta(j) for j in (2, 3, 4, 5, 6, 8, 9, 12)}
    for trial in range(1500):
        q = rng.choice([3, 3, 3, 5, 7, 9, 11, 13, 15, 21, 27])
        j = rng.choice([2, 2, 3, 4, 6, 5, 8, 9, 12])
        m = rng.randint(2, 9 if j <= 4 else 5)
        x = rng.randint(1, m - 1)
        K, X = j * m, j * x
        # d odd, coprime to q; sometimes a multiple of the clock so that the word is integral
        D = (1 << K) - q ** X
        while True:
            d = rng.choice([1, -1, rng.randrange(-999, 1000, 2)])
            if d % 2 and gcd(abs(d), q) == 1: break
        if rng.random() < 0.3:
            d = D * rng.choice([1, -1, 3, -5]) if gcd(abs(D), q) == 1 else d
        w = [1] * X + [0] * (K - X); rng.shuffle(w)
        # the cycle (rational): y_i = d c_{rot^i w} / D
        Y = [d * carry(rot(w, i), q) for i in range(K)]  # numerators over D
        # check the step map on numerators: 2 Y_{i+1} = q^{w_i} Y_i + w_i d D
        for i in range(K):
            assert 2 * Y[(i + 1) % K] == q ** w[i] * Y[i] + w[i] * d * D, "one-step law failed"
        # general window equation for ANY window [r, r+m): 2^m Y_{r+m} = q^{x_r} Y_r + d D c_r
        for r in range(K):
            win = (w + w)[r:r + m]
            xr = sum(win); cr = carry(win, q)
            assert (1 << m) * Y[(r + m) % K] == q ** xr * Y[r] + d * D * cr, "window equation failed"
        n_tests += 1
        Z = ZZ[j]
        # clock factorisation prod_k (2^m - zeta^k q^x) = 2^K - q^X in Z[zeta_j]
        prod = Z.const(1)
        for k in range(j):
            prod = Z.mul(prod, Z.add(Z.const(1 << m), Z.scal(-q ** x, Z.zpow(k))))
        assert Z.iszero(Z.add(prod, Z.const(-D))), "clock factorisation failed"
        # prod_{e|j} Phi_e(2^m, q^x) (homogenised) = D
        pr = 1
        for e in range(1, j + 1):
            if j % e == 0:
                co = phi_coeffs(e); deg = len(co) - 1
                pr *= sum(co[t] * (2 ** m) ** t * (q ** x) ** (deg - t) for t in range(deg + 1))
        assert pr == D, "cyclotomic product failed"
        for r in fair_cuts(w, j):
            n_fair += 1
            ns = [Y[(r + i * m) % K] for i in range(j)]
            cs = [carry((w + w)[r + i * m: r + (i + 1) * m], q) for i in range(j)]
            # circulant share equations
            for i in range(j):
                assert (1 << m) * ns[(i + 1) % j] == q ** x * ns[i] + d * D * cs[i]
            # c_{rot^r w} = sum_i q^{(j-1-i)x} 2^{im} c_i  (j=2: c_w = q^x c_0 + 2^m c_1)
            cw = carry(rot(w, r), q)
            assert cw == sum(q ** ((j - 1 - i) * x) * (1 << (i * m)) * cs[i] for i in range(j)), "c_w decomposition failed"
            # DFT: N_k (2^m zeta^k - q^x) = d C_k with N_k = sum zeta^{-ik} n_i (numerators: times D on the right)
            for k in range(j):
                Nk = Z.const(0); Ck = Z.const(0)
                for i in range(j):
                    Nk = Z.add(Nk, Z.scal(ns[i], Z.zpow(-i * k)))
                    Ck = Z.add(Ck, Z.scal(cs[i], Z.zpow(-i * k)))
                lhs = Z.mul(Nk, Z.add(Z.scal(1 << m, Z.zpow(k)), Z.const(-q ** x)))
                rhs = Z.scal(d * D, Ck)
                if not Z.iszero(Z.add(lhs, Z.scal(-1, rhs))):
                    n_bad += 1
            # k = 0 and k = j/2 in plain integers
            assert ((1 << m) - q ** x) * sum(ns) == d * D * sum(cs)
            if j % 2 == 0:
                a = sum((-1) ** i * ns[i] for i in range(j)); ac = sum((-1) ** i * cs[i] for i in range(j))
                assert ((1 << m) + q ** x) * a == -d * D * ac, "alternating identity failed"
            if j == 2:
                c0, c1 = cs
                # residues: c_w = 2^m (c0 + c1) mod (2^m - q^x); c_w = 2^m (c1 - c0) mod (2^m + q^x)
                A = (1 << m) - q ** x; B = (1 << m) + q ** x
                assert (cw - (1 << m) * (c0 + c1)) % A == 0
                assert (cw - (1 << m) * (c1 - c0)) % B == 0
                assert gcd(abs(A), abs(B)) == 1, "factors not coprime"
                # CRT equivalence of the integrality criterion (for THIS word)
                lhs_int = (d * cw) % D == 0
                rhs_int = (d * (c0 + c1)) % A == 0 and (d * (c0 - c1)) % B == 0
                assert lhs_int == rhs_int, "CRT form failed"
    P(f"    random words tested: {n_tests}; fair splits found and checked: {n_fair}; DFT identity failures: {n_bad}")
    P("    -> share equation, clock factorisation (both forms), c_w decomposition, DFT identities,")
    P("       k=0 / k=j/2 integer forms, residue signs and the j=2 CRT equivalence: all HOLD (0 failures).")
    P("    Remark: eigenvalues of 2^m P - q^x I are 2^m zeta^k - q^x = zeta^k (2^m - zeta^{-k} q^x);")
    P("       note 1.2 calls '2^m - zeta^{-k} q^x' the eigenvalues (off by the unit zeta^k; harmless).")

    # ---------- A2: fair 2-split existence (discrete IVT) and the census K <= 18 ----------
    P("\nA2. Discrete IVT: h(r) = F(r+m) - F(r) - X/2 has h(r+m) = -h(r), unit steps -> zero. Exhaustive check.")
    Kmax = 18
    rows = {}; pairs = 0; ivt_fail = 0; crit_mismatch = 0
    for K in range(2, Kmax + 1):
        for X in range(1, K):
            g = gcd(K, X)
            if g == 1: continue
            neck = prim_necklaces(K, X)
            for j in range(2, g + 1):
                if g % j: continue
                m = K // j; x = X // j
                ok = 0
                for w in neck:
                    fc = fair_cuts(w, j)
                    has = bool(fc)
                    ok += has
                    pairs += 1
                    if j == 2:
                        # IVT argument itself: sign change of h on [r, r+m]
                        F = [0]
                        for b in w + w: F.append(F[-1] + b)
                        h = [F[r + m] - F[r] - X // 2 for r in range(K)]
                        assert all(h[(r + m) % K] == -h[r] for r in range(K))
                        assert all(abs(h[(r + 1) % K] - h[r]) <= 1 for r in range(K))
                        if not has: ivt_fail += 1
                    # criterion (a): G(t) = K F(t) - t X constant on the coset {r + i m}
                    F = [0]
                    for b in w + w: F.append(F[-1] + b)
                    crit_a = any(len({K * F[r + i * m] - (r + i * m) * X for i in range(j)}) == 1 for r in range(m))
                    if crit_a != has: crit_mismatch += 1
                    # criterion (b): X = j, runs of t cyclic gaps in ((t-1)m, (t+1)m)
                    if X == j:
                        ones = [i for i in range(K) if w[i]]
                        gaps = [(ones[(t + 1) % j] - ones[t]) % K for t in range(j)]
                        gg = gaps + gaps
                        crit_b = all((t - 1) * m < sum(gg[s:s + t]) < (t + 1) * m for t in range(1, j) for s in range(j))
                        if crit_b != has: crit_mismatch += 1
                        # criterion (c): j = 3: every gap <= 2m - 1
                        if j == 3:
                            crit_c = max(gaps) <= 2 * m - 1
                            if crit_c != has: crit_mismatch += 1
                rows[(K, X, j)] = (len(neck), ok)
    P(f"    (necklace, j) pairs with K <= {Kmax}: {pairs} (note: 23,201); IVT failures (j=2 without fair split): {ivt_fail};")
    P(f"    criterion (a)/(b)/(c) mismatches against brute force: {crit_mismatch}")
    census = {(12, 6, 3): (47, 75), (15, 6, 3): (204, 333), (18, 6, 3): (621, 1026), (18, 9, 3): (1592, 2700),
              (16, 4, 4): (42, 112), (16, 8, 4): (255, 800), (18, 6, 6): (107, 1026), (16, 8, 8): (30, 800), (18, 9, 9): (56, 2700)}
    P("    note's census rows vs mine:")
    allok = True
    for key, (a, b) in census.items():
        mine = rows[key]
        okk = mine == (b, a)
        allok &= okk
        P(f"      (K,X,j)={key}: note {a} of {b}; mine {mine[1]} of {mine[0]} {'OK' if okk else 'MISMATCH'}")
    P(f"    all census rows agree: {allok}")
    P("    three-bead law m(m-1) of 3 C(m,2) primitive necklaces of shape (3m,3):")
    for mm in range(2, 9):
        n_, ok_ = rows[(3 * mm, 3, 3)] if 3 * mm <= Kmax else (None, None)
        if n_ is None:
            neck = prim_necklaces(3 * mm, 3); n_ = len(neck); ok_ = sum(1 for w in neck if fair_cuts(w, 3))
        P(f"      m={mm}: {ok_} of {n_}; m(m-1)={mm * (mm - 1)}, 3C(m,2)={3 * comb(mm, 2)} {'OK' if (ok_, n_) == (mm * (mm - 1), 3 * comb(mm, 2)) else 'MISMATCH'}")
    # ---------- A3: expected number of fair cuts, exact second moment ----------
    P("\nA3. E[N_j] = m C(m,x)^j / C(jm,jx): brute-force mean over ALL words of five shapes (independent code).")
    for (K, X, j) in [(12, 6, 3), (12, 4, 4), (15, 6, 3), (16, 8, 4), (18, 9, 3)]:
        m = K // j; x = X // j; tot = 0; nw = 0; atl = 0
        for ones in combinations(range(K), X):
            w = [0] * K
            for p in ones: w[p] = 1
            n = len(fair_cuts(w, j)); tot += n; nw += 1; atl += (n > 0)
        ex = Fraction(m * comb(m, x) ** j, comb(K, X))
        P(f"      shape ({K},{X}) j={j}: mean {Fraction(tot, nw)} = {float(Fraction(tot, nw)):.6f}; formula {float(ex):.6f} {'OK' if Fraction(tot, nw) == ex else 'MISMATCH'}; P(N>=1)={atl / nw:.4f}")
    P("    Stirling constants: sqrt(3)/(2 pi rho(1-rho)) at rho=1/2: %.4f; at rho=log_3 2: %.4f; at rho=12/19: %.4f"
      % (sqrt(3) / (2 * pi * 0.25), sqrt(3) / (2 * pi * (log(2) / log(3)) * (1 - log(2) / log(3))), sqrt(3) / (2 * pi * (12 / 19) * (7 / 19))))
    P("    (note says '1.1847 at rho = log_3 2'; 1.1847 is the value at rho = 12/19, the exact-log_3 2 value is 1.1839)")

    def lc(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    import numpy as np
    P("    exact second moment E[N_3^2] = (m/C(3m,3x)) sum_delta sum_s C(delta,s)^3 C(m-delta,x-s)^3 (log-space, my code):")
    prev = None
    for m in [100, 300, 1000, 3000, 10000, 30000]:
        x = m // 2; LC = lc(3 * m, 3 * x)
        EN = m * exp(3 * lc(m, x) - LC)
        lg = np.array([lgamma(t + 1) for t in range(m + 1)])
        tot = 0.0
        for dlt in range(m):
            s_lo = max(0, x - (m - dlt)); s_hi = min(dlt, x)
            if s_lo > s_hi: continue
            s = np.arange(s_lo, s_hi + 1)
            la = lg[dlt] - lg[s] - lg[dlt - s]
            lb = lg[m - dlt] - lg[x - s] - lg[m - dlt - (x - s)]
            tot += np.exp(3 * la + 3 * lb - LC).sum()
        EN2 = m * tot
        slope = (EN2 - prev[1]) / (log(m) - log(prev[0])) if prev else float('nan')
        P(f"      m={m}: E[N]={EN:.4f} E[N^2]={EN2:.4f} E[N^2]/log m={EN2 / log(m):.4f} slope dE[N^2]/dlog m={slope:.3f} bound P>=E^2/E[N^2]={EN * EN / EN2:.4f} (x log m: {EN * EN / EN2 * log(m):.3f})")
        prev = (m, EN2)
    P("    -> the second-moment bound values 0.249/0.180/0.141 at m=10^2/10^3/10^4 are reproduced; E[N^2] grows with")
    P("       slope ~0.8 per unit log m (not 0.93: E[N^2]/log m is still falling), so 'E[N^2] ~ 0.93 log m' is loose.")
    P("       The O(log m) growth is not proved in the note (it is an exact computation to m = 3e4 + a local-CLT heuristic).")
    # ---------- A4: Monte Carlo (own seed, own code) ----------
    P("\nA4. Monte Carlo P(N_3 >= 1), rho = 1/2, own generator (numpy):")
    rng2 = np.random.default_rng(20260930)
    for m, Tn in [(1000, 4000), (10000, 4000), (100000, 3000), (300000, 1500)]:
        Kk = 3 * m; Xx = 3 * (m // 2); hits = 0; tot = 0
        for _ in range(Tn):
            w = np.zeros(Kk, dtype=np.int64); w[rng2.choice(Kk, Xx, replace=False)] = 1
            F = np.concatenate(([0], np.cumsum(np.concatenate((w, w)))))
            r = np.arange(m)
            ok = (F[r + m] - F[r] == Xx // 3) & (F[r + 2 * m] - F[r + m] == Xx // 3) & (F[r + 3 * m] - F[r + 2 * m] == Xx // 3)
            n = int(ok.sum()); hits += (n > 0); tot += n
        p = hits / Tn
        P(f"      m={m}: P(N>=1)={p:.3f} +- {sqrt(p * (1 - p) / Tn):.3f}; E[N]={tot / Tn:.3f}; P log m = {p * log(m):.2f}; E[N|N>=1]={tot / max(hits, 1):.2f}")
    # ---------- A5: fair Eliahou identity on real cycles ----------
    P("\nA5. Fair Eliahou identity n'/n = sqrt(P_0/P_1) on real integer cycles with a fair 2-split:")
    def cycle_from(y0, q, d, maxlen=10 ** 5):
        ys = [y0]; y = T(y0, q, d)
        while y != y0:
            ys.append(y); y = T(y, q, d)
            if len(ys) > maxlen: return None
        return ys
    for (q, d, y0) in [(3, 7, 5), (3, 11, 1), (3, 11, 13), (3, -17, 65), (3, -17, 73), (5, -9, 7), (3, 25, 7), (3, -29, 109)]:
        ys = cycle_from(y0, q, d); w = [y % 2 for y in ys]; K = len(ys); X = sum(w); m = K // 2
        if K % 2 or X % 2: continue
        nmin = min(ys)
        for r in fair_cuts(w, 2):
            n0 = ys[r]; n1 = ys[(r + m) % K]
            P0 = Fraction(1); P1 = Fraction(1)
            for i in range(m):
                y = ys[(r + i) % K]
                if y % 2: P0 *= (1 + Fraction(d, q * y))
                y = ys[(r + m + i) % K]
                if y % 2: P1 *= (1 + Fraction(d, q * y))
            ok = Fraction(n1, n0) ** 2 == P0 / P1
            if q * nmin > abs(d):
                bound = X * abs(d) / (2 * (q * nmin - abs(d)))
                btxt = f"|log(n'/n)|={abs(log(n1 / n0)):.4f} <= bound {bound:.4f}: {abs(log(n1 / n0)) <= bound}"
            else:
                btxt = f"|log(n'/n)|={abs(log(n1 / n0)):.4f}; bound N/A (needs q n_min > |d|, here {q}*{nmin} <= {abs(d)})"
            P(f"      {q}x{d:+d} cycle min {nmin} (K,X)=({K},{X}) cut r={r}: n={n0}, n'={n1}; (n'/n)^2 == P0/P1: {ok}; {btxt}")
    P("    Remark: the bound needs the implicit hypothesis q n_min > |d| (so that |d/(qy)| < 1); the identity itself needs nothing.")
    P("       The general bound X|d|/(2(q n_min - |d|)) gives X/(6 n_min) for 3x+1, while the note's 'relative precision")
    P("       X/(12 n_min)' uses that all factors 1 + 1/(3y) exceed 1 (one-sided); both are correct.")
    # ---------- A6: the Belaga-Mignotte 3x+17021 cycle ----------
    P("\nA6. The 3x+17021 cycle through 5: shape, fair antipodal cuts, alternating sum, fair 4-split.")
    ys = cycle_from(5, 3, 17021, maxlen=10 ** 5)
    w = [y % 2 for y in ys]; K = len(ys); X = sum(w)
    P(f"      length K={K}, odd steps X={X}, gcd={gcd(K, X)}, least element {min(ys)}, primitive word: {is_primitive(w)}")
    fc2 = fair_cuts(w, 2); fc4 = fair_cuts(w, 4)
    m = K // 2; x = X // 2
    r = fc2[0]; n0 = ys[r]; n1 = ys[(r + m) % K]
    c0 = carry((w + w)[r:r + m], 3); c1 = carry((w + w)[r + m:r + 2 * m], 3)
    alt = n0 - n1
    P(f"      fair 2-cuts: {len(fc2)} of {m} (note: 79 of 1070), first at r={fc2[0]} (cut elements {n0}, {n1}); fair 4-cuts: {len(fc4)} (note: none)")
    P(f"      alternating sum n0 - n1 = {alt}; -d(c0 - c1)/(2^m + 3^x) = {Fraction(-17021 * (c0 - c1), (1 << m) + 3 ** x)} -> {'OK' if alt == Fraction(-17021 * (c0 - c1), (1 << m) + 3 ** x) else 'MISMATCH'}")
    P(f"      sum: (2^m - 3^x)(n0 + n1) == d (c0 + c1): {((1 << m) - 3 ** x) * (n0 + n1) == 17021 * (c0 + c1)}")
    P(f"      2^m - 3^x has {len(str(abs((1 << m) - 3 ** x)))} digits, 2^m + 3^x has {len(str((1 << m) + 3 ** x))} digits (m={m}, x={x})")
    # ---------- A7: census of small 3x+d cycles with gcd(K,X) > 1 ----------
    P("\nA7. Primitive 3x+d cycles, |d| <= 100, d odd prime to 3, least element <= 100|d|, gcd(K,X) > 1 (own search):")
    hits = []
    for ad in range(1, 101, 2):
        if ad % 3 == 0: continue
        for d in (ad, -ad):
            seen = set()
            for y0 in range(1, 100 * ad + 1):
                if y0 in seen: continue
                y = y0; steps = 0
                while True:
                    y = T(y, 3, d); steps += 1
                    if y < y0 or steps > 5000: break
                    if y == y0:
                        ys = cycle_from(y0, 3, d); seen.update(ys)
                        w = [t % 2 for t in ys]; K = len(w); X = sum(w)
                        if gcd(K, X) > 1 and gcd(gcd(*ys), ad) == 1:
                            hits.append((d, y0, K, X, gcd(K, X)))
                        break
    P(f"      total {len(hits)} (note: 49); first ten: {hits[:10]}")
    ys = cycle_from(131, 3, 13); w = [y % 2 for y in ys]
    P(f"      3x+13 cycle through 131: shape ({len(w)},{sum(w)}), fair 3-cuts: {fair_cuts(w, 3)} (note: none)")
    # ---------- A8: the 30 necklaces of shape (11,7) ----------
    P("\nA8. Shape (11,7), clock 2^11 - 3^7 = -139: necklaces, run classes, residues of carries mod 139, rotation law.")
    K, X = 11, 7; Dc = 2 ** K - 3 ** X
    words = []
    for ones in combinations(range(K), X):
        w = [0] * K
        for p in ones: w[p] = 1
        words.append(w)
    necks = {}
    for w in words:
        necks.setdefault(canon(w), []).append(w)
    def nruns(w):
        # number of maximal runs of ones cyclically
        return sum(1 for i in range(K) if w[i] == 1 and w[i - 1] == 0)
    runclass = {}
    minima = []; leastabs = []
    integral = []
    for c, ws in necks.items():
        runclass[nruns(list(c))] = runclass.get(nruns(list(c)), 0) + 1
        cs = [carry(w, 3) for w in ws]
        minima.append(min(Fraction(cc, Dc) for cc in cs))          # true minimum (most negative element)
        leastabs.append(max(Fraction(cc, Dc) for cc in cs))        # element of least absolute value = -c_min/139
        if all(cc % 139 == 0 for cc in cs): integral.append((c, max(Fraction(cc, Dc) for cc in cs)))
    hist = {}
    for w in words:
        hist[carry(w, 3) % 139] = hist.get(carry(w, 3) % 139, 0) + 1
    P(f"      necklaces: {len(necks)} (all primitive: {all(is_primitive(list(c)) for c in necks)}); words: {len(words)}")
    P(f"      run classes (m-cycles): {dict(sorted(runclass.items()))} (note: 1/9/15/5)")
    P(f"      elements of least absolute value (-c_min/139, the note's 'minima') range [{float(min(leastabs)):.2f}, {float(max(leastabs)):.2f}] (note: [-27.10, -14.81]);")
    P(f"      true minima (-c_max/139) range [{float(min(minima)):.2f}, {float(max(minima)):.2f}]; integral necklaces: {[(''.join(map(str, c)), float(v)) for c, v in integral]}")
    P(f"      residue classes hit mod 139: {len(hist)} of 139 (note 130); class 0 multiplicity {hist.get(0)} (note 11); max nonzero multiplicity {max(v for k, v in hist.items() if k)} (note 5)")
    # rotation law c_{rot w} = (q^{w_0} c_w + w_0 D)/2
    okrot = all(2 * carry(rot(w, 1), 3) == 3 ** w[0] * carry(w, 3) + w[0] * Dc for w in words)
    P(f"      rotation law c_(rot w) = (3^(w_0) c_w + w_0 D)/2 on all 330 words: {okrot}; hence c -> 3^(w_0) c / 2 mod 139: consistent")
    P("\nSECTION A VERDICT: THM-4515 (1)-(4) HOLD; the census, the 17021-cycle facts and the (11,7) orderings HOLD;")
    P("  (5) E[N_j] exact formula and Stirling asymptotics HOLD; the 'PROVED lower bound of order 1/log m' is proved")
    P("  only for the m values computed (E[N_3^2] = O(log m) for all m is not proved in the note).")

# ================================================================================================
# SECTION B code
# ================================================================================================
"""Section B of the audit: THM-4516 (perfect-power clocks are Fermat-Catalan identities)."""
import sys, time
from math import gcd, log, log2
from itertools import combinations
import numpy as np
from sympy import primerange, isprime, primitive_root

def iroot(n, r):
    """floor of the r-th root of n >= 0 (exact integer Newton, converging from above)"""
    if n < 2: return n
    x = 1 << ((n.bit_length() + r - 1) // r)
    while True:
        y = ((r - 1) * x + n // pow(x, r - 1)) // r
        if y >= x: return x
        x = y

def perfect_power_exact(n):
    """(m, e) with n = m^e, e maximal (e = 1 if n is not a perfect power); exact, all prime exponents tested"""
    if n < 4: return (n, 1)
    e = 1; M = n
    Rmax = int(log(n) / log(3)) + 2
    for r in primerange(2, Rmax + 1):
        while True:
            m = iroot(M, r)
            if m ** r == M and m > 1:
                M = m; e *= r
            else:
                break
    return (M, e)

class PowerSieve:
    """necessary-condition sieve for |2^K - q^X| = m^r over all primes r <= Rmax, then exact confirmation"""
    def __init__(self, Rmax, log):
        self.log = log
        self.rs = list(primerange(2, Rmax + 1))
        sched = lambda r: 24 if r == 2 else 16 if r == 3 else 12 if r == 5 else 10 if r == 7 else 6 if r <= 30 else 4 if r <= 100 else 3
        self.Lr = {}
        self.isres = {}
        allL = set()
        for r in self.rs:
            need = sched(r); ls = []
            l = r + 1
            while len(ls) < need:
                if isprime(l): ls.append(l)
                l += r
            self.Lr[r] = ls
            for l in ls:
                allL.add(l)
                g = primitive_root(l); h = pow(g, r, l)
                tab = np.zeros(l, dtype=bool); tab[0] = True
                v = 1
                for _ in range((l - 1) // r):
                    tab[v] = True; v = v * h % l
                self.isres[(l, r)] = tab
        self.allL = sorted(allL)
        self.small = list(primerange(2, 200))
        log(f"    sieve: {len(self.rs)} prime exponents r <= {Rmax}, {len(self.allL)} auxiliary primes (max {max(self.allL)}), small-factor filter primes < 200")

    def run(self, qs, Kmax, Xmin=2, qblock=100):
        rs = self.rs
        Ks, Xs = [], []
        for K in range(Xmin, Kmax + 1):
            for X in range(Xmin, K + 1):
                Ks.append(K); Xs.append(X)
        Kp = np.array(Ks, dtype=np.int64); Xp = np.array(Xs, dtype=np.int64); Pn = len(Kp)
        # per-l tables of 2^K mod l (and mod l^2 for the small filter)
        pow2 = {}
        for l in set(self.allL) | {p * p for p in self.small}:
            t = np.ones(Kmax + 1, dtype=np.int64)
            for K in range(1, Kmax + 1): t[K] = t[K - 1] * 2 % l
            pow2[l] = t
        found = []; units = []
        ln2 = log(2.0); l3 = log(3.0)
        qs = list(qs)
        for b0 in range(0, len(qs), qblock):
            qb = np.array(qs[b0:b0 + qblock], dtype=np.int64); nq = len(qb)
            lnq = np.log(qb.astype(float))
            diff = Kp[None, :] * ln2 - Xp[None, :] * lnq[:, None]          # log(2^K / q^X)
            sgn = np.where(diff > 0, 1, -1).astype(np.int64)
            near = np.abs(diff) < 1e-6
            for (qi, pi) in zip(*np.nonzero(near)):
                q = int(qb[qi]); K = int(Kp[pi]); X = int(Xp[pi])
                sgn[qi, pi] = 1 if (1 << K) > q ** X else -1
            alive = np.ones((nq, Pn), dtype=bool)
            # ---- small-factor filter: a prime l dividing N exactly once kills every exponent
            for p in self.small:
                L2 = p * p
                tq = np.ones((nq, Kmax + 1), dtype=np.int64); qm = qb % L2
                for X in range(1, Kmax + 1): tq[:, X] = tq[:, X - 1] * qm % L2
                D2 = (pow2[L2][Kp][None, :] - tq[:, Xp]) % L2
                alive &= ~((D2 % p == 0) & (D2 != 0))
            qi_s, pi_s = np.nonzero(alive)
            K_s = Kp[pi_s]; X_s = Xp[pi_s]; sg_s = sgn[qi_s, pi_s]
            logN_upper = np.maximum(K_s * ln2, X_s * lnq[qi_s])              # log N <= this
            self.log(f"      q-block {qs[b0]}..{qs[min(b0 + qblock, len(qs)) - 1]}: {nq * Pn} pairs, {len(qi_s)} survive the small-factor filter")
            # per-l tables q^X mod l for this block
            tql = {}
            for l in self.allL:
                t = np.ones((nq, Kmax + 1), dtype=np.int64); qm = qb % l
                for X in range(1, Kmax + 1): t[:, X] = t[:, X - 1] * qm % l
                tql[l] = t
            for r in rs:
                cand = np.nonzero(logN_upper >= r * l3 - 1e-9)[0]      # m >= 3 needs N >= 3^r (N = 1 handled below)
                cand = np.concatenate((cand, np.nonzero(logN_upper < r * l3 - 1e-9)[0][:0]))
                for l in self.Lr[r]:
                    if len(cand) == 0: break
                    D = (pow2[l][K_s[cand]] - tql[l][qi_s[cand], X_s[cand]]) % l
                    if r == 2: D = (sg_s[cand] * D) % l
                    cand = cand[self.isres[(l, r)][D]]
                for idx in cand:
                    q = int(qb[qi_s[idx]]); K = int(K_s[idx]); X = int(X_s[idx])
                    g = (1 << K) - q ** X; N = abs(g)
                    if N == 1:
                        units.append((q, K, X, g)); continue
                    m = iroot(N, r)
                    if m ** r == N:
                        found.append((q, K, X, 1 if g > 0 else -1, N, r))
            # N = 1 candidates that fail the size test for every r: check the near-ties exactly
            for (qi, pi) in zip(*np.nonzero(near)):
                q = int(qb[qi]); K = int(Kp[pi]); X = int(Xp[pi]); g = (1 << K) - q ** X
                if abs(g) == 1 and (q, K, X, g) not in units: units.append((q, K, X, g))
        # consolidate: one entry per (q,K,X) with the maximal exponent
        res = {}
        for (q, K, X, s, N, r) in found:
            if (q, K, X) not in res:
                M, e = perfect_power_exact(N)
                res[(q, K, X)] = (s, M, e)
        return res, sorted(set(units))

def brute_census(qs, Kmax, Xmin=2):
    res = {}; units = []
    for q in qs:
        for K in range(Xmin, Kmax + 1):
            for X in range(Xmin, K + 1):
                g = (1 << K) - q ** X; N = abs(g)
                if N == 1: units.append((q, K, X, g)); continue
                M, e = perfect_power_exact(N)
                if e > 1: res[(q, K, X)] = (1 if g > 0 else -1, M, e)
    return res, sorted(set(units))

def classify(res):
    fam = []; other = []
    for (q, K, X), (s, M, e) in sorted(res.items()):
        if X == 2 and q == (1 << (K - 2)) + 1: fam.append((q, K, X, s, M, e))
        else: other.append((q, K, X, s, M, e))
    return fam, other

def cycle_from(y0, q, d, maxlen=10 ** 5):
    ys = [y0]; y = T(y0, q, d)
    while y != y0:
        ys.append(y); y = T(y, q, d)
        if len(ys) > maxlen: return None
    return ys

def all_necklaces(K, X):
    s = set()
    for ones in combinations(range(K), X):
        w = [0] * K
        for p in ones: w[p] = 1
        s.add(canon(w))
    return sorted(s)

def sectionB(out):
    P = out
    P("=" * 78)
    P("SECTION B: THM-4516 perfect-power clocks (independent census: residue sieve + exact roots, no gmpy2)")
    P("=" * 78)
    # ---------- B1: dictionary ----------
    P("\nB1. Dictionary. 2^K - q^X = +-m^r with q = p^s odd, K >= 1: m is odd (2^K - q^X is odd), p does not divide m")
    P("    (else p | 2^K), so gcd(m, 2p) = 1 and (2, p, m) with exponents (K, sX, r) is a primitive solution of")
    P("    x^a + y^b = z^c; Darmon-Granville finiteness applies when 1/K + 1/(sX) + 1/r < 1. By THM-4484(1) (clock | d)")
    P("    the shape (K,X) is free for T_{q, +-m^r} since the clock equals +-d. CONSISTENT (the converse is NOT claimed")
    P("    correctly in section 2.5 / THM-4516(5); see B5).")
    # ---------- B2: census ----------
    P("\nB2. Census of |2^K - q^X| in {1} u {m^r : r >= 2}, X >= 2 (own sieve; every prime exponent r <= log_3 N tested).")
    t0 = time.time()
    # validation of the sieve against brute force on a small window
    qs_small = list(range(3, 42, 2)); Kv = 100
    sv = PowerSieve(int(Kv * log(41) / log(3)) + 1, P)
    res_s, units_s = sv.run(qs_small, Kv, qblock=20)
    res_b, units_b = brute_census(qs_small, Kv)
    P(f"    validation (odd q <= 41, 2 <= X <= K <= {Kv}): sieve {len(res_s)} clocks / brute force {len(res_b)}; identical: {res_s == res_b and units_s == units_b}")
    P(f"      found: {[(q, K, X, s * M ** e, M, e) for (q, K, X), (s, M, e) in sorted(res_b.items())]}; units: {units_b}")
    # (i) odd q <= 201, 2 <= X <= K <= 400
    Kc = 400; qs = list(range(3, 202, 2))
    sv = PowerSieve(int(Kc * log(201) / log(3)) + 1, P)
    res1, units1 = sv.run(qs, Kc, qblock=50)
    fam, other = classify(res1)
    P(f"    (i) odd q <= 201, 2 <= X <= K <= {Kc}: {len(res1)} perfect-power clocks + {len(units1)} unit clocks [{time.time() - t0:.0f} s]")
    P(f"        units |2^K - q^X| = 1 with X >= 2: {units1} (Gersonides/Catalan: only (3,3,2) expected)")
    P(f"        Pythagorean family q = 2^(K-2)+1, X = 2: {[(q, K) for (q, K, X, s, M, e) in fam]}")
    P(f"        all others (q, K, X, sign, m, r): {[(q, K, X, s, M, e) for (q, K, X, s, M, e) in other]}")
    exp_other = {(3, 5, 4): (-1, 7, 2), (7, 9, 3): (1, 13, 2), (13, 9, 2): (1, 7, 3), (71, 7, 2): (-1, 17, 3)}
    got_other = {(q, K, X): (s, M, e) for (q, K, X, s, M, e) in other}
    P(f"        matches the note's list exactly (four Fermat-Catalan readings, family q=5,9,17,33,65,129): {got_other == exp_other and [q for (q, K, X, s, M, e) in fam] == [5, 9, 17, 33, 65, 129]}")
    # (ii) q = 3, K <= 3000
    t1 = time.time()
    Kc = 3000
    sv = PowerSieve(Kc + 1, P)
    res2, units2 = sv.run([3], Kc, qblock=1)
    P(f"    (ii) q = 3, 2 <= X <= K <= 3000: clocks {[(q, K, X, s * M ** e if M ** e < 10 ** 12 else ('sign', s, 'm^r with', len(str(M)), 'digits'), M, e) for (q, K, X), (s, M, e) in sorted(res2.items())]}; units {units2} [{time.time() - t1:.0f} s]")
    P(f"        nothing beyond (3,3,2) [-1] and (3,5,4) [-7^2]: {set(res2) == {(3, 5, 4)} and units2 == [(3, 3, 2, -1)]}")
    # (iii) odd q <= 10^4, K <= 60
    t1 = time.time()
    Kc = 60; qs = list(range(3, 10001, 2))
    sv = PowerSieve(int(Kc * log(10000) / log(3)) + 1, P)
    res3, units3 = sv.run(qs, Kc, qblock=1000)
    fam3, other3 = classify(res3)
    P(f"    (iii) odd q <= 10^4, 2 <= X <= K <= 60: {len(res3)} clocks = {len(fam3)} family members (K = {min(K for (q, K, X, s, M, e) in fam3)}..{max(K for (q, K, X, s, M, e) in fam3)}) + others {[(q, K, X, s, M, e) for (q, K, X, s, M, e) in other3]}; units {units3} [{time.time() - t1:.0f} s]")
    P(f"        note claims 17 = 13 family (K=3..15) + the four readings; with the unit clock (3,3,2) counted as the K = 3 family member")
    P(f"        (2^3 + 1^2 = 3^2), mine: {len(res3) + len(units3)} = {len(fam3) + len(units3)} family + {len(other3)} others: {len(res3) + len(units3) == 17 and len(fam3) + len(units3) == 13 and {(q, K, X): (s, M, e) for (q, K, X, s, M, e) in other3} == exp_other}")
    # (iv) odd q <= 10^5, K <= 40
    t1 = time.time()
    Kc = 40; qs = list(range(3, 100001, 2))
    sv = PowerSieve(int(Kc * log(100000) / log(3)) + 1, P)
    res4, units4 = sv.run(qs, Kc, qblock=2500)
    fam4, other4 = classify(res4)
    P(f"    (iv) odd q <= 10^5, 2 <= X <= K <= 40: {len(res4)} clocks = {len(fam4)} family members (K = {min(K for (q, K, X, s, M, e) in fam4)}..{max(K for (q, K, X, s, M, e) in fam4)}) + others {[(q, K, X, s, M, e) for (q, K, X, s, M, e) in other4]}; units {units4} [{time.time() - t1:.0f} s]")
    P(f"        note claims 20 = 16 family (K=3..18) + the four readings: {len(res4) + len(units4) == 20 and len(fam4) + len(units4) == 16}")
    P("    Method remark: the session's power_clocks.py identifies the exponent with gmpy2.iroot(n, r) for r <= 64 only and")
    P("      DROPS a perfect power whose least prime exponent exceeds 64 (power_data returns None); the wide script keeps it")
    P("      with r = None. My sieve tests every prime r <= log_3 N, so the census statement is now covered for all r.")
    # ---------- B3: free cycles by direct orbit computation ----------
    P("\nB3. Free cycles of the perfect-power clocks, every necklace of the shape iterated directly:")
    for (q, d, K, X, idn) in [(3, -49, 5, 4, "2^5+7^2=3^4"), (9, -49, 5, 2, "2^5+7^2=9^2"), (7, 169, 9, 3, "7^3+13^2=2^9"),
                              (13, 343, 9, 2, "7^3+13^2=2^9"), (71, -4913, 7, 2, "2^7+17^3=71^2"), (5, -9, 4, 2, "2^4+3^2=5^2"),
                              (17, -225, 6, 2, "2^6+15^2=17^2"), (3, 125, 7, 1, "2^7=3+5^3")]:
        D = (1 << K) - q ** X
        assert d % D == 0
        rows = []
        for wn in all_necklaces(K, X):
            y0 = d * carry(list(wn), q) // D
            ys = cycle_from(y0, q, d); assert ys is not None
            assert [y % 2 for y in ys] == list(wn)[:len(ys)]
            g = gcd(gcd(*ys), abs(d))
            rows.append((''.join(map(str, wn)), min(ys), len(ys), 'prim' if (len(ys) == K and g == 1) else f'period {len(ys)}, gcd {g}'))
        P(f"      {q}x{d:+d} [{idn}] shape ({K},{X}) clock {D}: {rows}")
    P("      3x-49 cycle from 65: " + str(cycle_from(65, 3, -49)) + "; 3x+125 cycle from 1: " + str(cycle_from(1, 3, 125)))
    P("      13x+343 fourth necklace: least element 21 = 7 x 3, i.e. 7 x (cycle of 13x+49 through 3): " + str(cycle_from(3, 13, 49)) + " x 7 = " + str(cycle_from(21, 13, 343)))
    P("      -> note 2.2's phrase '21 = 7 x (3x+49 ...) scaled' should read '7 x (13x+49 cycle through 3)'.")
    # ---------- B4: Eisenstein ----------
    P("\nB4. Eisenstein factorisation of 2^9 - 7^3:")
    Z3 = Zzeta(3)
    lhs = Z3.add(Z3.const(8), Z3.scal(-7, Z3.zpow(1)))
    sq = Z3.mul(Z3.add(Z3.const(3), Z3.scal(-1, Z3.zpow(1))), Z3.add(Z3.const(3), Z3.scal(-1, Z3.zpow(1))))
    lhs2 = Z3.add(Z3.const(8), Z3.scal(-7, Z3.zpow(2)))
    sq2 = Z3.mul(Z3.add(Z3.const(3), Z3.scal(-1, Z3.zpow(2))), Z3.add(Z3.const(3), Z3.scal(-1, Z3.zpow(2))))
    P(f"      2^9 - 7^3 = {2 ** 9 - 7 ** 3} = 13^2; (2^3 - 7) Phi_3(8,7) = {(8 - 7) * (64 + 56 + 49)}; 8 - 7 zeta == (3 - zeta)^2 mod Phi_3: {Z3.iszero(Z3.add(lhs, Z3.scal(-1, sq)))};"
      f" 8 - 7 zeta^2 == (3 - zeta^2)^2: {Z3.iszero(Z3.add(lhs2, Z3.scal(-1, sq2)))}; N(3 - zeta) = 9 + 3 + 1 = {9 + 3 + 1}")
    P("      13 = 1 mod 3 splits in Z[zeta_3]; (3 - zeta)(3 - zeta^2) = 13, so the twisted clocks are squares of the two primes above 13: HOLDS.")
    # ---------- B5: the Beal-shadow overreach ----------
    P("\nB5. Beal shadow. THM-4516(5)/note 2.5 conclude: 'Beal implies no map py +- m^r (r >= 3) has a free mixed shape (K,X)")
    P("    with K >= 3, sX >= 3, X >= 2'. Freeness needs only (2^K - p^X) | m^r (THM-4484), NOT 2^K - p^X = +-m^r. Counterexamples:")
    for (q, d, K, X) in [(3, 125, 5, 3), (3, 47 ** 3, 7, 4)]:
        D = (1 << K) - q ** X
        rows = []
        for wn in all_necklaces(K, X):
            y0 = d * carry(list(wn), q) // D
            ys = cycle_from(y0, q, d)
            rows.append((''.join(map(str, wn)), min(ys), len(ys), all(isinstance(y, int) for y in ys)))
        M, e = perfect_power_exact(d)
        P(f"      {q}y+{d} = {q}y+{M}^{e}: clock 2^{K} - {q}^{X} = {D} divides d, so shape ({K},{X}) is free (K={K}, sX={X}, r={e}); necklaces -> integer cycles {rows}; the clock {D} is not +-m^r, no Beal solution involved")
    P("    -> the sentence is FALSE as stated; the true consequence of Beal is only: no perfect-power CLOCK 2^K - p^{sX} = +-m^r with K, sX, r >= 3.")
    # ---------- B6: the (4,3,17) clock shadow ----------
    P("\nB6. The (4,3,17) 'clock shadow is empty' sentence. A clock with exponent multiset {4,3,17} means 2^K = q^X +- m^r")
    P("    with {K, sX, r} = {4, 3, 17}. The note's reason ('no term can be a pure power of two while the other two are powers")
    P("    of one odd base') is not an argument: a coprime solution of x^4 + y^3 = z^17 has exactly one even variable, and")
    P("    the clock shadow is the sub-case where that variable is a power of 2 (e.g. x = 2: 16 = z^17 - y^3). Finite checks:")
    # sum cases (both odd terms positive) are finite
    sols = []
    for (K, a, b) in [(4, 3, 17), (4, 17, 3), (3, 4, 17), (3, 17, 4), (17, 3, 4), (17, 4, 3)]:
        N = 1 << K; m = 1
        while m ** b < N:
            rest = N - m ** b
            if rest > 0:
                r_ = iroot(rest, a)
                if r_ ** a == rest and rest % 2 == 1 and m % 2 == 1: sols.append((K, a, b, r_, m))
            m += 2
    P(f"      2^K = u^a + v^b with u, v odd, {{K,a,b}} = {{4,3,17}}: solutions {sols} (finite search complete)")
    # difference cases: bounded search only
    cnt = 0; sols = []
    B = 3000
    for u in range(1, B, 2):
        for (K, a, b) in [(4, 3, 17), (4, 17, 3), (3, 4, 17), (3, 17, 4), (17, 3, 4), (17, 4, 3)]:
            N = 1 << K
            for val in (u ** a + N, u ** a - N):
                if val > 1:
                    r_ = iroot(val, b)
                    if r_ ** b == val and r_ % 2 == 1: sols.append((K, a, b, u, r_, val))
    P(f"      2^K = u^a - v^b or v^b - u^a with u odd < {B}: solutions {sols} (bounded search; NOT a proof of emptiness)")
    P("    -> 'Its clock shadow is empty in any case' is UNSUPPORTED; emptiness of the shadow is exactly as open as the (4,3,17)")
    P("       equation restricted to solutions with a power-of-two variable (plus the K = 4a readings for x = 2^a).")
    P("\nSECTION B VERDICT: dictionary, census (now covered for all exponents r), free cycles, Eisenstein square: HOLD;")
    P("  the Beal-shadow consequence in THM-4516(5)/title and the 'clock shadow is empty' sentence: FAIL as stated.")

# ================================================================================================
# SECTION C code
# ================================================================================================
"""Section C of the audit: THM-4517 (no root-uniform positive density; sheet cap; basin densities)."""
import sys, os, time, subprocess, tempfile, shutil
from math import gcd
from fractions import Fraction

C_MINUS = r'''
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
/* audit sieve, 3x-1 on [1,N]: T(n)=n/2 (even), (3n-1)/2 (odd). For each n record the first orbit value
   whose odd part lies in {1} u {5,7} u {17,25,37,55,41,61,91} as 1 + idx*64 + v_2; memoise on the first
   value below n (the orbit of n reaches it before any ray element, so the first ray element is inherited). */
static const uint32_t ODD[10]={1,5,7,17,25,37,55,41,61,91};
static const int CYC[10]={1,2,2,3,3,3,3,3,3,3};
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10);
  uint16_t*code=(uint16_t*)calloc(N+1,sizeof(uint16_t)); if(!code){fprintf(stderr,"alloc\n");return 1;}
  uint8_t lut[92]; memset(lut,0xff,92); for(int i=0;i<10;i++) lut[ODD[i]]=(uint8_t)i;
  u128 maxv=0;
  for(uint64_t n=1;n<=N;n++){
    u128 v=n;
    for(;;){
      uint64_t k=0; u128 o=v; while(!(o&1)){o>>=1;k++;}
      if(o<=91 && lut[(int)o]!=0xff){ code[n]=(uint16_t)(1+lut[(int)o]*64+(k>63?63:k)); break; }
      if(v<n){ code[n]=code[(uint64_t)v]; break; }
      v=(v&1)?(3*v-1)/2:(v>>1); if(v>maxv) maxv=v;
    }
  }
  int S=0; while((1ULL<<(S+1))<=N) S++;
  uint64_t cum[4]={0,0,0,0}; uint64_t unset=0;
  static uint64_t spec[641]; memset(spec,0,sizeof spec);
  static uint64_t spec29[641]; memset(spec29,0,sizeof spec29);
  uint64_t cum29[4]={0,0,0,0};
  printf("N=%llu maxorbit_bits=%d\n",(unsigned long long)N,(int)(maxv>>64?64+63-__builtin_clzll((uint64_t)(maxv>>64)):63-__builtin_clzll((uint64_t)maxv)));
  printf("dyadic rows: s  count1 count2 count3  dens1 dens2 dens3\n");
  for(int s=0;s<=S;s++){
    uint64_t lo=1ULL<<s, hi=(1ULL<<(s+1))-1; if(hi>N) hi=N; uint64_t b[4]={0,0,0,0};
    for(uint64_t n=lo;n<=hi;n++){ uint16_t c=code[n]; if(!c){unset++;continue;} int id=(c-1)/64; b[CYC[id]]++; spec[c]++; if(n<=(1ULL<<29)) spec29[c]++; }
    uint64_t w=hi-lo+1;
    printf("[2^%d,2^%d) %llu %llu %llu  %.7f %.7f %.7f\n",s,s+1,(unsigned long long)b[1],(unsigned long long)b[2],(unsigned long long)b[3],(double)b[1]/w,(double)b[2]/w,(double)b[3]/w);
    for(int i=1;i<4;i++){ cum[i]+=b[i]; if(hi<=(1ULL<<29)) cum29[i]+=b[i]; }
  }
  printf("unset=%llu\n",(unsigned long long)unset);
  printf("cumulative [1,2^29]: %llu %llu %llu  %.7f %.7f %.7f\n",(unsigned long long)cum29[1],(unsigned long long)cum29[2],(unsigned long long)cum29[3],(double)cum29[1]/(1ULL<<29),(double)cum29[2]/(1ULL<<29),(double)cum29[3]/(1ULL<<29));
  printf("cumulative [1,N]: %llu %llu %llu  %.8f %.8f %.8f\n",(unsigned long long)cum[1],(unsigned long long)cum[2],(unsigned long long)cum[3],(double)cum[1]/N,(double)cum[2]/N,(double)cum[3]/N);
  printf("ray-entry spectrum on [1,2^29] (odd,k,density>1e-4):\n");
  for(int c=1;c<641;c++) if(spec29[c] && (double)spec29[c]/(1ULL<<29)>1e-4) printf("  odd=%u k=%d dens=%.6f\n",ODD[(c-1)/64],(c-1)%64,(double)spec29[c]/(1ULL<<29));
  printf("ray-entry spectrum on [1,N] (odd,k,density>1e-4):\n");
  for(int c=1;c<641;c++) if(spec[c] && (double)spec[c]/N>1e-4) printf("  odd=%u k=%d dens=%.6f\n",ODD[(c-1)/64],(c-1)%64,(double)spec[c]/N);
  free(code); return 0;
}
'''

C_PLUS = r'''
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
/* audit sieve, 3x+1 (T-form) on [1,N]: entry index i such that the first power of two on the orbit is 2^(2i-1)
   (entered from a_i = (4^i-1)/3). Powers of two themselves: code 255. Even-exponent first hit (impossible for
   non-powers of two): code 254. Memoise on the first value below n. */
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10);
  uint8_t*code=(uint8_t*)calloc(N+1,1); if(!code){fprintf(stderr,"alloc\n");return 1;}
  for(uint64_t n=1;n<=N;n++){
    u128 v=n;
    for(;;){
      if((v&(v-1))==0){ int j=0; u128 t=v; while(t>1){t>>=1;j++;} if(v==(u128)n) code[n]=255; else if(j&1) code[n]=(uint8_t)((j+1)/2); else code[n]=254; break; }
      if(v<n){ code[n]=code[(uint64_t)v]; break; }
      v=(v&1)?(3*v+1)/2:(v>>1);
    }
  }
  int S=0; while((1ULL<<(S+1))<=N) S++;
  static uint64_t cnt[256]; static uint64_t tot[256]; memset(tot,0,sizeof tot);
  printf("N=%llu\n",(unsigned long long)N);
  printf("dyadic rows: density of entry index i = 2..14 (first trunk element 2^(2i-1), from a_i = (4^i-1)/3)\n");
  for(int s=20;s<=S;s++){
    uint64_t lo=1ULL<<s, hi=(1ULL<<(s+1))-1; if(hi>N) hi=N; memset(cnt,0,sizeof cnt);
    for(uint64_t n=lo;n<=hi;n++) cnt[code[n]]++;
    uint64_t w=hi-lo+1; printf("[2^%d,2^%d):",s,s+1);
    for(int i=2;i<=14;i++) printf(" %.5f",(double)cnt[i]/w);
    printf("  pow2=%llu anomalies=%llu unset=%llu\n",(unsigned long long)cnt[255],(unsigned long long)cnt[254],(unsigned long long)cnt[0]);
  }
  for(uint64_t n=1;n<=N;n++) tot[code[n]]++;
  printf("cumulative exact counts on [1,N] by entry index i (a_i, count, density):\n");
  for(int i=2;i<=40;i++) if(tot[i]) { unsigned long long a=((1ULL<<(2*i))-1)/3; printf("  i=%d a_i=%llu count=%llu dens=%.9f\n",i,a,(unsigned long long)tot[i],(double)tot[i]/N); }
  printf("  powers of two: %llu; anomalies (even-exponent first hit): %llu; unset: %llu\n",(unsigned long long)tot[255],(unsigned long long)tot[254],(unsigned long long)tot[0]);
  free(code); return 0;
}
'''

def have_gcc():
    try:
        subprocess.run(['gcc', '--version'], capture_output=True, check=True); return True
    except Exception:
        return False

def compile_and_run(src, name, arg, log):
    d = tempfile.mkdtemp(prefix='collatz_audit_')
    cpath = os.path.join(d, name + '.c'); exe = os.path.join(d, name + ('.exe' if os.name == 'nt' else ''))
    with open(cpath, 'w') as f: f.write(src)
    r = subprocess.run(['gcc', '-O3', '-o', exe, cpath], capture_output=True, text=True)
    if r.returncode != 0:
        log("      gcc failed: " + r.stderr[:500]); return None
    t0 = time.time()
    r = subprocess.run([exe, str(arg)], capture_output=True, text=True)
    log(f"      [{name} {arg}: {time.time() - t0:.0f} s]")
    shutil.rmtree(d, ignore_errors=True)
    return r.stdout

def py_minus(N):
    """pure-python fallback (small N): basin index 1/2/3 of 3x-1 for n <= N"""
    odd = {1: 1, 5: 2, 7: 2, 17: 3, 25: 3, 37: 3, 55: 3, 41: 3, 61: 3, 91: 3}
    code = [0] * (N + 1)
    for n in range(1, N + 1):
        v = n
        while True:
            o = v
            while o % 2 == 0: o //= 2
            if o in odd: code[n] = odd[o]; break
            if v < n: code[n] = code[v]; break
            v = (3 * v - 1) // 2 if v % 2 else v // 2
    return code

def sectionC(out):
    P = out
    P("=" * 78)
    P("SECTION C: THM-4517 no root-uniform positive density; sheet cap; FINITE-EXACT basin densities")
    P("=" * 78)
    # ---------- C1: the proposition ----------
    P("\nC1. Proposition 3.2 re-derived.")
    def Tp(n): return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
    a = lambda i: (4 ** i - 1) // 3
    P(f"    a_i = (4^i-1)/3: a_2..a_8 = {[a(i) for i in range(2, 9)]}; T(a_i) = 2^(2i-1): {all(Tp(a(i)) == 2 ** (2 * i - 1) for i in range(2, 40))}")
    P(f"    a_i = 0 mod 3 iff 3 | i (i <= 200): {all((a(i) % 3 == 0) == (i % 3 == 0) for i in range(1, 201))}; ord_9(4) = {next(k for k in range(1, 10) if pow(4, k, 9) == 1)}")
    P(f"    a_i (i >= 2) is never a power of two (odd and > 1): {all(a(i) % 2 == 1 and a(i) > 1 for i in range(2, 60))}; the orbit after a_i is the trunk down to 1, so no orbit contains two distinct a_i: disjointness HOLDS")
    # multiples of 3 have no odd preimage
    P(f"    odd preimage of n exists iff n = 2 mod 3 (n <= 10^5): {all(((2 * n - 1) % 3 == 0 and ((2 * n - 1) // 3) % 2 == 1) == (n % 3 == 2) for n in range(1, 100001))}")
    # decomposition B(a) = {a} u B(2a) u B((2a-1)/3): check on n <= 2^18 for a few roots via entry sets
    N0 = 1 << 18
    reach = {}
    def orbit_set(n, cap=1 << 40):
        s = []; v = n; seen = 0
        while v != 1 and seen < 5000:
            s.append(v); v = Tp(v); seen += 1
        s.append(1); return s
    import collections
    B = collections.defaultdict(set)
    for n in range(1, N0 + 1):
        for v in orbit_set(n): B[v].add(n)
    ok = True
    for aa in [5, 85, 341, 11, 23, 53, 113, 8, 32, 14]:
        rhs = {aa} | B[2 * aa] | (B[(2 * aa - 1) // 3] if aa % 3 == 2 else set())
        ok &= (B[aa] == rhs)
    P(f"    B(a) = {{a}} u B(2a) u B((2a-1)/3) (last term iff a = 2 mod 3), checked on n <= 2^18 for ten roots: {ok}")
    P(f"    5 in B(16)? {5 in B[16]} (T(5) = 8 < 16). So 'B(4^i) = union_(i' >= i) B(a_i')' is off by one: B(2^(2i-1)) = trunk-above u union_(i' >= i) B(a_i'),")
    P(f"      and B(4^i) = trunk-above u union_(i' >= i+1) B(a_i'). Check: B(16) \\ {{16,32,64,...}} == union_(i'>=3) B(a_i') on n <= 2^18: {B[16] - {2 ** k for k in range(4, 19)} == set().union(*[B[a(i)] for i in range(3, 30)])};"
      f" B(8) \\ {{8,16,32,...}} == union_(i'>=2) B(a_i'): {B[8] - {2 ** k for k in range(3, 19)} == set().union(*[B[a(i)] for i in range(2, 30)])}")
    P(f"    B(5) u B(32) u {{1,2,4,8,16}} == B(1) on n <= 2^18 (this is the set reaching 1, i.e. all n only under Collatz): {B[5] | B[32] | {1, 2, 4, 8, 16} == B[1]}; B(5) and B(32) disjoint: {not (B[5] & B[32])}")
    P("    Superadditivity of lower density on disjoint sets: liminf (A u B)/X >= liminf A/X + liminf B/X: standard. Proposition (iii) HOLDS.")
    P("    (iv) 'e_i sum to 1 iff almost every orbit reaches 1': the IF direction is proved (finite unions give lower density -> 1);")
    P("      the ONLY-IF direction needs dens B(2^(2M+1)) -> 0 as M -> oo (tightness of the entry index), which is not proved: OVERCLAIM.")
    # ---------- C2: the sheet cap ----------
    P("\nC2. Sheet cap: a sheet-blind proof of 'every positive cycle basin has lower density >= c' transfers to 3x-1, whose three")
    P("    known positive cycles have disjoint basins, so 3c <= 1 and c <= liminf of each; if the densities exist, c <= min: HOLDS.")
    P("    The value quoted for the cap (0.3248) is the [1,2^29] count of the {5,7,10} basin; that basin is RISING (0.32505 on [2^31,2^32)),")
    P("    so the cap, if the limit exists, is >= 0.3250 (HYP-9165's own status text says so; its title still says 0.3248).")
    # ---------- C3: finite-exact densities, own sieves ----------
    P("\nC3. FINITE-EXACT densities with my own sieves (C, unsigned __int128, memoised on the first value below n).")
    if have_gcc():
        NN = 1 << 30
        outm = compile_and_run(C_MINUS, 'audit_minus', NN, P)
        if outm:
            lines = outm.strip().split('\n')
            P("    3x-1 basins of {1} / {5,7,10} / {17,...,91}:")
            for ln in lines:
                if ln.startswith('[2^2') or ln.startswith('[2^1') or ln.startswith('cumulative') or ln.startswith('unset') or ln.startswith('N='):
                    P("      " + ln)
            # compare with the note
            row = [ln for ln in lines if ln.startswith('[2^28,2^29)')][0].split()
            P(f"      note [2^28,2^29): 87721223 87208132 93506101 (dens 0.326787 0.324876 0.348337); mine: {row[1]} {row[2]} {row[3]} -> {'EXACT MATCH' if row[1:4] == ['87721223', '87208132', '93506101'] else 'MISMATCH'}")
            cum29 = [ln for ln in lines if ln.startswith('cumulative [1,2^29]')][0]
            P(f"      note cumulative [1,2^29]: 0.3268569 0.3247598 0.3483833; mine: {cum29.split(':')[1]}")
            for s in (29,):
                r = [ln for ln in lines if ln.startswith(f'[2^{s},2^{s + 1})')][0].split()
                P(f"      note [2^29,2^30) (2-bit sieve): 0.3267479 0.3249565 0.3482957; mine: {r[4]} {r[5]} {r[6]}")
            spec_i = lines.index([ln for ln in lines if ln.startswith('ray-entry spectrum on [1,2^29]')][0])
            spec_j = lines.index([ln for ln in lines if ln.startswith('ray-entry spectrum on [1,N]')][0])
            spec = {}
            for ln in lines[spec_i + 1:spec_j]:
                t = ln.split(); spec[(int(t[0][4:]), int(t[1][2:]))] = float(t[2][5:])
            top = sorted(spec.items(), key=lambda kv: -kv[1])[:10]
            P("      top-10 first-ray-element densities on [1,2^29] (odd, k): " + str(top))
            note_top = [((7, 2), 0.1984), ((61, 2), 0.1914), ((1, 6), 0.1878), ((1, 4), 0.1336), ((7, 8), 0.0526), ((17, 1), 0.0460), ((25, 2), 0.0433), ((37, 4), 0.0335), ((5, 5), 0.0282), ((5, 7), 0.0264)]
            P(f"      note's top-10: {note_top}; agree to 4 decimals: {all(abs(spec.get(k, 0) - v) < 5e-5 for k, v in note_top)}")
            P(f"      ray of 1 entered at k = {sorted(k for (o, k) in spec if o == 1)} (never at k = 2, 8, 14: {all(k % 6 != 2 for (o, k) in spec if o == 1)})")
            P("      third basin, dyadic rows from 2^20 (note: 'flat to four decimals from 2^20 on'):")
            vals = []
            for s in range(20, 30):
                r = [ln for ln in lines if ln.startswith(f'[2^{s},2^{s + 1})')][0].split(); vals.append(float(r[6]))
            P(f"        {vals}; range {min(vals):.5f}..{max(vals):.5f}: flat to THREE decimals (0.348), not four.")
        outp = compile_and_run(C_PLUS, 'audit_plus', NN, P)
        if outp:
            lines = outp.strip().split('\n')
            P("    3x+1 trunk-entry basins e_i = dens B(a_i):")
            for ln in lines:
                if ln.startswith('[2^29') or ln.startswith('[2^28') or ln.startswith('  i=') or ln.startswith('  powers') or ln.startswith('N='):
                    P("      " + ln)
            cnts = {}
            for ln in lines:
                if ln.startswith('  i='):
                    t = ln.split(); cnts[int(t[0][2:])] = int(t[2][6:])
            P(f"      exact counts for 3 | i (doubling rays of a_i): i=3: {cnts.get(3)} (21*2^k <= 2^30: 26), i=6: {cnts.get(6)} (20), i=9: {cnts.get(9)} (14), i=12: {cnts.get(12)} (8)")
            P("      -> 'e_3 = e_6 = e_9 = e_12 = 0 exactly' is a 7-decimal rounding: the counts are the ray lengths, not zero (density 0 in the limit).")
            P(f"      note top range [2^29,2^30): 0.93794 0.02366 0.03780 for i = 2, 4, 5; mine: {[ln for ln in lines if ln.startswith('[2^29')][0]}")
            P(f"      cumulative e_2, e_4, e_5, e_7, e_8 (note: 0.9379633 0.0236466 0.0377892 0.0000805 0.0004854): {[round(cnts.get(i, 0) / NN, 7) for i in (2, 4, 5, 7, 8)]}")
            P(f"      e_5 > e_4 (341 owns a larger basin than 85): {cnts.get(5, 0) > cnts.get(4, 0)}")
    else:
        P("    gcc not available: pure-python fallback on [1, 2^20] for the 3x-1 basins only")
        N = 1 << 20; code = py_minus(N)
        c = [code[1:].count(b) for b in (1, 2, 3)]
        P(f"      3x-1 basins on [1,2^20]: {c} -> {[x / N for x in c]}")
    # ---------- C4: Krasikov-Lagarias ----------
    P("\nC4. Krasikov-Lagarias (2003) is CITED: for every a not divisible by 3, #{n <= X : a on the orbit of n} >= X^0.84 for all")
    P("    X >= X_0(a). The note states it in this form (root-uniform in type, not in the threshold): consistent; not verified here.")
    P("\nSECTION C VERDICT: (1) disjoint basins, (iii) no root-uniform positive proportion, (2) sheet cap: HOLD; (iv) has an index slip")
    P("  and an unproved 'only if'; (3) FINITE-EXACT tables reproduced exactly, with two wording corrections (e_3 'exactly 0', 'four decimals').")

# ================================================================================================
# SECTION D code
# ================================================================================================
"""Section D of the audit: HYP-9165 and the status-label / overclaim audit of the note."""

def sectionD(out):
    P = out
    P("=" * 78)
    P("SECTION D: HYP-9165 and the status-label audit")
    P("=" * 78)
    P("\nD1. HYP-9165 claims (verbatim scope): the natural densities of the three 3x-1 basins and of B(5) on 3x+1 exist,")
    P("    with values 0.3269 / 0.3248 / 0.3484 and 0.938; the strict ordering dens B({5,7,10}) < dens B({1}) < dens B({17,...});")
    P("    hence a sheet-blind cycle-uniform lower bound is capped at 0.3248 rather than 1/3. Status OPEN with FINITE-EXACT support.")
    # the note's own 2^32 rows
    rows = {28: (0.3267870, 0.3248756, 0.3483374), 29: (0.3267479, 0.3249565, 0.3482957), 30: (0.3267403, 0.3250124, 0.3482472), 31: (0.3267500, 0.3250474, 0.3482025)}
    inc2 = [rows[s + 1][1] - rows[s][1] for s in (28, 29, 30)]
    inc3 = [rows[s + 1][2] - rows[s][2] for s in (28, 29, 30)]
    P(f"    session's dyadic rows [2^s,2^(s+1)) s=28..31: increments of basin {{5,7,10}} per doubling {['%.1e' % v for v in inc2]},")
    P(f"      of basin {{17,...}} {['%.1e' % v for v in inc3]}; basin {{1}} spread {max(r[0] for r in rows.values()) - min(r[0] for r in rows.values()):.1e}.")
    P("    The status label OPEN is right: existence of the natural density of any Collatz basin is open on both sheets; the")
    P("    evidence supports the ORDERING to 2^32 (gaps 1.7e-3 and 2.2e-2) but the numerical VALUE 0.3248 in the title is stale:")
    P("    the {5,7,10} basin is still rising (0.32499 cumulative, 0.32505 in the top row); the hypothesis body already says")
    P("    'between 0.3250 and the limit'. A geometric extrapolation of the decaying increments (ratio ~0.7) would put the limit")
    P("    near 0.3251-0.3252; this is a heuristic, not part of the hypothesis.")
    P("    What would refute it: (a) a fourth positive 3x-1 cycle or a divergent 3x-1 orbit with positive upper density (the")
    P("    three basins then need not exhaust density 1; note the three cumulative densities sum to 1 - 2e-8 at 2^32 only because")
    P("    every n <= 2^32 was found to reach one of the three cycles); (b) a proof that some basin has different lower and")
    P("    upper densities; (c) a crossing of the {1} and {5,7,10} curves at larger N (the gap 1.7e-3 is closing at ~5e-5 per")
    P("    doubling with the increments themselves decaying, so a crossing would need the decay to stop).")
    P("    Internal inconsistency: HYP-9165's title says 'the sheet cap ... is 0.3248', its status text says 'between 0.3250 and")
    P("    the limit of the rising basin'; note section 3.1 and section 4 use 0.3248. The cap should be stated as 'the limit of")
    P("    the {5,7,10} basin density, >= 0.3250 on the evidence'.")

    P("\nD2. Status labels and overclaims found in the note and theorem files:")
    items = [
        ("THM-4516(5), note 2.5, THM-4516 title", "'Beal implies no map py +- m^r (r >= 3) has a free mixed shape (K,X) with K >= 3, sX >= 3, X >= 2'",
         "FALSE: freeness needs only (2^K - p^X) | m^r; 3y+125 has the free shape (5,3) with clock 5 | 125 (section B5). Correct statement: Beal implies no perfect-power CLOCK with K, sX, r >= 3."),
        ("note 2.5, THM-4516(5)", "'Its clock shadow is empty in any case: no term of x^4 + y^3 = z^17 can be a pure power of two while the other two are powers of one odd base'",
         "UNSUPPORTED: no argument is given and none is elementary; the shadow is the set of solutions with a power-of-two variable (section B6). Should be marked OPEN/UNVERIFIED like the statement itself."),
        ("note 3.2(iv), THM-4517", "'e_i ... sum to 1 if and only if almost every orbit reaches 1'",
         "Only 'if' is proved; 'only if' needs dens B(2^(2M+1)) -> 0 (tightness), not proved. Replace 'if and only if' by 'if' (converse open)."),
        ("note 3.2(iv)", "'B(4^i) = union_{i' >= i} B(a_{i'})'", "index slip: a_i maps to 2^(2i-1) < 4^i, so B(4^i) = trunk-above u union_{i' >= i+1} B(a_{i'}) (section C1)."),
        ("note 3.2(iv)", "'B(5) and B(32) partition the positive integers minus the powers of two'",
         "equivalent to the Collatz conjecture as stated (B(5) u B(32) u {1,2,4,8,16} is the set of n reaching 1, and 32, 64, ... lie in B(32)); replace by 'partition B(1) minus {1,2,4,8,16}'."),
        ("note 3.3, THM-4517(3)", "'e_3 = e_6 = e_9 = e_12 = 0 exactly'",
         "the counts are the doubling-ray lengths (26, 20, 14, 8 up to 2^30), i.e. densities ~1e-8 that round to 0.0000000; 'exactly' should be 'to seven decimals (ray lengths 26, 20, 14, 8)'. THM-4517's 'exact check of (1)' is a check that these basins are rays."),
        ("note 3.3, THM-4517(3)", "'the third [basin] is flat to four decimals from 2^20 on' / 'flat at 0.3483 to four decimals'",
         "the dyadic rows range 0.34804..0.34854: flat to THREE decimals (0.348)."),
        ("note 1.8, THM-4515(5)", "'PROVED lower bound of order 1/log m'",
         "the second-moment bound is exact for each m computed (<= 10^4; I extend to 3e4); the growth E[N_3^2] = O(log m) for all m is not proved in the note (local-CLT heuristic), so 'of order 1/log m' is PROVED only for the computed m. The slope is ~0.81 per unit log m, not 0.93."),
        ("note 1.8", "'E[N_3] -> sqrt(3)/(2 pi rho(1-rho)) (... 1.1847 at rho = log_3 2)'", "1.1847 is the value at rho = 12/19 (the shape used); at rho = log_3 2 exactly it is 1.1839."),
        ("note 1.2", "'The j twisted clocks 2^m - zeta^{-k} q^x are the eigenvalues'", "the eigenvalues of 2^m P - q^x I are 2^m zeta^k - q^x = zeta^k (2^m - zeta^{-k} q^x); equal up to a unit (harmless)."),
        ("note 1.5", "'|log(n'/n)| <= X|d|/(2(q n_min - |d|))'", "needs the implicit hypothesis q n_min > |d| (fails for 3x+11 through 1 and 3x+25 through 7, where the right side is negative); add 'for q n_min > |d|'."),
        ("note 1.9", "'minimum -14.81', 'minima range over [-27.10, -14.81]', 'minimum element exactly -17'", "these are the elements of least absolute value (-c_min/139, the maximum of a negative cycle); the true minima range over [-237.0, -72.9]. Say 'least-absolute-value element'."),
        ("note 2.2 table", "'21 = 7 x (3x+49 ...) scaled'", "should read '21 = 7 x 3, the 13x+49 cycle through 3 scaled by 7'."),
        ("THM-4516 census method (power_clocks.py)", "gmpy2.iroot loop r <= 64; power_data returns None otherwise",
         "a perfect power whose least prime exponent exceeds 64 would have been dropped from the q <= 201 / q = 3 censuses; my sieve covers every prime r <= log_3 N and finds nothing new, so the census statement HOLDS (method gap closed here)."),
        ("note 2.5 / THM-4516(5)", "'smallest open Beal signature (3,5,7)' alongside 'among (3,4,n) only n = 4, 5 are solved'",
         "internally consistent only if 'smallest' is by largest exponent (or if (3,4,6), (3,4,7), ... count as solved by reduction); under lexicographic or reciprocal-sum order (3,4,6) or (3,4,17) would precede (3,5,7). CITED content not checked (no browsing); the ordering sense should be stated."),
    ]
    for i, (where, phrase, verdict) in enumerate(items, 1):
        P(f"    D2.{i} [{where}] {phrase}\n        -> {verdict}")
    P("\nD3. Labels that are correct as used: PROVED for THM-4515(1)-(4) and THM-4517(1)-(2)(iii); FINITE-EXACT for the censuses and")
    P("    tables (all reproduced); CITED for Darmon-Granville, the ten Fermat-Catalan solutions, Krasikov-Lagarias, Barina; EMPIRICAL")
    P("    for the ~2.2/log m law; UNVERIFIED for the (4,3,17) statement; OPEN for HYP-9165. 'Collatz OPEN' is stated throughout.")

# ================================================================================================
# main
# ================================================================================================
def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join('05-knowledge', 'results', 'collatz_necklace_20260929_audit.out')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    f = open(out_path, 'w', encoding='utf-8', newline=chr(10))
    def P(s=''):
        s = str(s); print(s, flush=True); f.write(s + chr(10)); f.flush()
    import numpy, sympy
    try:
        gccv = subprocess.run(['gcc', '--version'], capture_output=True, text=True).stdout.splitlines()[0]
    except Exception:
        gccv = 'not available'
    P("Independent audit of collatz-necklace-20260929 (THM-4515/4516/4517, HYP-9165) -- auditor's own code, 2026-09-29")
    P(f"python {platform.python_version()}, numpy {numpy.__version__}, sympy {sympy.__version__}, gcc: {gccv}, platform {platform.platform()}")
    t0 = time.time()
    for sec in (sectionA, sectionB, sectionC, sectionD):
        t1 = time.time(); sec(P); P(f"[{sec.__name__}: {time.time() - t1:.0f} s]" + chr(10))
    P("=" * 78)
    P("SECTION E: SUMMARY OF VERDICTS")
    P("=" * 78)
    P("THM-4515: HOLDS (all identities exact on 1500 random words and 1206 fair splits; census K <= 18 reproduced; 17021-cycle facts")
    P("          reproduced; E[N_j] formula exact; second-moment bounds reproduced). Two wording corrections (1/log m 'PROVED' scope, 1.1847).")
    P("THM-4516: HOLDS WITH CORRECTION (census reproduced with a sieve that covers every exponent r, closing the r <= 64 gap of the")
    P("          session's method; free cycles and the Eisenstein square verified). FAILS as stated: the Beal-shadow consequence in (5)")
    P("          and in the title, and the 'clock shadow is empty' sentence.")
    P("THM-4517: HOLDS WITH CORRECTION ((1), (iii), (2) proved; FINITE-EXACT tables reproduced to the last digit; (iv) has an index slip,")
    P("          a Collatz-equivalent phrasing of the B(5)/B(32) partition and an unproved 'only if'; 'e_3 = 0 exactly' and 'four decimals').")
    P("HYP-9165: status OPEN is right; the title's 0.3248 is stale (rising basin, >= 0.3250 on the evidence).")
    P("OVERALL: SOUND WITH CORRECTIONS.")
    P(f"[total {time.time() - t0:.0f} s]")
    f.close()

if __name__ == '__main__':
    main()
