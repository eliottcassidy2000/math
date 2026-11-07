# audit_H1: independent checks of THM-4900(6) and HYP-9490 Lemmas A-B: cover polynomial by subset convolution vs the odd-cycle set-partition (Grinberg-Stanley) form at even and odd n, sw_1 = 2c3, tower formula vs brute-force H of explicit towers, U_(d,r) three ways (L-polynomials, exponential walk W_r, Taylor series of M^r e^y).
import math
import random
import subprocess
import sys
from fractions import Fraction
from audit_tournament_lib import (rand_tour, classes, transitive, gate, c3, cover_poly, ham_cycles)

HP = sys.argv[1] if len(sys.argv) > 1 else "./hp_count"
rng = random.Random(4900)
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def upper(A):
    n = len(A)
    return "".join("1" if A[i][j] else "0" for i in range(n) for j in range(i + 1, n))


def brute_H(tours):
    inp = "".join(f"{len(A)} {upper(A) if len(A) > 1 else '0'}\n" for A in tours)
    out = subprocess.run([HP], input=inp, capture_output=True, text=True, check=True).stdout.split()
    return [int(x) for x in out]


# ---------- odd-cycle set-partition form of P_T(q) ----------
def cycle_form_coeffs(A):
    """sum over set partitions of V into singletons (weight q) and odd blocks B, |B| >= 3, of weight 2*HC(T[B])*q;
    returns coefficient list in q."""
    n = len(A)
    size = 1 << n
    wB = {}
    for S in range(1, size):
        k = bin(S).count("1")
        if k >= 3 and k % 2 == 1:
            verts = [v for v in range(n) if (S >> v) & 1]
            sub = [[A[a][b] for b in verts] for a in verts]
            hc = ham_cycles(sub)
            if hc:
                wB[S] = 2 * hc
        elif k == 1:
            wB[S] = 1
    # DP over subsets: f[S] = polynomial (dict exp->coef); fix lowest vertex in block
    f = [None] * size
    f[0] = {0: 1}
    for S in range(1, size):
        low = S & -S
        rest = S ^ low
        acc = {}
        B = rest
        while True:
            blk = B | low
            if blk in wB and f[S ^ blk] is not None:
                for e, c in f[S ^ blk].items():
                    acc[e + 1] = acc.get(e + 1, 0) + c * wB[blk]
            if B == 0:
                break
            B = (B - 1) & rest
        f[S] = acc
    co = [0] * (n + 1)
    for e, c in f[size - 1].items():
        co[e] += c
    return co


# 1. cover polynomial structure and sw_1 = 2 c3
tested = {}
pool = []
for n in range(1, 8):
    cl = classes(n)
    pool += [(n, A) for A in cl]
for n, cnt in ((8, 120), (9, 25), (10, 8)):
    pool += [(n, rand_tour(n, rng)) for _ in range(cnt)]
gs_checked = 0
for n, A in pool:
    co = cover_poly(A)
    check(co[n] == 1, "not monic")
    check(all(co[k] == 0 for k in range(n + 1) if (n - k) % 2 == 1), "parity")
    check(co[0] == 0, "constant term")
    if n >= 2:
        check(co[n - 2] == 2 * c3(A), f"sw_1 != 2c3 at n={n}")
    if n <= 8:
        check(cycle_form_coeffs(A) == co, f"odd-cycle set-partition form fails at n={n}")
        gs_checked += 1
    tested[n] = tested.get(n, 0) + 1
print("(6) P_T monic, parity n, P_T(0)=0, coefficient of q^(n-2) = 2 c3:", tested)
print(f"(6) P_T(q) = sum over partitions into singletons and odd cycles of prod (2 HC q or q) [Grinberg-Stanley form], "
      f"incl. EVEN n: {gs_checked} tournaments of orders 1..8, mismatches {sum(1 for f in fails if 'set-partition' in f)}")

# ---------- universal weights U_(d,r) three ways ----------
def L_poly(co):
    """L f(q) = q[f(q+1) + q f(q) - f(q-1)] on coefficient lists."""
    deg = len(co) - 1
    out = [0] * (deg + 3)
    for e, c in enumerate(co):
        if c == 0:
            continue
        out[e + 2] += c                         # q * q * q^e
        for i in range(1, e + 1, 2):            # q[(q+1)^e - (q-1)^e] = 2 sum_{i odd} C(e,i) q^(e-i+1)
            out[e - i + 1] += 2 * math.comb(e, i) * c
    return out


def U_poly(d, r):
    co = [0] * d + [1]
    for _ in range(r):
        co = L_poly(co)
    return sum(co)


def W_walk(rmax):
    """M^r e^y = sum_k W_r(k) e^(ky), M e^(ky) = k e^((k+1)y) + k^2 e^(ky) - k e^((k-1)y)."""
    Ws = [{1: 1}]
    for r in range(rmax):
        W = Ws[-1]
        new = {}
        for k, c in W.items():
            if k == 0 or c == 0:
                continue
            new[k + 1] = new.get(k + 1, 0) + k * c
            new[k] = new.get(k, 0) + k * k * c
            new[k - 1] = new.get(k - 1, 0) - k * c
        # keep k = 0 entries (constant term) since they matter for d = 0
        for k, c in W.items():
            if k == 0:
                new[0] = new.get(0, 0) + 0  # M(1) = 0
        Ws.append(new)
    return Ws


def U_walk(Ws, d, r):
    return sum(c * (k ** d if (k, d) != (0, 0) else 1) for k, c in Ws[r].items())


def taylor_M_power(r, K):
    """Taylor coefficients a_0..a_(K-2r) of M^r e^y at 0, exact."""
    a = [Fraction(1, math.factorial(m)) for m in range(K + 1)]
    sinh = [Fraction(1, math.factorial(i)) if i % 2 == 1 else Fraction(0) for i in range(K + 1)]
    for _ in range(r):
        L = len(a) - 1
        d1 = [(m + 1) * a[m + 1] for m in range(L)]
        d2 = [(m + 2) * (m + 1) * a[m + 2] for m in range(L - 1)]
        new = []
        for m in range(L - 1):
            conv = sum(sinh[i] * d1[m - i] for i in range(1, m + 1))
            new.append(d2[m] + 2 * conv)
        a = new
    return a


Ws = W_walk(260)
mism = 0
for r in range(0, 21):
    tay = taylor_M_power(r, 12 + 2 * r)
    for d in range(0, 13):
        up = U_poly(d, r)
        uw = U_walk(Ws, d, r)
        ut = tay[d] * math.factorial(d)
        if not (up == uw == ut):
            mism += 1
check(mism == 0, "U three ways")
print(f"Lemma A: (L^r q^d)(1) = sum_k W_r(k) k^d = d^d/dy^d [M^r e^y](0) for r <= 20, d <= 12: mismatches {mism}")

# H(Gamma T) = P(1) + P(2) - P(0): r = 1 of Lemma A
for n in range(1, 7):
    for A in classes(n)[:20]:
        co = cover_poly(A)
        Pv = lambda q: sum(c * q ** e for e, c in enumerate(co))
        check(brute_H([gate(A)])[0] == Pv(1) + Pv(2) - Pv(0), "H(Gamma T) = P(1)+P(2)-P(0)")

# ---------- tower formula vs brute force ----------
batch, expect = [], []
for n in range(0, 8):
    for r in range(0, 8):
        if n + 2 * r > 19 or n + 2 * r < 1:
            continue
        # U_(n,r) = H(Gamma^r TT_n) by brute force (TT_0 empty: Gamma(empty) = TT_2)
        X = transitive(n) if n > 0 else None
        if X is None:
            if r == 0:
                continue
            X = transitive(2)
            rr = r - 1
        else:
            rr = r
        for _ in range(rr):
            X = gate(X)
        batch.append(X)
        expect.append(("U", n, r, U_walk(Ws, n, r)))
cores = []
for n in range(1, 8):
    for t in range(4):
        cores.append(rand_tour(n, rng))
cores += classes(5)
for A in cores:
    n = len(A)
    co = cover_poly(A)
    for r in range(0, 7):
        if n + 2 * r > 19:
            break
        X = A
        for _ in range(r):
            X = gate(X, rng)
        pred = sum(co[n - 2 * j] * U_walk(Ws, n - 2 * j, r) for j in range(0, n // 2 + 1))
        batch.append(X)
        expect.append(("T", n, r, pred))
got = brute_H(batch)
bad = [(e, g) for e, g in zip(expect, got) if e[3] != g]
check(not bad, f"tower formula / U brute force mismatches: {bad[:3]}")
print(f"(6) brute-force H of explicit towers: {sum(1 for e in expect if e[0]=='U')} values U_(d,r) = H(Gamma^r TT_d) and "
      f"{sum(1 for e in expect if e[0]=='T')} values H(Gamma^r T) = sum_j sw_j U_(n-2j,r) (total order <= 19): mismatches {len(bad)}")

# ---------- Lemma B ----------
v1 = v2 = v3 = 0
Utab = [[U_walk(Ws, d, rr) for d in range(0, 43)] for rr in range(0, 201)]
for r in range(1, 201):
    for e in range(0, 41):
        U = lambda d, rr: Utab[rr][d]
        a, b, c_ = U(e, r), U(e + 2, r), U(e + 2, r - 1)
        if not (a > 0):
            v1 += 1
        if not (a <= b):
            v1 += 1
        if not (a >= c_):
            v2 += 1
        rho = Fraction(U(e, r - 1), U(e + 2, r - 1))
        ratio = Fraction(a, c_)
        if not (1 <= ratio <= 1 + 2 ** e * rho):
            v3 += 1
check(v1 == v2 == v3 == 0, "Lemma B")
print(f"Lemma B on e <= 40, r <= 200: violations of positivity/B1 {v1}, B2 {v2}, B3 {v3}")

# K1 tower and one-defect counts
k1 = [U_walk(Ws, 1, k) for k in range(9)]
print("H(Gamma^k K1), k = 0..8:", k1)
onedef = [k1[k + 1] - (2 * k + 3) * k1[k] for k in range(8)]
print("one-defect counts m(Gamma^k K1) = H(Gamma^(k+1) K1) - (2k+3) H(Gamma^k K1), k = 0..7:", onedef)
print("FAILURES:", len(fails))
