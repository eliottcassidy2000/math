#!/usr/bin/env python3
"""Two-anchor theory of the residual first-reset-2 branch (mac-mini-2026-10-07-twoanchor).

U(x) = (3x+1)/2^v2(3x+1) on odd x; Terras T(x) = x/2 or (3x+1)/2.
Residual source: n = 2^K t - 1, t odd, K >= 2, with first reset letter 2 (v2(3^K t - 1) = 1).
Run end x = 2*3^(K-1) t - 1 = 1 + 2^j u (u odd, j >= 3); the source then has J = floor((j-1)/2) letters 2.
Deletion child h_D = (n+1)/2^D - 1, run end y_D = (x+1)/3^D - 1.

Parts:
 A  run transparency at the anchors c_a = 1/(2^a - 3) (a = 1: -1, a = 2: +1), random exact tests
 B  the +1 barrier: exact limit chains (x* = 1, y* = (2-3^D)/3^D) for D <= DMAX never absorb; debts k_inf(D)
 C  universal residual state: after a long two-run the D = 3, 4 chain is exactly (3, 1 - 27)
 D  doubly-uniform collisions (uniform in K and J): two explicit word pairs, algebra + random sources
 E  Mersenne line: K = 1 + 2^m q families certified by K -> K-3 uniformly in m
All checks are exact integer / rational arithmetic.
"""
import random, sys
from fractions import Fraction as Fr

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000

def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
def U(x):
    y = 3*x + 1; return y >> v2(y)
def letter(x): return v2(3*x + 1)
def F(word, x):
    for a in word: x = (3*x + 1) / Fr(2**a)
    return x
CHECKS = 0
def ok(cond, msg):
    global CHECKS
    if not cond: raise AssertionError(msg)
    CHECKS += 1

# ---------------- A: run transparency ----------------
def part_A(rnd):
    for a in range(1, 9):
        c = Fr(1, 2**a - 3)
        ok(3*c + 1 == 2**a * c, "anchor fixed point")
        for _ in range(300):
            k = rnd.randint(1, 12); L = rnd.randint(1, 6)
            # v with U-word a^L: v - c has v2 >= a*L + 1 ; build v = c + 2^(aL+1) * s  (as an integer: need c 2-adic)
            M = 2**(a*L + 1 + 40)
            cm = (pow(2**a - 3, -1, M)) % M          # c mod M
            s = rnd.getrandbits(30) | 1
            v = (cm + 2**(a*L + 1) * s) % M + M * rnd.getrandbits(8)
            if v % 2 == 0: continue
            # u with u - c = 3^k (v - c)   (mod M; then lift to an integer congruent mod M)
            u = (cm + pow(3, k) * (v - cm)) % M + M * rnd.getrandbits(8)
            if u % 2 == 0: continue
            uu, vv = u, v
            for i in range(L):
                ok(letter(vv) == a and letter(uu) == a, "common letter a in the run")
                uu, vv = U(uu), U(vv)
            # relation persists modulo the remaining precision
            ok((uu - cm - pow(3, k) * (vv - cm)) % 2**(40) == 0, "anchored relation persists")
    print("A  run transparency: anchors c_a, a=1..8, exact random runs ok")

# ---------------- B: the +1 barrier ----------------
def limit_chain(D):
    a = 2 - 3**D; d = D; s = 0; oddy = 0
    while d > 0:
        if a & 1: a = (a + 3**(d-1)) // 2; d -= 1; oddy += 1
        else: a //= 2
        s += 1
    N = a; land = s; y = N
    while y not in (0, 1, 2):
        if y & 1: y = (3*y + 1) // 2; oddy += 1
        else: y //= 2
        s += 1
    if y == 0: return dict(D=D, land=land, N=N, kind='to0')
    xs = 1 if s % 2 == 0 else 2
    k = D + (s + 1)//2 - oddy
    return dict(D=D, land=land, N=N, entry=s, inphase=(y == xs), k=k)
def part_B():
    recs = [limit_chain(D) for D in range(1, DMAX + 1)]
    for r in recs:
        ok(r['N'] >= 0, "limit child lands on a nonnegative integer")
        if r.get('kind') == 'to0': ok(r['D'] in (1, 2), "only D=1,2 go to 0")
        else: ok(not (r['inphase'] and r['k'] == 0), "barrier: no absorption at +1")
    inph = [r for r in recs if r.get('inphase')]
    kmin = min(inph, key=lambda r: r['k'])
    ok(kmin['k'] == 3 and kmin['D'] in (3, 4), "minimal in-phase debt 3 at D=3,4")
    small = [(r['D'], r.get('k'), r.get('inphase')) for r in recs[:20]]
    print(f"B  +1 barrier: D<= {DMAX}: none absorbs; to0 D=1,2; in-phase {len(inph)}, out-of-phase {DMAX-2-len(inph)};"
          f" min debt {kmin['k']} at D={kmin['D']}; max landing N {max(r['N'] for r in recs)}")
    print("   first debts (D, k_inf, inphase):", small)
    return recs

# ---------------- C: universal residual state ----------------
def residual_source(rnd, K, j, bits=60):
    M = 1 << (j - 1)
    t0 = pow(3, -(K-1), M)
    while True:
        t = t0 + M * rnd.getrandbits(bits)
        if t & 1 and v2(3**(K-1)*t - 1) == j - 1: return t
def part_C(rnd):
    for D, kexp in ((3, 3), (4, 3), (5, 6), (6, 6)):
        for j in (16, 25, 40, 64):
            for _ in range(40):
                K = rnd.randint(D + 1, 60)
                t = residual_source(rnd, K, j)
                x = 2*3**(K-1)*t - 1; y = (x + 1)//3**D - 1
                ok(v2(x - 1) == j and v2(3**K*t - 1) == 1, "residual source")
                J = (j - 1)//2
                u, v, k = x, y, D
                for s in range(2*J):
                    k += (u & 1) - (v & 1); u, v = T(u), T(v)
                ok(k == kexp and u - 1 == 3**k * (v - 1), "state (k, 1-3^k) at run end")
    print("C  universal residual state: D=3,4 -> (3, 1-27), D=5,6 -> (6, 1-729) at the end of every two-run tested")

# ---------------- D: doubly-uniform collisions ----------------
PAIRS = {1: ((1, 10, 1), (1, 1, 2, 2, 3, 3)), 2: ((10, 2), (3, 1, 1, 3, 4))}
def part_D(rnd):
    # algebra: clearing segment and anchor
    X = Fr(rnd.getrandbits(40)); Y = (X + 1)/27 - 1
    ok(F((2, 2, 2), X) - 1 == 27*(F((4, 1, 1), Y) - 1), "clearing identity F222(x)-1 = 27(F411(y)-1)")
    for r, (bu, bv) in PAIRS.items():
        ok(len(bv) - len(bu) == 3 and sum(bu) == sum(bv), "length/total")
        ok(F(bu, Fr(1)) == F(bv, Fr(1)), "collision at +1")
        Z = Fr(rnd.getrandbits(50)); W = 1 + 27*(Z - 1)
        ok(F(bu, W) == F(bv, Z), "affine identity from the anchored state")
    # native classes: w' (Y_end = 1 + 2^r w') for which both post-run words are actual
    classes = {}
    for r, (bu, bv) in PAIRS.items():
        S = sum(bu); mod = 1 << (S + 2)
        good = []
        for wp in range(1, mod, 2):
            Ye = 1 + (1 << r) * wp; Xe = 1 + 27 * (Ye - 1)
            def actual(z, word):
                for a in word:
                    if letter(z) != a: return False
                    z = U(z)
                return True
            # last letters are 'uncapped': require all but the last exactly, last at least as given and equal endpoints
            if actual(Xe, bu[:-1]) and actual(Ye, bv[:-1]):
                zx = Xe
                for a in bu[:-1]: zx = U(zx)
                zy = Ye
                for a in bv[:-1]: zy = U(zy)
                if v2(3*zx+1) >= bu[-1] and v2(3*zy+1) >= bv[-1] and (3*zx+1)*2**bv[-1] == (3*zy+1)*2**bu[-1]:
                    good.append(wp)
        classes[r] = (mod, good)
    # random sources across K, J
    cnt = 0
    for r, (bu, bv) in PAIRS.items():
        mod, good = classes[r]
        ok(len(good) >= 1, "nonempty native class")
        for _ in range(120):
            K = rnd.choice([4, 5, 6, 9, 17, 40, 101, 333, 1001]); J = rnd.choice([3, 4, 5, 8, 13, 30, 77, 150])
            j = 2*J + r
            wp = rnd.choice(good) + mod * rnd.getrandbits(20)
            # u = 3^(3-J) w'  (mod 2^big): choose t with run end x = 1 + 2^j u
            Mb = 1 << (j + 64)
            u = (pow(3, 3 - J, Mb) * wp) % Mb if J >= 3 else None
            t = (pow(3, -(K-1), Mb) * ((1 + (1 << j) * u + 1) // 2)) % Mb
            if t % 2 == 0: t += Mb
            n = (1 << K)*t - 1; h = (1 << (K-3))*t - 1
            x = 2*3**(K-1)*t - 1
            ok(v2(x - 1) == j and v2(3**K*t - 1) == 1, "residual source with prescribed J")
            # source word 1^(K-1) 2^J bu ; child word 1^(K-4) 4 1 1 2^(J-3) bv  (last letters uncapped)
            src = [1]*(K-1) + [2]*J + list(bu); chd = [1]*(K-4) + [4, 1, 1] + [2]*(J-3) + list(bv)
            zn, zh = n, h
            for a in src[:-1]: ok(letter(zn) == a, "source word actual"); zn = U(zn)
            for a in chd[:-1]: ok(letter(zh) == a, "child word actual"); zh = U(zh)
            ok(U(zn) == U(zh) and h < n, "common odd endpoint, smaller child")
            cnt += 1
    print(f"D  doubly-uniform collisions: 2 explicit pairs, native classes w' mod 2^(S+2): "
          f"r=1 {len(classes[1][1])} of {classes[1][0]//2}, r=2 {len(classes[2][1])} of {classes[2][0]//2}; {cnt} random sources (K up to 1001, J up to 150) certified")
    return classes

# ---------------- E: Mersenne line ----------------
def part_E():
    # odd K with m = v2(K-1) >= 4 (J >= 3): direct check of merge 2^K - 1 ~ 2^(K-3) - 1 via the anchored chain
    hits = {}
    for K in range(17, 4000, 2):
        m = v2(K - 1)
        if m < 4: continue
        x = 2*3**(K-1) - 1; y = (x + 1)//27 - 1
        j = v2(x - 1); J = (j - 1)//2
        ok(j == 3 + m, "Mersenne j = 3 + v2(K-1)")
        u, v, k = x, y, 3
        dep = None
        for s in range(2*J + 400):
            if u == v and k == 0: dep = s - 2*J; break
            k += (u & 1) - (v & 1); u, v = T(u), T(v)
        if dep is not None: hits[K] = (m, dep)
    tot = sum(1 for K in range(17, 4000, 2) if v2(K - 1) >= 4)
    print(f"E  Mersenne: of {tot} odd K<4000 with v2(K-1)>=4, {len(hits)} merge with 2^(K-3)-1 within post-run depth 400;"
          f" examples {sorted(hits.items())[:6]}")
    return hits

if __name__ == "__main__":
    rnd = random.Random(20261007)
    part_A(rnd); part_B(); part_C(rnd); part_D(rnd); part_E()
    print(f"ALL CHECKS PASSED ({CHECKS} assertions)")
