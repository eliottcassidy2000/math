#!/usr/bin/env python3
"""death-star-2026-07-16-S20 (HYP-7021/THM-890): THE RELATION-LATTICE RATE THEOREM referee.

THEOREM (paper derivation, this session): for owner e with others f = (f1..f5) + stationary,
the a-twisted owner sum T_e(a) = sum_j u_e(j) w^{aj}  (w = e(1/7); = the owner part of
S(ae)) satisfies the EXACT finite identity

  T_e(a) = 7e * SUM_{(c,c0)} Ghat(c,c0) * SUM_{k in Z_{7e}^5 : SUM k_i f_i + (c0+a)e = 0 mod 7e}
           PROD_i hhat_{c_i}(k_i)

with hhat_c(k) = [k = c mod 7] * (1/e) * (1 - e(-k/7)) / (1 - e(-k/(7e)))  (k != 0 mod 7e),
hhat_c(0) = [c = 0], and Ghat = the Z_7^6 DFT of the dipole rule
G(sigma, r) = g(sigma, r) - g(sigma, r-1), g = [ {0} u {sigma_i} u {r} = Z7 \\ {s} ].

MAIN TERM (k = 0 forces c = 0, c0 = -a):  7e * Ghat(0, -a) = e * (1 - w^a) * mhat*_s(a)
— the independent limit EXACTLY (grid marginals exactly uniform; no correction).
ERROR = the nonzero terms = SMALL ADDITIVE RELATIONS  sum k_i f_i = l*e (mod 7e) with
7-adic type constraints, weighted by the folded-decay products.

Referee parts:
  R1: hhat closed form vs direct DFT (unit test).
  R2: Ghat(0,-a) vs (1/7)(1-w^a)-normalized mhat* (constant check incl. sign/units).
  R3: THE IDENTITY: T_e(a) computed directly (endpoint words) vs main term + enumerated
      small-relation terms (|k_i| <= K); the residual must drop as K grows (and match the
      tail bound scale). Run on incoherent + mildly-coherent cores.
  R4: the coherence meter: cores WITH a planted small relation (f1 + f2 - f3 = e) must
      show the predicted relation-term dominating the deviation.
"""
from fractions import Fraction as Fr
from math import gcd, pi
from itertools import product as iproduct
import cmath, sys, time

def e2(x):  # e(x)
    return cmath.exp(2j * pi * x)

W7 = [e2(Fr(c, 7)) for c in range(7)]

# ---------------- owner words (exact, klein-verbatim R_s logic, per-owner) ----------------

def owner_word(E, s, e):
    """u_e on Z_{7e}: net sign at j/(7e) attributed to owner e = min-owner convention."""
    bps = sorted(set(Fr(k, 7 * f) for f in E if f > 0 for k in range(7 * f)) | {Fr(0), Fr(1)})
    word = {}
    prev_in = None
    events = []
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((f * mid % 1) * 7) for f in E)
        cur = (len(occ) == 6) and (s not in occ)
        if prev_in is None:
            first = cur
        else:
            if cur != prev_in:
                events.append((bps[i], +1 if cur else -1))
        prev_in = cur
    # wraparound: compare first cell and last cell
    if prev_in is not None and prev_in != first:
        events.append((Fr(1), +1 if first else -1))
    for p, sg in events:
        p = p % 1
        os_ = [f for f in E if f > 0 and (p * 7 * f).denominator == 1]
        if not os_: continue
        own = min(os_)
        if own != e: continue
        j = int(p * 7 * e) % (7 * e)
        word[j] = word.get(j, 0) + sg
    return word

def T_direct(E, s, e, a):
    w = owner_word(E, s, e)
    return sum(v * W7[(a * j) % 7] for j, v in w.items())

# ---------------- the identity's ingredients ----------------

def hhat(c, k, e):
    """Z_{7e}-DFT coefficient of y -> w^{c*floor(y/e)} at frequency k (exact formula)."""
    P = 7 * e
    k %= P
    if k % 7 != c % 7:
        return 0j
    if k == 0:
        return 1.0 + 0j if c % 7 == 0 else 0j
    num = 1 - e2(Fr(-k, 7))
    den = 1 - e2(Fr(-k, P))
    return (num / den) / e

def hhat_direct(c, k, e):
    P = 7 * e
    return sum(W7[(c * (y // e)) % 7] * e2(Fr(-k * y, P)) for y in range(P)) / P

def g_pattern(sigma, r, s):
    occ = {0} | set(sigma) | {r % 7}
    return 1 if (len(occ) == 6 and (s % 7) not in occ) else 0

_GV = {}
def gv_table(s):
    if s in _GV:
        return _GV[s]
    tab = []
    for sigma in iproduct(range(7), repeat=5):
        row = [g_pattern(sigma, r, s) for r in range(7)]
        tab.append((sigma, row))
    _GV[s] = tab
    return tab

def Ghat_one(s, c, c0):
    """Single DFT coefficient of G at character (c, c0) — on demand, 7^6 ops."""
    z = 0j
    for sigma, row in gv_table(s):
        base = sum(ci * si for ci, si in zip(c, sigma)) % 7
        for r in range(7):
            G = row[r] - row[(r - 1) % 7]
            if G:
                z += G * W7[(-base - c0 * r) % 7]
    return z / 7 ** 6

def mstar(s, a):
    A = 360 / 16807; B = 120 / 16807
    tot = 0j
    for c in range(7):
        if c == s % 7: continue
        Ac = 0.0 if c == 0 else A
        tot += (Ac + B) * W7[(a * c) % 7]
    return tot

if __name__ == "__main__":
    t0 = time.time()
    print("R1: hhat closed form vs direct DFT (e = 5, all c, k sample)")
    ok = True
    for c in range(7):
        for k in [0, 1, 2, 5, 7, 12, 33]:
            d = abs(hhat(c, k, 5) - hhat_direct(c, k, 5))
            if d > 1e-9: ok = False; print(f"  FAIL c={c} k={k} diff={d:.2e}")
    print(f"  {'PASS' if ok else 'FAIL'}")
    sys.stdout.flush()

    print("R2: main-term constant: 7e*Ghat(0,-a) vs e*(1-w^a)*mhat*(a)")
    s = 3
    for a in range(1, 7):
        c0 = (-a) % 7
        gh = Ghat_one(s, (0, 0, 0, 0, 0), c0)
        e = 10
        main = 7 * e * gh
        pred = e * (1 - W7[a % 7]) * mstar(s, a)
        rr = main / pred if abs(pred) > 0 else float('nan')
        print(f"  a={a}: 7e*Ghat = {main:.6f}   e(1-w^a)m* = {pred:.6f}   "
              f"main/pred = {rr:.6f}")
    sys.stdout.flush()

    print("\nR3: THE IDENTITY on real cores — T_e(a) vs main + small-relation terms")
    def relation_terms(E, s, e, a, K):
        """Enumerate k with |k_i| <= K, (k,c0) != (0,-a)-main, sum k_i f_i + (c0+a)e = 0 mod 7e;
        weights Ghat(c(k),c0)*prod hhat. c_i = k_i mod 7 forced; c0 free 0..6."""
        others = [f for f in E if f > 0 and f != e]
        P = 7 * e
        tot = 0j
        cnt = 0
        rng = range(-K, K + 1)
        for k in iproduct(rng, repeat=5):
            kf = sum(ki * fi for ki, fi in zip(k, others))
            for c0 in range(7):
                if (kf + (c0 + a) * e) % P != 0:
                    continue
                if all(ki == 0 for ki in k) and c0 == (-a) % 7:
                    continue  # main term
                c = tuple(ki % 7 for ki in k)
                gh = Ghat_one(s, c, c0)
                if abs(gh) < 1e-14:
                    continue
                wgt = 1.0 + 0j
                for ki in k:
                    wgt *= hhat(ki % P, ki % P, e) if False else hhat(ki % 7, ki % P, e)
                tot += 7 * e * gh * wgt
                cnt += 1
        return tot, cnt
    for (Etag, E) in [("incoherent", [0, 100, 138, 191, 284, 413, 588]),
                      ("planted f1+f2-f3=e", [0, 97, 140, 137, 284, 413, 100])]:
        E = sorted(set(E))
        e = 100 if 100 in E else max(E)
        s = 3
        for a in [1, 2]:
            td = T_direct(E, s, e, a)
            main = e * (1 - W7[a % 7]) * mstar(s, a)
            for K in [0, 2, 4]:
                rel, cnt = relation_terms(E, s, e, a, K)
                resid = td - main - rel
                print(f"  {Etag} e={e} a={a} K={K}: |T|={abs(td):.3f} |main|={abs(main):.3f} "
                      f"|rel({cnt})|={abs(rel):.3f} |resid|={abs(resid):.3f}")
            sys.stdout.flush()
    print(f"[{time.time()-t0:.1f}s]")
