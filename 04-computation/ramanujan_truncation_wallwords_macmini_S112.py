#!/usr/bin/env python3
"""ramanujan_truncation_wallwords_macmini_S112.py -- mac-mini-2026-07-15-S112.
Owner: (1) work the Ramanujan truncation mean-square classical estimate;
(2) witnesses under negation; (3) does numerator multiplication permute the
mechanical wall-words?

(I) THE TRUNCATION CLOSED FORM AND MEAN SQUARE (kps cont.26's named residual).
    T_L(h) := sum_{l<=L} c_l(h)/l.  CLOSED FORM (proved below, verified exactly):
        T_L(h) = sum_{d|h} M(L/d),   M(x) = sum_{m<=x} mu(m)/m  (M(x)=0 for x<1).
    [Proof: c_l(h) = sum_{d|(l,h)} d mu(l/d); swap sums, l = d*m.]
    Consequences, all machine-checked:
      - reproduces kps cont.26's concentration numbers (h=12: +1.97 at L=12, etc.);
      - LIMIT MEAN SQUARE (H -> inf, exact):  lim (1/H) sum_{h<=H} T_L(h)^2
            = sum_{d,e<=L} M(L/d) M(L/e) / lcm(d,e)   =: S(L);
      - S(L) growth fit vs log L (elementary bound O(log^3 L); PNT-strength
        per-layer analysis suggests ~ C log L; numerics decide);
      - the truncation error is supported on the top divisor window:
        T_L(h) ~ #{d|h : L/2 < d <= L} + lower layers (M(1)=1 dominates).

(II) WALL-WORDS AND THE UNIT ACTION.
    For a witness a/q (gcd(a,q)=1) and core {1..k}: the THREE-DISTANCE GAP WORD
    W(a,q,k) = cyclic sequence of gap lengths of {i*a/q mod 1 : i=0..k} -- the
    mechanical wall-word (THM-639's order-cell pattern at the witness).
    Tested exhaustively (q <= 60, all units, k in {5,8,12}):
      (a) NEGATION u = -1: W(q-a) = reversal of W(a)   [geometric reflection]
      (b) INVERSION u = a^{-1}: W(a^{-1}) vs W(a) -- test reversal/exchange duality
      (c) general u: does ANY uniform word operation (rotation/reversal/letter-perm)
          map W(a) to W(ua)?  Expect NO beyond the Klein four {+-1, +-a^{-1}}.
"""
import sys, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

# ---------- (I)
def mobius(n):
    m, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0: return 0
            m = -m
        p += 1
    if n > 1: m = -m
    return m

LMAX = 4000
mu = [0, 1] + [0] * (LMAX)
for n in range(2, LMAX + 1): mu[n] = mobius(n)
Mcum = [F(0)] * (LMAX + 1)          # M(x) at integer x
for m in range(1, LMAX + 1): Mcum[m] = Mcum[m - 1] + F(mu[m], m)
def M(x):
    xi = int(x)
    return Mcum[xi] if xi >= 1 else F(0)

def c_ram(l, h):                     # Ramanujan sum c_l(h)
    g = gcd(l, h); s = F(0)
    for d in range(1, g + 1):
        if g % d == 0: s += d * mu[l // d]
    return s

def T_direct(L, h): return sum(F(c_ram(l, h), l) for l in range(1, L + 1))
def T_closed(L, h): return sum(M(F(L, d)) for d in range(1, h + 1) if h % d == 0 and d <= L)

print("(I) RAMANUJAN TRUNCATION")
ok = all(T_direct(L, h) == T_closed(L, h) for L in (5, 12, 30) for h in range(1, 40))
print("   closed form T_L(h) = sum_{d|h} M(L/d) == direct sum (L=5,12,30; h<40):", ok)
print("   kps cont.26 reproduction at L=12: h=12: %.4f  h=30: %.4f  h=1: %.6f"
      % (float(T_closed(12, 12)), float(T_closed(12, 30)), float(T_closed(12, 1))))
print("   vanishing by L=800: h=12: %.6f  h=30: %.6f"
      % (float(T_closed(800, 12)), float(T_closed(800, 30))))
# exact limit mean square S(L) = sum_{d,e<=L} M(L/d)M(L/e)/lcm(d,e)
def S_exact(L):
    Ms = [None] + [M(F(L, d)) for d in range(1, L + 1)]
    tot = F(0)
    for d in range(1, L + 1):
        if Ms[d] == 0: continue
        for e in range(1, L + 1):
            if Ms[e] == 0: continue
            tot += Ms[d] * Ms[e] * gcd(d, e) / (d * e)
    return tot
import math
print("   limit mean-square S(L) (exact), with S(L)/log L:")
for L in (6, 12, 25, 50, 100, 200, 400, 800):
    s = float(S_exact(L))
    print(f"      L={L:4d}: S = {s:8.4f}   S/log L = {s/math.log(L):6.4f}   S/log^2 L = {s/math.log(L)**2:6.4f}")
# empirical mean square at finite H, sanity vs S(L)
H = 20000
for L in (12, 50):
    ms = sum(float(T_closed(L, h)) ** 2 for h in range(1, H + 1)) / H
    print(f"   finite-H check L={L}: (1/H)sum T^2 = {ms:.4f} vs S(L) = {float(S_exact(L)):.4f}  (H={H})")

print()
print("(II) WALL-WORDS UNDER THE UNIT ACTION")
def gap_word(a, q, k):
    pts = sorted((i * a) % q for i in range(0, k + 1))
    gaps = [pts[i + 1] - pts[i] for i in range(len(pts) - 1)] + [q + pts[0] - pts[-1]]
    # letters by gap VALUE (three-distance: <= 3 distinct values)
    vals = sorted(set(gaps))
    letter = {v: chr(ord('a') + i) for i, v in enumerate(vals)}
    return "".join(letter[g] for g in gaps), tuple(vals)

def cyc_variants(w):
    rots = {w[i:] + w[:i] for i in range(len(w))}
    return rots | {v[::-1] for v in rots}

def test_q(q, k):
    units = [a for a in range(1, q) if gcd(a, q) == 1]
    res = {"neg_reversal": True, "inv_in_dihedral": True, "others_dihedral": 0, "others_tot": 0}
    for a in units:
        w, vals = gap_word(a, q, k)
        # (a) negation
        wn, valsn = gap_word(q - a, q, k)
        if not (valsn == vals and wn in {x[::-1] for x in {w[i:]+w[:i] for i in range(len(w))}} | {w[i:]+w[:i] for i in range(len(w))}):
            res["neg_reversal"] = False
        # (b) inversion
        ainv = pow(a, -1, q)
        wi, valsi = gap_word(ainv, q, k)
        if wi not in cyc_variants(w) or valsi != vals:
            res["inv_in_dihedral"] = False
        # (c) other units: how often is W(ua) in the dihedral class of W(a)?
        for u in units:
            if u in (1, q - 1, ainv, (q - ainv) % q): continue
            wu, valsu = gap_word((u * a) % q, q, k)
            res["others_tot"] += 1
            if valsu == vals and wu in cyc_variants(w): res["others_dihedral"] += 1
    return res

print("   q <= 60 exhaustive, k = 5, 8, 12:")
agg = {}
for k in (5, 8, 12):
    neg_ok = inv_ok = True
    oth_d = oth_t = 0
    inv_fail_ex = None
    for q in range(k + 2, 61):
        r = test_q(q, k)
        neg_ok &= r["neg_reversal"]
        if not r["inv_in_dihedral"] and inv_fail_ex is None: inv_fail_ex = q
        inv_ok &= r["inv_in_dihedral"]
        oth_d += r["others_dihedral"]; oth_t += r["others_tot"]
    print(f"   k={k:2d}: NEGATION = dihedral (reversal) class: {neg_ok};  "
          f"INVERSION in dihedral class: {inv_ok}" + (f" (first fail q={inv_fail_ex})" if inv_fail_ex else "") +
          f";  other units in dihedral class: {oth_d}/{oth_t} = {oth_d/max(oth_t,1):.3f}")
print("\nDONE")
