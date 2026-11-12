#!/usr/bin/env python3
"""Skeptic check 6: the CRT form of (B2), exactly (integers and Fractions), for both parities of w.

For every P_u = (a, b) (|a|, |b| < 2^s1, b odd) and P_v = (a', b') (|a'|, |b'| < 2^s2, b' odd):
  direction (k, k') = (X, Xi)/j, j = +-gcd, k' >= 1.
 U-side (P_u fixed):  e = 2^ds k' a' - k a
   (i)   k b' - k' b = 2^s2 e
   (ii)  |e| < |k| + 2^ds k'
   (iii) e/m = -sgn(k) a/(2^ds k') - b inv(2^s1;|k|)/|k|  (mod 1),  m = 2^ds k' |k|
   (iv)  e determines (a', b'); #(a', b') with direction (k, k') <= 6M/m + 1,  M = max(|k|, k')   [strongest K]
   (v)   || a/(2^ds k') + theta_k || < 3M/(2^ds |k| k'),  theta_k = sgn(k) b inv(2^s1;|k|)/|k|
 V-side (P_v fixed):  e = k a - 2^ds k' a'
   (i')  k' b = k b' + 2^s2 e
   (ii') |e| < |k| + 2^ds k'
   (iii') e/(|k| k') = -sgn(k) [2^ds a'/k + theta'_(k')] (mod 1),  theta'_(k') = b' inv(2^s2; k')/k'
   (iv') #(a, b) with direction (k, k') <= 6M/(|k| k') + 1
   (v')  || 2^ds a'/k + theta'_(k') || < 3M/(|k| k')
"""
import sys, math
from fractions import Fraction
from collections import defaultdict

def dist_int(fr):
    f = fr - math.floor(fr)
    return min(f, 1 - f)

def inv_mod(x, m):
    if m == 1:
        return 0
    return pow(x, -1, m)

def main():
    out = []
    def P(*a):
        s = " ".join(str(t) for t in a)
        print(s, flush=True)
        out.append(s)
    wlist = [int(t) for t in sys.argv[1].split(",")]
    for w in wlist:
        s1 = (w + 1) // 2; s2 = w // 2; ds = s1 - s2
        viol = defaultdict(int)
        cntU = defaultdict(int); cntV = defaultdict(int)
        eU = defaultdict(set); eV = defaultdict(set)
        npairs = 0
        maxratU = Fraction(0); maxratV = Fraction(0)
        for a in range(-(1 << s1) + 1, 1 << s1):
            for b in range(-(1 << s1) + 1, 1 << s1, 2):
                for ap in range(-(1 << s2) + 1, 1 << s2):
                    X = ap * (1 << s1) + b
                    for bp in range(-(1 << s2) + 1, 1 << s2, 2):
                        Xi = a * (1 << s2) + bp
                        g = math.gcd(X, Xi)
                        j = g if Xi > 0 else -g
                        k = X // j; kp = Xi // j
                        assert kp >= 1 and k % 2 and kp % 2 and math.gcd(k, kp) == 1
                        npairs += 1
                        M = max(abs(k), kp)
                        sg = 1 if k > 0 else -1
                        # U side
                        e = (1 << ds) * kp * ap - k * a
                        if k * bp - kp * b != (1 << s2) * e:
                            viol["i"] += 1
                        if not abs(e) < abs(k) + (1 << ds) * kp:
                            viol["ii"] += 1
                        m = (1 << ds) * kp * abs(k)
                        iv = inv_mod(1 << s1, abs(k))
                        # (iii): e/m + sgn(k) a/(2^ds k') + b inv/|k| integer
                        val = Fraction(e, m) + Fraction(sg * a, (1 << ds) * kp) + Fraction(b * iv, abs(k))
                        if val.denominator != 1:
                            viol["iii"] += 1
                        theta = Fraction(sg * b * iv, abs(k))
                        if not dist_int(Fraction(a, (1 << ds) * kp) + theta) < Fraction(3 * M, (1 << ds) * abs(k) * kp):
                            viol["v"] += 1
                        keyU = (a, b, k, kp)
                        cntU[keyU] += 1
                        if e in eU[keyU]:
                            viol["iv-e-not-injective"] += 1
                        eU[keyU].add(e)
                        # V side
                        e2 = k * a - (1 << ds) * kp * ap
                        if kp * b != k * bp + (1 << s2) * e2:
                            viol["i'"] += 1
                        if not abs(e2) < abs(k) + (1 << ds) * kp:
                            viol["ii'"] += 1
                        m2 = abs(k) * kp
                        iv2 = inv_mod(1 << s2, kp)
                        thetap = Fraction(bp * iv2, kp)
                        val2 = Fraction(e2, m2) + sg * (Fraction((1 << ds) * ap, k) + thetap)
                        if val2.denominator != 1:
                            viol["iii'"] += 1
                        if not dist_int(Fraction((1 << ds) * ap, k) + thetap) < Fraction(3 * M, abs(k) * kp):
                            viol["v'"] += 1
                        keyV = (ap, bp, k, kp)
                        cntV[keyV] += 1
                        if e2 in eV[keyV]:
                            viol["iv'-e-not-injective"] += 1
                        eV[keyV].add(e2)
        for (a, b, k, kp), c in cntU.items():
            M = max(abs(k), kp); m = (1 << ds) * kp * abs(k)
            r = Fraction(c) / (Fraction(6 * M, m) + 1)
            if r > 1:
                viol["iv"] += 1
            if r > maxratU:
                maxratU = r
        for (ap, bp, k, kp), c in cntV.items():
            M = max(abs(k), kp); m2 = abs(k) * kp
            r = Fraction(c) / (Fraction(6 * M, m2) + 1)
            if r > 1:
                viol["iv'"] += 1
            if r > maxratV:
                maxratV = r
        P("w=%d (s1=%d s2=%d ds=%d): (P_u, P_v) pairs %d ; violations %s ; max count/(6M/m+1): U-side %.4f, V-side %.4f ; "
          "max #points per direction: U %d, V %d"
          % (w, s1, s2, ds, npairs, dict(viol) if viol else "{}", float(maxratU), float(maxratV),
             max(cntU.values()), max(cntV.values())))
    with open(sys.argv[2], "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
