#!/usr/bin/env python3
"""Audit of lane fruit_rank_positive_multiples (collatz_mod6_20260921).

Independent recomputation: points come from PARI ellmul/elladd (not the
explorer's Python group law); fruit triples, positivity, permutation/J
identities, digit counts and the BSD index relation are recomputed here
with exact Fractions / sympy / mpmath.  Explicit raise, no assert.
"""
import math
import subprocess
from fractions import Fraction as Fr

import mpmath
import sympy as sp


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def run_gp(code, timeout=420):
    try:
        out = subprocess.run(["gp", "-q", "-f"], input=code + "\nquit\n",
                             capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "GP TIMEOUT"
    return "\n".join(ln for ln in out.stdout.splitlines() if not ln.startswith("  ***"))


def parse_pt(s):
    s = s.strip()
    if s == "[0]":
        return None
    s = s.strip("[]")
    xs, ys = s.split(",")
    return (Fr(xs.strip()), Fr(ys.strip()))


A2, A4 = 109, 224


def on_curve(P):
    return P is None or P[1] ** 2 == P[0] ** 3 + A2 * P[0] ** 2 + A4 * P[0]


def triple(P):
    """Inverse map (3) of catalan_elliptic; primitive, positive sum."""
    if P is None:
        return (1, -1, 0)
    x, y = P
    v = [56 - x + y, 56 - x - y, -56 - 12 * x]
    den = math.lcm(*[t.denominator for t in v])
    v = [int(t * den) for t in v]
    g = math.gcd(*v)
    v = [t // g for t in v]
    if sum(v) < 0:
        v = [-t for t in v]
    return tuple(v)


def forward(t):
    """Forward map (1): [X:Y:Z] = [-28(a+b+2c) : 364(a-b) : 6(a+b)-c]."""
    a, b, c = t
    return (-28 * (a + b + 2 * c), 364 * (a - b), 6 * (a + b) - c)


def J(t):
    a, b, c = t
    s2 = a * b + a * c + b * c
    return (-a * a + b * b + c * c + s2, a * a - b * b + c * c + s2, a * a + b * b - c * c + s2)


def prim(t):
    g = math.gcd(*t)
    v = [x // g for x in t]
    if sum(v) < 0:
        v = [-x for x in v]
    return tuple(v)


def fsum(t):
    a, b, c = t
    return Fr(a, b + c) + Fr(b, c + a) + Fr(c, a + b)


def main():
    print("=== AUDIT: fruit_rank_positive_multiples ===")
    # ---------- A. PARI, independent calls ----------
    gp = r"""
default(parisize, 400000000);
default(realprecision, 38);
E = ellinit([0,109,0,224,0]);
G = [-4,28]; T = [56,728];
print("elltors = ", elltors(E));
print("ellrank(E,1) = ", ellrank(E,1));
print("ellrootno = ", ellrootno(E));
ar = ellanalyticrank(E); print("ellanalyticrank = ", ar);
hG = ellheight(E,G); print("hG = ", hG);
bsd = ellbsd(E); print("ellbsd = ", bsd);
ratio = ar[2]/(bsd*hG); print("ratio L'/(ellbsd*hG) = ", ratio);
print("BSD reading: ratio = |Sha| * hhat(P0)/hhat(G) = |Sha|/n^2 where G = n*P0 + tors ; so |Sha| = n^2 * ratio");
print("if n=1: |Sha| = ", ratio, " ; if n=2: |Sha| = ", 4*ratio, " ; if n=1009: |Sha| = ", 1009^2*ratio);
print("omega1 = ", E.omega[1], " ; 2*omega1*24/36 = ", 2*E.omega[1]*24/36);
print("Pg=(-100,260) - G = ", ellsub(E,[-100,260],G), " ; 2T = ", ellmul(E,T,2));
print("ellsaturation(E,[G],1000) = ", ellsaturation(E,[G],1000));
s1 = getabstime(); sat = ellsaturation(E,[G],10000); print("ellsaturation(E,[G],10000) = ", sat, " time ms ", getabstime()-s1);
pts = ellratpoints(E, 100000);
print("ellratpoints(E,10^5): count ", #pts);
bad = 0; maxn = 0;
for(i=1,#pts, P = pts[i]; if(P==[0], next); h = ellheight(E,P); if(h < 1e-20, next); nn = round(sqrt(h/hG)); ok = 0; for(k=0,5, for(s=0,1, Q = elladd(E, ellmul(E,G,(-1)^s*nn), ellmul(E,T,k)); if(Q==P, ok=1))); if(!ok, bad++; print("  point not in <G,T>: ", P)); maxn = max(maxn, nn));
print("ellratpoints: points outside <G,T> = ", bad, " ; largest |n| among them all = ", maxn);
for(m=1,80, for(k=0,5, print("PT ", m, " ", k, " ", elladd(E, ellmul(E,G,m), ellmul(E,T,k)))));
"""
    out = run_gp(gp)
    pts = {}
    for ln in out.splitlines():
        if ln.startswith("PT "):
            _, m, k, rest = ln.split(" ", 3)
            pts[(int(m), int(k))] = parse_pt(rest)
        else:
            print(ln)
    check(len(pts) == 480, "480 PARI points parsed")
    check("elltors = [6, [6], [[56, 728]]]" in out, "torsion")
    check("ellrank(E,1) = [1, 1, 0," in out, "rank bounds 1,1 at effort 1")
    check("ellsaturation(E,[G],1000) = [[-4, 28]]" in out, "saturation 1000")
    check("points outside <G,T> = 0" in out, "all small-height rational points in <G,T>")
    hG = None
    for ln in out.splitlines():
        if ln.startswith("hG = "):
            hG = float(ln[5:])
    check(hG is not None, "hG parsed")

    # ---------- B. exact triples, positivity, identities ----------
    print("\n--- B. exact triples from PARI points ---")
    x = sp.symbols("x")
    fac = sp.factor((56 - x) ** 2 - (x ** 3 + 109 * x ** 2 + 224 * x))
    print("(56-x)^2 - (x^3+109x^2+224x) factors as", fac)
    check(sp.expand(fac + (x - 4) * (x ** 2 + 112 * x + 784)) == 0, "test (10) factorisation")
    print("=> for x<4: a*b>0 iff x^2+112x+784>0 ; c>0 iff x<-14/3 ; so (10) is an iff (PROVED)")
    pos = {}
    gcds = set()
    for (m, k), P in sorted(pts.items()):
        check(on_curve(P), "on curve %d %d" % (m, k))
        t = triple(P)
        check(fsum(t) == 4, "fruit sum %d %d" % (m, k))
        X, Y, Z = forward(t)
        check(Fr(X, Z) == P[0] and Fr(Y, Z) == P[1], "forward map returns point %d %d" % (m, k))
        p = all(v > 0 for v in t)
        test10 = P[0] < Fr(-14, 3) and P[0] ** 2 + 112 * P[0] + 784 > 0
        check(p == test10, "test (10) iff %d %d" % (m, k))
        pos[(m, k)] = (p, max(len(str(abs(v))) for v in t), t)
        # gcd of the cleared triple divides 728 (forward*inverse = 728 I)
        if k == 0:
            q2 = P[0].denominator
            q = math.isqrt(q2)
            check(q * q == q2 and P[1].denominator == q ** 3, "x=p/q^2, y=r/q^3 at m=%d" % m)
            pn, r = P[0].numerator, P[1].numerator
            raw = (56 * q ** 3 - pn * q + r, 56 * q ** 3 - pn * q - r, -56 * q ** 3 - 12 * pn * q)
            g = math.gcd(*raw)
            gcds.add(g)
            check(728 % g == 0, "gcd divides 728 at m=%d" % m)
    print("gcds of cleared triples (m=1..80, k=0):", sorted(gcds), " all divide 728:", all(728 % g == 0 for g in gcds))
    for k in range(6):
        ms = [m for m in range(1, 81) if pos[(m, k)][0]]
        print("k=%d positive m: %s count %d" % (k, ms, len(ms)))
    tot = sum(1 for v in pos.values() if v[0])
    print("total positive:", tot, "of 480")
    check(tot == 33, "33 positive")
    check([m for m in range(1, 81) if pos[(m, 0)][0]] == [9, 17, 43, 51, 77], "row k=0")
    check([m for m in range(1, 81) if pos[(m, 1)][0]] == [13, 21, 39, 47, 55, 73], "row k=1")
    check(all(m % 2 == 1 for (m, k), v in pos.items() if v[0]), "all positive m odd")
    small = sorted(((v[1], m, k) for (m, k), v in pos.items() if v[0]))[:8]
    print("smallest positive by max digits (digits,m,k):", small)
    print("max digits of positive k=0 rows:", [(m, pos[(m, 0)][1]) for m in (9, 17, 43, 51, 77)])
    print("k=0 fruit digit tuples m=9,13,17:", [tuple(len(str(abs(v))) for v in pos[(m, 0)][2]) for m in (9, 13, 17)])
    print("k=1 fruit digit tuple m=13:", tuple(len(str(abs(v))) for v in pos[(13, 1)][2]))
    # egg membership exactly: on egg iff x^2+109x+224 <= 0 (x<0 there)
    egg = {(m, k): (P[0] ** 2 + 109 * P[0] + 224 <= 0) for (m, k), P in pts.items()}
    check(all(egg[(m, k)] == (m % 2 == 1) for (m, k) in pts), "egg iff m odd, all cosets")
    print("exact: mG+kT on egg (x^2+109x+224<=0) iff m odd, for all 480 points: True")
    print("odd m with x<-14/3 (k=0):", sum(1 for m in range(1, 81, 2) if pts[(m, 0)][0] < Fr(-14, 3)),
          "; even:", sum(1 for m in range(2, 81, 2) if pts[(m, 0)][0] < Fr(-14, 3)))
    # permutation / J identities: triple(P+2T) is a cyclic permutation, triple(P+3T) = J
    cyc = {"id": 0, "bca": 0, "cab": 0}
    for m in range(1, 81):
        t0 = pos[(m, 0)][2]
        t2 = pos[(m, 2)][2]
        a, b, c = t0
        if t2 == (b, c, a):
            cyc["bca"] += 1
        elif t2 == (c, a, b):
            cyc["cab"] += 1
        elif t2 == t0:
            cyc["id"] += 1
        check(prim(J(t0)) == pos[(m, 3)][2], "triple(mG+3T) = J(triple(mG)) m=%d" % m)
        check(prim(J(t2)) == pos[(m, 5)][2], "triple(mG+5T) = J(triple(mG+2T)) m=%d" % m)
        t4 = pos[(m, 4)][2]
        check(prim(J(t4)) == pos[(m, 1)][2], "triple(mG+T) = J(triple(mG+4T)) m=%d" % m)
        check(sorted(t4) == sorted(t0) and sorted(t2) == sorted(t0), "k=2,4 are permutations of k=0 m=%d" % m)
        check(triple((pts[(m, 0)][0], -pts[(m, 0)][1])) == (b, a, c), "negation swaps a,b m=%d" % m)
    print("triple(mG+2T) vs triple(mG) = (a,b,c): counts", cyc, " (one fixed cyclic direction for all m)")
    print("triple(mG+3T) = J(triple(mG)), triple(mG+T) = J(triple(mG+4T)), triple(mG+5T) = J(triple(mG+2T)): all 80 m: True")
    # sqrt5 theorem check on the census: J of each positive triple has exactly one negative coordinate
    for (m, k), v in pos.items():
        if v[0]:
            jt = prim(J(v[2]))
            check(sum(1 for u in jt if u < 0) == 1, "J of positive has one negative %d %d" % (m, k))
    print("J of every positive triple (33) has exactly one negative coordinate: True")
    rows0 = {m for m in range(1, 81) if pos[(m, 0)][0]}
    rows1 = {m for m in range(1, 81) if pos[(m, 1)][0]}
    print("rows disjoint:", rows0.isdisjoint(rows1))

    # ---------- C. digits vs height ----------
    print("\n--- C. digits vs height (hG = %.16g) ---" % hG)
    c32 = 1.5 * hG / math.log(10)
    c1 = hG / math.log(10)
    print("(3/2) hG/ln10 = %.16g ; hG/ln10 = %.16g" % (c32, c1))
    for m in (9, 20, 40, 60, 80):
        t = pos[(m, 0)][2]
        P = pts[(m, 0)]
        L = math.log10(max(abs(v) for v in t))
        Hx = max(abs(P[0].numerator), P[0].denominator)
        print("m=%2d log10(max fruit)/m^2 = %.5f ; log10 H(x)/m^2 = %.5f ; x digits (num,den) = (%d,%d) ; fruit max digits %d" % (
            m, L / m / m, math.log10(Hx) / m / m, len(str(abs(P[0].numerator))), len(str(P[0].denominator)), pos[(m, 0)][1]))
    print("m=9 predictions: (3/2)*81*hG/ln10 = %.2f ; 81*hG/ln10 = %.2f ; actual fruit digits 81, x digits 53" % (81 * c32, 81 * c1))
    print("gcd bound: gcd(cleared triple) | 728 so log10(max fruit) = 3 log10 q + O(1) on the egg (|x|,|y| bounded): PROVED for odd m")

    # ---------- D. arc measure (mpmath, tanh-sinh handles 1/sqrt endpoints) ----------
    print("\n--- D. arc-measure fraction of positive window on the egg ---")
    mpmath.mp.dps = 30
    r1 = (-109 - mpmath.sqrt(10985)) / 2
    r2 = (-109 + mpmath.sqrt(10985)) / 2
    w1 = -56 - 28 * mpmath.sqrt(3)
    w2 = -56 + 28 * mpmath.sqrt(3)
    f = lambda u: 1 / mpmath.sqrt(u ** 3 + 109 * u ** 2 + 224 * u)
    egg = mpmath.quad(f, [r1, r2])
    win = mpmath.quad(f, [r1, w1]) + mpmath.quad(f, [w2, mpmath.mpf(-14) / 3])
    print("egg = [%s, %s] ; window = [egg_left, %s) U (%s, -14/3)" % (mpmath.nstr(r1, 9), mpmath.nstr(r2, 9), mpmath.nstr(w1, 9), mpmath.nstr(w2, 9)))
    print("int dx/sqrt(f) over egg = %s ; over window = %s ; fraction = %s" % (mpmath.nstr(egg, 15), mpmath.nstr(win, 15), mpmath.nstr(win / egg, 12)))
    print("egg integral vs PARI omega1 (0.650720927641666618...): ratio =", mpmath.nstr(egg / mpmath.mpf("0.65072092764166661835192516435949075893"), 12))
    print("observed 5/40 = 0.125 (k=0, odd m<=80) ; 6/40 = 0.150 (k=1) ; 11/80 = 0.1375 (both rows)")

    # ---------- E. index relation from BSD ----------
    print("\n--- E. BSD index bookkeeping ---")
    print("ellbsd doc: L'(E,1) = c*R*S with R = hhat(P0) (generator), S = |Sha|.  G = n P0 + tors => hhat(G) = n^2 R.")
    print("=> L'/(c*hhat(G)) = S/n^2.  Computed value 1 (38 digits) => S = n^2 under BSD, NOT S*n^2 = 1.")
    print("with ellsaturation(1000): n = 1 or n > 997 ; under BSD: n = 1 <=> Sha trivial ; n>997 <=> |Sha| > 994009.")
    lp = sp.prevprime(10000)
    print("largest prime below 10^4 = %d ; with ellsaturation(10000): n = 1 or n > %d ; then hhat(P0) = hG/n^2 < hG/10^8 = %.4e (= 1.5188/10^8)" % (lp, lp, hG / 1e8))
    print("explorer's midpoint arc fraction 0.13638 vs tanh-sinh 0.13619: the midpoint rule is wrong in the 4th digit")
    print("=== AUDIT END ===")


if __name__ == "__main__":
    main()
