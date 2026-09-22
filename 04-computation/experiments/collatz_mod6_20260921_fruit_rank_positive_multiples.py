#!/usr/bin/env python3
"""collatz_mod6_20260921 -- lane fruit_rank_positive_multiples.

Fruit curve E: y^2 = x^3 + 109 x^2 + 224 x  (Bremner--Macleod N=4).
(1) PARI: torsion, rank (2-descent bounds), analytic rank, height of G=(-4,28),
    saturation of <G>, BSD consistency.
(2) Exact multiples mG (1<=m<=80) and torsion translates mG+kT (0<=k<=5),
    fruit triples via the projective inverse map of
    catalan_elliptic_20260921_elliptic.md eq. (3), exact verification of
    a/(b+c)+b/(c+a)+c/(a+b)=4, positivity census, digit counts vs height.
(3) Structural statement (E(Q) = Z x Z/6) and the sparse positive subset.
Audit 2026-09-22: mpmath quadrature for the arc measure, saturation to 10^4,
ellratpoints membership check added.

No Collatz map. Exact Fractions throughout; PARI only for rank/height data.
Explicit raise, no assert.
"""
import math
import subprocess
import sys
from fractions import Fraction as Fr

A2, A4 = 109, 224  # y^2 = x^3 + A2 x^2 + A4 x
INF = None


def on_curve(P):
    if P is INF:
        return True
    x, y = P
    return y * y == x ** 3 + A2 * x * x + A4 * x


def neg(P):
    return INF if P is INF else (P[0], -P[1])


def add(P, Q):
    if P is INF:
        return Q
    if Q is INF:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 + y2 == 0:
            return INF
        lam = (3 * x1 * x1 + 2 * A2 * x1 + A4) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - A2 - x1 - x2
    y3 = lam * (x1 - x3) - y1
    return (x3, y3)


def mul(m, P):
    R = INF
    Q = P
    while m > 0:
        if m & 1:
            R = add(R, Q)
        Q = add(Q, Q)
        m >>= 1
    return R


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def fruit_triple(P):
    """Projective inverse map (3): [a:b:c] = [56Z-X+Y : 56Z-X-Y : -56Z-12X].
    Returns primitive integer triple with positive coordinate sum (sum = 14(4-x)
    times a positive scalar; sum is zero only at x=4, which is torsion)."""
    if P is INF:
        return (1, -1, 0)
    x, y = P
    a, b, c = 56 - x + y, 56 - x - y, -56 - 12 * x
    den = math.lcm(a.denominator, b.denominator, c.denominator)
    a, b, c = int(a * den), int(b * den), int(c * den)
    g = math.gcd(math.gcd(abs(a), abs(b)), abs(c))
    a, b, c = a // g, b // g, c // g
    s = a + b + c
    if s < 0:
        a, b, c = -a, -b, -c
    return (a, b, c)


def fruit_sum(t):
    a, b, c = t
    if 0 in (b + c, c + a, a + b):
        return None
    return Fr(a, b + c) + Fr(b, c + a) + Fr(c, a + b)


def F4(t):
    a, b, c = t
    s = a + b + c
    return s ** 3 - 6 * s * (a * b + a * c + b * c) + 7 * a * b * c


def positive(t):
    return all(v > 0 for v in t)


def positivity_test(P):
    """catalan_elliptic (10): x < -14/3 and x^2+112x+784 > 0."""
    x = P[0]
    return x < Fr(-14, 3) and x * x + 112 * x + 784 > 0


def digits(n):
    return len(str(abs(n)))


def run_gp(code):
    out = subprocess.run(["gp", "-q", "-f"], input=code + "\nquit\n",
                         capture_output=True, text=True, timeout=400)
    lines = [ln for ln in out.stdout.splitlines() if not ln.startswith("  ***")]
    return "\n".join(lines)


def main():
    print("=== collatz_mod6_20260921 lane fruit_rank_positive_multiples ===")
    print("E: y^2 = x^3 + 109x^2 + 224x ; G = (-4,28) ; T = (56,728)")

    # ---------- Section 1: PARI data ----------
    print("\n--- S1: PARI 2.17 data (verbatim outputs) ---")
    gp_code = r"""
default(parisize, 200000000);
default(realprecision, 38);
E = ellinit([0,109,0,224,0]);
print("disc = ", E.disc, " = ", factor(E.disc));
gr = ellglobalred(E);
print("conductor = ", gr[1], " = ", factor(gr[1]));
print("minimal model change [u,r,s,t] = ", gr[2]);
print("local data (kodaira, cp) = ", gr[5]);
print("kodaira codes 4+nu = I_nu: ", [k[2] for k in gr[5]], " -> I_nu with nu = ", [k[2]-4 for k in gr[5]]);
print("minimal discriminant = ", ellminimalmodel(E).disc, " = ", factor(ellminimalmodel(E).disc));
print("tamagawa product = ", elltamagawa(E));
print("elltors = ", elltors(E));
G = [-4,28]; T = [56,728];
print("ellisoncurve(G) = ", ellisoncurve(E,G));
print("ellorder(G) = ", ellorder(E,G));
print("ellheight(G) = ", ellheight(E,G));
print("ellheight(2G) = ", ellheight(E, ellmul(E,G,2)));
print("ellheight(9G) = ", ellheight(E, ellmul(E,G,9)));
print("ellheight(9G)/81 = ", ellheight(E, ellmul(E,G,9))/81);
print("ellheight(G+T) = ", ellheight(E, elladd(E,G,[56,728])));
print("ellheight(T) = ", ellheight(E,[56,728]));
r = ellrank(E);
print("ellrank(E) = ", r);
print("ellrank(E,,[G]) = ", ellrank(E,,[G]));
print("ellsaturation(E,[G],1000) = ", ellsaturation(E,[G],1000));
print("ellsaturation(E,[[-100,260]],1000) = ", ellsaturation(E,[[-100,260]],1000));
print("ellsaturation(E,[G],10000) = ", ellsaturation(E,[G],10000));
pts = ellratpoints(E, 100000);
bad = 0; maxn = 0;
for(i=1,#pts, P = pts[i]; if(P==[0], next); h = ellheight(E,P); if(h < 1e-20, next); nn = round(sqrt(h/ellheight(E,G))); ok = 0; for(k=0,5, for(s=0,1, Q = elladd(E, ellmul(E,G,(-1)^s*nn), ellmul(E,T,k)); if(Q==P, ok=1))); if(!ok, bad++; print("  point not in <G,T>: ", P)); maxn = max(maxn, nn));
print("ellratpoints(E,10^5): ", #pts, " points with naive x-height <= 10^5 ; outside <G,T>: ", bad, " ; largest |n| among them: ", maxn);
print("ellrootno = ", ellrootno(E));
print("ellanalyticrank = ", ellanalyticrank(E));
print("L(E,1) = ", elllseries(E,1));
print("L'(E,1) = ", ellL1(E,1));
print("E.omega = ", E.omega);
print("ellbsd(E) = ", ellbsd(E));
print("Sha_an = L'(E,1)/(ellbsd(E)*ellheight(G)) = ", ellL1(E,1)/(ellbsd(E)*ellheight(E,G)));
print("9G = ", ellmul(E,G,9));
print("ellmul(E,G,80) x-numerator digits = ", #Str(numerator(ellmul(E,G,80)[1])));
"""
    gp_out = run_gp(gp_code)
    print(gp_out)
    check("elltors = [6, [6], [[56, 728]]]" in gp_out, "torsion Z/6 by T")
    check("ellrank(E) = [1, 1, 0, [[-100, 260]]]" in gp_out, "ellrank bounds 1,1")
    check("ellsaturation(E,[G],1000) = [[-4, 28]]" in gp_out, "G saturated to 1000")
    check("ellsaturation(E,[G],10000) = [[-4, 28]]" in gp_out, "G saturated to 10000")
    check("outside <G,T>: 0" in gp_out, "all rational points of naive height <= 10^5 lie in <G,T>")
    check("ellorder(G) = 0" in gp_out, "G infinite order")
    hG = None
    for ln in gp_out.splitlines():
        if ln.startswith("ellheight(G) = "):
            hG = float(ln.split("= ")[1])
    check(hG is not None, "parsed hG")

    # ---------- Section 2: exact multiples ----------
    print("\n--- S2: exact group law, torsion, PARI generator vs G ---")
    G = (Fr(-4), Fr(28))
    T = (Fr(56), Fr(728))
    check(on_curve(G) and on_curve(T), "G,T on curve")
    tors = [mul(k, T) for k in range(6)]
    print("kT, k=0..5:", [None if P is INF else (int(P[0]), int(P[1])) for P in tors])
    check(mul(6, T) is INF and mul(3, T) == (Fr(0), Fr(0)) and mul(2, T) == (Fr(4), Fr(52)),
          "torsion orders")
    # torsion points are exactly the zero-denominator fruit points (catalan_elliptic sec 4)
    for k in range(6):
        t = fruit_triple(tors[k])
        print("  kT fruit triple k=%d: %s  pair-sum-zero: %s" % (k, t, 0 in (t[0] + t[1], t[1] + t[2], t[0] + t[2])))
        check(0 in (t[0] + t[1], t[1] + t[2], t[0] + t[2]), "torsion = zero pair sum")
    Pg = (Fr(-100), Fr(260))
    check(on_curve(Pg), "PARI generator on curve")
    rel = None
    for s in (1, -1):
        for k in range(6):
            cand = add(mul(1, (G[0], s * G[1])), tors[k])
            if cand == Pg:
                rel = (s, k)
    print("PARI generator (-100,260) = %s" % ("(%+dG) + %dT" % rel if rel else "NOT in G+torsion"))
    check(rel is not None, "PARI generator differs from G by torsion/sign")

    print("\n--- S2: multiples mG, 1<=m<=80 (exact) ---")
    print("m | x(mG) digits(num,den) | fruit digits (a,b,c) | positive? | test(10) | log10(maxcoord)/m^2 | F4==0 | sum==4")
    mults = {}
    P = INF
    pos_m = []
    ratios = []
    for m in range(1, 81):
        P = add(P, G)
        mults[m] = P
        check(on_curve(P), "mG on curve m=%d" % m)
        t = fruit_triple(P)
        fs = fruit_sum(t)
        check(fs == 4, "fruit sum 4 at m=%d" % m)
        check(F4(t) == 0, "F4 zero at m=%d" % m)
        pos = positive(t)
        check(pos == positivity_test(P), "positivity test (10) agrees m=%d" % m)
        if pos:
            pos_m.append(m)
        dmax = max(digits(v) for v in t)
        ratio = math.log10(max(abs(v) for v in t)) / (m * m)
        ratios.append(ratio)
        xd = (digits(P[0].numerator), digits(P[0].denominator))
        if m <= 20 or pos or m == 80:
            print("%2d | %s | %s | %s | %s | %.4f | %s | %s" % (
                m, xd, tuple(digits(v) for v in t), pos, positivity_test(P), ratio, F4(t) == 0, fs == 4))
    print("all-positive m in 1..80:", pos_m)
    for m in pos_m:
        kids = [(c, (c in pos_m) if c <= 80 else "out of range") for c in (3 * m - 1, 3 * m, 3 * m + 1)]
        print("  prime_shells tree children of %dG (3m-1,3m,3m+1) positive?: %s" % (m, kids))
    print("count positive:", len(pos_m), " parity of positive m:", sorted(set(m % 2 for m in pos_m)))
    print("odd multiples with x<-14/3 (bounded component) count:", sum(1 for m in range(1, 81, 2) if mults[m][0] < Fr(-14, 3)))
    print("even multiples with x<-14/3 count:", sum(1 for m in range(2, 81, 2) if mults[m][0] < Fr(-14, 3)))
    t9 = fruit_triple(mults[9])
    print("9G fruit triple:", t9)
    print("9G digits:", tuple(digits(v) for v in t9))
    print("9G point x:", mults[9][0], "\n     y:", mults[9][1])
    A_rep = 154476802108746166441951315019919837485664325669565431700026634898253202035277999
    B_rep = 36875131794129999827197811565225474825492979968971970996283137471637224634055579
    C_rep = 4373612677928697257861252602371390152816537558161613618621437993378423467772036
    check(sorted(t9) == sorted((A_rep, B_rep, C_rep)), "9G is the repaired triple")
    print("9G == repaired input triple (A, B/10, C/10):", sorted(t9) == sorted((A_rep, B_rep, C_rep)))

    # digit growth vs height
    print("\n--- S2: digit growth vs canonical height ---")
    print("hhat(G) (PARI, Cremona normalisation) =", hG)
    print("hhat(G)/log(10) =", hG / math.log(10))
    print("(3/2) hhat(G)/log(10) =", 1.5 * hG / math.log(10))
    for m in (9, 20, 40, 60, 80):
        t = fruit_triple(mults[m])
        L = math.log10(max(abs(v) for v in t))
        xnum = math.log10(abs(mults[m][0].numerator)) if mults[m][0] != 0 else 0.0
        xden = math.log10(mults[m][0].denominator)
        print("m=%2d: log10(max fruit coord)/m^2 = %.5f ; log10|x num|/m^2 = %.5f ; log10 x den/m^2 = %.5f ; log10 H(x)/m^2 = %.5f" % (
            m, L / m / m, xnum / m / m, xden / m / m, max(xnum, xden) / m / m))
    print("prediction: fruit digits ~ (3/2) m^2 hhat(G)/log10 ; x-height digits ~ m^2 hhat(G)/log10")
    print("m=9 predicted fruit digits (3/2)*81*hG/ln10 = %.2f ; actual max = %d" % (1.5 * 81 * hG / math.log(10), max(digits(v) for v in t9)))
    print("m=9 predicted with factor 1: 81*hG/ln10 = %.2f" % (81 * hG / math.log(10)))

    # ---------- torsion translates ----------
    print("\n--- S2: torsion translates mG + kT, 1<=m<=80, 0<=k<=5 ---")
    print("(convention: sign s=+1 only; -(mG+kT) gives the same triple with a,b swapped)")
    census = {}
    pos_list = []
    for k in range(6):
        for m in range(1, 81):
            P = add(mults[m], tors[k])
            t = fruit_triple(P)
            check(fruit_sum(t) == 4, "translate fruit sum m=%d k=%d" % (m, k))
            pos = positive(t)
            check(pos == positivity_test(P), "translate test (10) m=%d k=%d" % (m, k))
            census[(m, k)] = pos
            if pos:
                pos_list.append((m, k, max(digits(v) for v in t)))
    for k in range(6):
        ms = [m for m in range(1, 81) if census[(m, k)]]
        print("k=%d: positive m = %s (count %d)" % (k, ms, len(ms)))
    print("total positive (m,k) in 1..80 x 0..5:", len(pos_list), "of", 480)
    print("smallest positive by digit count:", sorted(pos_list, key=lambda z: z[2])[:8])
    # even multiples: which component?
    print("x(2G) =", mults[2][0], " x(2G+3T)=", add(mults[2], tors[3])[0])
    print("x(G+3T) =", add(G, tors[3])[0], " (catalan_elliptic: U-translate of G is (-56,-392))")
    # component bookkeeping: bounded real component = x in [x1,x2] with x^3+109x^2+224x>=0 negative roots
    disc_roots = "roots of x^2+109x+224: (-109 +- sqrt(10985))/2"
    print("egg component: x between", disc_roots)
    r1 = (-109 - math.sqrt(10985)) / 2
    r2 = (-109 + math.sqrt(10985)) / 2
    print("numerically: [%.6f, %.6f]" % (r1, r2))
    print("G on egg (x=-4 in interval):", r1 <= -4 <= r2, " ; 2G x=%s on egg: %s" % (mults[2][0], r1 <= float(mults[2][0]) <= r2))
    odd_on_egg = all(r1 <= float(mults[m][0]) <= r2 for m in range(1, 81, 2))
    even_off_egg = all(float(mults[m][0]) > r2 for m in range(2, 81, 2))
    print("all odd mG on egg:", odd_on_egg, " all even mG on unbounded component:", even_off_egg)
    # positive window inside egg: x < -14/3 and x^2+112x+784>0 -> x < -56-28sqrt3 or x > -56+28sqrt3
    w1 = -56 - 28 * math.sqrt(3)
    w2 = -56 + 28 * math.sqrt(3)
    print("positive window (10): x in [egg_left, %.6f) U (%.6f, -14/3)" % (w1, w2))
    print("egg_left = %.6f ; window lengths %.6f and %.6f ; egg length %.6f" % (r1, w1 - r1, -14 / 3 - w2, r2 - r1))
    # Haar-measure proportion of positive arc on egg (elliptic-log measure dx/|y|).
    # Audit fix: the earlier midpoint rule mishandled the 1/sqrt endpoint
    # singularities (gave 0.13638); tanh-sinh quadrature (mpmath) is used instead
    # and its egg integral is checked against PARI's real period omega1.
    import mpmath
    mpmath.mp.dps = 30
    R1 = (-109 - mpmath.sqrt(10985)) / 2
    R2 = (-109 + mpmath.sqrt(10985)) / 2
    W1 = -56 - 28 * mpmath.sqrt(3)
    W2 = -56 + 28 * mpmath.sqrt(3)
    fq = lambda u: 1 / mpmath.sqrt(u ** 3 + 109 * u ** 2 + 224 * u)
    egg = mpmath.re(mpmath.quad(fq, [R1, R2]))
    win = mpmath.re(mpmath.quad(fq, [R1, W1]) + mpmath.quad(fq, [W2, mpmath.mpf(-14) / 3]))
    omega1 = None
    for ln in gp_out.splitlines():
        if ln.startswith("E.omega = ["):
            omega1 = mpmath.mpf(ln.split("[")[1].split(",")[0])
    check(omega1 is not None and abs(egg - omega1) < mpmath.mpf(10) ** (-12), "egg integral equals PARI omega1")
    print("egg integral int dx/sqrt(f) = %s = PARI omega1 to 1e-12 ; window integral = %s" % (mpmath.nstr(egg, 15), mpmath.nstr(win, 15)))
    print("arc-measure fraction of egg that is positive: %.5f (limit density of positive odd m in each coset row, Weyl equidistribution)" % float(win / egg))
    print("observed fraction among odd m<=80: k=0 row %d/40 = %.4f ; k=1 row %d/40 = %.4f ; both rows %d/80 = %.4f" % (
        len(pos_m), len(pos_m) / 40, len(pos_list) // 3 - len(pos_m), (len(pos_list) // 3 - len(pos_m)) / 40, len(pos_list) // 3, len(pos_list) / 3 / 80))

    print("\n--- S3: structural statement ---")
    print("E(Q) = Z*G (+) Z/6*T  [rank 1 PROVED by ellrank bounds; G generates mod torsion: saturation at p<=1000 FINITE-EXACT, full CITED Bremner-Macleod Rem 2.2]")
    print("positive solutions among {mG+kT : 1<=m<=80, 0<=k<=5} = sparse subset listed above; no Collatz map used.")
    print("cited canon/session labels: THM-3341, THM-3756, THM-3620, THM-4139, THM-4146; session collatz-mod6-20260917; wave 20260921")
    print("=== END ===")


if __name__ == "__main__":
    main()
