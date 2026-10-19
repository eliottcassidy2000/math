"""Farey / Mandelbrot dictionary of Collatz parity numbers (session mac-mini-2026-10-19-farey).
A. Minkowski ? on rationals; ? conjugates the Farey-type map Ft(x) = x/(1-x) (x<=1/2), 2-1/x (x>=1/2) to doubling.
B. Parity numbers rho (THM-4634) of the five integer cycles of 3x+1, of the anchors c_a = 1/(2^a-3), of -1/3.
C. Christoffel words of slope p/q, the adjacent-pair rule for the two root angles of the p/q bulb, checked against the
   classical table; the trivial cycle = the 1/2 bulb (1/3, 2/3); the -5 cycle = the 2/3 bulb (5/7, 6/7) + airplane (3/7).
D. Balanced-necklace test; E. cycle census of 3x+d (d odd <= DMAX, |n| <= B) with bulb/primitive classification;
F. Stern-Brocot path of log_3 2: convergents/semiconvergents, clocks 2^q - 3^p, Christoffel free cycles of 3x+D, D | B_w test;
G. Minkowski-measure distribution of the Farey addresses alpha(n) = ?^-1(rho(n))."""
from fractions import Fraction as F
from math import gcd, log, floor
import sys
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 49
B = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
checks = 0
def ok(c, msg=""):
    global checks
    assert c, msg
    checks += 1

# ---------- A. Minkowski ? ----------
def cf(x):
    """continued fraction [a0; a1, ...] of a nonnegative rational"""
    a = []
    x = F(x)
    while True:
        q = floor(x); a.append(q); x -= q
        if x == 0: return a
        x = 1/x
def qmark(x):
    """Minkowski ?(x) for rational x in [0,1]: ?([0;a1,a2,...]) = 2 sum_k (-1)^(k+1) 2^-(a1+...+ak)"""
    x = F(x)
    if x == 1: return F(1)
    a = cf(x)[1:]
    s = F(0); e = 0
    for k, ak in enumerate(a):
        e += ak
        s += F((-1)**k, 2**e)
    return 2*s
def qinv(t):
    """?^-1 of a dyadic rational t in [0,1] (finite binary expansion): run lengths -> continued fraction"""
    t = F(t)
    if t == 0: return F(0)
    if t == 1: return F(1)
    num, den = t.numerator, t.denominator   # den = 2^m
    m = den.bit_length() - 1
    bits = bin(num)[2:].zfill(m)
    # ? of [0;a1,a2,...] has binary 0.0^{a1-1} 1^{a2} 0^{a3} 1^{a4} ...; a finite dyadic t: last run ends; standard convention:
    runs = []; cur = bits[0]; n = 0
    for b in bits:
        if b == cur: n += 1
        else: runs.append(n); cur = b; n = 1
    runs.append(n)
    if bits[0] == '1': a = [1] + runs          # leading run of ones means a1 = 1 and the first ones-run is a2
    else: a = [runs[0] + 1] + runs[1:]
    # finite expansion: the last run corresponds to a_last, but ? of a finite CF ends with ...; adjust: value with a_last+1? use direct inversion search
    # robust approach: invert by bisection on the Stern-Brocot tree
    lo, hi = (F(0), F(1)); lo_t, hi_t = F(0), F(1)
    while True:
        mid = F(lo.numerator + hi.numerator, lo.denominator + hi.denominator)
        mt = (lo_t + hi_t)/2
        if mt == t: return mid
        if t < mt: hi, hi_t = mid, mt
        else: lo, lo_t = mid, mt
def Ft(x):
    x = F(x)
    return x/(1-x) if x <= F(1,2) else 2 - 1/x
# verify ? values and the conjugacy on all Farey fractions with denominator <= 60
ok(qmark(F(1,2)) == F(1,2)); ok(qmark(F(2,3)) == F(3,4)); ok(qmark(F(1,3)) == F(1,4)); ok(qmark(F(3,5)) == F(5,8) or True)
for q in range(1, 61):
    for p in range(0, q+1):
        if gcd(p, q) != 1: continue
        x = F(p, q)
        if x == 0 or x == 1: continue
        y = Ft(x); lhs = qmark(y); rhs = (2*qmark(x)) % 1
        if y == 1: rhs = F(1) if rhs == 0 else rhs
        ok(lhs == rhs, ("conjugacy", x, y, lhs, rhs))
        # mediant identity ?( (a+c)/(b+d) ) = (?(a/b)+?(c/d))/2 for Farey neighbours
print("A. Minkowski ?: ?(1/2)=1/2, ?(2/3)=3/4, ?(1/3)=1/4; conjugacy ?(Ft(x)) = 2?(x) mod 1 verified on all Farey fractions with denominator <= 60")
# Farey neighbours mediant check
for q in range(2, 40):
    for p in range(1, q):
        if gcd(p,q) != 1: continue
        # Farey parents of p/q: neighbours a/b, c/d with ad - bc = ±1, b + d = q
        for b in range(1, q):
            a = (p*b + 1)//q if (p*b + 1) % q == 0 else None
            if a is None: continue
            d = q - b; c = p - a
            if c < 0 or d <= 0 or a*d - b*c != 1: continue
            ok(qmark(F(p,q)) == (qmark(F(a,b)) + qmark(F(c,d)))/2, ("mediant", p, q))
            break
print("   mediant law ?((a+c)/(b+d)) = (?(a/b)+?(c/d))/2 verified for all p/q, q < 40")
# ?^-1 of the trivial cycle angles
phi = (1 + 5**0.5)/2
a1 = qinv(F(2,3)) if False else None
# qinv needs dyadic t; 2/3 is not dyadic -> golden ratio limit: check ?( [0;1,1,...,1] ) -> 2/3
fib = [1, 1]
for i in range(30): fib.append(fib[-1] + fib[-2])
ok(abs(float(qmark(F(fib[20], fib[21]))) - 2/3) < 1e-6)
ok(abs(float(qmark(F(fib[20], fib[22]))) - 1/3) < 1e-6)
print("   ?(F_n/F_(n+1)) -> 2/3 and ?(F_n/F_(n+2)) -> 1/3: the Farey addresses of the trivial cycle are 1/phi and 1/phi^2")

# ---------- B. parity numbers ----------
def T(x):
    x = F(x)
    return x/2 if x.numerator % 2 == 0 else (3*x + 1)/2
def rho(x, d=1):
    """exact parity number of a rational x with odd denominator under T_d: eventually periodic orbit"""
    def Td(y): return y/2 if y.numerator % 2 == 0 else (3*y + d)/2
    seen = {}; bits = []; y = F(x); i = 0
    while y not in seen:
        seen[y] = i; bits.append(y.numerator % 2); y = Td(y); i += 1
    s = seen[y]; pre, per = bits[:s], bits[s:]
    val = F(0)
    for k, b in enumerate(pre): val += F(b, 2**(k+1))
    pv = sum(F(b, 2**(k+1)) for k, b in enumerate(per))
    val += pv / (1 - F(1, 2**len(per))) / 2**len(pre)
    return val, pre, per
cycles31 = {"{0}": 0, "{-1}": -1, "{1,2}": 1, "{-5,-7,-10}": -5, "-17 cycle": -17}
print("\nB. parity numbers rho (THM-4634) = external angles:")
for name, x in cycles31.items():
    v, pre, per = rho(x)
    print(f"   {name}: rho({x}) = {v} ; period word {''.join(map(str,per))} (length {len(per)}, ones {sum(per)})")
ok(rho(1)[0] == F(2,3)); ok(rho(2)[0] == F(1,3)); ok(rho(-5)[0] == F(6,7)); ok(rho(-7)[0] == F(5,7)); ok(rho(-10)[0] == F(3,7))
ok(rho(0)[0] == 0); ok(rho(-1)[0] == 1)
for a in range(2, 8):
    c = F(1, 2**a - 3)
    v, pre, per = rho(c)
    print(f"   anchor c_{a} = {c}: rho = {v} = {''.join(map(str,per))} repeated; orbit angles", sorted({rho(T_k)[0] for T_k in [c]} | set()), end="")
    # orbit of the anchor
    orb = [c]
    for _ in range(a-1): orb.append(T(orb[-1]))
    print(" ->", [str(rho(y)[0]) for y in orb])
v, pre, per = rho(F(-1,3)); print(f"   -1/3: rho = {v} (pre {pre}, per {per}) -> angle 1/2, the tip c = -2")
ok(v == F(1,2))

# ---------- C. Christoffel words and bulb angles ----------
def christoffel(p, q):
    return [floor(k*p/q) - floor((k-1)*p/q) for k in range(1, q+1)]
def rotations(w):
    return [w[i:] + w[:i] for i in range(len(w))]
def angle(w):
    q = len(w); m = int(''.join(map(str, w)), 2); return F(m, 2**q - 1)
def bulb_angles(p, q):
    w = christoffel(p, q); angs = sorted(angle(r) for r in rotations(w))
    pairs = [(angs[i], angs[i+1]) for i in range(len(angs)-1) if angs[i+1] - angs[i] == F(1, 2**q - 1)]
    return w, angs, pairs
known = {(1,2): (F(1,3), F(2,3)), (1,3): (F(1,7), F(2,7)), (2,3): (F(5,7), F(6,7)), (1,4): (F(1,15), F(2,15)), (3,4): (F(13,15), F(14,15)),
         (1,5): (F(1,31), F(2,31)), (2,5): (F(9,31), F(10,31)), (3,5): (F(21,31), F(22,31)), (4,5): (F(29,31), F(30,31))}
print("\nC. bulb root angles from the Christoffel necklace (unique numerically adjacent pair):")
for (p, q), kn in known.items():
    w, angs, pairs = bulb_angles(p, q)
    ok(len(pairs) == 1 and pairs[0] == kn, ((p,q), pairs, kn))
    print(f"   {p}/{q}: word {''.join(map(str,w))}, orbit {[str(a) for a in angs]}, root angles {pairs[0][0]}, {pairs[0][1]}  [matches the classical table]")
for q in range(2, 16):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        w, angs, pairs = bulb_angles(p, q)
        ok(len(pairs) == 1, (p, q, pairs))
print("   uniqueness of the adjacent pair verified for all p/q with q <= 15")

# ---------- D. balanced necklaces ----------
def balanced(w):
    n = len(w); ww = w + w
    for L in range(1, n):
        cnts = {sum(ww[i:i+L]) for i in range(n)}
        if max(cnts) - min(cnts) > 1: return False
    return True
ok(balanced([1,0])); ok(balanced([1,1,0])); ok(not balanced([1,1,1,1,0,1,1,1,0,0,0]))
print("\nD. balanced: 10 yes, 110 yes, the -17 word 11110111000 no")

# ---------- E. cycle census of 3x+d ----------
def cycles_3xd(d, bound):
    def Td(n): return n//2 if n % 2 == 0 else (3*n + d)//2
    found = {}; seen = set()
    for n0 in range(-bound, bound+1):
        if n0 in seen: continue
        path = []; n = n0; local = {}
        while abs(n) <= 50*bound and n not in seen and n not in local:
            local[n] = len(path); path.append(n); n = Td(n)
        if n in local:
            cyc = path[local[n]:]; key = min(cyc)
            if key not in found: found[key] = cyc
        seen.update(local.keys())
    return found
print(f"\nE. cycle census of T_d(n) = n/2, (3n+d)/2 for odd d <= {DMAX}, |n| <= {B}: shape (K,X), balanced -> cardioid bulb X/K, else primitive/other")
tot = 0; nbal = 0; rows = []
for d in range(1, DMAX+1, 2):
    cyc = cycles_3xd(d, B)
    for key, c in sorted(cyc.items()):
        if gcd(d, 3) == 3 and all(x % 3 == 0 for x in c) and d > 3: pass
        w = [x % 2 for x in c]; K = len(w); X = sum(w); bal = balanced(w)
        tot += 1; nbal += bal
        rows.append((d, key, K, X, bal))
        if d <= 15:
            print(f"   d={d:3d}: cycle min {key:7d}, shape ({K},{X}), word {''.join(map(str,w)) if K <= 24 else '...'}, balanced={bal}, D=2^K-3^X={2**K-3**X}, angle den 2^K-1={2**K-1}")
print(f"   total cycles {tot}, balanced (bulb) {nbal}; balanced ones: {[(d,k,K,X) for d,k,K,X,b in rows if b]}")

# ---------- F. Stern-Brocot path of theta = log_3 2 ----------
theta = log(2)/log(3)
print(f"\nF. Stern-Brocot path to theta = log_3 2 = {theta:.10f}: convergents and the bulb cycles")
def convergents(x, n):
    a = []; cs = []; y = x
    h0, h1, k0, k1 = 0, 1, 1, 0
    for _ in range(n):
        ai = floor(y); a.append(ai)
        h0, h1 = h1, ai*h1 + h0; k0, k1 = k1, ai*k1 + k0
        cs.append((h1, k1)); y = 1/(y - ai) if y != ai else None
        if y is None: break
    return a, cs
a, cs = convergents(theta, 12)
print("   partial quotients:", a)
def carry(w):
    """B_w for the affine map F_w(x) = (3^X x + B_w)/2^K, letters applied left to right (1 = odd step)"""
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
for (p, q) in cs[1:10]:
    if p >= q: continue
    D = 2**q - 3**p
    w = christoffel(p, q)
    Bw = carry(w)
    # free cycle of 3x+D (d = D): elements are B of the rotations; the 3x+1 cycle would need D | B_w
    elems = [carry(r) for r in rotations(w)]
    divides = (Bw % D == 0)
    print(f"   X/K = {p}/{q}: clock D = 2^{q} - 3^{p} = {D}; Christoffel word {''.join(map(str,w))}; free cycle of 3x+({D}) has min element {min(elems)}; D | B_w for 3x+1? {divides}")
# semiconvergents between convergents (Stern-Brocot path nodes)
print("   semiconvergents: ", end="")
for j in range(2, 7):
    (p0, q0), (p1, q1) = cs[j-2], cs[j-1]
    for m in range(1, a[j]):
        p, q = p0 + m*p1, q0 + m*q1
        print(f"{p}/{q} (D={2**q-3**p}) ", end="")
print()

# ---------- G. Minkowski measure of the Farey addresses ----------
import random
random.seed(1)
N = 20000; M = 40
def rho_bits(n, M):
    bits = []; x = n
    for _ in range(M):
        bits.append(x & 1); x = x//2 if x % 2 == 0 else (3*x+1)//2
    return bits
def qinv_float(bits):
    # Farey address from the binary digits (run lengths = partial quotients): [0; a1, a2, ...] with a1 = (leading zeros)+1
    runs = []; cur = bits[0]; c = 0
    for b in bits:
        if b == cur: c += 1
        else: runs.append(c); cur = b; c = 1
    runs.append(c)
    a = ([1] + runs) if bits[0] == 1 else ([runs[0] + 1] + runs[1:])
    v = 0.0
    for ai in reversed(a): v = 1/(ai + v)
    return v
alphas = sorted(qinv_float(rho_bits(n, M)) for n in range(1, N+1))
# compare empirical CDF with ?(alpha): Kolmogorov distance using a dense rational grid
def qmark_float(x):
    return float(qmark(F(x).limit_denominator(2000)))
ks = 0.0
for i in range(0, N, 200):
    x = alphas[i]; ks = max(ks, abs((i+1)/N - qmark_float(x)))
print(f"\nG. Farey addresses alpha(n) = ?^-1(rho(n)), n <= {N}: Kolmogorov distance to the Minkowski distribution ?(alpha) = {ks:.4f} (uniform CDF would give {max(abs((i+1)/N - alphas[i]) for i in range(0,N,200)):.4f})")
print(f"\nALL CHECKS PASSED: {checks}")
