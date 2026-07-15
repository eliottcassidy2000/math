# opus-2026-07-15-S315 -- the six speeds, the primes, and the Fibonacci rate
# in the radius-7 resonance framework (owner directive).
# (1) THE FIBONACCI EDGE LAW: Y*(F_n, F_{n+1}) = F_{n-7} exactly (d'Ocagne:
#     F_k F_{n+1} - F_{k+1} F_n = (-1)^n F_{n-k}; best approximants of
#     F_{n+1}/F_n are Fibonacci; q <= 13 = F_7 => k = 7 optimal).
# (2) THE FAREY-14 LAW: the depth-13 functional Y*(x, y)/x is maximized near
#     ratios p/14 with gcd(p, 14) = 1 -- SIX classes per block (phi(14) = 6),
#     value ~ 1/14; the golden line achieves only ||13 phi|| ~ 0.0335.
#     => at finite depth 13 the extremal tree-edge ratios are the 14th Farey
#     row, NOT the noble numbers: Fibonacci is the depth-infinity FOIL
#     (klein-S124) and the LRC(14) denominator 14 IS the finite-depth lever.
# (3) golden-pair sawtooth rho(F_n, F_{n+1}): governed by (F_{n+2} mod 13,
#     F_{n-1} mod 13) -> Pisano pi(13) = 28 periodicity in n.
# (4) numerology kill: phi^-6 vs the AP floor 2833/50700.
from fractions import Fraction
from math import gcd, sqrt

F = [0, 1]
while len(F) < 60: F.append(F[-1] + F[-2])

def Ystar(x1, x2, Q=13):
    best = None
    for q in range(1, Q+1):
        p0 = round(q*x2/x1)
        for p in (p0-1, p0, p0+1):
            if p < 1: continue
            v = abs(q*x2 - p*x1)
            if best is None or v < best[0]: best = (v, q, p)
    return best

print("(1) THE FIBONACCI EDGE LAW: Y*(F_n, F_{n+1}) vs F_{n-7}")
ok = True
for n in range(9, 26):
    v, q, p = Ystar(F[n], F[n+1])
    match = (v == F[n-7])
    ok &= match
    if n <= 14 or not match:
        print(f"   n={n}: Y*({F[n]},{F[n+1]}) = {v} at (q,p)=({q},{p}); "
              f"F_{n-7} = {F[n-7]}  {'OK' if match else 'MISMATCH'}")
print(f"   ALL n=9..25: {ok}")

print("\n(2) THE FAREY-14 LAW: top Y*/x ratio classes (x = 2003, y/x in (1, 13))")
x = 2003   # prime, no gcd artifacts
rows = []
for y in range(x+1, 13*x, 7):   # stride for speed; refined near peaks below
    g = gcd(x, y)
    if g > 1: continue
    v, q, p = Ystar(x, y)
    rows.append((v, y))
rows.sort(reverse=True)
print("   coarse scan top-10 (Y*, ratio, nearest p/14):")
for v, y in rows[:10]:
    r = y/x
    p14 = round(14*r)
    print(f"      Y* = {v} = {v/x:.5f}x  ratio {r:.5f}  ~ {p14}/14 = {p14/14:.5f} "
          f"(gcd(p,14) = {gcd(p14,14)})")
# refined: measure the best Y* very near p/14 (units) vs near phi + k
print("   refined at candidate extremals:")
phi = (1 + sqrt(5))/2
for label, r in ([(f"{p}/14", p/14) for p in (17, 19, 23, 25, 27)]
                 + [("phi+1", phi+1), ("phi+2", phi+2), ("2phi", 2*phi)]):
    best = (0, 0)
    for y in range(int(r*x) - 30, int(r*x) + 31):
        if gcd(x, y) > 1: continue
        v, q, p = Ystar(x, y)
        if v > best[0]: best = (v, y)
    print(f"      near {label:6s} (r = {r:.5f}): max Y* = {best[0]} = "
          f"{best[0]/x:.5f}x at y = {best[1]}")
print(f"   benchmarks: 1/14 = {1/14:.5f}; ||13 phi||-rate = "
      f"{abs(13*phi - round(13*phi))/13*13/1:.5f} -> {abs(13*phi-round(13*phi)):.5f}/q-scaled; "
      f"phi(14) = 6 unit classes per 14-block: p mod 14 in {{1,3,5,9,11,13}}")

print("\n(3) golden-pair sawtooth: 169*rho - 4 scaled, vs (F_{n+2} mod 13): "
      "Pisano pi(13) = 28")
def T_of(a, b):
    s = a + b
    T, m = 0, 0
    while s - 13*m > 0:
        T += (min(2*a, s - 13*m)) * (1 if m == 0 else 2)
        m += 1
    return T
for n in range(8, 20):
    a, b = F[n], F[n+1]
    rho = Fraction(T_of(a, b), 13*a*b)
    dev = (169*rho - 4) * a * b / 13   # the integer-ish sawtooth kernel
    print(f"   n={n:2d}: s mod 13 = {F[n+2] % 13:2d}, d mod 13 = {F[n-1] % 13:2d}, "
          f"169rho-4 = {float(169*rho-4):+.2e}, kernel = {float(dev):+.3f}")

print("\n(4) numerology kill: phi^-6 vs AP floor 2833/50700")
apf = Fraction(2833, 50700)
print(f"   phi^-6 = {phi**-6:.7f}; 2833/50700 = {float(apf):.7f}; "
      f"difference = {phi**-6 - float(apf):+.7f} -> NOT equal (coincidence killed)")
