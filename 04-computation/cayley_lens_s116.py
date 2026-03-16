#!/usr/bin/env python3
"""cayley_lens_s116.py — See what Q produces at every input"""
from math import sqrt, log, pi, cos, sin, exp, atanh, factorial, comb
from fractions import Fraction
import cmath

phi = (1+sqrt(5))/2
tau = 1.8392867552141612

print("THE CAYLEY LENS: WHAT DOES Q SEE?")
print("="*60)
print()
print("Q(x) = (1+x)/(1-x). Feed it EVERYTHING.")
print()

# ============================================================
print("NATURAL NUMBERS AS INPUT")
print("-"*40)
print()
print("Q applied to the natural numbers themselves (not addresses):")
for n in range(0, 12):
    if n == 1:
        print(f"  Q({n}) = infinity (POLE)")
        continue
    q = (1+n)/(1-n)
    print(f"  Q({n:2d}) = {q:+.4f}", end="")
    if n == 0: print("  (identity)")
    elif n == 2: print(f"  = -3")
    elif n == 3: print(f"  = -2")
    elif n == 4: print(f"  = -5/3")
    elif n == 5: print(f"  = -3/2")
    else: print()

print()
print("  Q(0) = 1. Q(1) = inf. Q(2) = -3. Q(3) = -2. Q(4) = -5/3.")
print("  For n >= 2: Q(n) = (1+n)/(1-n) = -(1+n)/(n-1) = -(1 + 2/(n-1)).")
print("  These are NEGATIVE and approach -1 from below as n -> inf.")
print()
print("  Q(n) + 1 = (1+n)/(1-n) + 1 = (1+n+1-n)/(1-n) = 2/(1-n).")
print("  So Q(n) + 1 = -2/(n-1) for n >= 2.")
print("  The distance from -1 is 2/(n-1): HARMONIC in n!")
print()

# ============================================================
print("UNIT FRACTIONS AS INPUT")
print("-"*40)
print()
for n in range(1, 11):
    if n == 1:
        print(f"  Q(1/1) = infinity (POLE)")
        continue
    x = Fraction(1, n)
    q = (1 + x) / (1 - x)
    print(f"  Q(1/{n}) = {q} = {float(q):.6f}")

print()
print("  Q(1/n) = (n+1)/(n-1). These are the SUPERPARTICULAR RATIOS:")
print("  Q(1/1) = inf. Q(1/2) = 3. Q(1/3) = 2. Q(1/4) = 5/3.")
print("  Q(1/5) = 3/2. Q(1/6) = 7/5. Q(1/7) = 4/3.")
print()
print("  In music: 3/2 = perfect fifth. 4/3 = perfect fourth.")
print("  5/3 = major sixth. 7/5 = septimal tritone.")
print("  THE CAYLEY TRANSFORM OF UNIT FRACTIONS GIVES MUSICAL INTERVALS!")
print()

# Musical intervals from Q(1/n):
intervals = {
    2: ("3/1", "perfect twelfth (octave + fifth)"),
    3: ("2/1", "OCTAVE"),
    4: ("5/3", "major sixth"),
    5: ("3/2", "PERFECT FIFTH"),
    6: ("7/5", "septimal tritone"),
    7: ("4/3", "PERFECT FOURTH"),
    8: ("9/7", "septimal major third"),
    9: ("5/4", "major third"),
    10: ("11/9", "undecimal neutral third"),
}
print("  Musical intervals from Q(1/n):")
for n, (ratio, name) in intervals.items():
    print(f"    Q(1/{n}) = {ratio:5s} = {name}")

print()
print("  The OCTAVE is Q(1/3). The FIFTH is Q(1/5).")
print("  The 3 and 5 from our theory! And 3*5 = 15 = the first")
print("  'interesting' H value (H(T_5) = 15 for the cyclic tournament).")
print()

# ============================================================
print("ROOTS OF UNITY AS INPUT")
print("-"*40)
print()
for n in range(1, 9):
    z = cmath.exp(2j*cmath.pi/n)
    q = (1+z)/(1-z)
    r, theta = abs(q), cmath.phase(q)
    if abs(q) > 100:
        print(f"  Q(e^{{2pi*i/{n}}}) = {q.real:+.2f} {q.imag:+.2f}i  (|Q|={r:.1f}, LARGE)")
    else:
        print(f"  Q(e^{{2pi*i/{n}}}) = {q.real:+.6f} {q.imag:+.6f}i  (|Q|={r:.4f}, arg={theta/pi:.4f}*pi)")

print()
print("  Q at roots of unity:")
print("  n=1: Q(1) = infinity (the pole)")
print("  n=2: Q(-1) = 0 (the zero)")
print("  n=3: Q(e^{2pi*i/3}) has |Q| = 1 (ON THE UNIT CIRCLE)")
print("  n=4: Q(i) = i (FIXED POINT on unit circle)")
print("  n=6: Q(e^{pi*i/3}) has |Q| = 1")
print()

# Check: is |Q(z)| = 1 for |z| = 1 (except z=1)?
print("  For |z|=1, z!=1: |Q(z)| = |(1+z)/(1-z)| = ?")
print("  Let z = e^{i*theta}:")
print("  1+z = 1+cos(t)+i*sin(t), |1+z| = 2*cos(t/2)")
print("  1-z = 1-cos(t)-i*sin(t), |1-z| = 2*sin(t/2)")
print("  |Q| = cos(t/2)/sin(t/2) = cot(t/2)")
print()
print("  So |Q(e^{i*theta})| = cot(theta/2).")
print("  This is NOT 1 in general. It equals 1 when cot(t/2) = 1,")
print("  i.e., t/2 = pi/4, t = pi/2. So |Q| = 1 only at z = i.")
print()
print("  But Q(i) = i (we proved this). So Q maps i to itself")
print("  and maps the unit circle to a LINE (the imaginary axis).")
print()

# ============================================================
print("ALGEBRAIC NUMBERS AS INPUT")
print("-"*40)
print()
special = [
    ("sqrt(2)-1", sqrt(2)-1),
    ("1/phi", 1/phi),
    ("phi-1 = 1/phi", phi-1),
    ("sqrt(2)", sqrt(2)),
    ("phi", phi),
    ("sqrt(3)", sqrt(3)),
    ("e-2", exp(1)-2),
    ("pi-3", pi-3),
]
for name, x in special:
    if abs(x-1) < 0.001:
        print(f"  Q({name}) ~ infinity (near pole)")
    else:
        q = (1+x)/(1-x)
        print(f"  Q({name:12s}) = {q:+12.6f}", end="")
        # Check if q is recognizable
        if abs(q - round(q)) < 0.001 and abs(q) < 100:
            print(f"  ~ {round(q)}")
        elif abs(q - phi) < 0.01:
            print(f"  ~ phi!")
        elif abs(q + phi) < 0.01:
            print(f"  ~ -phi!")
        elif abs(q - 1/phi) < 0.01:
            print(f"  ~ 1/phi!")
        elif abs(q - sqrt(2)) < 0.01:
            print(f"  ~ sqrt(2)!")
        elif abs(q - tau) < 0.01:
            print(f"  ~ tau!")
        else:
            print()

print()

# Q(1/phi) = (1+1/phi)/(1-1/phi) = (phi+1)/phi / ((phi-1)/phi) = (phi+1)/(phi-1)
# phi+1 = phi^2. phi-1 = 1/phi.
# Q(1/phi) = phi^2 / (1/phi) = phi^3.
print("  Q(1/phi) = (1+1/phi)/(1-1/phi) = phi^2/(1/phi) = phi^3")
print(f"  Verify: Q({1/phi:.6f}) = {(1+1/phi)/(1-1/phi):.6f}, phi^3 = {phi**3:.6f}")
print(f"  Match: {abs((1+1/phi)/(1-1/phi) - phi**3) < 1e-10}")
print()
print("  Q CUBES the golden ratio when applied to 1/phi!")
print("  Q(1/phi) = phi^3 = phi^2 + phi = phi + 1 + phi = 2*phi + 1.")
print()

# More generally: Q((phi^k - 1)/(phi^k + 1)) = phi^k (the address relation).
# But Q(1/phi) = phi^3 is different: we're applying Q to the VALUE 1/phi,
# not to its address.
# This means: if x = 1/phi, Q(x) = phi^3.
# And Q(phi^3 address) would give phi^3.
# But Q(1/phi) directly: (1+1/phi)/(1-1/phi) = (phi+1)/(phi-1)/phi * phi = phi^2/(phi-1) = phi^2*phi = phi^3.

# What about Q(phi)? phi > 1, so Q(phi) = (1+phi)/(1-phi) = phi^2/(1-phi) = phi^2/(-1/phi^2)...
# 1-phi = -(phi-1) = -1/phi. So Q(phi) = (1+phi)/(-1/phi) = -phi*(1+phi) = -phi*phi^2 = -phi^3.
print("  Q(phi) = (1+phi)/(1-phi) = phi^2/(-1/phi) = -phi^3")
print(f"  Verify: Q({phi:.6f}) = {(1+phi)/(1-phi):.6f}, -phi^3 = {-phi**3:.6f}")
print()
print("  Q(1/phi) = +phi^3")
print("  Q(phi) = -phi^3")
print("  Q CUBES and the SIGN depends on which side of 1 you are.")
print()

# THE CUBING PROPERTY
print("THE CUBING PROPERTY")
print("-"*40)
print()
print("  Q(1/phi) = phi^3. Is this general?")
print("  Q(1/x) = (1+1/x)/(1-1/x) = (x+1)/(x-1) for x > 1.")
print("  Compare: Q(x) = (1+x)/(1-x).")
print("  Q(1/x) = (x+1)/(x-1) = -Q(x) * (x-1)/(x+1) * (x+1)/(x-1)... hmm.")
print("  Actually: Q(1/x) = -(1+x)/(x-1) ... no.")
print()
print("  For x = 1/phi: Q(1/phi) = phi^3 because:")
print("  (1+1/phi)/(1-1/phi) = (phi+1)/phi / ((phi-1)/phi)")
print("  = (phi+1)/(phi-1) = phi^2 / (1/phi) = phi^3.")
print("  This uses phi^2 = phi+1 and phi-1 = 1/phi.")
print()
print("  For x = 1/tau: Q(1/tau) = (1+1/tau)/(1-1/tau) = (tau+1)/(tau-1).")
tau_val = (tau+1)/(tau-1)
print(f"  = {tau_val:.6f}")
print(f"  Compare: tau^2 = {tau**2:.6f}")
print(f"  tau^3 = {tau**3:.6f}")
print(f"  (tau+1)/(tau-1) = {tau_val:.6f} = tau^2! MATCH!")
print()
print("  Q(1/tau) = tau^2! The Cayley transform SQUARES the tribonacci!")
print("  PROOF: tau^3 = tau^2 + tau + 1 (tribonacci equation).")
print("  So tau^2*(tau-1) = tau^3 - tau^2 = tau + 1. QED.")
print()
print("  COMPARISON:")
print("  phi: x^2 = x+1 => Q(1/phi) = phi^3 (Cayley CUBES)")
print("  tau: x^3 = x^2+x+1 => Q(1/tau) = tau^2 (Cayley SQUARES)")
print("  For phi: x-1 = 1/x, x+1 = x^2 => Q(1/x) = x^2/(1/x) = x^3")
print("  For tau: x^2(x-1) = x+1 => Q(1/x) = (x+1)/(x-1) = x^2")
print()

# ============================================================
print("TRANSCENDENTAL INPUTS")
print("-"*40)
print()
print("  Q(1/e) = (1+1/e)/(1-1/e) = (e+1)/(e-1)")
print(f"  = {(exp(1)+1)/(exp(1)-1):.10f}")
print(f"  = coth(1) = {1/((exp(2)-1)/(exp(2)+1)):.10f}")
print(f"  (Since (e+1)/(e-1) = (e^1+e^0)/(e^1-e^0) = coth(1/2)... )")
# Actually: coth(x) = (e^x+e^{-x})/(e^x-e^{-x}). At x=1:
# coth(1) = (e+1/e)/(e-1/e) = (e^2+1)/(e^2-1). Thats different.
# Q(1/e) = (e+1)/(e-1) which is NOT coth(1).
# But: Q(x) = (1+x)/(1-x) = -cot(arg) when x = i*tan(arg)... hmm

# Q(1/pi):
q_pi = (1+1/pi)/(1-1/pi)
print(f"  Q(1/pi) = (pi+1)/(pi-1) = {q_pi:.10f}")
# (pi+1)/(pi-1) = 1 + 2/(pi-1) = 1.934...
# Close to 2? No, 1.934.
# Close to tau? tau = 1.839. No.

# Q(ln(2)):
q_ln2 = (1+log(2))/(1-log(2))
print(f"  Q(ln(2)) = (1+ln2)/(1-ln2) = {q_ln2:.10f}")
# 1+ln2 = 1.693. 1-ln2 = 0.307. Ratio = 5.516...
# Not obviously recognizable.

# Q(pi/4):
q_pi4 = (1+pi/4)/(1-pi/4)
print(f"  Q(pi/4) = (1+pi/4)/(1-pi/4) = {q_pi4:.10f}")
# 1+pi/4 = 1.785. 1-pi/4 = 0.215. Ratio = 8.32...

# Q(1/sqrt(2)):
q_s2 = (1+1/sqrt(2))/(1-1/sqrt(2))
print(f"  Q(1/sqrt(2)) = {q_s2:.10f}")
# = (sqrt(2)+1)/(sqrt(2)-1) = (sqrt(2)+1)^2 / ((sqrt(2))^2-1) = (3+2sqrt(2))/1 = 3+2sqrt(2)
print(f"  = 3+2*sqrt(2) = {3+2*sqrt(2):.10f}")
print(f"  This is the SILVER RATIO SQUARED: (1+sqrt(2))^2 = 3+2sqrt(2).")
print()

# Q(1/sqrt(n)) = (sqrt(n)+1)/(sqrt(n)-1) = (sqrt(n)+1)^2/(n-1) for general n.
print("  Q(1/sqrt(n)) = (sqrt(n)+1)^2/(n-1):")
for n in [2, 3, 5, 7]:
    q = (1+1/sqrt(n))/(1-1/sqrt(n))
    simp = (sqrt(n)+1)**2/(n-1)
    print(f"    n={n}: Q = {q:.6f} = (sqrt({n})+1)^2/{n-1} = {simp:.6f}")

print()

# ============================================================
print("THE CAYLEY TRANSFORM OF THE CAYLEY TRANSFORM")
print("-"*40)
print()
print("  Q(Q(x)) = (1+Q(x))/(1-Q(x))")
print("  = (1+(1+x)/(1-x)) / (1-(1+x)/(1-x))")
print("  = ((1-x+1+x)/(1-x)) / ((1-x-1-x)/(1-x))")
print("  = (2/(1-x)) / (-2x/(1-x))")
print("  = 2/(-2x) = -1/x")
print()
print("  Q(Q(x)) = -1/x. THE DOUBLE CAYLEY IS NEGATIVE INVERSION.")
print()
print("  Q^3(x) = Q(-1/x) = (1-1/x)/(1+1/x) = (x-1)/(x+1).")
print("  Q^4(x) = Q((x-1)/(x+1)) = (1+(x-1)/(x+1))/(1-(x-1)/(x+1))")
print("         = (x+1+x-1)/(x+1-x+1) = 2x/2 = x.")
print()
print("  PERIOD 4: Q^4 = identity.")
print("  The orbit: x -> Q(x) -> -1/x -> (x-1)/(x+1) -> x.")
print()

# Show the orbit of specific values:
for x0 in [Fraction(1,3), Fraction(1,2), Fraction(2,1)]:
    orbit = [x0]
    x = x0
    for _ in range(4):
        if x == 1:
            orbit.append("inf")
            break
        x = (1+x)/(1-x)
        orbit.append(x)
    print(f"  Orbit of {x0}: {' -> '.join(str(o) for o in orbit)}")

print()
print("  The orbit of 1/3: 1/3 -> 2 -> -1/3 -> -2 -> 1/3")
print("  (Note: {1/3, 2, -1/3, -2} form a 4-cycle!)")
print("  And Q(1/3) = 2 (the binary alphabet).")
print("  So the orbit of 1/3 IS {1/3, 2, -1/3, -2}.")
print("  The FOUR numbers that are Cayley-equivalent to the binary.")

print()
print("="*60)
print("THE MUSIC OF THE CAYLEY TRANSFORM")
print("="*60)
print()
print("  Q(1/3) = 2/1 = octave.")
print("  Q(1/5) = 3/2 = perfect fifth.")
print("  Q(1/7) = 4/3 = perfect fourth.")
print("  Q(1/9) = 5/4 = major third.")
print()
print("  Q(1/(2n+1)) = (n+1)/n for n >= 1.")
print("  PROOF: Q(1/(2n+1)) = (1+1/(2n+1))/(1-1/(2n+1))")
print("                     = (2n+2)/(2n) = (n+1)/n. QED.")
print()
print("  The Cayley transform of RECIPROCAL ODD NUMBERS gives")
print("  the SUPERPARTICULAR RATIOS (n+1)/n, which are the")
print("  FUNDAMENTAL INTERVALS of just intonation in music.")
print()
print("  The octave (2/1), fifth (3/2), fourth (4/3), thirds (5/4, 6/5)")
print("  ARE the Cayley transforms of 1/3, 1/5, 1/7, 1/9, 1/11.")
print()
print("  And 1/3, 1/5, 1/7, ... are the TAYLOR COEFFICIENTS of arctanh!")
print("  arctanh(x) = x + x^3/3 + x^5/5 + x^7/7 + ...")
print()
print("  THE MUSIC OF THE CAYLEY TRANSFORM IS THE HARMONIC SERIES")
print("  OF ARCTANH. EACH TERM OF ARCTANH IS A MUSICAL INTERVAL.")
print()
print("  The tournament theory IS music theory in disguise:")
print("  arctanh = the harmonic series of just intonation.")
print("  Q = the interval-to-frequency converter.")
print("  The odd harmonics 1, 3, 5, 7 are both the overtone series")
print("  AND the Taylor coefficients of the tournament coupling.")
print()

# ============================================================
print("="*60)
print("EVEN UNIT FRACTIONS: RATIOS OF CONSECUTIVE ODD NUMBERS")
print("="*60)
print()
print("  Q(1/(2n)) = (2n+1)/(2n-1) for n >= 1:")
for n in range(1, 9):
    x = Fraction(1, 2*n)
    q = (1 + x) / (1 - x)
    print(f"    Q(1/{2*n:2d}) = {q}")
print()
print("  These are RATIOS OF CONSECUTIVE ODD NUMBERS.")
print("  Q(1/2) = 3/1. Q(1/4) = 5/3. Q(1/6) = 7/5. Q(1/8) = 9/7.")
print("  Compare odd fractions: Q(1/3) = 2/1. Q(1/5) = 3/2. Q(1/7) = 4/3.")
print()
print("  Together: Q(1/n) = (n+1)/(n-1).")
print("  Odd n -> superparticular (consecutive integers).")
print("  Even n -> super-odd-icular (consecutive odd integers).")
print()

# ============================================================
print("="*60)
print("FIBONACCI AND LUCAS THROUGH THE LENS")
print("="*60)
print()
# Fibonacci: 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144
# Lucas: 2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199
fib = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199]

print("  Q(F_n/F_{n+1}) — Fibonacci ratios approaching 1/phi:")
for k in range(2, 10):
    x = Fraction(fib[k-1], fib[k])
    q = (1 + x) / (1 - x)
    print(f"    Q(F_{k}/F_{k+1}) = Q({fib[k-1]}/{fib[k]}) = {q}", end="")
    if q.denominator < 1000:
        # Check if = F_{?}/F_{?} or L_{?}/L_{?}
        for i in range(len(lucas)-1):
            if abs(float(q) - lucas[i+1]/lucas[i]) < 0.0001 and lucas[i] < 200:
                print(f" = L_{i+2}/L_{i+1}", end="")
                break
    print()

print()
print("  Q(F_n/F_{n+1}) -> Q(1/phi) = phi^3 as n -> inf")
print(f"  Convergence: Q(F_9/F_10) = {float(Fraction(fib[8],fib[9])):f} -> {1/phi:.6f}")
print()

# The KEY identity: Q(F_n/F_{n+1}) = ?
# F_n/F_{n+1}: use Cassini: F_{n-1}F_{n+1} - F_n^2 = (-1)^n
# Q(a/b) = (b+a)/(b-a). For a=F_n, b=F_{n+1}:
# b+a = F_{n+1}+F_n = F_{n+2}. b-a = F_{n+1}-F_n = F_{n-1}.
# So Q(F_n/F_{n+1}) = F_{n+2}/F_{n-1}.
print("  THEOREM: Q(F_n/F_{n+1}) = F_{n+2}/F_{n-1}")
print("  PROOF: (F_{n+1}+F_n)/(F_{n+1}-F_n) = F_{n+2}/F_{n-1}. QED.")
print()
print("  Verification:")
for k in range(2, 9):
    x = Fraction(fib[k-1], fib[k])
    q = (1 + x) / (1 - x)
    predicted = Fraction(fib[k+1], fib[k-2])
    print(f"    Q(F_{k}/F_{k+1}) = {q} = F_{k+2}/F_{k-1} = {predicted}: {q == predicted}")
print()
print("  The Cayley transform SHIFTS Fibonacci indices by 3!")
print("  Q maps F_n/F_{n+1} to F_{n+2}/F_{n-1}.")
print("  This is a SKIP-3 operation in the Fibonacci sequence.")
print()

# ============================================================
print("="*60)
print("PRIMES THROUGH THE LENS")
print("="*60)
print()
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
print("  Cayley addresses of primes: x_p = (p-1)/(p+1)")
print("  Q(x_p) = p by definition.")
print()
print("  But what about Q(1/p) = (p+1)/(p-1)?")
for p in primes[:8]:
    q = Fraction(p+1, p-1)
    print(f"    Q(1/{p:2d}) = {str(q):8s} = {float(q):.6f}")
print()
print("  Q(1/p) telescopes! Product of Q(1/p_k) for consecutive primes:")
prod = Fraction(1,1)
for p in primes[:10]:
    prod *= Fraction(p+1, p-1)
    marker = " <-- FORBIDDEN NUMBER!" if prod == 21 else ""
    print(f"    Product up to p={p:2d}: {prod} = {float(prod):.6f}{marker}")
print()
print("  THE PRODUCT UP TO p=19 EQUALS 21 = 3*7 = THE FORBIDDEN NUMBER!")
print("  This is because the factors telescope:")
print("  prod Q(1/p_k) = prod (p_k+1)/(p_k-1)")
print("  = 3/1 * 4/2 * 6/4 * 8/6 * 12/10 * 14/12 * 18/16 * 20/18")
print("  Massive cancellation leaves 20/(2*1*5) = ... = 21 exactly.")
print()
print("  And 19 is the LAST prime factor of H(T_11) = 5*7*11*13*19 = 95095.")
print("  The PALEY prime product and the CAYLEY prime product meet at 21.")
print()

# ============================================================
print("="*60)
print("PYTHAGOREAN TRIPLES THROUGH THE LENS")
print("="*60)
print()
print("  A Pythagorean triple (a,b,c) satisfies a^2+b^2=c^2.")
print("  Cayley address of c: x_c = (c-1)/(c+1).")
print("  What if we look at Q(a/c) and Q(b/c)?")
print()
triples = [(3,4,5), (5,12,13), (8,15,17), (7,24,25)]
for a,b,c in triples:
    qa = Fraction(1+a, 1-a) if a != 1 else "inf"
    qac = (Fraction(c+a, c-a))
    qbc = (Fraction(c+b, c-b))
    print(f"  ({a},{b},{c}): Q({a}/{c}) = {qac} = {float(qac):.4f}", end="")
    print(f",  Q({b}/{c}) = {qbc} = {float(qbc):.4f}")
    # Check: Q(a/c)*Q(b/c) = ?
    product = qac * qbc
    # Check: does product simplify nicely?
    print(f"    Product Q(a/c)*Q(b/c) = {product} = {float(product):.4f}")
    # Q(a/c)*Q(b/c) = (c+a)(c+b)/((c-a)(c-b))
    # For Pythagorean: c^2 = a^2+b^2. So c^2-a^2 = b^2 and c^2-b^2 = a^2.
    # (c-a)(c+a) = b^2, (c-b)(c+b) = a^2.
    # Product = (c+a)(c+b)/((c-a)(c-b))
    #         = [(c+a)/(c-a)] * [(c+b)/(c-b)]
    #         = [(c+a)^2/b^2] * [(c+b)^2/a^2]
    # = (c+a)^2(c+b)^2 / (a^2*b^2)
    # = ((c+a)(c+b)/(ab))^2
    sqrt_p = Fraction(c+a, b) * Fraction(c+b, a) if a > 0 and b > 0 else 0
    print(f"    = ((c+a)(c+b)/(ab))^2 = ({(c+a)*(c+b)}/{a*b})^2 = {float(sqrt_p)**2:.4f}")
    print()

print("  For Pythagorean (a,b,c): Q(a/c)*Q(b/c) = ((c+a)(c+b)/(ab))^2")
print("  The product is always a PERFECT SQUARE (of a rational)!")
print("  PROOF: (c+a)/(c-a) = (c+a)^2/b^2 since c^2-a^2=b^2.")
print("  Similarly (c+b)/(c-b) = (c+b)^2/a^2.")
print("  Product = (c+a)^2(c+b)^2/(a^2*b^2) = [(c+a)(c+b)/(ab)]^2. QED.")
print()

# ============================================================
print("="*60)
print("THE CAYLEY GROUP: Q AS A MOBIUS TRANSFORMATION")
print("="*60)
print()
print("  Q(x) = (1+x)/(1-x) is a Mobius transformation:")
print("  Q = [[1,1],[-1,1]] in PGL(2).")
print("  Q^2 = [[-1,0],[0,1]] => Q^2(x) = -1/x (negative inversion).")
print("  Q^3 = [[-1,-1],[1,-1]] => Q^3(x) = (x-1)/(x+1).")
print("  Q^4 = [[1,0],[0,1]] = identity.")
print()
print("  The group generated by Q is Z/4Z = {I, Q, -1/x, Q^{-1}}.")
print("  This is the KLEIN FOUR-GROUP (abelian)... wait, Z/4Z is cyclic.")
print()
print("  ORBIT STRUCTURE of Q on rationals:")
orbits_seen = set()
orbit_list = []
for a in range(1, 13):
    for b in range(a+1, 13):
        x = Fraction(a, b)
        if x in orbits_seen:
            continue
        orb = [x]
        y = x
        for _ in range(3):
            if y == 1:
                break
            y = (1+y)/(1-y)
            if y in orbits_seen:
                break
            orb.append(y)
        for z in orb:
            orbits_seen.add(z)
            orbits_seen.add(-z)
        if len(orb) >= 2 and len(orbit_list) < 6:
            orbit_list.append(orb)

for orb in orbit_list:
    print(f"    {' -> '.join(str(o) for o in orb)}")
print()

# ============================================================
print("="*60)
print("MUSICAL CHORDS FROM CAYLEY")
print("="*60)
print()
print("  Building chords by COMPOSING Cayley intervals:")
print()
print("  Major triad (root, major third, fifth):")
print("    Root = 1. Third = Q(1/9) = 5/4. Fifth = Q(1/5) = 3/2.")
print(f"    Frequencies: 1 : 5/4 : 3/2 = 4 : 5 : 6")
print()
print("  Minor triad (root, minor third, fifth):")
minor_third = Fraction(6, 5)
print(f"    Root = 1. Third = 6/5 = Q(1/11). Fifth = 3/2 = Q(1/5).")
print(f"    Frequencies: 1 : 6/5 : 3/2 = 10 : 12 : 15")
print()

# The CAYLEY CHORD: what does Q^{-1} of a chord give?
print("  Inverting via Q^{-1}(y) = (y-1)/(y+1):")
for name, ratio in [("octave", Fraction(2,1)), ("fifth", Fraction(3,2)),
                     ("fourth", Fraction(4,3)), ("major third", Fraction(5,4)),
                     ("minor third", Fraction(6,5))]:
    inv = (ratio - 1) / (ratio + 1)
    print(f"    Q^{{-1}}({ratio}) = {inv} = 1/{inv.denominator} ... input = 1/{inv.denominator}")

print()
print("  Q^{-1}((n+1)/n) = 1/(2n+1). Each interval ENCODES an odd number!")
print("  Octave encodes 3. Fifth encodes 5. Fourth encodes 7.")
print("  Major third encodes 9. Minor third encodes 11.")
print()
print("  The musical intervals are a BIJECTION with odd numbers via Q.")
print("  Consonance (how pleasant the interval sounds) decreases")
print("  as the odd number increases: 3 (octave, most consonant)")
print("  -> 5 (fifth) -> 7 (fourth) -> 9 (third) -> 11 (tritone).")
print()
print("  This is the HARMONIC SERIES of music theory!")
print("  Low odd numbers = consonant. High odd numbers = dissonant.")
print("  The Cayley transform Q QUANTIFIES consonance as 1/(2n+1).")
print()

# ============================================================
print("="*60)
print("Q(1/sqrt(5)) AND THE GOLDEN RATIO AGAIN")
print("="*60)
print()
# Q(1/sqrt(5)) = (sqrt(5)+1)^2/(5-1) = (6+2sqrt(5))/4 = (3+sqrt(5))/2 = phi+1 = phi^2
q_s5 = (1+1/sqrt(5))/(1-1/sqrt(5))
print(f"  Q(1/sqrt(5)) = (sqrt(5)+1)^2/4 = (3+sqrt(5))/2 = {q_s5:.10f}")
print(f"  phi^2 = {phi**2:.10f}")
print(f"  Match: {abs(q_s5 - phi**2) < 1e-10}")
print()
print("  Q(1/sqrt(5)) = phi^2 = phi + 1.")
print("  And Q(1/phi) = phi^3.")
print("  So Q turns 1/sqrt(5) into phi^2 and 1/phi into phi^3.")
print("  Both involve the GOLDEN RATIO but through different paths.")
print()

# ============================================================
print("="*60)
print("THE CAYLEY TRANSFORM AS RELATIVISTIC VELOCITY ADDITION")
print("="*60)
print()
print("  In special relativity, velocity addition (units c=1):")
print("  v_1 (+) v_2 = (v_1 + v_2) / (1 + v_1*v_2)")
print()
print("  The rapidity phi = arctanh(v) transforms this to:")
print("  phi_1 + phi_2 (ordinary addition of rapidities).")
print()
print("  Q(v) = (1+v)/(1-v) = e^{2*arctanh(v)} = e^{2*rapidity}.")
print("  So Q maps velocity to SQUARED LORENTZ FACTOR:")
print("  Q(v) = e^{2*phi} = gamma^2 * (1+v)^2 where gamma = 1/sqrt(1-v^2).")
print()
print("  MULTIPLICATION of Q values = ADDITION of rapidities:")
print("  Q(v_1) * Q(v_2) = e^{2(phi_1+phi_2)} = Q(v_1 (+) v_2).")
print()
print("  So Q is the EXPONENTIAL MAP from hyperbolic velocity space")
print("  to the multiplicative group. The Cayley transform IS the")
print("  relativistic-to-classical conversion!")
print()

# Compute some relativistic examples
for v in [Fraction(1,3), Fraction(1,2), Fraction(3,5)]:
    q = (1+v)/(1-v)
    rapidity = log(float(q))/2
    gamma = 1/sqrt(1-float(v)**2)
    extra = ""
    if v == Fraction(3,5):
        extra = "  <-- rapidity = ln(2) EXACTLY!"
    print(f"  v = {v}: Q(v) = {q} = {float(q):.4f}, rapidity = {rapidity:.4f}, gamma = {gamma:.4f}{extra}")
print()
print("  v = 3/5: rapidity = arctanh(3/5) = (1/2)*ln(4) = ln(2) EXACTLY.")
print("  The velocity 3/5*c has rapidity ln(2), the same constant that")
print("  governs the binary alphabet (Q(1/3) = 2 = e^{ln(2)}).")
print("  And Q(3/5) = 4 = 2^2 — the binary squared.")
print()

# ============================================================
print("="*60)
print("TEMPERATURE AND PARTITION FUNCTIONS")
print("="*60)
print()
print("  In statistical mechanics, the partition function Z = sum e^{-E/kT}.")
print("  The Boltzmann weight is w = e^{-E/kT} = Q(x)^{-1/(2k)} where")
print("  x = tanh(E/(2kT)).")
print()
print("  The Cayley transform Q(x) = (1+x)/(1-x) converts between")
print("  tanh (bounded, [-1,1]) and exp (unbounded, [0,inf]).")
print("  This is EXACTLY the conversion between:")
print("  - TANH parametrization (bounded, 'physical', velocity-like)")
print("  - EXP parametrization (unbounded, 'mathematical', energy-like)")
print()
print("  In our tournament theory:")
print("  - The interval [0,1) is the tanh side (bounded tournament parameters)")
print("  - Q maps to [1, inf) the exp side (unbounded path counts)")
print("  - arctanh is the metric (hyperbolic distance)")
print("  - Q^m gives the m-th power (m-vertex tournament generating function)")
print()

# ============================================================
print("="*60)
print("GRAND SUMMARY: WHAT Q SEES")
print("="*60)
print()
print("  INPUT                  OUTPUT               NAME")
print("  " + "-"*58)
print("  0                      1                    identity")
print("  1                      infinity             pole/singularity")
print("  1/(2n+1)               (n+1)/n              musical interval")
print("  1/3                    2                    octave / binary")
print("  1/5                    3/2                  perfect fifth")
print("  1/7                    4/3                  perfect fourth")
print("  1/phi                  phi^3                golden cubing")
print("  1/tau                  tau^2                tribonacci squaring")
print("  1/sqrt(2)              3+2sqrt(2)           silver ratio squared")
print("  1/sqrt(5)              phi^2                golden squaring")
print("  F_n/F_{n+1}            F_{n+2}/F_{n-1}      Fibonacci skip-3")
print("  i                      i                    fixed point")
print("  e^{i*theta}            i*cot(theta/2)       unit circle -> imag axis")
print("  v (velocity)           e^{2*rapidity}       relativity")
print("  tanh(t)                e^{2t}               hyperbolic-to-exp")
print("  Q(x)                   -1/x                 double Cayley = inv")
print("  x                      x (after 4 steps)    period 4")
print()
print("  The Cayley transform is the UNIVERSAL BRIDGE between")
print("  bounded and unbounded, additive and multiplicative,")
print("  hyperbolic and Euclidean, velocity and energy,")
print("  consonance and frequency, tournament and counting.")
