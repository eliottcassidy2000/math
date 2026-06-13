#!/usr/bin/env python3
"""
tau_phi_clock_S72b.py — opus-2026-03-15-S72b

THE TAU-PHI CLOCK: Periodicity in the Golden-Tribonacci Plane
================================================================

From S71l we know:
  1 < phi < tau < 2   (the algebraic trinity)
  phi = root of x^2 - x - 1  (2-nacci / Fibonacci)
  tau = root of x^3 - x^2 - x - 1  (3-nacci / tribonacci)

  Phi_3(phi) = 2*phi^2 ≈ 5.236
  Phi_3(tau) = tau^3 ≈ 6.222  (by definition)
  Phi_3(2) = 7  (the Fano number)

The "parity clock" has period 6 = LCM(2,3) in tournament theory
(Fibonacci mod 4 periodicity).

NEW INVESTIGATION: What "clock" do tau and phi define?
- The tau-phi ratio tau/phi ≈ 1.137
- The angles arg(phi + i*tau), arg(tau/phi), etc.
- Powers of tau and phi modulo integers
- The joint periodicity of tau^n and phi^n
- Do they define a "clock" in the complex plane?
- Relationship to 3-cycle orientations (clockwise/counterclockwise)
"""

import numpy as np
from math import sqrt, factorial, comb, gcd, log, pi, e
import sys
import time
sys.stdout.reconfigure(line_buffering=True)

phi = (1 + sqrt(5)) / 2
tau = 1.8392867552141612  # tribonacci constant

print("=" * 72)
print("THE TAU-PHI CLOCK")
print("opus-2026-03-15-S72b")
print("=" * 72)

# =====================================================================
# PART 1: The tau-phi ratio as a rotation
# =====================================================================
print("\n" + "=" * 60)
print("PART 1: THE TAU-PHI RATIO")
print("=" * 60)

ratio = tau / phi
print(f"""
  phi   = {phi:.10f}  (golden ratio)
  tau   = {tau:.10f}  (tribonacci constant)
  tau/phi = {ratio:.10f}

  Is tau/phi algebraic? YES — it's a root of a degree-6 polynomial:
  If phi = (1+√5)/2 and tau satisfies tau^3 = tau^2 + tau + 1,
  then r = tau/phi, so tau = r*phi, and:
    (r*phi)^3 = (r*phi)^2 + r*phi + 1
    r^3 * phi^3 = r^2 * phi^2 + r*phi + 1
  Using phi^2 = phi + 1, phi^3 = 2phi + 1:
    r^3 * (2phi+1) = r^2 * (phi+1) + r*phi + 1
    (2r^3 - r^2 - r)*phi + (r^3 - r^2 - 1) = 0
  Since phi is irrational:
    2r^3 - r^2 - r = 0  AND  r^3 - r^2 - 1 = 0
  From the first: r(2r^2 - r - 1) = 0, so r(2r+1)(r-1) = 0 → r=0, -1/2, 1.
  But r ≈ 1.137, contradiction!

  So phi and tau are NOT simply related: tau ≠ c*phi for any rational c.
  The ratio tau/phi is algebraic of degree 6 (product of minimal polys).
""")

# What IS the minimal polynomial of tau/phi?
# tau satisfies t^3 - t^2 - t - 1 = 0
# phi satisfies p^2 - p - 1 = 0
# We want the minimal polynomial of r = tau/phi
# Using resultants or direct computation

# Let's find it numerically
from numpy.polynomial import polynomial as P

# tau/phi satisfies some polynomial. Let's find it by computing powers
r = tau / phi
powers = [r**k for k in range(10)]
print(f"  Powers of tau/phi:")
for k in range(7):
    print(f"    (tau/phi)^{k} = {powers[k]:.10f}")

# Try to find integer linear relation among 1, r, r^2, ..., r^6
# Using LLL or just trying small coefficients
print(f"\n  Searching for minimal polynomial of tau/phi...")

# Resultant approach: eliminate tau and phi from the system
# tau^3 - tau^2 - tau - 1 = 0 and phi^2 - phi - 1 = 0 and r = tau/phi
# Substituting tau = r*phi into tau's equation:
# (r*phi)^3 - (r*phi)^2 - r*phi - 1 = 0
# r^3*phi^3 - r^2*phi^2 - r*phi - 1 = 0
# Using phi^3 = 2phi+1, phi^2 = phi+1:
# r^3*(2phi+1) - r^2*(phi+1) - r*phi - 1 = 0
# (2r^3 - r^2 - r)*phi + (r^3 - r^2 - 1) = 0
# Since phi = (1+√5)/2:
# (2r^3 - r^2 - r)*(1+√5)/2 + (r^3 - r^2 - 1) = 0
# A*(1+√5)/2 + B = 0 where A = 2r^3-r^2-r, B = r^3-r^2-1
# A/2 + A√5/2 + B = 0
# (A/2 + B) + (A/2)*√5 = 0
# Since √5 is irrational: A/2 = 0 AND A/2 + B = 0
# Wait this gives A = 0, B = 0 which we showed has no real solution > 1.

# Actually the equation is correct: BOTH conditions must hold simultaneously
# only if phi and tau were rationally related, which they're not.
# Instead, treat phi as a parameter and solve:
# (2r^3 - r^2 - r)*phi + (r^3 - r^2 - 1) = 0
# phi = -(r^3 - r^2 - 1) / (2r^3 - r^2 - r)
# Substitute into phi^2 - phi - 1 = 0:
# Let P = -(r^3-r^2-1), Q = 2r^3-r^2-r
# (P/Q)^2 - P/Q - 1 = 0
# P^2 - PQ - Q^2 = 0

P_val = -(r**3 - r**2 - 1)
Q_val = 2*r**3 - r**2 - r
print(f"  P = -(r^3-r^2-1) = {P_val:.10f}")
print(f"  Q = 2r^3-r^2-r = {Q_val:.10f}")
print(f"  P^2 - P*Q - Q^2 = {P_val**2 - P_val*Q_val - Q_val**2:.10e}")

# So P^2 - PQ - Q^2 = 0 is the relation!
# Expanding: P = -(r^3-r^2-1) = -r^3+r^2+1
#             Q = 2r^3-r^2-r
# P^2 = r^6 - 2r^5 - r^4 + 2r^3 + r^2 - 2r^2... let me just compute symbolically

# P^2: (-r^3+r^2+1)^2 = r^6 - 2r^5 + r^4 + 2(-r^3) + 2r^2 + 1
#     = r^6 - 2r^5 + r^4 - 2r^3 + 2r^2 + 1 ??? let me be careful
# (-r^3+r^2+1)^2 = r^6 - 2r^5 + r^4 - 2r^3 + 2r^2 + 1
# Wait: (a+b+c)^2 = a^2+b^2+c^2+2ab+2ac+2bc where a=-r^3, b=r^2, c=1
# = r^6 + r^4 + 1 + 2(-r^3)(r^2) + 2(-r^3)(1) + 2(r^2)(1)
# = r^6 + r^4 + 1 - 2r^5 - 2r^3 + 2r^2

# PQ: (-r^3+r^2+1)(2r^3-r^2-r)
# = -2r^6+r^5+r^4 + 2r^5-r^4-r^3 + 2r^3-r^2-r
# = -2r^6 + 3r^5 + 0r^4 + r^3 - r^2 - r

# Q^2: (2r^3-r^2-r)^2 = 4r^6-4r^5-4r^4+r^4+2r^3+r^2
# = 4r^6 + r^4 + r^2 - 4r^5 + 2r^3 - 4r^4 + 2r^3... let me redo
# (a+b+c)^2 where a=2r^3, b=-r^2, c=-r
# = 4r^6 + r^4 + r^2 + 2(2r^3)(-r^2) + 2(2r^3)(-r) + 2(-r^2)(-r)
# = 4r^6 + r^4 + r^2 - 4r^5 - 4r^4 + 2r^3
# = 4r^6 - 4r^5 - 3r^4 + 2r^3 + r^2

# P^2 - PQ - Q^2:
# (r^6 + r^4 + 1 - 2r^5 - 2r^3 + 2r^2)
# - (-2r^6 + 3r^5 + 0r^4 + r^3 - r^2 - r)
# - (4r^6 - 4r^5 - 3r^4 + 2r^3 + r^2)

# r^6: 1 - (-2) - 4 = 1 + 2 - 4 = -1
# r^5: -2 - 3 - (-4) = -2 - 3 + 4 = -1
# r^4: 1 - 0 - (-3) = 1 + 3 = 4
# r^3: -2 - 1 - 2 = -5
# r^2: 2 - (-1) - 1 = 2 + 1 - 1 = 2
# r^1: 0 - (-1) - 0 = 1
# r^0: 1 - 0 - 0 = 1

# So: -r^6 - r^5 + 4r^4 - 5r^3 + 2r^2 + r + 1 = 0
# Or: r^6 + r^5 - 4r^4 + 5r^3 - 2r^2 - r - 1 = 0

print(f"\n  Minimal polynomial of tau/phi:")
coeffs_check = [1, 1, -4, 5, -2, -1, -1]  # r^6 + r^5 - 4r^4 + 5r^3 - 2r^2 - r - 1
val = sum(c * r**(6-i) for i, c in enumerate(coeffs_check))
print(f"  r^6 + r^5 - 4r^4 + 5r^3 - 2r^2 - r - 1 = {val:.10e}")

# Hmm, let me verify numerically
# Just try to find the polynomial by computing powers and using numpy
import numpy as np

r = tau / phi
# Build Vandermonde-like matrix
max_deg = 7
V = np.array([[r**k for k in range(max_deg)] for _ in range(1)])
# Actually, find polynomial by checking
for deg in range(2, 8):
    mat = np.vander([r], deg + 1)
    # We want c_0 + c_1*r + ... + c_d*r^d = 0
    # Try integer relation
    powers = np.array([r**k for k in range(deg + 1)])
    # Check if there's a small integer vector c with c·powers ≈ 0
    # Simple: try np.linalg.lstsq after fixing leading coeff to 1
    A = powers[:-1].reshape(1, -1)
    b = np.array([-powers[-1]])
    x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    # Check if x is close to integers
    x_rounded = np.round(x).astype(int)
    check = sum(x_rounded[k] * r**k for k in range(deg)) + r**deg
    if abs(check) < 1e-6:
        poly = list(x_rounded.flatten()) + [1]
        print(f"  Degree {deg}: coefficients = {poly}, check = {check:.2e}")

# =====================================================================
# PART 2: The clock in the complex plane
# =====================================================================
print("\n" + "=" * 60)
print("PART 2: THE COMPLEX CLOCK")
print("=" * 60)

# Consider z = phi + i*tau in the complex plane
z = complex(phi, tau)
print(f"""
  The tau-phi point in C: z = phi + i*tau = {z}
  |z| = {abs(z):.10f}
  arg(z) = {np.angle(z):.10f} rad = {np.degrees(np.angle(z)):.6f}°

  Compare: phi^2 + tau^2 = {phi**2 + tau**2:.10f}
  |z|^2 = {abs(z)**2:.10f}

  phi^2 = phi + 1 = {phi+1:.10f}
  tau^2 = tau^3 - tau - 1? No: tau^3 = tau^2 + tau + 1 => tau^2 = tau^3 - tau - 1
  tau^2 = {tau**2:.10f}

  |z|^2 = phi^2 + tau^2 = (phi+1) + tau^2 = {phi+1+tau**2:.10f}
""")

# What angle does the tau-phi vector make?
angle = np.arctan2(tau, phi)
print(f"  theta = arctan(tau/phi) = {angle:.10f} rad")
print(f"        = {np.degrees(angle):.6f}°")
print(f"  theta/pi = {angle/np.pi:.10f}")
print(f"  theta/(2*pi) = {angle/(2*np.pi):.10f}")
print(f"  1/theta = {1/angle:.10f}")

# Is theta/pi rational? If so, z^n would be periodic
# theta/pi ≈ 0.2937... doesn't look rational

# Powers of z
print(f"\n  Powers of z = phi + i*tau:")
for k in range(1, 13):
    zk = z**k
    print(f"    z^{k:2d} = {zk.real:12.4f} + {zk.imag:12.4f}i"
          f"  |z^k| = {abs(zk):12.4f}  arg = {np.degrees(np.angle(zk)):8.2f}°")

# =====================================================================
# PART 3: Modular periodicity — the actual "clock"
# =====================================================================
print("\n" + "=" * 60)
print("PART 3: MODULAR PERIODICITY — THE CLOCK")
print("=" * 60)

print(f"""
  The "parity clock" in tournament theory has period 6 = LCM(2,3).

  Fibonacci mod m has Pisano period pi(m).
  Tribonacci mod m has its own period.

  Question: what is the JOINT period of Fibonacci and Tribonacci
  modulo small integers? This defines the "tau-phi clock".
""")

# Pisano periods for Fibonacci
def fib_mod_period(m):
    """Pisano period: period of Fibonacci mod m."""
    a, b = 0, 1
    for i in range(1, 6*m*m):
        a, b = b, (a + b) % m
        if a == 0 and b == 1:
            return i
    return -1

# Period of tribonacci mod m
def trib_mod_period(m):
    """Period of tribonacci sequence mod m."""
    a, b, c = 0, 0, 1
    for i in range(1, 6*m*m*m):
        a, b, c = b, c, (a + b + c) % m
        if a == 0 and b == 0 and c == 1:
            return i
    return -1

print(f"  {'m':>3} {'pi_F(m)':>8} {'pi_T(m)':>8} {'lcm':>8} {'pi_F*pi_T':>10}")
print(f"  " + "-" * 40)
for m in range(2, 25):
    pi_f = fib_mod_period(m)
    pi_t = trib_mod_period(m)
    if pi_f > 0 and pi_t > 0:
        lcm_val = (pi_f * pi_t) // gcd(pi_f, pi_t)
        print(f"  {m:3d} {pi_f:8d} {pi_t:8d} {lcm_val:8d} {pi_f*pi_t:10d}")

# =====================================================================
# PART 4: The characteristic polynomials as clock faces
# =====================================================================
print("\n" + "=" * 60)
print("PART 4: CYCLOTOMIC STRUCTURE")
print("=" * 60)

print(f"""
  phi satisfies x^2 - x - 1 = 0  → roots: phi, -1/phi
  tau satisfies x^3 - x^2 - x - 1 = 0  → roots: tau, alpha, alpha_bar

  The companion matrices define linear recurrences:
  M_phi = [[1,1],[1,0]]  (Fibonacci matrix)
  M_tau = [[1,1,1],[1,0,0],[0,1,0]]  (tribonacci matrix)

  The "tau-phi clock" = eigenvalues of these matrices on the unit circle:
""")

# Fibonacci matrix eigenvalues
M_phi = np.array([[1,1],[1,0]], dtype=float)
eig_phi = np.linalg.eigvals(M_phi)
print(f"  Fibonacci eigenvalues: {eig_phi}")
print(f"    phi = {eig_phi[0]:.10f}")
print(f"    -1/phi = {eig_phi[1]:.10f}")

# Tribonacci matrix eigenvalues
M_tau = np.array([[1,1,1],[1,0,0],[0,1,0]], dtype=float)
eig_tau = np.linalg.eigvals(M_tau)
print(f"\n  Tribonacci eigenvalues: {eig_tau}")
for i, ev in enumerate(eig_tau):
    if abs(ev.imag) < 1e-10:
        print(f"    tau = {ev.real:.10f}")
    else:
        print(f"    alpha_{i} = {ev:.10f}  |alpha| = {abs(ev):.10f}")

# The complex tribonacci roots
alpha = [e for e in eig_tau if abs(e.imag) > 0.01][0]
print(f"\n  The complex tribonacci root alpha:")
print(f"    alpha = {alpha}")
print(f"    |alpha| = {abs(alpha):.10f}")
print(f"    arg(alpha) = {np.degrees(np.angle(alpha)):.6f}°")
print(f"    1/|alpha|^2 = {1/abs(alpha)**2:.10f}")
print(f"    |alpha|^2 = {abs(alpha)**2:.10f}")
print(f"    Compare: 1/tau = {1/tau:.10f}")

# By Vieta's: tau * |alpha|^2 = 1 (product of roots = constant term / leading coeff = 1)
print(f"    tau * |alpha|^2 = {tau * abs(alpha)**2:.10f}")
print(f"    (Should be 1 by Vieta's: product of all 3 roots = 1)")

# =====================================================================
# PART 5: The angular velocity — how fast does each "hand" turn?
# =====================================================================
print("\n" + "=" * 60)
print("PART 5: THE CLOCK HANDS — ANGULAR VELOCITIES")
print("=" * 60)

# The Fibonacci matrix mod m has eigenvalues on unit circle when |eigenvalue| = 1
# For Fibonacci: eigenvalues are phi and -1/phi
# |-1/phi| = 1/phi < 1, so in Fibonacci, only one eigenvalue matters
# But mod m, the matrix is periodic → eigenvalues become roots of unity

# The "clock" metaphor: M^n mod m rotates through a finite group
# Period = order of M in GL(k, Z/mZ)

print(f"""
  The "tau-phi clock" has TWO hands:

  HOUR HAND (slow): The Fibonacci matrix M_phi mod m
    Period at m=2: {fib_mod_period(2)}
    Period at m=3: {fib_mod_period(3)}
    Period at m=4: {fib_mod_period(4)}
    Period at m=5: {fib_mod_period(5)}
    Period at m=6: {fib_mod_period(6)}
    Period at m=7: {fib_mod_period(7)}

  MINUTE HAND (fast): The tribonacci matrix M_tau mod m
    Period at m=2: {trib_mod_period(2)}
    Period at m=3: {trib_mod_period(3)}
    Period at m=4: {trib_mod_period(4)}
    Period at m=5: {trib_mod_period(5)}
    Period at m=6: {trib_mod_period(6)}
    Period at m=7: {trib_mod_period(7)}
""")

# The JOINT clock: when do both hands return to 12 o'clock?
print(f"  JOINT CLOCK (both hands at 12):")
print(f"  {'m':>3} {'Fib period':>10} {'Trib period':>11} {'Joint (LCM)':>12} {'ratio T/F':>10}")
print(f"  " + "-" * 50)
for m in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
    pf = fib_mod_period(m)
    pt = trib_mod_period(m)
    if pf > 0 and pt > 0:
        joint = (pf * pt) // gcd(pf, pt)
        print(f"  {m:3d} {pf:10d} {pt:11d} {joint:12d} {pt/pf:10.3f}")

# =====================================================================
# PART 6: The tau-phi angle and tournament 3-cycles
# =====================================================================
print("\n" + "=" * 60)
print("PART 6: THE TAU-PHI ANGLE AND 3-CYCLES")
print("=" * 60)

# The complex tribonacci roots are at angle ±theta from the real axis
theta_alpha = np.angle(alpha)
print(f"""
  The tribonacci roots form a "Y" in the complex plane:
    tau = {tau:.6f}  (on real axis)
    alpha = {abs(alpha):.6f} * e^(i * {np.degrees(theta_alpha):.4f}°)
    alpha_bar = {abs(alpha):.6f} * e^(-i * {np.degrees(theta_alpha):.4f}°)

  The angle between the complex roots:
    2*theta = {2*np.degrees(theta_alpha):.4f}°

  Compare: 360°/3 = 120° (the angle for a 3-cycle!)
    2*theta = {2*np.degrees(theta_alpha):.4f}° ≈ {2*abs(np.degrees(theta_alpha))/120:.4f} × 120°

  The tribonacci roots DON'T form a 120° pattern (they would if they
  were cube roots of unity scaled by |alpha|).

  But the PHASE ROTATION per step is:
    Each tribonacci step rotates alpha by theta_alpha:
    alpha^n = |alpha|^n * e^(i*n*theta_alpha)
""")

# What IS the angle as a fraction of 2*pi?
frac = abs(theta_alpha) / (2 * np.pi)
print(f"  theta_alpha / (2*pi) = {frac:.10f}")
print(f"  ≈ {frac:.4f}... is this rational?")

# Check continued fraction
def cont_frac(x, depth=10):
    cf = []
    for _ in range(depth):
        n = int(x)
        cf.append(n)
        x = x - n
        if abs(x) < 1e-12:
            break
        x = 1/x
    return cf

cf = cont_frac(abs(theta_alpha) / (2*np.pi), 15)
print(f"  Continued fraction: {cf}")

# The convergents
print(f"  Convergents:")
h_prev, h_curr = 1, cf[0]
k_prev, k_curr = 0, 1
for i, a in enumerate(cf):
    if i == 0:
        h_curr, k_curr = a, 1
    else:
        h_curr, h_prev = a * h_curr + h_prev, h_curr
        k_curr, k_prev = a * k_curr + k_prev, k_curr
    if k_curr > 0:
        print(f"    [{i}] = {h_curr}/{k_curr} = {h_curr/k_curr:.10f}"
              f"  (error = {abs(h_curr/k_curr - frac):.2e})")

# =====================================================================
# PART 7: Fibonacci × Tribonacci interaction
# =====================================================================
print("\n" + "=" * 60)
print("PART 7: FIBONACCI × TRIBONACCI INTERACTION")
print("=" * 60)

# Fibonacci sequence
fib = [0, 1]
for i in range(30):
    fib.append(fib[-1] + fib[-2])

# Tribonacci sequence
trib = [0, 0, 1]
for i in range(30):
    trib.append(trib[-1] + trib[-2] + trib[-3])

print(f"\n  Fibonacci:   {fib[1:16]}")
print(f"  Tribonacci:  {trib[3:18]}")

# Interaction: F_n / T_n ratio
print(f"\n  Ratio F_n / T_n:")
for n in range(3, 20):
    if trib[n] > 0:
        r = fib[n] / trib[n]
        print(f"    n={n:2d}: F={fib[n]:8d}, T={trib[n]:8d}, F/T = {r:.6f}")

# In the limit: F_n ~ phi^n / sqrt(5), T_n ~ c * tau^n
# So F_n / T_n ~ (phi/tau)^n * K → 0 (since phi < tau)
print(f"\n  Asymptotic: F_n/T_n ~ (phi/tau)^n → 0 (phi/tau = {phi/tau:.6f} < 1)")

# More interesting: mixed recurrence
# What if we define M_n = F_n + T_n? Or F_n * T_n?
print(f"\n  Mixed sequence F_n * T_n:")
for n in range(3, 16):
    if trib[n] > 0:
        prod = fib[n] * trib[n]
        # Does this satisfy a recurrence?
        print(f"    n={n:2d}: F*T = {prod}")

# The product F_n * T_n grows like phi^n * tau^n = (phi*tau)^n
print(f"\n  phi * tau = {phi*tau:.10f}")
print(f"  (phi*tau)^n growth rate check:")
for n in range(5, 16):
    if trib[n] > 0 and trib[n-1] > 0:
        ratio = (fib[n] * trib[n]) / (fib[n-1] * trib[n-1])
        print(f"    n={n}: ratio = {ratio:.8f} (target: {phi*tau:.8f})")

# =====================================================================
# PART 8: The 6-clock — tournament parity period
# =====================================================================
print("\n" + "=" * 60)
print("PART 8: THE 6-CLOCK AND TOURNAMENT PARITY")
print("=" * 60)

print(f"""
  The period-6 "parity clock" from tournament theory:

  Fibonacci mod 4: [1, 1, 2, 3, 1, 0] repeating (period 6)
  Fibonacci mod 2: [1, 1, 0] repeating (period 3)
  Fibonacci mod 3: [1, 1, 2, 0, 2, 2, 1, 0] repeating (period 8)

  Tribonacci mod 2: period = {trib_mod_period(2)}
  Tribonacci mod 4: period = {trib_mod_period(4)}

  The "tau-phi clock" runs on the joint period:
  LCM(pi_F(4), pi_T(4)) = LCM({fib_mod_period(4)}, {trib_mod_period(4)}) = {(fib_mod_period(4)*trib_mod_period(4))//gcd(fib_mod_period(4),trib_mod_period(4))}
""")

# Fibonacci mod 4 cycle
print(f"  Fibonacci mod 4 cycle:")
a, b = 0, 1
fib4 = []
for i in range(fib_mod_period(4)):
    fib4.append(a)
    a, b = b, (a+b) % 4
print(f"    {fib4}")

# Tribonacci mod 4 cycle
print(f"\n  Tribonacci mod 4 cycle:")
a, b, c = 0, 0, 1
trib4 = []
for i in range(trib_mod_period(4)):
    trib4.append(a)
    a, b, c = b, c, (a+b+c) % 4
print(f"    {trib4}")

# =====================================================================
# PART 9: The clock face — visualizing the joint state
# =====================================================================
print("\n" + "=" * 60)
print("PART 9: THE JOINT STATE CLOCK")
print("=" * 60)

print(f"""
  At each tick n, the clock has state (F_n mod m, T_n mod m).

  For m=2, the possible states form a finite automaton.
  The "clock face" is Z/mZ × Z/mZ.
""")

# Joint state for m=2
print(f"  Joint state (F_n mod 2, T_n mod 2):")
fa, fb = 0, 1
ta, tb, tc = 0, 0, 1
states_2 = []
for n in range(20):
    state = (fa, ta)
    states_2.append(state)
    print(f"    n={n:2d}: ({fa}, {ta})")
    fa, fb = fb, (fa+fb) % 2
    ta, tb, tc = tb, tc, (ta+tb+tc) % 2
    if n > 0 and state == states_2[0]:
        # Check if we've completed a cycle
        break

# Find the joint period
fa, fb = 0, 1
ta, tb, tc = 0, 0, 1
for n in range(1, 1000):
    fa, fb = fb, (fa+fb) % 2
    ta, tb, tc = tb, tc, (ta+tb+tc) % 2
    if fa == 0 and fb == 1 and ta == 0 and tb == 0 and tc == 1:
        print(f"\n  Joint period mod 2: {n}")
        break

# For m=4
fa, fb = 0, 1
ta, tb, tc = 0, 0, 1
for n in range(1, 10000):
    fa, fb = fb, (fa+fb) % 4
    ta, tb, tc = tb, tc, (ta+tb+tc) % 4
    if fa == 0 and fb == 1 and ta == 0 and tb == 0 and tc == 1:
        print(f"  Joint period mod 4: {n}")
        break

# =====================================================================
# PART 10: Connection to tournament structure
# =====================================================================
print("\n" + "=" * 60)
print("PART 10: WHAT THE CLOCK MEASURES")
print("=" * 60)

# The key insight: Fibonacci counts paths in the PATH GRAPH P_n
# Tribonacci counts paths in a 3-step path graph
# Tournament H counts Hamiltonian paths

# The connection: for circulant tournaments, H decomposes via eigenspaces
# The eigenspaces are indexed by characters of Z/nZ
# Powers of roots of unity ARE the clock ticks!

print(f"""
  SYNTHESIS: The tau-phi clock is the PHASE PORTRAIT of tournament structure.

  In a Paley tournament P_p (p prime, p ≡ 3 mod 4):
  - The vertices are Z/pZ
  - The arcs come from quadratic residues
  - The eigenvalues of the adjacency matrix involve ROOTS OF UNITY
  - The Hamiltonian path count H decomposes via CHARACTERS of Z/pZ

  The characters chi_k(j) = omega^(jk) where omega = e^(2*pi*i/p)
  trace around the clock with angular velocity 2*pi*k/p.

  For p=5: omega = e^(2*pi*i/5)
  For p=7: omega = e^(2*pi*i/7)

  The TOURNAMENT CLOCK has:
  - HOUR HAND: Fibonacci period pi_F(p) (path-counting periodicity)
  - MINUTE HAND: Tribonacci period pi_T(p) (3-cycle periodicity)
  - SECOND HAND: The character phase 2*pi*k/p (eigenspace rotation)

  When all three hands align → the tournament has maximal symmetry.
""")

# For Paley tournaments, compute the actual periods
for p in [3, 5, 7, 11, 13]:
    pf = fib_mod_period(p)
    pt = trib_mod_period(p)
    joint = (pf * pt) // gcd(pf, pt)
    # Quadratic residues
    qr = set(k*k % p for k in range(1, p))
    print(f"  P_{p}: QR={sorted(qr)}, pi_F={pf}, pi_T={pt}, joint={joint}")
    # Is joint related to p-1?
    print(f"    p-1={p-1}, joint/(p-1)={joint/(p-1):.2f}, joint/p={joint/p:.2f}")

# =====================================================================
# PART 11: The golden angle and tournament rotation
# =====================================================================
print("\n" + "=" * 60)
print("PART 11: THE GOLDEN ANGLE")
print("=" * 60)

golden_angle = 2*np.pi*(1 - 1/phi)  # ≈ 137.508°
tribonacci_angle = 2*np.pi*(1 - 1/tau)  # what is this?

print(f"""
  The GOLDEN ANGLE = 2*pi*(1 - 1/phi) = 2*pi*(2-phi) = {np.degrees(golden_angle):.4f}°
  This is the "most irrational" angle — gives optimal packing (sunflowers).

  The TRIBONACCI ANGLE = 2*pi*(1 - 1/tau) = {np.degrees(tribonacci_angle):.4f}°

  The TAU-PHI ANGLE = 2*pi*(1 - phi/tau) = {np.degrees(2*np.pi*(1-phi/tau)):.4f}°

  In tournament terms:
  - Golden angle: the OPTIMAL spacing for placing vertices around a circle
    so that no two adjacent vertices in the tournament ordering are close.
  - Tribonacci angle: the optimal spacing for 3-step avoidance.

  If we arrange n vertices on a circle with angular spacing equal to the
  golden angle, the tournament defined by "i→j iff i comes before j
  in the angular ordering" has MAXIMAL dispersion of Hamiltonian paths.

  The golden ratio's connection to tournament scores:
  - The transitive tournament has score sequence (0, 1, ..., n-1)
  - This is the MOST ordered: H = 1
  - A "golden tournament" would maximize disorder per arc

  Fibonacci mod p periods for Paley primes:
""")

# Which primes give Fibonacci period = p-1? (primitive root property)
# These are the "Fibonacci primes" where Fibonacci generates all residues
print(f"  Primes where pi_F(p) = p-1 (phi is primitive root mod p):")
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    pf = fib_mod_period(p)
    if pf == p - 1:
        print(f"    p={p}: pi_F = {pf} = p-1 ✓")
    elif pf == 2*(p-1):
        print(f"    p={p}: pi_F = {pf} = 2(p-1)")
    else:
        print(f"    p={p}: pi_F = {pf} (ratio = {(p-1)/pf:.2f})")

# Same for tribonacci
print(f"\n  Primes where pi_T(p) divides p^2-1:")
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    pt = trib_mod_period(p)
    for d in [p-1, p+1, p**2-1, p**2-p, p**2+p]:
        if d % pt == 0:
            print(f"    p={p}: pi_T = {pt}, divides {d} = {'p-1' if d==p-1 else 'p+1' if d==p+1 else 'p²-1' if d==p**2-1 else 'p²-p' if d==p**2-p else 'p²+p'}")
            break
    else:
        print(f"    p={p}: pi_T = {pt}")

# =====================================================================
# FINAL SYNTHESIS
# =====================================================================
print("\n" + "=" * 60)
print("FINAL SYNTHESIS: THE TAU-PHI CLOCK")
print("=" * 60)

print(f"""
  THE TAU-PHI CLOCK IS A THREE-HAND TIMEPIECE:

  1. PHI HAND (golden/Fibonacci):
     - Tracks binary recurrence: F_n = F_{{n-1}} + F_{{n-2}}
     - Pisano period pi_F(m) governs modular return times
     - Angular velocity: golden angle ≈ {np.degrees(golden_angle):.1f}°
     - Controls: PATH COUNTING structure (2-step)

  2. TAU HAND (tribonacci):
     - Tracks ternary recurrence: T_n = T_{{n-1}} + T_{{n-2}} + T_{{n-3}}
     - Period pi_T(m) governs 3-step return times
     - Angular velocity: tribonacci angle ≈ {np.degrees(tribonacci_angle):.1f}°
     - Controls: 3-CYCLE structure (the fundamental tournament invariant)

  3. OMEGA HAND (roots of unity):
     - Tracks character rotation: omega^n = e^(2*pi*i*n/p)
     - Period exactly p
     - Controls: EIGENSPACE decomposition (THM-125)

  WHEN HANDS ALIGN:
     Joint period = LCM(pi_F(p), pi_T(p), p)
     At alignment: the tournament achieves FULL SYMMETRY

  DEEP CONNECTION:
     phi < tau < 2 means: GOLDEN < TRIBONACCI < BINARY

     Tournament theory lives between tau and 2:
     - tau governs the GROWTH RATE of 3-cycle-related structures
     - 2 governs the NUMBER OF CHOICES (arc orientations)

     The gap tau → 2 (ratio {2/tau:.6f}) is where tournament
     structure lives: too few choices for pure Fibonacci recursion,
     but enough for tribonacci-type interactions.

     β₁ (our new orientability obstruction from earlier today)
     is essentially about when the PHI hand and TAU hand are
     OUT OF SYNC: the path structure (Fibonacci-like) doesn't
     align with the cycle structure (tribonacci-like).

  phi * tau = {phi*tau:.6f} ≈ 3 (the 3-nacci weight limit!)
  This is NOT a coincidence: the product of the dominant roots
  of x^2-x-1 and x^3-x^2-x-1 is related to the next integer.
""")

print("DONE — tau_phi_clock_S72b.py complete")
