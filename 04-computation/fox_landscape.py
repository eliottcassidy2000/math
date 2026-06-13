#!/usr/bin/env python3
"""fox_landscape.py — The fox sees the whole landscape, not one equation.

Session: kind-pasteur-2026-03-20-S12

THE HEDGEHOG says: "2^10 = 10^3 + 4! is the equation."
THE FOX says: "That's one point on a surface. What's the surface?"

The surface is S(y, k) = 2^y - y^k.
The tower equation asks: when is S(y, k) = b! for some b?
The near-miss asks: when is S(y, k) CLOSE to a factorial?

But the fox goes further: maybe the factorials are just one family
of landmarks on the surface. What OTHER landmarks are there?
"""

from math import factorial, comb, log, log2, sqrt
from collections import defaultdict

def surface_scan():
    """Scan the surface S(y,k) = 2^y - y^k for ALL factorial hits and near-misses."""
    print("=" * 70)
    print("THE FOX'S LANDSCAPE: S(y,k) = 2^y - y^k")
    print("=" * 70)

    # Find ALL (y, k, b) with 2^y - y^k = b! (exact)
    # and ALL near-misses within 1%
    print(f"\n  EXACT factorial gaps: 2^y = y^k + b!")
    exact = []
    near = []

    for k in range(1, 12):
        for y in range(1, 80):
            try:
                gap = 2**y - y**k
            except OverflowError:
                break
            if gap <= 0:
                continue
            # Check if gap is a factorial
            for b in range(1, 25):
                fb = factorial(b)
                if fb == gap:
                    exact.append((y, k, b, gap))
                    print(f"    y={y:3d}, k={k:2d}: 2^{y} = {y}^{k} + {b}! "
                          f"({2**y} = {y**k} + {fb})")
                elif fb > gap * 2:
                    break
                elif abs(fb - gap) < gap * 0.005 and gap > 100:
                    pct = abs(fb - gap) / gap * 100
                    near.append((y, k, b, gap, fb, pct))

    print(f"\n  Near-misses (within 0.5% and gap > 100):")
    for y, k, b, gap, fb, pct in sorted(near, key=lambda t: t[5])[:20]:
        sign = '+' if fb > gap else '-'
        diff = abs(fb - gap)
        print(f"    y={y:3d}, k={k:2d}: 2^{y} - {y}^{k} = {gap}, "
              f"nearest {b}! = {fb}, off by {sign}{diff} ({pct:.3f}%)")


def the_exact_spine():
    """The exact solutions of 2^y = y^k form the SPINE of the landscape."""
    print(f"\n{'='*70}")
    print(f"THE SPINE: Exact solutions of 2^y = y^k")
    print(f"{'='*70}")

    # 2^y = y^k iff y*ln2 = k*ln(y) iff k = y*ln2/ln(y)
    # Integer solutions: y = 2^m where m | 2^m, i.e., m is a power of 2
    # y=2: k=2 (2^2 = 4 = 2^2)
    # y=4: k=2 (2^4 = 16 = 4^2)
    # y=16: k=4 (2^16 = 65536 = 16^4)
    # y=256: k=32 (2^256 = 256^32)

    print(f"\n  Exact solutions (y, k) where 2^y = y^k:")
    for m in range(1, 20):
        y = 2**m
        if m > 0 and (2**m) % m == 0:
            k = (2**m) // m
            # Verify: y^k should equal 2^y
            # y^k = (2^m)^k = 2^(mk) = 2^(2^m) = 2^y. Check: mk = y?
            # mk = m * (2^m/m) = 2^m = y. Yes!
            print(f"    m={m:3d}: y = 2^{m} = {y}, k = 2^{m}/{m} = {k}")
            if y <= 1000:
                print(f"           Verify: {y}^{k} = 2^{y} = {2**y}")

    print(f"""
  The spine is the TETRATION staircase: y = 2, 4, 16, 256, 65536, ...
  These are 2^^j = 2, 4, 16, 65536, ... (tower of 2s, j times)

  Between spine points, 2^y - y^k oscillates. The tower equation
  2^y = y^3 + 4! lives between spine points y=4 (k=2) and y=16 (k=4).

  At k=3 (between spine exponents 2 and 4): the "midpoint" behavior
  produces a gap that happens to be 4! at y=10.
""")


def between_spine_points():
    """What happens between spine points y=4 and y=16?"""
    print(f"\n{'='*70}")
    print(f"BETWEEN THE SPINE: y=4 (k=2) to y=16 (k=4)")
    print(f"{'='*70}")

    # For each integer k from 2 to 4, find where 2^y crosses y^k
    # and what the gap looks like
    for k in [2, 3, 4]:
        print(f"\n  k={k}: 2^y vs y^{k}")
        crossing_y = None
        for y in range(1, 30):
            gap = 2**y - y**k
            if gap >= 0 and crossing_y is None:
                crossing_y = y
            if abs(y - 10) <= 5 or (crossing_y and abs(y - crossing_y) <= 3):
                print(f"    y={y:3d}: 2^{y} - {y}^{k} = {gap:>12d}", end="")
                # Check factorials
                for b in range(1, 20):
                    if factorial(b) == gap:
                        print(f"  = {b}!", end="")
                        break
                    elif factorial(b) == abs(gap):
                        print(f"  = -{b}!", end="")
                        break
                print()


def the_real_equation():
    """What equation GENERATES both the exact hit and near-miss?"""
    print(f"\n{'='*70}")
    print(f"THE FOX'S QUESTION: What generates BOTH y=10 and y=11?")
    print(f"{'='*70}")

    # At y=10: 2^10 - 10^3 = 24 = 4!
    # At y=11: 2^11 - 11^3 = 717 = 6! - 3

    # What if the REAL equation is:
    # 2^y - y^3 = Gamma(f(y)) for some smooth function f?
    # f(10) = 5 (since Gamma(5) = 4! = 24)
    # f(11) should give Gamma(f(11)) = 717
    # Gamma(6.something) ≈ 717? Gamma(6) = 120, Gamma(7) = 720.
    # So f(11) ≈ 6.996 (since Gamma(6.996) ≈ 717)

    # Using Stirling: Gamma(z+1) ≈ z^z * e^{-z} * sqrt(2*pi*z)
    # At z=6: Gamma(7) = 720. 717/720 = 0.9958. So z ≈ 6 - epsilon.

    import math
    # Inverse gamma: find z such that Gamma(z+1) = 717
    # Newton's method starting from z=6
    z = 6.0
    for _ in range(20):
        val = math.gamma(z + 1)
        # Derivative: psi(z+1) * Gamma(z+1)
        dval = val * (math.log(z + 0.5))  # rough approximation
        z = z - (val - 717) / dval
    print(f"\n  Gamma^{{-1}}(717) ≈ {z:.6f}")
    print(f"  So Gamma({z:.6f} + 1) ≈ {math.gamma(z+1):.4f}")

    # The function f(y) that generates both points:
    # f(10) = 5 (Gamma(5) = 24)
    # f(11) ≈ 6.996 (Gamma(6.996) ≈ 717)
    f10 = 5
    f11 = z  # approximately 6.996

    print(f"\n  f(10) = {f10}")
    print(f"  f(11) = {f11:.6f}")
    print(f"  f(11) - f(10) = {f11 - f10:.6f}")
    print(f"  Almost exactly 2! The slope is ≈ 2.")

    # If f(y) = 2y - 15 (linear with slope 2):
    # f(10) = 20 - 15 = 5 ✓
    # f(11) = 22 - 15 = 7 → Gamma(7+1) = 5040 ≠ 717
    # That doesn't work because Gamma is not linear.

    # What if the REAL relation is:
    # 2^y = y^3 + Gamma(2y - 15)  ?
    # y=10: Gamma(5) = 24 ✓
    # y=11: Gamma(7) = 720. But actual gap is 717 = 720 - 3. Off by 3.

    print(f"\n  Testing: 2^y = y^3 + Gamma(2y - 15)")
    for y in range(8, 15):
        arg = 2*y - 15
        if arg > 0:
            gamma_val = math.gamma(arg)
            predicted = y**3 + gamma_val
            actual = 2**y
            error = actual - predicted
            print(f"    y={y}: 2^{y}={actual}, y^3+Gamma({arg})={y**3}+{gamma_val:.0f}"
                  f"={predicted:.0f}, error={error:.0f}")

    # The "3" error at y=11 might not be the min_cycle coincidence.
    # Let's look at it differently.

    print(f"\n  The error 3 at y=11:")
    print(f"  2^11 - 11^3 - 6! = 2048 - 1331 - 720 = {2048-1331-720}")
    print(f"  = -3")
    print(f"")
    print(f"  What IS 3 here?")
    print(f"  3 = 11 - 8 = y - 2^3")
    print(f"  3 = min(y-1, y-2, ...) such that 2^y = y^3 + (2y-15)! - something")
    print(f"  Or simply: 717 = 3 * 239, where 239 is prime.")
    print(f"")
    print(f"  The fox says: maybe 3 is just 3. Not everything is meaningful.")

    # But let's look at the FULL landscape of errors
    print(f"\n  Full error landscape: 2^y - y^3 - Gamma(2y-15)")
    for y in range(8, 20):
        arg = 2*y - 15
        if arg > 0:
            try:
                gamma_val = math.gamma(arg)
                error = 2**y - y**3 - gamma_val
                print(f"    y={y}: arg={arg}, error = {error:.1f}")
            except OverflowError:
                print(f"    y={y}: overflow")


def what_the_fox_sees():
    """The fox's synthesis: the landscape has many features, not just one equation."""
    print(f"\n{'='*70}")
    print(f"WHAT THE FOX SEES")
    print(f"{'='*70}")

    print(f"""
  The hedgehog says: "2^10 = 10^3 + 4! is THE identity."

  The fox says: "That's one point on a vast surface.
  The surface S(y,k) = 2^y - y^k has:

  1. A SPINE of exact zeros at y = 2, 4, 16, 256, ...
     (the tetration staircase, where 2^y = y^k exactly)

  2. FACTORIAL CONTOURS where S = b! (our tower equation)
     Only one non-trivial hit at (y=10, k=3, b=4).

  3. NEAR-MISS CONTOURS where S ≈ b!
     At (y=11, k=3): S = 717 ≈ 720 = 6! (off by 3)
     At (y=1, k=any): S = 1 = 1! (trivial)

  4. OTHER LANDMARKS we haven't looked for:
     - When is S a perfect power?
     - When is S a binomial coefficient?
     - When is S a Fibonacci number?
     - When is S a primorial?

  The hedgehog-equation 2^y = y^(x+1) + (x+2)! is the FACTORIAL CONTOUR
  of the surface. It's beautiful but it's one curve on a rich landscape."

  Let me look for other landmarks:
""")

    # Check: when is 2^y - y^3 a perfect power?
    print(f"  Perfect powers on the k=3 slice:")
    for y in range(1, 30):
        gap = 2**y - y**3
        if gap <= 0:
            continue
        # Check if gap = a^b for small b
        for b in range(2, 20):
            a_float = gap ** (1.0/b)
            a = round(a_float)
            if a > 1 and a**b == gap:
                print(f"    y={y}: 2^{y} - {y}^3 = {gap} = {a}^{b}")
                break

    # Check: when is gap a binomial coefficient?
    print(f"\n  Binomial coefficients on the k=3 slice:")
    binom_set = set()
    for n in range(1, 50):
        for r in range(n+1):
            binom_set.add(comb(n, r))

    for y in range(1, 25):
        gap = 2**y - y**3
        if gap > 0 and gap in binom_set:
            # Find which C(n,r)
            for n in range(1, 50):
                for r in range(n+1):
                    if comb(n, r) == gap and r > 0 and r < n:
                        print(f"    y={y}: 2^{y} - {y}^3 = {gap} = C({n},{r})")

    # Check: Fibonacci numbers
    print(f"\n  Fibonacci numbers on the k=3 slice:")
    fibs = set()
    a, b = 1, 1
    while a < 10**8:
        fibs.add(a)
        a, b = b, a+b

    for y in range(1, 25):
        gap = 2**y - y**3
        if gap > 0 and gap in fibs:
            print(f"    y={y}: 2^{y} - {y}^3 = {gap} (Fibonacci)")

    # Check: triangular numbers
    print(f"\n  Triangular numbers on the k=3 slice:")
    for y in range(1, 25):
        gap = 2**y - y**3
        if gap <= 0:
            continue
        # T_n = n(n+1)/2, so gap = n(n+1)/2, n = (-1 + sqrt(1+8*gap))/2
        disc = 1 + 8*gap
        sq = int(sqrt(disc))
        if sq*sq == disc and (sq - 1) % 2 == 0:
            n = (sq - 1) // 2
            print(f"    y={y}: 2^{y} - {y}^3 = {gap} = T({n}) = {n}*{n+1}/2")


if __name__ == "__main__":
    surface_scan()
    the_exact_spine()
    between_spine_points()
    the_real_equation()
    what_the_fox_sees()
