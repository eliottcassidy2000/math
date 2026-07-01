#!/usr/bin/env python3
"""
gauss_sqrt21_aut21_bridge_klein.py  --  klein-2026-07-01-S84

The ι-odd certificate next step + the |Aut|=21 <-> sqrt21 convergence (owner task 5, building on opus-S26
HYP-3818 and my HYP-3814/3817).

opus-S26: OPEN-Q-108's odd index = the iota-ODD Gauss certificate i*sqrt7 (heptagon, p=7=3mod4); the two
order-2 structures = Gal(Q(sqrt-3,sqrt-7))=Z2xZ2; g(3)*g(7) = i*sqrt3 * i*sqrt7 = -sqrt21 = the iota-EVEN
cross-term (cup H^1_odd (x) H^1_odd -> H^2_even). sqrt21 = the OPEN-Q-108 residual.

MY CONVERGENCE (S82/S83): the flip-rank OBSTRUCTION at n=7 = the Paley heptagon, |Aut| = 21 = 3*7 (argmax,
HYP-3817); its CAYLEY spectrum (HYP-3814) sits at cos(theta) = -(p-1)/(p+1) = -3/4, encoding g(7) = i*sqrt7.
So 21 = 3*7 appears as BOTH |Aut(Paley_7)| (the covering obstruction) AND norm(sqrt21) (the LRC certificate
residual). The SHARED prime 3 = Eisenstein/covering (Phi6); the 7 = heptagon/odd index (14=2*7).

This script: (1) verifies the Gauss-sum iota-parity g(3)=i*sqrt3, g(7)=i*sqrt7 (odd, imaginary),
g(21)=+-sqrt21 (even, real); (2) verifies the cross-term g(3)g(7)=-sqrt21; (3) tabulates the 21=3*7 bridge.
"""
import cmath, math

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def gauss_quadratic(p):
    """g(p) = sum_{a=1}^{p-1} (a|p) * exp(2 pi i a / p)."""
    z = cmath.exp(2j*math.pi/p)
    return sum(legendre(a, p) * z**a for a in range(1, p))

def jacobi_symbol_char(a, m, factors):
    """Real quadratic-ish character mod m = product of Legendre symbols over prime factors."""
    r = 1
    for p in factors:
        r *= legendre(a, p)
    return r

def gauss_composite(m, factors):
    z = cmath.exp(2j*math.pi/m)
    return sum(jacobi_symbol_char(a, m, factors) * z**a for a in range(1, m) if math.gcd(a, m) == 1)

if __name__ == "__main__":
    print("="*74)
    print("(1) Gauss-sum iota-parity  (p = 3 mod 4 => g = i*sqrt(p), IMAGINARY = iota-ODD)")
    print("="*74)
    for p in (3, 7, 11, 19):
        g = gauss_quadratic(p)
        print(f"  g({p}) = {g:.5f}   |g|^2={abs(g)**2:.3f} (=p={p});  i*sqrt({p}) = {1j*math.sqrt(p):.5f}  "
              f"=> {'i*sqrt(p) (iota-ODD, p=3mod4)' if p%4==3 else 'sqrt(p) (iota-EVEN)'}")

    print("\n" + "="*74)
    print("(2) THE CROSS-TERM  g(3)*g(7) = i*sqrt3 * i*sqrt7 = -sqrt21  (iota-ODD (x) iota-ODD -> iota-EVEN)")
    print("="*74)
    g3, g7 = gauss_quadratic(3), gauss_quadratic(7)
    prod = g3*g7
    print(f"  g(3)*g(7) = {prod:.5f}   (imag ~ 0 => REAL);  -sqrt21 = {-math.sqrt(21):.5f}")
    print(f"  sqrt21 = {math.sqrt(21):.5f} = sqrt(3*7);  in Q(sqrt-3,sqrt-7): sqrt21 = (i*sqrt3)(i*sqrt7)/(-1)")
    print(f"  g(21) composite (char mod 21 = (.|3)(.|7)): {gauss_composite(21,[3,7]):.5f}  (real, ~ +-sqrt21*unit)")
    print(f"  Gal(Q(sqrt-3,sqrt-7)) = Z2 x Z2; fixed fields Q(sqrt-3),Q(sqrt-7),Q(sqrt21); sqrt21 needs BOTH primes")

    print("\n" + "="*74)
    print("(3) THE 21 = 3*7 BRIDGE  (flip-rank obstruction  <->  LRC sqrt21 certificate)")
    print("="*74)
    print("  |Aut(Paley_7)| = 21 = 3*7  = argmax|Aut| = the n=7 flip-rank OBSTRUCTION (HYP-3817, opus-S15)")
    print("  norm(sqrt21)   = 21 = 3*7  = the OPEN-Q-108 residual cross-term (opus-S26 HYP-3818)")
    print("  Paley_7 CAYLEY spectrum: cos(theta) = -(7-1)/(7+1) = -3/4  <=>  eig(S)=+-i*sqrt7 = g(7) (HYP-3814)")
    print("  SHARED prime 3 = Eisenstein/covering (Phi6(14)=183=3*61);  7 = heptagon/odd index (14=2*7)")
    print("  => the covering obstruction (|Aut|=21) and the certificate residual (sqrt21) are the SAME 3*7 arithmetic:")
    print("     3 (iota-odd i*sqrt3, covering) x 7 (iota-odd i*sqrt7, odd index) -> 21 (iota-even sqrt21, residual).")
    print("  The flip-rank obstruction is a COVERING witness (fixed-point/atom, HYP-3817);")
    print("  the sqrt21 residual is its MOMENT/certificate shadow -- 'reach for a covering or a moment, not a transform.'")
