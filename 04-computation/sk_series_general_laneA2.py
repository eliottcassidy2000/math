"""
Lane A2 :  the general closed form for

    S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ((kn+1) 64^n),   k a positive integer.

Everything below is reproducible with mpmath alone.  Notation

    a_n  = C(2n,n)C(4n,2n)/64^n = (1/4)_n (3/4)_n / (n!)^2
    f(z) = 2F1(1/4,3/4;1;z)   (so sum a_n z^n = f(z))
    G(s) = int_0^1 z^{s-1} f(z) dz = sum_n a_n/(n+s),      S(k) = (1/k) G(1/k)
    L    = log(1+sqrt2),   varpi = int_0^1 dv/sqrt(1-v^4) = Gamma(1/4)^2/(4 sqrt(2 pi))

Contents
  1. numerics(...)              two independent >=160-digit routes for S(k)
  2. thomae_form(k)             THEOREM A2-1  S(k)=(8 sqrt2/(3 pi k)) 3F2(1-1/k,1,1;5/4,7/4;1)
  3. sector_form(k)             THEOREM A2-2  sector / surface-period representation
  4. ellipticK_form(k)          THEOREM A2-3  master elliptic-K integral
  5. recursion_check()          THEOREM A2-4  (s-1)^2 G(s-1) = (s-1/4)(s-3/4) G(s) - 1/(pi sqrt2)
  6. k4_analysis()              S(4) reduced to a CM elliptic-curve Eichler integral (j = 66^3)
  7. pslq_bounded(...)          the bounded-coefficient PSLQ protocol used for the refutations
"""
from mpmath import (mp, mpf, mpc, sqrt, pi, log, atan, asin, asinh, quad, hyp2f1,
                    hyp3f2, cos, sin, tan, catalan, gamma, pslq, nstr, findpoly, polylog)

# ----------------------------------------------------------------- 1. numerics
def S_routeA(k):
    """S(k) = int_0^1 2F1(1/4,3/4;1;t^k) dt   (integrable log endpoint singularity)."""
    kk = mpf(k)
    return quad(lambda t: hyp2f1(mpf(1)/4, mpf(3)/4, 1, t**kk), [0, 1])

def S_routeB(k):
    """S(k) = (2/pi) int_0^{pi/2} 2F1(1/2, 2/k; 1+2/k; cos 2th) dth  (preamble (v)).
       Integrand is ANALYTIC on [0,pi/2] because c-a-b = 1/2."""
    kk = mpf(k)
    g = lambda th: hyp2f1(mpf(1)/2, 2/kk, 1+2/kk, cos(2*th))
    return (2/pi) * quad(g, [0, pi/4, pi/2])

# --------------------------------------------------------- 2. THEOREM A2-1
def thomae_form(k):
    """S(k) = (8 sqrt2/(3 pi k)) * 3F2(1-1/k, 1, 1; 5/4, 7/4; 1).
       Obtained from S(k)=3F2(1/4,3/4,1/k;1,1+1/k;1) by Thomae singling out c=1/k:
         3F2(a,b,c;e,f;1) = G(e)G(f)G(s)/(G(c)G(s+a)G(s+b)) 3F2(e-c,f-c,s;s+a,s+b;1)
       with s=1, G(1+1/k)/G(1/k)=1/k, G(5/4)G(7/4)=3 pi sqrt2/16.
       Equivalent term form uses (5/4)_n (7/4)_n = 4^{-n} (5/2)_{2n}."""
    kk = mpf(k)
    return (8*sqrt(2)/(3*pi*kk)) * hyp3f2(1-1/kk, 1, 1, mpf(5)/4, mpf(7)/4, 1)

def thomae_euler_double_integral(k):
    """Euler form of THEOREM A2-1 (v-integral outer so the only singular point is an endpoint):
         S(k) = (8 sqrt2/(pi k)) int_0^1 int_0^{pi/2}
                  sin th cos^2 th (cos^2 2th + v^2 sin^2 2th)^{1/k-1} dth dv."""
    kk = mpf(k); A = 1 - 1/kk
    def th_int(v):
        g = lambda th: sin(th)*cos(th)**2*(cos(2*th)**2 + v**2*sin(2*th)**2)**(-A)
        return quad(g, [0, pi/4, pi/2])
    return (8*sqrt(2)/(pi*kk)) * quad(th_int, [0, 1])

# --------------------------------------------------------- 3. THEOREM A2-2
def sector_form(k):
    """S(k) = (2k/pi) int int_{Sigma_k} (1 + Im(zeta^k))^{-1/2} dA,
       Sigma_k = { |zeta|<1, |arg zeta| < pi/(2k) }.   In polar zeta = r e^{i psi}."""
    kk = mpf(k)
    inner = lambda r: quad(lambda ps: r*(1 + r**kk*sin(kk*ps))**mpf('-0.5'),
                           [-pi/(2*kk), 0, pi/(2*kk)])
    return (2*kk/pi) * quad(inner, [0, 1])

# --------------------------------------------------------- 4. THEOREM A2-3
def ellipticK_form(k):
    """S(k) = (16/(sqrt2 pi k)) int_0^1 K(kappa) kappa^{4/k-1} (2-kappa^2)^{-1/2-2/k} dkappa.
       (mpmath's ellipk takes the PARAMETER m = kappa^2.)"""
    from mpmath import ellipk
    kk = mpf(k)
    g = lambda c: ellipk(c**2)*c**(4/kk-1)*(2-c**2)**(-mpf(1)/2-2/kk)
    return (16/(sqrt(2)*pi*kk)) * quad(g, [0, 1])

# --------------------------------------------------------- 5. THEOREM A2-4
def G(s):
    return quad(lambda z: z**(mpf(s)-1)*hyp2f1(mpf(1)/4, mpf(3)/4, 1, z), [0, 1])

def recursion_residual(s):
    """(s-1)^2 G(s-1) - (s-1/4)(s-3/4) G(s) + 1/(pi sqrt2)  == 0 for Re s > 1.
       Derived by integrating z^{s-1} * [ z(1-z)f'' + (1-2z)f' - (3/16) f ] = 0;
       the boundary term is z^{s-1}(1-z)[z f' - (s-1) f], which tends to
       1/(pi sqrt2) at z=1 (because f ~ -(1/(pi sqrt2)) log(1-z)) and to 0 at z=0."""
    s = mpf(s)
    return (s-1)**2*G(s-1) - (s-mpf(1)/4)*(s-mpf(3)/4)*G(s) + 1/(pi*sqrt(2))

def shifted_family(mmax=8):
    """Corollary: sum_n a_n/(n+m) = sqrt2 * r_m / pi with r_m in Q, r_1 = 8/3 and
         r_{m+1} = ((m)^2 r_m + 1/2) / ((m+3/4)(m+1/4)).
       (equivalently G(m) = Gamma(m)^2/(Gamma(m+1/4)Gamma(m+3/4)) only at m=1.)"""
    from fractions import Fraction as F
    r = F(8, 3); out = [r]
    for m in range(1, mmax):
        s = F(m+1)
        r = ((s-1)**2*r + F(1, 2)) / ((s-F(1, 4))*(s-F(3, 4)))
        out.append(r)
    return out

# --------------------------------------------------------- 6. the k = 4 object
def k4_pieces():
    """S(4) = (2/pi)(I1+I2) with the two 'halves' of preamble (vii) after th = arcsin v:
         I1 = int_0^{pi/2} th / sqrt(1+sin^2 th) dth        (the arcsin half)
         I2 = int_0^{pi/2} asinh(sin th)/sqrt(1+sin^2 th) dth (the arcsinh half)
       Both integrands are analytic on [0,pi/2] -> 400 digits are cheap.
       S(4) also equals (2 sqrt2/pi) int_0^1 K(kappa)/(2-kappa^2) dkappa
       and (2/pi) int_0^{pi/2} log(sqrt2 cos th + sqrt(cos 2th))/sqrt(cos 2th) dth."""
    I1 = quad(lambda th: th/sqrt(1+sin(th)**2), [0, pi/2])
    I2 = quad(lambda th: asinh(sin(th))/sqrt(1+sin(th)**2), [0, pi/2])
    return I1, I2

def k4_curve_j():
    """I1+I2 = int_{iL}^{pi/2} z dz / sqrt(2-cos^2 z),  L = log(1+sqrt2).
       W = e^{2iz} turns this into -(1/2) int log(W) dW/Y on  Y^2 = W(6W-W^2-1),
       an elliptic curve with 2-torsion {0, 3+2 sqrt2, 3-2 sqrt2}.  Its j-invariant:"""
    lam = (3-2*sqrt(2))**2                       # = 17 - 12 sqrt2
    return 256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2)   # = 287496 = 66^3  (CM, disc -16)

# --------------------------------------------------------- 7. PSLQ protocol
def pslq_bounded(target, basis, digits, maxcoeff):
    """Bounded-coefficient PSLQ.  With n numbers at D digits, PSLQ finds SPURIOUS
       relations of size ~10^{D/n}; choosing maxcoeff well below that makes a
       `None` meaningful.  ALWAYS run pslq(basis) first: a hit with a zero leading
       coefficient is an internal identity of the basis, not a hit on the target."""
    mp.dps = digits
    internal = pslq([+b for b in basis], maxcoeff=maxcoeff, maxsteps=10**6,
                    tol=mpf(10)**(-(digits-10)))
    rel = pslq([+target]+[+b for b in basis], maxcoeff=maxcoeff, maxsteps=10**6,
               tol=mpf(10)**(-(digits-10)))
    return internal, rel


if __name__ == '__main__':
    import sys, time
    mp.dps = 160
    print("== S(k), two routes, dps=160 ==")
    for k in range(1, 13):
        a, b = S_routeA(k), S_routeB(k)
        print("k=%2d %s   |A-B|=%s" % (k, mp.nstr(b, 60), mp.nstr(abs(a-b), 5)))
        sys.stdout.flush()
    mp.dps = 50
    print("\n== THEOREM A2-1 (Thomae) ==")
    for k in [1, 2, 3, 4, 5, 6, 7, 10]:
        print("k=%2d residual %s" % (k, mp.nstr(S_routeB(k)-thomae_form(k), 5)))
    print("\n== THEOREM A2-3 (elliptic K) ==")
    for k in [1, 2, 3, 4, 5]:
        print("k=%2d residual %s" % (k, mp.nstr(S_routeB(k)-ellipticK_form(k), 5)))
    print("\n== THEOREM A2-4 (contiguous recursion) ==")
    for s in ['1.5', '1.7', '2', '2.2', '7/3']:
        ss = mpf(s) if '/' not in s else mpf(7)/3
        print("  s=%-6s residual %s" % (s, mp.nstr(recursion_residual(ss), 5)))
    print("  r_m =", shifted_family())
    print("\n== k=4 ==")
    I1, I2 = k4_pieces()
    print("  I1 =", mp.nstr(I1, 40)); print("  I2 =", mp.nstr(I2, 40))
    print("  (2/pi)(I1+I2) - S(4) =", mp.nstr((2/pi)*(I1+I2)-S_routeB(4), 5))
    print("  j-invariant of the k=4 curve =", mp.nstr(k4_curve_j(), 12), " (66^3 =", 66**3, ")")
