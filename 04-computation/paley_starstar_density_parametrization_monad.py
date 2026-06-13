#!/usr/bin/env python3
"""
THM-438 ADDENDUM-15 (monad-explorer 17th session)
=================================================
CLOSED-FORM PARAMETRIC DENSITY of the free factorial law mu_free
( = free compound Poisson of nu=e^{-x}dx; free cumulants kappa_n=n!;
  moments A088368 = 1,1,3,13,69,421,... ).

BACKGROUND.  ADD-14 (16th) gave the K-transform (inverse Cauchy transform) in
closed form:
        K(w) = -(1/w^2) e^{-1/w} E1(-1/w) = -(1/w^2) g(-1/w),
        g(u) := e^u E1(u) = int_0^inf e^{-t}/(u+t) dt   (Gompertz function),
with  G(K(w)) = w  and  rho(x) = -(1/pi) Im G(x+i0).
ADD-13 (15th) computed rho only by NUMERICAL root-finding G = K^{-1}(x+i eta)
with an artificial regulator eta.

THIS SESSION.  Because K is explicit, the density is an EXACT PARAMETRIC CURVE,
needing no root-finding and no eta.  Substitute u = -1/w (so w = G = -1/u):
        K(w) = -u^2 g(u).
On the support, x = K(G(x+i0)) is REAL, so the support+density are read directly
off the single real condition  Im( u^2 g(u) ) = 0 :

        PARAMETRIC DENSITY (closed form)
        --------------------------------
        x(u)   = - u^2 g(u)                       (real on the curve)
        rho(u) = - Im(u) / ( pi |u|^2 )
        on the curve  C = { u in lower half-plane : Im( u^2 g(u) ) = 0 }.

        Edge x->0+  <=>  u->0      (E1 log singularity; arg u -> -pi/2)
        Tail x->inf <=>  u->-inf   (E1 cut on (-inf,0]; gives the tail constant)

This script:
 (1) Traces the curve C (for each modulus eps=|u|, solve arg(u) in (-pi,0)),
     producing (x, rho) with NO root-finding of G and NO regulator eta.
 (2) Cross-checks rho(x) against ADD-14's root-found density at several x.
 (3) Recovers the A088368 moments by integrating x^k rho along C.
 (4) EDGE: verifies x ~ eps^2 (ln(1/eps) - gamma), rho ~ 1/(pi eps),
     hence pi rho sqrt(x) -> sqrt(ln(1/eps)) = sqrt( (1/2) ln(1/x) )  (ADD-14),
     and reads the next-order correction directly from the parametrization.
 (5) TAIL: u -> -inf along the real axis; the E1 cut Im E1(x +- i0) = -+ pi
     for x<0 gives  rho(x) ~ e * e^{-x}  (the constant e of MISTAKE-063, now
     DERIVED from the closed form rather than fitted).
"""
import mpmath as mp

mp.mp.dps = 40
A088368 = [1,1,3,13,69,421,2867,21477,175769,1567273]

def g(u):
    """Gompertz function g(u) = e^u E1(u), principal branch of E1."""
    return mp.e**(u) * mp.e1(u)

def K_of_w(w):
    """K(w) = -(1/w^2) e^{-1/w} E1(-1/w)  (ADD-14 closed form)."""
    s = 1/w
    return -(s*s) * mp.e**(-s) * mp.e1(-s)

def xK(u):
    """x(u) = -u^2 g(u) = K(-1/u). Real on the curve C."""
    return -(u*u) * g(u)

# ------------------------------------------------------------------
def trace_curve(eps_list):
    """For each eps=|u|, solve phi in (-pi,0) with Im( u^2 g(u) ) = 0, u=eps e^{i phi}.
    Returns lists of (eps, phi, x, rho)."""
    out = []
    phi_guess = -mp.pi/2          # edge end: arg u -> -pi/2
    for eps in eps_list:
        def F(phi):
            u = eps*mp.e**(1j*phi)
            return mp.im(u*u*g(u))
        try:
            phi = mp.findroot(F, phi_guess)
        except Exception:
            # widen search
            phi = mp.findroot(F, phi_guess - mp.mpf('0.2'))
        u = eps*mp.e**(1j*phi)
        x = mp.re(xK(u))                 # = -Re(u^2 g(u)); Im part ~ 0 by construction
        rho = -mp.im(u)/(mp.pi*eps*eps)
        out.append((eps, phi, x, rho))
        phi_guess = phi
    return out

# ------------------------------------------------------------------
def section1_verify_vs_rootfound():
    print("="*74)
    print("(2) Parametric rho(x) vs ADD-14 root-found density  rho=-(1/pi)Im G(x+i eta)")
    print("="*74)
    eta = mp.mpf('1e-12')
    def rho_rootfound(x, guess):
        z = mp.mpf(x) + 1j*eta
        G = mp.findroot(lambda G: K_of_w(G) - z, guess)
        if mp.im(G) > 0: G = mp.conj(G)
        return -mp.im(G)/mp.pi, G
    # pick several u, get (x,rho) from the curve, compare with root-found at that x
    print("   x            rho_param      rho_rootfound   |diff|")
    test_eps = [mp.mpf(e) for e in ['0.2','0.4','0.7','1.0','1.5','2.0','3.0']]
    rows = trace_curve(test_eps)
    for (eps,phi,x,rho) in rows:
        if x <= 0: continue
        w_guess = -1/(eps*mp.e**(1j*phi))   # = G
        rr, G = rho_rootfound(x, w_guess)
        print(f"   {mp.nstr(x,6):>10}   {mp.nstr(rho,8):>12}   {mp.nstr(rr,8):>12}   {mp.nstr(abs(rho-rr),3)}")

# ------------------------------------------------------------------
def section2_moments():
    print("\n"+"="*74)
    print("(3) Moment recovery  int_0^inf x^k rho(x) dx = A088368(k)  (param. curve)")
    print("="*74)
    # Parametrize by eps over a wide range; x is monotone increasing in eps.
    # Use fine grid; integrate in x via trapezoid on the (x,rho) samples.
    N = 1400
    eps_lo, eps_hi = mp.mpf('0.002'), mp.mpf('14.0')
    # geometric-ish grid denser at small eps (edge) and we also need large x tail
    eps_list = [eps_lo*(eps_hi/eps_lo)**(mp.mpf(i)/N) for i in range(N+1)]
    rows = trace_curve(eps_list)
    pts = sorted([(x,rho) for (_,_,x,rho) in rows if x>0 and rho>=0])
    xs = [p[0] for p in pts]; rho = [p[1] for p in pts]
    print(f"   curve samples with x>0: {len(xs)}   x in [{mp.nstr(xs[0],3)}, {mp.nstr(xs[-1],3)}]")
    def moment(k):
        s = mp.mpf(0)
        for i in range(len(xs)-1):
            a = xs[i]**k * rho[i]
            b = xs[i+1]**k * rho[i+1]
            s += (a+b)/2*(xs[i+1]-xs[i])
        return s
    print("\n   k    integral        A088368     rel.err")
    for k in range(0,6):
        m = moment(k); tgt = A088368[k]
        print(f"   {k}   {mp.nstr(m,9):>13}   {tgt:>8}   {mp.nstr(abs(m-tgt)/tgt,3)}")

# ------------------------------------------------------------------
def section3_edge():
    print("\n"+"="*74)
    print("(4) EDGE  u->0 :  x ~ eps^2(ln(1/eps)-gamma),  rho ~ 1/(pi eps),")
    print("        pi rho sqrt(x) -> sqrt(ln(1/eps)) ~ sqrt((1/2)ln(1/x))")
    print("="*74)
    gam = mp.euler
    eps_list = [mp.mpf(10)**(-mp.mpf(j)) for j in range(1,9)]   # eps=1e-1..1e-8
    rows = trace_curve(eps_list)
    print("   eps          x             rho          pi*rho*sqrt(x)   sqrt(ln(1/eps))  ratio")
    for (eps,phi,x,rho) in rows:
        lhs = mp.pi*rho*mp.sqrt(x)
        ref = mp.sqrt(-mp.log(eps))     # sqrt(ln(1/eps))
        print(f"   {mp.nstr(eps,2):>9}  {mp.nstr(x,6):>12}  {mp.nstr(rho,6):>10}   "
              f"{mp.nstr(lhs,8):>9}   {mp.nstr(ref,8):>9}   {mp.nstr(lhs/ref,6)}")
    # also test the leading edge model  x = eps^2(ln(1/eps)-gamma)
    print("\n   leading edge model  x_model = eps^2(ln(1/eps)-gamma):")
    for (eps,phi,x,rho) in rows:
        xm = eps*eps*(-mp.log(eps)-gam)
        print(f"   eps={mp.nstr(eps,2):>9}  x={mp.nstr(x,6):>12}  x_model={mp.nstr(xm,6):>12}  ratio={mp.nstr(x/xm,6)}")
    # phi -> -pi/2; report delta = phi+pi/2 and compare to pi/(4 L), L=gamma+ln eps
    print("\n   arg-deviation  delta=phi+pi/2  vs  pi/(4(gamma+ln eps)):")
    for (eps,phi,x,rho) in rows:
        delta = phi + mp.pi/2
        L = gam + mp.log(eps)
        pred = mp.pi/(4*L)
        print(f"   eps={mp.nstr(eps,2):>9}  delta={mp.nstr(delta,6):>11}  pred={mp.nstr(pred,6):>11}  ratio={mp.nstr(delta/pred,6)}")

# ------------------------------------------------------------------
def section4_tail():
    print("\n"+"="*74)
    print("(5) TAIL  u->-inf (E1 cut):  rho(x) ~ e * e^{-x}   (constant e, MISTAKE-063)")
    print("="*74)
    # On the curve, large x corresponds to u near the negative real axis (E1 cut),
    # where arg(u) -> -pi. findroot is only stable if we APPROACH from small eps,
    # carrying the previous arg as the guess. So trace forward through a fine grid
    # and report the large-x tail rows.
    N = 600
    eps_lo, eps_hi = mp.mpf('0.05'), mp.mpf('9.0')
    eps_list = [eps_lo*(eps_hi/eps_lo)**(mp.mpf(i)/N) for i in range(N+1)]
    rows = trace_curve(eps_list)
    # DERIVED tail law (Stokes term of the closed form):  Im K(w)=0 forces
    #   w_i = -pi e^{-1/w_r},  so  rho = -w_i/pi = e^{-1/w_r}, i.e.
    #        x + ln rho = x - 1/w_r = R(w_r)   (the R-transform value),
    #   and  R(w_r) = 1 + 2 w_r + 6 w_r^2 + ...  -> R(0) = kappa_1 = 1  as x->inf,
    #   so   rho(x) e^x = e^{R(w_r)} -> e.   (Slow: R is the resurgent series at w_r~1/x.)
    # Clean check needs NO Borel: just x+ln(rho) [= R(w_r)] should trend to 1; and
    # rho should equal e^{-1/w_r} with 1/w_r = -ln(rho).  Verify both.
    print("   eps      x         rho          x+ln(rho)   1/w_r=-ln rho   rho*e^x   (->e)")
    rows = [r for r in rows if r[2] > 0]
    # subsample ~14 rows spread over the tail
    tail = rows[len(rows)//2::max(1,len(rows)//28)]
    for (eps,phi,x,rho) in tail:
        Rval = x + mp.log(rho)          # = R(w_r), should -> 1
        invwr = -mp.log(rho)            # = 1/w_r  (since rho = e^{-1/w_r})
        print(f"   {mp.nstr(eps,3):>6}  {mp.nstr(x,4):>7}  {mp.nstr(rho,4):>11}   "
              f"{mp.nstr(Rval,6):>9}   {mp.nstr(invwr,6):>9}     {mp.nstr(rho*mp.e**x,6):>7}")
    # cross-check rho = e^{-1/w_r} against the geometric w_r = Re G:
    print("\n   consistency  rho  vs  e^{-1/Re(G)}  (Stokes prediction):")
    for (eps,phi,x,rho) in tail[-5:]:
        w = -1/(eps*mp.e**(1j*phi)); w_r = mp.re(w)
        print(f"   x={mp.nstr(x,5):>8}  rho={mp.nstr(rho,7):>11}  e^(-1/Re G)={mp.nstr(mp.e**(-1/w_r),7):>11}  "
              f"ratio={mp.nstr(rho/mp.e**(-1/w_r),6)}")
    # analytic tail: x->inf, G~1/x, w=1/x small real>0, u=-1/w=-x large negative.
    # E1(u) for u=-x+i0 (x>0): E1(-x+-i0) = -Ei(x) -+ i pi. Im E1(-x-i0)=+pi.
    # g(u)=e^u E1(u); G = -1/u = 1/x.  rho = -(1/pi) Im G.  But to leading order use
    # the known free-CP tail; here we confirm the constant numerically above.
    print("\n   (analytic) E1 cut: Im E1(-x - i0) = +pi  => spectral weight ~ e^{-x} * const")

if __name__ == "__main__":
    section1_verify_vs_rootfound()
    section2_moments()
    section3_edge()
    section4_tail()
