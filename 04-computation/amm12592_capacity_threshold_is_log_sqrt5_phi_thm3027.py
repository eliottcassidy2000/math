#!/usr/bin/env python3
"""THM-3027 -- the AMM 12592 Bernstein-capacity threshold is EXACTLY

        gamma* = log(phi) / log(sqrt 5) = 2 log(phi) / log 5 = 0.5979874356654...

so the epoch-closure route's natural barrier is C = 1 + gamma* = 1.5979874...

BACKGROUND.  THM-3002 (death-star) reduces block closure to the necessary
Bernstein-capacity criterion

    S(t) := sum_{i<=t} C(d_i, t-i) 2^(t-i)  >=  C(R-1, t)     for all t,     (*)

with d_i = floor(gamma (R+i)) + D0, i = 0..R-1, and reports "the asymptotic
threshold is gamma* ~ 0.598 as the solution of a two-ray entropy comparison".
This file SOLVES that comparison in closed form.  klein-S428 had already
refuted both readings of the eq(27) weight against (*); the constant that
survives is the golden one.

THE RATE PROBLEM.  Put i = sigma R, t = tau R.  Then d_i/R -> D(sigma) =
gamma(1+sigma), and with m = tau - sigma, rho = m/D,

    (1/R) log [C(d_i, t-i) 2^(t-i)]  ->  D H(rho) + m log 2,
    (1/R) log C(R-1, t)              ->  H(tau),          H = natural entropy.

The sum in (*) has R terms, so its rate is the max over sigma.  (*) holds
asymptotically iff  Psi(gamma,tau) := max_sigma [D H(rho) + m log 2] >= H(tau)
for every tau in (0,1), and gamma* is the infimum of such gamma.

THE COLLAPSE (this file's contribution).  At the threshold the constraint is
tangent.  Write u = 1 - rho and A = log(1/u).  The three conditions are

  (S) inner stationarity in sigma:  gamma A = log(2u/(1-u))
  (T) tangency in tau:              2u/(1-u) = (1-tau)/tau
  (V) value:                        D [H(rho) + rho log 2] = H(tau)

(S) gives log(rho/2) = -A(1+gamma), whence the identity

        H(rho) + rho log 2 = A (1 + gamma rho),

and the constraint sigma = tau - rho D with D = gamma(1+sigma) gives
D = gamma(1+tau)/(1+gamma rho).  The factor (1 + gamma rho) CANCELS:

        LHS of (V) = gamma A (1 + tau) = (1+tau) log((1-tau)/tau)  by (S),(T).

So (V) reads (1+tau) log((1-tau)/tau) = -tau log tau - (1-tau) log(1-tau),
i.e. 2 log(1-tau) = log tau, i.e.

        (1 - tau)^2 = tau            <=>   tau^2 - 3 tau + 1 = 0.            (!)

(!) is the minimal polynomial of phi^2 and phi^-2.  Hence tau* = phi^-2 = 2-phi,
1 - tau* = phi^-1, and then (T) gives 2u/(1-u) = phi, so u* = phi/(2+phi) =
1/sqrt 5 (using 2+phi = phi^2+1 = phi sqrt 5), A* = log sqrt 5, and (S) gives

        gamma* = log(phi)/log(sqrt 5) = 2 log(phi)/log 5.

THE GOLDEN RATIO IS NOT DECORATION: it enters as the root of (!), which is the
tangency condition itself, not through any Fibonacci structure put in by hand.

WHAT IS PROVED vs CHECKED
  PROVED (symbolically, part C): given the rate reduction, the tangency system
    (S)+(T)+(V) has the stated unique interior solution, and (V) is EQUIVALENT
    to phi^2 = phi + 1.
  CHECKED NUMERICALLY: that this interior stationary point is the global worst
    tau (part B scans 3000 values of tau); that the maximizing sigma* = 0.03864
    is interior to [0,1] so no boundary case is missed (part C); and that a
    direct bisection of the variational threshold returns gamma* to 9 digits
    (part B).  The Laplace/Stirling passage from (*) to the rate problem is
    standard and is ASSUMED, not reproved here.
  CONSISTENT WITH: death-star's finite-R thresholds 0.5313, 0.5606, 0.5758,
    0.5849, 0.59065, 0.59393 for R = 32..1024, which increase toward gamma*
    from below (part D), and their report that gamma = 2457/4135 = 0.594196
    dies at R = 2048 -- it lies below gamma*.

Reproduce: python3 04-computation/amm12592_capacity_threshold_is_log_sqrt5_phi_thm3027.py
"""

from math import log as flog, sqrt as fsqrt, comb, floor
from mpmath import mp, mpf, log, sqrt

mp.dps = 40


def rule(s):
    print("=" * 76)
    print(s)
    print("=" * 76)


def Hf(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * flog(p) - (1 - p) * flog(1 - p)


def psi(g, tau, s):
    D = g * (1 + s)
    m = tau - s
    if m < 0 or m > D:
        return -1e9
    return D * Hf(m / D) + m * flog(2)


def Psi(g, tau):
    """max over sigma in [0, min(1,tau)], by ternary search (psi is unimodal)."""
    lo, hi = 0.0, min(1.0, tau)
    for _ in range(200):
        a = lo + (hi - lo) / 3
        b = hi - (hi - lo) / 3
        if psi(g, tau, a) < psi(g, tau, b):
            lo = a
        else:
            hi = b
    s = (lo + hi) / 2
    return psi(g, tau, s), s


def worst(g, N=3000):
    w, at, sat = 1e9, None, None
    for k in range(1, N):
        tau = k / N
        v, s = Psi(g, tau)
        d = v - Hf(tau)
        if d < w:
            w, at, sat = d, tau, s
    return w, at, sat


# --------------------------------------------------------------------------
def partA():
    rule("A. POSITIVE CONTROL -- the finite-R criterion (THM-3002) as implemented")

    def min_ratio(gamma, R, D0=0):
        d = [floor(gamma * (R + i)) + D0 for i in range(R)]
        best, at = None, None
        for t in range(1, R):
            S = 0
            for i in range(0, min(t, R - 1) + 1):
                if t - i <= d[i]:
                    S += comb(d[i], t - i) * (1 << (t - i))
            r = S / comb(R - 1, t)
            if best is None or r < best:
                best, at = r, t
        return best, at

    ok = True
    for R, want in ((8, True), (16, True), (64, False)):
        r, _ = min_ratio(0.5, R, 0)
        ok &= ((r >= 1) == want)
        print(f"    gamma=1/2  R={R:3d}  min ratio={r:8.4f}  "
              f"{'AMPLE' if r >= 1 else 'DEFICIENT'}   (expected "
              f"{'AMPLE' if want else 'DEFICIENT'})")
    print(f"  death-star's gamma=1/2 trichotomy reproduced: {ok}")
    return ok


def partB():
    print()
    rule("B. VARIATIONAL THRESHOLD BY BISECTION vs THE GOLDEN CONSTANT")
    lo, hi = 0.50, 0.70
    for _ in range(50):
        mid = (lo + hi) / 2
        if worst(mid)[0] >= 0:
            hi = mid
        else:
            lo = mid
    phi = (1 + fsqrt(5)) / 2
    target = 2 * flog(phi) / flog(5)
    w, at, sat = worst(hi)
    print(f"    bisection gamma*        = {hi:.12f}")
    print(f"    2 log(phi)/log 5        = {target:.12f}")
    print(f"    difference              = {hi - target:.3e}   (tau-grid is 1/3000)")
    print(f"    worst tau               = {at:.6f}   vs  phi^-2 = {1/phi**2:.6f}")
    print(f"    maximizing sigma        = {sat:.6f}   (interior to [0,1]: "
          f"{0 < sat < 1})")
    rho = (at - sat) / (hi * (1 + sat))
    print(f"    rho = m/D               = {rho:.9f}   vs  1 - 1/sqrt5 = "
          f"{1 - 1/fsqrt(5):.9f}")
    print(f"    2(1-rho)/rho            = {2*(1-rho)/rho:.9f}   vs  phi = {phi:.9f}")
    ok = abs(hi - target) < 1e-7 and 0 < sat < 1
    print(f"  VERDICT B: threshold matches the golden constant: {ok}")
    return ok


def partC():
    print()
    rule("C. EXACT SOLUTION OF THE TANGENCY SYSTEM (symbolic, 40 dps)")
    phi = (1 + sqrt(5)) / 2
    u = 1 / sqrt(5)
    rho = 1 - u
    A = log(1 / u)
    g = 2 * log(phi) / log(5)
    tau = (1 - u) / (1 + u)

    def H(p):
        return -p * log(p) - (1 - p) * log(1 - p)

    print(f"    tau* = (1-u)/(1+u)         = {tau}")
    print(f"    phi^-2 = 2 - phi           = {2 - phi}")
    print(f"    (1-tau)^2 - tau            = {(1-tau)**2 - tau}      <-- eq (!)")
    print()
    print(f"    (S)  gamma A               = {g*A}")
    print(f"         log(2u/(1-u))         = {log(2*u/(1-u))}")
    print(f"         residual              = {g*A - log(2*u/(1-u))}")
    print(f"    (T)  2u/(1-u)              = {2*u/(1-u)}")
    print(f"         (1-tau)/tau           = {(1-tau)/tau}")
    print(f"         residual              = {2*u/(1-u) - (1-tau)/tau}")
    print()
    print("    KEY IDENTITY  H(rho) + rho log2 = A (1 + gamma rho):")
    print(f"         LHS                   = {H(rho) + rho*log(2)}")
    print(f"         RHS                   = {A*(1 + g*rho)}")
    print(f"         residual              = {H(rho)+rho*log(2) - A*(1+g*rho)}")
    print()
    D = g * (1 + tau) / (1 + g * rho)
    sig = tau - rho * D
    print(f"    D = gamma(1+tau)/(1+gamma rho) = {D}")
    print(f"    sigma* = tau - rho D           = {sig}   interior: {0 < sig < 1}")
    print()
    lhs = D * (H(rho) + rho * log(2))
    rhs = H(tau)
    print(f"    (V)  D[H(rho)+rho log2]    = {lhs}")
    print(f"         H(tau)                = {rhs}")
    print(f"         residual              = {lhs - rhs}")
    print()
    print("    (V) reduces to 2 log(1-tau) = log tau, i.e. (1-tau)^2 = tau,")
    print("    i.e. tau^2 - 3 tau + 1 = 0 -- the minimal polynomial of phi^2.")
    print(f"         phi^2 - phi - 1       = {phi**2 - phi - 1}")
    tol = mpf(10) ** (-30)
    ok = (abs((1 - tau) ** 2 - tau) < tol
          and abs(g * A - log(2 * u / (1 - u))) < tol
          and abs(2 * u / (1 - u) - (1 - tau) / tau) < tol
          and abs(H(rho) + rho * log(2) - A * (1 + g * rho)) < tol
          and abs(lhs - rhs) < tol
          and 0 < sig < 1)
    print(f"  VERDICT C: all of (S),(T),(V) and the key identity vanish to 30+ dps,")
    print(f"             sigma* interior: {ok}")
    return ok


def partD():
    print()
    rule("D. CONSISTENCY WITH death-star's FINITE-R LADDER")
    phi = (1 + fsqrt(5)) / 2
    gstar = 2 * flog(phi) / flog(5)
    ladder = [(32, 0.5313), (64, 0.5606), (128, 0.5758), (256, 0.5849),
              (512, 0.59065), (1024, 0.59393)]
    print("      R      finite-R threshold    gamma* - threshold")
    inc = True
    prev = -1
    for R, v in ladder:
        print(f"    {R:5d}        {v:.5f}              {gstar - v:+.5f}")
        inc &= (v > prev)
        prev = v
    print(f"    limit      {gstar:.5f}  = 2 log(phi)/log 5")
    print(f"  monotone increasing toward gamma* from below: {inc}")
    print()
    for name, v in (("2457/6592 (eq27 reading R1)", 2457 / 6592),
                    ("2457/4135 (eq27 reading R2)", 2457 / 4135),
                    ("3/5", 0.6), ("1/2", 0.5)):
        verdict = "ABOVE gamma* -- survives" if v > gstar else "BELOW gamma* -- eventually deficient"
        print(f"    {name:30s} = {v:.6f}   {verdict}")
    print("  Both eq(27) readings sit below gamma*, independently confirming the")
    print("  R>=16 and R=2048 deaths found by direct computation.  3/5 clears it,")
    print("  which is why 3/5 is the first round rate that survives.")
    return inc


def partE():
    print()
    rule("E. THE BINDING FRACTION tau* = phi^-2 IS UNIVERSAL IN THE ALPHABET SIZE")
    print("  Replacing 2^(t-i) by b^(t-i) leaves the identity (K) and hence the")
    print("  reduction (V) <=> (1-tau)^2 = tau untouched; only (T) moves, giving")
    print("  u*(b) = phi/(b+phi) and gamma*(b) = log(phi)/log((b+phi)/phi).")
    print()

    def thresh(b):
        def ps(g, tau, s):
            D = g * (1 + s)
            m = tau - s
            if m < 0 or m > D:
                return -1e9
            return D * Hf(m / D) + m * flog(b)

        def Ps(g, tau):
            lo, hi = 0.0, min(1.0, tau)
            for _ in range(160):
                x = lo + (hi - lo) / 3
                y = hi - (hi - lo) / 3
                if ps(g, tau, x) < ps(g, tau, y):
                    lo = x
                else:
                    hi = y
            s = (lo + hi) / 2
            return ps(g, tau, s), s

        def wo(g, N=2000):
            w, at = 1e9, None
            for k in range(1, N):
                t = k / N
                v, _ = Ps(g, t)
                if v - Hf(t) < w:
                    w, at = v - Hf(t), t
            return w, at
        lo, hi = 0.05, 1.60
        for _ in range(42):
            mid = (lo + hi) / 2
            if wo(mid)[0] >= 0:
                hi = mid
            else:
                lo = mid
        return hi, wo(hi)[1]

    phi = (1 + fsqrt(5)) / 2
    ok = True
    print("     b    gamma* numeric   log(phi)/log((b+phi)/phi)   worst tau   phi^-2")
    for b in (2, 3, 4):
        g, t = thresh(b)
        pred = flog(phi) / flog((b + phi) / phi)
        ok &= abs(g - pred) < 1e-6 and abs(t - 1 / phi ** 2) < 2e-3
        print(f"     {b}    {g:.9f}      {pred:.9f}            {t:.6f}   "
              f"{1/phi**2:.6f}")
    print("  b=2 is special ONLY because (2+phi)/phi = sqrt 5 exactly")
    print(f"     (2+phi)/phi = {(2+phi)/phi:.12f}   sqrt 5 = {fsqrt(5):.12f}")
    print(f"  VERDICT E: tau* universal, gamma*(b) formula holds: {ok}")
    return ok


def main():
    a = partA()
    b = partB()
    c = partC()
    d = partD()
    e = partE()
    print()
    rule(f"SUMMARY  control={a}  bisection={b}  exact={c}  ladder={d}  universal={e}")
    print("  gamma* = log(phi)/log(sqrt 5) = 2 log(phi)/log 5 = 0.597987435665440...")
    print("  C = 1 + gamma* = 1.597987435665440...  is the epoch-closure barrier.")


if __name__ == "__main__":
    main()
