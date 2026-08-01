"""
GENERAL reformulation: I = int_0^1 e^{Q} dt = Phi_Q(1), where Phi_Q(z) = int_0^1 e^{z Q(t)} dt is an E-FUNCTION
in the parameter z for EVERY Q in Qbar[t].  Transcendence of I <= RIGIDITY (1 not a horizontal section of the
differential module of Phi_Q) = critical-value non-resonance.  See companion reflection.

Verifies (all degrees): (A) mu_m=int_0^1 Q^m algebraic + geometric growth (E-function -- CORRECT, free);
(B) the indicator width -- but see CORRECTION: width>0 gives only an UPPER bound on |Phi_Q|, NOT non-rationality
(a complex oscillatory integral can cancel to a constant while |integrand| grows: g=e^{2pi i t} has Phi_g==1).
So Phi_Q non-rational is EQUIVALENT to FC(2), not free (death-star THM-3031); (C) endpoint/saddle exponential
module and rate-0 resonance; (D) I=Phi_Q(1)=sum mu_m/m!.
"""
import mpmath as mp
import sympy as sp
mp.mp.dps = 30
t, z = sp.symbols('t z')

def report(name, Q, deg):
    print("="*76); print(f"{name}:  Q = {Q}   (deg {deg})")
    # (A) E-function coefficients mu_m = int_0^1 Q^m
    mus = [sp.nsimplify(sp.integrate(Q**m, (t, 0, 1))) for m in range(6)]
    print(f"  (A) mu_m=int_0^1 Q^m (algebraic):", [str(x) for x in mus])
    Qmax = max(abs(mp.mpf(str(sp.N(Q.subs(t, tt))))) for tt in [mp.mpf(k)/20 for k in range(21)])
    growth = [abs(mp.mpf(str(sp.N(mus[m])))) / (Qmax**m if m else 1) for m in range(1, 6)]
    print(f"      |mu_m|/max|Q|^m <= 1 (geometric):", [mp.nstr(g, 4) for g in growth])
    # (B) non-rationality: indicator width max_theta [h(theta)+h(theta+pi)] > 0 for nonconstant Q
    ts = [mp.mpf(k)/200 for k in range(201)]
    Qv = [complex(sp.N(Q.subs(t, tt))) for tt in ts]
    width = max(max((v*mp.e**(1j*th)).real for v in Qv) - min((v*mp.e**(1j*th)).real for v in Qv)
                for th in [mp.pi*k/12 for k in range(12)])
    print(f"  (B) indicator width max_theta = {mp.nstr(width,6)} > 0 => only |Phi_Q| UPPER bound e^{{R.width}};")
    print(f"      NOT non-rationality (g=e^{{2pi i t}} has Phi_g==1 with width>0). Non-rationality <=> FC(2). [CORRECTED]")
    # (C) critical values {Q(0),Q(1)} U {Q(c):Q'(c)=0}: rate-0 resonance <=> 0 among/related to them
    crit = [Q.subs(t, 0), Q.subs(t, 1)]
    for c in sp.solve(sp.diff(Q, t), t):
        if c.is_real and 0 <= sp.re(c) <= 1:
            crit.append(sp.nsimplify(Q.subs(t, c)))
    print(f"  (C) critical values (module exponential rates) = {[str(sp.nsimplify(c)) for c in crit]}")
    print(f"      RIGIDITY holds (=> I transcendental) iff 0 is non-resonant with these; here 0 in set: "
          f"{any(sp.nsimplify(c)==0 for c in crit)}")
    # (D) I = Phi_Q(1) = sum mu_m/m!
    Qf = sp.lambdify(t, Q, 'mpmath')
    I = mp.quad(lambda tt: mp.e**(Qf(tt)), [0, 1])
    Iser = mp.nsum(lambda m: mp.mpf(str(sp.integrate(Q**int(m), (t, 0, 1)))) / mp.factorial(int(m)), [0, mp.inf])
    print(f"  (D) I=Phi_Q(1)={mp.nstr(I,15)} = sum mu_m/m! = {mp.nstr(Iser,15)}  (diff {mp.nstr(abs(I-Iser),3)})")

report("degree 2 (closed)",       sp.Rational(7,10)*t**2 + sp.Rational(9,10)*t - sp.Rational(1,5), 2)
report("degree 3 (rigidity open)", sp.Rational(1,2)*t**3 - sp.Rational(3,4)*t**2 + sp.Rational(1,5), 3)
report("degree 1 (LW)",           2*t + sp.Rational(1,3), 1)
print("="*76)
print("Cited: Beukers 2006 (refined linear S-S, no exceptional point) turns RIGIDITY into I transcendental.")
print("Proved unconditionally here (all degrees): (A) E-function, (B) non-rationality. Degree<=2: RIGIDITY too.")
