"""
collatz_rational_approximants_20261013.py -- the Collatz transfer function as a rational function on P^1 x P^1 restricted
to a winding curve, and its rational approximants along the convergents of theta = log2 3 (S29; owner theme: rational
functions on the Riemann sphere, the constant function infinity).

The scalar law of THM-4629 for the character (alpha, beta) mod m has resolvent R(s) = 1/(1 - x - y) with
x = (1/2) zeta^alpha e^(-s), y = (1/2) zeta^(alpha+beta) e^(s(theta-1)): the rational function 1/(1 - x - y) on P^1 x P^1
(polar line x + y = 1) restricted to the curve Gamma_theta.  For a rational theta' = P/Q the curve is algebraic and
parametrised by w = e^(-s/Q): the restriction is the rational function of w
    R_(P/Q)(w) = w^(P-Q) / p(w),   p(w) = w^(P-Q) - (1/2) zeta^alpha w^P - (1/2) zeta^(alpha+beta),
of degree P on the Riemann sphere, with poles at the P roots of p.  Statement: if Q = alpha and P = -beta (mod m) then
w_0 = zeta^(-1) is a root on the unit circle (an exactly undamped mode), and the exact law moves it to w_0 e^(-w_THM/Q).
Each root w of p is a mode with s = -Q log w: frequency f = -Q arg(w)/(2 pi) (mod Q) and decay per octave e^(Re s)
= |w|^(-Q)... note e^(-s) = w^Q so Re s = -Q ln|w|.
Part 1: for the three parities and the convergents (5,8), (12,19), (41,65), (53,84), (94,149): all roots of p, the ones
within 1.3 of the unit circle (in |w|^Q), their (f, decay), against the exact roots of THM-4629 at the same f (continued
with the true theta).  Part 2: the unfolding: the approximant's single pole at w_0 stands for every multiple (jQ, jP);
the exact law gives the lines f = jQ/m with defect j theta_f and the decays, against the measured census decays.
Part 3: the Moebius dictionary: for each near-return (a, b) = (P - Q halvings? no: a = P halvings... the word has b = Q odd
steps and a = P - Q halvings, total P steps), the multiplier lambda = 3^b/2^(a+b) = 3^Q/2^P, ln lambda = m theta_f ln 2,
and the decay formula exp(-4 pi^2 theta_f^2/(2-theta)^3) = exp(-4 pi^2 (ln lambda)^2/(m^2 (ln 2)^2 (2-theta)^3)) ... printed.
Part 4: the constant function infinity: |R_theta(2 pi i Q/m)| ~ 1/(pi |theta_f|) along the convergents (unbounded on the
line Re s = 0 though finite at every point), against the approximants' poles.
Part 5: Codex's Conway constant: C_* = max phi_N(z)/(2(1 + Phi(z))) recomputed.
"""
import numpy as np
from mpmath import mp, mpf, mpc, exp, findroot, log, pi as mpi, sqrt, erf, quad, nstr

mp.dps = 30
TH = float(log(3) / log(2)); L = log(3) / log(2) - 1
chars = {"tau": (1, 0), "omega": (0, 1), "C": (1, 1)}
convs = [(5, 8), (12, 19), (41, 65), (53, 84), (94, 149)]

def exact_root(zeta_a, zeta_ab, f):
    s0 = 2j * mpi * f
    w = findroot(lambda w: zeta_a * exp(-(s0 + w)) + zeta_ab * exp((s0 + w) * L) - 2, mpc(-0.01, 0))
    return float(w.real), float(w.imag) / (2 * np.pi)

print("# Part 1: poles of the rational approximants R_(P/Q)(w) on the Riemann sphere, against the exact lines")
for name, (a, b) in chars.items():
    za, zab = (-1) ** a, (-1) ** (a + b)
    for Q, P in convs:
        hit = (Q % 2 == a) and (P % 2 == (-b) % 2)
        coeffs = np.zeros(P + 1, dtype=complex)        # polynomial in w: w^(P-Q) - (za/2) w^P - zab/2
        coeffs[P] = -za / 2; coeffs[P - Q] = 1.0; coeffs[0] = -zab / 2
        roots = np.polynomial.polynomial.polyroots(coeffs)
        near = [w for w in roots if abs(np.log(abs(w)) * Q) < 1.3]
        near.sort(key=lambda w: abs(np.log(abs(w)) * Q))
        print(f"\n## {name} parity, convergent (Q,P) = ({Q},{P}), exact root at w = zeta^-1 = -1 expected: {'YES' if hit else 'no'} | p(-1) = {np.polynomial.polynomial.polyval(-1.0, coeffs):.2e}")
        print("   approximant poles with |Re s| < 1.3:  f (mod Q, in [-Q/2, Q/2])   decay/oct   |  exact law at that f: f_eff   decay/oct")
        for w in near[:8]:
            f = -Q * np.angle(w) / (2 * np.pi); res = -Q * np.log(abs(w))
            # the exact law at the lattice frequency nearest to f for this character: f must be = alpha/2 mod 1 (m = 2)
            f_lat = round(f - a / 2) + a / 2 if abs(f) > 1e-9 else f
            try:
                rw, drift = exact_root(za, zab, abs(f_lat) if f_lat != 0 else 0.5)
                ex = f"{abs(f_lat) + drift:9.4f}   {np.exp(rw):.4f}"
            except Exception:
                ex = "   (no continuation)"
            print(f"      {f:+9.4f}                     {np.exp(res):.4f}   |  {ex}")

print("\n# Part 2: unfolding of the Mercator pole: the exact lines at the multiples j (53, 84), defect j theta_f, and their decays")
tf = (53 * TH - 84) / 2
for j in (1, 3, 5, 7, 9):
    Q, P = 53 * j, 84 * j
    rw, drift = exact_root(-1, -1, Q / 2)
    print(f"   j={j}: (Q,P)=({Q},{P}) f = {Q/2:6.1f} -> {Q/2 + drift:9.4f}, theta_f = {j*tf:+.5f}, decay/oct {np.exp(rw):.4f}, quadratic law {np.exp(-4*np.pi**2*(j*tf)**2/(2-TH)**3):.4f}")
print("   (measured decays from collatz_parity_census_m2_m3_20261013.out: 26.5: 0.9997, 79.5: 0.9892, 132.5: 0.9744 fit of the census, 185.5: weak)")

print("\n# Part 3: Moebius dictionary: the word with b = Q odd steps and a = P - Q halvings is the affine map z -> (3^b z + c)/2^P with multiplier 3^Q/2^P")
for Q, P in [(5, 8), (12, 19), (17, 27), (29, 46), (41, 65), (53, 84), (94, 149), (147, 233), (306, 485), (359, 569), (665, 1054)]:
    lam = mpf(3) ** Q / mpf(2) ** P
    lnlam = float(log(lam))
    d = Q * TH - P
    print(f"   (Q,P)=({Q:4d},{P:5d}): multiplier 3^Q/2^P = {float(lam):.6f}, ln multiplier = {lnlam:+.6f} = (Q theta - P) ln 2 = {d*np.log(2):+.6f}; parabolicity defect |ln lambda| = {abs(lnlam):.2e}; fixed point of the word x_w = c_w/(2^P - 3^Q) with 2^P - 3^Q = {nstr(mpf(2)**P - mpf(3)**Q, 4)}")

print("\n# Part 4: the exact transfer function on Re s = 0 at the convergents' frequencies: |R| = 1/|1 - x - y| ~ 1/(pi |theta_f|) (unbounded along the convergents, finite everywhere)")
for Q, P in [(5, 8), (29, 46), (53, 84), (159, 252), (665, 1054)]:
    tf = (Q * TH - P) / 2
    s = 2j * np.pi * Q / 2
    val = 1 / (1 - 0.5 * (-1) * np.exp(-s) - 0.5 * (-1) * np.exp(s * (TH - 1)))
    print(f"   tau parity, f = Q/2 = {Q/2:6.1f}: theta_f = {tf:+.5f}, |R(2 pi i f)| = {abs(val):10.2f}, 1/(pi |theta_f|) = {1/(np.pi*abs(tf)):10.2f}")

print("\n# Part 5: Codex's Conway constant C_* = max_z phi_N(z)/(2 (1 + Phi(z))) recomputed")
phi = lambda z: exp(-z * z / 2) / sqrt(2 * mpi)
Phi = lambda z: (1 + erf(z / sqrt(2))) / 2
F = lambda z: z * (1 + Phi(z)) + phi(z)
zs = findroot(F, mpf(-0.3))
Cs = phi(zs) / (2 * (1 + Phi(zs)))
print(f"   z_* = {zs}, C_* = {Cs}, x_* = 1 + Phi(z_*) = {1 + Phi(zs)}; Codex: z_* = -0.2760298047981433, C_* = 0.13801490239907166, x_* = 1.3912625871346345")
print(f"   measured A_k sqrt(k) at k = 10, 20, 29: 0.1584, 0.1461, 0.1438 vs C_* + c/k with c = {(0.1438 - float(Cs)) * 29:.3f} (k = 29), {(0.1461 - float(Cs)) * 20:.3f} (k = 20), {(0.1584 - float(Cs)) * 10:.3f} (k = 10)")


print("\n# Part 4b: the maximum of |R| on the line Re s = 0 near each convergent (tau parity): location and height against")
print("#   f_peak = Q/2 + theta_f/(2 - theta) (the drifted line frequency of THM-4629), height (2 - theta)^2/(2 pi^2 theta_f^2),")
print("#   the pole formula 2/((2 - theta) |Re w_THM|), and the proved line bound theta^2/(4 theta_f^2)")
def Rabs(f):
    s = 2j * np.pi * f
    return 1 / abs(1 + 0.5 * np.exp(-s) + 0.5 * np.exp(s * (TH - 1)))
for Q, P in [(5, 8), (29, 46), (53, 84), (159, 252), (665, 1054)]:
    tf = (Q * TH - P) / 2
    grid = np.linspace(Q / 2 - 0.08, Q / 2 + 0.08, 400001)
    vals = np.array([Rabs(f) for f in grid])
    i = int(np.argmax(vals)); fpk = grid[i]; hpk = vals[i]
    rw, drift = exact_root(-1, -1, Q / 2)
    print(f"   (Q,P)=({Q},{P}): theta_f = {tf:+.6f}; |R(Q/2)| = {Rabs(Q/2):9.2f}; max on line {hpk:10.2f} at f = {fpk:.5f} (Q/2 {fpk - Q/2:+.5f});"
          f" predicted f_peak shift theta_f/(2-theta) = {tf/(2-TH):+.5f}, THM-4629 drift {drift:+.5f}; height (2-theta)^2/(2 pi^2 theta_f^2) = {(2-TH)**2/(2*np.pi**2*tf**2):10.2f};"
          f" 2/((2-theta)|Re w|) = {2/((2-TH)*abs(rw)):10.2f}; bound theta^2/(4 theta_f^2) = {TH**2/(4*tf**2):10.2f}")

print("\n# Part 4c: the proved line bound |R(2 pi i f)| <= pi^2/(phi^2 + psi^2) <= theta^2/(4 theta_f(f)^2) on a dense grid (tau parity, m = 2):")
print("#   phi = principal 2 pi (1/2 - f), psi = principal 2 pi (1/2 + f (theta - 1)), theta_f(f) := (psi + (theta - 1) phi)/(2 pi) = (Q theta - P)/2 for the integers Q, P (Q odd, P even) fixed by f")
def princ(t): return (t + np.pi) % (2 * np.pi) - np.pi
fs = np.linspace(0.0005, 400, 4000001)
phi_ = princ(2 * np.pi * (0.5 - fs)); psi_ = princ(2 * np.pi * (0.5 + fs * (TH - 1)))
R_ = 1 / np.abs(1 + 0.5 * np.exp(-2j * np.pi * fs) + 0.5 * np.exp(2j * np.pi * fs * (TH - 1)))
b1 = np.pi ** 2 / (phi_ ** 2 + psi_ ** 2); tf_ = (psi_ + (TH - 1) * phi_) / (2 * np.pi); b2 = TH ** 2 / (4 * tf_ ** 2)
print(f"   grid of {len(fs)} frequencies in (0, 400]: max |R| = {R_.max():.2f} at f = {fs[np.argmax(R_)]:.4f}; violations of |R| <= pi^2/(phi^2+psi^2): {int(np.sum(R_ > b1 * (1 + 1e-9)))}; of pi^2/(phi^2+psi^2) <= theta^2/(4 theta_f^2): {int(np.sum(b1 > b2 * (1 + 1e-9)))}")
print(f"   min over the grid of (bound theta^2/(4 theta_f^2)) / |R| = {(b2 / R_).min():.3f} (at f = {fs[np.argmin(b2 / R_)]:.4f})")


print("\n# Part 4d: |R| at the THM-4629 drifted frequency Q/2 + drift, and a two-stage local search (window 0.01 at spacing 1e-7, then 2e-5 at 1e-10) around it")
for Q, P in [(5, 8), (29, 46), (53, 84), (159, 252), (665, 1054)]:
    rw, drift = exact_root(-1, -1, Q / 2)
    f0 = Q / 2 + drift
    loc1 = np.linspace(f0 - 0.01, f0 + 0.01, 200001)
    lv1 = np.array([Rabs(f) for f in loc1]); f1 = loc1[int(np.argmax(lv1))]
    loc = np.linspace(f1 - 2e-5, f1 + 2e-5, 400001)
    lv = np.array([Rabs(f) for f in loc]); i = int(np.argmax(lv))
    print(f"   (Q,P)=({Q},{P}): |R(Q/2 + drift)| = {Rabs(f0):12.2f}; local max {lv[i]:12.2f} at f - (Q/2 + drift) = {loc[i] - f0:+.2e}; 2/((2-theta)|Re w|) = {2/((2-TH)*abs(rw)):12.2f}; ratio local max / pole formula = {lv[i]/(2/((2-TH)*abs(rw))):.4f}")
print("done")
