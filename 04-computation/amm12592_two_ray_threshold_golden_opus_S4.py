"""
AMM 12592: death-star's two-ray entropy threshold gamma* has closed form 2 ln(phi)/ln 5.

Two-ray comparison (THM-3002, lanes C2/F2):
    max_y [ gamma(1+y) H((x-y)/(gamma(1+y))) + (x-y) ln2 ] >= H(x)   for all x in (0,1),
H = natural-log binary entropy. gamma* = smallest gamma for which this holds.

Result:  gamma* = 2 ln(phi)/ln 5 = log_5(phi^2),  phi=(1+sqrt5)/2,  binding x*=(3-sqrt5)/2=1/phi^2.
         C* = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654401497...
Reproduces death-star's amm12592_archimedean_threshold_asymptotic.out (gamma*~0.597987, delta~0.6180=1/phi).
"""
import mpmath as mp
mp.mp.dps = 50
ln = mp.log; sqrt = mp.sqrt; L = ln(2)
phi = (1 + sqrt(5)) / 2

def H(x):
    if x <= 0 or x >= 1: return mp.mpf(0)
    return -x*ln(x) - (1-x)*ln(1-x)

def inner(g, x, iters=200):
    """max over y of the two-ray bracket (golden-section search on valid y-domain)."""
    ylo = max(mp.mpf('-0.999'), (x-g)/(1+g)); yhi = x
    def val(y):
        cap = g*(1+y)
        if cap <= 0: return mp.mpf('-inf')
        fr = (x-y)/cap
        if fr < 0 or fr > 1: return mp.mpf('-inf')
        return cap*H(fr) + (x-y)*L
    a, b = ylo, yhi
    for _ in range(iters):
        m1 = a + (b-a)/3; m2 = b - (b-a)/3
        if val(m1) < val(m2): a = m1
        else: b = m2
    return val((a+b)/2)

# --- 1. raw threshold (independent cross-check; coarse dps for speed) ---
def slack(g):
    return min(inner(g, mp.mpf(i)/120, iters=70) - H(mp.mpf(i)/120) for i in range(1, 120))
_dps = mp.mp.dps; mp.mp.dps = 25
gamma_num = mp.findroot(slack, mp.mpf('0.6'))
mp.mp.dps = _dps

# --- 2. closed form ---
gamma_cf = 2*ln(phi)/ln(5)
xstar = (3 - sqrt(5))/2

# --- 3. the two Lagrangian branches at x* (both must equal gamma*) ---
def P(x): return -2*x*ln(x) - (1-x)*ln(1-x) + (1+x)*ln(1+x)
def gB(x): return ln((1-x)/x) / ln((1+x)/(1-x))                     # (B'')
def gC(x): return (1+x)*H(x) / ((1+x)*P(x) - 2*x*H(x))             # (C'')

print("gamma* (raw two-ray optimizer) =", mp.nstr(gamma_num, 40))
print("gamma* = 2 ln(phi)/ln 5        =", mp.nstr(gamma_cf, 40))
print("  |raw - closed form|          =", mp.nstr(abs(gamma_num - gamma_cf), 4))
print()
print("binding x* = (3-sqrt5)/2 = 1/phi^2 =", mp.nstr(xstar, 30))
print("  death-star delta = 1-x* = 1/phi  =", mp.nstr(1-xstar, 30), "(their 0.6180)")
print("  y-opt frac = 2x*/(1+x*) = 1-1/sqrt5 =", mp.nstr(2*xstar/(1+xstar), 30))
print()
print("branch B'' at x* =", mp.nstr(gB(xstar), 40))
print("branch C'' at x* =", mp.nstr(gC(xstar), 40))
print("both = 2 ln phi/ln5? ",
      abs(gB(xstar)-gamma_cf) < mp.mpf(10)**-45 and abs(gC(xstar)-gamma_cf) < mp.mpf(10)**-45)
print("raw slack at (x*, gamma*) =", mp.nstr(inner(gamma_cf, xstar) - H(xstar), 6), "(~0 confirms saddle)")
print()
print("structural identities at x*:")
print("  ln x*   = -2 ln phi :", mp.nstr(ln(xstar)+2*ln(phi), 4))
print("  ln(1-x*)= -ln phi   :", mp.nstr(ln(1-xstar)+ln(phi), 4))
print("  ln(1+x*)= ln5/2-lnphi:", mp.nstr(ln(1+xstar)-(ln(5)/2-ln(phi)), 4))
print("  H(x*) = (1+x*)ln phi :", mp.nstr(H(xstar)-(1+xstar)*ln(phi), 4))
print()
print("=> C* = 1 + log_5(phi^2) = log_5(5 phi^2) = log_5((15+5 sqrt5)/2)")
print("     =", mp.nstr(1+gamma_cf, 40))
print("   log_5(5 phi^2) =", mp.nstr(ln(5*phi**2)/ln(5), 40), "(same)")
print("   8/5 = 1.6 ruled out; answer is transcendental.")
