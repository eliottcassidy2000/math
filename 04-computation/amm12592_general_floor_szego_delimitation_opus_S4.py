"""
AMM 12592 general-class floor: the analytic/Szego rigidity method is confined to EXACTLY gamma=0.

g(p)=sum p^m W_m(p) is analytic on D_gamma={|p|(|p|+|1-p|)^gamma<1}; with the p<->1-p mirror it continues to
D_gamma U D'_gamma = {min(|p|,|1-p|)(|p|+|1-p|)^gamma<1}. The two-circle point e^{i pi/3} (|p|=|1-p|=1) that
drives the gamma=0 proof has domain value 2^gamma: ON the boundary at gamma=0, OUTSIDE for gamma>0. And the
radius of convergence R(gamma) (min|p| on the singular set, on the negative axis r(1+2r)^gamma=1) is 1 at
gamma=0 and <1 for gamma>0 => coefficients unbounded => Szego inapplicable. See the reflection of the same name.
"""
import mpmath as mp
mp.mp.dps = 25

def domain_value(p, gamma):
    a, b = abs(p), abs(1-p)
    return min(a, b)*(a+b)**gamma

w = mp.exp(1j*mp.pi/3)   # e^{i pi/3} = |p|=1 cap |p-1|=1
print("two-circle point e^{i pi/3}: |p|=%.4f, |1-p|=%.4f" % (abs(w), abs(1-w)))
for g in [mp.mpf('0'), mp.mpf('0.1'), mp.mpf('0.3'), mp.mpf('0.59799'), mp.mpf('1')]:
    v = domain_value(w, g)
    where = ("ON boundary -> rigidity applies" if abs(v-1) < 1e-12
             else "INSIDE" if v < 1 else "OUTSIDE -> no analytic control")
    print("  gamma=%.5f: domain value = %s  (%s)" % (float(g), mp.nstr(v, 8), where))

print("\nradius of convergence R(gamma)  (r(1+2r)^gamma = 1 on the negative axis):")
def R(gamma):
    if gamma == 0: return mp.mpf(1)
    return mp.findroot(lambda r: r*(1+2*r)**gamma - 1, mp.mpf('0.9'))
for g in [mp.mpf('0'), mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('1')]:
    print("  gamma=%.2f: R = %s   (<1 for gamma>0 => coefficients unbounded => Szego inapplicable)"
          % (float(g), mp.nstr(R(g), 8)))

print("\nCONCLUSION: Carlson-Szego-two-circle proves EXACTLY gamma=0 impossible (C*>1, general) and,")
print("provably, nothing more. A quantitative general floor must be archimedean, not analytic-rigidity.")

# resummation sanity: von Neumann sum_{label 1} p^h q^t = 1/2 for all p (each term ->0 as p->0)
def vn(p, terms=4000):
    q = 1-p
    return sum((p*p+q*q)**k*(p*q) for k in range(terms))
for p in [mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('0.01')]:
    print("  von Neumann sum at p=%.2f: %s  (=1/2; resummation is essential -- truncation/limits invalid)"
          % (float(p), mp.nstr(vn(p), 12)))
