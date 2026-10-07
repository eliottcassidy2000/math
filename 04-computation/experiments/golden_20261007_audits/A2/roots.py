# Linear (Haar) model for log-periodic modes: psi(u) = 1/2 psi(u+alpha) + 1/2 psi(u-1), alpha=log2(3/2).
# Modes e^{z u}, z = mu + 2 pi i nu, roots of e^{z alpha} + e^{-z} = 2. Decay per octave = mu (<0).
import cmath, math
L3 = math.log2(3); al = L3 - 1
def g(z): return cmath.exp(z*al) + cmath.exp(-z) - 2
def dg(z): return al*cmath.exp(z*al) - cmath.exp(-z)
def newton(z):
    for _ in range(200):
        dz = g(z)/dg(z); z -= dz
        if abs(dz) < 1e-14: break
    return z
def theta(f):
    x = f*L3; return x - round(x)
obs = {53:0.993, 41:0.917, 12:0.90, 65:0.885, 24:0.82, 29:0.82, 17:0.78, 36:0.72, 5:None,
       106:None, 147:None, 159:None, 200:None, 94:None, 265:None, 253:None, 306:None}
print(" f    theta    gauss(552th^2)  exact_mu   exact_nu   obs_rate(HYP)")
for f in sorted(obs):
    th = theta(f)
    # start Newton from the small-theta approximation
    z0 = complex(-552.7*th*th, 2*math.pi*(f + th/(1-al)))
    z = newton(z0)
    o = obs[f]
    print(f"{f:4d} {th:+.5f}  {552.7*th*th:8.4f}     {z.real:+.4f}   {z.imag/(2*math.pi):9.4f}   {(-math.log(o)) if o else float('nan'):.4f}")
# also: odd-step variance per octave
mu = (al - 1)/2  # mean log change per step
var = ((al+1)/2)**2
print("odd-step variance per octave (renewal):", (1/abs(mu))**2 * 0.25 * (1/abs(mu)) * 0 + ( (1/(1-L3/2)) **2 ) * 0.25 * (1/(1-L3/2)) )
