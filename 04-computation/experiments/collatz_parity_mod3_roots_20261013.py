"""mod-3 refined transfer roots for the three step-count parities.
Class means m_r(u), r mod 3, u = log2 n: even n = r go to n/2 = 2r mod 3; odd n go to (3n+1)/2 = 2 mod 3 at u + theta - 1.
tau:   m_r = -1/2 m_{2r}(u-1) - 1/2 m_2(u+theta-1)  ->  a^2 + a b - 2 b = 4
omega: m_r = +1/2 m_{2r}(u-1) - 1/2 m_2(u+theta-1)  ->  a^2 - a b - 2 b = 4
C:     m_r = -1/2 m_{2r}(u-1) + 1/2 m_2(u+theta-1)  ->  a^2 - a b + 2 b = 4
with a = e^{-s}, b = e^{s(theta-1)}, s = 2 pi i f + w (classes 1 and 2 form the closed subsystem; class 0 is fed by it).
Scalar model for comparison: tau a + b = -2; omega a - b = 2; C -a + b = 2."""
import numpy as np
from mpmath import mp, mpc, exp, findroot, log, pi as mpi
mp.dps = 30
L = log(3)/log(2) - 1
def eqs(kind):
    if kind == 'tau':
        return (lambda s: exp(-s) + exp(s*L) + 2,            lambda s: exp(-2*s) + exp(-s)*exp(s*L) - 2*exp(s*L) - 4)
    if kind == 'omega':
        return (lambda s: exp(-s) - exp(s*L) - 2,            lambda s: exp(-2*s) - exp(-s)*exp(s*L) - 2*exp(s*L) - 4)
    if kind == 'C':
        return (lambda s: -exp(-s) + exp(s*L) - 2,           lambda s: exp(-2*s) - exp(-s)*exp(s*L) + 2*exp(s*L) - 4)
preds = {'tau': [2.5, 9.5, 14.5, 26.5, 332.5], 'omega': [1, 6, 18, 47, 153], 'C': [1.5, 3.5, 8.5, 20.5, 73.5, 179.5]}
print(' kind    f     scalar: f_eff   factor/oct   |  mod-3: f_eff   factor/oct')
for kind, fs in preds.items():
    e1, e3 = eqs(kind)
    for f in fs:
        s0 = 2j*mpi*f
        out = []
        for e in (e1, e3):
            try:
                w = findroot(lambda w: e(s0 + w), mpc(-0.01, 0.0))
                out.append((f + float(w.imag)/(2*np.pi), float(exp(w.real))))
            except Exception as ex:
                out.append((float('nan'), float('nan')))
        print(f' {kind:5s} {f:7.1f}   {out[0][0]:9.4f}   {out[0][1]:.4f}      |   {out[1][0]:9.4f}   {out[1][1]:.4f}')
