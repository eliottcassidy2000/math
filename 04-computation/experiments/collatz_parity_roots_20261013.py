"""exact decay roots of the signed transfer psi(u) = -1/2 psi(u + log2(3/2)) - 1/2 psi(u - 1):
mode e^{(2 pi i f + w) u} needs exp(2 pi i theta_f + w log2(3/2)) + exp(-w) = 2 with theta_f = f log2 3 - nearest integer,
f a half-integer (sign flip per step); Newton continuation from w(0) = 0.  Re w = decay per octave (natural log), Im w/2pi = drift."""
import numpy as np
from mpmath import mp, mpc, exp, findroot, log
mp.dps = 30
L = log(3)/log(2) - 1   # log2(3/2)
def root(theta):
    w = mpc(0)
    for t in np.linspace(0, float(theta), 60)[1:]:
        F = lambda w: exp(2j*np.pi*t + w*L) + exp(-w) - 2
        w = findroot(F, w)
    return w
th = log(3)/log(2)
print(' f (cycles/oct)   q*theta   theta_f      Re w (per octave)   factor/octave   Im w/2pi')
for q in (1, 3, 5, 7, 9, 19, 29, 53, 665):
    f = q/2
    qt = float(q*th)
    theta_f = qt/2 - round(qt/2)
    w = root(theta_f)
    print(f'{f:8.1f}   {qt:12.5f}   {theta_f:+.5f}   {float(w.real):+.5f}          {float(exp(w.real)):.4f}         {float(w.imag)/(2*np.pi):+.5f}')
