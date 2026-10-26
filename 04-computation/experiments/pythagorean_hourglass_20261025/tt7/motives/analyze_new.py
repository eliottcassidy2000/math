# Separate the four new seven-vertex motives inside the cheap triple products Y02, Y05, Y09, Y18.
# 1. Weighted least squares (weight q^-5) of each Y on the base columns 1, q, q^2, q^3, a84, q a84, a86, a46, |b|^2;
#    distinct irreducible motives are asymptotically orthogonal in this inner product, so the coefficients should be
#    close to small rationals.  2. Gram matrix of the residuals N_i, normalised by q^5 and the prime count.
#    3. Moments of the normalised residuals.  mac-mini-2026-10-25-pyth.
import os, sys, subprocess, tempfile
import numpy as np
rows = []
for fn in sorted(f for f in os.listdir('.') if f.startswith('cheapY_') and f.endswith('.txt')):
    rows += [list(map(int, l.split())) for l in open(fn)]
rows = sorted({r[0]: r for r in rows}.values())
qmin = int(sys.argv[1]) if len(sys.argv) > 1 else 100
rows = [r for r in rows if r[0] >= qmin]
qs = np.array([r[0] for r in rows], dtype=float)
Y = np.array([[r[k] for k in (1, 3, 5, 6)] for r in rows], dtype=float)     # Y02, Y05, Y09, Y18
P = int(qs.max()) + 2
gpf = os.path.join(tempfile.gettempdir(), 'forms_an.gp')
open(gpf, 'w').write(f'default(parisize,1000000000);\nf=mfeigenbasis(mfinit([8,4,1],0))[1];g=mfeigenbasis(mfinit([8,6,1],0))[1];h=mfeigenbasis(mfinit([4,6,1],0))[1];b=mfeigenbasis(mfinit([32,3,-4],0))[1];\n'
                     f'print(mfcoefs(f,{P}));print(mfcoefs(g,{P}));print(mfcoefs(h,{P}));c=mfcoefs(b,{P});print(vector({P}+1,n,norml2(subst(lift(c[n]),y,I))));\nquit;\n')
out = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()
f84, f86, f46, B2 = (eval(l.replace('I', '1j')) for l in out[-4:])
qi = qs.astype(int)
base = np.array([np.ones_like(qs), qs, qs ** 2, qs ** 3, [f84[q] for q in qi], qs * np.array([f84[q] for q in qi]),
                 [f86[q] for q in qi], [f46[q] for q in qi], [round(abs(B2[q])) for q in qi]], dtype=float).T
w = qs ** -2.5
Bw = base * w[:, None]; Yw = Y * w[:, None]
coef, *_ = np.linalg.lstsq(Bw, Yw, rcond=None)
names = ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46', '|b|^2']
print(f"{len(qs)} primes in [{int(qs.min())}, {int(qs.max())}]")
for i, nm in enumerate(['Y02', 'Y05', 'Y09', 'Y18']):
    print(f"  {nm}: " + ", ".join(f"{n}: {c:+.3f}" for n, c in zip(names, coef[:, i])))
R = (Yw - Bw @ coef)
Q = R.T @ R / len(qs)
np.set_printoptions(precision=3, suppress=True)
print("Gram matrix of the normalised new parts (Y02, Y05, Y09, Y18):\n", Q)
print("eigenvalues:", np.linalg.eigvalsh(Q))
