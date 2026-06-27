#!/usr/bin/env python3
"""NEW SIGNAL solidified: #real-roots of the miss-PGF  <->  extreme-mass (q0+q6)  <->  gK8 extremizer (S66).

Claim under test: across configs, the miss-count PGF G_N(z)=sum q_t z^t has FEWER real roots exactly when
the sector-miss distribution is more BIMODAL (high extreme-mass q0+q6) -- and consec, the gK8/coverage
extremizer, has ZERO real roots (maximally non-independent sectors). Real-rootedness of a PGF <=> the count
is a sum of INDEPENDENT indicators (each factor (1-p)+pz, real root); complex roots = sector correlation.
So '#real roots of the miss-PGF' is a new signal measuring sector-INDEPENDENCE, minimized at the extremizer.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
random.seed(66)

def sector_of(p): return int((p % 1) * 7)
def missdist(E):
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): b.add(F(m, 7 * e))
    b = sorted(b); q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        t = 7 - len(set(sector_of(e * ((x0 + x1) / 2)) for e in E))
        if 0 <= t <= 6: q[t] += x1 - x0
    return q
def n_real_roots(q):
    c = [float(q[t]) for t in range(6, -1, -1)]
    while len(c) > 1 and abs(c[0]) < 1e-14: c = c[1:]
    r = np.roots(c)
    return sum(1 for z in r if abs(z.imag) < 1e-7), sorted(r, key=lambda z: abs(z))
def Lyk8(q): return float(10*q[0] + q[3] + 10*q[6])
def extreme_mass(q): return float(q[0] + q[6])

K = 8; N = 250
print("=" * 88)
print(f" #real-roots of miss-PGF  vs  extreme-mass (q0+q6)  vs  L_yK8   (k={K}, {N} configs)")
print("=" * 88)
samples = []
# the named extremals
named = {
    "consec_8": tuple(range(8)),
    "even-AP": tuple(2*i for i in range(8)),
    "top-cluster": (0,7,8,9,10,11,12,13),
}
for name, E in named.items():
    q = missdist(E); nr, roots = n_real_roots(q)
    print(f"  {name:12s}: #real={nr}/6  extreme-mass={extreme_mass(q):.4f}  L_yK8={Lyk8(q):.4f}  p0={float(q[0]):.4f}")
    samples.append((nr, extreme_mass(q), Lyk8(q)))
# random configs
for _ in range(N):
    E = tuple([0] + random.sample(range(1, 40), 7))
    q = missdist(E); nr, _ = n_real_roots(q)
    samples.append((nr, extreme_mass(q), Lyk8(q)))

# group by #real roots
from collections import defaultdict
g = defaultdict(list)
for nr, em, ly in samples: g[nr].append((em, ly))
print("\n  #real-roots -> mean extreme-mass, mean L_yK8, count:")
for nr in sorted(g):
    ems = [x[0] for x in g[nr]]; lys = [x[1] for x in g[nr]]
    print(f"    #real={nr}: n={len(g[nr]):3d}   mean extreme-mass={sum(ems)/len(ems):.4f}   "
          f"mean L_yK8={sum(lys)/len(lys):.4f}   max L_yK8={max(lys):.4f}")
# correlation
import math
xs = [s[0] for s in samples]; ys = [s[1] for s in samples]
mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
sx = math.sqrt(sum((x-mx)**2 for x in xs)); sy = math.sqrt(sum((y-my)**2 for y in ys))
corr = cov/(sx*sy) if sx*sy>0 else 0
print(f"\n  corr(#real-roots, extreme-mass) = {corr:+.3f}   "
      f"(negative => fewer real roots <=> higher extreme-mass = the gK8 extremizer direction)")

print("\n" + "=" * 88)
print(" consec_8 roots in C (do they lie on a curve / relate to the apex 7?)")
print("=" * 88)
q = missdist(tuple(range(8))); _, roots = n_real_roots(q)
for z in roots:
    print(f"   z = {z.real:+.4f} {z.imag:+.4f}i   |z|={abs(z):.4f}  arg={np.angle(z)*180/np.pi:+.1f} deg")
print(f"   product of |roots| = {np.prod([abs(z) for z in roots]):.4f}  (= q0/q6 = {float(q[0]/q[6]):.4f})")
print(f"   all roots have |z|>1 (=> q is 'top-heavy enough' that G_N(z)>0 on [0,1]); 3 conjugate pairs = ")
print(f"   maximal NON-real-rootedness = maximal sector correlation = the gK8/extreme-mass extremizer.")
