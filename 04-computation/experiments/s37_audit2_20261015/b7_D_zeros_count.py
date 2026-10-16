#!/usr/bin/env python3
"""S37 audit 2, script b7: zeros of D(s) = 1 + 2^(-s) + (3/2)^s/3 in the strip 0 < Im s < 3000.

 (1) Argument-principle counts (winding number of D along rectangle boundaries, sampled adaptively so that every
     phase increment is below 0.3 rad, minimum |D| on each boundary reported) for the strips
     1/2 < Re s < 1, 0.65 < Re s < 1, 0.85 < Re s < 1 and 1 < Re s < 3.01 (D has no zeros with Re s >= 3:
     |(3/2)^s|/3 > 1 + 2^(-Re s) there), split into rectangles of height 10.
 (2) Every zero located (Newton from the rectangle centres plus a scan; each rectangle's count matched), listed and
     compared with the 28 zeros of collatz_mercator_waveform_prediction_20261015.json.
 (3) Simplicity: D and D' have no common zero (closed-form argument), and min |D'| over the zeros found.
"""
import sys, json, math, cmath
import numpy as np
import mpmath as mp

PRED = sys.argv[1]
L2, L15 = math.log(2), math.log(1.5)


def D(s):
    return 1 + np.exp(-s * L2) + np.exp(s * L15) / 3


def Dp(s):
    return -L2 * np.exp(-s * L2) + L15 * np.exp(s * L15) / 3


def seg_winding(a, b, n0=64):
    # phase change of D along the segment a -> b with adaptive refinement
    ts = np.linspace(0, 1, n0 + 1)
    total = 0.0
    minabs = np.inf
    stack = [(0.0, 1.0)]
    pts = a + (b - a) * ts
    vals = D(pts)
    minabs = np.abs(vals).min()
    dph = np.angle(vals[1:] / vals[:-1])
    for i in range(n0):
        if abs(dph[i]) < 0.3:
            total += dph[i]
        else:
            # refine this piece recursively
            w, m = seg_winding(pts[i], pts[i + 1], 16)
            total += w
            minabs = min(minabs, m)
    return total, minabs


def rect_count(s0, s1, t0, t1):
    c = [complex(s0, t0), complex(s1, t0), complex(s1, t1), complex(s0, t1)]
    tot = 0.0
    mn = np.inf
    for k in range(4):
        a, b = c[k], c[(k + 1) % 4]
        n0 = max(64, int(abs(b - a) * 40))
        w, m = seg_winding(a, b, n0)
        tot += w
        mn = min(mn, m)
    return tot / (2 * math.pi), mn


strips = [(0.5, 1.0), (0.65, 1.0), (0.85, 1.0), (1.0, 3.01)]
H = 10.0
counts = {}
rects = {}
for (s0, s1) in strips:
    tot = 0
    worst_int = 0.0
    minb = np.inf
    per = []
    for k in range(300):
        t0, t1 = k * H + (1e-9 if k == 0 else 0.0), (k + 1) * H
        w, mn = rect_count(s0, s1, t0, t1)
        worst_int = max(worst_int, abs(w - round(w)))
        minb = min(minb, mn)
        tot += round(w)
        per.append(round(w))
    counts[(s0, s1)] = tot
    rects[(s0, s1)] = per
    print(f'strip {s0} < Re s < {s1}, 0 < Im s < 3000: {tot} zeros (max distance of a winding number from an integer '
          f'{worst_int:.1e}; min |D| on the rectangle boundaries {minb:.3e})')

# locate zeros: Newton from a dense set of starts, keep distinct roots, then check against the counts
mp.mp.dps = 30
Dm = lambda s: 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3
found = []
for t in np.arange(0.5, 3000.0, 0.5):
    for sg in (0.55, 0.7, 0.85, 0.95, 1.0, 1.3, 1.8, 2.5):
        s = complex(sg, t)
        for _ in range(60):
            ds = D(s) / Dp(s)
            s -= ds
            if abs(ds) < 1e-14:
                break
        if abs(D(s)) < 1e-10 and 0 < s.imag < 3000 and 0.5 < s.real < 3.01:
            if all(abs(s - z) > 1e-7 for z in found):
                found.append(s)
found.sort(key=lambda z: z.imag)
for (s0, s1) in strips:
    n = sum(1 for z in found if s0 < z.real < s1)
    print(f'  zeros located by Newton in {s0} < Re < {s1}: {n} (argument principle {counts[(s0, s1)]})')
# per-rectangle agreement in the main strip
mism = 0
for (s0, s1) in strips:
    for k in range(300):
        n = sum(1 for z in found if s0 < z.real < s1 and k * H < z.imag <= (k + 1) * H)
        mism += n != rects[(s0, s1)][k]
print(f'  per-rectangle mismatches between located zeros and winding numbers: {mism}')
main = [z for z in found if 0.5 < z.real < 1.0]
print('zeros with 1/2 < Re < 1, 0 < Im < 3000 (Re, Im, cycles per octave):')
for z in sorted(main, key=lambda z: -z.real):
    print(f'  {z.real:.6f} + {z.imag:10.4f} i   {z.imag * L2 / (2 * math.pi):8.3f}')
hi = [z for z in found if z.real > 1.0]
print(f'zeros with Re > 1 and 0 < Im < 3000: {len(hi)}; first few: ' +
      ', '.join(f'{z.real:.4f}+{z.imag:.3f}i' for z in hi[:6]) + f'; max Re {max(z.real for z in hi):.4f}')
P = json.load(open(PRED))
pz = [complex(float(a), float(b)) for a, b in P['zeros']]
miss = [z for z in pz if min(abs(z - w) for w in main) > 1e-6]
extra = [w for w in main if w.real > 0.65 and min(abs(z - w) for z in pz) > 1e-6]
print(f'prediction JSON: {len(pz)} zeros; not among the located zeros: {len(miss)}; located zeros with Re > 0.65 '
      f'missing from the JSON: {len(extra)}')
print(f'numbers of zeros with Re > 0.85 / 0.65 / 0.5 (and < 1): {sum(z.real > 0.85 for z in main)} / '
      f'{sum(z.real > 0.65 for z in main)} / {len(main)}')
print(f'min |D\'| over the located zeros: {min(abs(Dp(z)) for z in found):.4f}; '
      f'a common zero of D and D\' would need 2^-s = -1/(1 + ln2/ln1.5) and (3/2)^s/3 = (ln2/ln1.5) 2^-s, '
      f'i.e. Re s = {math.log(1 + L2 / L15) / L2:.4f} and Re s = {math.log(3 * (L2 / L15) / (1 + L2 / L15)) / L15:.4f}: incompatible')
