"""Auditor: HYP-9370 mechanism section.  (1) truncated Gaussian exp((n-1) z^2 / n) to degree floor((n-1)/2) in w:
smallest-root |z| and arg w; (2) sign structure of the exact P_T: negative coefficients and positive real roots
(w > 0, i.e. real zeros z) on the complete n = 8, 10 censuses and the samples.  Usage: python szego_check.py"""
import os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for n in (8, 10, 12, 14, 16, 100):
    m = (n - 1) // 2; c = (n - 1) / n
    co = [c ** k / math.factorial(k) for k in range(m + 1)]
    r = np.roots(co[::-1]); k = np.argmin(abs(r))
    print(f'n={n}: truncated exp(c z^2), degree {m} in w: smallest |z| = {math.sqrt(abs(r[k])):.4f}, arg w = {abs(np.angle(r[k])):.3f}; '
          f'sqrt(0.2785 m / c) = {math.sqrt(0.27846 * m / c):.4f}')
for tag, n in (('ar8_classes', 8), ('ar10_classes', 10), ('ar12_samples', 12), ('ar14_samples', 14), ('ar16_samples', 16),
               ('ar10_classes_completed', 11)):
    fn = os.path.join(HERE, tag + '.rows')
    if not os.path.exists(fn):
        if tag == 'ar8_classes':
            continue
        print(tag, 'missing'); continue
    negcoef = 0; posroot = 0; tot = 0; allpos_low = 0
    for l in open(fn):
        t = l.split()
        N = list(map(int, t[5:]))
        co = [N[2 * k] / n ** (2 * k) for k in range((n - 1) // 2 + 1)]
        while len(co) > 1 and co[-1] == 0: co.pop()
        tot += 1
        if min(co) < 0: negcoef += 1
        r = np.roots(co[::-1])
        if any(abs(x.imag) < 1e-9 * max(1, abs(x)) and x.real > 0 for x in r): posroot += 1
    print(f'{tag}: {tot} tournaments; P_T has a negative coefficient in {negcoef}; a positive real root w > 0 in {posroot}')
