#!/usr/bin/env python3
"""Where does the resonant Fourier coefficient mu_hat_h(2^s) come from?  Joint law of (Y_h mod 3^h, total cost A)
by the recursion with a cost coordinate (A <= AMAX), then the contribution of each cost A to mu_hat_h(t) at the
resonant t = 2^s and at a generic unit t.  Also the coefficient by the LAST valuation a_h.
Run: python 04-computation/experiments/collatz_five_mirrors_costsplit_20260929.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

if __name__ == "__main__":
    h = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    AMAX = 3 * h + 20
    mod = 3 ** h
    # joint[y, A]
    joint = np.zeros((1, AMAX + 1))
    joint[0, 0] = 1.0
    for lev in range(1, h + 1):
        m = 3 ** lev
        new = np.zeros((m, AMAX + 1))
        y = np.arange(joint.shape[0])
        base = (3 * y + 1) % m
        for a in range(1, AMAX + 1):
            inv = pow(2, -a, m)
            z = (inv * base) % m
            w = 2.0 ** (-a)
            # new[z, A + a] += w * joint[y, A]
            np.add.at(new[:, a:], (z, slice(None)), 0)  # no-op to keep shapes explicit
            shifted = w * joint[:, : AMAX + 1 - a]
            np.add.at(new, (z[:, None], np.arange(a, AMAX + 1)[None, :]), shifted)
        joint = new
    mu = joint.sum(axis=1)
    print(f"h={h}: total mass {mu.sum():.6f} (cost truncation at A <= {AMAX}); P(A) mean = {(joint.sum(axis=0) * np.arange(AMAX + 1)).sum():.3f}")
    mh = np.fft.fft(mu)
    amp = np.abs(mh)
    units = [t for t in range(1, mod) if t % 3]
    s_star = int(np.argmax([amp[pow(2, s, mod)] for s in range(0, 60)]))
    t_star = pow(2, s_star, mod)
    print(f"   resonant t = 2^{s_star} mod 3^h = {t_star}: |mu_hat| = {amp[t_star]:.5f}; generic t = 7: {amp[7]:.5f}; max over units {max(amp[t] for t in units):.5f}")
    phase = np.exp(-2j * math.pi * t_star * np.arange(mod) / mod)
    contrib = phase @ joint  # per cost A
    tot = contrib.sum()
    print(f"   sum over A of contributions = {abs(tot):.5f} (check)")
    order = np.argsort(-np.abs(contrib))[:12]
    print("   largest |contribution| by total cost A: " + ", ".join(f"A={A}: {abs(contrib[A]):.4f} (mass {joint[:, A].sum():.4f}, phase {np.angle(contrib[A]):+.2f})" for A in order))
    # coherence: |sum| vs sum of |.| over A
    print(f"   sum_A |contribution| = {np.abs(contrib).sum():.4f}; coherent fraction |sum|/sum|.| = {abs(tot) / np.abs(contrib).sum():.3f}")
    # by last valuation a_h: Y_h = 2^-a (3 Y_(h-1) + 1): contribution of each a
    pm = 3 ** (h - 1)
    prev = joint  # not needed; recompute mu_(h-1) by reduction
    mu_prev = mu.reshape(3, pm).sum(axis=0)  # consistency: reduction mod 3^(h-1)
    base = (3 * np.arange(pm) + 1) % mod
    parts = []
    for a in range(1, 40):
        inv = pow(2, -a, mod)
        z = (inv * base) % mod
        parts.append((a, (2.0 ** (-a)) * np.exp(-2j * math.pi * t_star * z / mod) @ mu_prev))
    print("   by last valuation a: " + ", ".join(f"a={a}: {abs(c):.4f}@{np.angle(c):+.2f}" for a, c in parts[:12]))
    print(f"   sum over a = {abs(sum(c for a, c in parts)):.5f} (equals |mu_hat| up to the a > 39 tail)")
    print("DONE")
