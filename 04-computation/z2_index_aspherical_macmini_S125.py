#!/usr/bin/env python3
"""
THM-1385: ind(aspherical free Z/2-space) = 1, and the explicit sharpness witness.
================================================================================
THEOREM. X aspherical, Z/2 acting freely => ind(X) = 1.
  ind>=1: the cover is nontrivial. ind<=1: a Z/2-map S^2 -> X would descend to
  RP^2 -> M = X/Z2 with pi_1(RP^2)=Z/2 -> Gamma -> Z/2 onto, giving a nontrivial
  gamma with gamma^2=1; but Gamma acts freely on a contractible universal cover,
  so Gamma is TORSION-FREE. Contradiction.
CONSEQUENCE. An odd map X -> R^m must vanish only for m=1; for m>=2 zero-free odd
maps exist -- so on T^k, Borsuk-Ulam yields ONE equation, never k.
SHARPNESS (checked here). For tau(x)=x+v, v=w/2, w not in 2Z^k: pick i with w_i odd,
a=e_i; f(x)=exp(2*pi*i*<a,x>) in R^2 is odd and never zero.
"""
import numpy as np

def witness(k, w, i):
    v = w / 2.0
    a = np.zeros(k); a[i] = 1.0
    f = lambda p: np.array([np.cos(2*np.pi*(a@p)), np.sin(2*np.pi*(a@p))])
    return v, a, f

rng = np.random.default_rng(0)
print("sharpness: odd zero-free map T^k -> R^2 for each k")
allok = True
for k in (1, 2, 3, 5, 8, 13):
    w = (rng.integers(0, 4, size=k) * 2).astype(float)
    i = int(rng.integers(0, k)); w[i] = 2*int(rng.integers(0, 3)) + 1   # one odd coord
    v, a, f = witness(k, w, i)
    ok = True
    for _ in range(200):
        x = rng.random(k)
        ok &= np.allclose(f(x + v), -f(x), atol=1e-12)      # ODD
        ok &= np.linalg.norm(f(x)) > 0.5                     # NO ZERO
    ok &= not np.allclose(v % 1.0, 0)                        # involution is free
    allok &= ok
    print(f"  k={k:2d}: w_i odd at i={i}; odd & zero-free & free: {ok}")
print(f"ALL k verified: {allok}")
print()
print("=> ind(T^k) < 2 for every k; with ind >= 1 this gives ind = 1 EXACTLY.")
print("   Contrast ind(S^k) = k: S^k is NOT aspherical for k >= 2, which is the")
print("   entire reason BU works there and collapses on the torus.")
