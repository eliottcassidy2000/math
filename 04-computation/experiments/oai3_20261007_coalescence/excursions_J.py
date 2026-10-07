#!/usr/bin/env python3
"""Number J of excursions (returns of k to 0) before absorption, start (0,1) (y vs y+1); heuristic constant
c1 ~ E[J] * P(excursion > t) t^(1/2) with time ~ 2 x flips: P(X > t) ~ sqrt(4/(pi t))."""
import random, math
from coal_check import step
rng = random.Random(5); paths = 4000; cap = 200000; Js = []; unabs = 0; flips = 0; steps = 0
for _ in range(paths):
    k, N = 0, 1; J = 0; t = 0
    while t < cap:
        k2, N2 = step(k, N, rng.getrandbits(1)); t += 1
        if N & 1: flips += 1
        steps += 1
        if k != 0 and k2 == 0: J += 1
        k, N = k2, N2
        if k == 0 and N == 0: break
    else:
        unabs += 1
    Js.append(J)
EJ = sum(Js) / len(Js)
print(f"start (0,1): E[J] = {EJ:.3f} (unabsorbed by {cap}: {unabs}/{paths}); P(J=1) = {Js.count(1)/paths:.3f}; flip share of steps = {flips/steps:.4f}")
print(f"heuristic c1 = E[J] * sqrt(4/pi) = {EJ*math.sqrt(4/math.pi):.2f}  (measured sqrt(T) q(T) ~ 11.0)")
