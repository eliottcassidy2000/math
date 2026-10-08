#!/usr/bin/env python3
"""Robustness of one-step Lamperti transience: an adversary chooses the coupling lag b in {0,...,4} at every step
(b = 0 freezes the walk, but the adversary must move at least every 10 steps), trying to keep the debt walk near the
origin: it picks the lag minimizing E[|x + xi|^2_A] - |x|^2_A ... no: it picks the lag MAXIMIZING E[V(x + xi)] with
V = |x|_A^(-alpha) (pulling the walk back).  Type A configuration on Z_5, A = Q^-1, Q = [[7,-3,1],[-3,7,-3],[1,-3,7]].
Prediction: even the adversarial walk is transient (V(x_t) a supermartingale far out): |x_t| grows like sqrt t and the
number of returns to the origin stays bounded."""
import random, math
Q = [[7, -3, 1], [-3, 7, -3], [1, -3, 7]]
def inv3(M):
    a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
    det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    return [[(e*i - f*h)/det, (c*h - b*i)/det, (b*f - c*e)/det],
            [(f*g - d*i)/det, (a*i - c*g)/det, (c*d - a*f)/det],
            [(d*h - e*g)/det, (b*g - a*h)/det, (a*e - b*d)/det]]
A = inv3(Q)
def nA(x): return math.sqrt(sum(x[i] * A[i][j] * x[j] for i in range(3) for j in range(3)))
v = [(0, 0, 0), (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
def incs(b): return [tuple(v[(j + b) % 5][k] - v[j][k] for k in range(3)) for j in range(5)]
LAGS = {b: incs(b) for b in range(1, 5)}
alpha = 0.15
def EV(x, b):
    out = 0.0
    for xi in LAGS[b]:
        y = [x[k] + xi[k] for k in range(3)]
        r = nA(y); out += (r ** -alpha if r > 0 else 1e9) / 5
    return out
rnd = random.Random(2)
for mode in ('adversarial', 'uniform'):
    returns = []; finals = []
    for _ in range(300):
        x = [0, 0, 0]; ret = 0; idle = 0
        for t in range(4000):
            if mode == 'adversarial':
                b = max(range(1, 5), key=lambda bb: EV(x, bb)) if any(x) else rnd.randrange(1, 5)
            else:
                b = rnd.randrange(1, 5)
            xi = LAGS[b][rnd.randrange(5)]
            x = [x[k] + xi[k] for k in range(3)]
            if not any(x): ret += 1
        returns.append(ret); finals.append(nA(x))
    print(f"{mode:12s}: mean returns to the origin in 4000 steps {sum(returns)/len(returns):.2f} (max {max(returns)}); "
          f"median |x_T|_A {sorted(finals)[len(finals)//2]:.1f}; fraction with no return after step 200 ...", flush=True)
