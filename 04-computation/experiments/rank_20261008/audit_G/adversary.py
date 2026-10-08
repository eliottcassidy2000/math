#!/usr/bin/env python3
"""Audit G, THM-4609 statement 3 (adversarial lags), stronger adversaries than adversarial_lag.py.
Debt walk of the Z_5 AP configuration (non-unit positions 2,3,4): at each step the adversary picks a lag b in 1..4
(lag 0 = freeze is not allowed: freezing cannot help return), the step is uniform over the 5 cyclic roots.
Adversaries: uniform; RADIAL-E (maximise xhat^T C_b xhat / tr C_b, Euclidean); RADIAL-Q (same in Q_AP^-1 metric);
GREEDY-V (maximise E[V(x + xi)], V = |x|_{Q^-1}^-alpha, alpha = 0.03 < alpha_max(Q_AP) = 0.0388);
NEAREST (maximise P(next step hits 0), else radial).  Control: the rank-2 configuration (non-units at 3, 4) with
RADIAL-E, which should be recurrent (many returns).
Output: mean number of returns to 0 per walk in [0, T), and fraction of walks with a return in (T/2, T]."""
import random, math
rnd = random.Random(2026)
Q = [[6, -2, 0], [-2, 6, -1], [0, -1, 5]]       # basis (end, middle, end) = positions (2, 3, 4)
def inv3(M):
    a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
    det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
    return [[(e*i - f*h)/det, (c*h - b*i)/det, (b*f - c*e)/det],
            [(f*g - d*i)/det, (a*i - c*g)/det, (c*d - a*f)/det],
            [(d*h - e*g)/det, (b*g - a*h)/det, (a*e - b*d)/det]]
A = inv3(Q)
def config(nonunits, d=5):
    rho = len(nonunits); v = []
    for j in range(d):
        w = [0] * rho
        if j in nonunits: w[nonunits.index(j)] = 1
        v.append(tuple(w))
    return v
def roots(v, b, d=5):
    return [tuple(x - y for x, y in zip(v[(j + b) % d], v[j])) for j in range(d)]
def quad(x, M): return sum(x[i] * M[i][j] * x[j] for i in range(len(x)) for j in range(len(x)))
def run(nonunits, mode, T=10000, walks=100, alpha=0.03):
    v = config(nonunits); rho = len(nonunits)
    R = {b: roots(v, b) for b in range(1, 5)}
    I = [[float(i == j) for j in range(rho)] for i in range(rho)]
    Cov = {b: [[sum(w[i] * w[j] for w in R[b]) / 5 for j in range(rho)] for i in range(rho)] for b in R}
    rets = []; late = 0
    for _ in range(walks):
        x = [0] * rho; nret = 0; lr = False
        for t in range(T):
            if mode == 'uniform' or not any(x):
                b = rnd.randrange(1, 5)
            elif mode in ('radialE', 'radialQ'):
                if mode == 'radialE' or rho != 3:
                    y = x; G = I
                else:
                    y = [sum(A[i][j] * x[j] for j in range(3)) for i in range(3)]   # metric A: <x, C x>_A-ish
                    G = A
                # radial share of variance in the metric G: (x^T G C G x)/(x^T G x) / tr(C G)
                def share(b):
                    C = Cov[b]
                    Gx = [sum(G[i][j] * x[j] for j in range(rho)) for i in range(rho)]
                    num = sum(Gx[i] * C[i][j] * Gx[j] for i in range(rho) for j in range(rho))
                    den = sum(x[i] * Gx[i] for i in range(rho))
                    tr = sum(C[i][j] * G[j][i] for i in range(rho) for j in range(rho))
                    return num / den / tr
                b = max(range(1, 5), key=share)
            elif mode == 'greedyV':
                def EV(b):
                    s = 0.0
                    for w in R[b]:
                        y = [x[k] + w[k] for k in range(rho)]
                        q = quad(y, A) if rho == 3 else sum(t * t for t in y)
                        s += (q ** (-alpha / 2) if q > 0 else 1e9)
                    return s
                b = max(range(1, 5), key=EV)
            elif mode == 'nearest':
                hits = {b: sum(1 for w in R[b] if all(x[k] + w[k] == 0 for k in range(rho))) for b in range(1, 5)}
                mh = max(hits.values())
                if mh > 0: b = max(hits, key=hits.get)
                else:
                    def share(b):
                        C = Cov[b]; num = sum(x[i] * C[i][j] * x[j] for i in range(rho) for j in range(rho))
                        return num / sum(t * t for t in x) / sum(C[i][i] for i in range(rho))
                    b = max(range(1, 5), key=share)
            w = R[b][rnd.randrange(5)]
            x = [x[k] + w[k] for k in range(rho)]
            if not any(x):
                nret += 1
                if t >= T // 2: lr = True
        rets.append(nret); late += lr
    return sum(rets) / len(rets), late / walks
for mode in ('uniform', 'radialE', 'radialQ', 'greedyV', 'nearest'):
    m, l = run([2, 3, 4], mode)
    print(f"rank 3 (AP, Z_5 example), {mode:8s}: mean returns in 10000 steps {m:.2f}; P(return in second half) {l:.3f}", flush=True)
m, l = run([3, 4], 'radialE'); print(f"rank 2 control (non-units 3,4), radialE: mean returns {m:.2f}; P(return in second half) {l:.3f}")
m, l = run([3, 4], 'uniform'); print(f"rank 2 control (non-units 3,4), uniform: mean returns {m:.2f}; P(return in second half) {l:.3f}")
