"""Lovasz theta of the unit-circle Cayley graphs G_q = Cay(F_(q^2), mu_(q+1)) by the abelian-Cayley LP,
giving lower bounds kappa(q) = chi(G_q) >= q^2 / theta(G_q) and upper bounds alpha(G_q) <= theta(G_q).

For a Cayley graph on a finite abelian group with symmetric connection set S,
    theta(G) = max { sum_x f(x) : f(0) = 1, f = 0 on S, f symmetric, fhat >= 0 },
a linear program in the orbit variables (the Fourier transform fhat is real for symmetric f).
Uses the field construction gq() of logic4_20261012/gq_automorphisms.py (same session folder).
Run from the repository root: python3 -B <script>
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
src = (ROOT/'04-computation/experiments/logic4_20261012/gq_automorphisms.py').read_text()
exec(src.split('qs=[q for q')[0])  # defines gq(q)

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def theta(q):
    adj, p, e, mu, elems, idx, add = gq(q); n = len(adj)
    E = np.array(elems)
    S = set(adj[0])
    dots = (E@E.T) % p
    C = np.cos(2*np.pi*dots/p)
    neg = {i: idx[tuple((-np.array(elems[i])) % p)] for i in range(n)}
    reps = sorted({min(i, neg[i]) for i in range(n)}); rid = {r: k for k, r in enumerate(reps)}
    m = len(reps)
    A = np.zeros((n, m)); cvec = np.zeros(m)
    for i in range(n):
        A[:, rid[min(i, neg[i])]] += C[:, i]; cvec[rid[min(i, neg[i])]] += 1
    bounds = [(0, None)]*m; bounds[rid[0]] = (1, 1)
    for s_ in S:
        bounds[rid[min(s_, neg[s_])]] = (0, 0)
    res = linprog(-cvec, A_ub=-A, b_ub=np.zeros(n), bounds=bounds, method='highs')
    check(res.status == 0, 'LP solved')
    return n, len(S), -res.fun


out = {}
for q in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25):
    n, deg, th = theta(q)
    out[q] = dict(n=n, degree=deg, theta=th, chi_lower=int(np.ceil(n/th-1e-9)), alpha_upper=th)
# calibration against the known values of THM-4624: kappa(2) = 4 (K_4), kappa(3) = 3 (K_3 x K_3), kappa(8) = 4, kappa(9) = 3
check(abs(out[2]['theta']-1) < 1e-9 and abs(out[3]['theta']-3) < 1e-9 and abs(out[8]['theta']-16) < 1e-9 and abs(out[9]['theta']-27) < 1e-9, 'theta on the frame cases')
check(out[13]['chi_lower'] == 4 and out[17]['chi_lower'] == 4 and out[23]['chi_lower'] == 5, 'LP lower bounds kappa(13), kappa(17) >= 4 and kappa(23) >= 5')
out['checks'] = CHECKS
print(json.dumps(out, indent=1))
