#!/usr/bin/env python3
"""Does R+ (natural order = best rotation order, THM-4666 1c) separate better by L(1, chi_-p) together with small-
conductor twists than by L(1, chi_-p) alone? (S40, HYP-9287 mechanism)

Features for each prime in the Farey-sweep census with L >= 2.4 and p > 479: L1 = L(1, chi_-p), the ratios
L(1, chi_-8p)/L1 and L(1, chi_-5p)/L1 (imaginary quadratic twists, class numbers from reduced forms), and log p. For each
feature set the script reports the best leave-one-out accuracy of a threshold rule (1 feature) or of a logistic
regression (several features).
Usage: python paley_rplus_twist_separation_20261017.py
Session opus-2026-10-17-S40.
"""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from paley_class_number_ranking_20261017 import class_number


def L_imag(D):
    return math.pi * class_number(D) / math.sqrt(-D)


def logistic_loo(X, y, iters=3000, lr=0.5, lam=1e-3):
    n = len(y); correct = 0
    Xs = (X - X.mean(0)) / X.std(0)
    Xs = np.hstack([np.ones((n, 1)), Xs])
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        w = np.zeros(Xs.shape[1])
        A, b = Xs[m], y[m]
        for _ in range(iters):
            z = 1 / (1 + np.exp(-A @ w))
            w -= lr * (A.T @ (z - b) / len(b) + lam * w)
        pred = (1 / (1 + np.exp(-Xs[i] @ w))) > 0.5
        correct += int(pred == y[i])
    return correct / n


def main():
    d = json.load(open(os.path.join(HERE, '..', '..', '05-knowledge', 'results', 'paley_rotation_farey_sweep_480_20000_20261017.json')))['range']
    rows = []
    for r in d:
        p = r['p']
        L1 = r['L']
        L8 = L_imag(-8 * p)
        L5 = L_imag(-5 * p) if p % 5 else float('nan')
        rows.append((p, L1, L8 / L1, L5 / L1, math.log(p), 0 if r['beaten'] else 1))
    R = np.array([r for r in rows if not math.isnan(r[3])])
    y = R[:, 5].astype(int)
    print('primes %d (R+ members %d)' % (len(y), y.sum()))
    # one-feature threshold on L1 (best split, in-sample upper bound)
    best = max(((np.mean((R[:, 1] >= t) == (y == 1)), t) for t in np.unique(R[:, 1])))
    print('threshold on L1 alone: best accuracy %.3f at L1 >= %.3f' % best)
    for name, cols in [('L1, log p', [1, 4]), ('L1, L8/L1', [1, 2]), ('L1, L8/L1, log p', [1, 2, 4]),
                       ('L1, L8/L1, L5/L1, log p', [1, 2, 3, 4])]:
        acc = logistic_loo(R[:, cols], y)
        print('logistic regression on %-26s leave-one-out accuracy %.3f' % (name, acc), flush=True)


if __name__ == '__main__':
    main()
