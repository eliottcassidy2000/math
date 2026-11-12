#!/usr/bin/env python3
"""Verifier: numerical check of the two inequalities used in the proof that F(W) > 0 for every proper nonempty W
(F(W) = V U - (2^k - 1)(Wd + sigma eps)), all W with k <= 16:
  sigma eps = 0:  F >= sum over runs (2^(k-r) + 2^r - 2) >= 1;
  sigma eps = 1, W proper:  F >= 2^(r_1) + 2^(r_s) + 2^(k - r_1 - r_s) - 3 >= 3  (r_1, r_s = lengths of the runs at 1 and at k);
  single run of length r not containing both ends: F = 2^(k-r) + 2^r - 2 exactly."""
bad = 0
n = 0
for k in range(1, 17):
    full = (1 << k) - 1
    for mask in range(1, full):
        W = [p for p in range(1, k + 1) if (mask >> (p - 1)) & 1]
        Ws = set(W)
        sig = 1 if 1 in Ws else 0
        eps = 1 if k in Ws else 0
        V = mask
        U = sum(1 << (k - p) for p in W)
        Wd = sum(1 for p in W if p + 1 in Ws) + sum(1 << (e - f - 1) for i, e in enumerate(W) for f in W[:i])
        F = V * U - full * (Wd + sig * eps)
        runs = []
        r = 0
        for p in range(1, k + 2):
            if p in Ws:
                r += 1
            elif r:
                runs.append(r)
                r = 0
        n += 1
        if sig * eps == 0:
            lb = sum((1 << (k - r)) + (1 << r) - 2 for r in runs)
            if F < lb or lb < 1:
                bad += 1
            if len(runs) == 1 and F != lb:
                bad += 1
        else:
            r1, rs = runs[0], runs[-1]
            lb = (1 << r1) + (1 << rs) + (1 << (k - r1 - rs)) - 3
            if F < lb or lb < 3:
                bad += 1
print(f"proper nonempty W, k <= 16: {n} sets, violations of the inequalities of the proof: {bad}")
