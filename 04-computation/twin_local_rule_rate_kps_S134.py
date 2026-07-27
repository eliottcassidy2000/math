#!/usr/bin/env python3
"""HYP-9025 prediction 2: the consecutive-rule rate (gap in K) is
predicted by the same Cramer model as THM-2447 PART3.
Per dyadic window: fit lambda_j on log(N_j(g)/w(g)) ~ -lambda W(g);
model rate = sum_{g in K} w(g) e^{-lambda W(g)} / sum_g w(g) e^{-lambda W(g)}.
Compare with observed rate."""
import numpy as np
from collections import Counter

LIMIT = 100_000_000
GMAX = 400
sieve = np.ones(LIMIT + 3, dtype=bool)
sieve[:2] = False
for p in range(2, int((LIMIT + 2) ** 0.5) + 1):
    if sieve[p]:
        sieve[p * p:: p] = False
mid = (np.where(sieve[:-2] & sieve[2:])[0] + 1).astype(np.int64)
K = (mid[mid >= 6] // 6)
Kset = set(K.tolist())
gaps = np.diff(K)
kv = K[:-1]

PR = [p for p in range(5, 200) if all(p % q for q in range(2, int(p ** .5) + 1))]
def w_of(g):
    w = 1.0
    for p in PR:
        if g % p == 0: w *= (p - 2) / (p - 4)
        elif (9 * g * g - 1) % p == 0: w *= (p - 3) / (p - 4)
    return w
W = np.zeros(GMAX + 1)
wv = np.zeros(GMAX + 1)
acc = 0.0
for g in range(1, GMAX + 1):
    wv[g] = w_of(g)
    acc += wv[g]
    W[g] = acc
inK = np.array([1 if g in Kset else 0 for g in range(GMAX + 1)])

print("window        observed_rate  model_rate  rel.err")
for j in range(17, 24):
    lo, hi = 2 ** j, 2 ** (j + 1)
    sel = (kv >= lo) & (kv < hi)
    gw = gaps[sel]
    obs = float(np.mean([1 if int(g) in Kset else 0 for g in gw]))
    Nj = Counter(int(g) for g in gw if g <= GMAX)
    xs, ys = [], []
    for g in range(1, GMAX + 1):
        if Nj.get(g, 0) >= 80:
            xs.append(W[g]); ys.append(np.log(Nj[g] / wv[g]))
    lam, b = np.polyfit(np.array(xs), np.array(ys), 1)
    mass = wv[1:] * np.exp(lam * W[1:])
    model = float((mass * inK[1:]).sum() / mass.sum())
    print(f"[2^{j},2^{j+1})   {obs:.4f}        {model:.4f}      {abs(model-obs)/obs:.2%}")
