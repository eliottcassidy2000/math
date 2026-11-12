#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: (a) the leading term H = 2 3^n S_n E_n + O(n 2^n) on all words, n <= 9 (own H);
(b) distribution function of H / (2 3^17) at x = 2^-j and the median, from c_17 (the lane's array, shown identical
to the verifier's independent census); (c) twisted rate of the primitive characters modulo 8 (power iteration with
re-projection, floating point).

usage: python v9_limit_law.py LANE_STATE_DIR
"""
import json
import os
import sys

import numpy as np

from v1_core import FIN, all_alpha

HERE = os.path.dirname(os.path.abspath(__file__))
ldir = sys.argv[1]

# ---- (a)
exp = {2: 0.6250, 3: 0.5972, 4: 0.5781, 5: 0.6062, 6: 0.5940, 7: 0.5864, 8: 0.5654, 9: 0.5561}
oka = True
for n in range(2, 10):
    Hn = (all_alpha(n) @ FIN).astype(np.float64)
    idx = np.arange(4 ** n)
    a = np.stack([(idx >> (2 * (n - 1 - k))) & 1 for k in range(n)], axis=1)
    b = np.stack([(idx >> (2 * (n - 1 - k) + 1)) & 1 for k in range(n)], axis=1)
    S = np.zeros(4 ** n)
    E = np.zeros(4 ** n)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i != j:
                S += a[:, i - 1] * b[:, j - 1] * 2.0 ** (abs(i - j) - 1) / 3.0 ** max(i, j)
                E += (1 - a[:, i - 1]) * (1 - b[:, j - 1]) * 2.0 ** (abs(i - j) - 1) / 3.0 ** (n + 1 - min(i, j))
    dev = np.abs(Hn - 2 * 3.0 ** n * S * E)
    r = dev.max() / (n * 2 ** n)
    oka &= abs(r - exp[n]) < 6e-5
    print(f"(a) n = {n}: max |H - 2 3^n S_n E_n| / (n 2^n) = {r:.4f} (lane {exp[n]:.4f})")
print("(a) leading-term table reproduced:", oka)

# ---- (b)
st = np.load(os.path.join(ldir, "census_state_b26_L17.npz"))
c17 = st["last"].astype(np.int64)
cs = np.cumsum(c17)
tot = 4 ** 17
lane = {2: 0.97941, 3: 0.85323, 4: 0.63500, 5: 0.41843, 6: 0.25213, 7: 0.14174, 8: 0.07569, 9: 0.03871, 10: 0.01908}
okb = True
vals = {}
for j in range(2, 11):
    thr = 2 * 3 ** 17 / 2 ** j
    k = int(np.floor((thr - 1) / 2))            # largest index with 2k + 1 <= thr
    if 2 * k + 1 >= thr:
        k -= 1
    p = cs[k] / tot
    vals[j] = p
    okb &= abs(p - lane[j]) < 6e-6
print("(b) P(H < 2 3^17 / 2^j), j = 2..10: " + " ".join(f"{vals[j]:.5f}" for j in range(2, 11)), "| reproduced:", okb)
print("(b) local exponents: " + " ".join(f"{np.log2(vals[j] / vals[j + 1]):.3f}" for j in range(2, 10)))
med = 2 * int(np.searchsorted(cs, tot / 2)) + 1
print(f"(b) median of H over the words of length 17: {med} = {med / (2 * 3 ** 17):.5f} * 2 * 3^17 (lane: 0.04114)")

# ---- (c)
d = json.load(open(os.path.join(HERE, "..", "minrep17.json")))
r17 = d["dimension"]
A = [np.array(d["A"][l], dtype=np.int64) for l in range(4)]
x0 = np.array(d["x0"], dtype=np.int64)
m = 8
pw = m ** np.arange(r17, dtype=np.int64)
key = lambda V: (V % m) @ pw
states = (x0 % m)[None, :]
keys = np.array([int(key(x0))], dtype=np.int64)
frontier = states
while len(frontier):
    new = np.concatenate([(frontier @ A[l]) % m for l in range(4)], axis=0)
    nk, ui = np.unique(key(new), return_index=True)
    new = new[ui]
    mask = ~np.isin(nk, keys)
    frontier = new[mask]
    keys = np.concatenate([keys, nk[mask]])
    states = np.concatenate([states, frontier], axis=0)
order = np.argsort(keys)
keys, states = keys[order], states[order]
N = len(keys)
nxt = [np.searchsorted(keys, key((states @ A[l]) % m)) for l in range(4)]
t1 = np.searchsorted(keys, key((states + 2) % m))
t2 = t1[t1]
t3 = t2[t1]
s0 = int(np.searchsorted(keys, key(x0 % m)))
for jchar, expect in ((1, 0.6175), (2, 0.7327856)):
    w = 1j ** jchar

    def proj(mu):
        return (mu + np.conj(w) * mu[t1] + np.conj(w) ** 2 * mu[t2] + np.conj(w) ** 3 * mu[t3]) / 4
    mu = np.zeros(N, complex)
    mu[s0] = 1.0
    mu = proj(mu)
    mu /= np.linalg.norm(mu)
    logs = []
    for n in range(300):
        q = np.zeros(N, complex)
        for l in range(4):
            q += np.bincount(nxt[l], weights=mu.real, minlength=N) + 1j * np.bincount(nxt[l], weights=mu.imag, minlength=N)
        q = proj(q)
        nrm = np.linalg.norm(q)
        logs.append(np.log(nrm / 4))
        mu = q / nrm
    rate = float(np.exp(np.mean(logs[100:300])))
    print(f"(c) modulus 8, orbit {N} states, component where x -> x + 2(1,..,1) acts by i^{jchar}: rate {rate:.6f} (lane: {expect})")
