#!/usr/bin/env python3
"""Adversarial verifier: refinement of the primitive twisted rate modulo 8 (power iteration, 700 steps, windows)."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
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
t1 = np.searchsorted(keys, key((states + 2) % m)); t2 = t1[t1]; t3 = t2[t1]
s0 = int(np.searchsorted(keys, key(x0 % m)))
w = 1j
proj = lambda mu: (mu + np.conj(w) * mu[t1] + np.conj(w) ** 2 * mu[t2] + np.conj(w) ** 3 * mu[t3]) / 4
char = 1j ** (-(states[:, 0] % 8 - 1) // 2 % 4)          # readout character of H (odd): i^(-(H-1)/2)
mu = np.zeros(N, complex); mu[s0] = 1.0; mu = proj(mu); mu /= np.linalg.norm(mu)
logs = []; cabs = []
scale = 0.0
for n in range(700):
    q = np.zeros(N, complex)
    for l in range(4):
        q += np.bincount(nxt[l], weights=mu.real, minlength=N) + 1j * np.bincount(nxt[l], weights=mu.imag, minlength=N)
    q = proj(q)
    nrm = np.linalg.norm(q)
    logs.append(np.log(nrm / 4)); scale += np.log(nrm / 4)
    mu = q / nrm
    cabs.append(scale + np.log(abs((mu * char).sum()) + 1e-300))
logs = np.array(logs); cabs = np.array(cabs)
for a, b in ((100, 300), (300, 500), (500, 700)):
    print(f"steps {a}..{b}: rate from the l2 norm {np.exp(logs[a:b].mean()):.6f}; from the twisted readout {np.exp((cabs[b - 1] - cabs[a - 1]) / (b - a)):.6f}")
