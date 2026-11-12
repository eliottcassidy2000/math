#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: the mixing rate of H modulo 3 (floating point).

The chain x -> x A_l (mod 3) on the orbit O of x(empty) commutes with the translation t(x) = x + (1,...,1).  The
push-forward operator P on complex measures on O therefore preserves the eigenspaces of t; on the eigenspace where
t acts by a primitive cube root of unity the spectral radius of P/4 is the rate at which the classes of H modulo 3
equidistribute.  Power iteration WITH re-projection on that eigenspace at every step (without it the rounding
errors seed the stationary component, which then dominates: the first version of this check returned 1.000000).
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "..", "minrep17.json")))
r = d["dimension"]
A = [np.array(d["A"][l], dtype=np.int64) for l in range(4)]
x0 = np.array(d["x0"], dtype=np.int64)
m = 3
pw = m ** np.arange(r, dtype=np.int64)
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
t1 = np.searchsorted(keys, key((states + 1) % m))
t2 = t1[t1]
s0 = int(np.searchsorted(keys, key(x0 % m)))
om = np.exp(2j * np.pi / 3)
print("orbit:", N, "states")
for sign in (1, -1):
    w = om ** sign

    def proj(mu):
        # component with mu(t x) = w mu(x)
        return (mu + np.conj(w) * mu[t1] + np.conj(w) ** 2 * mu[t2]) / 3

    mu = np.zeros(N, complex)
    mu[s0] = 1.0
    mu = proj(mu)
    assert np.allclose(mu[t1], w * mu)
    mu /= np.linalg.norm(mu)
    logs = []
    for n in range(1, 601):
        q = np.zeros(N, complex)
        for l in range(4):
            q += np.bincount(nxt[l], weights=mu.real, minlength=N) + 1j * np.bincount(nxt[l], weights=mu.imag, minlength=N)
        q = proj(q)
        nrm = np.linalg.norm(q)
        logs.append(np.log(nrm / 4))
        mu = q / nrm
    logs = np.array(logs)
    print(f"eigenspace t = omega^{sign:+d}: rate from steps 201..400: {np.exp(logs[200:400].mean()):.7f}, "
          f"from steps 401..600: {np.exp(logs[400:600].mean()):.7f}   (lane: 0.7706329)")
    # readout: sum_x mu(x) omega^(-sign * H(x)) is the twisted sum c_n up to normalisation; it must not vanish identically
    c = (mu * om ** (-sign * (states[:, 0] % 3))).sum()
    print(f"   |twisted readout of the normalised iterate| = {abs(c):.3e}")
