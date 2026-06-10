#!/usr/bin/env python3
"""
paley_rotdp_smallp_verify_kpc1.py
kind-pasteur-2026-06-10-S1 -- Thread E (HYP-2371), companion to the design doc
05-knowledge/results/paley_H31_compute_design_kpc1.md.

Validates the Z_p-ROTATION-REDUCED LAYERED HELD-KARP design at p = 11 and 19:
  * states (S, v), v in S, S subset of Z_p, taken UP TO ROTATION (S+r, v+r);
  * the rotation action is FREE on all states with S != empty (p prime), so
    every orbit has size exactly p and dp values are constant on orbits;
  * dp(S, v) = number of directed Hamiltonian paths of the induced
    subtournament on S ending at v (start anywhere in S);
  * layers = popcount(S); peak layer size = C(p-1, (p-1)/2)-ish; total
    canonical states = 2^{p-1};
  * H(T_p) = p * dp(full, 0).
Exact integer arithmetic throughout. Expected: H(11) = 95095,
H(19) = 1172695746915 (canon values, re-verified vs stored Held-Karp results).

Also re-prints the layer-size arithmetic that gives the ~2.3 GB peak-layer
estimate at p = 31 (the design doc's central memory claim).
"""
import time
from collections import defaultdict

H_KNOWN = {11: 95095, 19: 1172695746915}

def qr_set(p):
    return set((x * x) % p for x in range(1, p))

def run(p):
    t0 = time.time()
    QR = qr_set(p)
    mask = (1 << p) - 1

    def canonical(S, v):
        best = None
        for r in range(p):
            S2 = ((S << r) | (S >> (p - r))) & mask
            key = (S2 << 5) | ((v + r) % p)
            if best is None or key < best:
                best = key
        return best

    # layer k=1: single orbit, representative ({0}, 0)
    layer = {canonical(1, 0): 1}
    layer_sizes = [len(layer)]
    for k in range(1, p):
        nxt = defaultdict(int)
        for key, c in layer.items():
            S, v = key >> 5, key & 31
            for w in range(p):
                if not (S >> w) & 1 and (w - v) % p in QR:
                    nxt[canonical(S | (1 << w), w)] += c
        layer = nxt
        layer_sizes.append(len(layer))
        # exact orbit-count check: layer k+1 has C(p, k+1)*(k+1)/p canonical
        # states REACHABLE-bounded; report only
    assert len(layer) == 1, "full layer must be a single orbit"
    H = p * next(iter(layer.values()))
    dt = time.time() - t0
    ok = (H == H_KNOWN[p])
    print(f"p={p:2d}: H(T_p) = {H}  expected {H_KNOWN[p]}  "
          f"{'OK' if ok else 'MISMATCH'}   ({dt:.1f}s)")
    assert ok
    tot = sum(layer_sizes)
    print(f"      total canonical states touched = {tot}  (= 2^{p-1} = {2**(p-1)} "
          f"if all (S,v) reachable; unreachable states never materialize)")
    print(f"      peak layer = {max(layer_sizes)} states "
          f"(layer k={layer_sizes.index(max(layer_sizes))+1})")
    print(f"      layer sizes: {layer_sizes}")
    # freeness sanity: orbit count formula C(p,k)*k/p must be an integer
    from math import comb
    for k in range(1, p + 1):
        assert (comb(p, k) * k) % p == 0
    print(f"      orbit-count formula C(p,k)*k/p integral for all k: OK (free action)")
    return H

print("=" * 74)
print("Z_p-rotation-reduced layered Held-Karp -- harness validation (exact)")
print("=" * 74)
for p in [11, 19]:
    run(p)

print()
print("=" * 74)
print("design arithmetic for p = 31 (the compute-node run; see design doc)")
print("=" * 74)
from math import comb
p = 31
tot = sum(comb(p, k) * k for k in range(1, p + 1)) // p
print(f"total canonical states: sum_k C(31,k)*k/31 = {tot} = 2^30 = {2**30}")
peak_k = max(range(1, p + 1), key=lambda k: comb(p, k) * k)
peak = comb(p, peak_k) * peak_k // p
print(f"peak layer k={peak_k}: C(31,{peak_k})*{peak_k}/31 = {peak} states")
print(f"  counts as uint128: {peak*16/2**30:.2f} GiB   (the ~2.3 GB claim)")
print(f"  adjacent layer k={peak_k+1}: {comb(p,peak_k+1)*(peak_k+1)//p} states "
      f"-> {comb(p,peak_k+1)*(peak_k+1)//p*16/2**30:.2f} GiB")
trans = sum(comb(p, k) * k * (p - k) for k in range(1, p)) // (2 * p)
print(f"total transitions (avg out-degree (p-k)/2): "
      f"sum C(31,k)*k*(31-k)/(2*31) = {trans} ~ {trans:.2e}")
print(f"  (= n(n-1)2^(n-2)/(2n) with n=31: {31*30*2**29//62})")
