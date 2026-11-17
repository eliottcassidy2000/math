#!/usr/bin/env python3
"""da_test.py: (Q-2) midpoint values v = h(A)h(D) + B_t(A)P_D(2) (cores of order <= 7 with every root, strong right sources
of order <= 8) against the bases-<= 9 near-join bitmap nj9.bits, by dyadic window, compared with the gap density of odd
numbers; (Q-3) the same-shape obstruction phi(D) = h(D) + P_D(2) vs h(D'') + 2^k P_{D''}(2); (Q-9) one-vertex-core values."""
import numpy as np
L = 8388608
bits = np.frombuffer(open('nj9.bits', 'rb').read(), dtype=np.uint8)
def member(v):  # v: numpy int64 array
    return (bits[v >> 3] >> (v & 7)) & 1
pairs = set(); 
for line in open('pairs7.txt'):
    n, h, B, e = map(int, line.split()); pairs.add((h, B))
pairs = np.array(sorted(pairs), dtype=np.int64); print(f"distinct rooted pairs (h, B) from orders <= 7: {len(pairs)}")
hp_all = []; hp_strong = []
for line in open('hp8.txt'):
    n, s, h, P = map(int, line.split()); hp_all.append((h, P))
    if s: hp_strong.append((n, h, P))
print(f"classes of orders <= 8: {len(hp_all)}; strong: {len(hp_strong)}")
vals = set()
for n, hD, P in hp_strong:
    v = pairs[:, 0] * hD + pairs[:, 1] * P; v = v[v <= L]; vals.update(v.tolist())
vals = np.array(sorted(vals), dtype=np.int64); print(f"distinct midpoint values <= 2^23: {len(vals)}; smallest 12: {vals[:12].tolist()}")
odd = np.arange(69, L, 2, dtype=np.int64); odd_in = member(odd)
first_miss = int(odd[odd_in == 0][0]); print(f"first odd miss of the bases-<= 9 near-join family above 67: {first_miss}")
print("window j: [midpoint values, misses, miss rate] vs [odd numbers gap density]  (ratio observed/expected)")
for j in range(14, 23):
    lo, hi = 1 << j, min(1 << (j+1), L)
    w = vals[(vals >= lo) & (vals < hi)]; wm = len(w) - int(member(w).sum()) if len(w) else 0
    o = odd[(odd >= lo) & (odd < hi)]; gap = 1 - member(o).mean()
    exp = len(w) * gap
    print(f"  j={j}: {len(w)} values, {wm} missed ({(wm/len(w) if len(w) else 0):.4f}) vs gap density {gap:.4f}; expected {exp:.1f}; ratio {(wm/exp if exp else float('nan')):.2f}")
above = vals[vals > first_miss]; print(f"midpoint values above the reach {first_miss}: {len(above)}, of which missed by the bases-<= 9 family: {len(above) - int(member(above).sum())}")
# Q-3 obstruction
hpA = np.array(hp_all, dtype=np.int64)
forms = set()
for h, P in hp_all:
    for k in range(1, 24):
        v = h + (P << k)
        if v > L: break
        forms.add(v)
cnt = 0; hit = []; strong7 = [(n, h, P) for (n, h, P) in hp_strong if n <= 7]
for n, h, P in strong7:
    if h + P in forms: cnt += 1; hit.append((n, h, P, h + P))
print(f"(Q-3) strong D of orders 3-7: {len(strong7)}; phi(D) = h + P of the form h'' + 2^k P'' (D'' any class of order <= 8, k >= 1): {cnt}; cases: {hit[:12]}")
print("      C_3 (3, 12) -> 15 of that form?", 15 in forms, "; strong order 4 (5, 32) -> 37?", 37 in forms)
phi = sorted(set(h + P for n, h, P in hp_strong))
print(f"(Q-9) one-vertex-core values phi(D) over strong D <= 8: {len(phi)} distinct; 241 in: {241 in phi}; 493: {493 in phi}; 997: {997 in phi}; 2005: {2005 in phi}; smallest 20: {phi[:20]}")
