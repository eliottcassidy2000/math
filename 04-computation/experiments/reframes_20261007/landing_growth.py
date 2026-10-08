#!/usr/bin/env python3
"""How large can the +1-barrier landing integer N_D be?  w = v + 1 evolves by w -> 3w/2 (v odd) or (w+1)/2 (v even),
w_0 = 2/3^D; landing after exactly D odd steps.  Record max N_D per window of D, and the D attaining records."""
import math
def landing(D):
    a, m, s = 2 - 3**D, D, 0
    while m > 0:
        if a & 1: a = (a + 3**(m-1)) // 2; m -= 1
        else: a //= 2
        s += 1
    return a, s
rec = 0; recs = []
win = {}
DMAX = 6000
for D in range(3, DMAX + 1):
    N, s = landing(D)
    if N > rec: rec = N; recs.append((D, N))
    w = D // 500
    win[w] = max(win.get(w, 0), N)
print("record holders (D, N_D):", recs)
print("max N_D per window of 500 D:", [(500*w, win[w]) for w in sorted(win)])
print("records N/D:", [(D, round(N / D, 3)) for D, N in recs[-8:]])
