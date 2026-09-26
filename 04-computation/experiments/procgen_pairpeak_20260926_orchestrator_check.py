#!/usr/bin/env python3
"""Orchestrator spot-check of lane `pairpeak`, Lemma 1 (coupling of a single
pairing flip), written from the note's statement; the lane's scripts were not read.
For odd v put x_1 = (v-1)/2 (the flipped image) and y_1 = T(v) = (3v+1)/2.
Claim: if the Collatz word of y_1 begins 1^r 0 0, the orbits of x_1 and y_1 merge
at time r+3 (x_(r+3) = y_(r+3)); and over odd v this event has frequency 1/2.
"""
def T(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2

def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)

hits = total = merged_other = 0
for v in range(3, 400001, 2):
    x, y = (v - 1) // 2, (3 * v + 1) // 2
    # word of y_1
    w, z = [], y
    for _ in range(64):
        w.append(z % 2); z = T(z)
    r = 0
    while r < len(w) and w[r] == 1:
        r += 1
    total += 1
    if r + 2 <= len(w) and w[r] == 0 and w[r + 1] == 0:
        hits += 1
        xs, ys = x, y
        for t in range(1, r + 3):          # advance from time 1 to time r+3
            xs, ys = T(xs), T(ys)
        assert xs == ys, (v, r, xs, ys)
check(True, f"for all odd 3 <= v <= 400001 whose T(v)-word begins 1^r00, the flipped and unflipped orbits coincide at time r+3 ({hits} cases)")
freq = hits / total
check(abs(freq - 0.5) < 0.01, f"the event '1^r 00' has frequency {freq:.4f} ~ 1/2 over odd v (re-merge probability 1/2 per flip)")
