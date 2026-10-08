#!/usr/bin/env python3
"""audit_C item 6: HYP-9240 'Nature of the problem' remark.

For D = 1..DMAX: the rational Terras orbit of y* = (2 - 3^D)/3^D = -1 + 2*3^-D until its denominator clears.
Checks:
  (1) exactly D odd steps clear the denominator (s_0 = time of the D-th odd step); landing N_D = T^(s_0)(y*) >= 0;
  (2) the archimedean bound 1 <= N_D + 1 <= (3/2)^D (in w = y + 1: odd step w -> 3w/2, even step w -> (w + 1)/2),
      and how far the true landings are from it;
  (3) the clearing word (parity bits at times 0..s_0-1) equals the Terras parity vector of the 2-adic integer
      -1 + 2z, z = 3^-D mod 2^(s_0+1) (integers mod 2^(s_0+1) only), and N_D = (2 - 3^D + B_word)/2^(s_0) with
      B_word accumulated from the word alone (B -> 3B + 2^t at an odd step at time t): N_D is determined by D and the
      first s_0 binary digits of 3^-D;
  (4) ranges of s_0/D and max N_D (HYP-9240: N_D <= 880 for D <= 3000).
"""
import sys, math
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
P3 = [1]
for i in range(DMAX + 1):
    P3.append(P3[-1] * 3)

maxN = (0, 0)
worst = (0.0, 0)
s0D = []
bad = 0
for D in range(1, DMAX + 1):
    a, d = 2 - P3[D], D            # y = a / 3^d, gcd(a, 3) = 1 throughout
    word = []
    while d > 0:
        if a & 1:
            word.append(1)
            a = (a + P3[d - 1]) >> 1
            d -= 1
        else:
            word.append(0)
            a >>= 1
    N = a
    s0 = len(word)
    assert sum(word) == D and N >= 0
    # (2)
    if not (1 <= N + 1 and (N + 1) * 2 ** D <= P3[D]):
        bad += 1
        print("archimedean bound fails at D =", D, N)
    if D >= 10:
        r = math.log(N + 1) / (D * math.log(1.5))
        if r > worst[0]:
            worst = (r, D)
    # (3) word from the binary digits of 3^-D only
    M = 1 << (s0 + 1)
    x = (2 * pow(3, -D, M) - 1) % M
    w2 = []
    for t in range(s0):
        b = x & 1
        w2.append(b)
        x = ((3 * x + 1) >> 1) if b else (x >> 1)
    assert w2 == word, D
    B = 0
    for t, b in enumerate(word):
        if b:
            B = 3 * B + (1 << t)
    num = 2 - P3[D] + B
    assert num % (1 << s0) == 0 and num >> s0 == N, D
    s0D.append(s0 / D)
    if N > maxN[0]:
        maxN = (N, D)
print(f"D = 1..{DMAX}: denominator cleared by exactly D odd steps, N_D >= 0: all")
print(f"archimedean bound 1 <= N_D + 1 <= (3/2)^D: failures {bad}")
print(f"max over 10 <= D <= {DMAX} of log(N_D+1)/log((3/2)^D) = {worst[0]:.4f} (at D = {worst[1]})")
print("clearing word = parity vector of -1 + 2*(3^-D mod 2^(s0+1)), and N_D = (2 - 3^D + B_word)/2^s0: all D")
print(f"max N_D = {maxN[0]} at D = {maxN[1]};  s_0/D in [{min(s0D):.3f}, {max(s0D):.3f}] overall, "
      f"[{min(s0D[99:]):.3f}, {max(s0D[99:]):.3f}] for D >= 100")
