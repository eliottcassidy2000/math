#!/usr/bin/env python3
"""
h_spectrum_n8_exhaustive_boxeph_S152b.py  (HYP-8220 stage 2)

COMPLETE n=8 h-spectrum by augmentation from the 456 canonical n=7 reps
(every 8-tournament minus a vertex is a 7-tournament => 456 x 128 = 58,368
candidates cover all 6880 classes).  Tests:
  * min strong h(8) = 41 = 2F(8)-1?  (Leonardo/Fibonacci strong-floor law,
    exact 3,5,9,15,25 at n=3..7)
  * do 7, 21, 63 stay missing at n=8?  (the 7*3^k gap tower)
  * full spectrum + missing odds + strong holes.
boxeph-2026-07-20-S152.
"""
N7, N8 = 7, 8

def hpaths(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c or not (S >> v) & 1: continue
            m = adj[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                dp[S | (1 << w)][w] += c
                m &= m - 1
    return sum(dp[(1 << n) - 1][v] for v in range(n))

def is_strong(adj, n):
    full = (1 << n) - 1
    for s in range(n):
        seen = 1 << s
        stack = [s]
        while stack:
            u = stack.pop()
            m = adj[u] & ~seen
            while m:
                w = (m & -m).bit_length() - 1
                seen |= 1 << w
                stack.append(w)
                m &= m - 1
        if seen != full: return False
    return True

pairs7 = [(i, j) for i in range(N7) for j in range(i + 1, N7)]
reps = [int(l) for l in open("05-knowledge/results/n7_class_reps_boxeph_S152.txt")]
print("n=7 reps loaded: %d" % len(reps))
h_all, h_strong = {}, {}
done = 0
for r in reps:
    adj7 = [0] * N7
    for k, (i, j) in enumerate(pairs7):
        if (r >> k) & 1: adj7[i] |= 1 << j
        else: adj7[j] |= 1 << i
    for pat in range(1 << N7):
        adj = [adj7[v] | ((1 << N7) if not (pat >> v) & 1 else 0) for v in range(N7)]
        adj.append(pat)
        h = hpaths(adj, N8)
        h_all[h] = h_all.get(h, 0) + 1
        if is_strong(adj, N8):
            h_strong[h] = h_strong.get(h, 0) + 1
    done += 1
    if done % 50 == 0: print("  ...%d/456 reps done" % done, flush=True)

print("\nCOMPLETE n=8 h-spectrum: %d distinct values" % len(h_all))
sa = sorted(h_all)
print("first 40:", sa[:40])
print("last 10:", sa[-10:])
print("\nh=7: %s  h=21: %s  h=63: %s  (gap tower 7*3^k)" %
      (7 in h_all, 21 in h_all, 63 in h_all))
print("h=189: %s  h=105=3*5*7: %s" % (189 in h_all, 105 in h_all))
ms = min(h_strong)
print("\nMIN strong h(8) = %d  (Leonardo law predicts 41 = 2F(8)-1): %s"
      % (ms, ms == 41))
missing = [v for v in range(1, sa[-1] + 1, 2) if v not in h_all]
print("missing odd values at n=8 (up to max): %s%s"
      % (missing[:30], " ..." if len(missing) > 30 else ""))
strong_holes = [v for v in range(ms, max(h_strong) + 1, 2) if v not in h_strong]
print("strong-spectrum holes above the floor: %s%s"
      % (strong_holes[:30], " ..." if len(strong_holes) > 30 else ""))
print("DONE.")
