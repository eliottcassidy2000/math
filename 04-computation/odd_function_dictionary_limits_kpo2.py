#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
odd_function_dictionary_limits_kpo2.py
kind-pasteur-2026-06-10-S2 / Thread A / HYP-2378 follow-up.

QUESTION left open by odd_function_dictionary_kpo2.py part (c):
the rotation (carousel) circulant family S = {1..(n-1)/2} has cluster
generators a_{2k} = [x^{2k-1}] tanh(x)  (PROVED there), so the NAIVE
linked-cluster formula predicts R_rot -> exp(tanh 1) ~ 2.14169.  But exact
R_rot(n) already exceeds that at n=9 and climbs through 2.487 at n=15.
So the exponential formula -- proved for quasirandom/DRT families in
THM-438 -- must fail beyond them.  WHERE does R_rot head?

This script:
 1. extends the EXACT R_rot(n) sequence to n = 17, 19 via a fixed-start
    Held-Karp (paths starting at vertex 0; H = n * count, valid because Z_n
    acts freely on directed Ham paths -- LEM-003), exact integers;
 2. cross-validates the fixed-start trick against full Held-Karp at n=11,13;
 3. estimates R at n = 21, 31, 51 (and 101 if time allows) by UNBIASED
    sequential importance sampling (SIS), validated against the exact values;
    Paley p = 31, 43 as a positive control (should keep heading to e).
Estimates are CLEARLY LABELED estimates (mean +- 2*stderr); all exact counts
are pure-int.  Fresh code (MISTAKE-001 guard).
"""

import sys, math, random, time

PASSED, FAILED = [], []
def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    (PASSED if ok else FAILED).append(name)
    print("  [%s] %s%s" % (tag, name, ((" -- " + detail) if detail else "")))

def sigma_rotation(n):
    h = (n - 1) // 2
    return [0] + [1] * h + [-1] * h

def sigma_paley(p):
    assert p % 4 == 3
    sgl = [0] * p
    for d in range(1, p):
        sgl[d] = 1 if pow(d, (p - 1) // 2, p) == 1 else -1
    return sgl

def adj_from_sigma(n, sgl):
    return [[1 if (i != j and sgl[(j - i) % n] == 1) else 0 for j in range(n)]
            for i in range(n)]

def ham_path_count_full(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                for w in range(n):
                    if not (mask >> w) & 1 and av[w]:
                        dp[mask | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])

def ham_paths_start0_times_n(n, sgl):
    """H = n * #(directed Ham paths starting at vertex 0).
    Valid for circulant tournaments: Z_n <= Aut acts freely on Ham paths
    (LEM-003), and the n rotations of a path realize each start vertex once.
    DP state: (mask over vertices 1..n-1, last vertex).  Exact ints."""
    m = n - 1
    # adjacency bitmasks among vertices 1..n-1 (index v-1)
    succ = []
    for v in range(1, n):
        bits = 0
        for w in range(1, n):
            if w != v and sgl[(w - v) % n] == 1:
                bits |= 1 << (w - 1)
        succ.append(bits)
    size = 1 << m
    dp = [0] * (size * m)
    for v in range(1, n):
        if sgl[v % n] == 1:           # arc 0 -> v
            dp[(1 << (v - 1)) * m + (v - 1)] = 1
    for mask in range(size):
        base = mask * m
        for v in range(m):
            c = dp[base + v]
            if c:
                avail = succ[v] & ~mask
                while avail:
                    b = avail & (-avail)
                    w = b.bit_length() - 1
                    dp[(mask | b) * m + w] += c
                    avail ^= b
    full = (size - 1) * m
    return n * sum(dp[full + v] for v in range(m))

def sis_R(n, sgl, samples, seed):
    """Unbiased SIS estimator of H, reported as R = H * 2^(n-1)/n!.
    Start vertex 0 (wlog: vertex-transitive => H = n*H_0), weight = product of
    branch counts.  Returns (R_mean, R_2stderr)."""
    rng = random.Random(seed)
    out = [[w for w in range(n) if w != v and sgl[(w - v) % n] == 1]
           for v in range(n)]
    tot = 0.0
    tot2 = 0.0
    for _ in range(samples):
        used = 1
        cur = 0
        w = float(n)               # start-vertex factor (H = n * H_0)
        for _ in range(n - 1):
            cands = [u for u in out[cur] if not (used >> u) & 1]
            k = len(cands)
            if k == 0:
                w = 0.0
                break
            w *= k
            cur = cands[int(rng.random() * k)]
            used |= 1 << cur
        tot += w
        tot2 += w * w
    mean = tot / samples
    var = max(tot2 / samples - mean * mean, 0.0)
    se = math.sqrt(var / samples)
    scale = 2.0 ** (n - 1) / math.factorial(n)
    return mean * scale, 2 * se * scale

print("=" * 78)
print("WHERE DOES R_rot HEAD?  exact extension + SIS estimates")
print("=" * 78)
tanh1 = math.tanh(1.0)
print("reference: e = %.5f ; naive linked-cluster exp(tanh 1) = %.5f"
      % (math.e, math.exp(tanh1)))
print()

print("1. fixed-start trick validation (exact, pure-int)")
for n, target in [(11, 93027), (13, 3711175)]:
    h_full = ham_path_count_full(adj_from_sigma(n, sigma_rotation(n)))
    h_fast = ham_paths_start0_times_n(n, sigma_rotation(n))
    check("n=%d: full HK = fixed-start*n = %d" % (n, target),
          h_full == h_fast == target)

print()
print("2. EXACT R_rot at n = 15, 17, 19  (fixed-start Held-Karp, exact ints)")
exact_R = {}
for n in [15, 17, 19]:
    t0 = time.time()
    H = ham_paths_start0_times_n(n, sigma_rotation(n))
    R = H * 2 ** (n - 1) / math.factorial(n)
    exact_R[n] = R
    print("   n=%2d  H = %18d   R = %.5f   (%.1f s)" % (n, H, R, time.time() - t0))
check("exact R_rot(15) reproduces 2.48658 (consistency with first script)",
      abs(exact_R[15] - 2.48658) < 5e-6)
check("R_rot still EXCEEDS exp(tanh 1) at n = 17, 19 (naive prediction stays wrong)",
      exact_R[17] > math.exp(tanh1) and exact_R[19] > math.exp(tanh1),
      "R(17)=%.5f R(19)=%.5f" % (exact_R[17], exact_R[19]))

print()
print("3. SIS validation against exact values (mean +- 2*stderr)")
ok_sis = True
for n, target in [(13, 3711175 * 2 ** 12 / math.factorial(13)),
                  (15, None), (17, None), (19, None)]:
    tgt = exact_R.get(n, target)
    if tgt is None:
        continue
    Rm, Re = sis_R(n, sigma_rotation(n), 40000, 4200 + n)
    ok = abs(Rm - tgt) < max(Re, 3e-2)
    ok_sis = ok_sis and ok
    print("   n=%2d  SIS R = %.4f +- %.4f   exact %.5f   %s"
          % (n, Rm, Re, tgt, "ok" if ok else "MISMATCH"))
check("SIS estimator agrees with exact R within error bars (n=13,15,17,19)", ok_sis)

print()
print("4. SIS estimates at larger n (rotation family + Paley control)")
t_start = time.time()
results = []
plan = [("rot", 21, 80000), ("rot", 25, 60000), ("rot", 31, 50000),
        ("rot", 41, 30000), ("rot", 51, 25000),
        ("paley", 31, 40000), ("paley", 43, 25000)]
if True:
    for fam, n, samples in plan:
        sgl = sigma_rotation(n) if fam == "rot" else sigma_paley(n)
        t0 = time.time()
        Rm, Re = sis_R(n, sgl, samples, 91000 + n)
        results.append((fam, n, Rm, Re))
        print("   %-5s n=%3d  R = %.4f +- %.4f   (%d samples, %.1f s)"
              % (fam, n, Rm, Re, samples, time.time() - t0))
# optional big-n probes if we are still fast
if time.time() - t_start < 240:
    for fam, n, samples in [("rot", 75, 12000), ("rot", 101, 8000)]:
        sgl = sigma_rotation(n)
        t0 = time.time()
        Rm, Re = sis_R(n, sgl, samples, 91000 + n)
        results.append((fam, n, Rm, Re))
        print("   %-5s n=%3d  R = %.4f +- %.4f   (%d samples, %.1f s)"
              % (fam, n, Rm, Re, samples, time.time() - t0))

rot_seq = [(n, Rm, Re) for fam, n, Rm, Re in results if fam == "rot"]
pal_seq = [(n, Rm, Re) for fam, n, Rm, Re in results if fam == "paley"]
print()
print("   exact R_rot:  5..19 = 2.00000 2.22222 2.30476 2.38646 2.44113 2.48658"
      " %.5f %.5f" % (exact_R[17], exact_R[19]))
print("   reference: e = %.5f, exp(tanh 1) = %.5f" % (math.e, math.exp(tanh1)))
print()
print("VERDICT GUIDANCE (for the .out reader):")
print(" - if R_rot flattens near some constant c with exp(tanh1) < c < e:")
print("   the limit exists but is NOT given by the naive single-run exponential;")
print(" - if R_rot keeps climbing past e at n ~ 100: no finite limit at the")
print("   1/n-corrections scale we can see; Alon's polynomial window is live.")
print()
print("=" * 78)
print("SUMMARY: %d checks passed, %d failed" % (len(PASSED), len(FAILED)))
for f in FAILED:
    print("  FAILED: %s" % f)
print("=" * 78)
sys.exit(1 if FAILED else 0)
