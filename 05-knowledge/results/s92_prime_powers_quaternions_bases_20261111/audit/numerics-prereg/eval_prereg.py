# Auditor numerics-prereg: independent evaluation of the registered tests P1-P4.
#  (a) from MY census (my_win_28_27.json, my_cum_26.json, my_cum_28.json), exact rationals -> floats
#  (b) from the STORED session outputs (window_prereg_28_27.out, cum_26.out, cum_28.out), as printed
# and an assessment of the discriminating power of the registered tolerances.
import sys
import json
import math
from aud_common import divisors, remainder, groups, L_m
from analyze import load, parse_out

S = "C:/tmp/wt-s92/04-computation/experiments/s92_prime_powers_20261111/"
P = json.load(open(S + "prereg_window_28_27.json"))
MODS = P["moduli"]
tol = {int(k): v for k, v in P["tol"].items()}
flat = {int(k): v for k, v in P["flat"].items()}

print("registered moduli:", MODS)
print("registered tol:", tol)
print("registered flat:", flat)

# ---------- (a) my census
J, T, mods = load("my_win_28_27.json")
print()
print("=== (a) from my own census of the window [2^27, 2^28): triples", T)
p1_ok = True
p2_ok = True
p1_cases = 0
p3_hits = 0
p3_n = 0
rows = []
for m in MODS:
    for sh in range(2):
        G = groups(m, *mods[m], sh, False)
        csig = -1 if sh == 0 else 0
        Tsh = sum(mods[m][0][sh])
        worst_dev = 0.0
        worst_flat = 0.0
        for g0 in divisors(m):
            mean, mn, mx, phi_meas, lim, ncls = G[g0]
            pred = P["phi"]["%d,%d,%d" % (sh, m, g0)]
            worst_dev = max(worst_dev, abs(mean - pred))
            worst_flat = max(worst_flat, mx - mean, mean - mn)
            if m in (9, 27, 25) and g0 > 1:
                p3_n += 1
                if abs(mean - pred) < abs(mean - lim):
                    p3_hits += 1
        p1_cases += 1
        ok1 = worst_dev <= tol[m]
        ok2 = worst_flat <= flat[m]
        p1_ok = p1_ok and ok1
        p2_ok = p2_ok and ok2
        rows.append((m, sh, worst_dev, tol[m], ok1, worst_flat, flat[m], ok2))
for r in rows:
    print("  m = %3d sheet %d: max |group mean - registered prediction| = %.4f (tol %.3f, %s, used %.0f%% of tol);  max class deviation from group mean = %.4f (flat %.3f, %s, used %.0f%%)" % (
        r[0], r[1], r[2], r[3], "ok" if r[4] else "FAIL", 100 * r[2] / r[3], r[5], r[6], "ok" if r[7] else "FAIL", 100 * r[5] / r[6]))
print("P1 (16 cases): %s;  P2 (16 cases): %s;  P3: %d of %d -> %s" % (
    "CONFIRMED" if p1_ok else "FAILED", "CONFIRMED" if p2_ok else "FAILED", p3_hits, p3_n, "CONFIRMED" if p3_hits == p3_n else "FAILED"))

# P4
J26, T26, m26 = load("my_cum_26.json")
J28, T28, m28 = load("my_cum_28.json")
cnt = 0
print("P4 (remainder of (B) at 2^28 against 2^26), from my census; noise = sqrt(m/#T) of the same sheet:")
for m in (9, 27, 81, 25, 125, 49):
    for sh in range(2):
        r26 = remainder(m, *m26[m], sh, False)
        r28 = remainder(m, *m28[m], sh, False)
        n26 = math.sqrt(m / T26[sh])
        n28 = math.sqrt(m / T28[sh])
        smaller = r28 < r26
        cnt += smaller
        z = (r28 - r26) / math.sqrt(n26 ** 2 + n28 ** 2)
        print("  m = %3d sheet %d: 2^26 %.5f -> 2^28 %.5f  change %+.5f  (%s)   noise 2^26 %.5f, 2^28 %.5f; change / combined noise = %+.1f" % (
            m, sh, r26, r28, r28 - r26, "smaller" if smaller else "NOT smaller", n26, n28, z))
print("P4: smaller in %d of 12 pairs (registered threshold: at least 8) -> %s" % (cnt, "CONFIRMED" if cnt >= 8 else "FAILED"))
ptail = sum(math.comb(12, k) for k in range(8, 13)) / 2 ** 12
print("   (for reference: P(at least 8 of 12) under independent fair coin flips = %.4f)" % ptail)

# ---------- (b) from stored outputs
print()
print("=== (b) from the stored session outputs as printed")
W = parse_out(S + "window_prereg_28_27.out")
p1b = True
p2b = True
p3b = 0
p3bn = 0
for (m, sh), d in W["groups"].items():
    for g0, (mean, mn, mx, rest) in d.items():
        mean = float(mean)
        mn = float(mn)
        mx = float(mx)
        apri = float(rest[0])
        lim = float(rest[2])
        if abs(mean - apri) > tol[m]:
            p1b = False
            print("  P1 violation", m, sh, g0)
        if max(mx - mean, mean - mn) > flat[m]:
            p2b = False
            print("  P2 violation", m, sh, g0)
        if m in (9, 27, 25) and g0 > 1:
            p3bn += 1
            if abs(mean - apri) < abs(mean - lim):
                p3b += 1
print("P1: %s; P2: %s; P3: %d of %d" % ("CONFIRMED" if p1b else "FAILED", "CONFIRMED" if p2b else "FAILED", p3b, p3bn))
C26 = parse_out(S + "cum_26.out")
C28 = parse_out(S + "cum_28.out")
cnt = 0
grew = []
for m in (9, 27, 81, 25, 125, 49):
    for sh in range(2):
        a = float(C26["table"][m][sh])
        b = float(C28["table"][m][sh])
        if b < a:
            cnt += 1
        else:
            grew.append("%d sheet %d" % (m, sh))
print("P4 from cum_26.out / cum_28.out: smaller in %d of 12; not smaller at: %s" % (cnt, ", ".join(grew)))

# ---------- discriminating power of P1
print()
print("=== discriminating power of the registered P1 tolerances (would other laws pass P1 on this window?)")
print("    alternatives evaluated against the MEASURED group means of the window:")
print("    FLAT = every class factor 1;  LIMIT = L_m(g0) of THM-5134 (D);  PRIME-LEVEL = factor depends only on the")
print("    primes dividing g0 (the measured mean over all classes sharing the same radical of g0, i.e. the law of THM-5128 (C)).")
for m in MODS:
    for sh in range(2):
        G = groups(m, *mods[m], sh, False)
        dev_flat = 0.0
        dev_lim = 0.0
        dev_rad = 0.0
        # prime-level alternative: pool the classes with the same set of primes dividing g0
        pools = {}
        for g0 in divisors(m):
            rad = tuple(p for p in (3, 5, 7) if g0 % p == 0)
            mean, mn, mx, phi_meas, lim, ncls = G[g0]
            pools.setdefault(rad, []).append((mean, ncls))
        poolmean = {rad: sum(a * n for a, n in v) / sum(n for a, n in v) for rad, v in pools.items()}
        for g0 in divisors(m):
            mean, mn, mx, phi_meas, lim, ncls = G[g0]
            rad = tuple(p for p in (3, 5, 7) if g0 % p == 0)
            dev_flat = max(dev_flat, abs(mean - 1.0))
            dev_lim = max(dev_lim, abs(mean - lim))
            dev_rad = max(dev_rad, abs(mean - poolmean[rad]))
        # also: does the registered prediction itself separate FLAT from the law within the tolerance?
        predmax = max(abs(P["phi"]["%d,%d,%d" % (sh, m, g0)] - 1.0) for g0 in divisors(m))
        print("  m = %3d sheet %d: tol %.3f | max |measured - 1| = %.4f (%s) | max |measured - limit| = %.4f (%s) | max |measured - prime-level mean| = %.4f (%s) | max |registered prediction - 1| = %.4f (%s)" % (
            m, sh, tol[m], dev_flat, "FLAT would pass" if dev_flat <= tol[m] else "flat rejected",
            dev_lim, "LIMIT would pass" if dev_lim <= tol[m] else "limit rejected",
            dev_rad, "PRIME-LEVEL would pass" if dev_rad <= tol[m] else "prime-level rejected",
            predmax, "prediction within tol of FLAT everywhere: P1 cannot reject flat data" if predmax <= tol[m] else "ok"))
