# Auditor numerics-prereg: what was already in the repository about the registered window [2^27, 2^28) and about
# the cumulative census at 2^28 BEFORE the registration commit 73c8437aac (2026-11-11 11:29:15)?
#  (1) S91 audit, recompute-5128/a6_window_27_28.out (commit 49e15115c6, 10:53:11): measured excess E of the
#      distinguished class modulo the primes 3, 5, 7 on the SAME window.  With class factors f (p * share):
#      f_other = p / (p + E), f_distinguished = p (1 + E) / (p + E).
#      f_other is exactly the g0 = 1 group mean registered in P1 at m = p^k; f_distinguished is the class-weighted
#      mean of the groups g0 > 1.
#  (2) S91, first_order_exact_28.out (commit 6cf358fd9c, 07:59:52): remainder p max|R_small| of the first-order law
#      at primes for X = 2^26 and 2^28: the prime-level analogue of the registered P4.
import re
import json
from aud_common import divisors, groups, remainder
from analyze import load

A6 = "C:/tmp/wt-s92/05-knowledge/results/s91_f1_local_densities_20261110/audit/recompute-5128/a6_window_27_28.out"
FO = "C:/tmp/wt-s92/04-computation/experiments/s91_f1_local_densities_20261110/first_order_exact_28.out"
PRE = "C:/tmp/wt-s92/04-computation/experiments/s92_prime_powers_20261111/prereg_window_28_27.json"

E = {}
for ln in open(A6):
    parts = ln.replace("|", " ").split()
    if len(parts) >= 11 and parts[0] in ("3", "5", "7"):
        # columns: p, f_pred s0, f_pred s1, f_meas s0, f_meas s1, E_pred s0, s1, E_claimed s0, s1, E_meas s0, s1, ...
        p = int(parts[0])
        E[p] = (float(parts[9]), float(parts[10]))
print("(1) measured excess E of the distinguished class on [2^27, 2^28), from the S91 audit output (before registration):", E)
P = json.load(open(PRE))
J, T, mods = load("my_win_28_27.json")
print("    m    sheet | g0 = 1 mean: from the S91 audit number p/(p+E) | measured in the registered run | registered prediction ||"
      " mean over the classes with g0 > 1: from the audit p(1+E)/(p+E) | measured | from the registered predictions")
for m, p in ((9, 3), (27, 3), (81, 3), (25, 5), (125, 5), (49, 7)):
    for sh in range(2):
        G = groups(m, *mods[m], sh, False)
        e = E[p][sh]
        f_other = p / (p + e)
        f_dist = p * (1 + e) / (p + e)
        num = sum(G[g0][0] * G[g0][5] for g0 in divisors(m) if g0 > 1)
        den = sum(G[g0][5] for g0 in divisors(m) if g0 > 1)
        nump = sum(P["phi"]["%d,%d,%d" % (sh, m, g0)] * G[g0][5] for g0 in divisors(m) if g0 > 1)
        print("  %4d    %d   |  %.4f | %.4f | %.4f  ||  %.4f | %.4f | %.4f" % (
            m, sh, f_other, G[1][0], P["phi"]["%d,%d,1" % (sh, m)], f_dist, num / den, nump / den))
print("    -> the g0 = 1 group means of the registered window at 9, 27, 81, 25, 125, 49 (12 of the 68 registered group values) and the")
print("       class-weighted mean of the other groups were fixed by numbers published before the registration.")

print()
print("(2) prime-level remainder p max|R_small| at 2^26 -> 2^28 (S91 first_order_exact_28.out, before registration) against the")
print("    registered P4 pairs (remainder of (B) at the prime powers, my census)")
rows = {}
X = None
for ln in open(FO):
    mo = re.match(r"X = 2\^(\d+)", ln)
    if mo:
        X = int(mo.group(1))
        continue
    parts = ln.split()
    if len(parts) == 7 and parts[0] in ("0", "1") and parts[1] in ("3", "5", "7"):
        rows[(X, int(parts[0]), int(parts[1]))] = float(parts[5])
J26, T26, m26 = load("my_cum_26.json")
J28, T28, m28 = load("my_cum_28.json")
agree = 0
for p, ms in ((3, (9, 27, 81)), (5, (25, 125)), (7, (49,))):
    for sh in range(2):
        a, b = rows[(26, sh, p)], rows[(28, sh, p)]
        cells = []
        for m in ms:
            r26 = remainder(m, *m26[m], sh, False)
            r28 = remainder(m, *m28[m], sh, False)
            same = (b < a) == (r28 < r26)
            agree += same
            cells.append("m = %d: %.5f -> %.5f (%s)%s" % (m, r26, r28, "smaller" if r28 < r26 else "not smaller", "" if same else " [direction differs from the prime]"))
        print("    p = %d sheet %d: %.5f -> %.5f (%s)   |  %s" % (p, sh, a, b, "smaller" if b < a else "not smaller", ";  ".join(cells)))
print("    the direction known at the prime before the registration equals the direction at the prime power in %d of the 12 registered pairs;" % agree)
npred = sum(1 for p, ms in ((3, (9, 27, 81)), (5, (25, 125)), (7, (49,))) for sh in range(2) for m in ms if rows[(28, sh, p)] < rows[(26, sh, p)])
print("    copying the prime-level directions would have predicted 'smaller' for %d of the 12 pairs (registered threshold: 8; outcome: 8)." % npred)
