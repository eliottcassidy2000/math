#!/usr/bin/env python3
"""
lrc14_nonisolated_truedisc_sweep_opus_S271.py
=============================================
opus-2026-07-13-S271.  THE NON-ISOLATED RESIDUAL vs THE TRUE-DISC CERTIFICATE.

State after kps-S128 + klein-S289: the far-element (ISOLATED) regime of Route B is CLOSED
(exact-Q Dedekind certificates, both extremal rays, all {1..11,a,b}); the crude arc-count bound
r < 3*sqrt2*v*|G'| is FALSE in general (klein-S289 census) and the open residual is the
NON-ISOLATED class: (i) compact / mid-band bodies (no far element at all, e.g. the compressed
near-dilate 2{1..12}u{13}), (ii) far-CLUSTER bodies ({1,90..101}: large diameter but the top is
a block, v/max(W) ~ 1).

klein-S289 only showed the CRUDE bound fails there.  Nobody has asked whether the TRUE disc_v
(THM-731, now exact-Q decidable per body via kps's THM-732) certifies at some peel.  S270 found
YES for 2{1..12}u{13} (8/13 peels).  This sweep answers it for the whole accessible class:

  for each non-isolated covering body, all 13 peels: M731(v) = 6|G'_{~v}| - sqrt(6 disc_v)
  (L >= M731/7; certifies iff > 0).  Report: certified? best peel; margin; isolation ratio.

PERSPECTIVE FRAME: a peel certifies iff that runner's clock resolves the peeled-away good set.
Isolated far elements are SWEEPERS (klein's classifier = the frozen-fan/sweeper dichotomy,
opus-S270).  The question here: in a body with NO sweeper, is some perspective still fine
enough?  If YES uniformly -> the residual = proving the Dedekind/exposure collapse (kps
HYP-6495) at the certifying peels.  If some body fails ALL 13 peels under TRUE disc -> THM-731
itself (the Cauchy-Schwarz step) is insufficient there = a new, sharper hard exemplar.
"""
import numpy as np
import math
from itertools import combinations

NG = 1 << 21
THR = 1.0 / 14.0

def good_indicator(S):
    t = np.arange(NG, dtype=np.float64) / NG
    g = np.ones(NG, dtype=np.float64)
    for w in S:
        frac = (w * t) % 1.0
        g *= (np.minimum(frac, 1.0 - frac) >= THR)
    return g

def analyze_body(body):
    """all-peel true-disc certificates; returns rows + true L"""
    Lfull = good_indicator(body).mean()
    rows = []
    for v in sorted(body):
        rest = [w for w in body if w != v]
        g = good_indicator(rest)
        Gp = g.mean()
        # interval count r on grid
        r = int(np.sum((g[1:] > 0.5) & (g[:-1] < 0.5)) + (1 if (g[0] > 0.5 and g[-1] < 0.5) else 0))
        G = np.fft.rfft(g)
        A = np.fft.irfft(G * np.conj(G), n=NG) / NG
        idx = np.round((np.arange(v) / v) * NG).astype(np.int64) % NG
        disc = float(A[idx].mean() - A.mean())
        M731 = 6 * Gp - math.sqrt(max(6 * disc, 0.0))
        M732 = 6 * Gp - math.sqrt(2) * r / v
        rows.append((v, r, Gp, disc, M731, M732))
    return Lfull, rows

def is_covering(body):
    return all(any(w % q == 0 for w in body) for q in range(2, 15))

def isolation(body):
    s = sorted(body)
    return s[-1] / s[-2]

# ---------------- the battery ----------------

NAMED = [
    ("compressed 2{1..12}u{13}",   [2,4,6,8,10,12,14,16,18,20,22,24,13]),
    ("M141 extremal M=3/37",       [2,3,5,8,9,11,12,13,14,15,17,20,23]),
    ("klein cex far-cluster",      [1] + list(range(90, 102))),
    ("klein cex {1,2,3,50..59}",   [1,2,3] + list(range(50, 60))),
    ("klein limitation {1,3..14}", [1] + list(range(3, 15))),
    ("consecutive {2..14}",        list(range(2, 15))),
]

# random primitive DC bodies, Vmax<=26, non-isolated (max/secondmax <= 1.5)
rng = np.random.default_rng(20260713)
pool = list(range(1, 27))
extra = []
seen = set()
tries = 0
while len(extra) < 12 and tries < 200000:
    tries += 1
    body = sorted(rng.choice(pool, size=13, replace=False).tolist())
    key = tuple(body)
    if key in seen: continue
    seen.add(key)
    if not is_covering(body): continue
    if math.gcd(*body) if len(body)==2 else np.gcd.reduce(body) != 1: continue
    if isolation(body) > 1.5: continue
    extra.append(("random DC Vmax<=26 #%d" % (len(extra)+1), body))

print("=" * 108)
print("NON-ISOLATED SWEEP: true-disc THM-731 all-peel certificates.  L >= M731/7 ; certifies iff M731>0")
print("(crude M732 shown for contrast; klein-S289: crude fails on non-isolated).  NG=2^21")
print("=" * 108)

summary = []
for name, body in NAMED + extra:
    L, rows = analyze_body(body)
    ncert = sum(1 for (_, _, _, _, M731, _) in rows if M731 > 0)
    ncert732 = sum(1 for (_, _, _, _, _, M732) in rows if M732 > 0)
    best = max(rows, key=lambda x: x[4])
    iso = isolation(body)
    summary.append((name, body, L, ncert, ncert732, best, iso))
    print("\n%-28s  iso=%.2f  true L=%.6f   cert-731 at %d/13 peels (crude at %d/13)"
          % (name, iso, L, ncert, ncert732))
    if ncert == 0:
        print("   *** FAILS ALL PEELS under TRUE disc -- new hard exemplar for THM-731 ***")
        for v, r, Gp, disc, M731, M732 in rows:
            print("     peel %4d  r=%3d  |G'|=%.5f  disc=%.3e  M731=%+.4f" % (v, r, Gp, disc, M731))
    else:
        v, r, Gp, disc, M731, M732 = best
        print("   best peel v=%d (r=%d, |G'|=%.5f, disc=%.3e): M731=%.5f -> L >= %.6f  [true %.6f]"
              % (v, r, Gp, disc, M731, M731 / 7, L))

print("\n" + "=" * 108)
print("SUMMARY   (iso = max/secondmax; klein classifier: crude fires iff isolated)")
print("%-28s %6s %10s %12s %12s %10s" % ("body", "iso", "true L", "731 peels", "732 peels", "best M731"))
for name, body, L, ncert, ncert732, best, iso in summary:
    print("%-28s %6.2f %10.6f %12s %12s %10.5f"
          % (name, iso, L, "%d/13" % ncert, "%d/13" % ncert732, best[4]))

nfail = sum(1 for s in summary if s[3] == 0)
print("\nbodies failing ALL peels under TRUE disc: %d / %d" % (nfail, len(summary)))
print("done.")
