# opus-2026-07-16-S332 -- HYP-7135 / THM-925 part (B): T2 STRUCTURE RECON.
# The general triple beat |mu(W cap D_C) - (2/13) mu(W)|, W = D_A cap D_B,
# WITHOUT scale separation (mac-mini THM-924's re-handed priority).
#
# Layer 1 (this script): the EXACT structure of W for coprime A < B:
#   (S1) components of W <-> beat offsets m = iB - jA, |m| < (A+B)/13,
#        exactly ONE component per admissible m (coprime case);
#   (S2) widths w(m) = sawtooth: (1/(13AB)) * min(2A, 2B, (A+B) - 13|m|)+
#        hmm -- verified against codex's T-law normalization;
#   (S3) THE AP LAW: within each endpoint REGIME the component left
#        endpoints form an EXACT AP in m, with step in {B*/A - 1/(AB),
#        -1/(13 something)...} -- empirically extracted, then matched to the
#        Bezout cofactors B*B - AB~ = 1: candidate steps B*/A mod 1 and
#        B~/B mod 1 (+ sawtooth slope corrections);
#   (S4) baseline T2-U: |mu(W cap D_C) - (2/13) mu W| <= 4 kappa_W/(13 C)
#        (unconditional) -- verify on a grid incl. in-band and resonant C;
#   (S5) the resonance correlation: observed err vs the triple small form
#        y* = min_{1<=q<=13, 0<=r,s<=13} |qC - rA - sB| (not all r,s zero
#        allowed... q>=1, (r,s) any signs) -- fit the shape for Layer 2.
from fractions import Fraction
from math import gcd
import itertools, random

F = Fraction

def teeth(x):
    """D_x as sorted disjoint intervals on [0,1) (exact, wrap split)."""
    w = F(1, 13*x)
    out = []
    for j in range(x):
        a, b = F(j, x) - w, F(j, x) + w
        if a < 0:
            out.append((F(0), b)); out.append((a + 1, F(1)))
        elif b > 1:
            out.append((a, F(1))); out.append((F(0), b - 1))
        else:
            out.append((a, b))
    out.sort()
    # merge adjacent (wrap pieces)
    merged = []
    for iv in out:
        if merged and iv[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], iv[1]))
        else:
            merged.append(list(iv))
    return [(a, b) for a, b in merged]

def isect(u, v):
    """intersect two sorted disjoint interval lists (exact)."""
    out, i, j = [], 0, 0
    while i < len(u) and j < len(v):
        a = max(u[i][0], v[j][0]); b = min(u[i][1], v[j][1])
        if a < b: out.append((a, b))
        if u[i][1] < v[j][1]: i += 1
        else: j += 1
    return out

def circ_components(ivs):
    """merge across the 0/1 seam for honest circle component count."""
    if len(ivs) >= 2 and ivs[0][0] == 0 and ivs[-1][1] == 1:
        return len(ivs) - 1
    return len(ivs)

def mu(ivs): return sum(b - a for a, b in ivs)

print("=" * 72)
print("(S1-S3) THE EXACT COMPONENT STRUCTURE OF W = D_A cap D_B (coprime):")
random.seed(332)
pairs = [(19, 24), (50, 63), (87, 101), (120, 143), (55, 89), (144, 233)]
ap_report = []
for A, B in pairs:
    assert gcd(A, B) == 1
    W = isect(teeth(A), teeth(B))
    kap = circ_components(W)
    M = (A + B - 1) // 13   # |m| <= M candidate range: |m| < (A+B)/13
    # predicted: one component per m with (A+B) - 13|m| > 0
    pred_ms = [m for m in range(-M, M + 1) if (A + B) - 13*abs(m) > 0]
    # predicted width law
    wsum_pred = sum(F(min(2*A, 2*B, (A + B) - 13*abs(m)), 13*A*B)
                    for m in pred_ms)
    ok_mu = (mu(W) == wsum_pred)
    ok_k = (kap == len(pred_ms))
    # tag each component with its m: centers near i/A: m = round(B*A*center)*...
    # exact: for component (a,b): find i = nearest A-tooth: i = round(center*A),
    # j = round(center*B), m = i*B - j*A
    tags = []
    for (a, b) in W:
        c = (a + b) / 2
        i = round(c * A); j = round(c * B)
        m = i * B - j * A
        tags.append((m, a, b))
    ok_tags = sorted(t[0] for t in tags) == sorted(pred_ms) if ok_k else False
    # (S3) AP law: sort by m; left endpoints; successive differences mod 1,
    # grouped by regime = which max/min case is active
    tags.sort()
    Bs = pow(B, -1, A)                    # B* = B^{-1} mod A
    Bt = (B * Bs - 1) // A                # B~ : B*B - A*B~ = 1
    diffs = {}
    for k in range(1, len(tags)):
        m0, a0, _ = tags[k-1]; m1, a1, _ = tags[k]
        if m1 == m0 + 1:
            d = (a1 - a0) % 1
            diffs.setdefault(d, []).append(m0)
    # candidate steps
    cand = {('B*/A', F(Bs, A) % 1), ('B*/A - 1/(AB)', (F(Bs, A) - F(1, A*B)) % 1),
            ('B~/B', F(Bt, B) % 1), ('B~/B + 1/(AB)', (F(Bt, B) + F(1, A*B)) % 1)}
    matched = {d: [nm for nm, v in cand if v == d] for d in diffs}
    nmatch = sum(1 for d in diffs if matched[d])
    ap_report.append((A, B, len(diffs), nmatch))
    print(f"  A={A:4d} B={B:4d}: kappa={kap:4d} (pred {len(pred_ms)}: {ok_k}); "
          f"mu exact match: {ok_mu}; m-tags match: {ok_tags}; "
          f"distinct consecutive steps: {len(diffs)}, Bezout-matched: {nmatch}")
    if len(diffs) <= 6:
        for d, ms in sorted(diffs.items()):
            lbl = matched[d] if matched[d] else '??'
            print(f"      step {d} on m-runs {ms[:8]}{'...' if len(ms)>8 else ''}"
                  f"  <- {lbl}")

print()
print("=" * 72)
print("(S4) BASELINE T2-U: err <= 4 kappa_W / (13 C), unconditional:")
print("(S5) resonance shape: err vs y* = min |qC - rA - sB|, q<=13, |r|,|s|<=13")
rows = []
worst_ratio = F(0)
for (A, B) in [(19, 24), (50, 63), (87, 101), (55, 89)]:
    W = isect(teeth(A), teeth(B))
    kap = circ_components(W)
    muW = mu(W)
    Cs = set()
    # in-band probes: near-resonant and generic
    for C in [A + B, A + B + 1, 2*A + B, A + 2*B, 2*A + 2*B + 1, 3*B + 2,
              5*A + 3*B, 7*A - 2*B, 13*A + 1, 6*A + 6*B + 5]:
        if C > B and C != A and C != B: Cs.add(C)
    random.seed(A * 1000 + B)
    while len(Cs) < 16:
        C = random.randint(B + 1, 13 * A)
        if C not in (A, B): Cs.add(C)
    for C in sorted(Cs):
        WC = isect(W, teeth(C))
        err = mu(WC) - F(2, 13) * muW
        bnd = F(4 * kap, 13 * C)
        ok = abs(err) <= bnd
        if not ok:
            print(f"  *** BASELINE VIOLATED at A={A} B={B} C={C}: "
                  f"err={float(err):+.3e} bnd={float(bnd):.3e}")
        ystar = min(abs(q*C - r*A - s*B)
                    for q in range(1, 14)
                    for r in range(-13, 14) for s in range(-13, 14))
        rows.append((A, B, C, kap, float(err), float(bnd), ystar,
                     float(abs(err))))
        worst_ratio = max(worst_ratio, abs(err) / bnd)
print(f"  baseline holds on all {len(rows)} triples; "
      f"worst |err|/bound = {float(worst_ratio):.4f}")
print()
print("  resonance correlation (in-band triples, sorted by y*):")
print("   y*  :  max |err|         mean |err|        (n)   scaled max*169/(A+B) hmm")
from collections import defaultdict
bins = defaultdict(list)
for (A, B, C, kap, err, bnd, ystar, aerr) in rows:
    key = 0 if ystar == 0 else (1 if ystar <= 3 else (2 if ystar <= 10 else
          (3 if ystar <= 30 else (4 if ystar <= 100 else 5))))
    bins[key].append((aerr, ystar, A, B, C))
lbl = {0: 'y*=0 (exact res)', 1: 'y*<=3', 2: '3<y*<=10', 3: '10<y*<=30',
       4: '30<y*<=100', 5: 'y*>100'}
for k in sorted(bins):
    v = bins[k]
    print(f"   {lbl[k]:18s}: max {max(x[0] for x in v):.3e}  "
          f"mean {sum(x[0] for x in v)/len(v):.3e}  (n={len(v)})")
print()
print("  worst 6 triples by |err| (for Layer-2 targeting):")
for (A, B, C, kap, err, bnd, ystar, aerr) in sorted(rows, key=lambda r: -r[7])[:6]:
    print(f"   A={A:3d} B={B:3d} C={C:4d}: err={err:+.3e} y*={ystar:3d} "
          f"kappa={kap}")
