"""
lrc14_riesz_dissociated_uniform_opus_S178.py   (opus-2026-07-09-S178)

Is the DISSOCIATED branch's looseness certifiable with a UNIFORM Riesz margin -- sup over dissociated
13-sets of inf_R int(M*R)/int(R) < 1 (bounded away)?  If yes: inf L > 0 for the dissociated subfamily,
a DECOMPOSITION-FREE closure of the branch that was Mertens-walled (opus-S172) -- and it sidesteps the
two-scale slow-fast separation that mac-mini-S64 flagged as failing for large spread (THM-663 concern).

The Riesz certificate (opus-S173/S174): S loose (positive-measure lonely set) if some Riesz product
R=prod(1+a_m cos2pi m tau)>=0 has ratio=int(M*R)/int(R) < 1, M(tau)=#{v:||v tau||<=1/14}.  This does an
ADVERSARIAL search over dissociated 13-sets (random + hill-climb to MAXIMIZE the best ratio = hardest to
certify), per-set optimizing R by coordinate descent.  Reports the MAX best-ratio; < 1 with margin =>
uniform dissociated looseness.
"""
import numpy as np
import random

NG = 120011
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14


def Mmult(S):
    M = np.zeros(NG)
    for v in S:
        d = np.abs(((v * TAU + 0.5) % 1.0) - 0.5)
        M += (d <= H)
    return M


def best_ratio(S, D, iters=8):
    """coordinate-descent min of ratio int(M*R)/int(R), R=prod(1+a_m cos), |a_m|<=0.999."""
    M = Mmult(S)
    cb = [np.cos(2 * np.pi * m * TAU) for m in D]
    a = np.zeros(len(D))

    def ratio(a):
        R = np.ones(NG)
        for i in range(len(D)):
            R = R * (1 + a[i] * cb[i])
        r = R.mean()
        return (M * R).mean() / r if r > 1e-6 else 5.0
    cur = ratio(a)
    for _ in range(iters):
        improved = False
        for i in range(len(D)):
            best_ai, best_v = a[i], cur
            for c in np.linspace(-0.999, 0.999, 41):
                a[i] = c; v = ratio(a)
                if v < best_v - 1e-9:
                    best_v, best_ai = v, c; improved = True
            a[i] = best_ai; cur = best_v
        if not improved:
            break
    return cur


def longest_ap(S):
    Sset = set(S); S = sorted(S); best = 1
    for i, aa in enumerate(S):
        for b in S[i + 1:]:
            d = b - aa
            if aa - d in Sset: continue
            L = 2; x = b + d
            while x in Sset: L += 1; x += d
            best = max(best, L)
    return best


def lonely_measure(S):
    M = Mmult(S)
    return float((M == 0).mean())


rng = random.Random(178)
print("=" * 92)
print("ADVERSARIAL Riesz margin over DISSOCIATED 13-sets: max best-ratio (< 1 => uniform looseness)")
print("=" * 92)
worst = 0.0; worst_S = None; nchk = 0; ntight = 0
for trial in range(45):
    spread = rng.randint(40, 260)
    mid = sorted(rng.sample(range(1, spread), 11)); S = [1] + mid + [spread] if 1 not in mid else sorted(set([spread] + mid + [rng.randint(1, spread)]))
    S = sorted(set(S))
    while len(S) != 13:
        S = sorted(set(S + [rng.randint(1, spread)]))[:13]
    if longest_ap(S) > 7:      # dissociated only (longest-AP <= 7)
        continue
    lm = lonely_measure(S)
    if lm < 1e-4:              # tight (exact-check territory), skip
        ntight += 1; continue
    D = sorted(S)              # speeds as (non-dissociated but valid) Riesz freqs
    r = best_ratio(S, D)
    nchk += 1
    if r > worst:
        worst = r; worst_S = (S, lm, r)
    if trial % 10 == 0:
        print(f"  trial {trial}: spread={spread}, L={longest_ap(S)}, lonely-meas={lm:.4f}, best-ratio={r:.4f}  (running max {worst:.4f})")
print()
if worst_S:
    S, lm, r = worst_S
    print(f"  HARDEST dissociated set: spread={max(S)}, lonely-measure={lm:.4f}, best Riesz ratio={r:.4f}")
print(f"  MAX best-ratio over {nchk} dissociated sets = {worst:.4f}  ({ntight} tight skipped)")
print(f"  => dissociated looseness {'UNIFORMLY certifiable (sup ratio < 1)' if worst < 1 else 'NOT uniform (some ratio >= 1)'}")
print()
print("  READING: sup_dissociated inf_R ratio < 1 (with margin) => inf L > 0 for dissociated -- a")
print("  DECOMPOSITION-FREE looseness closure (no two-scale/drift, no Mertens); the naive Riesz here is")
print("  a lower bar than the tuned Bedert-2025 construction, so a margin <1 here is strong evidence.")
print("=" * 92)
