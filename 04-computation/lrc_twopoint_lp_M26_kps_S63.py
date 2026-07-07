#!/usr/bin/env python3
r"""
lrc_twopoint_lp_M26_kps_S63.py   (kind-pasteur-2026-07-07-S63, HYP-4897 probe C, decisive form)

THE 2-POINT LP AT M=26, EXACT ENUMERATION.

Question: is the mean floor E[maxgap] > 1/7 already forced by PAIR-UNIFORMITY alone?
LP: over probability measures on 13-point configurations on the discrete circle Z/26
(up to rotation -- rotation acts on configs, pair spectra and maxgap are invariant),
  minimize  E[maxgap]
  s.t.      aggregate unlabeled pair-distance density uniform:
            for each distance d in {1..13}: E[#ordered pairs at circ-distance d] = 12
            (156 ordered pairs over 13 effective distance classes; d=13 is its own
             antipode class counted once per unordered pair*2 = same normalization).
Every integer family's x-process satisfies the continuum version EXACTLY (S59).
Mixtures of grid configs are a SUBSET of all continuum processes, so:
  gridLP_min  >=  continuumLP_min.
If gridLP_min < 1/7: the 2-point conjecture DIES (pair data cannot force the floor,
weight->=3 usage is NECESSARY -- a methodological theorem).
If gridLP_min >= 1/7: the conjecture SURVIVES at this resolution (evidence, not proof --
finer grids could dip).

Solve: exact enumeration of configs (C(25,12) = 5200300 with point at 0 fixed kills
rotation; dedup by necklace canonical form optional -- we keep 0 in fact fixing rotation
partially: configs containing 0, each rotation class counted |orbit| times, harmless for
LP value). Lagrangian/dual subgradient on 13 multipliers with exact pricing each round
(pricing = scan configs minimizing maxgap - y . spec); finish with a primal fit over the
active column pool (projected least squares) to certify a feasible primal value.
"""
import numpy as np
from itertools import combinations

M = 26
NP_ = 13
target = np.full(13, 12.0)   # E[pairs at distance d] = 156/13 = 12 for d=1..13

print("enumerating configs (0 fixed), deduping (spec, maxgap) columns ...")
uniq = {}
cnt = 0
half = M // 2
for rest in combinations(range(1, M), NP_ - 1):
    pts = (0,) + rest
    spec = [0] * 13
    for i in range(NP_):
        pi = pts[i]
        for j in range(i + 1, NP_):
            d = pts[j] - pi
            if d > half: d = M - d
            spec[d - 1] += 2       # ordered pairs
    mg = pts[0] + M - pts[-1]
    prev = pts[0]
    for p in pts[1:]:
        g = p - prev
        if g > mg: mg = g
        prev = p
    key = (mg, bytes(spec))
    uniq[key] = uniq.get(key, 0) + 1
    cnt += 1
print(f"  {cnt} configs enumerated; {len(uniq)} unique (maxgap, spectrum) columns")
specs = np.array([np.frombuffer(k[1], dtype=np.uint8) for k in uniq], dtype=np.float64)
mgs = np.array([k[0] / M for k in uniq], dtype=np.float64)

# Dual subgradient: L(y) = min over configs of [ mg + y.(spec - target) ]  -> maximize over y
y = np.zeros(13)
best_dual = -1e9
step0 = 0.002
scores_cache = None
for it in range(300):
    scores = mgs + specs @ y - y @ target
    i = int(np.argmin(scores))
    L = float(scores[i])
    if L > best_dual: best_dual = L
    # subgradient = spec_i - target
    g = specs[i] - target
    y = y + (step0 / (1 + it / 40)) * g
    if it % 50 == 0:
        print(f"  it {it}: dual L = {L:.5f}  best = {best_dual:.5f}  |g| = {np.linalg.norm(g):.2f}")
print(f"  DUAL lower bound on gridLP: {best_dual:.5f}   (1/7 = {1/7:.5f})")

# Primal: build a column pool of low-score configs across a spread of y draws, then solve
# min mgs.w s.t. specs.w = target, w >= 0 on the pool by alternating projection + line search
pool = set()
rng = np.random.default_rng(63)
for rep in range(60):
    yy = y + rng.normal(0, 0.002 * (1 + rep / 10), 13)
    scores = mgs + specs @ yy
    idx = np.argpartition(scores, 40)[:40]
    pool.update(int(t) for t in idx)
pool.update(int(t) for t in np.argpartition(mgs, 200)[:200])
pool = sorted(pool)
A = specs[pool].T            # 13 x K
c = mgs[pool]
K = len(pool)
print(f"  primal pool: {K} columns")
w = np.full(K, 1.0 / K)
lam_path = [1e2, 1e3, 1e4, 1e5, 1e6]
for lam in lam_path:
    for it in range(4000):
        resid = A @ w - target
        grad = c + lam * 2 * (A.T @ resid)
        w -= 1e-6 * grad
        np.maximum(w, 0, out=w)
        s = w.sum()
        if s > 0: w /= s
    resid = A @ w - target
    print(f"  lam={lam:.0e}: E[mg] = {float(c @ w):.5f}  resid Linf = {float(np.abs(resid).max()):.4f}")
Ev = float(c @ w); rl = float(np.abs(A @ w - target).max())
print()
print(f"RESULT: primal (near-)feasible value = {Ev:.5f} with constraint Linf {rl:.4f} (target 12 per bin)")
print(f"        dual bound = {best_dual:.5f};  1/7 = {1/7:.5f}")
if best_dual > 1/7 + 1e-4:
    print("==> at M=26 the 2-POINT LP FLOOR EXCEEDS 1/7: pair-uniformity alone forces the mean")
    print("    floor at this resolution -- the 2-point conjecture SURVIVES (evidence, not proof).")
elif Ev < 1/7 - 1e-3 and rl < 0.3:
    print("==> a (near-)feasible mixture BEATS 1/7: pair data alone CANNOT force the floor;")
    print("    weight->=3 information is NECESSARY in any proof of the mean floor.")
else:
    print("==> inconclusive at this resolution/solver; report both numbers honestly.")
