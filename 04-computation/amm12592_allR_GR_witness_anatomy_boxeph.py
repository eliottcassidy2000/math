#!/usr/bin/env python3
"""ANGLE B3 step 4 -- y-space anatomy of winning witnesses: the cancellation law.

For each verified witness (R = 8..512), map rows to gamma_i = (Delta_i-(p-q))/2
and measure, per x-coefficient t of G_R:
    S_t := sum_i |[x^{t-i}] gamma_i|     (gross mass the construction routes)
    |g_t|                                 (net mass demanded)
    cancellation factor  rho_t = S_t / max(1,|g_t|)   (>= 1; log2 reported)
Also per-row: gamma_i(1) in {0,-1} (marker word), gamma_i(0) in {0,1},
max cell saturation |y|/box, count of nonzero cells.

If max_t log2(rho_t) grows ~linearly in R, the winning constructions use
exponentially deep cancellation (carries are heavy-duty); if it stays bounded,
an explicit bounded-cancellation scheme is plausible.  This decides the SHAPE
an all-R explicit construction must have.
"""
import json, os
from math import comb
from fractions import Fraction

WT = "/tmp/math-wt-boxeph-multifront"
CP = WT + "/04-computation"
OUT = WT + "/05-knowledge/results/amm12592_allR_GR_witness_anatomy_boxeph.out"

log_lines = []
def log(s=""):
    print(s); log_lines.append(s)
def flush():
    with open(OUT, "w") as f: f.write("\n".join(log_lines) + "\n")

def ptrim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a
def qpow(m): return [((-1)**k)*comb(m,k) for k in range(m+1)]
def bern_to_poly(cells, d):
    U = [0]*(d+1)
    U[0] = cells[d]
    for k in range(d-1, -1, -1):
        prev = 0
        for i in range(d+1):
            cur = U[i]; U[i] = cur - prev; prev = cur
        U[d-k] += cells[k]
    return U  # length d+1, untrimmed
def ballot(d):
    return [comb(d-1,k) - (comb(d-1,k-1) if k >= 1 else 0) for k in range(d+1)]
def Em(m): return ptrim([-1] + [1]*m)
def Gpoly(R):
    t = [a - b for a, b in
         zip(qpow(R-1) + [0]*R, (Em(R-1) + [0]*R))][:R]
    assert all(v % 2 == 0 for v in t)
    return [v//2 for v in t]

FILES = [
    ("amm12592_floor_witness_R8_ruleB_boxeph.json", "R8 ruleB"),
    ("amm12592_floor_witness_R16_ruleB_boxeph.json", "R16 ruleB"),
    ("amm12592_witness_R128_ruleB_D0_0_boxeph.json", "R128 ruleB D0=0"),
    ("amm12592_witness_R256_ruleA_D0_1_boxeph.json", "R256 ruleA D0=1"),
    ("amm12592_witness_R512_ruleA_D0_8_boxeph.json", "R512 ruleA D0=8"),
]
COMBINED = "amm12592_floor_witnesses_R8_R16_R32.json"

def anatomy(R, prof, blocks, label):
    G = Gpoly(R)
    # per-row gamma polys and stats
    rows = []
    maxsat_num, maxsat_den = 0, 1     # track max |y|/box as fraction
    nz_total = 0
    markers = []
    for i in range(R-1):
        d = prof[i]; b = ballot(d)
        y = []
        for t in range(d+1):
            diff = blocks[i][t] - b[t]
            assert diff % 2 == 0
            y.append(diff//2)
        # saturation
        for k, v in enumerate(y):
            if v == 0: continue
            box = comb(d-1, k) if v < 0 else comb(d-1, k-1) if k >= 1 else 0
            if box == 0 and v != 0:
                box = 1  # shouldn't happen for feasible
            if abs(v) * maxsat_den > maxsat_num * box:
                maxsat_num, maxsat_den = abs(v), box
        nz_total += sum(1 for v in y if v)
        gp = bern_to_poly(y, d)
        rows.append(gp)
        markers.append(sum(gp))
    # cancellation per coefficient
    T = R - 1 + max(prof[:R-1])
    worst = (0, 0)  # (log2rho approx as bitdiff, t)
    bitdiffs = []
    for t in range(len(G)):
        S = 0
        for i, gp in enumerate(rows):
            s = t - i
            if 0 <= s < len(gp):
                S += abs(gp[s])
        g = abs(G[t]) if t < len(G) else 0
        bd = (S.bit_length() - max(1, g).bit_length())
        bitdiffs.append(bd)
        if bd > worst[0]: worst = (bd, t)
    word = "".join("-" if m == -1 else ("0" if m == 0 else "?") for m in markers)
    n_minus = sum(1 for m in markers if m == -1)
    log("%s:" % label)
    log("  marker word (gamma_i(1)): %s%s" % (word[:100], "..." if len(word) > 100 else ""))
    log("  #(-1 rows) = %d (need R/2-1 = %d), all in {0,-1}: %s" %
        (n_minus, R//2 - 1, all(m in (0,-1) for m in markers)))
    log("  max cell saturation |y|/box = %d/%d ~ %.3g, nonzero cells = %d / %d" %
        (maxsat_num, maxsat_den, maxsat_num/maxsat_den if maxsat_den else 0.0,
         nz_total, sum(p+1 for p in prof[:R-1])))
    log("  cancellation: max_t log2(S_t/|g_t|) ~ %d bits at t = %d; mean bitdiff = %.1f"
        % (worst[0], worst[1], sum(bitdiffs)/len(bitdiffs)))
    # profile of bitdiff along t (coarse)
    step = max(1, len(bitdiffs)//16)
    prof_s = " ".join("%d" % bitdiffs[t] for t in range(0, len(bitdiffs), step))
    log("  bitdiff profile (t = 0, %d, ...): %s" % (step, prof_s))
    log("")
    flush()

comb_data = json.load(open(os.path.join(CP, COMBINED)))
for w in comb_data:
    anatomy(w["R"], w["profile"], w["blocks"], "combined R=%d" % w["R"])
for fn, label in FILES:
    try:
        d = json.load(open(os.path.join(CP, fn)))
        anatomy(d["R"], d["profile"], d["blocks"], label)
    except FileNotFoundError:
        log("[SKIP] %s" % fn)

log("interpretation: bitdiff ~ log2 of gross/net mass; linear-in-R growth =")
log("deep cancellation is intrinsic to winning constructions; bounded = shallow.")
flush()
