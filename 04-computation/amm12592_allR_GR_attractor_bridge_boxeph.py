#!/usr/bin/env python3
"""ANGLE B3 step 6 -- the attractor bridge: sigma_i = E_{R-2-i} + 2 T_{i+1}.

Delta-space residuals: sigma_{-1} = q^{R-1}, p sigma_i = sigma_{i-1} - Delta_i.
y-space tail targets:  T_0 = G_R,      x T_{i+1} = T_i - gamma_i.
Claim (derived from the master equation + (p-q)(1+...+x^i) = E_i + 2x^{i+1}):
    sigma_i = E_{R-2-i} + 2 T_{i+1}      for -1 <= i <= R-2,
with the convention E_{-1} = -1 ... actually E_m = -1+x+...+x^m, E_{-1} = -1?
E_{-1} means empty positive part: E_{-1} = -1.  We verify on witnesses.
So Rule A's attractor sigma_i ~ "E-shape" is EXACTLY the statement T_{i+1} ~ 0:
the y-space tail is HALF the deviation of sigma_i from the E-attractor.
"""
import json, os
from math import comb

WT = "/tmp/math-wt-boxeph-multifront"
CP = WT + "/04-computation"
OUT = WT + "/05-knowledge/results/amm12592_allR_GR_attractor_bridge_boxeph.out"
log_lines = []
def log(s=""):
    print(s); log_lines.append(s)
def flush():
    with open(OUT, "w") as f: f.write("\n".join(log_lines) + "\n")

def ptrim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a
def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i,v in enumerate(a): r[i] += v
    for i,v in enumerate(b): r[i] += v
    return ptrim(r)
def psub(a, b): return padd(a, [-v for v in b])
def pscale(a, c): return [c*v for v in a]
def qpow(m): return [((-1)**k)*comb(m,k) for k in range(m+1)]
def Em(m): return ptrim([-1] + [1]*m) if m >= 0 else ([-1] if m == -1 else [])
def bern_to_poly(cells, d):
    U = [0]*(d+1); U[0] = cells[d]
    for k in range(d-1, -1, -1):
        prev = 0
        for i in range(d+1):
            cur = U[i]; U[i] = cur - prev; prev = cur
        U[d-k] += cells[k]
    return ptrim(U)
def ballot(d):
    return [comb(d-1,k) - (comb(d-1,k-1) if k >= 1 else 0) for k in range(d+1)]
def Gpoly(R):
    tg = psub(qpow(R-1), Em(R-1))
    assert all(v % 2 == 0 for v in tg)
    return [v//2 for v in tg]

ok_all = True
for fn, R in (("amm12592_floor_witness_R16_ruleB_boxeph.json", 16),
              ("amm12592_witness_R128_ruleB_D0_0_boxeph.json", 128)):
    d = json.load(open(os.path.join(CP, fn)))
    prof, blocks = d["profile"], d["blocks"]
    # Delta residuals
    sig = qpow(R-1); sigs = []
    for i in range(R):
        tpoly = psub(sig, bern_to_poly(blocks[i], prof[i]))
        assert (not tpoly) or tpoly[0] == 0
        sig = tpoly[1:] if tpoly else []
        sigs.append(ptrim(sig))
    # y tails
    T = Gpoly(R); Ts = []
    for i in range(R-1):
        b = ballot(prof[i])
        y = [(blocks[i][k] - b[k])//2 for k in range(prof[i]+1)]
        g = bern_to_poly(y, prof[i])
        tpoly = psub(T, g)
        assert (not tpoly) or tpoly[0] == 0
        T = tpoly[1:] if tpoly else []
        Ts.append(ptrim(T))          # this is T_{i+1}
    ok = True
    for i in range(R-1):             # i = 0..R-2: sigma_i vs E_{R-2-i} + 2T_{i+1}
        lhs = sigs[i]
        rhs = padd(Em(R-2-i), pscale(Ts[i], 2))
        if lhs != rhs: ok = False; break
    ok_all = ok_all and ok
    log("  [%s] sigma_i = E_{R-2-i} + 2 T_{i+1} on %s (i = 0..%d)"
        % ("PASS" if ok else "FAIL", fn, R-2))
log("  => Rule A attractor (sigma ~ E-shape) == y-tail smallness, exactly.")
flush()
