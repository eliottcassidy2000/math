"""Rule C: hybrid bottom-exact / Bernstein-clamp attractor rule.

row i: ideal := sigma_{i-1} - x*E_{R-2-i};  T := trunc_{d_i}(ideal)
  1. bottom pass (cells s = d..d-G): coefficient-exact greedy w/ box+parity clamp;
  2. rest := T - partial;  cells k < d-G := Bernstein clamp of rest's cells.
G = guard depth (gfun(d, R, i)).
"""
import json, time
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

GS = (5979874356654402, 10**16)
def prof(R): return [GS[0]*(R+i)//GS[1] for i in range(R)]
def Em(m): return [-1] + [1]*m if m >= 0 else []

def clampfix(v, cap):
    w = min(cap, max(-cap, v))
    if (w - cap) % 2:
        w = w - 1 if w > 0 else w + 1
    return w

def emit_hybrid(ideal, d, G, qpows):
    tgt = ideal[:d+1] + [0]*max(0, d+1-len(ideal))
    delta = [0]*(d+1)
    partial = [0]*(d+1)
    G = min(G, d)
    for j in range(G+1):
        s = d - j
        want = tgt[j] - partial[j]
        w = clampfix(want, comb(d, s))
        delta[s] = w
        if w:
            qs = qpows[s]; base = d - s
            for t in range(min(len(qs), d - base + 1)):
                partial[base + t] += w * qs[t]
    rest = [tgt[j] - partial[j] for j in range(d+1)]
    cr = poly_to_bern(rest, d)
    for k in range(d - G):
        delta[k] = clampfix(cr[k], comb(d, k))
    return delta

def run_ruleC(R, gfun, profile=None, qpows=None):
    pr = profile if profile is not None else prof(R)
    if qpows is None: qpows = [qpow(s) for s in range(pr[-1]+1)]
    sig = qpow(R-1)
    blocks = []
    for i in range(R):
        d = pr[i]
        ideal = psub(sig, pshift(Em(R-2-i), 1))
        delta = emit_hybrid(ideal, d, gfun(d, R, i), qpows)
        D = bern_to_poly(delta, d)
        t = psub(sig, D)
        if t and t[0] != 0:
            return None, f"die row {i}: const |.|~2^{abs(t[0]).bit_length()}"
        sig = t[1:] if t else []
        blocks.append(delta)
    if sig == []:
        ok = all(admissible(blocks[i], pr[i]) for i in range(R)) and \
             epoch_sum(R, blocks, pr) == qpow(R-1)
        return blocks, f"CLOSED verified={ok}"
    return None, f"final residual L1~2^{sum(abs(v) for v in sig).bit_length()}"

if __name__ == "__main__":
    GFUNS = {
        "G=d/2":   lambda d, R, i: d//2,
        "G=d/3":   lambda d, R, i: d//3,
        "G=2d/3":  lambda d, R, i: (2*d)//3,
        "G=d-R/4": lambda d, R, i: max(0, d - R//4),
        "G=i+1":   lambda d, R, i: i+1,
    }
    for name, gf in GFUNS.items():
        line = [name]
        for R in (32, 64, 128, 256):
            t0 = time.time()
            sol, msg = run_ruleC(R, gf)
            line.append(f"R{R}:{'CLOSED' if sol else msg}")
        print("  ".join(line), flush=True)
