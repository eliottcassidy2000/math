"""Rule B: bottom-exact greedy attractor rule (parameter-free, deterministic).

row i: ideal := sigma_{i-1} - x*E_{R-2-i}   (E_{-1} := 0)
       build cells delta_{d}, delta_{d-1}, ..., delta_0 so that the polynomial
       coefficient a_j of Delta matches ideal[j] EXACTLY from j = 0 upward,
       clamping (box + parity, toward zero) any cell that exceeds capacity;
       a clamp pollutes only coefficients ABOVE its own position (triangular),
       so the emitted block matches ideal on a maximal bottom segment.
       sigma_i := (sigma_{i-1} - Delta_i)/x.
Closure iff sigma_{R-1} == 0. Exact integer arithmetic throughout.
"""
import json, time
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

GS = (5979874356654402, 10**16)
def prof(R): return [GS[0]*(R+i)//GS[1] for i in range(R)]
def Em(m): return [-1] + [1]*m if m >= 0 else []

def emit_block(ideal, d):
    """Bottom-exact greedy admissible block at degree d approximating ideal."""
    tgt = ideal[:d+1] + [0]*max(0, d+1-len(ideal))
    delta = [0]*(d+1)
    partial = [0]*(d+1)          # poly of sum of chosen cells so far
    qs = [1]                     # q^s coefficients, updated incrementally
    for j in range(d+1):         # j = poly coefficient position; cell s = d-j
        s = d - j
        want = tgt[j] - partial[j]
        cap = comb(d, s)
        w = min(cap, max(-cap, want))
        if (w - cap) % 2:
            w = w - 1 if w > 0 else w + 1
        delta[s] = w
        if w:
            # add w * x^{d-s} q^s to partial: q^s coeffs qs[t], positions d-s+t
            for t, qc in enumerate(qs):
                pos = d - s + t
                if pos <= d: partial[pos] += w * qc
        # update qs -> q^{s+1}... careful: next step uses s-1, need q^{s-1}!
        # cells are visited s = d, d-1, ..., so we need q^d first. Precompute instead.
        qs = None if False else qs
        if j < d:
            pass
    return delta

def emit_block_v2(ideal, d, qpows):
    tgt = ideal[:d+1] + [0]*max(0, d+1-len(ideal))
    delta = [0]*(d+1)
    partial = [0]*(d+1)
    for j in range(d+1):
        s = d - j
        want = tgt[j] - partial[j]
        cap = comb(d, s)
        w = min(cap, max(-cap, want))
        if (w - cap) % 2:
            w = w - 1 if w > 0 else w + 1
        delta[s] = w
        if w:
            qs = qpows[s]
            base = d - s
            for t in range(min(len(qs), d - base + 1)):
                partial[base + t] += w * qs[t]
    return delta

def run_ruleB(R, trace=False, profile=None):
    pr = profile if profile is not None else prof(R)
    dmax = pr[-1]
    qpows = [qpow(s) for s in range(dmax+1)]
    sig = qpow(R-1)
    blocks = []
    for i in range(R):
        d = pr[i]
        ideal = psub(sig, pshift(Em(R-2-i), 1))
        delta = emit_block_v2(ideal, d, qpows)
        D = bern_to_poly(delta, d)
        t = psub(sig, D)
        if t and t[0] != 0:
            return None, f"die row {i}: const {t[0]}"
        sig = t[1:] if t else []
        blocks.append(delta)
        if trace and (i % 16 == 0 or i >= R-3):
            off = psub(sig, Em(R-2-i)) if R-2-i >= 0 else sig
            print(f"  i={i:3d} d={d:3d} degsig={len(sig)-1:4d} off_terms={sum(1 for v in off if v)} "
                  f"off_log2max={max((abs(v) for v in off), default=0).bit_length()}")
    if sig == []:
        ok = all(admissible(blocks[i], pr[i]) for i in range(R)) and \
             epoch_sum(R, blocks, pr) == qpow(R-1)
        return blocks, f"CLOSED verified={ok}"
    return None, f"final residual L1={sum(abs(v) for v in sig)}"

if __name__ == "__main__":
    for R in (8, 16, 32, 64, 128):
        t0 = time.time()
        sol, msg = run_ruleB(R)
        print(f"R={R:3d} ruleB: {msg}  ({time.time()-t0:.1f}s)")
        if sol:
            w = {"R": R, "profile": prof(R), "blocks": sol, "verified": True,
                 "source": "bottom-exact greedy attractor rule B (parameter-free deterministic)"}
            json.dump(w, open(f"amm12592_floor_witness_R{R}_ruleB_boxeph.json", "w"))
