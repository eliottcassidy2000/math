"""Parameter-free deterministic construction rule for the gamma* floor epoch:

  row i:  ideal_i := sigma_{i-1} - x*E_{R-2-i}     (E_{-1} := 0)
          Delta_i := AdmClamp_{d_i}( trunc_{d_i}( ideal_i ) )
          sigma_i := (sigma_{i-1} - Delta_i)/x
Closure iff sigma_{R-1} = 0. AdmClamp: cellwise box clamp + parity fix (variant).
"""
import sys, json, os
from math import comb
from amm12592_allR_family_toolbox_boxeph import *

GS = (5979874356654402, 10**16)
def prof(R): return [GS[0]*(R+i)//GS[1] for i in range(R)]
def Em(m): return [-1] + [1]*m if m >= 0 else []

def adm_clamp(P, d, pf='tozero'):
    """P poly deg<=d -> admissible cells at degree d."""
    cells = poly_to_bern(P + [0]*max(0, 0), d) if len(P) <= d+1 else None
    out = []
    for k, v in enumerate(cells):
        cap = comb(d, k)
        w = min(cap, max(-cap, v))
        if (w - cap) % 2:
            if pf == 'tozero':
                w = w - 1 if w > 0 else w + 1
            elif pf == 'toward':
                if v > w and w + 1 <= cap: w += 1
                elif v < w and w - 1 >= -cap: w -= 1
                else: w = w - 1 if w > 0 else w + 1
            elif pf == 'away':
                if w + 1 <= cap and w >= 0: w += 1
                elif w - 1 >= -cap: w -= 1
        out.append(w)
    return out

def run_rule(R, pf='tozero', trace=False):
    pr = prof(R)
    sig = qpow(R-1)
    blocks = []
    for i in range(R):
        d = pr[i]
        m = R-2-i
        ideal = psub(sig, pshift(Em(m), 1))
        T = ideal[:d+1]
        delta = adm_clamp(T, d, pf)
        D = bern_to_poly(delta, d)
        t = psub(sig, D)
        if t and t[0] != 0:
            return None, f"die row {i}: sigma-Delta not divisible by x (const {t[0]})", blocks
        sig = t[1:] if t else []
        blocks.append(delta)
        if trace:
            L1 = sum(abs(v) for v in sig)
            off = psub(sig, Em(R-2-i)) if R-2-i >= 0 else sig
            offL1 = sum(abs(v) for v in off)
            print(f"  i={i:3d} d={d:3d} degsig={len(sig)-1:4d} L1={L1} offtrack_L1={offL1}")
    if sig == []:
        # verify fully
        ok = all(admissible(blocks[i], pr[i]) for i in range(R)) and \
             epoch_sum(R, blocks, pr) == qpow(R-1)
        return blocks, f"CLOSED verified={ok}", blocks
    return None, f"final residual nonzero L1={sum(abs(v) for v in sig)}", blocks

if __name__ == "__main__":
    for pf in ('tozero', 'toward', 'away'):
        for R in (8, 16, 32):
            sol, msg, _ = run_rule(R, pf, trace=(R==8 and pf=='tozero'))
            print(f"R={R:3d} pf={pf:7s}: {msg}")
