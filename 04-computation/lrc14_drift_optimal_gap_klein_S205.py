#!/usr/bin/env python3
"""
klein-2026-07-09-S205 (CORRECTED): where does the drift-absorbed embed first fire, and WHICH gap wins?

Earlier version was flawed: it read gaps off one ruler and quoted a threshold for another (the teeth
c_i = frac(e_i*j/Vmax) depend on Vmax). Redo properly: sweep Vmax; at each Vmax test EVERY j and EVERY
tooth-free gap against the two conditions

  COARSE:  g > 1/7 + 2*spread/Vmax             (drift bounded by spread/Vmax, i.e. phi<=1)
  SHARP :  g > 1/7 + 2*spread*phi/Vmax         (phi-proportional drift -- what LRCDriftEmbed proves)

Record the least Vmax where each fires (V0), and for SHARP the WINNING gap's (a, g, phi) -- to see
whether the anchor gap (a=0, the observer's non-drifting tooth c_0=0) is the one that wins.
SOUNDNESS: the constructed tau=(j+phi)/Vmax must give minReach>=1/14 (the Lean conclusion).
"""
import numpy as np

def gaps_at(e, Vmax, j):
    c = np.unique(np.array([(ei * j) % Vmax for ei in e], dtype=float) / Vmax)
    cyc = np.append(c, 1.0)            # c[0]=0 (anchor, e_0=0); close the circle at 1
    return [(float(cyc[k]), float(cyc[k+1] - cyc[k])) for k in range(len(cyc) - 1)]

def nearInt(x):
    f = x - np.floor(x); return np.minimum(f, 1.0 - f)
def minReach(v, tau):
    return float(np.min(nearInt(np.array(v, float) * tau)))

def first_fire(e, spread, sharp: bool, Vhi):
    for Vmax in range(spread + 1, Vhi):
        v = [Vmax - ei for ei in e]
        for j in range(1, Vmax):
            for (a, g) in gaps_at(e, Vmax, j):
                phi = a + g / 2
                thr = 1/7 + (2*spread*phi/Vmax if sharp else 2*spread/Vmax)
                if g > thr:
                    tau = (j + phi) / Vmax
                    ok = minReach(v, tau) >= 1/14 - 1e-12
                    return Vmax, j, a, g, phi, ok
    return None

CLUSTERS = {
    "tightAP e={0..12}":   list(range(13)),
    "7-struct (M-128)":    [0,7,14,21,26,29,37,44,51,58,67,75,82],
    "dissociated":         [0,1,3,7,12,20,30,44,65,80,96,112,129],
    "near-AP 3*{0..10}+2": [0,3,6,9,12,15,18,21,24,27,30,31,32],
}

print("Least Vmax at which the drift-absorbed embed FIRES (some j, some free gap).")
print("hard regime is Vmax <= 7*spread/6 = 1.167*spread (j=1 fails there).\n")
print(f"{'cluster':>22} {'spread':>7} {'V0_sharp':>9} {'/spread':>8} {'V0_coarse':>10} {'/spread':>8} {'gain':>5} {'win a':>7} {'win phi':>8} {'sound':>6}")
for name, e in CLUSTERS.items():
    e = sorted(set(e)); spread = max(e)
    rs = first_fire(e, spread, True, 8*spread+2)
    rc = first_fire(e, spread, False, 12*spread+2)
    if rs is None or rc is None:
        print(f"{name:>22} {spread:>7}  (no fire in range)"); continue
    Vs, js, a_s, g_s, phi_s, ok_s = rs
    Vc = rc[0]
    print(f"{name:>22} {spread:>7} {Vs:>9} {Vs/spread:>8.3f} {Vc:>10} {Vc/spread:>8.3f} "
          f"{Vc/Vs:>5.2f} {a_s:>7.4f} {phi_s:>8.4f} {str(ok_s):>6}")

print("\nREADING: V0_sharp/spread ~ the explicit finite-check threshold. 'win a' = where the winning gap")
print("starts (a=0 would mean the ANCHOR gap wins, giving the 14x drift floor phi~1/14).")
print("If V0_sharp/spread > 1.167, the drift embed does NOT reach the hard regime -- the residual")
print("hembed corner is Vmax in (spread, V0_sharp], a BOUNDED finite check once spread is bounded.")
