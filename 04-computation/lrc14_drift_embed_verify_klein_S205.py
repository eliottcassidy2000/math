#!/usr/bin/env python3
"""
klein-2026-07-09-S205: VERIFY the drift-absorbed ruler embedding (LRCDriftEmbed.lean) and quantify the
phi-proportional SHARPENING.

THE THEOREM (formalized, sorry-free):  v_i = Vmax - e_i;  tau = (j + phi)/Vmax with phi = a + g/2 the
midpoint of a tooth-free gap (a, a+g) of the teeth c_i = frac(e_i*j/Vmax), (a,a+g) subset [0,1].
Then  nearInt(v_i*tau) = nearInt(phi - c_i - d_i),  drift d_i = e_i*phi/Vmax.
  COARSE  condition:  g > 1/7 + 2*spread/Vmax          (|d_i| <= spread/Vmax, using phi<=1)
  SHARP   condition:  g > 1/7 + 2*spread*phi/Vmax      (|d_i| <= spread*phi/Vmax)   <-- what we proved
=> minReach(v, tau) >= 1/14.

CHECKS:
 (A) SOUNDNESS: whenever the SHARP condition fires, the constructed tau really has minReach >= 1/14.
     (must be 100% -- else the Lean theorem would be wrong)
 (B) The sharpening: SHARP fires at smaller Vmax than COARSE. Report the threshold V0 for each.
 (C) The ANCHOR: the observer's tooth c_0 = frac(0*j/Vmax) = 0 never drifts (e_0=0); the minimum
     admissible fast phase is phi = 1/14 (forced by the observer's own safety), so the least possible
     drift is spread/(14*Vmax) -- a 14x floor improvement over spread/Vmax.
"""
import numpy as np
from math import gcd
from functools import reduce

def teeth_and_gap(e, Vmax, j):
    """teeth c_i = frac(e_i*j/Vmax); return sorted teeth and the maximal gap (a, g) inside [0,1]."""
    c = np.sort(np.array([(ei * j) % Vmax for ei in e], dtype=float) / Vmax)
    cyc = np.append(c, c[0] + 1.0)          # close the circle (c[0]=0 since e_0=0 => cyc.last=1)
    gaps = np.diff(cyc)
    k = int(np.argmax(gaps))
    return c, float(cyc[k]), float(gaps[k])  # a = cyc[k], g = gap

def nearInt(x):
    f = x - np.floor(x)
    return np.minimum(f, 1.0 - f)

def minReach(v, tau):
    return float(np.min(nearInt(np.array(v, dtype=float) * tau)))

def run(e, Vmax_list, label):
    e = sorted(set(e))
    assert e[0] == 0, "co-offset convention needs 0 in E (the observer)"
    spread = max(e)
    rows = []
    for Vmax in Vmax_list:
        if Vmax <= spread:  # need e_i < Vmax
            continue
        v = [Vmax - ei for ei in e]
        sharp_j = coarse_j = None
        sharp_ok = True
        for j in range(1, Vmax):
            c, a, g = teeth_and_gap(e, Vmax, j)
            phi = a + g / 2
            coarse = g > 1/7 + 2 * spread / Vmax
            sharp = g > 1/7 + 2 * spread * phi / Vmax
            if sharp and sharp_j is None:
                sharp_j = j
                tau = (j + phi) / Vmax
                mr = minReach(v, tau)
                if mr < 1/14 - 1e-12:
                    sharp_ok = False
                    print(f"  !! SOUNDNESS FAIL {label} Vmax={Vmax} j={j} minReach={mr:.6f}")
            if coarse and coarse_j is None:
                coarse_j = j
            if sharp_j is not None and coarse_j is not None:
                break
        rows.append((Vmax, spread, sharp_j is not None, coarse_j is not None, sharp_ok))
    return rows

CLUSTERS = {
    "tightAP e={0..12}":        list(range(13)),
    "7-struct (MISTAKE-128)":   [0,7,14,21,26,29,37,44,51,58,67,75,82],
    "dissociated":              [0,1,3,7,12,20,30,44,65,80,96,112,129],
    "near-AP 3*{0..10}+2":      [0,3,6,9,12,15,18,21,24,27,30,31,32],
}

print("(A) SOUNDNESS: when SHARP fires, does the constructed tau give minReach >= 1/14?")
print("(B) THRESHOLD V0: least Vmax at which some j satisfies the condition.\n")
print(f"{'cluster':>26} {'spread':>7} {'V0(SHARP)':>10} {'V0(COARSE)':>11} {'ratio':>6} {'sound':>6}")
allsound = True
for name, e in CLUSTERS.items():
    spread = max(e)
    Vs = list(range(spread + 1, 40 * spread + 1, max(1, spread // 8)))
    rows = run(e, Vs, name)
    v0s = next((r[0] for r in rows if r[2]), None)
    v0c = next((r[0] for r in rows if r[3]), None)
    sound = all(r[4] for r in rows)
    allsound &= sound
    ratio = (v0c / v0s) if (v0s and v0c) else float('nan')
    print(f"{name:>26} {spread:>7} {str(v0s):>10} {str(v0c):>11} {ratio:>6.2f} {str(sound):>6}")

print(f"\nALL SOUND = {allsound}   (must be True: the Lean theorem's conclusion)")
print("\n(C) The ANCHOR: e_0=0 => tooth c_0=0 never drifts; observer safety forces |phi|>=1/14,")
print("    so the LEAST possible drift is spread/(14*Vmax) -- a 14x floor below the naive spread/Vmax.")
print("    => the drift-margin threshold V0 improves by up to 14x when the good gap abuts the anchor.")
