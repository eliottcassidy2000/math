"""
lrc14_Vstar_crossover_across_families_kps.py  (kind-pasteur-2026-06-27-S31ag)

Answers mac-mini-S61 handoff A: confirm the V* crossover (= the unified constant
D ~ V* ~ paper-wall) is worst-case across covering families, not just {1..12,14V}.

For each covering family parameterized by apex V, compute the LARGEST direct lonely
arc L_arc(V) = max interval of {x: ||s_i x|| >= 1/14 for all i}. It is bounded below
(~l_core) for V <= V*, then decays ~1/V. Find V* = crossover where L_arc drops below
the bounded-core floor, and confirm V* and the floor across families. Pins D = V*.
"""
import sys

def nrm(x):
    f = x - int(x)
    if f < 0: f += 1.0
    return min(f, 1 - f)

def largest_lonely_arc(S, thr=1/14):
    """float: largest interval of {x in [0,1): ||s x||>=thr all s in S}."""
    S = [s for s in S if s != 0]
    bps = [0.0, 1.0]
    for s in S:
        inv = 1.0 / s
        # lonely arcs of s: [(k+thr)/s, (k+1-thr)/s]; collect endpoints
        k = 0
        while k <= s:
            lo = (k + thr) * inv
            hi = (k + 1 - thr) * inv
            if 0 <= lo <= 1: bps.append(lo)
            if 0 <= hi <= 1: bps.append(hi)
            k += 1
    bps = sorted(set(round(b, 12) for b in bps))
    best = 0.0; cur0 = None
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        mid = 0.5 * (a + b)
        lonely = True
        for s in S:
            if nrm(s * mid) < thr - 1e-12:
                lonely = False; break
        if lonely:
            if cur0 is None: cur0 = a
            cure = b
        else:
            if cur0 is not None:
                best = max(best, cure - cur0); cur0 = None
    if cur0 is not None: best = max(best, cure - cur0)
    return best

FAMILIES = {
    "aliasing {1..12,14V}":      lambda V: list(range(1, 13)) + [14 * V],
    "top-bal {1..11,13,14V}":    lambda V: list(range(1, 12)) + [13, 14 * V],
    "wide-doublet{1..11,14V,+1}":lambda V: list(range(1, 12)) + [14 * V, 14 * V + 1],
    "cap8-core {1,5,7,8,9,..,14V}": lambda V: [1,5,7,8,9,10,11,12,13,2,3,4] + [14*V],
}

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    FLOOR = 0.0045   # bounded-core long-arc floor (~ THM-575/V* atlas)
    print(f"V* crossover (L_arc drops below floor {FLOOR}) across covering families:")
    print(f"{'family':30} {'V*':>6} {'L_arc@2':>9} {'L_arc@V*':>9} {'L_arc@2V*':>10}")
    for name, fam in FAMILIES.items():
        Vstar = None; larc2 = None; larcStar = None; larc2Star = None
        prev_above = True
        for V in range(2, 360):
            S = fam(V)
            if len(set(s % 1 for s in S)) and len(S) != 13:
                # keep families at 13 speeds
                S = sorted(set(S))
                if len(S) != 13: continue
            L = largest_lonely_arc(sorted(set(S)))
            if V == 2: larc2 = L
            if L < FLOOR and Vstar is None and V > 14:
                Vstar = V; larcStar = L
            if Vstar is not None and V == min(2 * Vstar, 359):
                larc2Star = L
                break
        print(f"{name:30} {str(Vstar):>6} {larc2:>9.5f} "
              f"{(larcStar if larcStar else 0):>9.5f} {(larc2Star if larc2Star else 0):>10.5f}")
    print()
    print("Expected (mac-mini-S61): V* ~ 200-234 across families = D = paper wall = 1/(14*l_core).")
    print("Bounded-core floor l_core ~ 0.005-0.006; V* worst-case confirms one unified constant.")
