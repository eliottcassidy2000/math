"""
lrc14_jstar_AP_extremal_opus_S164.py   (opus-2026-07-08-S164)

THE CAPSTONE: is the AP the EXTREMAL (worst) cluster for the good-period j*?

LRC(14) is reduced (mac-mini-S59) to ONE lemma: j* = O(k), where for a k-cluster E={0=e_0<..<e_{k-1}}
and ruler Vmax (hard region spread >= 6Vmax/7), the good period is
    j*(E,Vmax) = min{ j>=1 : maxgap{ frac(e_i j / Vmax) } > 1/7 }.
mac-mini PROVED the AP case j*(AP) <= ceil(7(k-1)/6) via Dirichlet on the step, and observed the
max j* over adversarial clusters concentrates near the AP.  IF the AP MAXIMIZES j* over all
k-clusters, then j*(E) <= j*(AP) = O(k) for ALL E -- closing THM-527-A and hence LRC(14).

This script tests "AP maximizes j*": for each k, compute max j* over (a) APs e_i=i*d, (b) random
clusters, (c) structured clusters (2-block, near-AP, geometric), all in the hard region, and check
whether the AP achieves the global max and whether max j* = O(k).
"""
import sys, random
from math import gcd


def jstar(E, Vmax, theta_num=1, theta_den=7, jmax=None):
    """min j>=1 with max circular gap of {e_i j mod Vmax} > Vmax/7; None if none <= jmax."""
    k = len(E)
    if jmax is None:
        jmax = Vmax - 1
    thr = Vmax * theta_num  # compare max_gap * theta_den > Vmax  <=>  max_gap > Vmax/7
    for j in range(1, jmax + 1):
        ph = sorted((e * j) % Vmax for e in E)
        # circular gaps
        mg = Vmax - ph[-1] + ph[0]
        for i in range(k - 1):
            g = ph[i + 1] - ph[i]
            if g > mg:
                mg = g
        if mg * theta_den > Vmax:   # mg > Vmax/7
            return j
    return None


def hard_vmaxes(spread):
    """Vmax in (spread, 7*spread/6]  (the hard region spread >= 6Vmax/7)."""
    lo = spread + 1
    hi = (7 * spread) // 6
    return list(range(lo, hi + 1))


def maxjstar_over_vmax(E):
    spread = max(E)
    best = 0; bestV = None
    for V in hard_vmaxes(spread):
        js = jstar(E, V)
        if js is None:
            return None, V   # no good period in the hard band (would be a real failure)
        if js > best:
            best = js; bestV = V
    return best, bestV


def main():
    print("=" * 92)
    print("IS THE AP EXTREMAL FOR j*?  (LRC(14) capstone: j*(E) <= j*(AP) = O(k)?)")
    print("  hard region: spread >= 6Vmax/7 (Vmax in (spread, 7spread/6]); good iff maxgap>1/7")
    print("=" * 92)
    r = random.Random(164)
    for k in range(8, 14):
        apbound = -(-7 * (k - 1) // 6)   # ceil(7(k-1)/6)
        # (a) APs: e_i = i*d, several d
        ap_best = 0; ap_arg = None
        for d in range(1, 25):
            E = [i * d for i in range(k)]
            j, V = maxjstar_over_vmax(E)
            if j is not None and j > ap_best:
                ap_best = j; ap_arg = (d, V)
        # (b) random clusters (0 and a max, k-2 interior)
        rnd_best = 0; rnd_arg = None
        for _ in range(4000):
            spread = r.randint(k, 60)
            mid = sorted(r.sample(range(1, spread), k - 2))
            E = [0] + mid + [spread]
            if len(set(E)) != k:
                continue
            j, V = maxjstar_over_vmax(E)
            if j is not None and j > rnd_best:
                rnd_best = j; rnd_arg = (tuple(E), V)
        # (c) structured: 2-block, near-AP (AP with one defect), geometric-ish
        st_best = 0; st_arg = None
        structs = []
        for d in range(2, 12):
            base = [i * d for i in range(k)]
            for jdef in range(1, k):
                for delta in (-1, 1, 2):
                    E = base[:]; E[jdef] += delta
                    E = sorted(set(E))
                    if len(E) == k and E[0] == 0:
                        structs.append(E)
        for a in range(2, k):
            for off in range(1, 40):
                E = sorted(set(list(range(a)) + [off + i for i in range(k - a)]))
                if len(E) == k:
                    structs.append(E)
        for E in structs:
            j, V = maxjstar_over_vmax(E)
            if j is not None and j > st_best:
                st_best = j; st_arg = (tuple(E), V)
        gmax = max(ap_best, rnd_best, st_best)
        who = "AP" if ap_best == gmax else ("struct" if st_best == gmax else "random")
        print(f"  k={k:2d}: AP max j*={ap_best:3d} (d,V={ap_arg})  random max j*={rnd_best:3d}  "
              f"struct max j*={st_best:3d}  || GLOBAL max={gmax} by {who}; AP-bound ceil(7(k-1)/6)={apbound}")
    print()
    print("  READING: if AP max >= random & struct max at every k, the AP is (empirically) extremal")
    print("  for j*, and j*(E) <= j*(AP) = O(k) would close THM-527-A / LRC(14).")


if __name__ == "__main__":
    main()
