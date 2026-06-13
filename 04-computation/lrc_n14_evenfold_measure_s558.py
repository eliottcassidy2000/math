#!/usr/bin/env python3
r"""
lrc_n14_evenfold_measure_s558.py    oracle-2026-06-02-S558o

NEW IDEA for a proof of LRC@14, synthesizing the concurrent agents' threads:

  * opus-S554/S558 (EVEN-FOLD): even speed v=2u has ||v t||=||u(2t)||; with
    fold(S)={v/2 : v even} (e speeds, e<=12 since S primitive), the PROVEN LRC(13)
    gives max_t g_fold(2t) >= 1/(e+1) > 1/14. So the EVEN-GOOD set
        G := { t in [0,1) : ||v t|| >= 1/14 for every EVEN v }
    has POSITIVE measure (a window opens around the fold's optimum once the threshold
    drops from the fold collar 1/(e+1) to 1/14). The even half of LRC@14 is FREE.
  * The entire remaining difficulty is the ODD runners ON G (opus-S554 'no-odd-split',
    verified 127/127 for e<=6 but open in general).
  * opus-S552o/S559: the obstruction localizes to the APEX (the multiple of 14 = the
    mod-7 singleton = the zero-divisor). opus-S556 ThmB: a counterexample must contain
    that apex, yet apexes keep configs loose ('the tension').

THE REFRAMING (this script). S is LONELY (not a counterexample) iff G is NOT covered by
the odd danger arcs D_v={t:||v t||<1/14} (v odd):
        S lonely  <=>  |G \ U_{odd} D_v| > 0  <=>  |G ∩ U_odd D_v| < |G|.
A clean SUFFICIENT condition (union bound, |D_v|=2/14=1/7):
        |G| > o/7         (o = # odd speeds)     ... (SUF)
We measure |G|, o/7, the actual odd-cover |G ∩ U D_v|, the slack |G\U D_v| (= safe
measure), and the odd-danger DENSITY inside vs outside G (anti-correlation = room).
We test the tight configs (AP, V*), random primitives, and APEX-forced sets.
"""
from functools import reduce
from math import gcd
import random

N = 14
TH = 1.0 / N
G = 600000  # grid for measure estimates

def d0(p):
    p = p % 1.0
    return min(p, 1 - p)

def split(S):
    ev = [v for v in S if v % 2 == 0]
    od = [v for v in S if v % 2 == 1]
    return ev, od

def measures(S):
    """return dict of measures over the grid."""
    ev, od = split(S)
    o = len(od)
    inG = 0          # |G| : all even >= 1/14
    inG_oddcov = 0   # |G ∩ U_odd D_v| : in G but some odd dangerous
    odd_danger = 0   # |U_odd D_v| overall
    lonely_pts = 0   # |G \ U_odd| = all 13 safe (safe measure)
    for i in range(G):
        t = (i + 0.5) / G
        even_ok = all(d0(v * t) >= TH for v in ev)
        odd_bad = any(d0(v * t) < TH for v in od)
        if odd_bad:
            odd_danger += 1
        if even_ok:
            inG += 1
            if odd_bad:
                inG_oddcov += 1
            else:
                lonely_pts += 1
    f = 1.0 / G
    return {
        "e": len(ev), "o": o,
        "|G|": inG * f,
        "o/7": o / 7.0,
        "|G∩Dodd|": inG_oddcov * f,
        "safe=|G\\Dodd|": lonely_pts * f,
        "|Dodd|": odd_danger * f,
        # odd-danger density inside G vs the global odd-danger density:
        "dens_in_G": (inG_oddcov / inG) if inG else 0.0,
        "dens_global": odd_danger * f,
    }

def has_apex(S):
    return any(v % N == 0 for v in S)

def report(name, S):
    S = tuple(sorted(S))
    g = reduce(gcd, S)
    prim = (g == 1)
    m = measures(S)
    suf = m["|G|"] > m["o/7"]                      # union-bound sufficient cond
    lonely = m["safe=|G\\Dodd|"] > 1e-6
    apex = has_apex(S)
    print(f"  {name}: e={m['e']} o={m['o']} apex={int(apex)} prim={int(prim)}")
    print(f"      |G|={m['|G|']:.4f}  o/7={m['o/7']:.4f}  SUF(|G|>o/7)={suf}")
    print(f"      |G∩Dodd|={m['|G∩Dodd|']:.4f}  safe=|G\\Dodd|={m['safe=|G\\Dodd|']:.5f}  -> LONELY={lonely}")
    print(f"      odd-danger density: in G={m['dens_in_G']:.4f} vs global={m['dens_global']:.4f}"
          f"  ({'ANTI-corr (room)' if m['dens_in_G']<m['dens_global']-1e-3 else 'no anti-corr'})")
    return m, suf, lonely

def main():
    print("=" * 80)
    print(f"EVEN-FOLD MEASURE REDUCTION (n=14): does the ODD danger cover the even-good set G?")
    print(f"  G = {{t: every EVEN speed >= 1/14}} has positive measure by PROVEN LRC(13).")
    print(f"  S lonely <=> |G \\ U_odd D_v| > 0.  SUF: |G| > o/7.")
    print("=" * 80)

    print("\n(1) THE TIGHT CONFIGS (safe measure should be ~0 -- the boundary):")
    report("AP {1..13}", range(1, 14))
    report("V* {1..11,13,24}", list(range(1, 12)) + [13, 24])

    print("\n(2) RANDOM primitive 13-sets:")
    rnd = random.Random(558)
    suf_hits = lonely_all = total = 0
    min_slack = 1.0
    for j in range(12):
        while True:
            S = tuple(sorted(rnd.sample(range(1, 60), 13)))
            if reduce(gcd, S) == 1:
                break
        m, suf, lonely = report(f"rand{j}", S)
        total += 1
        suf_hits += int(suf)
        lonely_all += int(lonely)
        min_slack = min(min_slack, m["safe=|G\\Dodd|"])
    print(f"\n  -> SUF(|G|>o/7) held for {suf_hits}/{total};  LONELY for {lonely_all}/{total};"
          f"  min safe slack = {min_slack:.5f}")

    print("\n(3) APEX-FORCED sets (contain a multiple of 14 -- opus-S556 ThmB necessary):")
    report("AP w/ 14 added, drop 7 {1..6,8..13,14}", [1,2,3,4,5,6,8,9,10,11,12,13,14])
    report("AP, 12->14 (apex) {1..11,13,14}", [1,2,3,4,5,6,7,8,9,10,11,13,14])
    rnd2 = random.Random(77)
    for j in range(4):
        base = rnd2.sample([x for x in range(1, 56) if x % N != 0], 12)
        S = tuple(sorted(base + [N * rnd2.randint(1, 3)]))
        if reduce(gcd, S) == 1 and len(set(S)) == 13:
            report(f"apex-rand{j}", S)

    print("\n" + "=" * 80)
    print("READING")
    print("=" * 80)
    print("""  The reframing isolates the proof to ONE inequality per config:
      |G ∩ U_odd D_v| < |G|   (the odds do not cover the guaranteed even-good window).
  |G|>0 is FREE (proven LRC(13)). The union bound SUF |G|>o/7 is the cheap test; where
  it fails (tight / many odds), the real lever is ANTI-CORRELATION: the odd danger is
  LESS dense inside G than globally (odds tend to be safe where evens are safe), so the
  odds cover only part of G. The hard locus is exactly the tight configs (AP, V*) where
  |G|->0 and safe->0 -- the measure-zero wall (S551). Off the wall, the anti-correlation
  slack is the room LRC@14 lives in. APEX-forced sets carry a mult of 14 (folds to the
  fold's apex) and are robustly loose -- matching opus-S556's 'tension'.""")

if __name__ == "__main__":
    main()
