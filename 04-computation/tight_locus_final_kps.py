"""
tight_locus_final_kps.py  (kind-pasteur, 2026-06-17)

FINAL consolidated census of the primitive tight locus for LRC(14).
Fast: precompute danger arcs, L-check = 'do closed arcs cover [0,1]?'.
Covers: exhaustive {1..18},{1..20}; k=1 repl to 300; k=2 repl to 60;
structurally-targeted families; max-min exact for every tight config.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import time, json

def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0: A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1: A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else: A.append((lo, hi))
    return A

_ARC = {}
def arcs(v):
    a = _ARC.get(v)
    if a is None:
        a = danger_arcs(v); _ARC[v] = a
    return a

def is_tight(S):
    """True iff closed danger arcs cover [0,1] (L=0)."""
    A = []
    for v in S: A.extend(arcs(v))
    A.sort()
    cl = ch = None; covered = Fr(0)
    for a, b in A:
        if b <= a: continue
        if ch is None: cl, ch = a, b
        elif a <= ch:
            if b > ch: ch = b
        else: covered += ch - cl; cl, ch = a, b
    if ch is not None: covered += ch - cl
    return covered == 1

def L_exact(S):
    A = []
    for v in S: A.extend(arcs(v))
    A.sort(); cl = ch = None; tot = Fr(0)
    for a, b in A:
        if b <= a: continue
        if ch is None: cl, ch = a, b
        elif a <= ch:
            if b > ch: ch = b
        else: tot += ch - cl; cl, ch = a, b
    if ch is not None: tot += ch - cl
    return 1 - tot

def gcd_list(S): return reduce(gcd, S)
def lcm2(a, b): return a*b//gcd(a, b)
def lcm_list(S): return reduce(lcm2, S)

def frac_norm(x):
    f = x - x.__floor__()
    return f if f <= Fr(1, 2) else 1 - f

def max_min_exact(S):
    S = list(S); cands = {Fr(0)}
    for v in S:
        w = Fr(1, 14*v)
        for k in range(v+1):
            c = Fr(k, v)
            for t in (c-w, c+w, c):
                if 0 <= t <= 1: cands.add(t)
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            a, b = S[i], S[j]
            for p in range(a+1):
                for q in range(b+1):
                    if a != b:
                        t = Fr(p-q, a-b)
                        if 0 <= t <= 1: cands.add(t)
                    t = Fr(p+q, a+b)
                    if 0 <= t <= 1: cands.add(t)
    best = Fr(0); arg = Fr(0)
    for t in cands:
        m = min(frac_norm(v*t) for v in S)
        if m > best: best = m; arg = t
    return best, arg

def main():
    t0 = time.time()
    tight = {}  # tuple -> source
    def consider(S, src):
        S = tuple(sorted(S))
        if len(set(S)) != 13 or any(x <= 0 for x in S): return
        if S in tight: return
        if gcd_list(S) != 1: return
        if is_tight(S):
            tight[S] = src

    AP = list(range(1, 14))

    # exhaustive {1..18}
    print("exhaustive {1..18}...", flush=True)
    for S in combinations(range(1, 19), 13):
        consider(S, "exh<=18")
    # exhaustive {1..20}
    print("exhaustive {1..20}...", flush=True)
    for S in combinations(range(1, 21), 13):
        consider(S, "exh<=20")
    # exhaustive {1..22}
    print("exhaustive {1..22}...", flush=True)
    for S in combinations(range(1, 23), 13):
        consider(S, "exh<=22")
    print(f"  after exhaustive<=22: {len(tight)} tight, t={time.time()-t0:.1f}s", flush=True)

    # k=1 replacements to 300
    print("k=1 repl to 300...", flush=True)
    for j in range(13):
        kept = [AP[i] for i in range(13) if i != j]; ks = set(kept)
        for w in range(1, 301):
            if w in ks: continue
            consider(kept+[w], "k1<=300")
    # k=2 replacements to 60
    print("k=2 repl to 60...", flush=True)
    for di, dj in combinations(range(13), 2):
        kept = [AP[i] for i in range(13) if i not in (di, dj)]; ks = set(kept)
        pool = [x for x in range(1, 61) if x not in ks]
        for a, b in combinations(pool, 2):
            consider(kept+[a, b], "k2<=60")
    print(f"  after k1,k2: {len(tight)} tight, t={time.time()-t0:.1f}s", flush=True)

    # fix {1..10}, vary top-3 to 200
    print("fix {1..10} vary top-3 to 200...", flush=True)
    base10 = list(range(1, 11))
    for trip in combinations(range(11, 201), 3):
        consider(base10+list(trip), "fix10-top3<=200")
    print(f"  after fix10: {len(tight)} tight, t={time.time()-t0:.1f}s", flush=True)

    # structurally-targeted far-outlier families (cheap)
    for W in range(14, 3001):
        consider(list(range(1, 12))+[13, W], "spor-family")
        consider(list(range(1, 13))+[W], "drop13-family")

    # report
    items = sorted(tight.items(), key=lambda kv: (max(kv[0]), lcm_list(kv[0])))
    print("\n===== PRIMITIVE TIGHT LOCUS (consolidated) =====", flush=True)
    out_list = []
    for S, src in items:
        mm, arg = max_min_exact(S)
        ce = mm < Fr(1, 14)
        print(f"  S={list(S)} lcm={lcm_list(S)} max={max(S)} max-min={mm}={float(mm):.7f} "
              f"{'<<<COUNTEREXAMPLE' if ce else '(=1/14 LRC-equality)' if mm==Fr(1,14) else ''} [{src}]",
              flush=True)
        out_list.append({"S": list(S), "lcm": lcm_list(S), "max": max(S),
                         "maxmin_num": mm.numerator, "maxmin_den": mm.denominator,
                         "is_counterexample": ce, "source": src})
    print(f"\nTotal primitive tight configs: {len(items)}", flush=True)
    if items:
        print(f"MAX ENTRY among tight: {max(max(S) for S in tight)}", flush=True)
        print(f"MAX LCM among tight:   {max(lcm_list(S) for S in tight)}", flush=True)
    ce_list = [d for d in out_list if d["is_counterexample"]]
    print(f"COUNTEREXAMPLES (max-min<1/14): {len(ce_list)}", flush=True)

    out = {"tight": out_list, "n_tight": len(items),
           "max_entry": max((max(S) for S in tight), default=None),
           "max_lcm": max((lcm_list(S) for S in tight), default=None),
           "n_counterexamples": len(ce_list),
           "elapsed_sec": time.time()-t0}
    with open("05-knowledge/results/tight_locus_final_kps.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nElapsed {time.time()-t0:.1f}s. Saved.", flush=True)

if __name__ == "__main__":
    main()
