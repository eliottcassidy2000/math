#!/usr/bin/env python3
r"""
lrc_64_classes_tight_realizable_s577b.py    oracle-2026-06-03-S577o

Refuted the naive model (s577): the 64 self-converse classes are NOT the 2^6 antipodal
tie-resolutions of the regular 14-gon (those give 33 round classes, 15 self-converse).
So the worry-set is SPREAD across the self-converse boundary, not concentrated at the AP.

THE REFRAMED QUESTION (census-backed, S576: only the AP is tight at n=14). Of the 64
self-converse round classes, which are realized at TIGHTNESS (M=1/n) by an actual speed
set, vs only by STRICTLY lonely configs (M>1/n)? If only the REGULAR class is
tight-realizable, the worry collapses 64 -> 1 (and the AP is lonely).

We work from the SPEED side: for a curated + sampled set of configs we compute M(S) exactly
(pinch) and, at the optimal time t* (perturbed off the tie-boundary by +/- eps to get an
OPEN time), the half-turn runner tournament's canonical class. We tabulate, per round class:
min M over realizations, whether self-converse, and which classes ever reach M=1/14.
"""
import importlib.util
from fractions import Fraction as Fr
from functools import reduce
from math import gcd
from pathlib import Path
import random

def load_s574():
    spec = importlib.util.spec_from_file_location(
        "s574", Path("04-computation/lrc_round_count_m89_s574.py").resolve())
    M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
    return M
MOD = load_s574()
canon = MOD.canon

N = 14
TH = Fr(1, N)

def circ(r, C):
    r %= C
    return min(r, C - r)

def pair_sums(S):
    return sorted({a + b for i, a in enumerate(S) for b in S[i + 1:]})

def M_exact(S):
    best = Fr(0); bt = None
    for C in pair_sums(S):
        for m in range(1, C):
            md = min(circ(v * m, C) for v in S)
            val = Fr(md, C)
            if val > best:
                best, bt = val, (m, C)
    return best, bt

def half_turn_class(S, m, C, sign):
    """tournament class at open time t = m/C + sign*eps (eps -> 0+). Positions v_i*t mod 1;
    order by exact position, ties broken by sign*(v_i) (the eps-perturbation direction)."""
    # exact base position numerators p_i = (v_i*m) mod C  over C ; perturb by sign*eps*v_i
    items = []
    for v in S:
        p = (v * m) % C
        items.append((p, sign * v, v))   # sort key: (base position, eps-direction)
    # sort around circle by (p, perturbation)
    order = sorted(range(len(S)), key=lambda i: (items[i][0], items[i][1]))
    k = len(S)
    # positions on circle in sorted order; half-turn tournament:
    # runner a beats b if (pos_b - pos_a) in (0, 1/2) going forward.
    # Use the sorted ranks with the C-grid + eps to decide; compute via fractional pos.
    # fractional position f_i = p_i/C + sign*eps*v_i ; compare circular distance to 1/2.
    adj = [[0] * k for _ in range(k)]
    for ia in range(k):
        for ib in range(ia + 1, k):
            pa, ea, _ = items[order[ia]]; pb, eb, _ = items[order[ib]]
            # circular forward gap from a to b (a earlier in order): (pb-pa) mod C, tie via eps
            diff = (pb - pa) % C
            half = Fr(C, 2)
            if diff < half:
                a_beats_b = True
            elif diff > half:
                a_beats_b = False
            else:                            # antipodal tie -> resolved by eps direction
                a_beats_b = (ea < eb)        # smaller perturbation is "behind"
            ra, rb = order[ia], order[ib]
            if a_beats_b: adj[ra][rb] = 1
            else: adj[rb][ra] = 1
    return canon(adj, k)

def opp_canon(c, k):
    # we only stored canon; recompute self-converse by building? store adj instead
    return None

def main():
    print("=" * 84)
    print("LRC n=14: which of the 64 self-converse classes are TIGHT-realizable (M=1/14)?")
    print("=" * 84)

    # self-converse set for reference
    reps = {}
    for d in MOD.valid_dvectors(N - 1):
        a = MOD.build_adj(d, N - 1); reps.setdefault(canon(a, N - 1), a)
    def opp(a):
        m = len(a); return [[0 if i==j else a[j][i] for j in range(m)] for i in range(m)]
    sc = {c for c, a in reps.items() if c == canon(opp(a), N - 1)}
    print(f"  round classes={len(reps)}, self-converse={len(sc)}")

    # curated tight / near-tight + random
    AP = tuple(range(1, 14))
    Vstar = tuple(sorted(list(range(1, 12)) + [13, 24]))
    configs = {AP: "AP", Vstar: "V*"}
    rnd = random.Random(577)
    # AP single-coordinate perturbations (near-tight hunters) + random
    for d in range(14, 40):
        c = tuple(sorted(list(range(1, 13)) + [d]))
        if reduce(gcd, c) == 1 and len(set(c)) == 13: configs[c] = f"AP12+{d}"
    while len([k for k in configs]) < 400:
        c = tuple(sorted(rnd.sample(range(1, 60), 13)))
        if reduce(gcd, c) == 1: configs[c] = "rand"

    per_class_minM = {}     # class -> (minM, example, tag)
    tight_classes = set()
    ap_classes = set()
    for S, tag in configs.items():
        M, (m, C) = M_exact(S)
        for sign in (+1, -1):
            try:
                cls = half_turn_class(S, m, C, sign)
            except Exception:
                continue
            cur = per_class_minM.get(cls)
            if cur is None or M < cur[0]:
                per_class_minM[cls] = (M, S, tag)
            if M == TH:
                tight_classes.add(cls)
                if tag == "AP": ap_classes.add(cls)

    nclasses = len(per_class_minM)
    n_sc_seen = sum(1 for c in per_class_minM if c in sc)
    print(f"\n  distinct optimal-time classes seen: {nclasses}  (self-converse among them: {n_sc_seen})")
    print(f"  classes that are TIGHT-realizable (some realization has M=1/14): {len(tight_classes)}")
    print(f"  the AP's optimal-time class(es): {len(ap_classes)}  self-converse? "
          f"{all(c in sc for c in ap_classes)}")
    # are all tight classes self-converse? are they the AP class?
    print(f"  all tight classes self-converse? {all(c in sc for c in tight_classes)}")
    print(f"  tight classes == AP classes? {tight_classes == ap_classes}")

    # min M over non-AP-class realizations
    others = [(M, tag, S) for c, (M, S, tag) in per_class_minM.items() if c not in ap_classes]
    if others:
        mo = min(others)
        print(f"\n  min M among classes OTHER than the AP's: {mo[0]} = {float(mo[0]):.4f} "
              f"(> 1/14={float(TH):.4f}? {mo[0] > TH})  via {mo[2]} [{mo[1]}]")
    print(f"\n  tight (M=1/14) examples by class:")
    shown = 0
    for c, (M, S, tag) in per_class_minM.items():
        if M == TH and shown < 6:
            print(f"     class {str(c)[:34]}...  M=1/14  via {tag} {S}  self-conv={c in sc}")
            shown += 1

    print("\n" + "=" * 84)
    print("READING")
    print("=" * 84)
    print("""  If the only TIGHT-realizable optimal-time class is the AP's (regular) one, the worry
  collapses 64 -> 1: every other self-converse class is realized only by configs with
  M>1/14 (strictly lonely). The proof target becomes: (a) the regular class is lonely
  (AP, t=1/14, proven); (b) every NON-regular self-converse class forces M>1/14 -- a
  realization-independent strict-loneliness statement. The data here tests (b) on a
  bounded sample; the gap to a proof is making it realization-independent + unbounded.""")

if __name__ == "__main__":
    main()
