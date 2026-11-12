#!/usr/bin/env python3
"""Independent verifier (part D/E) - compare the exact collision census with the bounds used in the proof of (D).

Inputs: de_collisions.json (exact counts per (t', s, max run)), de_gcdsum.json (S_w, per-pair lattice sums).
For every n and every (t', s) with t' >= 1, s >= 0, w = n - t' - s >= 1:
   E(t', s)   = number of ORDERED colliding pairs of n-bit words with these parameters (exact)
   L3 <= L2 <= L1 = per-pair lattice-line bounds summed over all odd window pairs (see de_gcdsum.py)
   B(t', s)   = 4 * 4^w + 4 * 2^max(t'-1, s) * S_w            (the bound written in the proof of (D))
and checks E <= L3 <= L2 <= L1 <= B.  Then the totals of (D).
"""
import sys, json, math

col = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'de_collisions.json'))
gs = json.load(open(sys.argv[2] if len(sys.argv) > 2 else 'de_gcdsum.json'))
nlist = [int(a) for a in sys.argv[3:]] if len(sys.argv) > 3 else list(range(8, 23))

S = {int(w): gs[w]['S'] for w in gs}
wavail = max(S)
print("S_w available for w <= %d" % wavail)
bad = 0
for n in nlist:
    ent = col[str(n)]
    cnt = {tuple(int(a) for a in k.split(',')): v for k, v in ent['cnt'].items()}
    Rn = math.ceil(math.log2(n)) + 1
    t = -(-n // 4)
    tot_pairs = 0
    worst = {'E/L3': (0, None), 'E/L1': (0, None), 'E/B': (0, None), 'L1/B': (0, None)}
    rows = []
    ncheck = 0
    for tp in range(1, n):
        for s in range(0, n - tp):
            w = n - tp - s
            c = cnt.get((tp, s))
            E = 2 * sum(c) if c else 0
            tot_pairs += E
            if w > wavail:
                # no S_w: the first term alone already exceeds the number of all pairs with these parameters
                allpairs = 2 ** (tp - 1) * 2 ** s * (4 ** (w - 1) if w >= 2 else 2)
                assert 4 * 4 ** w >= allpairs >= E
                continue
            cb = gs[str(w)]['combos']['%d,%d' % (tp, s)]
            L1, L2, L3 = cb['L1'], cb['L2'], cb['L3']
            B = 4 * 4 ** w + 4 * 2 ** max(tp - 1, s) * S[w]
            ok = (E <= L3 <= L2) and (L2 <= L1 * (1 + 1e-12)) and (L1 <= B * (1 + 1e-12))
            ncheck += 1
            if not ok:
                bad += 1
                print("VIOLATION n=%d t'=%d s=%d w=%d E=%d L3=%d L2=%d L1=%.1f B=%.1f" % (n, tp, s, w, E, L3, L2, L1, B))
            for key, val in (('E/L3', E / L3), ('E/L1', E / L1), ('E/B', E / B), ('L1/B', L1 / B)):
                if val > worst[key][0]:
                    worst[key] = (val, (tp, s, w))
            rows.append((tp, s, w, E, L3, L2, L1, B))
    assert tot_pairs == 2 * ent['unordered_pairs']
    print("\nn=%d  (R=%d, t=%d): ordered colliding pairs=%d ; (t',s) cells checked=%d ; worst ratios: %s"
          % (n, Rn, t, tot_pairs, ncheck, {k: ('%.4f' % v[0], v[1]) for k, v in worst.items()}))
    # run-free sets, several R
    for R in sorted(set([2, 3, 4, 5, 6, 7, 8, Rn])):
        ER = {}
        for (tp, s), c in cnt.items():
            e = 2 * sum(c[:R + 1])
            if e:
                ER[(tp, s)] = e
        mx = max((abs(s - tp) for (tp, s) in ER), default=-1)
        tot = sum(ER.values())
        totD = sum(e for (tp, s), e in ER.items() if tp >= t)
        ncell = len([1 for (tp, s) in ER if tp >= t])
        flag = 'OK' if mx <= R + 3 else 'VIOLATION'
        if mx > R + 3:
            bad += 1
        print("   R=%d%s: ordered colliding pairs in W_R: %d ; max|s-t'|=%d (<= R+3=%d %s) ; with t'>=t: %d pairs in %d cells "
              "(E/2^n = %.5f)" % (R, '*' if R == Rn else ' ', tot, mx, R + 3, flag, totD, ncell, totD / 2 ** n))
    # totals of (D) with the theorem's R and t
    ED = 0
    BD = 0.0
    L1D = 0.0
    L3D = 0
    ncells = 0
    t1 = 0.0
    t2 = 0.0
    for tp in range(t, n):
        for s in range(max(0, tp - Rn - 3), min(n - 1 - tp, tp + Rn + 3) + 1):
            w = n - tp - s
            c = cnt.get((tp, s))
            ED += 2 * sum(c[:Rn + 1]) if c else 0
            ncells += 1
            if w <= wavail:
                cb = gs[str(w)]['combos']['%d,%d' % (tp, s)]
                L1D += cb['L1']
                L3D += cb['L3']
                t1 += 4 * 4 ** w
                t2 += 4 * 2 ** max(tp - 1, s) * S[w]
                BD += 4 * 4 ** w + 4 * 2 ** max(tp - 1, s) * S[w]
            else:
                BD = float('nan')
    # pairs outside the cell range (must be none)
    outside = sum(2 * sum(c[:Rn + 1]) for (tp, s), c in cnt.items() if tp >= t and abs(s - tp) > Rn + 3)
    print("   (D) totals, R=%d, t=%d: E=%d (outside |s-t'|<=R+3: %d) ; cells=%d (<= n(2R+7)=%d) ; sum L3=%d ; sum L1=%.0f ; "
          "sum B=%.0f = %.0f (4*4^w terms) + %.0f (gcd terms)" % (Rn, t, ED, outside, ncells, n * (2 * Rn + 7), L3D, L1D, BD, t1, t2))
    print("       E/2^n=%.5f  sumL3/2^n=%.2f  sumL1/2^n=%.2f  sumB/2^n=%.2f  [4^w part %.2f, gcd part %.2f]  "
          "sumB/(n^(9/2) ln n 2^n)=%.3e  n^(9/2) ln n=%.3e"
          % (ED / 2 ** n, L3D / 2 ** n, L1D / 2 ** n, BD / 2 ** n, t1 / 2 ** n, t2 / 2 ** n,
             BD / (n ** 4.5 * math.log(n) * 2 ** n), n ** 4.5 * math.log(n)))
    if n in (12, 16, 22):
        print("   table (t', s, w, E, L3, L2, L1, B) for cells with E > 0 or t' >= t and |s-t'| <= 2:")
        for (tp, s, w, E, L3, L2, L1, B) in rows:
            if E > 0 or (tp >= t and abs(s - tp) <= 2):
                print("     t'=%2d s=%2d w=%2d  E=%8d  L3=%12d  L2=%12d  L1=%14.1f  B=%16.1f  E/L3=%.4f  L1/B=%.4f"
                      % (tp, s, w, E, L3, L2, L1, B, E / L3, L1 / B))
print("\nTOTAL VIOLATIONS:", bad)
