#!/usr/bin/env python3
"""
lrc14_Dm_cell_sharpening_klein_S214.py

Sharpen boxeph-S3's D_m gap (measured <=4% vs needed 5.7%): the CELL structure.

KEY STRUCTURE: within a maximal run of j where the widest gap is bounded by the
same teeth pair (e_L, e_R) and no wrap occurs, boxeph's canonical phase is
LINEAR in j:
    phi_j = alpha + beta j,  beta = [e_L/(1-r) + e_R/(1+r)]/(2V),  r = s/V,
so the killer position m*tau_j = m(1+beta)j/V + m*alpha/V is a per-cell
ARITHMETIC PROGRESSION: the kill comb is per-cell Denjoy-Koksma territory.

MEASURE (on boxeph-bank-style instances, exact):
 [1] #cells over Good; cell-length distribution; per-cell linearity check.
 [2] per-killer, per-cell: e_cell = kill_cell - len_cell/7 (SIGNED);
     D_m = |sum e_cell|; Sum|e_cell| (the absolute per-cell route);
     the boundary-cancellation ratio D_m / Sum|e_cell|.
 [3] verdict per killer: D_m vs 0.057|Good| (needed) and what a provable
     per-cell bound requires: if D_m ~ 2#cells the route is vacuous unless
     #cells small; if signs cancel across cells (ratio << 1), the provable
     object is the TELESCOPED comb (continuity of m*tau_j within same-gap
     runs) -- the sharpening.
"""

from fractions import Fraction as F
from math import gcd
import itertools as it

ONE14 = F(1, 14)


def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def is_primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return g == 1


def good_periods_with_pairs(E, V):
    """boxeph's pipeline + the widest-gap PAIR identity per j.
    Returns {j: (tau, (eL, eR))}."""
    s = max(E)
    out = {}
    for j in range(V):
        pos = sorted(((e * j) % V, e) for e in set(E))
        n = len(pos)
        best, ustart, pair = -1, 0, None
        for i in range(n):
            a, eL = pos[i]
            b, eR = (pos[(i + 1) % n][0] + (V if i == n - 1 else 0),
                     pos[(i + 1) % n][1])
            if b - a > best:
                best, ustart, pair = b - a, a, (eL, eR)
        if 7 * best <= V:
            continue
        a, g, r = F(ustart, V), F(best, V), F(s, V)
        lo = (a + ONE14) / (1 - r)
        hi = (a + g - ONE14) / (1 + r)
        if lo >= hi:
            continue
        phi = (lo + hi) / 2
        tau = (j + phi) / V
        ok = True
        for e in set(E):
            v = V - e
            x = (v * tau) % 1
            if min(x, 1 - x) < ONE14:
                ok = False
                break
        if ok:
            out[j] = (tau, pair)
    return out


def analyze(name, killers, E, V):
    S = sorted(killers + [V - e for e in set(E)])
    good = good_periods_with_pairs(E, V)
    js = sorted(good)
    nG = len(js)
    print(f"\n=== {name}: V={V}, |Good|={nG}, covering={is_covering(S)}, "
          f"primitive={is_primitive(S)} ===")
    if nG == 0:
        return
    # [1] cells: maximal runs of consecutive-in-Good j with same pair AND
    # linear phi continuation (detect slope breaks exactly)
    cells = []
    cur = [js[0]]
    for prev, j in zip(js, js[1:]):
        same_pair = good[j][1] == good[prev][1]
        if same_pair and len(cur) >= 2:
            # check linearity: phi_j - phi_prev == previous step?
            step_prev = good[cur[-1]][0] * V - cur[-1] - (good[cur[-2]][0] * V - cur[-2])
            step_now = (good[j][0] * V - j) - (good[prev][0] * V - prev)
            lin = (j - prev == cur[-1] - cur[-2]) and step_now == step_prev
        else:
            lin = same_pair
        if same_pair and lin:
            cur.append(j)
        else:
            cells.append(cur)
            cur = [j]
    cells.append(cur)
    lens = sorted(len(c) for c in cells)
    print(f"  [1] #cells={len(cells)}; len dist: min={lens[0]}, med={lens[len(lens)//2]}, "
          f"max={lens[-1]}; #len1={sum(1 for c in cells if len(c)==1)}")

    need = 0.0571 * nG
    print(f"  [2] per-killer cell ledger (needed D_m <= {need:.1f}):")
    for m in killers:
        tot_kill = 0
        e_cells = []
        for cell in cells:
            kill = 0
            for j in cell:
                x = (m * good[j][0]) % 1
                if min(x, 1 - x) < ONE14:
                    kill += 1
            tot_kill += kill
            e_cells.append(kill - len(cell) / 7)
        D = abs(tot_kill - nG / 7)
        sabs = sum(abs(e) for e in e_cells)
        zone = 'P' if m <= 13 else ('T' if m <= V / 14 else ('M' if m < 9 * V / 14 else 'L'))
        print(f"    m={m:5d} ({zone}, gcd={gcd(m, V):3d}): kill={tot_kill:4d} "
              f"D_m={D:6.1f} ({D/nG*100:4.1f}%)  Sum|e_cell|={sabs:6.1f} "
              f"cancel={D/sabs if sabs else 0:5.2f}  "
              f"{'NEEDED-OK' if D <= need else 'over'}")


def find_covering_variant(base_killers, E, V, scan):
    for deltas in it.product(range(0, 8), repeat=len(base_killers)):
        K = [k + d for k, d in zip(base_killers, deltas)]
        if len(set(K)) < len(K):
            continue
        S = sorted(K + [V - e for e in set(E)])
        if len(set(S)) == 13 and is_covering(S) and is_primitive(S):
            return K
    return None


def main():
    print("=" * 78)
    print("D_m cell sharpening (klein-S214), on boxeph-S3-style banks")
    print("=" * 78)

    # bank-2 style: V=842, k~=10 cluster, 3 mid-band killers
    V = 842
    E = [0, 11, 37, 68, 105, 133, 160, 191, 224, 260]
    K = find_covering_variant([V // 2 - 1, V // 3 + 1, V // 5 + 2], E, V, 8)
    if K:
        analyze("bank2-style 3M", K, E, V)

    # bank-4 style: k~=8, V=1006, 2P+3M
    V = 1006
    E8 = [0, 17, 55, 96, 141, 190, 243, 300]
    K = None
    for cand in ([12, 13, 251, 402, 503], [12, 13, 253, 405, 505]):
        K2 = find_covering_variant(cand, E8, V, 8)
        if K2:
            K = K2
            break
    if K:
        analyze("bank4-style 2P+3M", K, E8, V)

    print("""
READING: if cancel << 1 (cross-cell cancellation) and D_m <= need with room,
the provable route is the TELESCOPED per-cell DK (boundary terms cancel along
same-gap runs); if #len1 cells dominate, the cell route needs the coarser
'continuity of m*tau_j' argument instead. Either way the ledger shows exactly
what a proof must control.""")
    print("DONE.")


if __name__ == '__main__':
    main()
