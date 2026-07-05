#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S53 -- HYP-4109: the l=3 FLOOR-LEVEL CLOSURE at beta = 2/25.

Ledger: l=1 (S52, floor 14/169), l=2 (S52, floor 2/25, 600k domain), l>=7
(opus-S78, floor 3/19).  THIS: l=3 complete on its chain domain.

CHAIN DOMAIN at beta = 2/25 (outside closed by citation fee/window lemmas;
sorted lifted values w1 <= w2 <= w3):
  ANCHOR: w1 <= A3 = 3600/13 ~ 276.9   (3-killer fee over 9 unlifted, cite 1/10)
  MIDDLE: w2 <= (1100/51) w1 ~ 21.57x  (2-killer fee over base u {w1}, cite 1/11)
  TOP:    w3 <= 24 w2                  (1-killer window, cite 1/12 = LRC(12))

THE PYRAMID (the engineering core):
  A witness t = a/q clears a lifted coordinate c+13k by RESIDUE CLASS of k mod
  period(q) = q/gcd(13a mod q, q).  At q = 13u the orbit {x, x+13u', ...} covers
  ALL 13-lifts of a residue: the coordinate is cleared for EVERY height iff
  dist_13(a*c mod 13) >= u+1 (window [u+1, 12-u] mod 13).  So:
  PLANE KILL (fixed C, perm, k_min; all k_mid, k_max at once): find q = 13u
    (u in 2..5), a with: base9 and w_min clear at FULL q-precision AND
    a*p_mid, a*p_max mod 13 in [u+1, 12-u].  O(q) per try, kills O(10^7) cells.
  ROW KILL (fixed k_mid; all k_max): same with only p_max orbit-reduced; w_mid
    checked at full precision.  General-q witnesses also clear by k_max classes
    with the uncleared-set trick.
  CELL: adaptive per-set scan (set-specific denominators), then exact M.
Zero cells may end < 2/25.  Everything is exact integer arithmetic.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations, permutations
import sys, time

sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(m=""):
    print(m, flush=True)

BETA = F(2, 25)
AP = list(range(1, 13))
R1, R2 = 24, F(1100, 51)
A3 = F(3600, 13)

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

def dq(x, q):
    x %= q
    return min(x, q - x)

def sieve_ok(W):
    return all(any(v % m == 0 for v in W) for m in range(2, 13))

def ok(x, q):                      # margin >= 2/25 at value-residue x mod q
    return dq(x, q) * 25 >= 2 * q

QGEN = [q for q in range(8, 42) if q % 13] + [13 * u for u in range(2, 22)] + [25, 50]

stats = dict(planes=0, plane_killed=0, rows=0, row_killed=0, cells=0,
             cell_class=0, cell_scan=0, cell_exact=0, sieved=0, nonprim=0)
sub, unresolved = [], []
t_last = [T0]
nC = 0

for C in combinations(range(1, 13), 3):
    nC += 1
    base9 = [v for v in AP if v not in C]
    forced = [r for r in C if r >= 7]
    # -------- per-C: the 13u plane/row-kill witness families (u = 2..5) -----
    # for each u, viable a's: base9 clear at full 13u precision
    plane_wit = {}   # u -> list of a with base9 clear mod 13u
    for u in (2, 3, 4, 5):
        q = 13 * u
        plane_wit[u] = [a for a in range(1, q // 2 + 1) if gcd(a, q) == 1
                        and all(ok(a * v, q) for v in base9)]
    # general per-C witnesses (base9 clear), used at row/cell level
    cg = []
    for q in QGEN:
        if any(v % q == 0 for v in base9):
            continue
        cg.extend((q, a) for a in range(1, q // 2 + 1)
                  if gcd(a, q) == 1 and all(ok(a * v, q) for v in base9))
    for pmin, pmid, pmax in permutations(C):
        for ka in range(1, int(float(A3) - pmin) // 13 + 2):
            wa = pmin + 13 * ka
            if wa > A3:
                break
            stats['planes'] += 1
            # ---- PLANE KILL: q = 13u, a: wa clear full-q; pmid,pmax windows
            plane_dead = False
            for u in (2, 3, 4, 5):
                q = 13 * u
                for a in plane_wit[u]:
                    if wa % q == 0 or not ok(a * wa, q):
                        continue
                    if dq(a * pmid, 13) >= u + 1 and dq(a * pmax, 13) >= u + 1:
                        plane_dead = True     # every (k_mid,k_max) cleared
                        break
                if plane_dead:
                    break
            if plane_dead:
                stats['plane_killed'] += 1
                continue
            wag = [(q, a) for (q, a) in cg if wa % q and ok(a * wa, q)]
            wag_u = {u: [a for a in plane_wit[u] if wa % (13 * u) and ok(a * wa, 13 * u)]
                     for u in (2, 3, 4, 5)}
            for kb in range(1, int(float(R2) * wa) // 13 + 2):
                wb = pmid + 13 * kb
                if wb < wa:
                    continue
                if wb > R2 * wa:
                    break
                stats['rows'] += 1
                # ---- ROW KILL: q = 13u with wb clear full-q; pmax window
                row_dead = False
                for u in (2, 3, 4, 5):
                    q = 13 * u
                    for a in wag_u[u]:
                        if wb % q == 0 or not ok(a * wb, q):
                            continue
                        if dq(a * pmax, 13) >= u + 1:
                            row_dead = True
                            break
                    if row_dead:
                        break
                if row_dead:
                    stats['row_killed'] += 1
                    continue
                # ---- segmented k_max clearing over the finite range
                kclo = max(1, -(-(wb - pmax) // 13))
                kchi = (24 * wb - pmax) // 13
                if kchi < kclo:
                    continue
                need = [r for r in forced if wa % r and wb % r]
                if need:
                    # k_max must satisfy 13k == -pmax mod r for all r in need
                    Rr, Mm = 0, 1
                    bad = False
                    for r in need:
                        rr = (-pmax * pow(13 % r, -1, r)) % r
                        g2 = gcd(Mm, r)
                        if (rr - Rr) % g2:
                            bad = True
                            break
                        t = ((rr - Rr) // g2 * pow(Mm // g2, -1, r // g2)) % (r // g2)
                        Rr, Mm = Rr + Mm * t, Mm // g2 * r
                        Rr %= Mm
                    if bad:
                        continue
                    first = kclo + ((Rr - kclo) % Mm)
                    cand = list(range(first, kchi + 1, Mm))
                else:
                    cand = list(range(kclo, kchi + 1))
                if not cand:
                    continue
                stats['cells'] += len(cand)
                unc = set(range(len(cand)))
                for q, a in wag:
                    if not unc:
                        break
                    if wb % q == 0 or not ok(a * wb, q):
                        continue
                    step = (13 * a) % q
                    b0 = (pmax * a) % q
                    if step == 0:
                        if ok(b0, q):
                            unc.clear()
                        continue
                    g = gcd(step, q)
                    period = q // g
                    goodt = {t for t in range(period) if ok(b0 + step * t, q)}
                    if not goodt:
                        continue
                    rm = [i for i in unc if (cand[i] % period) in goodt]
                    for i in rm:
                        unc.discard(i)
                        stats['cell_class'] += 1
                # ---- leftover cells
                for i in sorted(unc):
                    kc = cand[i]
                    wc = pmax + 13 * kc
                    W = sorted(base9 + [wa, wb, wc])
                    if not sieve_ok(W):
                        stats['sieved'] += 1
                        continue
                    if reduce(gcd, W) != 1:
                        stats['nonprim'] += 1
                        continue
                    found = None
                    for q in ([wc + 1, wc - 1, wb + wc, wc - wb, wa + wb, wb + 1,
                               wa + wc, wc - wa, wb - wa, wa + 1] +
                              [13 * u for u in range(2, 40)] + [25, 50]):
                        if q < 8 or any(v % q == 0 for v in W):
                            continue
                        for a in range(1, q // 2 + 1):
                            if gcd(a, q) != 1:
                                continue
                            if all(ok(a * v, q) for v in W):
                                found = (q, a)
                                break
                        if found:
                            break
                    if found:
                        stats['cell_scan'] += 1
                        continue
                    stats['cell_exact'] += 1
                    if max(W) <= 3000:
                        M = M_exact(W)
                        if M < BETA:
                            sub.append((M, tuple(W), C, (ka, kb, kc)))
                            log(f"  << SUB-2/25: M={M} W={list(W)}")
                    else:
                        unresolved.append((tuple(W), C))
                        log(f"  ?? UNRESOLVED big: W={list(W)}")
            if time.time() - t_last[0] > 90:
                t_last[0] = time.time()
                log(f"  ...C {nC}/220 {C}  {stats}  t={time.time()-T0:.0f}s")

log(f"\nl=3 sweep complete: {stats}")
sub.sort()
log(f"sub-2/25 sets: {len(sub)}; unresolved: {len(unresolved)}")
for M, W, C, ks in sub[:20]:
    log(f"   M = {str(M):>10}  W = {list(W)}")
log("l=3 RESULT: " + ("FLOOR >= 2/25 on the FULL chain domain (zero sub-2/25, zero unresolved)"
                      if not sub and not unresolved else "see lists above"))
log(f"[t = {time.time()-T0:.0f}s]")
