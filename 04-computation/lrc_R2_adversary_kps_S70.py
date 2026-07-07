#!/usr/bin/env python3
r"""
lrc_R2_adversary_kps_S70.py   (kind-pasteur-2026-07-07-S70, HYP-5077)

R2 STRESS TEST: are spread APs the GLOBAL PA_2-minimizer? Test HARD before any proof
(fleet discipline: Schur / 2-7 floor / E[maxgap]-AP-min were all refuted by strong
adversaries; a weak adversary that only tries APs would trivially "confirm" R2).

PA_2(E) = meas{x: max(gap@0(x), gap@1/2(x)) > 1/7}.  Two questions:
  Q1 (the one that matters for A'): does ANY k-family have PA_2 < T_k?  (If yes the
     2-anchor route is insufficient -- need a 3rd anchor or the full mu.)
  Q2 (R2 as stated): does ANY NON-AP family beat the spread-AP min (0.7202 at k=8)?

Adversary: structured families (spread/dilated/defected APs, two-block, geometric, Sidon,
near-AP single/double perturbations, mult-of-7 heavy, parity-interlaced) + random + a
LOCAL-DESCENT PA_2-minimizer from many seeds (the move-one-speed neighborhood), for
k=8 (tightest T_k=0.6185) and k=13.
"""
import random, math
from fractions import Fraction as F

def PA2(E, res=16000):
    E = sorted(set(E)); c = 0
    for r in range(res):
        x = (r + .5) / res
        ph = sorted((e * x) % 1.0 for e in E); n = len(ph)
        # gap@0 and gap@1/2
        def gapc(a):
            for i in range(n-1):
                if ph[i] <= a < ph[i+1]: return ph[i+1]-ph[i]
            return ph[0] + 1 - ph[-1]
        if max(gapc(0.0), gapc(0.5)) > 1/7: c += 1
    return c / res

Tk = {8: 0.6185, 13: 0.0565}
spreadAP_min = {8: 0.7202, 13: 0.3187}

rng = random.Random(70)
def structured(k, H):
    fams = {}
    # spread / dilated / defected APs
    for d in (1,2,3,5,7):
        fams[f"AP a=1 d={d}"] = [1 + d*j for j in range(k)]
        fams[f"AP a=5 d={d}"] = [5 + d*j for j in range(k)]
    fams["AP+outlier"] = [1+2*j for j in range(k-1)] + [1+2*(k-1)+50]
    fams["AP-defect"] = [1,2,3,4,5,6,7,9][:k] if k<=8 else list(range(1,k))+[k+3]
    fams["two-block"] = list(range(1,k//2+1)) + list(range(H-k+k//2+1, H+1))
    fams["geometric"] = [2**j for j in range(k)]
    fams["Sidon-greedy"] = ([1,2,4,8,13,21,31,45,66,81,97,123,148])[:k]
    fams["primes"] = ([2,3,5,7,11,13,17,19,23,29,31,37,41])[:k]
    fams["mult7-heavy"] = [7*j for j in range(1,k+1)]
    fams["parity-interlace"] = sorted([2*j for j in range(1,k//2+2)] + [2*j-1 for j in range(1,k-k//2)])[:k]
    # near-AP single/double perturbations of the d=2 dip family
    base = [5 + 2*j for j in range(k)]
    for _ in range(40):
        E = base[:]
        for _ in range(rng.randrange(1,4)):
            E[rng.randrange(k)] += rng.choice([-2,-1,1,2,3,50])
        E = sorted(set(e for e in E if e > 0))
        if len(E) == k: fams[f"nearAP-{_}"] = E
    return fams

def descend(k, H, iters=200, rs=0):
    rng2 = random.Random(rs)
    cur = sorted(rng2.sample(range(1, H), k))
    cv = PA2(cur, 6000)
    for _ in range(iters):
        i = rng2.randrange(k)
        cand = cur[:]; cand[i] += rng2.choice([-3,-2,-1,1,2,3])
        cand = sorted(set(e for e in cand if 0 < e < H))
        if len(cand) != k: continue
        v = PA2(cand, 6000)
        if v < cv - 1e-6: cur, cv = cand, v
        elif rng2.random() < 0.1:  # occasional sideways
            cur = sorted(rng2.sample(range(1,H), k))
            cv = PA2(cur, 6000)
    return cv, cur

for k in (8, 13):
    T = Tk[k]; sa = spreadAP_min[k]; H = 8*k
    print("=" * 90)
    print(f"R2 STRESS TEST k={k}: T_k={T:.4f}, spread-AP min={sa:.4f}")
    print("=" * 90)
    best_all = (9.9, None, None)
    below_Tk = []
    beat_sa = []
    for nm, E in structured(k, H).items():
        v = PA2(E)
        if v < best_all[0]: best_all = (v, nm, E)
        if v < T - 0.003: below_Tk.append((v, nm, E))
        if v < sa - 0.005: beat_sa.append((v, nm, E))
    # descent
    dmin = (9.9, None)
    for rs in range(25):
        v, E = descend(k, H, 160, rs)
        if v < dmin[0]: dmin = (v, E)
        if v < best_all[0]: best_all = (v, "descent", E)
        if v < T - 0.003: below_Tk.append((v, "descent", E))
        if v < sa - 0.005: beat_sa.append((v, "descent", E))
    print(f"  descent min PA_2 = {dmin[0]:.4f} at {dmin[1]}")
    print(f"  GLOBAL min found = {best_all[0]:.4f} ({best_all[1]}: {best_all[2]})")
    print(f"  Q1 -- any family below T_k={T:.4f}? {'YES *** (2-anchor route FAILS)' if below_Tk else 'NO (2-anchor route holds)'}")
    if below_Tk:
        for v,nm,E in sorted(below_Tk)[:3]: print(f"      {v:.4f} {nm}: {E}")
    print(f"  Q2 -- any NON-AP beats spread-AP min {sa:.4f}? {'YES (R2 as-stated false)' if beat_sa else 'NO (R2 survives: spread AP is the min)'}")
    if beat_sa:
        for v,nm,E in sorted(beat_sa)[:5]:
            st = sorted(E); diffs = [st[i+1]-st[i] for i in range(len(st)-1)]
            isAP = len(set(diffs)) == 1
            print(f"      {v:.4f} {nm}: {E}  (AP? {isAP})")
    print()
print("DONE.")
