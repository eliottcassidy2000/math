#!/usr/bin/env python3
"""ADVERSARIAL RE-VERIFICATION (independent, fresh code) of LEM-003 / HYP-2369 claims
(session kind-pasteur-2026-06-10-S1, thread B-universal-freeness).

Claims under test:
  B1: |Aut(T)| divides H(T) for ALL 2^10 labeled n=5 and 2^15 labeled n=6 tournaments;
      0 failures; 0 cases with H=0; all H odd (Redei); Burnside sums
      sum_T |Aut(T)| = 120*12 (n=5) and 720*56 (n=6).
  B2: orbits of the Aut-action on directed Ham paths have size EXACTLY |Aut|
      (freeness) -- checked here for EVERY mask with |Aut|>1 at n=5 and n=6
      (worker claims 184 and 3248 such masks); n=6 (|Aut|,H,#orbits) profile
      claimed = {(3,3,1),(3,15,5),(3,27,9),(3,33,11),(3,45,15),(5,15,3),(9,9,1)};
      |Aut| values at n=6 claimed {1,3,5,9}; max H(6)=45.
  B3: cycle caveat: C3 has 1 Ham cycle, |Aut|=3, 3 does not divide 1;
      RQ5 (i->i+1,i+2 mod 5) has exactly 2 Ham cycles, both rotation-fixed,
      orbit sizes [1,1], |Aut|=5, 5 does not divide 2; same tournament's PATH
      action is free: H=15, 3 orbits of size 5.

METHODS -- deliberately DIFFERENT from the worker's aut_divides_H_freeness_kpc1.py
(which used Held-Karp bitmask DP for H and, per its claims, per-mask scans for Aut):
  * H(T) for ALL masks simultaneously: each of the n! vertex sequences fixes the
    orientation of its n-1 consecutive pairs; it contributes +1 to exactly the
    2^(C(n,2)-(n-1)) masks agreeing with that pattern, enumerated by submask
    iteration over the free pairs.  No DP.
  * |Aut(T)| for ALL masks simultaneously: each sigma in S_n induces a signed
    permutation of the C(n,2) pair-bits; T is fixed by sigma iff every cycle of
    that signed permutation has even flip parity, and then the fixed masks are
    exactly 2^(#cycles) explicitly enumerable patterns.  Increments aut[T] only
    on actual fixed points (total increments = n! * #iso-classes, the Burnside
    identity, which doubles as a structural cross-check).
  * Independent spot-check of Aut by direct arc-preservation testing on a random
    sample of masks.
  * Freeness: explicit orbit construction on Ham-path vertex sequences for EVERY
    mask with |Aut|>1 (superset of the worker's sample).
All arithmetic exact integers.  Convention: bit k of mask T corresponds to
pairs[k]=(i,j), i<j (combinations order); bit=1 means arc i->j, bit=0 arc j->i.
"""

import itertools
import random
import sys
import time

A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56}  # tournament iso classes


def analyze(n):
    pairs = list(itertools.combinations(range(n), 2))
    NB = len(pairs)
    idx = {p: k for k, p in enumerate(pairs)}
    NMASK = 1 << NB
    FULL = NMASK - 1
    perms = list(itertools.permutations(range(n)))

    # ---------- H for all masks: pattern accumulation ----------
    H = [0] * NMASK
    perm_req = []  # (pairmask, pattern) per vertex sequence
    for p in perms:
        pm = 0
        pv = 0
        for a in range(n - 1):
            u, v = p[a], p[a + 1]
            if u < v:
                k = idx[(u, v)]
                pm |= 1 << k
                pv |= 1 << k
            else:
                k = idx[(v, u)]
                pm |= 1 << k
        perm_req.append((pm, pv))
        free = FULL ^ pm
        s = free
        while True:
            H[pv | s] += 1
            if s == 0:
                break
            s = (s - 1) & free

    # ---------- Aut for all masks: signed-cycle fixed-point enumeration ----------
    aut = [0] * NMASK
    autlist = {}  # mask -> list of fixing perms (only filled where some non-id fixes... we store all)
    total_fix = 0
    for sig in perms:
        tgt = [0] * NB
        flip = [0] * NB
        for k, (i, j) in enumerate(pairs):
            si, sj = sig[i], sig[j]
            if si < sj:
                tgt[k] = idx[(si, sj)]
            else:
                tgt[k] = idx[(sj, si)]
                flip[k] = 1
        seen = [False] * NB
        cycles = []
        consistent = True
        for k0 in range(NB):
            if seen[k0]:
                continue
            cyc = []
            fl = []
            k = k0
            while not seen[k]:
                seen[k] = True
                cyc.append(k)
                fl.append(flip[k])
                k = tgt[k]
            if sum(fl) % 2 == 1:
                consistent = False
                break
            cycles.append((cyc, fl))
        if not consistent:
            continue
        # enumerate all fixed masks: 2 starting values per cycle, propagate
        masks = [0]
        for cyc, fl in cycles:
            pats = []
            for start in (0, 1):
                m = 0
                val = start
                for t in range(len(cyc)):
                    if val:
                        m |= 1 << cyc[t]
                    val ^= fl[t]
                pats.append(m)
            masks = [m | q for m in masks for q in pats]
        for m in masks:
            aut[m] += 1
            autlist.setdefault(m, []).append(sig)
        total_fix += len(masks)

    # Burnside structural identity
    burnside_classes = total_fix // len(perms)
    assert total_fix == sum(aut), "internal accounting"

    # ---------- divisibility / parity / H=0 ----------
    div_fail = 0
    h0 = 0
    h_even = 0
    for T in range(NMASK):
        if H[T] == 0:
            h0 += 1
        if H[T] % 2 == 0:
            h_even += 1
        if H[T] % aut[T] != 0:
            div_fail += 1
    maxH = max(H)
    aut_values = sorted(set(aut))

    # ---------- independent direct Aut spot-check ----------
    rng = random.Random(20260610)
    sample = [rng.randrange(NMASK) for _ in range(min(60, NMASK))]
    for T in sample:
        cnt = 0
        for sig in perms:
            good = True
            for k, (i, j) in enumerate(pairs):
                b = (T >> k) & 1
                u, v = (i, j) if b else (j, i)  # arc u->v in T
                su, sv = sig[u], sig[v]
                if su < sv:
                    if not ((T >> idx[(su, sv)]) & 1):
                        good = False
                        break
                else:
                    if (T >> idx[(sv, su)]) & 1:
                        good = False
                        break
            if good:
                cnt += 1
        assert cnt == aut[T], f"Aut mismatch at n={n} T={T}: direct {cnt} vs cycle-method {aut[T]}"

    # ---------- freeness: explicit orbits for EVERY |Aut|>1 mask ----------
    big = [T for T in range(NMASK) if aut[T] > 1]
    free_fail = 0
    profile = {}
    for T in big:
        paths = [p for p, (pm, pv) in zip(perms, perm_req) if (T & pm) == pv]
        assert len(paths) == H[T], "H cross-check"
        A = autlist[T]
        assert len(A) == aut[T]
        pset = set(paths)
        seenp = set()
        norb = 0
        ok = True
        for p in paths:
            if p in seenp:
                continue
            orb = set()
            for sig in A:
                q = tuple(sig[v] for v in p)
                if q not in pset:
                    ok = False  # image of a Ham path must be a Ham path
                orb.add(q)
            seenp |= orb
            norb += 1
            if len(orb) != len(A):
                ok = False
        if not ok or norb * aut[T] != H[T]:
            free_fail += 1
        key = (aut[T], H[T], H[T] // aut[T])
        profile[key] = profile.get(key, 0) + 1

    return {
        "n": n, "NMASK": NMASK, "div_fail": div_fail, "h0": h0, "h_even": h_even,
        "sum_aut": total_fix, "burnside_classes": burnside_classes,
        "n_big": len(big), "free_fail": free_fail, "profile": profile,
        "maxH": maxH, "aut_values": aut_values,
    }


def report(r):
    print(f"--- n={r['n']} : {r['NMASK']} labeled tournaments ---")
    print(f"  divisibility |Aut| | H failures : {r['div_fail']}")
    print(f"  H=0 cases                      : {r['h0']}")
    print(f"  H even cases (Redei says 0)    : {r['h_even']}")
    exp = A000568[r['n']] * len(list(itertools.permutations(range(r['n']))))
    print(f"  Burnside sum |Aut| = {r['sum_aut']}  (expected n!*A000568 = {exp}, "
          f"classes = {r['burnside_classes']} vs A000568 {A000568[r['n']]})")
    print(f"  masks with |Aut|>1             : {r['n_big']}")
    print(f"  freeness failures (all |Aut|>1 masks, explicit orbits): {r['free_fail']}")
    print(f"  |Aut| values occurring         : {r['aut_values']}")
    print(f"  max H                          : {r['maxH']}")
    print(f"  (|Aut|,H,#orbits) profile over |Aut|>1 masks (with multiplicities):")
    for k in sorted(r['profile']):
        print(f"      {k} x {r['profile'][k]}")
    ok = (r['div_fail'] == 0 and r['h0'] == 0 and r['h_even'] == 0
          and r['sum_aut'] == exp and r['free_fail'] == 0)
    print(f"  PART VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


# ============================ PART 1+2: n=3..6 ============================
t0 = time.time()
print("=" * 72)
print("PART 1+2: exhaustive |Aut| | H + freeness, all labeled tournaments n=3..6")
print("=" * 72)
all_ok = True
results = {}
for n in (3, 4, 5, 6):
    r = analyze(n)
    results[n] = r
    all_ok &= report(r)

# worker-claim cross-checks
print()
print("worker-claim cross-checks:")
checks = [
    ("n=5 masks with |Aut|>1 == 184", results[5]["n_big"] == 184),
    ("n=6 masks with |Aut|>1 == 3248", results[6]["n_big"] == 3248),
    ("n=6 |Aut| values == {1,3,5,9}", results[6]["aut_values"] == [1, 3, 5, 9]),
    ("n=6 max H == 45", results[6]["maxH"] == 45),
    ("n=6 profile keys == {(3,3,1),(3,15,5),(3,27,9),(3,33,11),(3,45,15),(5,15,3),(9,9,1)}",
     set(results[6]["profile"].keys()) ==
     {(3, 3, 1), (3, 15, 5), (3, 27, 9), (3, 33, 11), (3, 45, 15), (5, 15, 3), (9, 9, 1)}),
    ("n=5 Burnside sum == 1440", results[5]["sum_aut"] == 1440),
    ("n=6 Burnside sum == 40320", results[6]["sum_aut"] == 40320),
]
for label, ok in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    all_ok &= ok

# ============================ PART 3: cycle caveat ============================
print()
print("=" * 72)
print("PART 3: cycle caveat (B3) -- C3 and RQ5, Ham CYCLES vs Ham PATHS")
print("=" * 72)


def digraph_aut(n, arcs):
    A = []
    for sig in itertools.permutations(range(n)):
        if all((sig[u], sig[v]) in arcs for (u, v) in arcs):
            A.append(sig)
    return A


def ham_cycles(n, arcs):
    found = set()
    for p in itertools.permutations(range(1, n)):
        seq = (0,) + p
        if all((seq[i], seq[(i + 1) % n]) in arcs for i in range(n)):
            found.add(frozenset((seq[i], seq[(i + 1) % n]) for i in range(n)))
    return sorted(found, key=sorted)


def ham_paths(n, arcs):
    return [p for p in itertools.permutations(range(n))
            if all((p[i], p[i + 1]) in arcs for i in range(n - 1))]


def cycle_orbits(cycles, A):
    cset = set(cycles)
    seen = set()
    sizes = []
    for c in cycles:
        if c in seen:
            continue
        orb = set()
        for sig in A:
            img = frozenset((sig[u], sig[v]) for (u, v) in c)
            assert img in cset, "automorphism image of a Ham cycle must be a Ham cycle"
            orb.add(img)
        seen |= orb
        sizes.append(len(orb))
    return sizes


# C3
arcs3 = {(0, 1), (1, 2), (2, 0)}
A3 = digraph_aut(3, arcs3)
C3cyc = ham_cycles(3, arcs3)
rot3 = (1, 2, 0)
fixed3 = [c for c in C3cyc if frozenset((rot3[u], rot3[v]) for (u, v) in c) == c]
print(f"C3: |Aut| = {len(A3)} (expect 3), #Ham cycles = {len(C3cyc)} (expect 1), "
      f"rotation-fixed cycles = {len(fixed3)}")
print(f"C3: 3 divides 1?  {len(C3cyc) % len(A3) == 0}  (expect False -> caveat real)")
c3_ok = (len(A3) == 3 and len(C3cyc) == 1 and len(fixed3) == 1
         and len(C3cyc) % len(A3) != 0)
print(f"  C3 caveat verdict: {'PASS' if c3_ok else 'FAIL'}")
all_ok &= c3_ok

# RQ5: i -> i+1, i+2 (mod 5)
arcs5 = set()
for i in range(5):
    arcs5.add((i, (i + 1) % 5))
    arcs5.add((i, (i + 2) % 5))
A5 = digraph_aut(5, arcs5)
RQ5cyc = ham_cycles(5, arcs5)
rot5 = (1, 2, 3, 4, 0)
fixed5 = [c for c in RQ5cyc if frozenset((rot5[u], rot5[v]) for (u, v) in c) == c]
sizes5 = cycle_orbits(RQ5cyc, A5)
plus1 = frozenset((i, (i + 1) % 5) for i in range(5))
plus2 = frozenset((i, (i + 2) % 5) for i in range(5))
both_are_steps = set(RQ5cyc) == {plus1, plus2}
print(f"RQ5: |Aut| = {len(A5)} (expect 5; rotations only: "
      f"{sorted(A5) == sorted(tuple((j + s) % 5 for j in range(5)) for s in range(5))})")
print(f"RQ5: #Ham cycles = {len(RQ5cyc)} (expect 2), they are the +1/+2 step cycles: {both_are_steps}")
print(f"RQ5: rotation-fixed cycles = {len(fixed5)} (expect 2), orbit sizes = {sorted(sizes5)} (expect [1,1])")
print(f"RQ5: 5 divides 2?  {len(RQ5cyc) % len(A5) == 0}  (expect False -> maximally non-free)")
# path action on same tournament
P5 = ham_paths(5, arcs5)
pset5 = set(P5)
seenp = set()
psizes = []
for p in P5:
    if p in seenp:
        continue
    orb = {tuple(sig[v] for v in p) for sig in A5}
    assert orb <= pset5
    seenp |= orb
    psizes.append(len(orb))
print(f"RQ5 PATHS: H = {len(P5)} (expect 15), orbit sizes = {sorted(psizes)} (expect [5,5,5])")
rq5_ok = (len(A5) == 5 and len(RQ5cyc) == 2 and both_are_steps and len(fixed5) == 2
          and sorted(sizes5) == [1, 1] and len(RQ5cyc) % len(A5) != 0
          and len(P5) == 15 and sorted(psizes) == [5, 5, 5])
print(f"  RQ5 caveat verdict: {'PASS' if rq5_ok else 'FAIL'}")
all_ok &= rq5_ok

# ============================ summary ============================
print()
print("=" * 72)
print(f"OVERALL: {'ALL CHECKS PASS' if all_ok else 'SOME CHECKS FAILED'}   "
      f"({time.time() - t0:.1f}s)")
print("=" * 72)
sys.exit(0 if all_ok else 1)
