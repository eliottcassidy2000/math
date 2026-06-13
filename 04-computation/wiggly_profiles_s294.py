#!/usr/bin/env python3
"""
wiggly_profiles_s294.py — How tilings within the same class differ in wiggly neighborhoods
opus-2026-03-24-S294

For each tiling T in class C, the WIGGLY PROFILE is:
  profile(T) = (class(T_A), class(T_B), ..., class(T_m))
  where T_α = tiling with tile α flipped

Two tilings T, T' in the same class C satisfy T' = σ(T) for some σ ∈ S_n.
Then: profile(T')_{σ(α)} = profile(T)_α

So the profiles are EQUIVARIANT under S_n — they transform together.
The distinct profiles within a class = orbits of Aut(C) on profiles.
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
    m_total = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines_raw = [l.strip() for l in (r.stdout or '').split('\n')
                 if len(l.strip()) == m_total and all(c in '01' for c in l.strip())]
    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m_total-1-k)): a[i][j] = 1
            else: a[j][i] = 1
        return a
    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if a[i][j] and a[j][k] and a[k][i]: c += 1
                    if a[i][k] and a[k][j] and a[j][i]: c += 1
        return c
    def H(a):
        dp = {}
        for v in range(n): dp[(1<<v,v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: dp[(S|(1<<u),u)] = dp.get((S|(1<<u),u),0) + val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(a):
        sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[sc[v]].append(v)
        gs = [sg[s] for s in sorted(sg.keys())]
        best = None; cnt = 0
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gp(gs2[1:]): yield list(p)+rest
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
            cnt += 1
            if cnt > 500000: break
        return best
    cls = []; hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines_raw):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc, 'comp_canon':cc2, 'canon':canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)
    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)
    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1
    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]; mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}
    m_tiles = comb(n-1, 2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2: tiles.append((x, y))
    tile_to_P = {}; base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1: base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k; break
    total = 2 ** m_tiles; V = mid
    tc = [0] * total
    for tb in range(total):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m_tiles - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits); cid = lk(adj)
        tc[tb] = mn[cid] if cid >= 0 else -1
    return V, mi, tc, m_tiles, total, tiles

def analyze_profiles(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: WIGGLY PROFILES WITHIN ISO CLASSES")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total, tiles = build_data(n)
    labels = [chr(65+i) if i < 26 else chr(65+i-26)+"'" for i in range(m)]
    print(f"  Built in {time.time()-t0:.1f}s. V={V}, m={m}, 2^m={total}")

    # For each tiling, compute its wiggly profile
    fiber_tilings = defaultdict(list)  # class → list of tiling indices
    for tb in range(total):
        if tc[tb] >= 0:
            fiber_tilings[tc[tb]].append(tb)

    # For each class, compute ALL profiles and find distinct patterns
    for mm in sorted(mi.keys()):
        tils = fiber_tilings[mm]
        if not tils: continue
        h = mi[mm]['H']; sc = mi[mm]['sc']

        # Compute profiles
        profiles = []
        for tb in tils:
            profile = []
            for ti in range(m):
                flipped = tb ^ (1 << (m - 1 - ti))
                profile.append(tc[flipped])
            profiles.append(tuple(profile))

        # Distinct profiles
        profile_counts = Counter(profiles)
        n_distinct = len(profile_counts)

        # The SORTED profile (class-wise, ignoring which tile goes where)
        sorted_profiles = Counter()
        for p in profiles:
            sorted_profiles[tuple(sorted(p))] += 1
        n_sorted_distinct = len(sorted_profiles)

        # The CLASS MULTISET of each profile (which classes appear, with multiplicity)
        class_multisets = Counter()
        for p in profiles:
            ms = tuple(sorted(Counter(p).items()))
            class_multisets[ms] += 1
        n_multiset_distinct = len(class_multisets)

        # Print for important classes
        if len(tils) <= 20 or mm < 3 or mm > V - 3 or h in [1, 3, 5] or n_distinct == 1:
            print(f"\n  Class {mm} (H={h}, {'SC' if sc else 'NS'}, fiber={len(tils)}):")
            print(f"    Distinct profiles: {n_distinct}")
            print(f"    Distinct sorted profiles: {n_sorted_distinct}")
            print(f"    Distinct class multisets: {n_multiset_distinct}")

            if n_distinct == 1:
                print(f"    *** ALL TILINGS HAVE IDENTICAL PROFILE! ***")
                print(f"    Profile: {profiles[0]}")
                print(f"    → classes visited: {sorted(set(profiles[0]))}")
                sl = sum(1 for x in profiles[0] if x == mm)
                print(f"    → self-loops: {sl}/{m}")

            elif n_distinct <= 5:
                for p, cnt in profile_counts.most_common():
                    classes_hit = sorted(set(p))
                    sl = sum(1 for x in p if x == mm)
                    print(f"    Profile ({cnt}×): classes={classes_hit}, self-loops={sl}")
                    if m <= 10:
                        print(f"      Detail: [{', '.join(str(x) for x in p)}]")

            else:
                # Show statistics
                sl_vals = [sum(1 for x in p if x == mm) for p in profiles]
                div_vals = [len(set(p)) for p in profiles]
                print(f"    Self-loops: min={min(sl_vals)}, max={max(sl_vals)}, "
                      f"mean={np.mean(sl_vals):.2f}")
                print(f"    Diversity: min={min(div_vals)}, max={max(div_vals)}, "
                      f"mean={np.mean(div_vals):.2f}")

                # Is self-loop count constant across tilings?
                if len(set(sl_vals)) == 1:
                    print(f"    → CONSTANT self-loop count: {sl_vals[0]}")
                else:
                    print(f"    → VARIABLE self-loop count: {sorted(Counter(sl_vals).items())}")

                # Is diversity constant?
                if len(set(div_vals)) == 1:
                    print(f"    → CONSTANT diversity: {div_vals[0]}")

    # SUMMARY: fraction of classes with constant profile/self-loops/diversity
    print(f"\n  SUMMARY:")
    constant_profile = 0
    constant_sl = 0
    constant_div = 0
    total_classes = 0

    for mm in sorted(mi.keys()):
        tils = fiber_tilings[mm]
        if not tils or len(tils) <= 1:
            constant_profile += 1; constant_sl += 1; constant_div += 1
            total_classes += 1
            continue

        profiles = []
        for tb in tils:
            profile = tuple(tc[tb ^ (1 << (m-1-ti))] for ti in range(m))
            profiles.append(profile)

        n_dist = len(set(profiles))
        sl_vals = set(sum(1 for x in p if x == mm) for p in profiles)
        div_vals = set(len(set(p)) for p in profiles)

        if n_dist == 1: constant_profile += 1
        if len(sl_vals) == 1: constant_sl += 1
        if len(div_vals) == 1: constant_div += 1
        total_classes += 1

    print(f"  Classes with constant profile: {constant_profile}/{total_classes} "
          f"({100*constant_profile/total_classes:.1f}%)")
    print(f"  Classes with constant self-loop count: {constant_sl}/{total_classes} "
          f"({100*constant_sl/total_classes:.1f}%)")
    print(f"  Classes with constant diversity: {constant_div}/{total_classes} "
          f"({100*constant_div/total_classes:.1f}%)")

if __name__ == '__main__':
    print("=" * 72)
    print("  WIGGLY PROFILES: HOW TILINGS IN THE SAME CLASS DIFFER")
    print("  opus-2026-03-24-S294")
    print("=" * 72)

    for n in [4, 5, 6]:
        analyze_profiles(n)

    print("\nDONE.")
