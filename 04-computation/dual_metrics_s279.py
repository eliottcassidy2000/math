#!/usr/bin/env python3
"""
dual_metrics_s279.py — Dual metrics on the merged metagraph
opus-2026-03-24-S279

TWO METRICS ON THE SAME NODES:
  1. Blue/black lines (complement pairs) → metagraph edges with weights
  2. Wiggly lines (tile flips) → wiggly edges with weights

CONJECTURE FOR #TILINGS:
  By orbit-stabilizer on Aut(T) acting on Hamiltonian paths:
  - No non-identity automorphism fixes a Hamiltonian path
  - So: #tilings(unmerged C) = H(C) / |Aut(C)|
  - For merged: #tilings = H/|Aut| (SC) or 2H/|Aut| (NS)

This script:
1. Verifies the #tilings formula at n=4,5,6
2. Computes blue/black line weights per metagraph edge
3. Computes wiggly line weights per metagraph edge
4. Checks if the two weight systems are proportional or related
5. Derives the tiling count from the dual metrics
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def build_full(n):
    """Build everything: classes, tilings, automorphisms, H."""
    m_tiles = comb(n-1, 2)
    m_total = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]

    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    tile_to_P = {}
    base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1:
            base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k
                    break

    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines = [l.strip() for l in (r.stdout or '').split('\n')
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
    def aut_size(a):
        cnt = 0
        for p in permutations(range(n)):
            if all(a[p[i]][p[j]] == a[i][j] for i in range(n) for j in range(n)):
                cnt += 1
        return cnt
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

    cls = []
    hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        aut = aut_size(adj) if n <= 7 else 1
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'adj':adj, 'score':s, 'c3':cc, 'H':h,
                    'sc':sc, 'comp_canon':cc2, 'canon':canon, 'aut':aut})
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
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'aut': cl0['aut']}

    return {
        'n': n, 'm_tiles': m_tiles, 'm_total': m_total,
        'tiles': tiles, 'tile_to_P': tile_to_P, 'base_P': base_P,
        'P': P, 'cls': cls, 'lk': lk, 'mn': mn, 'V': mid, 'mi': mi, 'b2a': b2a,
    }

def analyze_dual(data):
    """Compute both metrics and the tiling formula."""
    n = data['n']; V = data['V']; mi = data['mi']
    m_tiles = data['m_tiles']; tiles = data['tiles']
    lk = data['lk']; mn = data['mn']; b2a = data['b2a']
    P = data['P']; base_P = data['base_P']; tile_to_P = data['tile_to_P']
    m_total = data['m_total']

    print(f"\n{'#'*72}")
    print(f"  n = {n}: DUAL METRICS ON MERGED METAGRAPH (V={V})")
    print(f"{'#'*72}")

    nf = factorial(n)

    # 1. VERIFY THE TILING FORMULA
    print(f"\n  TILING FORMULA: #tilings = H/|Aut| × (2 if NS, 1 if SC)")

    # Count actual tilings per merged class
    actual_tilings = Counter()
    for tiling_bits in range(2 ** m_tiles):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k in base_P:
                pass  # base path: bit = 0
            else:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tiling_bits >> (m_tiles - 1 - ti)) & 1
                        if bit:
                            full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits)
        cid = lk(adj)
        if cid >= 0:
            actual_tilings[mn[cid]] += 1

    print(f"  {'Class':>6} {'H':>5} {'|Aut|':>6} {'SC?':>4} {'pred':>6} {'actual':>7} {'match':>6}")
    all_match = True
    for mm in sorted(mi.keys()):
        h = mi[mm]['H']
        aut = mi[mm]['aut']
        sc = mi[mm]['sc']
        predicted = h // aut if sc else 2 * h // aut
        actual = actual_tilings.get(mm, 0)
        match = "✓" if predicted == actual else "✗"
        if predicted != actual: all_match = False
        print(f"  {mm:>6} {h:>5} {aut:>6} {'SC' if sc else 'NS':>4} {predicted:>6} {actual:>7} {match:>6}")

    if all_match:
        print(f"\n  *** FORMULA VERIFIED: #tilings = H/|Aut| × (2 if NS, 1 if SC) ***")
    else:
        print(f"\n  Formula has mismatches — check automorphism sizes")

    # Check total
    total_predicted = sum(
        mi[mm]['H'] // mi[mm]['aut'] * (1 if mi[mm]['sc'] else 2)
        for mm in mi
    )
    total_actual = sum(actual_tilings.values())
    print(f"\n  Total tilings: predicted={total_predicted}, actual={total_actual}, "
          f"expected=2^{m_tiles}={2**m_tiles}")

    # 2. COMPUTE WIGGLY EDGE WEIGHTS (tile flips between classes)
    wiggly_weights = defaultdict(int)  # (ma, mb) → count
    wiggly_self = defaultdict(int)     # ma → self-loop count

    for tiling_bits in range(2 ** m_tiles):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k in base_P:
                pass
            else:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tiling_bits >> (m_tiles - 1 - ti)) & 1
                        if bit:
                            full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits)
        cid = lk(adj)
        if cid < 0: continue
        ma = mn[cid]

        # Flip each tile
        for ti in range(m_tiles):
            flipped = tiling_bits ^ (1 << (m_tiles - 1 - ti))
            full_bits_f = 0
            for k, (i, j) in enumerate(P):
                if k in base_P:
                    pass
                else:
                    for tj in range(m_tiles):
                        if tile_to_P.get(tj) == k:
                            bit = (flipped >> (m_tiles - 1 - tj)) & 1
                            if bit:
                                full_bits_f |= (1 << (m_total - 1 - k))
                            break
            adj_f = b2a(full_bits_f)
            cid_f = lk(adj_f)
            if cid_f < 0: continue
            mb = mn[cid_f]

            if ma == mb:
                wiggly_self[ma] += 1
            else:
                wiggly_weights[(min(ma,mb), max(ma,mb))] += 1

    # 3. COMPARE WITH BLUE/BLACK WEIGHTS (complement-based)
    # Blue/black lines connect tiling to its complement
    # Complement tiling = flip ALL tiles
    comp_weights = defaultdict(int)
    comp_self = defaultdict(int)

    for tiling_bits in range(2 ** m_tiles):
        comp_bits = ((1 << m_tiles) - 1) ^ tiling_bits  # flip all tiles

        full_bits_orig = 0
        full_bits_comp = 0
        for k, (i, j) in enumerate(P):
            if k in base_P:
                pass
            else:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit_o = (tiling_bits >> (m_tiles - 1 - ti)) & 1
                        bit_c = (comp_bits >> (m_tiles - 1 - ti)) & 1
                        if bit_o: full_bits_orig |= (1 << (m_total - 1 - k))
                        if bit_c: full_bits_comp |= (1 << (m_total - 1 - k))
                        break

        adj_o = b2a(full_bits_orig)
        adj_c = b2a(full_bits_comp)
        cid_o = lk(adj_o)
        cid_c = lk(adj_c)
        if cid_o < 0 or cid_c < 0: continue
        ma = mn[cid_o]
        mb = mn[cid_c]

        if ma == mb:
            comp_self[ma] += 1
        else:
            comp_weights[(min(ma,mb), max(ma,mb))] += 1

    # 4. DISPLAY DUAL METRICS
    all_edges = set(wiggly_weights.keys()) | set(comp_weights.keys())
    me_edges = set(e for e in all_edges if e[0] != e[1])

    print(f"\n  DUAL METRIC TABLE:")
    print(f"  {'Edge':>10} {'H_a':>4} {'H_b':>4} {'wiggly_w':>9} {'comp_w':>7} {'ratio':>7}")

    for (a, b) in sorted(me_edges):
        ww = wiggly_weights.get((a, b), 0)
        cw = comp_weights.get((a, b), 0)
        ratio = ww / cw if cw > 0 else float('inf')
        print(f"  ({a:>3},{b:>3}) {mi[a]['H']:>4} {mi[b]['H']:>4} {ww:>9} {cw:>7} {ratio:>7.2f}")

    # Self-loop comparison
    print(f"\n  SELF-LOOP COMPARISON:")
    print(f"  {'Node':>5} {'H':>4} {'SC':>3} {'wiggly_SL':>10} {'comp_SL':>8} {'#tilings':>9}")
    for mm in sorted(mi.keys()):
        ws = wiggly_self.get(mm, 0)
        cs = comp_self.get(mm, 0)
        nt = actual_tilings.get(mm, 0)
        sc_str = 'SC' if mi[mm]['sc'] else 'NS'
        print(f"  {mm:>5} {mi[mm]['H']:>4} {sc_str:>3} {ws:>10} {cs:>8} {nt:>9}")

    # 5. THE KEY RELATIONSHIP
    print(f"\n  KEY RELATIONSHIPS:")

    # Wiggly self-loops = #tilings × (avg neutral arcs per tiling)
    # comp self-loops = #tilings where comp is in same class = tilings of SC classes
    for mm in sorted(mi.keys()):
        nt = actual_tilings.get(mm, 0)
        ws = wiggly_self.get(mm, 0)
        cs = comp_self.get(mm, 0)
        if nt > 0:
            avg_neutral = ws / nt
            print(f"  Class {mm} (H={mi[mm]['H']}): "
                  f"avg_neutral={avg_neutral:.2f}/{m_tiles}, "
                  f"comp_self={'=tilings' if cs == nt else f'{cs}/{nt}'}")

    return all_match

if __name__ == '__main__':
    print("=" * 72)
    print("  DUAL METRICS: BLUE/BLACK + WIGGLY")
    print("  opus-2026-03-24-S279")
    print("=" * 72)

    for n in [4, 5, 6]:
        t0 = time.time()
        data = build_full(n)
        ok = analyze_dual(data)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    print(f"\n{'='*72}")
    print("  THE TILING FORMULA")
    print(f"{'='*72}")
    print("""
  THEOREM: #tilings(merged class C) = H(C)/|Aut(C)| × multiplier
    where multiplier = 1 if SC, 2 if NS

  PROOF:
  1. A tiling = a tournament containing the base path P_0.
  2. Each Hamiltonian path P in T_C gives one relabeling σ with σ(T_C) ⊃ P_0.
  3. Two paths give the same tiling iff they're related by Aut(T_C).
  4. No non-identity automorphism fixes a Hamiltonian path.
  5. So #tilings(unmerged C) = H(C)/|Aut(C)| by orbit-stabilizer.
  6. For merged: SC classes count once, NS pairs count twice.

  This formula connects three quantities:
    #tilings (from wiggly model)
    H (from tournament theory)
    |Aut| (from group theory)

  And it means: the WIGGLY METAGRAPH DETERMINES H/|Aut| for every class!
  (Just count how many tilings land in each merged node.)
""")

    print("DONE.")
