#!/usr/bin/env python3
r"""
metagraph_line_census_kps_S66.py   (kind-pasteur-2026-07-07-S66, HYP-4967)

THE COMPLETE LINE/NODE STRUCTURE of the merged tournament metagraph (owner directive:
formulas for blue/black lines and pure-blue/pure-black/mixed nodes at each n, the even/odd
multiplicity, and the tiling->iso-class fibration solved via the blue/black line-pair partner).

STRICT DEFINITIONS (CLAUDE.md, klein-S13 -- NOT the legacy opus-S211 class-level meaning):
  * A TILING = a binary string on the m = C(n-1,2) non-base-path arcs (tiles).  Fixing the
    base Hamiltonian path P0 = n->n-1->...->1, each tiling is a labeled tournament containing
    P0 as a directed path.  There are 2^m tilings.
  * grid involution g: tile (x,y) [x>y, x-y>=2] <-> (n-y+1, n-x+1).  A tiling is GRID-SYM
    (BLUE) iff invariant under g; else BLACK.  This g is the tiling-level transpose, so a blue
    tiling = a self-complementary (transpose-fixed) tournament.
  * A LINE = a {t, flip(t)} pair, flip = complement TILING (flip ALL m tiles; stays in the
    cube, NOT T^op).  isGridSym(flip(t)) = isGridSym(t), so lines are monochromatic.
  * BLUE LINE = grid-sym tiling; BLACK LINE = not.  #lines = 2^{m-1}.
  * NODE = iso class.  PURE-BLUE = all its tilings blue; PURE-BLACK = none; MIXED = both.

OUTPUTS:
 (1) per-n census: |V|=A000568, SC, and pure-blue / pure-black / mixed node counts; blue/black
     tiling & line counts; verify BLUE(n)=2^{e(n)} (square/pronic) and #tilings(C)*|Aut(C)|=H(C).
 (2) the even/odd multiplicity: blue-tiling multiplicity per class (ODD over SC, 0 over NS).
 (3) LINE-NODE-TYPE INCIDENCE: for each line {t,flip(t)}, the (type,type) of the two endpoint
     CLASSES -- which node-type pairs do blue lines connect? black lines?  (the combination law)
 (4) the fibration solved: #tilings per class, and the flip-partner map on classes.
 (5) formulas + residuals; conjecture seeds for the pure-blue/mixed split.
"""
from itertools import combinations, permutations
from math import comb

def tiles_of(n):
    """non-base-path arcs (x,y), x>y, x-y>=2, in a fixed order."""
    T = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            if x - y >= 2:
                T.append((x, y))
    return T

def grid_perm(n, T):
    """index permutation of the grid involution g:(x,y)->(n-y+1,n-x+1) on the tile list."""
    pos = {t: i for i, t in enumerate(T)}
    return [pos[(n - y + 1, n - x + 1)] for (x, y) in T]

def tournament_from_tiling(n, T, bits):
    """adjacency: base path k->k-1 fixed; tile (x,y) bit 0 => x->y, bit 1 => y->x."""
    A = [[0]*(n+1) for _ in range(n+1)]           # 1-indexed
    for k in range(2, n + 1):
        A[k][k-1] = 1                              # base arc k -> k-1
    for (x, y), b in zip(T, bits):
        if b == 0: A[x][y] = 1
        else:      A[y][x] = 1
    return A

def canon(n, A):
    """canonical form under relabeling: min over permutations of the arc-bitstring.
    n<=6 brute (720); returns a hashable tuple."""
    best = None
    verts = list(range(1, n + 1))
    for p in permutations(verts):
        # p[i] is the new label of vertex i+1
        key = 0
        for a in range(n):
            for b in range(n):
                if a != b and A[p[a]][p[b]]:
                    key |= 1 << (a * n + b)
        if best is None or key < best:
            best = key
    return best

def aut_size(n, A):
    cnt = 0
    for p in permutations(range(1, n + 1)):
        ok = True
        for a in range(1, n+1):
            for b in range(1, n+1):
                if a != b and A[a][b] != A[p[a-1]][p[b-1]]:
                    ok = False; break
            if not ok: break
        if ok: cnt += 1
    return cnt

def ham_paths(n, A):
    dp = {}
    for i in range(1, n+1): dp[(1 << (i-1), i)] = 1
    full = (1 << n) - 1
    for mask in range(1 << n):
        for last in range(1, n+1):
            c = dp.get((mask, last), 0)
            if not c: continue
            for nx in range(1, n+1):
                if not (mask >> (nx-1)) & 1 and A[last][nx]:
                    k = (mask | (1 << (nx-1)), nx)
                    dp[k] = dp.get(k, 0) + c
    return sum(dp.get((full, i), 0) for i in range(1, n+1))

def is_SC(n, A):
    """self-complementary: T iso to its arc-reversal."""
    B = [[A[b][a] for b in range(n+1)] for a in range(n+1)]
    return canon(n, A) == canon(n, B)

print("=" * 100)
print("METAGRAPH LINE/NODE CENSUS (strict grid-sym definition)")
print("=" * 100)
NS_MAX = 6
summary = {}
for n in range(3, NS_MAX + 1):
    T = tiles_of(n); m = len(T)
    g = grid_perm(n, T)
    tilings = 1 << m
    # per-tiling data
    cls_of = {}                 # tiling index -> canonical class key
    blue_of = {}                # tiling index -> bool grid-sym
    class_reps = {}             # class key -> a representative A
    for code in range(tilings):
        bits = [(code >> i) & 1 for i in range(m)]
        A = tournament_from_tiling(n, T, bits)
        ck = canon(n, A)
        cls_of[code] = ck
        blue = all(bits[i] == bits[g[i]] for i in range(m))
        blue_of[code] = blue
        if ck not in class_reps: class_reps[ck] = A
    classes = sorted(class_reps)
    V = len(classes)
    # per-class aggregation
    per = {}
    for ck in classes:
        per[ck] = {"til": 0, "blue": 0, "A": class_reps[ck]}
    for code in range(tilings):
        d = per[cls_of[code]]; d["til"] += 1
        if blue_of[code]: d["blue"] += 1
    # invariants
    SC = 0; pureB = pureK = mixed = 0
    blue_tilings = sum(1 for c in range(tilings) if blue_of[c])
    for ck in classes:
        d = per[ck]; A = d["A"]
        d["H"] = ham_paths(n, A); d["aut"] = aut_size(n, A); d["SC"] = is_SC(n, A)
        if d["SC"]: SC += 1
        if d["blue"] == d["til"]: pureB += 1; d["type"] = "pureB"
        elif d["blue"] == 0:     pureK += 1; d["type"] = "pureK"
        else:                     mixed += 1; d["type"] = "mixed"
    # line/node type incidence
    blue_lines = blue_tilings // 2
    black_lines = (tilings - blue_tilings) // 2
    # flip-partner map on classes; incidence of (type,type) over LINES (count each line once)
    seen = set(); inc_blue = {}; inc_black = {}
    for code in range(tilings):
        fc = code ^ (tilings - 1)     # complement tiling = flip all m bits
        key = min(code, fc)
        if key in seen: continue
        seen.add(key)
        t1, t2 = per[cls_of[code]]["type"], per[cls_of[fc]]["type"]
        pair = tuple(sorted((t1, t2)))
        if blue_of[code]: inc_blue[pair] = inc_blue.get(pair, 0) + 1
        else:             inc_black[pair] = inc_black.get(pair, 0) + 1
    # e(n) closed form = floor((n-1)^2 / 4) = A002620(n-1); = k^2 (odd n=2k+1) / k(k-1) (even n=2k)
    e = ((n - 1) ** 2) // 4
    summary[n] = dict(V=V, SC=SC, pureB=pureB, pureK=pureK, mixed=mixed,
                      m=m, e=e, blue_til=blue_tilings, blue_lines=blue_lines,
                      black_lines=black_lines, per=per, classes=classes,
                      inc_blue=inc_blue, inc_black=inc_black)
    print(f"\n--- n={n}: |V|={V} (A000568), m={m}, tilings=2^{m}={tilings} ---")
    print(f"  SC={SC}, pure-blue={pureB}, mixed={mixed}, pure-black={pureK}  (pure-black == NS == {V-SC}: {pureK==V-SC})")
    print(f"  BLUE tilings = {blue_tilings} = 2^{e} ? {blue_tilings == (1<<e)}   (e(n)=floor((n-1)^2/4)={e}, {'sq' if n%2 else 'pronic'})")
    exp_bl = (1 << (e-1)) if e >= 1 else 0
    exp_bk = (1 << (m-1)) - exp_bl if m >= 1 else 0
    print(f"  blue lines = {blue_lines} = 2^(e-1)={exp_bl}? {blue_lines==exp_bl}; black lines = {black_lines} = 2^(m-1)-2^(e-1)={exp_bk}? {black_lines==exp_bk}")
    # verify #tilings * |Aut| = H
    fib_ok = all(per[ck]["til"] * per[ck]["aut"] == per[ck]["H"] for ck in classes)
    print(f"  fibration  #tilings(C)*|Aut(C)| == H(C) for all C: {fib_ok}")
    # even/odd multiplicity
    odd_over_SC = all((per[ck]["blue"] % 2 == 1) for ck in classes if per[ck]["SC"])
    zero_over_NS = all(per[ck]["blue"] == 0 for ck in classes if not per[ck]["SC"])
    print(f"  even/odd: blue-multiplicity ODD over every SC class: {odd_over_SC}; ==0 over every NS class: {zero_over_NS}")
    print(f"  LINE-NODE-TYPE incidence (endpoints' class types), each line once:")
    print(f"    blue lines by (type,type): {dict(sorted(inc_blue.items()))}")
    print(f"    black lines by (type,type): {dict(sorted(inc_black.items()))}")

print()
print("=" * 100)
print("FORMULA TABLE + residual patterns")
print("=" * 100)
print(f"  {'n':>2} {'V':>4} {'SC':>3} {'pureB':>5} {'mixed':>5} {'pureK/NS':>8} {'e(n)':>4} {'blueLines':>9} {'blackLines':>11}")
for n in range(3, NS_MAX+1):
    s = summary[n]
    print(f"  {n:>2} {s['V']:>4} {s['SC']:>3} {s['pureB']:>5} {s['mixed']:>5} {s['pureK']:>8} {s['e']:>4} {s['blue_lines']:>9} {s['black_lines']:>11}")
print()
print("  sequences to fit:")
print(f"    SC(n):        {[summary[n]['SC'] for n in range(3,NS_MAX+1)]}")
print(f"    pure-blue(n): {[summary[n]['pureB'] for n in range(3,NS_MAX+1)]}  (+ n=7 = 4 from prior work; seq 2,1,3,2,4)")
print(f"    mixed(n):     {[summary[n]['mixed'] for n in range(3,NS_MAX+1)]}")
print(f"    pure-blue+mixed = SC: {[summary[n]['pureB']+summary[n]['mixed']==summary[n]['SC'] for n in range(3,NS_MAX+1)]}")

print()
print("=" * 100)
print("PURE-BLUE STRUCTURE: characterize pure-blue classes (H, |Aut|, tc=H/|Aut|, blue-mult)")
print("=" * 100)
for n in range(3, NS_MAX+1):
    s = summary[n]; per = s['per']
    print(f"  n={n} pure-blue classes:")
    for ck in s['classes']:
        d = per[ck]
        if d['type'] == 'pureB':
            print(f"    H={d['H']:>3} |Aut|={d['aut']:>2} tc={d['H']//d['aut']:>3} blue-mult={d['blue']:>3} (tc==1? {d['H']==d['aut']})")
    # is pure-blue EXACTLY the tc==1 (rigid) SC classes?
    tc1_SC = sum(1 for ck in s['classes'] if per[ck]['SC'] and per[ck]['H']==per[ck]['aut'])
    print(f"    => #SC classes with tc==1 (H==|Aut|): {tc1_SC}   vs pure-blue={s['pureB']}   match: {tc1_SC==s['pureB']}")

print()
print("=" * 100)
print("FLIP-PARTNER CLASS MAP: does the complement-tiling flip relate H / |Aut| / class of endpoints?")
print("=" * 100)
for n in range(3, NS_MAX+1):
    s = summary[n]; per = s['per']
    # rebuild flip-partner class pairs at the CLASS level (which class -> which under flip, per line)
    T = tiles_of(n); m = s['m']; tilings = 1 << m
    g = grid_perm(n, T)
    same_class = 0; diff_class = 0; H_preserved = 0; total = 0
    seen = set()
    # need cls_of again -- recompute compactly
    cls_of = {}
    for code in range(tilings):
        bits = [(code >> i) & 1 for i in range(m)]
        A = tournament_from_tiling(n, T, bits)
        cls_of[code] = canon(n, A)
    for code in range(tilings):
        fc = code ^ (tilings - 1)
        key = min(code, fc)
        if key in seen: continue
        seen.add(key); total += 1
        c1, c2 = cls_of[code], cls_of[fc]
        if c1 == c2: same_class += 1
        else: diff_class += 1
        if per[c1]['H'] == per[c2]['H']: H_preserved += 1
    print(f"  n={n}: {total} lines; flip stays in same class (self-loop): {same_class}, crosses: {diff_class};"
          f" H(endpoint1)==H(endpoint2): {H_preserved}/{total} ({'ALWAYS' if H_preserved==total else 'NOT always'})")
print("DONE.")
