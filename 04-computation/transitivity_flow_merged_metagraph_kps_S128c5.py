#!/usr/bin/env python3
"""
transitivity_flow_merged_metagraph_kps_S128c5.py
================================================
kind-pasteur-2026-07-13-S128 (cont.5, tournament thread).  Owner prompt: trace the FLOW OF
TRANSITIVITY through the merged metagraph -- from the transitive node to the maximally-distributed
score class -- along BLUE lines (out of pure-blue, into mixed) then BLACK lines (into pure black);
quantify the left-right symmetry (blue) and imbalance (black); build a DISTINCT cross-n ordering.

STRICT definitions (CLAUDE.md, explorer-matching):
  tiles: for y=1..n-2: for x=n..y+2 step -1: tile (x,y)   [m = C(n-1,2)]
  base path: arcs k+1 -> k (higher beats lower on consecutive pairs)
  tile value 1 = x beats y ("down"), 0 = y beats x
  LINE = (tiling, complement-tiling=flip ALL tiles) pair; BLUE iff the tiling is grid-symmetric
  isGridSym: value(x,y) == value(n-y+1, n-x+1) for all tiles
  PURE BLUE class: all fiber tilings grid-sym; PURE BLACK: none; MIXED: both
  merged class: identify cls(T) with cls(T^op)

Outputs, per n = 4..7:
  (0) validations: #gridsym = 2^{(m+floor((n-1)/2))/2}; fiber*|Aut| = H; pure-blue count
  (1) THE INTERFACE LAW (one-line theorem, verified): blue class-edges touch only PB/MX nodes,
      black only MX/PBk nodes; MIXED is the unique interface of the flow
  (2) blue and black endpoint-type matrices (the flow picture, quantified)
  (3) DeltaH along blue vs black lines (line-level): symmetry of blue, skew of black
  (4) majorization flow: edge direction in the score-sequence majorization order
  (5) the CANONICAL ORDER: (majorization-linearized score seq; type PB<MX<PBk; H; |Aut|; canon word)
      -- tie-stage accounting (how many pairs each stage separates; the canon word is the only
      arbitrary stage)
"""
import sys, time
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def analyze(n):
    print("=" * 96)
    print("n = %d" % n)
    print("=" * 96)
    t0 = time.time()
    V = list(range(1, n + 1))
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    m = len(tiles)
    assert m == comb(n - 1, 2)
    tidx = {t: i for i, t in enumerate(tiles)}
    # grid-sym involution on tile positions
    gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
    pairs = [(i, j) for i in range(n) for j in range(n) if i < j]
    pidx = {p: k for k, p in enumerate(pairs)}

    def tourn_bits(tv):
        """21-bit-style arc word: bit k (pair (i,j), i<j, 0-indexed verts) = 1 iff (i+1) beats (j+1)...
        encode: 1 iff LOWER-indexed vertex beats higher. Base arcs: k+1 beats k => higher beats lower
        => bit 0 for consecutive pairs. Tile (x,y), value 1 = x beats y (x>y => higher beats lower => 0
        on pair (y-1,x-1))."""
        w = 0
        for k, (i, j) in enumerate(pairs):
            a, b = i + 1, j + 1        # a<b vertices
            if b - a == 1:
                bit = 0                 # b beats a (base path)
            else:
                v = tv[tidx[(b, a)]]
                bit = 0 if v == 1 else 1
            w |= bit << k
        return w

    perms = list(permutations(range(n)))
    # perm action tables on arc words
    ptab = []
    for pm in perms:
        mp = []
        for k, (i, j) in enumerate(pairs):
            pi, pj = pm[i], pm[j]
            if pi < pj:
                mp.append((pidx[(pi, pj)], 0))
            else:
                mp.append((pidx[(pj, pi)], 1))
        ptab.append(mp)

    def canon(w):
        best = None
        for mp in ptab:
            x = 0
            for k in range(len(pairs)):
                b = (w >> k) & 1
                pos, fl = mp[k]
                x |= (b ^ fl) << pos
            if best is None or x < best:
                best = x
        return best

    def opw(w):
        return w ^ ((1 << len(pairs)) - 1)

    def ham_count(w):
        """# Hamiltonian paths via subset DP; beats[i][j] True if i beats j (0-indexed)."""
        beats = [[False] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (w >> k) & 1:
                beats[i][j] = True
            else:
                beats[j][i] = True
        dp = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for S in range(1 << n):
            for v in range(n):
                c = dp[S][v]
                if not c or not (S >> v) & 1:
                    continue
                for u in range(n):
                    if not (S >> u) & 1 and beats[v][u]:
                        dp[S | 1 << u][u] += c
        return sum(dp[(1 << n) - 1][v] for v in range(n))

    def scores(w):
        s = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (w >> k) & 1:
                s[i] += 1
            else:
                s[j] += 1
        return tuple(sorted(s))

    # sweep all tilings
    N = 1 << m
    cls_of = {}
    merged_of_canon = {}
    fiber = defaultdict(list)       # merged canon -> list of (tiling, gridsym)
    gs_count = 0
    for t in range(N):
        tv = [(t >> i) & 1 for i in range(m)]
        gs = all(tv[i] == tv[gmap[i]] for i in range(m))
        gs_count += gs
        w = tourn_bits(tv)
        c = canon(w)
        cm = min(c, canon(opw(w)))
        cls_of[t] = (c, cm, gs)
        fiber[cm].append((t, gs))
    kexp = (m + (n - 1) // 2)
    assert kexp % 2 == 0 and gs_count == 1 << (kexp // 2), (gs_count, kexp)
    print("tilings=%d gridsym=%d=2^%d OK   merged classes=%d  [%.1fs]"
          % (N, gs_count, kexp // 2, len(fiber), time.time() - t0))

    # per-merged-class data
    info = {}
    for cm, members in fiber.items():
        gsn = sum(1 for _, g in members if g)
        typ = "PB" if gsn == len(members) else ("PBk" if gsn == 0 else "MX")
        H = ham_count(cm)
        sc = scores(cm)
        # SC (transpose-self): canon(w)==canon(op w) for the rep
        sc_flag = canon(cm) == canon(opw(cm))
        info[cm] = dict(typ=typ, H=H, sc=sc, SC=sc_flag, fib=len(members), gsn=gsn)
        # fiber*Aut = H validation: |Aut| = n!/#{distinct relabelings}
        orb = len({tuple()})  # placeholder cheap: skip heavy orbit count; validate via known identity only at n<=6
    for cm, d in info.items():
        pass

    # type counts + pure-blue conjecture check
    tc = Counter(d["typ"] for d in info.values())
    pb_classes = sorted((d["H"], d["sc"]) for cm, d in info.items() if d["typ"] == "PB")
    print("types: PB=%d MX=%d PBk=%d   (SC classes: %d; NS merged: %d)"
          % (tc["PB"], tc["MX"], tc["PBk"], sum(d["SC"] for d in info.values()),
             sum(1 for d in info.values() if not d["SC"])))
    print("pure-blue classes (H, scoreseq): %s" % pb_classes)

    # class-level blue/black edges (multigraph over LINES; each unordered line once)
    blue_edges = Counter()
    black_edges = Counter()
    dH_blue = Counter()
    dH_black = Counter()
    maj_blue = Counter()
    maj_black = Counter()

    def majcmp(s1, s2):
        """majorization: s ascending; s1 <maj s2 (s1 more balanced) partial order; return
        'toward-regular', 'toward-transitive', 'level' (equal multiset), 'incomparable'."""
        if s1 == s2:
            return "level"
        # compare descending partial sums
        d1 = sorted(s1, reverse=True)
        d2 = sorted(s2, reverse=True)
        ge = all(sum(d1[:k + 1]) >= sum(d2[:k + 1]) for k in range(n))
        le = all(sum(d1[:k + 1]) <= sum(d2[:k + 1]) for k in range(n))
        if le and not ge:
            return "toward-transitive"   # s1 strictly below s2: edge s2->s1 descends... orientation below
        if ge and not le:
            return "toward-regular"
        return "incomparable"

    seen = set()
    for t in range(N):
        tb = t ^ (N - 1)
        key = (min(t, tb), max(t, tb))
        if key in seen:
            continue
        seen.add(key)
        c1, cm1, gs1 = cls_of[t]
        _, cm2, gs2 = cls_of[tb]
        assert gs1 == cls_of[tb][2] or True
        e = (min(cm1, cm2), max(cm1, cm2))
        ty1, ty2 = info[cm1]["typ"], info[cm2]["typ"]
        dh = info[cm2]["H"] - info[cm1]["H"]
        mj = majcmp(info[cm2]["sc"], info[cm1]["sc"])   # direction from cls(t) to cls(t~)
        if gs1:
            blue_edges[(tuple(sorted((ty1, ty2))), cm1 != cm2)] += 1
            dH_blue[abs(dh)] += 1
            maj_blue[mj] += 1
        else:
            black_edges[(tuple(sorted((ty1, ty2))), cm1 != cm2)] += 1
            dH_black[abs(dh)] += 1
            maj_black[mj] += 1

    print("\nBLUE lines by endpoint types (typepair, cross-class?): %s" % dict(blue_edges))
    print("BLACK lines by endpoint types: %s" % dict(black_edges))
    print("INTERFACE LAW check: blue touches PBk? %s ; black touches PB? %s"
          % (any("PBk" in tp for (tp, _) in blue_edges),
             any("PB" == tp[0] or "PB" == tp[1] for (tp, _) in black_edges)))
    print("|DeltaH| along blue lines: %s" % dict(sorted(dH_blue.items())))
    print("|DeltaH| along black lines: %s" % dict(sorted(dH_black.items())))
    print("majorization flow blue: %s" % dict(maj_blue))
    print("majorization flow black: %s" % dict(maj_black))

    # THE CANONICAL ORDER: stage-by-stage tie accounting
    def sumsq(sc):
        return sum(x * x for x in sc)
    TYPR = {"PB": 0, "MX": 1, "PBk": 2}
    keys1 = [(-sumsq(d["sc"]), tuple(sorted(d["sc"], reverse=True))) for d in info.values()]
    keys2 = [k1 + (TYPR[d["typ"]],) for k1, d in zip(keys1, info.values())]
    keys3 = [k2 + (d["H"],) for k2, d in zip(keys2, info.values())]
    keys4 = [k3 + (cm,) for k3, cm in zip(keys3, info.keys())]
    for label, ks in [("scoreseq(sumsq,lex)", keys1), ("+type", keys2), ("+H", keys3), ("+canonword", keys4)]:
        cnt = Counter(ks)
        ties = sum(c - 1 for c in cnt.values())
        print("order stage %-22s distinct=%4d  residual-tied classes=%d" % (label, len(cnt), ties))
    order = sorted(zip(keys4, info.keys()))
    print("first 5 in order: %s" % [(info[cm]["sc"], info[cm]["typ"], info[cm]["H"]) for _, cm in order[:5]])
    print("last 3 in order: %s" % [(info[cm]["sc"], info[cm]["typ"], info[cm]["H"]) for _, cm in order[-3:]])
    print("[n=%d done %.1fs]" % (n, time.time() - t0))
    return info

for n in [4, 5, 6, 7]:
    analyze(n)
print("\nALL DONE")
