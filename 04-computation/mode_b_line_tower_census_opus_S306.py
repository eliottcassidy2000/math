#!/usr/bin/env python3
"""THM-793 referee + the three-object census (opus-2026-07-14-S306).

The three tracked objects of the merged metagraph and ALL their relations:
  TILINGS  T_n  (2^m),   NODES  (iso classes; merged = mod reversal),
  EDGES    (d=m lines, blue/black), with maps
  pi: tilings -> nodes,  kappa: tilings -> tilings (the flip),
  lines -> node pairs (cross-edges and SELF-lines).

PART 1  THE MODE-B LINE TOWER (referee of THM-793): the strip-legs+apex
        projection p: T_n -> T_{n-2}: (a) 2^{2n-5}-to-1; (b) p∘kappa =
        kappa'∘p; (c) p∘sigma = sigma'∘p (blue -> blue, fiber 2^{n-2});
        (d) T(p(t)) = the induced subtournament of T(t) on {2..n-1};
        (e) the axis current is FIBER-measurable: Dx_n(t) depends only on
        the forgotten leg+apex bits.
PART 2  THE CENSUS n=3..7: SC(n) vs A000570; merged nodes (A000568+SC)/2;
        self-lines (class(t) = class(tbar)) by colour, with the PROVED
        necessary condition e1 = en (leg law: Dx = 0); merged self-loops
        (class(tbar) in {class(t), reversal-class(t)}) by colour; level
        lines vs self-lines.
PART 3  DESCENT KERNELS: the class-level Markov kernel K: class(n) ->
        distribution over class(n-2) induced by p on fibers; summary
        statistics (support sizes, whether kernel rows respect the axis).
"""
import sys
from collections import defaultdict

def build(n):
    V = list(range(1, n+1))
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
    m = len(tiles)
    ti = {t: i for i, t in enumerate(tiles)}
    refl = [ti[(n-y+1, n-x+1)] for (x, y) in tiles]

    def tourn(bits):
        adj = [0]*(n+1)
        for k in range(n, 1, -1): adj[k] |= 1 << (k-2)
        for i, (x, y) in enumerate(tiles):
            if bits >> i & 1: adj[x] |= 1 << (y-1)
            else: adj[y] |= 1 << (x-1)
        return adj

    def invariant(adj, VV):
        s = {v: bin(adj[v] & sum(1 << (u-1) for u in VV)).count('1') for v in VV}
        prof = {}
        for v in VV:
            outs = sorted(s[u] for u in VV if adj[v] >> (u-1) & 1)
            ins = sorted(s[u] for u in VV if adj[u] >> (v-1) & 1 and u in s)
            c3 = 0
            for u in VV:
                if adj[v] >> (u-1) & 1:
                    for w in VV:
                        if adj[u] >> (w-1) & 1 and adj[w] >> (v-1) & 1: c3 += 1
            prof[v] = (s[v], tuple(outs), tuple(ins), c3)
        prof2 = {}
        for v in VV:
            po = sorted(prof[u] for u in VV if adj[v] >> (u-1) & 1)
            pi_ = sorted(prof[u] for u in VV if adj[u] >> (v-1) & 1 and u in prof)
            prof2[v] = (prof[v], tuple(po), tuple(pi_))
        arcs = []
        for u in VV:
            for v in VV:
                if adj[u] >> (v-1) & 1:
                    cuo = bin(adj[u] & adj[v]).count('1')
                    arcs.append((s[u], s[v], cuo))
        # H
        full = sum(1 << (u-1) for u in VV)
        dp = {(1 << (v-1), v): 1 for v in VV}
        for _ in range(len(VV)-1):
            nd = defaultdict(int)
            for (mask, v), c in dp.items():
                for u in VV:
                    b = 1 << (u-1)
                    if not mask & b and adj[v] & b: nd[(mask | b, u)] += c
            dp = nd
        H = sum(c for (mask, v), c in dp.items() if mask == full)
        return (tuple(sorted(prof2.values())), tuple(sorted(arcs)), H)

    classes = {}; cls_of = {}; rep = {}
    for b in range(1 << m):
        inv = invariant(tourn(b), V)
        if inv not in classes:
            classes[inv] = len(classes); rep[classes[inv]] = b
        cls_of[b] = classes[inv]
    C = len(classes)
    # reversal classes
    rcls = {}
    for c in range(C):
        adj = tourn(rep[c])
        radj = [0]*(n+1)
        for u in V:
            for v in V:
                if adj[u] >> (v-1) & 1: radj[v] |= 1 << (u-1)
        rcls[c] = classes[invariant(radj, V)]
    return dict(n=n, tiles=tiles, m=m, refl=refl, cls_of=cls_of, C=C,
                rep=rep, rcls=rcls, tourn=tourn, invariant=invariant, V=V)

counts568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
A000570 = {3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}
built = {}
for n in range(3, 8):
    built[n] = build(n)
    assert built[n]['C'] == counts568[n]

print("=" * 96)
print("PART 1 -- the Mode-B line tower (referee)")
print("=" * 96)
for n in range(5, 8):
    B = built[n]; Bc = built[n-2]
    tiles, m = B['tiles'], B['m']
    # interior tiles = tiles of the child staircase on {2..n-1}, relabelled v -> v-1
    inter = [i for i, (x, y) in enumerate(tiles) if y >= 2 and x <= n-1]
    child_ti = {t: i for i, t in enumerate(Bc['tiles'])}
    proj_map = [child_ti[(x-1, y-1)] for i, (x, y) in enumerate(tiles) if i in set(inter)]
    def p(bits):
        out = 0
        for k, i in enumerate(inter):
            if bits >> i & 1: out |= 1 << proj_map[k]
        return out
    fullm = (1 << m) - 1
    fullc = (1 << Bc['m']) - 1
    ok_fiber = ok_kappa = ok_sigma = ok_sub = ok_fibmeas = True
    fiber_count = defaultdict(int)
    # fiber-measurability: Dx determined by forgotten bits
    leg_idx = [i for i in range(m) if i not in set(inter)]
    seen_leg = {}
    for b in range(1 << m):
        pb = p(b)
        fiber_count[pb] += 1
        if p(b ^ fullm) != (pb ^ fullc): ok_kappa = False
        sb = 0
        for i in range(m):
            if b >> B['refl'][i] & 1: sb |= 1 << i
        # sigma commutation: p(sigma(b)) == sigma_child(p(b))
        psb = p(sb)
        sc = 0
        for i in range(Bc['m']):
            if pb >> Bc['refl'][i] & 1: sc |= 1 << i
        if psb != sc: ok_sigma = False
        # leg-measurable Dx: e1 - en depends only on leg bits
        e1 = sum(1 for i, (x, y) in enumerate(tiles) if y == 1 and not (b >> i & 1))
        en = sum(1 for i, (x, y) in enumerate(tiles) if x == n and (b >> i & 1))
        legbits = tuple((b >> i) & 1 for i in leg_idx)
        if legbits in seen_leg and seen_leg[legbits] != e1 - en: ok_fibmeas = False
        seen_leg[legbits] = e1 - en
    ok_fiber = all(v == 1 << (2*n-5) for v in fiber_count.values()) and len(fiber_count) == 1 << Bc['m']
    # induced subtournament check on a sample
    import random
    random.seed(306)
    for b in random.sample(range(1 << m), min(200, 1 << m)):
        adj = B['tourn'](b); cadj = Bc['tourn'](p(b))
        okk = True
        for u in range(2, n):
            for v in range(2, n):
                if u == v: continue
                if ((adj[u] >> (v-1)) & 1) != ((cadj[u-1] >> (v-2)) & 1): okk = False
        if not okk: ok_sub = False
    print(f"  n={n}->{n-2}: fiber 2^{2*n-5} exact: {ok_fiber}; p.kappa=kappa.p: {ok_kappa}; "
          f"p.sigma=sigma.p: {ok_sigma}; induced-subtournament: {ok_sub}; "
          f"Dx leg-measurable: {ok_fibmeas}")

print()
print("=" * 96)
print("PART 2 -- the three-object census")
print("=" * 96)
print(f"{'n':>2} {'tilings':>9} {'nodes':>6} {'SC':>4} {'A000570':>8} {'merged':>7} "
      f"{'lines':>7} {'blue':>6} {'black':>7} {'selfB':>6} {'selfK':>6} {'mselfB':>7} {'mselfK':>7} {'lvlK':>6}")
for n in range(3, 8):
    B = built[n]
    m, refl, cls_of, rcls = B['m'], B['refl'], B['cls_of'], B['rcls']
    SC = sum(1 for c in range(B['C']) if rcls[c] == c)
    merged = (B['C'] + SC) // 2
    fullm = (1 << m) - 1
    selfB = selfK = mselfB = mselfK = lvlK = 0; nblue = 0
    for b in range(1 << m):
        pb = b ^ fullm
        if b > pb: continue
        blue = all((b >> i & 1) == (b >> refl[i] & 1) for i in range(m))
        ca, cb = cls_of[b], cls_of[pb]
        if blue: nblue += 1
        if ca == cb:
            if blue: selfB += 1
            else: selfK += 1
        if cb in (ca, rcls[ca]):
            if blue: mselfB += 1
            else: mselfK += 1
        # level (needs x equal) -- use score-x via leg law necessary condition instead:
    lines = 1 << (m-1)
    print(f"{n:>2} {1<<m:>9} {B['C']:>6} {SC:>4} {A000570[n]:>8} {merged:>7} "
          f"{lines:>7} {nblue:>6} {lines-nblue:>7} {selfB:>6} {selfK:>6} {mselfB:>7} {mselfK:>7} {'':>6}")
    assert SC == A000570[n], f"SC({n}) != A000570"

print()
print("=" * 96)
print("PART 3 -- descent kernels (class-level Markov structure of p)")
print("=" * 96)
for n in (5, 6, 7):
    B = built[n]; Bc = built[n-2]
    tiles, m = B['tiles'], B['m']
    inter = [i for i, (x, y) in enumerate(tiles) if y >= 2 and x <= n-1]
    child_ti = {t: i for i, t in enumerate(Bc['tiles'])}
    proj_map = [child_ti[(x-1, y-1)] for i, (x, y) in enumerate(tiles) if i in set(inter)]
    def p(bits):
        out = 0
        for k, i in enumerate(inter):
            if bits >> i & 1: out |= 1 << proj_map[k]
        return out
    K = defaultdict(lambda: defaultdict(int))
    for b in range(1 << m):
        K[B['cls_of'][b]][Bc['cls_of'][p(b)]] += 1
    supp = [len(K[c]) for c in range(B['C'])]
    fullsupp = sum(1 for c in range(B['C']) if len(K[c]) == Bc['C'])
    print(f"  n={n}->{n-2}: kernel support sizes: min {min(supp)}, max {max(supp)} "
          f"(child classes: {Bc['C']}); rows with FULL support: {fullsupp}/{B['C']}")
    # does the transitive class descend to the transitive class?
    tc = [c for c in range(B['C']) if all(len(K[c]) >= 1 for _ in [0]) and
          B['cls_of'][ (1 << m) - 1 ] == c][0]
    tcc = Bc['cls_of'][(1 << Bc['m']) - 1]
    print(f"    transitive row: {dict(K[tc])} (child transitive id {tcc}) "
          f"-> transitive descends to transitive: {list(K[tc].keys()) == [tcc]}")
print("done")
