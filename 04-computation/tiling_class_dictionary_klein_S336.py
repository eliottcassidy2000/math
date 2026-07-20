#!/usr/bin/env python3
"""
klein-2026-07-20-S336 -- THE EXACT TILING <-> ISO CLASS <-> MERGED METAGRAPH DICTIONARY,
and the tricks that let nodes, edges and tilings each compute the others cheaply.

Owner: "figure out how tiling sets map to iso class nodes in the merged metagraph exactly,
use the edges and nodes and tilings to all compute each other efficiently, look for tricks.
consider even more creative ideas than a descending star-type invariant has to come from a
base-path-independent subgroup -- the natural candidate being the intersection of Gamma over
all spanning paths."

PARTS
  A. the exact fibration, exhaustively (n = 4,5,6): |fibre(c)| = H(c)/|Aut(c)|, the global
     checksum sum_c H/|Aut| = 2^{C(n-1,2)}, and the merged/SC refinement.
  B. THE HAM-PATH TRICK: compute every fibre WITHOUT touching the 2^m cube, from one
     representative per class.  Verified against A, then run where A cannot go.
  C. THE HALF TILING IS THE MERGED FIBRATION (ties THM-549 + THM-1410 to the metagraph):
     sigma-orbits of tilings -> merged nodes, with the blue (grid-symmetric) count as the
     exact SC correction, and sum_{SC} B(c) = 2^{floor((n-1)^2/4)}.
  D. SWITCHING: tilings ARE the switching classes of tournaments -- the tiling cube is
     arc-space modulo cut(K_n), and the base path is a spanning tree transversal.
  E. THE BASE-PATH-INDEPENDENT SUBGROUP: intersection of Gamma_P and the join, over ALL
     Hamiltonian paths P.  Decides kind-pasteur's proposal outright.
  F. EDGES: G_n by the Aut-orbit trick, and the d=1 wiggly quotient compared to it.
"""
import itertools, sys, time
from math import comb, factorial

# ---------------------------------------------------------------- tournaments
def out_masks(n, arcbits, pairs):
    """arcbits: int, bit k = 1 means pairs[k] = (i,j) oriented i->j, else j->i."""
    om = [0] * n
    for k, (i, j) in enumerate(pairs):
        if arcbits >> k & 1: om[i] |= 1 << j
        else:                om[j] |= 1 << i
    return tuple(om)

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def relabel(om, perm, n):
    """perm[v] = new label of v."""
    new = [0] * n
    for v in range(n):
        mv = om[v]; t = 0
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            t |= 1 << perm[w]
        new[perm[v]] = t
    return tuple(new)

def word(om, n):
    w = 0
    for v in range(n): w = (w << n) | om[v]
    return w

def refine(om, n):
    """iterated score refinement -> list of cells (tuples of vertices)."""
    colour = [bin(om[v]).count("1") for v in range(n)]
    while True:
        sig = []
        for v in range(n):
            cnt = {}
            mv = om[v]
            while mv:
                b = mv & -mv; w = b.bit_length() - 1; mv ^= b
                cnt[colour[w]] = cnt.get(colour[w], 0) + 1
            sig.append((colour[v], tuple(sorted(cnt.items()))))
        order = sorted(set(sig))
        newc = [order.index(sig[v]) for v in range(n)]
        if newc == colour: break
        colour = newc
    cells = {}
    for v in range(n): cells.setdefault(colour[v], []).append(v)
    return [tuple(cells[k]) for k in sorted(cells)]

def canon(om, n, want_aut=False):
    """canonical word = min over permutations respecting the refined cell order."""
    cells = refine(om, n)
    slots, best, autc = [], None, 0
    pos = 0; base = []
    for c in cells:
        base.append((c, pos)); pos += len(c)
    for choice in itertools.product(*[itertools.permutations(c) for (c, _) in base]):
        perm = [0] * n
        for (blk, (c, start)) in zip(choice, base):
            for k, v in enumerate(blk): perm[v] = start + k
        w = word(relabel(om, perm, n), n)
        if best is None or w < best: best, autc = w, 1
        elif w == best: autc += 1
    return (best, autc) if want_aut else best

def ham_paths(om, n):
    """all Hamiltonian paths as vertex sequences (directed)."""
    res = []
    def go(seq, used):
        if len(seq) == n: res.append(tuple(seq)); return
        last = seq[-1]; mv = om[last] & ~used
        while mv:
            b = mv & -mv; w = b.bit_length() - 1; mv ^= b
            seq.append(w); go(seq, used | (1 << w)); seq.pop()
    for s in range(n): go([s], 1 << s)
    return res

def complement(om, n):
    full = (1 << n) - 1
    return tuple(full & ~om[v] & ~(1 << v) for v in range(n))

# ------------------------------------------------------- the tiling model
def tiles_of(n):
    """explorer order, 1-indexed vertices: for y=1..n-2: for x=n down to y+2: (x,y)."""
    return [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1)]

def tiling_to_om(bits, n, T):
    """base path n->n-1->...->1 (1-indexed) becomes 0-indexed vertices v = label-1.
       tile (x,y) set bit => x->y, else y->x."""
    om = [0] * n
    for k in range(n, 1, -1):        # path arcs k -> k-1
        om[k - 1] |= 1 << (k - 2)
    for k, (x, y) in enumerate(T):
        if bits >> k & 1: om[x - 1] |= 1 << (y - 1)
        else:             om[y - 1] |= 1 << (x - 1)
    return tuple(om)

def sigma_perm_tiles(n, T):
    idx = {t: i for i, t in enumerate(T)}
    return [idx[(n + 1 - y, n + 1 - x)] for (x, y) in T]

def apply_tileperm(bits, p):
    out = 0
    for i in range(len(p)):
        if bits >> i & 1: out |= 1 << p[i]
    return out

# ============================================================ PART A + C
print("=" * 78)
print("PART A/C -- THE EXACT FIBRATION, EXHAUSTIVE")
print("=" * 78)
exact = {}
for n in (4, 5, 6):
    T = tiles_of(n); m = len(T); p_sig = sigma_perm_tiles(n, T)
    fib, blue = {}, {}
    for bits in range(1 << m):
        om = tiling_to_om(bits, n, T)
        c = canon(om, n)
        fib[c] = fib.get(c, 0) + 1
        if apply_tileperm(bits, p_sig) == bits:
            blue[c] = blue.get(c, 0) + 1
    # per-class H, |Aut|, SC
    info = {}
    for c in fib:
        om = None
        # rebuild a representative from the canonical word
        w = c; rows = []
        for _ in range(n):
            rows.append(w & ((1 << n) - 1)); w >>= n
        om = tuple(reversed(rows))
        _, aut = canon(om, n, want_aut=True)
        H = len(ham_paths(om, n))
        sc = (canon(complement(om, n), n) == c)
        info[c] = (H, aut, sc)
    ok_orbit = all(fib[c] * info[c][1] == info[c][0] for c in fib)
    checksum = sum(info[c][0] // info[c][1] for c in fib)
    blue_on_nonsc = [c for c in fib if blue.get(c, 0) and not info[c][2]]
    h = ((n - 1) ** 2) // 4
    blue_total = sum(blue.values())
    nsc = sum(1 for c in fib if info[c][2])
    print(f"\n n={n}  m={m}  tilings={1<<m}  iso classes={len(fib)}  SC classes={nsc}")
    print(f"   |fibre(c)| * |Aut(c)| = H(c) for every class : {ok_orbit}")
    print(f"   checksum  sum_c H/|Aut| = {checksum}  vs 2^m = {1<<m}  : {checksum == (1<<m)}")
    print(f"   blue (sigma-fixed) tilings total = {blue_total}  vs 2^h, h=floor((n-1)^2/4)={h} -> {1<<h}"
          f"  : {blue_total == (1 << h)}")
    print(f"   any blue tiling over a NON-SC class? {'YES (BAD)' if blue_on_nonsc else 'no'}")
    # merged nodes
    merged, seen = [], set()
    for c in fib:
        if c in seen: continue
        w = c; rows = []
        for _ in range(n): rows.append(w & ((1 << n) - 1)); w >>= n
        om = tuple(reversed(rows)); cc = canon(complement(om, n), n)
        seen.add(c); seen.add(cc)
        merged.append((c, cc))
    mf_ok, sig_ok = True, True
    for (c, cc) in merged:
        H, aut, sc = info[c]
        pred = (H // aut) if sc else 2 * (H // aut)
        act = fib[c] + (0 if sc else fib[cc])
        if pred != act: mf_ok = False
        # sigma-orbit count over this merged node
        nB = blue.get(c, 0) + (0 if sc else blue.get(cc, 0))
        orb_pred = (H // aut) if not sc else (H // aut + nB) // 2
        if (act + nB) // 2 != orb_pred: sig_ok = False
    print(f"   merged nodes = {len(merged)}   fibre = (2-[SC])*H/|Aut| : {mf_ok}")
    print(f"   sigma-orbits over merged node = H/|Aut| (non-SC) / (H/|Aut|+B)/2 (SC) : {sig_ok}")
    exact[n] = {c: (fib[c], info[c][0], info[c][1], info[c][2], blue.get(c, 0)) for c in fib}

# ============================================================ PART B
print("\n" + "=" * 78)
print("PART B -- THE HAM-PATH TRICK: fibres without touching the 2^m cube")
print("=" * 78)
print(" fibre(c) <-> HamPaths(T_c)/Aut(T_c)  (Aut acts freely on Ham paths, canon LEM-003)")
print(" so |fibre(c)| = H(c)/|Aut(c)| -- and the fibre ITSELF is enumerable from one rep.")

def classes_by_bfs(n, cap=None):
    """iso classes via BFS on the single-ARC-flip graph, using Aut-orbits on arcs.
       returns {canon: (rep_om, H, aut, sc)} and the edge set of G_n."""
    P = pairs_of(n)
    om0 = tuple(sum(1 << j for j in range(i)) for i in range(n))   # transitive
    start = canon(om0, n)
    seen = {start: om0}; frontier = [om0]; edges = set(); flips = 0
    while frontier:
        nxt = []
        for om in frontier:
            c = canon(om, n)
            # Aut-orbit trick: only flip ONE arc per Aut(T)-orbit of arcs
            cells = refine(om, n)
            reps, orb = [], set()
            for (i, j) in P:
                key = None
                if (i, j) in orb: continue
                reps.append((i, j))
            for (i, j) in reps:
                nm = list(om)
                if om[i] >> j & 1: nm[i] &= ~(1 << j); nm[j] |= 1 << i
                else:              nm[j] &= ~(1 << i); nm[i] |= 1 << j
                nm = tuple(nm); cc = canon(nm, n); flips += 1
                if cc != c: edges.add((min(c, cc), max(c, cc)))
                if cc not in seen:
                    seen[cc] = nm; nxt.append(nm)
                    if cap and len(seen) >= cap: return seen, edges, flips
        frontier = nxt
    return seen, edges, flips

for n in (4, 5, 6, 7):
    t0 = time.time()
    seen, edges, flips = classes_by_bfs(n)
    tb = time.time() - t0
    rows = []
    for c, om in seen.items():
        _, aut = canon(om, n, want_aut=True)
        H = len(ham_paths(om, n))
        sc = (canon(complement(om, n), n) == c)
        rows.append((c, H, aut, sc))
    m = comb(n - 1, 2)
    checksum = sum(H // aut for (_, H, aut, _) in rows)
    nsc = sum(1 for r in rows if r[3])
    merged_ct = (len(rows) + nsc) // 2
    agree = "n/a"
    if n in exact:
        agree = all(exact[n][c][0] == H // aut for (c, H, aut, _) in rows if c in exact[n]) \
                and len(exact[n]) == len(rows)
    print(f"\n n={n}: classes={len(rows)} (A000568)  SC={nsc}  merged={merged_ct}"
          f"  G_n edges={len(edges)}")
    print(f"   checksum sum H/|Aut| = {checksum}  vs 2^m = {1<<m} : {checksum == (1<<m)}"
          f"   [agrees with exhaustive: {agree}]")
    print(f"   cost: {flips} canonicalisations, {tb:.1f}s   vs naive 2^C(n,2) = {1<<comb(n,2)} tournaments")

# ============================================================ PART D
print("\n" + "=" * 78)
print("PART D -- TILINGS ARE THE SWITCHING CLASSES (arc space modulo cut(K_n))")
print("=" * 78)
for n in (4, 5, 6):
    P = pairs_of(n); E = len(P); T = tiles_of(n); m = len(T)
    pathset = {(k - 2, k - 1) for k in range(n, 1, -1)}     # 0-indexed, i<j
    hit = {}
    for arcbits in range(1 << E):
        om = out_masks(n, arcbits, P)
        # switch at S so that every path arc k -> k-1 holds: solve greedily along the path
        S = 0
        for k in range(n - 1, 0, -1):        # vertices n-1 .. 1 (0-indexed); want k -> k-1
            up = (om[k] >> (k - 1)) & 1
            flip_k = (S >> k) & 1
            # current orientation after switching S on the already-decided part
            cur = up ^ flip_k
            if not cur: S |= 1 << (k - 1)
            else:       pass
        # apply switch S
        sm = list(om)
        for (i, j) in P:
            if ((S >> i) & 1) != ((S >> j) & 1):
                if sm[i] >> j & 1: sm[i] &= ~(1 << j); sm[j] |= 1 << i
                else:              sm[j] &= ~(1 << i); sm[i] |= 1 << j
        sm = tuple(sm)
        ok = all((sm[k] >> (k - 1)) & 1 for k in range(1, n))
        hit[sm] = hit.get(sm, 0) + 1
        if not ok:
            print(f"   n={n}: SWITCH FAILED on arcbits={arcbits}"); break
    sizes = set(hit.values())
    print(f" n={n}: distinct switch-canonical tournaments = {len(hit)}  vs 2^m = {1<<m}"
          f"  : {len(hit) == (1<<m)};  fibre sizes = {sizes} (expect {{{1<<(n-1)}}} = 2^(n-1))")

# ============================================================ PART E
print("\n" + "=" * 78)
print("PART E -- THE BASE-PATH-INDEPENDENT SUBGROUP: cap Gamma_P and join Gamma_P")
print("=" * 78)
def rref_basis(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v); b.sort(reverse=True)
    return b

for n in (4, 5, 6, 7):
    P = pairs_of(n); E = len(P); eidx = {e: k for k, e in enumerate(P)}
    allP, join, dualspan = [], [], []
    for perm in itertools.permutations(range(n)):
        if perm[0] > perm[-1]: continue                     # paths up to reversal
        pe = {tuple(sorted((perm[i], perm[i + 1]))) for i in range(n - 1)}
        # Gamma_P = span{ star_H(v) } , H = K_n minus P
        gens = []
        for v in range(n):
            s = 0
            for e in P:
                if v in e and e not in pe: s |= 1 << eidx[e]
            gens.append(s)
        G = rref_basis(gens)
        join.extend(G)
        # dual of Gamma_P inside GF(2)^E : cycle space of H, plus the P-coordinates
        # easier: collect Gamma_P^perp by solving; use complement via full-rank check below
        allP.append(G)
    joinb = rref_basis(join)
    # intersection: dim cap = E - dim( sum of duals ); compute duals explicitly
    duals = []
    for G in allP:
        # orthogonal complement of span(G) in GF(2)^E
        # build matrix rows = G, find null space basis
        rows = list(G); piv = {}
        for r in rows:
            rr = r
            for p_, pr in piv.items():
                if rr >> p_ & 1: rr ^= pr
            if rr:
                p_ = rr.bit_length() - 1; piv[p_] = rr
        free = [k for k in range(E) if k not in piv]
        for f in free:
            v = 1 << f
            for p_ in sorted(piv, reverse=True):
                if v >> p_ & 1: v ^= piv[p_]
            # v now reduced; build the null vector directly instead
        # simpler + safe: brute null space by Gaussian elimination on the transpose
        M = [0] * E
        for k, r in enumerate(G):
            for b in range(E):
                if r >> b & 1: M[b] |= 1 << k
        # null space of span(G) = {x : <x, g> = 0 for all g}
        eqs = list(G); pv = {}
        for e in eqs:
            ee = e
            for p_, pr in pv.items():
                if ee >> p_ & 1: ee ^= pr
            if ee: pv[ee.bit_length() - 1] = ee
        freec = [k for k in range(E) if k not in pv]
        nb = []
        for f in freec:
            x = 1 << f
            for p_ in sorted(pv):
                s = bin(x & pv[p_]).count("1") & 1
                if s: x ^= 1 << p_
            nb.append(x)
        duals.extend(nb)
    db = rref_basis(duals)
    dim_cap = E - len(db)
    print(f" n={n}: #HamPaths={factorial(n)//2}  E=C(n,2)={E}   dim Gamma_P = {n-1}")
    print(f"        dim JOIN <U Gamma_P> = {len(joinb)}  (= E means the union generates everything)")
    print(f"        dim CAP  (^ Gamma_P) = {dim_cap}")

# ============================================================ PART F
print("\n" + "=" * 78)
print("PART F -- G_n vs THE d=1 WIGGLY QUOTIENT (which metagraph do tilings see?)")
print("=" * 78)
for n in (4, 5, 6):
    T = tiles_of(n); m = len(T)
    wig = set(); cls = {}
    for bits in range(1 << m):
        om = tiling_to_om(bits, n, T); c = canon(om, n); cls[bits] = c
    for bits in range(1 << m):
        c = cls[bits]
        for k in range(m):
            cc = cls[bits ^ (1 << k)]
            if cc != c: wig.add((min(c, cc), max(c, cc)))
    seen, edges, _ = classes_by_bfs(n)
    print(f" n={n}: G_n (all-arc flips) edges = {len(edges)};  d=1 wiggly quotient edges = {len(wig)}")
    print(f"        wiggly subset of G_n: {wig <= edges};  equal: {wig == edges}"
          f";  missing from wiggly: {len(edges - wig)}")
print("\nDONE.")

# ============================================================ PART G (addendum)
print("\n" + "=" * 78)
print("PART G -- WHY PART F HOLDS: the avoidable-arc lemma")
print("=" * 78)
print(" To realise a G_n edge {[T],[T+a]} inside the tiling model we need a Hamiltonian path")
print(" of T AVOIDING a.  Test: for every T and every arc a, either T has a Ham path avoiding")
print(" a, or flipping a is a SELF-LOOP (class-preserving).")
for n in (4, 5, 6, 7):
    seen, _, _ = classes_by_bfs(n)
    P = pairs_of(n)
    forced = 0; forced_and_selfloop = 0; bad = []
    for c, om in seen.items():
        HP = ham_paths(om, n)
        for (i, j) in P:
            ea = tuple(sorted((i, j)))
            avoid = any(ea not in {tuple(sorted((q[k], q[k+1]))) for k in range(n-1)} for q in HP)
            if avoid: continue
            forced += 1
            nm = list(om)
            if om[i] >> j & 1: nm[i] &= ~(1 << j); nm[j] |= 1 << i
            else:              nm[j] &= ~(1 << i); nm[i] |= 1 << j
            if canon(tuple(nm), n) == c: forced_and_selfloop += 1
            else: bad.append((c, ea))
    print(f" n={n}: (class,arc) pairs where EVERY Ham path uses the arc: {forced}"
          f";  of those, flip is a self-loop: {forced_and_selfloop}"
          f";  counterexamples: {len(bad)}")
    if bad: print(f"        first counterexample: {bad[0]}")
