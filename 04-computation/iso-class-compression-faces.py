#!/usr/bin/env python3
"""
iso-class-compression-faces.py

Verifies and extends the "tiling-cube compression" observation:

  At n=4, the m = C(n-1,2) = 3 free tiles of the tiling model live on the cube
  Q_3 (8 tilings). The iso-class map f : Q_3 -> {T,+,-,S} (A000568(4)=4 classes)
  is NOT injective. But there is a 2-dimensional FACE of Q_3 on which f is a
  bijection. On that face the 4 classes carry the Klein four-group V_4 = (Z2)^2
  structure (XOR of two bits = "source-destroyed" and "sink-destroyed").
  Complementation acts on V_4 as the coordinate swap; SC classes = the diagonal
  {00,11}={T,S}; the NS complement-pair = the antidiagonal {10,01}={+,-}.

This script:
  (a) enumerates tilings on a fixed Hamiltonian base path, computes exact iso
      classes by brute-force canonical form over S_n, and verifies class counts
      against A000568 for n=4,5,6;
  (b) for n=4 prints both Cayley-style tables (3 generators a,b,c vs 2 gens x,y),
      checks the V_4 / XOR law, and verifies complement = coordinate swap + SC/NS;
  (c) for n=4,5,6 searches for the MINIMUM-dimension face of Q_m that surjects
      onto every iso class ("compression budget"), and compares to the
      information-theoretic bound ceil(log2 A000568) and to log2(n!/2^{n-1}).
"""

import itertools
from math import comb, log2, factorial

A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}

def pairs(n):
    """Canonical ordering of unordered vertex pairs (i,j), i<j."""
    return [(i, j) for i in range(n) for j in range(i+1, n)]

def base_path_pairs(n):
    """Consecutive pairs (k,k+1) = the fixed Hamiltonian base path arcs."""
    return [(k, k+1) for k in range(n-1)]

def tile_pairs(n):
    """Non-consecutive pairs = the m = C(n-1,2) flippable tiles."""
    bp = set(base_path_pairs(n))
    return [p for p in pairs(n) if p not in bp]

def tiling_bits(n, flips):
    """
    Bitmask over canonical pairs. Bit b(i,j)=1 means i beats j (i<j).
    Transitive base (higher beats lower) = all zeros. Each flipped tile sets its bit.
    """
    idx = {p: k for k, p in enumerate(pairs(n))}
    bits = 0
    for p in flips:
        bits |= (1 << idx[p])
    return bits

def scores(n, bits):
    P = pairs(n); idx = {p:k for k,p in enumerate(P)}
    s = [0]*n
    for (i,j) in P:
        if bits >> idx[(i,j)] & 1:   # i beats j
            s[i]+=1
        else:                         # j beats i
            s[j]+=1
    return tuple(sorted(s))

# ---- canonical form over all permutations (exact iso class) ----
def perm_tables(n):
    """For each permutation, precompute (target_bit_index, invert?) per canonical pair."""
    P = pairs(n); idx = {p:k for k,p in enumerate(P)}
    tables = []
    for perm in itertools.permutations(range(n)):
        row = []
        for (i,j) in P:
            a,b = perm[i], perm[j]
            if a < b:
                row.append((idx[(a,b)], False))
            else:
                row.append((idx[(b,a)], True))
        tables.append(row)
    return tables

def canonical(bits, tables, m_pairs):
    best = None
    for row in tables:
        v = 0
        for q,(tgt,inv) in enumerate(row):
            bit = (bits >> tgt) & 1
            if inv: bit ^= 1
            v |= (bit << q)
        if best is None or v < best:
            best = v
    return best

def complement_bits(bits, n):
    full = (1 << comb(n,2)) - 1
    return bits ^ full   # reverse every arc

def analyze(n, do_face_search=True):
    print(f"\n{'='*70}\n n = {n}   (m = C(n-1,2) = {comb(n-1,2)} tiles,  2^m = {2**comb(n-1,2)} tilings)\n{'='*70}")
    tiles = tile_pairs(n); m = len(tiles)
    tables = perm_tables(n)
    # iso class (canonical form) for every tiling
    class_of = {}        # flip-bitvector (over tiles) -> canonical int
    canon_set = {}       # canonical int -> small label index
    for mask in range(2**m):
        flips = [tiles[k] for k in range(m) if (mask>>k)&1]
        bits = tiling_bits(n, flips)
        c = canonical(bits, tables, comb(n,2))
        class_of[mask] = c
        if c not in canon_set:
            canon_set[c] = len(canon_set)
    nclasses = len(canon_set)
    print(f" iso classes found: {nclasses}   (A000568({n}) = {A000568[n]})   "
          f"{'OK' if nclasses==A000568[n] else 'MISMATCH!'}")
    # fiber sizes
    from collections import Counter
    fib = Counter(class_of.values())
    print(f" fiber sizes (tilings per class), sorted: {sorted(fib.values(), reverse=True)}")

    if n == 4:
        n4_tables(n, tiles, class_of, canon_set, tables)

    if do_face_search:
        face_search(n, tiles, class_of, canon_set, nclasses)

def n4_tables(n, tiles, class_of, canon_set, tables):
    """Reproduce the user's two tables and verify Klein/complement structure."""
    print("\n --- n=4 detailed structure ---")
    # identify which tile-flip gives which class by score sequence
    name = {(0,1,2,3):'T', (0,2,2,2):'+', (1,1,1,3):'-', (1,1,2,2):'S'}
    # map canonical int -> name via a representative's score sequence
    cint_to_name = {}
    for mask, c in class_of.items():
        flips=[tiles[k] for k in range(len(tiles)) if (mask>>k)&1]
        ss = scores(n, tiling_bits(n,flips))
        cint_to_name[c] = name[ss]
    # tiles in user's labeling: a=(4,2),b=(3,1),c=(4,1) in 1-indexed => 0-indexed (3,1),(2,0),(3,0)
    # find tile indices by their (i,j) 0-indexed identity
    def tidx(pair0):
        return tiles.index(pair0)
    a = tidx((1,3)); b = tidx((0,2)); c = tidx((0,3))   # 0-indexed (4,2)->(3,1)? careful
    # NOTE vertices are 0..3 representing labels 1..4; pair (i,j) i<j.
    # user arc a=(4,2): labels 4,2 -> 0-indexed (1,3). b=(3,1)->(0,2). c=(4,1)->(0,3).
    gens = {'a':a, 'b':b, 'c':c}
    def cls(mask):
        return cint_to_name[class_of[mask]]
    def flipmask(*tile_indices):
        m=0
        for t in tile_indices: m |= (1<<t)
        return m
    print("  single-tile flips:  flip a ->", cls(flipmask(a)),
          " flip b ->", cls(flipmask(b)), " flip c ->", cls(flipmask(c)),
          "  (expected +,-,S)")
    # Scheme 1 table over E,a,b,c  (entry = class of XOR of the two generators)
    labels1 = ['E','a','b','c']; gmask={'E':0,'a':flipmask(a),'b':flipmask(b),'c':flipmask(c)}
    print("\n  Scheme 1 (tiling model, generators a,b,c):")
    hdr = "    *  | " + "  ".join(f"{l}" for l in labels1)
    print(hdr)
    for r in labels1:
        row = [cls(gmask[r]^gmask[col]) for col in labels1]
        print(f"    {r}  | " + "  ".join(row))
    # Scheme 2 table over E,x,y  (x=a, y=b ; the c=0 face)
    labels2=['E','x','y']; g2={'E':0,'x':flipmask(a),'y':flipmask(b)}
    print("\n  Scheme 2 (compressed, x=a, y=b, tile c fixed=0):")
    print("    *  | " + "  ".join(labels2))
    for r in labels2:
        row=[cls(g2[r]^g2[col]) for col in labels2]
        print(f"    {r}  | " + "  ".join(row))
    # Klein / XOR check on the c=0 face: bijection 2 bits -> 4 classes?
    face = [0, flipmask(a), flipmask(b), flipmask(a,b)]
    facecls = [cls(mm) for mm in face]
    print("\n  c=0 face  (00,10,01,11) -> classes:", facecls,
          "  bijection:", len(set(facecls))==4)
    # partial score sequence of the frame (base path + tile c unflipped)
    frame_bits = tiling_bits(n, [])  # all tiles 0 incl c
    # only base path + c are "fixed"; partial scores counting base path + arc c only:
    # arc c=(4,1)->0-indexed (0,3); base path consecutive pairs.
    fixed_pairs = base_path_pairs(n) + [(0,3)]
    ps=[0]*n; idx={p:k for k,p in enumerate(pairs(n))}
    for (i,j) in fixed_pairs:
        if (frame_bits>>idx[(i,j)])&1: ps[i]+=1
        else: ps[j]+=1
    print("  partial score sequence of frame (base path + c):", tuple(sorted(ps)),
          "  (user said 0,1,1,2)")
    # complement action on the 4 classes
    print("\n  complement (reverse all arcs) action & SC/NS:")
    # need canonical of complement; recompute canonical with full tables
    reps={}
    for mask,c in class_of.items():
        reps.setdefault(c,mask)
    for c,nm in cint_to_name.items():
        mask=reps[c]; flips=[tiles[k] for k in range(len(tiles)) if (mask>>k)&1]
        bits=tiling_bits(n,flips); comp=complement_bits(bits,n)
        cc=canonical(comp,tables,comb(n,2))
        print(f"    {nm:2s} -> complement -> {cint_to_name[cc]:2s}   "
              f"{'SC (self-complementary)' if cc==c else 'NS'}")

def face_search(n, tiles, class_of, nclasses_canon_set=None, nclasses=None):
    """Min-dimension face of Q_m surjecting onto all iso classes."""
    m=len(tiles)
    classes_full=set(class_of.values())
    info_bound=0
    while (1<<info_bound) < len(classes_full): info_bound+=1
    print(f"\n --- compression / minimal covering face (n={n}) ---")
    print(f"   classes={len(classes_full)},  info bound ceil(log2)={info_bound},  "
          f"m={m},  redundancy budget log2(n!/2^(n-1))={log2(factorial(n)/2**(n-1)):.2f}")
    found=None
    for k in range(info_bound, m+1):
        ok_example=None
        for free in itertools.combinations(range(m), k):
            freeset=list(free)
            fixed=[i for i in range(m) if i not in free]
            # iterate over all fixings of the 'fixed' coords; group tilings
            # by fixed-pattern, check coverage of each coset
            # build cosets
            cover={}
            for mask in class_of:
                key=tuple((mask>>i)&1 for i in fixed)
                cover.setdefault(key,set()).add(class_of[mask])
            for key,cset in cover.items():
                if cset==classes_full:
                    ok_example=(freeset,fixed,key); break
            if ok_example: break
        if ok_example:
            found=(k,ok_example); break
    if found:
        k,(freeset,fixed,key)=found
        print(f"   MINIMUM covering face dimension = {k}  "
              f"(fix {m-k} tiles, free {k});  2^{k}={2**k} tilings cover all {len(classes_full)} classes")
        print(f"   tight to info bound: {k==info_bound}")
        print(f"   example: free tiles (indices) {freeset}, fixed tiles {fixed} at pattern {key}")
    else:
        print("   no covering face found (should not happen, k=m always works)")

if __name__=="__main__":
    for n in (4,5,6):
        analyze(n, do_face_search=True)
    # info-theoretic table for higher n (no enumeration)
    print(f"\n{'='*70}\n compression budget across n (info-theoretic)\n{'='*70}")
    print(f"  {'n':>2} {'classes A000568':>16} {'m=C(n-1,2)':>11} {'ceil log2':>10} {'budget=m-ceil':>13} {'log2(n!/2^(n-1))':>17}")
    for n in range(3,11):
        cl=A000568[n]; m=comb(n-1,2)
        ib=0
        while (1<<ib)<cl: ib+=1
        print(f"  {n:>2} {cl:>16} {m:>11} {ib:>10} {m-ib:>13} {log2(factorial(n)/2**(n-1)):>17.2f}")
