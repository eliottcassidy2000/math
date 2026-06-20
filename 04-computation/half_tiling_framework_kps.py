#!/usr/bin/env python3
"""
half_tiling_framework_kps.py
================================================================
The HALF-TILING framework (kind-pasteur-2026-06-20).

USER INSIGHT (verbatim seed):
  "Taking the mirror image of a tiling over the line y=x is equivalent to
   reversing all arcs in the tournament (including those in the fixed path).
   Thus we can keep tiles along the y=x line and on one side, discarding the
   other side -> HALF-TILINGS.  Sizes (from a tournament of size 2):
       0,1,2,4,6,9,12,16,20,25,30  (n = 2,3,4,...)
   EVEN n half-tiling = A+B-C   (A,B half of size n-1; C tiling of size n-2)
   ODD  n half-tiling = A+B-C+D-E-F+G (A,B ~ n-1; C,D ~ n-2; E,F ~ n-3; G ~ n-4)
   producing corners A,D,B; edges A+D-E, B+D-F, A+B-C; center A+B+G-E-F."

This script makes that framework precise and verifies every claim against the
canonical tiling definitions (01-canon/definitions.md):
  Base path P_0 = n -> n-1 -> ... -> 1.
  Tile (a,b), a >= b+2 ; bit 0 => a->b (forward), bit 1 => b->a (backward).
  m = C(n-1,2) tiles.  Pin grid (r,c): r=a-b-1, c=b.
  isGridSym (CLAUDE.md): invariant under (x,y)->(n+1-y, n+1-x).

We establish:
  (1) The converse-with-relabel map is a PURE coordinate involution
        rho:(a,b) -> (n+1-b, n+1-a)  -- NO GF(2) bit flip.   [structural thm]
  (2) rho has d = floor((n-1)/2) fixed cells (on anti-diagonal a+b=n+1) and
        (m-d)/2 transposed pairs; #orbits = half(n) = floor((n-1)^2/4).
  (3) #(grid-symmetric tilings) = 2^half(n); grid-sym fraction = 2^((d-m)/2)
        reproducing canon exponents 0,-1,-2,-4,-6,-9 for n=3..8.
  (4) Recursions:  even  h(n)=2h(n-1)-h(n-2)
                   odd   h(n)=2h(n-1)-2h(n-3)+h(n-4)
        == the user's A+B-C and A+B-C+D-E-F+G signed decompositions.
  (5) Fundamental-domain shape (odd -> k^2 'gnomon square', even -> k(k-1)).
  (6) grid-sym tilings <=> phi-self-converse tournaments (phi = i->n+1-i),
        and their H-spectrum (SC maximizers, Paley).
"""

import itertools
from math import comb, floor

# ----------------------------------------------------------------------
# Tiling / tournament machinery (canonical conventions)
# ----------------------------------------------------------------------

def tiles(n):
    """Tiles (a,b) with a>=b+2, in the explorer's enumeration order:
       for y=1..n-2: for x=n down to y+2: tile (x,y)."""
    out = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            out.append((x, y))
    return out

def rc(a, b):
    """pin-grid coordinates: r = a-b-1, c = b."""
    return (a-b-1, b)

def rho(a, b, n):
    """converse-with-relabel reflection on a tile: (a,b)->(n+1-b, n+1-a).
       Returns the image tile with upper vertex first."""
    u, v = n+1-b, n+1-a          # n+1-b > n+1-a since b<a
    return (u, v)

def tiling_to_tournament(t, tile_list, n):
    """t: dict tile->bit. Return adjacency: T[(u,v)]=1 if u->v."""
    T = {}
    # base path arcs k->k-1
    for k in range(n, 1, -1):
        T[(k, k-1)] = 1
        T[(k-1, k)] = 0
    for (a, b) in tile_list:
        bit = t[(a, b)]
        if bit == 0:               # forward a->b
            T[(a, b)] = 1; T[(b, a)] = 0
        else:                      # backward b->a
            T[(a, b)] = 0; T[(b, a)] = 1
    return T

def tournament_to_tiling(T, tile_list):
    """inverse: read off bit per tile."""
    t = {}
    for (a, b) in tile_list:
        t[(a, b)] = 0 if T[(a, b)] == 1 else 1
    return t

def converse_relabel(T, n):
    """Reverse all arcs then relabel i -> n+1-i.  Returns new adjacency."""
    R = {}
    for (u, v), val in T.items():
        # arc u->v with val ; converse gives v->u ; relabel phi(i)=n+1-i
        pu, pv = n+1-u, n+1-v
        # reversed arc is v->u i.e. relabeled (n+1-v)->(n+1-u)
        if val == 1:               # u->v present
            R[(pv, pu)] = 1
            R[(pu, pv)] = 0
    return R

def count_hamiltonian_paths(T, n):
    """Number of directed Hamiltonian paths via subset DP. Vertices 1..n."""
    verts = list(range(1, n+1))
    idx = {v: i for i, v in enumerate(verts)}
    full = (1 << n) - 1
    # dp[mask][last] = # paths covering 'mask' ending at 'last'
    from collections import defaultdict
    dp = [defaultdict(int) for _ in range(1 << n)]
    for v in verts:
        dp[1 << idx[v]][v] = 1
    for mask in range(1 << n):
        if not dp[mask]:
            continue
        for last, cnt in list(dp[mask].items()):
            for w in verts:
                b = 1 << idx[w]
                if mask & b:
                    continue
                if T.get((last, w), 0) == 1:
                    dp[mask | b][w] += cnt
    return sum(dp[full].values())

# ----------------------------------------------------------------------
# (1)(2)  involution structure
# ----------------------------------------------------------------------

def involution_structure(n):
    tl = tiles(n)
    tset = set(tl)
    m = len(tl)
    fixed, pairs, seen = [], [], set()
    for (a, b) in tl:
        img = rho(a, b, n)
        assert img in tset, f"rho image {img} of {(a,b)} not a tile (n={n})"
        if img == (a, b):
            fixed.append((a, b))
        else:
            key = frozenset([(a, b), img])
            if key not in seen:
                seen.add(key)
                pairs.append(((a, b), img))
    d = len(fixed)
    half = d + len(pairs)
    return tl, m, d, fixed, pairs, half

def verify_no_bitflip(n, exhaustive_upto=5, sample=2000):
    """Confirm rho is a PURE coordinate permutation (converse+relabel never
       flips a tile bit). Returns (#checked, #ok)."""
    import random
    tl = tiles(n)
    m = len(tl)
    rho_map = {p: rho(*p, n) for p in tl}
    checked = ok = 0
    if m <= exhaustive_upto*3 and 2**m <= 200000:
        iterator = itertools.product([0, 1], repeat=m)
    else:
        random.seed(12345)
        iterator = ([random.randint(0, 1) for _ in range(m)] for _ in range(sample))
    for bits in iterator:
        t = {p: bits[i] for i, p in enumerate(tl)}
        T = tiling_to_tournament(t, tl, n)
        R = converse_relabel(T, n)
        t2 = tournament_to_tiling(R, tl)            # tiling of converse-relabel
        # PURE permutation prediction: t2(P) == t(rho(P))
        pred = {p: t[rho_map[p]] for p in tl}
        checked += 1
        if pred == t2:
            ok += 1
    return checked, ok

# ----------------------------------------------------------------------
# (3) grid-sym counting + fraction exponents
# ----------------------------------------------------------------------

def grid_sym_count_exhaustive(n):
    """Count tilings fixed by rho (t==t∘rho) AND cross-check they are exactly
       the phi-self-converse tournaments."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    rho_map = {p: rho(*p, n) for p in tl}
    cnt_fixed = cnt_selfconv = 0
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        is_fixed = all(t[p] == t[rho_map[p]] for p in tl)
        if is_fixed:
            cnt_fixed += 1
        T = tiling_to_tournament(t, tl, n)
        R = converse_relabel(T, n)
        is_selfconv = all(T[k] == R.get(k, 0) for k in T)
        if is_selfconv:
            cnt_selfconv += 1
        assert is_fixed == is_selfconv, "fixed != phi-self-converse mismatch!"
    return cnt_fixed, cnt_selfconv, 2**half, half

# ----------------------------------------------------------------------
# (4) recursions
# ----------------------------------------------------------------------

def half_formula(n):
    return floor((n-1)**2 / 4) if n >= 1 else 0

def verify_recursions(nmax=40):
    h = [half_formula(n) for n in range(0, nmax+1)]
    even_ok = odd_ok = True
    rows = []
    for n in range(4, nmax+1):
        if n % 2 == 0:
            lhs, rhs = h[n], 2*h[n-1] - h[n-2]
            ok = (lhs == rhs)
            even_ok &= ok
            rows.append((n, 'even', lhs, rhs, ok, "2A-C  (A=h(n-1),C=h(n-2))"))
        else:
            lhs, rhs = h[n], 2*h[n-1] - 2*h[n-3] + h[n-4]
            ok = (lhs == rhs)
            odd_ok &= ok
            rows.append((n, 'odd', lhs, rhs, ok, "2A-2E+G (A=h(n-1),E=h(n-3),G=h(n-4))"))
    return rows, even_ok, odd_ok, h

# ----------------------------------------------------------------------
# (5) fundamental-domain shape (ascii)
# ----------------------------------------------------------------------

def fundamental_domain_cells(n):
    """Cells kept in the half-tiling: anti-diagonal a+b=n+1 plus the side
       a+b < n+1 (equivalently 2c+r < n).  Returns set of tiles."""
    tl = tiles(n)
    keep = set()
    for (a, b) in tl:
        s = a + b
        if s <= n+1:               # one side incl. the fixed anti-diagonal
            keep.add((a, b))
    return keep

def ascii_shape(n):
    """Draw the staircage in (r,c): row = r (skip), col = c (lower vertex).
       '#' kept cell, '.' discarded (mirror) cell, '=' fixed diagonal cell."""
    tl = tiles(n)
    keep = fundamental_domain_cells(n)
    maxr = max((rc(a, b)[0] for (a, b) in tl), default=0)
    maxc = max((rc(a, b)[1] for (a, b) in tl), default=0)
    cellset = {rc(a, b): (a, b) for (a, b) in tl}
    lines = []
    for r in range(maxr, 0, -1):
        row = []
        for c in range(1, maxc+1):
            if (r, c) in cellset:
                a, b = cellset[(r, c)]
                if a + b == n+1:
                    row.append('=')           # fixed (diagonal) cell
                elif (a, b) in keep:
                    row.append('#')            # kept
                else:
                    row.append('.')            # discarded mirror cell
            else:
                row.append(' ')
        lines.append("r=%d | %s" % (r, ''.join(row)))
    return lines

# ----------------------------------------------------------------------
# (6) H-spectrum of grid-sym (phi-self-converse) tilings
# ----------------------------------------------------------------------

def grid_sym_H_spectrum(n):
    tl, m, d, fixed, pairs, half = involution_structure(n)
    rho_map = {p: rho(*p, n) for p in tl}
    from collections import Counter
    spec = Counter()
    allspec = Counter()
    hmax_all = 0
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        T = tiling_to_tournament(t, tl, n)
        H = count_hamiltonian_paths(T, n)
        allspec[H] += 1
        hmax_all = max(hmax_all, H)
        if all(t[p] == t[rho_map[p]] for p in tl):
            spec[H] += 1
    return spec, allspec, hmax_all

# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    print("="*72)
    print("HALF-TILING FRAMEWORK  --  kind-pasteur 2026-06-20")
    print("="*72)

    # ---- size table ----
    print("\n[A] Half-tiling size half(n) = #orbits of rho = floor((n-1)^2/4)")
    print(f"{'n':>3} {'m=C(n-1,2)':>11} {'d=fix':>6} {'pairs':>6} {'half':>6} "
          f"{'floor((n-1)^2/4)':>17} {'odd k^2/even k(k-1)':>22}")
    seq = []
    for n in range(2, 13):
        tl, m, d, fixed, pairs, half = involution_structure(n)
        f = half_formula(n)
        if n % 2 == 1:
            k = (n-1)//2; note = f"k={k}, k^2={k*k}"
        else:
            k = n//2; note = f"k={k}, k(k-1)={k*(k-1)}"
        seq.append(half)
        flag = "OK" if half == f else "  <-- MISMATCH"
        print(f"{n:>3} {m:>11} {d:>6} {len(pairs):>6} {half:>6} {f:>17} {note:>22} {flag}")
    print("  sequence half(2..12) =", seq, " (expected 0,1,2,4,6,9,12,16,20,25,30)")

    # ---- (1) no bit flip ----
    print("\n[B] rho = converse+relabel is a PURE coordinate permutation (no bit flip)")
    for n in range(3, 8):
        chk, ok = verify_no_bitflip(n)
        mode = "exhaustive" if 2**len(tiles(n)) <= 200000 else "sampled"
        print(f"  n={n}: {ok}/{chk} tilings confirm  t_converse(P)=t(rho(P))  [{mode}]  "
              f"{'OK' if ok==chk else 'FAIL'}")

    # ---- (3) grid-sym count & fraction ----
    print("\n[C] grid-symmetric (phi-self-converse) tilings = 2^half(n)")
    print(f"{'n':>3} {'#fixed':>8} {'#selfconv':>10} {'2^half':>8} {'frac':>12} "
          f"{'log2 frac':>10} {'canon exp':>10}")
    canon_exp = {3:0, 4:-1, 5:-2, 6:-4, 7:-6, 8:-9}
    for n in range(3, 8):
        cf, cs, p2, half = grid_sym_count_exhaustive(n)
        m = len(tiles(n))
        frac = cf / 2**m
        log2 = (cf and (half - m))
        ce = canon_exp.get(n, None)
        print(f"{n:>3} {cf:>8} {cs:>10} {p2:>8} {frac:>12.6f} {log2:>10} "
              f"{str(ce):>10} {'OK' if (cf==cs==p2 and log2==ce) else 'CHECK'}")

    # ---- (4) recursions ----
    print("\n[D] Recursions for half(n)=floor((n-1)^2/4)")
    rows, even_ok, odd_ok, h = verify_recursions(40)
    for (n, par, lhs, rhs, ok, desc) in rows[:9]:
        print(f"  n={n:>2} [{par:>4}] half={lhs:>3} = {rhs:>3}  {desc}  {'OK' if ok else 'FAIL'}")
    print(f"  ... verified to n=40:  EVEN recursion (2A-C) all-OK = {even_ok};  "
          f"ODD recursion (2A-2E+G) all-OK = {odd_ok}")
    # also verify the user's literal signed 7-term sum for odd n (C,D cancel)
    print("  user's literal ODD A+B-C+D-E-F+G with sizes (n-1,n-1,n-2,n-2,n-3,n-3,n-4):")
    lit_ok = True
    for n in range(5, 30, 2):
        A=B=h[n-1]; C=D=h[n-2]; E=F=h[n-3]; G=h[n-4]
        val = A+B-C+D-E-F+G
        ok = (val == h[n]); lit_ok &= ok
    print(f"    all-OK to n=29 = {lit_ok}   (note: -C+D cancels -> 2A-2E+G)")
    print("  user's literal EVEN A+B-C with A,B=h(n-1), C=h(n-2):")
    lit_ok2 = all((h[n-1]+h[n-1]-h[n-2])==h[n] for n in range(4,30,2))
    print(f"    all-OK to n=28 = {lit_ok2}")

    # ---- (5) shapes ----
    print("\n[E] Fundamental-domain shape  ('=' fixed/diagonal, '#' kept, '.' discarded)")
    for n in [5, 6, 7, 8, 9]:
        keep = fundamental_domain_cells(n)
        print(f"  --- n={n}  (kept cells = half(n) = {half_formula(n)}) ---")
        for line in ascii_shape(n):
            print("    " + line)

    # ---- gnomon decomposition for odd n ----
    print("\n[F] Odd-n gnomon decomposition  k^2 = 1+3+5+...+(2k-1)")
    for n in range(3, 14, 2):
        k = (n-1)//2
        gn = [2*i-1 for i in range(1, k+1)]
        print(f"  n={n}: k={k}, half={k*k} = sum{gn} = {sum(gn)}  "
              f"({'OK' if sum(gn)==half_formula(n) else 'FAIL'})")

    # ---- (6) H-spectrum ----
    print("\n[G] H-spectrum of grid-symmetric (phi-self-converse) tilings")
    for n in range(3, 8):
        spec, allspec, hmax = grid_sym_H_spectrum(n)
        gmax = max(spec) if spec else 0
        print(f"  n={n}: H_max(all tilings)={hmax}  H_max(grid-sym)={gmax}  "
              f"{'== matches global max' if gmax==hmax else '!= (grid-sym misses global max)'}")
        print(f"        grid-sym H-spectrum = {dict(sorted(spec.items()))}")
    print("\nDONE.")

if __name__ == "__main__":
    main()
