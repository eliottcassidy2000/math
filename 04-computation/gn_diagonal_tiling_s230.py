#!/usr/bin/env python3
"""
gn_diagonal_tiling_s230.py — The y=x diagonal as the keystone
opus-2026-03-23-S230

THE INSIGHT: The staircase δ_{n-2} has a y=x diagonal. This diagonal
is the LINE OF COMPLEMENT SYMMETRY. A tournament is SC iff its tiling
is symmetric under y=x reflection (up to vertex relabeling).

For SC tournaments with anti-automorphism σ:
- Arcs at cells FIXED by σ → flipping preserves SC → BLUE edges
- Arcs at cells NOT fixed by σ → flipping breaks SC → BLACK edges

The simplest σ (vertex reversal i→n-1-i) fixes the ANTI-DIAGONAL
cells where r+c = n-3 (in 0-indexed staircase coordinates).

This explains:
1. H = 1 + 2^(n-2) for the first blue neighbor (anti-diagonal flip)
2. trans_SC_deg = floor((n-1)/2) (# anti-auto fixed arcs)
3. Blue edges ≈ anti-diagonal flips, Black ≈ off-diagonal flips

Let's VERIFY this tiling-level structure by computing, for each edge
in the merged metagraph, WHICH CELL of the staircase was flipped,
and checking if anti-diagonal cells produce blue edges.
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE y=x DIAGONAL: TILING-LEVEL BLUE/BLACK STRUCTURE")
print("  opus-2026-03-23-S230")
print("=" * 80)

# ============================================================================
# STAIRCASE COORDINATES
# ============================================================================

def staircase_position(i, j, n):
    """Map arc (i,j) with i<j to staircase coordinate (row, col).
    The staircase δ_{n-2} has rows r=0..n-3, with row r having n-2-r cells.
    Arc (i,j) maps to row r = j-i-1, col c = i.
    The anti-diagonal is where r+c = n-3, i.e., j-i-1+i = n-3, i.e., j = n-2."""
    r = j - i - 1  # row (0-indexed)
    c = i           # column (0-indexed)
    return r, c

def is_antidiagonal(r, c, n):
    """Check if cell (r,c) is on the anti-diagonal r+c = n-3."""
    return r + c == n - 3

def arc_range(i, j):
    """Range = j - i - 1 (distance between vertices minus 1)."""
    return j - i - 1

# ============================================================================
# BUILD + CLASSIFY EACH ARC FLIP BY STAIRCASE POSITION
# ============================================================================

for n in range(3, 8):
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]

    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0 = time.time()

    # Build classes (small n brute force)
    if n <= 6:
        AP = list(permutations(range(n)))
        def b2a(b):
            a=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(P):
                if b&(1<<k): a[i][j]=1
                else: a[j][i]=1
            return a
        def cn(a):
            best=None
            for p in AP:
                f=tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
                if best is None or f<best: best=f
            return best
        cm={}
        for b in range(1<<m):
            a=b2a(b); c=cn(a)
            if c not in cm: cm[c]={'bits':b,'adj':a,'size':0}
            cm[c]['size']+=1
        cls=list(cm.values())
        for i,c in enumerate(cls): c['cid']=i
        ny=False
    else:
        r=subprocess.run(['gentourng',str(n)],capture_output=True,text=True)
        ls=[l.strip() for l in (r.stdout or '').split('\n') if len(l.strip())==m and all(c in '01' for c in l.strip())]
        def b2a_n(b):
            a=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(P):
                if b&(1<<(m-1-k)): a[i][j]=1
                else: a[j][i]=1
            return a
        cls=[{'bits':int(l,2),'adj':b2a_n(int(l,2)),'cid':i} for i,l in enumerate(ls)]
        ny=True

    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c=0
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    if a[i][j] and a[j][k] and a[k][i]: c+=1
                    if a[i][k] and a[k][j] and a[j][i]: c+=1
        return c
    def Hf(a):
        dp={}
        for v in range(n): dp[(1<<v,v)]=1
        for S in range(1,1<<n):
            for v in range(n):
                if not(S&(1<<v)): continue
                val=dp.get((S,v),0)
                if val==0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: key=(S|(1<<u),u); dp[key]=dp.get(key,0)+val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(a):
        scores=[sum(a[i][j] for j in range(n)) for i in range(n)]
        sg=defaultdict(list)
        for v in range(n): sg[scores[v]].append(v)
        gs=[sg[s] for s in sorted(set(scores))]
        if all(len(g)==1 for g in gs):
            p=[g[0] for g in gs]
            return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        best=None
        def gen(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gen(gs2[1:]): yield list(p)+rest
        for p in gen(gs):
            f=tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f<best: best=f
        return best

    hm=defaultdict(list); cm2={}
    for cl in cls:
        a=cl['adj']; cl['score']=ss(a); cl['c3']=c3(a); cl['H']=Hf(a)
        cl['own_canon']=cp(a)
        comp=[[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon']=cp(comp); cl['sc']=(cl['own_canon']==cl['comp_canon'])
        hm[(cl['score'],cl['c3'],cl['H'])].append(cl['cid'])
        cm2[cl['own_canon']]=cl['cid']
    for cl in cls: cl['comp_cid']=cm2.get(cl['comp_canon'],-1)

    def lk(a):
        s=ss(a);c=c3(a);h=Hf(a)
        cids=hm.get((s,c,h))
        if not cids: return -1
        if len(cids)==1: return cids[0]
        return cm2.get(cp(a),-1)

    # Merged graph
    mn={}; mid=0
    for cl in cls:
        ci=cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci]=mid; mid+=1
        else:
            mn[ci]=mid
            comp=cl['comp_cid']
            if comp>=0 and comp!=ci: mn[comp]=mid
            mid+=1

    mi={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        mi[mm]={'H':cls[cids[0]]['H'],'sc':cls[cids[0]]['sc']}

    # For each arc flip, classify by staircase position AND blue/black
    # Count: how many blue vs black edges come from each staircase cell
    cell_blue = Counter()  # (r,c) → # blue lines from this cell
    cell_black = Counter() # (r,c) → # black lines from this cell
    cell_self = Counter()  # (r,c) → # self-loop lines from this cell
    cell_total = Counter() # (r,c) → total lines

    for cl in cls:
        b = cl['bits']
        for arc_idx in range(m):
            i_arc, j_arc = P[arc_idx]
            r, c = staircase_position(i_arc, j_arc, n)

            if ny:
                fb = b ^ (1 << (m-1-arc_idx))
            else:
                fb = b ^ (1 << arc_idx)
            if ny:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(P):
                    if fb&(1<<(m-1-k)): fa[i][j]=1
                    else: fa[j][i]=1
            else:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(P):
                    if fb&(1<<k): fa[i][j]=1
                    else: fa[j][i]=1

            nb = lk(fa)
            if nb == cl['cid']:
                cell_self[(r,c)] += 1
            elif nb >= 0:
                ma_nb = mn[nb]; ma_cl = mn[cl['cid']]
                if ma_nb == ma_cl:
                    cell_self[(r,c)] += 1  # collapsed in merged
                elif mi[ma_nb]['sc'] == mi[ma_cl]['sc']:
                    cell_blue[(r,c)] += 1
                else:
                    cell_black[(r,c)] += 1
            cell_total[(r,c)] += 1

    # Print staircase grid
    print(f"\n  STAIRCASE δ_{n-2} GRID ({n-2} rows):")
    print(f"  Cell format: B/K (blue/black line count)")
    print(f"  Anti-diagonal cells marked with [*]")
    print()

    max_row = n - 2
    for r in range(max_row):
        row_str = f"    r={r}: "
        for c in range(max_row - r):
            b = cell_blue.get((r,c), 0)
            k = cell_black.get((r,c), 0)
            s = cell_self.get((r,c), 0)
            ad = "[*]" if is_antidiagonal(r, c, n) else "   "
            row_str += f"{b:3d}/{k:3d}{ad} "
        print(row_str)

    # Summary: anti-diagonal vs off-diagonal
    antidiag_blue = sum(cell_blue.get((r,c),0) for r in range(max_row) for c in range(max_row-r) if is_antidiagonal(r,c,n))
    antidiag_black = sum(cell_black.get((r,c),0) for r in range(max_row) for c in range(max_row-r) if is_antidiagonal(r,c,n))
    offdiag_blue = sum(cell_blue.get((r,c),0) for r in range(max_row) for c in range(max_row-r) if not is_antidiagonal(r,c,n))
    offdiag_black = sum(cell_black.get((r,c),0) for r in range(max_row) for c in range(max_row-r) if not is_antidiagonal(r,c,n))

    total_blue = antidiag_blue + offdiag_blue
    total_black = antidiag_black + offdiag_black

    print(f"\n  ANTI-DIAGONAL vs OFF-DIAGONAL:")
    print(f"    Anti-diagonal: blue={antidiag_blue}, black={antidiag_black}, "
          f"blue_frac={antidiag_blue/(antidiag_blue+antidiag_black):.4f}" if antidiag_blue+antidiag_black>0 else "")
    print(f"    Off-diagonal:  blue={offdiag_blue}, black={offdiag_black}, "
          f"blue_frac={offdiag_blue/(offdiag_blue+offdiag_black):.4f}" if offdiag_blue+offdiag_black>0 else "")
    print(f"    Total:         blue={total_blue}, black={total_black}, "
          f"blue_frac={total_blue/(total_blue+total_black):.4f}" if total_blue+total_black>0 else "")

    # Blue fraction per row (distance from main diagonal)
    print(f"\n  BLUE FRACTION BY ROW (r = range - 1):")
    for r in range(max_row):
        row_blue = sum(cell_blue.get((r,c),0) for c in range(max_row-r))
        row_black = sum(cell_black.get((r,c),0) for c in range(max_row-r))
        row_total = row_blue + row_black
        frac = row_blue/row_total if row_total > 0 else 0
        n_cells = max_row - r
        print(f"    r={r} ({n_cells} cells): blue={row_blue:5d}, black={row_black:5d}, frac={frac:.4f}")

    # Blue fraction per column
    print(f"\n  BLUE FRACTION BY COLUMN (c = source vertex):")
    for c in range(max_row):
        col_blue = sum(cell_blue.get((r,c),0) for r in range(max_row-c))
        col_black = sum(cell_black.get((r,c),0) for r in range(max_row-c))
        col_total = col_blue + col_black
        frac = col_blue/col_total if col_total > 0 else 0
        n_cells = max_row - c
        print(f"    c={c} ({n_cells} cells): blue={col_blue:5d}, black={col_black:5d}, frac={frac:.4f}")

    # Check: is blue fraction UNIFORM across cells?
    fracs = []
    for r in range(max_row):
        for c in range(max_row-r):
            b = cell_blue.get((r,c),0)
            k = cell_black.get((r,c),0)
            if b+k > 0:
                fracs.append(b/(b+k))
    if fracs:
        print(f"\n  BLUE FRACTION UNIFORMITY:")
        print(f"    min={min(fracs):.4f}, max={max(fracs):.4f}, "
              f"std={((sum((f-sum(fracs)/len(fracs))**2 for f in fracs)/len(fracs))**0.5):.6f}")
        is_uniform = max(fracs) - min(fracs) < 0.01
        print(f"    UNIFORM: {'YES ✓' if is_uniform else 'NO ✗'}")

    elapsed = time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
