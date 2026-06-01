#!/usr/bin/env python3
"""
lrc_realizable_isoclasses_s523.py    oracle-2026-06-01-S523

"LRC at a particular n IS a tournament question, and the structure of the set of
possible isomorphism classes."  (User prompt.)

THM-381: LRC(n) <=> for every speed system the clock-movie t -> T_S(t) (tournament
on n vertices, observer 0 marked) makes vertex 0 a SOURCE.  Runner-runner arcs use
the half-turn comparator on circle positions.

A half-turn tournament of points on a circle is a ROUND tournament: each vertex's
out-set is a contiguous clockwise ARC.  So the movie can realize ONLY round
tournaments -- a strict, small subclass of A000568(n).  We:

 (a) generate ROUND tournament iso-classes directly for m=3..8 (clockwise-arc
     out-degree vectors, dedup, canon) -> the realizable sequence,
 (b) for m=3..6 do the FULL enumeration and verify
        round  ==  locally-transitive  ==  half-turn-realizable,
 (c) tabulate against A000568(m) and the S520 reachable-LONELY menu.
"""
from itertools import combinations, permutations, product
import random

# ---------- canon ----------
_PERMS = {}
def canon(adj):
    m = len(adj)
    perms = _PERMS.setdefault(m, list(permutations(range(m))))
    best = None
    for p in perms:
        flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or flat < best: best = flat
    return best

# ---------- direct ROUND generation ----------
def round_valid(d, m):
    """d-vector valid round tournament? check each pair directly (no adj)."""
    for i in range(m):
        di = d[i]
        for j in range(i+1, m):
            k = j - i
            if (k <= di) == ((m - k) <= d[j]):
                return False
    return True

def round_classes(m):
    """All round tournaments: out-set of position i = next d_i vertices clockwise.
    Validity checked from d directly; build adj + canon only for valid, distinct."""
    seen = set(); advs = set()
    for d in product(range(m), repeat=m):
        if not round_valid(d, m): continue
        adj = [[0]*m for _ in range(m)]
        for i in range(m):
            for k in range(1, d[i]+1):
                adj[i][(i+k) % m] = 1
        key = tuple(tuple(r) for r in adj)
        if key in advs: continue
        advs.add(key); seen.add(canon(adj))
    return seen

# ---------- full enumeration helpers (small m) ----------
def all_classes(m):
    pairs = list(combinations(range(m), 2)); seen = {}
    for bits in product((0,1), repeat=len(pairs)):
        adj = [[0]*m for _ in range(m)]
        for (i,j),b in zip(pairs,bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        c = canon(adj)
        if c not in seen: seen[c]=adj
    return seen

def locally_transitive(adj):
    m=len(adj)
    def trans(S):
        S=list(S); return sorted(sum(1 for b in S if a!=b and adj[a][b]) for a in S)==list(range(len(S)))
    for v in range(m):
        if not trans([w for w in range(m) if w!=v and adj[v][w]]): return False
        if not trans([w for w in range(m) if w!=v and adj[w][v]]): return False
    return True

def half_turn_realizable(adj, samples=4000):
    """is adj the half-turn tournament of SOME circle config? sample gap vectors."""
    m=len(adj); target=canon(adj)
    for _ in range(samples):
        pos=sorted(random.random() for _ in range(m))
        a=[[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                if i==j: continue
                dd=(pos[j]-pos[i])%1.0
                if 0<dd<0.5: a[i][j]=1
        if canon(a)==target: return True
    return False

A000568={1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}
MENU={4:1,5:2,6:6,7:6,8:'>=12'}  # reachable LONELY menu (S520), indexed by n; m=n-1

def main():
    print("LRC realizable iso-class structure (oracle-S523)\n")
    # (a) round sequence m=3..7
    rnd_count={}; rnd_sets={}
    for m in range(3,8):
        rc=round_classes(m); rnd_count[m]=len(rc); rnd_sets[m]=rc
        print(f" m={m}: ROUND iso-classes = {len(rc):4d}   (A000568({m})={A000568[m]})")
    print()
    print(" ROUND sequence m=3..7 :", [rnd_count[m] for m in range(3,8)])
    print(" A000568        m=3..7 :", [A000568[m] for m in range(3,8)])
    print()
    # (b) equivalence check on small m (cheap full enumeration)
    print("EQUIVALENCE CHECK (full enumeration):  round == locally-transitive == half-turn")
    for m in range(3,6):
        cls=all_classes(m)
        lt={c for c,a in cls.items() if locally_transitive(a)}
        rd=rnd_sets[m]
        ht={c for c,a in cls.items() if half_turn_realizable(a, samples=1500)}
        print(f"  m={m}: |all|={len(cls)} |round|={len(rd)} |loc.trans|={len(lt)} "
              f"round==loctrans:{lt==rd}  round==half-turn:{ht==rd}")
    print()
    # (c) the nested chain vs menu
    print("THE NESTED CHAIN for LRC(n)  [m=n-1 runners]:")
    print(f"  {'n':>2} {'m=n-1':>5} {'A000568(m)':>10} {'round(m)':>8} {'lonely-menu':>11}")
    for n in range(4,8):
        m=n-1
        print(f"  {n:>2} {m:>5} {A000568[m]:>10} {rnd_count[m]:>8} {str(MENU[n]):>11}")
    print("\n ALL tournaments ⊋ ROUND (realizable movie) ⊋ source-marked ⊋ LONELY menu.")
    print(" Extremal symmetric round tournament = regular m-gon R_m (S522, H=175 at m=7).")

if __name__=="__main__":
    main()
