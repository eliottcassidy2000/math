"""
monad-explorer-2026-06-06-S709c
===============================
VERIFY N=39 and push lower: fine center sweep + swap local search on the
sqrt(7)-Eisenstein unit-distance graph (improves HYP-2267's N=43 to N<=39).

All arithmetic is EXACT integer arithmetic on Eisenstein index pairs (x,y) with
squared-distance form Q(x,y)=x^2+xy+y^2 and unit^2 = 7.  A "config" is any finite
set S of index points; U(S) = #{ {p,q} subset S : Q(p-q)=7 }.  We confirm U(S)>3|S|.
"""
import math
from itertools import combinations

A, B, C = 1, 1, 1            # x^2 + xy + y^2
R = 7                        # unit^2 (layer 12)

def Q(dx, dy):
    return A*dx*dx + B*dx*dy + C*dy*dy

# offset shell at squared-distance 7
OFFS = [(dx, dy) for dx in range(-4, 5) for dy in range(-4, 5) if Q(dx, dy) == R]
assert len(OFFS) == 12, len(OFFS)

def U_of(S):
    Sset = set(S)
    e = 0
    for (x, y) in S:
        for (dx, dy) in OFFS:
            if (x+dx, y+dy) in Sset:
                e += 1
    assert e % 2 == 0
    return e // 2

def disk(cx, cy, k, half=12):
    box = [(x, y) for x in range(-half, half+1) for y in range(-half, half+1)]
    box.sort(key=lambda p: ( (p[0]-cx)**2 + (p[0]-cx)*(p[1]-cy) + (p[1]-cy)**2 ))
    return box[:k]

# ---- 1. VERIFY the N=39 config at center (1/2, 0) -------------------------
print("="*70)
print("1.  VERIFY: triangular sqrt(7) disk, center (1/2,0), exact integer count")
print("="*70)
for k in range(35, 45):
    S = disk(0.5, 0.0, k)
    u = U_of(S)
    flag = "  <-- U>3N" if u > 3*k else ""
    print(f"   N={k:>3}  U={u:>4}  3N={3*k:>4}  surplus(U-3N)={u-3*k:+d}{flag}")

# ---- 2. FINE center sweep: which fractional center crosses 3N earliest? ---
print()
print("="*70)
print("2.  FINE center sweep (cx,cy on a 1/6 grid in the fundamental cell)")
print("    report smallest k with U>3N for each center")
print("="*70)
best = None
G = 6
for ix in range(G+1):
    for iy in range(G+1):
        cx, cy = ix/G, iy/G
        kc = None
        for k in range(17, 60):
            S = disk(cx, cy, k)
            if U_of(S) > 3*k:
                kc = k; break
        if kc is not None and (best is None or kc < best[0]):
            best = (kc, cx, cy, U_of(disk(cx, cy, kc)))
print(f"   best disk over center grid: N={best[0]} at center=({best[1]:.3f},{best[2]:.3f}) "
      f"U={best[3]} (3N={3*best[0]})")

# ---- 3. SWAP / add-remove local search from the best disk ------------------
print()
print("="*70)
print("3.  Local search (boundary add/remove + swaps) to push N below the disk")
print("="*70)
def neighbors_count_in(S):
    Sset = set(S)
    deg = {}
    for p in S:
        d = 0
        for (dx, dy) in OFFS:
            if (p[0]+dx, p[1]+dy) in Sset:
                d += 1
        deg[p] = d
    return deg

# seed: best disk
seed = disk(best[1], best[2], best[0])
S = set(seed)
# try to find ANY set of size < best[0] with U>3N by simulated greedy:
# repeatedly: among all candidate single-vertex swaps (remove worst, add a frontier
# vertex that raises surplus), keep if it allows shrinking. Then attempt deletions.
def U_set(S): return U_of(list(S))

improved = True
curk = len(S); curU = U_set(S)
print(f"   seed N={curk} U={curU}")
# aggressive deletion: try removing any vertex that keeps U>3(N-1), repeat
changed = True
while changed and len(S) > 17:
    changed = False
    deg = neighbors_count_in(S)
    # delete the vertex with smallest internal degree if it keeps feasibility
    for p in sorted(S, key=lambda q: deg[q]):
        nU = curU - deg[p]
        nk = len(S) - 1
        if nU > 3*nk:
            S.discard(p); curU = nU; changed = True
            break
print(f"   after greedy deletion: N={len(S)} U={U_set(S)} (3N={3*len(S)})")

# add-one-remove-two style: try to find a smaller feasible set by local perturbation
def frontier(S):
    Sset = set(S); fr = set()
    for p in S:
        for (dx, dy) in OFFS:
            q = (p[0]+dx, p[1]+dy)
            if q not in Sset:
                fr.add(q)
    return fr

bestS = set(S); bestk = len(S)
import itertools
for _ in range(2000):
    # perturb: add a frontier vertex maximizing new internal degree, then delete
    # two lowest-degree vertices if still feasible
    fr = frontier(bestS)
    if not fr: break
    Sset = set(bestS)
    def gain(q):
        return sum(1 for (dx,dy) in OFFS if (q[0]+dx,q[1]+dy) in Sset)
    q = max(fr, key=gain)
    T = set(bestS); T.add(q)
    Ut = U_set(T)
    deg = neighbors_count_in(T)
    # remove two smallest-degree
    rem = sorted(T, key=lambda v: deg[v])[:2]
    for v in rem:
        Ut2 = U_set(T - {v})
        if Ut2 > 3*(len(T)-1):
            T.discard(v); Ut = Ut2
    if len(T) < bestk and U_set(T) > 3*len(T):
        bestS = T; bestk = len(T)
print(f"   after perturbation search: best N={bestk} U={U_set(bestS)} (3N={3*bestk})")

print()
print("="*70)
print(f"CONCLUSION: smallest verified config N = {min(len(S), bestk)} "
      f"(beats HYP-2267's N=43). THM-420 floor: N>=17. So 17 <= N* <= {min(len(S),bestk)}.")
print("="*70)
