"""
Apex tile <-> H=49 (=7^2), 75 in the tiling model (kind-pasteur S31f).

Owner lead: "balanced cuts supply the strong n=8 H-values; a single unbalanced weight
w=1 supplies 49 and 75 -- probably the APEX TILE in the tiling model."  49 = 7^2 ties
to the apex prime (E_7 / {7,21} / LRC-14).  Test it directly: in the fixed-base-path
tiling model (base path n-1 -> ... -> 0, tiles = non-consecutive arcs (x,y) x>y+1),
the APEX TILE = (n-1, 0) (the longest-range arc, the geometric apex of the staircase).

For each tile configuration we compute H (Hamiltonian-path count) and record, per
achievable H value, whether the apex tile is flipped and the tile-weight (#flipped
tiles).  Question: are H=49 and H=75 reachable ONLY with the apex tile flipped?
"""
import sys
from itertools import combinations
from collections import defaultdict

def hcount(n, out):
    # out[v] = bitmask of vertices v beats. H = # Hamiltonian paths (directed).
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    full = (1<<n)-1
    for mask in range(1<<n):
        for v in range(n):
            c = dp[mask][v]
            if not c: continue
            ov = out[v] & ~mask
            w = ov
            while w:
                u = (w & -w).bit_length()-1
                dp[mask|(1<<u)][u] += c
                w &= w-1
    return sum(dp[full][v] for v in range(n))

def analyze(n, targets=(49,75)):
    # base path: i -> i-1 for i=1..n-1  (vertex i beats i-1)
    base = [0]*n
    for i in range(1, n):
        base[i] |= (1 << (i-1))     # i beats i-1
    # tiles: arcs (x,y) with x > y+1.  default orientation x->y (x beats y).
    tiles = [(x,y) for y in range(n) for x in range(y+2, n)]
    apex = (n-1, 0)
    apex_idx = tiles.index(apex)
    T = len(tiles)                  # = C(n-1,2)
    print(f"n={n}: {T} tiles, apex tile = {apex} (index {apex_idx}), enumerating 2^{T}={1<<T} configs")
    # for each config (bitmask over tiles): bit set => tile flipped to y->x (y beats x)
    Hinfo = defaultdict(lambda: {"apex_yes":0, "apex_no":0, "weights":set(), "min_w":99, "count":0})
    for cfg in range(1<<T):
        out = base[:]
        for t in range(T):
            x,y = tiles[t]
            if (cfg>>t)&1:
                out[y] |= (1<<x)     # flipped: y beats x
            else:
                out[x] |= (1<<y)     # default: x beats y
        H = hcount(n, out)
        info = Hinfo[H]
        info["count"] += 1
        w = bin(cfg).count("1")
        info["weights"].add(w)
        if w < info["min_w"]: info["min_w"] = w
        if (cfg>>apex_idx)&1: info["apex_yes"] += 1
        else: info["apex_no"] += 1
    # report targets + summary
    print(f"  achievable H values: {len(Hinfo)}; sample max H = {max(Hinfo)}")
    for tH in sorted(targets):
        if tH in Hinfo:
            i = Hinfo[tH]
            apex_required = (i["apex_no"] == 0)
            print(f"  H={tH}: count={i['count']}  apex_flipped={i['apex_yes']}  apex_NOT={i['apex_no']}  "
                  f"=> APEX REQUIRED? {apex_required};  min tile-weight={i['min_w']}")
        else:
            print(f"  H={tH}: NOT achievable at n={n} (fixed base path)")
    # which H values REQUIRE the apex tile (apex_no==0)?
    apex_only = sorted(h for h,i in Hinfo.items() if i["apex_no"]==0 and i["count"]>0)
    print(f"  H values REQUIRING the apex tile flipped (apex_no==0): {apex_only[:40]}{' ...' if len(apex_only)>40 else ''}")

if __name__ == "__main__":
    analyze(7)
    if len(sys.argv) > 1 and sys.argv[1] == "n8":
        analyze(8)
