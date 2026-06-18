"""
lrc14_refute_residue-orbit_M4_ceiling_kps-S2-wf.py

FAST + EXHAUSTIVE absolute ceiling for the M4 band-depth orbit-majority map (mod 14),
n=5 runners.

M4: vertices = runners. For each unit a in (Z/14)* = {1,3,5,9,11,13}, residue
    r_i(a)=v_i*a mod 14, depth d_i(a)=min(r_i(a),14-r_i(a)) in {0..7}.
    Arc i->j iff #{a: d_i>d_j} > #{a: d_j>d_i}. Tie (equal counts) -> tie-break by raw
    speed (smaller speed -> arc OUT), matching method4_band_majority exactly.

KEY FACTS used to make this exhaustive and exact:
 * M4 depends ONLY on the residue VALUES (v_i mod 14), via the pairwise depth-majority
   outcome WIN[x][y] in {+1,-1,0}.  We precompute WIN for all 14x14 residue value pairs.
 * Loneliness on the grid only needs section 0 empty; it does NOT force distinct residues
   (SDR is strictly stronger). So two runners may share a residue. We therefore allow
   REPETITION: enumerate all residue MULTISETS of size 5 from {0..13} = C(18,5) = 8568.
   (Residue 0 included; including the parked section only ENLARGES the ceiling.)
 * The ONLY freedom beyond WIN is how tied pairs (WIN=0) get oriented. Tie-break is by raw
   speed; speeds map to residues arbitrarily, so EVERY strict total order of the 5 runners
   is achievable. Hence for each multiset we sweep all 5! = 120 runner orderings.

This is the MAXIMUM set M4 can ever produce over ANY runner-vertex assignment. If a
forbidden class is not hit here, it is unreachable by M4.

Exact integer arithmetic throughout.
"""
from math import gcd
from itertools import combinations, permutations, combinations_with_replacement

MOD = 14
UNITS = [a for a in range(1, MOD) if gcd(a, MOD) == 1]  # {1,3,5,9,11,13}
def depth(r):
    r %= MOD
    return min(r, MOD - r)

WIN = [[0]*MOD for _ in range(MOD)]
for x in range(MOD):
    for y in range(MOD):
        wi = wj = 0
        for a in UNITS:
            di = depth(x*a); dj = depth(y*a)
            if di > dj: wi += 1
            elif dj > di: wj += 1
        WIN[x][y] = 1 if wi > wj else (-1 if wj > wi else 0)

m = 5

def build_adj(res, order):
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            w = WIN[res[i]][res[j]]
            if w == 1: adj[i][j] = 1
            elif w == -1: adj[j][i] = 1
            else:
                if order[i] < order[j]: adj[i][j] = 1
                else: adj[j][i] = 1
    return adj

def canon(adj):
    best = None
    for p in permutations(range(m)):
        b = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or b < best: best = b
    return best
def H(adj):
    c = 0
    for p in permutations(range(m)):
        if all(adj[p[k]][p[k+1]] for k in range(m-1)): c += 1
    return c
def c3(adj):
    c = 0
    for a, b, cc in combinations(range(m), 3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c += 1
    return c
def score(adj):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))
def signature(adj):
    return (H(adj), c3(adj), score(adj))

FORB = {
    (9, 3, (1,1,2,3,3)),
    (13,4, (1,2,2,2,3)),
    (15,4, (1,2,2,2,3)),
    (15,5, (2,2,2,2,2)),
}

def main():
    realized_keys = set()
    realized_sig = {}
    forb_hits = []
    n_multisets = 0
    for res in combinations_with_replacement(range(MOD), m):
        n_multisets += 1
        tie = False
        for i in range(m):
            for j in range(i+1, m):
                if WIN[res[i]][res[j]] == 0:
                    tie = True; break
            if tie: break
        orders = permutations(range(m)) if tie else [tuple(range(m))]
        for order in orders:
            adj = build_adj(res, order)
            realized_keys.add(canon(adj))
            sg = signature(adj)
            if sg not in realized_sig:
                realized_sig[sg] = (res, order)
            if sg in FORB:
                forb_hits.append((res, order, sg))
    print(f"enumerated {n_multisets} residue multisets (C(18,5)=8568 expected)")
    print(f"M4 absolute ceiling: {len(realized_keys)} distinct iso KEYS")
    print(f"M4 absolute ceiling: {len(realized_sig)} distinct SIGNATURES (H,#3cyc,score)")
    print()
    print("ALL realized signatures (the M4 ceiling):")
    for sg in sorted(realized_sig):
        flag = "  <-- FORBIDDEN!" if sg in FORB else ""
        print(f"   {sg}{flag}")
    print()
    hit = sorted(set(s for *_, s in forb_hits))
    if forb_hits:
        print(f"*** FORBIDDEN CLASSES REALIZED: {hit} ***  -> CLAIM REFUTED")
        for w in forb_hits[:8]:
            print("    witness res(multiset), tie-order, sig:", w)
    else:
        print("NO forbidden class realized at the absolute residue ceiling.")
        print(f"Forbidden classes confirmed UNREACHABLE by M4: {sorted(FORB)}")
    print()
    print(f"Total signatures realized: {len(realized_sig)}  (claim: 8 reachable of 12 free)")

if __name__ == "__main__":
    main()
