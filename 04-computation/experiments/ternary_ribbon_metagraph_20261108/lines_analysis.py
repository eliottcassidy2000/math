#!/usr/bin/env python3
"""lines_analysis.py -- node-level statistics of the metagraph lines from lines_nodes_n{n}.json (n = 4..8):
checks of THM-643 (fibre sizes H/|Aut|, blue endpoints odd at self-converse nodes, black even there, black odd at
non-self-converse nodes), and the line counts between the new node classifications:
  S: self-converse (SC) / not (NS);  K: number of strong components;  V: eternal colour of the value H (THM-4930 lists,
  exact below 1131);  F: first creator (least order realising H) / not;  plus multiplicities between node pairs."""
import json, sys
from collections import Counter, defaultdict
BLACK = set([7, 21, 19, 23, 35, 39, 43, 47, 53, 55, 59, 63, 101, 103, 107, 109, 113, 119, 131, 163, 167, 181, 183, 185, 203, 205, 207, 209, 215,
             231, 235, 239, 251, 267, 271, 291, 299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359, 377, 415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523, 527, 563,
             671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983, 1003, 1009, 1067, 1087, 1115])
def colour(h):
    if h >= 1131: return 'unknown'
    return 'black' if h in BLACK else 'blue'
data = {n: json.load(open(f"lines_nodes_n{n}.json")) for n in range(4, 9)}
born = {}
for n in range(4, 9):
    for cid, c in data[n]["classes"].items():
        born.setdefault(c["H"], n)
for h in (1, 3, 5): born.setdefault(h, {1: 1, 3: 3, 5: 4}[h])
for n in range(5, 9):
    D = data[n]; cls = {int(k): v for k, v in D["classes"].items()}
    lines = D["lines"]
    ends = Counter(); blue_ends = Counter(); black_ends = Counter()
    pair_mult = Counter(); pair_blue = Counter()
    for c1, c2, b, cnt in lines:
        for c in (c1, c2):
            ends[c] += cnt
            (blue_ends if b else black_ends)[c] += cnt
        key = (min(c1, c2), max(c1, c2))
        pair_mult[key] += cnt
        if b: pair_blue[key] += cnt
    # THM-643 checks
    ok_fibre = all(ends[c] == cls[c]["H"] // cls[c]["aut"] * 1 for c in cls) if all(cls[c]["H"] % cls[c]["aut"] == 0 for c in cls) else False
    fib = {c: ends[c] for c in cls}
    sc = {c for c in cls if cls[c]["sc"]}
    parity_ok = all(blue_ends[c] % 2 == 1 and black_ends[c] % 2 == 0 for c in sc) and all(blue_ends[c] == 0 and black_ends[c] % 2 == 1 for c in cls if c not in sc)
    nblue = sum(cnt for c1, c2, b, cnt in lines if b); nlines = sum(cnt for *_, cnt in lines)
    print(f"\n=== n = {n}: {len(cls)} nodes ({len(sc)} self-converse), {nlines} lines ({nblue} blue); fibre sizes = H/|Aut|: {ok_fibre}; THM-643 parities: {parity_ok}")
    # S classification
    S = Counter()
    for c1, c2, b, cnt in lines:
        key = tuple(sorted(('SC' if c in sc else 'NS') for c in (c1, c2))); S[(key, 'blue' if b else 'black')] += cnt
    print("  S (self-converse): ", dict(S), f"| black SC-NS lines {S[(('NS','SC'),'black')]} vs #NS nodes {len(cls)-len(sc)}: parity law {S[(('NS','SC'),'black')] % 2 == (len(cls)-len(sc)) % 2}")
    # K classification (closed end always strong)
    K = Counter()
    for c1, c2, b, cnt in lines:
        K[(cls[c1]["k"], cls[c2]["k"], 'blue' if b else 'black')] += cnt
    print("  K (components of closed, open end):", sorted(K.items()))
    # V classification
    V = Counter()
    for c1, c2, b, cnt in lines:
        key = tuple(sorted((colour(cls[c1]["H"]), colour(cls[c2]["H"])))); V[(key, 'blue' if b else 'black')] += cnt
    vnodes = Counter(colour(cls[c]["H"]) for c in cls)
    print("  V (eternal colour of H): nodes", dict(vnodes), "| lines", dict(V))
    # F classification
    FC = {c for c in cls if born[cls[c]["H"]] == n}
    F = Counter()
    for c1, c2, b, cnt in lines:
        key = tuple(sorted(('FC' if c in FC else 'old') for c in (c1, c2))); F[(key, 'blue' if b else 'black')] += cnt
    print(f"  F (first creators): {len(FC)} first-creator nodes creating {len({cls[c]['H'] for c in FC})} values | lines", dict(F))
    # multiplicities
    distinct = [(k, v) for k, v in pair_mult.items() if k[0] != k[1]]
    mx = max(v for k, v in distinct); arg = [k for k, v in distinct if v == mx]
    selfl = sum(v for k, v in pair_mult.items() if k[0] == k[1])
    print(f"  multiplicity: max lines between two distinct nodes {mx} at pairs {[(cls[a]['H'], cls[b]['H'], cls[a]['sc'], cls[b]['sc']) for a, b in arg][:4]}; self-lines (both ends one node) {selfl} (blue {sum(v for k, v in pair_blue.items() if k[0] == k[1])})")
    bound_ok = all(v <= min(fib[a], fib[b]) for (a, b), v in distinct)
    nbrs = defaultdict(set)
    for (a, b), v in pair_mult.items(): nbrs[a].add(b); nbrs[b].add(a)
    print(f"  multiplicity <= min fibre size on every pair: {bound_ok}; nodes whose distinct neighbours = fibre size: {sum(1 for c in cls if len(nbrs[c] - {c}) + (1 if c in nbrs[c] else 0) == fib[c])} of {len(cls)}")
    # same-value lines and H relations
    same = sum(cnt for c1, c2, b, cnt in lines if cls[c1]["H"] == cls[c2]["H"])
    hc_lt = sum(cnt for c1, c2, b, cnt in lines if cls[c1]["H"] < cls[c2]["H"]); hc_gt = sum(cnt for c1, c2, b, cnt in lines if cls[c1]["H"] > cls[c2]["H"])
    print(f"  H: same-value lines {same}; closed end smaller {hc_lt}, larger {hc_gt}; open end transitive (H=1) {sum(cnt for c1,c2,b,cnt in lines if cls[c2]['H']==1)}")

# ---- merged metagraph (converse pairs identified, as in THM-4729/785/793) ----
print("\n==== merged nodes (converse pairs identified) ====")
for n in range(5, 9):
    D = data[n]; cls = {int(k): v for k, v in D["classes"].items()}
    if any(v.get("conv") is None for v in cls.values()): print(f"n={n}: no converse ids in this file"); continue
    rep = {c: min(c, cls[c]["conv"]) for c in cls}
    merged = set(rep.values())
    pair_mult = Counter(); pair_blue = Counter(); selfl = 0; selfb = 0; ends = Counter()
    for c1, c2, b, cnt in D["lines"]:
        a, d = rep[c1], rep[c2]
        key = (min(a, d), max(a, d)); pair_mult[key] += cnt
        if b: pair_blue[key] += cnt
        ends[a] += cnt; ends[d] += cnt
    distinct = {k: v for k, v in pair_mult.items() if k[0] != k[1]}
    mx = max(distinct.values()); arg = [k for k, v in distinct.items() if v == mx]
    mxb = max((v for k, v in pair_blue.items() if k[0] != k[1]), default=0)
    argb = [k for k, v in pair_blue.items() if v == mxb and k[0] != k[1]]
    self_lines = sum(v for k, v in pair_mult.items() if k[0] == k[1]); self_blue = sum(v for k, v in pair_blue.items() if k[0] == k[1])
    sc_nodes = sum(1 for c in merged if cls[c]["sc"])
    print(f"n={n}: {len(merged)} merged nodes ({sc_nodes} self-converse, {len(merged)-sc_nodes} converse pairs); max lines between two distinct merged nodes {mx} at {[(cls[a]['H'], cls[b]['H']) for a, b in arg][:3]}; max blue lines between distinct nodes {mxb} at {[(cls[a]['H'], cls[b]['H']) for a, b in argb][:3]}; self-lines {self_lines} (blue {self_blue}); fibre of a merged pair node = 2 tc: {all(ends[c] == (1 if cls[c]['sc'] else 2) * cls[c]['H'] // cls[c]['aut'] for c in merged)}")
    # distinct neighbours vs fibre
    nb = defaultdict(set)
    for (a, d), v in pair_mult.items(): nb[a].add(d); nb[d].add(a)
    deg = sorted(((len(nb[c]), ends[c]) for c in merged), reverse=True)[:3]
    print(f"      largest neighbourhoods (distinct neighbours, fibre size): {deg}")
